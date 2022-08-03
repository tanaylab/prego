#include "port.h"
BASE_CC_FILE
#include "GenomeSeqIntervSet.h"
#include "strutil.h"

void GenomeSeqIntervSet::read(ifstream &tab, int marg, int generic_format)
{
	string chrom;
	int from, to, strand;

	int lcount = 0;
	vector<string> fields;
	
	tab >> chrom;
	while(tab) {
		if(!generic_format) {
    		tab >> from >> to >> strand;
		} else {
			split_line(tab, fields, '\t');
			ASSERT(fields.size() >= 3, "Bad line in intervals table line " << lcount);
			from = atoi(fields[1].c_str());
			to = atoi(fields[2].c_str());
			strand = 1;
			lcount++;
		}
		int newid = m_Intervs.size();
		m_interv_map[pair<string, int>(chrom, from)] = newid;
		ASSERT(to < m_Genome.chrom_size(chrom), "Bad interval, " << chrom << " " << from << " " << to << " max larger than chrom size");
		ASSERT(from < to, "Bad interval, " << chrom << " " << from << " " << to << " min > max");
		ASSERT(from >= 0, "Bad interval, " << chrom << " " << from << " " << to << " min < 0");
		from -= marg;
		if(from < 0) {
			from = 0;
		}
		to += marg;
		if(to >= m_Genome.chrom_size(chrom)) {
			to = m_Genome.chrom_size(chrom);
		}
		add_locus(chrom, from, to, strand);
		tab >> chrom;
	}
}

GenomeSeqIntervSet::~GenomeSeqIntervSet()
{
	for(int i = 0; i < m_seq_repository.size(); i++) {
		delete m_seq_repository[i];
	}
}

void GenomeSeqIntervSet::write(ostream &tab)
{
	for(vector<GenomeSeqInterval>::const_iterator i = m_Intervs.begin(); i != m_Intervs.end(); i++) {
		tab << m_Genome.chrom_name(i->m_chrom_id) 
		<< "\t" << i->m_from 
		<< "\t" << i->m_to
		<< "\t" << i->m_strand << "\n";
	}
}

void GenomeSeqIntervSet::init_from_interv(const GenomeIntervSet &intervs)
{
	for(int id = 0; id < intervs.max_interv_id(); id++) {
			const GenomeInterval &itr = intervs.interv(id);
			add_locus(itr.m_chrom, itr.m_from, itr.m_to, itr.m_strand);
	}
	
}

void GenomeSeqIntervSet::add_locus(const string &chrom, int from, int to, int strand)
{
	int id = m_Intervs.size();
	m_Intervs.resize(id + 1);
	GenomeSeqInterval &interv = m_Intervs[id];
	interv.m_from = from;
	interv.m_to = to;
	interv.m_strand = strand;
	interv.m_chrom_id = m_Genome.chrom_code(chrom);

	if(m_project_align) {
		//figure out the aln coordinate
		string::const_iterator alnseq = m_project_align->get_seq(m_project_spid, chrom, 0);
		
		int aln_from = m_project_align->get_alnpos_by_refgenome(from);
		int aln_to = m_project_align->get_alnpos_by_refgenome(to);
		cerr << "mapped " << from << " " << to << " to aln " << aln_from << " " << aln_to << endl;

		//filter gaps?
		string *interv_seq = new string;
		if(m_project_seglen != -1) {
    		int maxpos = m_project_align->get_max_pos(chrom);
			
			int aln_center = m_project_align->get_alnpos_by_refgenome(int(from+to)/2);
			int count = 0;
			int half_size = int(m_project_seglen/2);
    		string::const_iterator src_i = alnseq + aln_center;
			while(count < half_size && src_i >= alnseq) {
    			if(*src_i != '-' && *src_i != 'N' && *src_i != 'n') {
    				count++;
    			}
    			src_i--;
			}
			count = 0;
			while(count < m_project_seglen && src_i < alnseq + maxpos) {
    			if(*src_i != '-' && *src_i != 'N' && *src_i != 'n') {
    				if(m_project_uppercase) {
                		interv_seq->push_back(toupper(*src_i));
    				} else {
                		interv_seq->push_back(*src_i);
    				}
    			}
    			count++;
    			src_i++;
			}
    	} else {
    		for(string::const_iterator src_i = alnseq + aln_from; src_i < alnseq + aln_to; src_i++) {
    			if(*src_i != '-' && *src_i != 'N' && *src_i != 'n') {
    				char tc = toupper(*src_i);
    				ASSERT(tc == 'A' || tc == 'C' || tc == 'G' || tc == 'T', "Bad character " << tc << " when initializing interval " << id << " chr "<< chrom << " " << from << " " << to);
    				if(m_project_uppercase) {
                		interv_seq->push_back(toupper(*src_i));
    				} else {
                		interv_seq->push_back(*src_i);
    				}
    			}
    		}
    	}
		m_seq_repository.push_back(interv_seq);
    	interv.m_seq = m_seq_repository.back()->begin();
    	interv.m_seq_end = m_seq_repository.back()->end();
    	
    	interv.m_chr_start = m_seq_repository.back()->begin();
    	interv.m_chr_end = m_seq_repository.back()->end();
	} else {
    	interv.m_seq = m_Genome.locus(interv.m_chrom_id, from);
    	interv.m_seq_end = interv.m_seq + (to - from);
    	
    	interv.m_chr_start = m_Genome.chrom_seq(chrom).begin();
    	interv.m_chr_end = m_Genome.chrom_seq(chrom).end();
	}
}

void GenomeSeqIntervSet::shift_intervs(int shift)
{
	ASSERT(m_seq_repository.size() == 0, "Cannot shift intervals in projected alignment settings");

	for(vector<GenomeSeqInterval>::iterator i = m_Intervs.begin(); i != m_Intervs.end(); i++) {
		if(i->m_seq  + shift > i->m_chr_start && i->m_chr_end > i->m_seq_end + shift) {
    		i->m_seq += shift;
    		i->m_seq_end += shift;
    		i->m_from += shift;
    		i->m_to += shift;
		}
	}
}
