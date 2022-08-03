#ifndef GENOMESEQINTERVSET_H_
#define GENOMESEQINTERVSET_H_

#include <vector>
#include <string>
#include <map>
#include "GenomeSequence.h"
#include "GenomeIntervSet.h"
#include "GenomesAlignment.h"

struct GenomeSeqInterval {
	
public:
	int m_chrom_id;
	
	int m_from;
	int m_to;
	int m_strand;
	
	string::const_iterator m_seq;
	string::const_iterator m_seq_end;
	
	string::const_iterator m_chr_start;
	string::const_iterator m_chr_end;
};

class GenomeSeqIntervSet {
	
public:

	GenomeSeqIntervSet(const GenomeSequence &gn) :
		m_Genome(gn),
		m_project_align(0),
		m_project_spid(-1)
	{}
	
	~GenomeSeqIntervSet();
	
	void set_alignment_projection(GenomesAlignment *align, int spid, int seglen, int uppercase = 0) {
		m_project_align = align;
		m_project_spid = spid;
		m_project_uppercase = uppercase;
		m_project_seglen = seglen;
	}
	
	const GenomeSequence &genome() {
		return(m_Genome);
	}
	
	void read(ifstream &tab, int marg = 0, int generic_format = 0);
	void write(ostream &tab);
	
	void init_from_interv(const GenomeIntervSet &intervs);

	void add_locus(const string &chrom, int from, int to, int strand);
	
	const GenomeSeqInterval &interv(int id) const {
		return(m_Intervs[id]);
	}

	int interv_start_coord(int interv_id) const {
		return(m_Intervs[interv_id].m_from);
	}
	int interv_end_coord(int interv_id) const {
		return(m_Intervs[interv_id].m_to);
	}
	string::const_iterator interv_start(int interv_id) const {
		return(m_Intervs[interv_id].m_seq);
	}
	string::const_iterator interv_end(int interv_id) const {
		return(m_Intervs[interv_id].m_seq_end);
	}
	const string &interv_chr(int interv_id) const {
		return(m_Genome.chrom_name(m_Intervs[interv_id].m_chrom_id));
	}
	string::const_iterator chr_start(int interv_id) const {
		return(m_Intervs[interv_id].m_chr_start);
	}
	string::const_iterator chr_end(int interv_id) const {
		return(m_Intervs[interv_id].m_chr_end);
	}
	int interv_strand(int interv_id) const {
		return(m_Intervs[interv_id].m_strand);
	}
	
	const vector<GenomeSeqInterval> &intervs() const {
		return(m_Intervs);
	}
	
	int max_interv_id() const {
		return(m_Intervs.size());
	}

	int interv_id(const string &chrom, int coord) const {
		pair<string, int> key(chrom, coord);
		map<pair<string, int>, int>::const_iterator i = m_interv_map.find(key);
		if(i == m_interv_map.end()) {
			return(-1);
		} else {
			return(i->second);
		}
	}
    void shift_intervs(int shift);
	
private:
	
	const GenomeSequence &m_Genome;
	
	vector<GenomeSeqInterval> m_Intervs;

	map<pair<string, int>, int> m_interv_map;
	
	GenomesAlignment *m_project_align;
	int m_project_spid;
	int m_project_uppercase;
	int m_project_seglen;
	
	vector<string *> m_seq_repository;
	
};
#endif /*GENOMESEQINTERVSET_H_*/
