#include "port.h"
BASE_CC_FILE
#include "GenomesAlignment.h"

GAlignsMafsCache::GAlignsMafsCache(const PhyloTree &phy, ifstream &key_tab, int max_miss, int spat_bin) :
	GenomesAlignment(phy),
	m_seq(phy.get_max_sp_node_id()),
	m_max_missing_align(max_miss),
	m_spat_bin(spat_bin),
	m_frag_to_genome(phy.get_max_node_id())
{
	m_base_ref_coord = -1;
	m_max_ref_coord = -1;

	m_max_global_aln_coord=-1;
	m_base_global_aln_coord=-1;
	m_use_ref_genome_coords=1;
//	m_base_aln_coor=1;

	m_cur_chrom = "";
	string chrom, fn;
	key_tab >> chrom;
	while(key_tab) {
		key_tab >> fn;
		m_chrom_maf_fn[chrom] = fn;
        m_chroms.push_back(chrom);

		key_tab >> chrom;
	}
}


void GAlignsMafsCache::set_focus_interval(const string &chr, int fr, int to)
{
	m_cur_chrom = chr;
	if (m_use_ref_genome_coords){
		m_base_ref_coord = fr;
		m_max_ref_coord = to;
		//Rcpp::Rcerr << "set interv " << m_base_ref_coord <<  " " << m_max_ref_coord << endl;
	} else {
		m_base_global_aln_coord = fr;
		m_max_global_aln_coord = to;
	}

	read_chrom(chr);
}

void GAlignsMafsCache::reset()
{
		for(unsigned int spi = 0; spi < m_seq.size(); spi++) {
			m_seq[spi].resize(0);
		}
    	m_frag_start.resize(0);
    	m_frag_end.resize(0);
    	m_alnpos_to_ref_genome_bins.resize(0);
    	m_ref_genome_bins.resize(0);
    	m_alnpos_frag.resize(0);
    	fill(m_frag_to_genome.begin(), m_frag_to_genome.end(), vector<AlnFrag>());
}

string::const_iterator GAlignsMafsCache::get_seq(int sp, const string &chrom, int fr)
{
	if(m_base_ref_coord != -1) {
    	ASSERT(chrom == m_cur_chrom,
			"Trying to acces alignment for chrom " << chrom << " when the focus was " << m_cur_chrom);
    } else {
    	if(chrom != m_cur_chrom) {
    		ASSERT(m_base_ref_coord == -1, "Tried to read seq from chrom "<< chrom
    				<< " when current chromo is " << m_cur_chrom
    				<< " and AlignMafCache is focused on interval - call set_base_ref_coord before");
    		//reinit chrom
    		read_chrom(chrom);
    		m_cur_chrom = chrom;
    	}
    }

	return(m_seq[sp].begin() + fr);
}

int GAlignsMafsCache::get_num_active_sp (MafParse& maf, int max_spi){
	int active_sp = 0;
	for(int spi = 0; spi < max_spi; spi++) {	//look for the species in the entry that is part of our phylo
		int sp = maf.get_species_id(m_phylo, spi);
		if(sp != -1) {
			active_sp++;
		}
	}
	return active_sp;
}

int GAlignsMafsCache::get_ref_sp_id(MafParse& maf , int max_spi){
	int ref_spid = -1;
	for(int spi = 0; spi < max_spi; spi++) {			//look for the id (in the entry) of the reference species
		const string &msp = maf.get_desc(spi);
		if(msp.find(m_phylo.get_name(0)) == 0) {
			ref_spid = spi;
		}
	}
	ASSERT(ref_spid != -1, "Cannot find ref species in maf, first entry " << maf.get_desc(0) << " " << maf.get_pos(0));
	return ref_spid;
}

void GAlignsMafsCache::read_chrom(const string &chrom)
{
	if (m_use_ref_genome_coords){
		read_chrom_ref_genome_coord(chrom);
	} else {
		read_chrom_global_aln_coord(chrom);
	}
}

void GAlignsMafsCache::read_chrom_ref_genome_coord(const string &chrom){
//	Rcpp::Rcerr << "read_chrom_ref_genome_coord" << endl;
	reset();

	ASSERT(m_chrom_maf_fn.find(chrom) != m_chrom_maf_fn.end(),
			"Cannot find maf fn from requested chrom " << chrom);

	ifstream maftab(m_chrom_maf_fn[chrom].c_str());

	ASSERT(maftab, "Cannot open maf file " << m_chrom_maf_fn[chrom]);

	MafParse maf(maftab);

	int max_sp = m_phylo.get_max_sp_node_id();

	m_max_chr_coord = 0;
	m_base_aln_coor = 0;

	//Rcpp::Rcerr << "will read chrom, base " << m_base_ref_coord << endl;
//initialize
	int mfi = 0;


	while(maf.next()) {	//get next maf parse entry (group of aligned sequences)


		int max_spi = maf.get_num_of_seq();		//the number of sequences present in this entry
		int active_sp = get_num_active_sp(maf , max_spi);

//		this lines are in comment in order to make global alignment coordinates the same for all max_missing_alignment values
//		instead the check is being made afterwards when the sequence is nucleotides (and not gaps) are replaced with N's if there are not enough species in the entry
//		if(m_max_missing_align < m_seq.size() - active_sp) {//if we have too few active speceis, give up on the entry
//			continue;
//		}
		int ref_spid = get_ref_sp_id(maf , max_spi);

		if(m_base_ref_coord != -1 && maf.get_pos(ref_spid) > m_max_ref_coord) {	//if the position in the ref is larger than our current focus interval - terminate
			if (m_max_chr_coord < m_max_ref_coord){
				int num_of_N = m_max_ref_coord - m_max_chr_coord;
				for (int sp=0; sp<max_sp; sp++){
					m_seq[sp].append(num_of_N , 'N');
				}
			}
			break;	//this means that we assume the maf to be sorted!
		}


		int base_chr_coord = maf.get_pos(ref_spid);	//this is the coordinate of the entry start in the reference genome coordinate system

		if(m_base_ref_coord != -1 && maf.get_max_pos(ref_spid) < m_base_ref_coord) {	//if the max position of the current entry is smaller than the begining of our windows, skip it
			m_base_aln_coor += (base_chr_coord - m_max_chr_coord);
			m_max_chr_coord = maf.get_max_pos(ref_spid);
			m_base_aln_coor += maf.get_seq(0).size();
			continue;
		}

		//Rcpp::Rcerr << "got maf entry, ref " << ref_spid << " pos ref " << maf.get_pos(ref_spid) << " gap " << -m_max_chr_coord-maf.get_pos(ref_spid) << endl;
		if(base_chr_coord > m_max_chr_coord) {		//if this entry starts after more than 0 character from the ending of the last entry, add padding
			int gap_size = base_chr_coord - m_max_chr_coord;
			int num_out_of_focus =  max(0 , m_base_ref_coord - m_max_chr_coord); // number of nucs before the interval
			num_out_of_focus = min(num_out_of_focus , gap_size ); // if the current entry starts before interval need bound number out_of_focus
			int num_in_focus  = gap_size - num_out_of_focus;
			for(int sp = 0; sp < max_sp; sp++) { // padding according to the number in focus
				m_seq[sp].append(num_in_focus, 'N');
			}
			m_base_aln_coor+= (num_out_of_focus); // enlarging base_aln_coor according to the number out of focus
		}

		int offset = 0;
		if(base_chr_coord < m_base_ref_coord) {	//if this entry starts before the beginning of the current focus interval
			string::const_iterator s = maf.get_seq(ref_spid).begin();
			while(base_chr_coord < m_base_ref_coord) {
				if(*s != '-') {
					base_chr_coord++;		//we advance the ref coordinate if this is not a gap
				}
				offset++;
				s++;
			}
			while (*s == '-'){
				s++;
				offset++;
			}
			m_base_aln_coor+=offset;
			//after finishing this, base_chr_coordinate equals m_base_ref_coord which is the begining of the focus interval
			//and the offest variale equals the number of characters (in alignment space) that we should skip in the begining of the entry
		}

		int suffix_offset = 0;
		if(m_max_ref_coord != -1 && maf.get_max_pos(ref_spid) > m_max_ref_coord) {	//if the ending of the focus interal is smaller than the ending of the current entry
    		string::const_reverse_iterator s = maf.get_seq(ref_spid).rbegin(); //note the use of a reverse iterator
    		int max_chr_coord = maf.get_max_pos(ref_spid);	//start from the max position
    		while(max_chr_coord > m_max_ref_coord) {
    			if(*s != '-') {
    				max_chr_coord--;
    			}
    			suffix_offset++;
    			s++;
    		}
    		//after finishing this, max_chr_coord equals m_max_ref_coord
    		//and suffix_offest euqals the number of of characters we should prune from the end of the max entry
		}

		ASSERT(base_chr_coord >= m_max_chr_coord,
			"Maf file should be sorted by the first species coordinate, but entry "
			<< base_chr_coord << " < cur max " << m_max_chr_coord);

		m_frag_start.push_back(m_seq[0].size());
		append_seq_by_maf_entry(maf, max_spi, active_sp, offset, suffix_offset);
		m_frag_end.push_back(m_seq[0].size());
	//Reverse all str for debugging
#ifdef REV_DBG
		for(int sp = 0; sp < m_seq.size(); sp++) {
			reverse(m_seq[sp].begin() + base_chr_coord, m_seq[sp].end());
		}
#endif

		m_max_chr_coord = maf.get_max_pos(ref_spid);
		mfi++;
	}
//	Rcpp::Rcerr << "done reading maf, seq at 0 internal len is " << m_seq[0].length() << endl;
	create_maps();
}

void GAlignsMafsCache::read_chrom_global_aln_coord(const string &chrom){
//	Rcpp::Rcerr << "read_chrom_global_aln_coord" << endl;
	reset();

	ASSERT(m_chrom_maf_fn.find(chrom) != m_chrom_maf_fn.end(),
			"Cannot find maf fn from requested chrom " << chrom);

	ifstream maftab(m_chrom_maf_fn[chrom].c_str());
	ASSERT(maftab, "Cannot open maf file " << m_chrom_maf_fn[chrom]);
	MafParse maf(maftab);
	int max_sp = m_phylo.get_max_sp_node_id();
	m_max_chr_coord = 0;
	if (m_base_global_aln_coord == -1){
		m_base_ref_coord = -1;
	}
	//Rcpp::Rcerr << "will read chrom, base " << m_base_ref_coord << endl;
//initialize
	int mfi = 0;
	m_curr_global_aln_coord = 0;
	bool is_seq_started = false;
	while(maf.next()) {	//get next maf parse entry (group of aligned sequences)

		int max_spi = maf.get_num_of_seq();		//the number of sequences present in this entry
		int active_sp = get_num_active_sp(maf , max_spi);

//		this lines are in comment in order to make global alignment coordinates the same for all max_missing_alignment values
//		instead the check is being made afterwards when the sequence is nucleotides (and not gaps) are replaced with N's if there are not enough species in the entry
//		if(m_max_missing_align < m_seq.size() - active_sp) {//if we have too few active speceis, give up on the entry
//			continue;
//		}
		int ref_spid = get_ref_sp_id(maf , max_spi);
		int base_chr_coord = maf.get_pos(ref_spid);	//this is the coordinate of the entry start in the reference genome coordinate system
		if(m_base_global_aln_coord != -1 && m_curr_global_aln_coord +base_chr_coord - m_max_chr_coord> m_max_global_aln_coord) {	//if the position in the ref is larger than our current focus interval - terminate
			if (m_curr_global_aln_coord < m_max_global_aln_coord){
				int num_of_N = m_max_global_aln_coord - m_curr_global_aln_coord;
				for (int sp=0; sp<max_sp; sp++){
					m_seq[sp].append(num_of_N , 'N');
				}
			}
			break;	//this means that we assume the maf to be sorted!
		}

		if(m_base_global_aln_coord != -1 && m_curr_global_aln_coord + base_chr_coord - m_max_chr_coord + maf.get_seq(0).size() < m_base_global_aln_coord) {	//if the max global position of the current entry is smaller than the begining of our windows, skip it
			m_curr_global_aln_coord += (base_chr_coord - m_max_chr_coord);
			m_max_chr_coord = maf.get_max_pos(ref_spid);
			m_curr_global_aln_coord += maf.get_seq(0).size();
			continue;
		}

		if(base_chr_coord > m_max_chr_coord) {		//if this entry starts after more than 0 character from the ending of the last entry, add padding
			int gap_size = base_chr_coord - m_max_chr_coord;
			int num_out_of_focus =  max(0 , m_base_global_aln_coord - m_curr_global_aln_coord); // number of nucs before the interval
			num_out_of_focus = min(num_out_of_focus , gap_size ); // if the current entry starts before interval need bound number out_of_focus
			int num_in_focus  = gap_size - num_out_of_focus;
			if (m_base_global_aln_coord!= -1 && !is_seq_started && num_out_of_focus > 0 && num_in_focus > 0){
				m_base_ref_coord = m_max_chr_coord + num_out_of_focus;
				is_seq_started = true;
			}
			for(int sp = 0; sp < max_sp; sp++) { // padding according to the number in focus
				m_seq[sp].append(num_in_focus, 'N');
			}
			m_curr_global_aln_coord += gap_size; // enlarging base_aln_coor according to the number out of focus
		}
		//Rcpp::Rcerr << "got maf entry, ref " << ref_spid << " pos ref " << maf.get_pos(ref_spid) << " gap " << -m_max_chr_coord-maf.get_pos(ref_spid) << endl;

		int offset = ( m_base_global_aln_coord == -1 ?  0 :  max(0 , m_base_global_aln_coord - m_curr_global_aln_coord) );
		if (m_base_global_aln_coord!= -1 && !is_seq_started){
			m_base_ref_coord = base_chr_coord;
			string::const_iterator s = maf.get_seq(ref_spid).begin();
			int k=0;
			while (k<offset ){ // first advancing in the entry until the first nucleotide
				if (*s!='-'){
					m_base_ref_coord++;
				}
				s++;
				k++;
				is_seq_started = true;
			}
		}
		int suffix_offset = ( m_max_global_aln_coord == -1 ? 0 : m_curr_global_aln_coord + maf.get_seq(0).size() - m_max_global_aln_coord );
		suffix_offset  = max(0 ,  suffix_offset);
		ASSERT(base_chr_coord >= m_max_chr_coord,
			"Maf file should be sorted by the first species coordinate, but entry "
			<< base_chr_coord << " < cur max " << m_max_chr_coord);

		m_frag_start.push_back(m_seq[0].size());
		append_seq_by_maf_entry(maf, max_spi, active_sp, offset, suffix_offset);
		m_frag_end.push_back(m_seq[0].size());
	//Reverse all str for debugging
#ifdef REV_DBG
		for(int sp = 0; sp < m_seq.size(); sp++) {
			reverse(m_seq[sp].begin() + base_chr_coord, m_seq[sp].end());
		}
#endif
		m_curr_global_aln_coord += (maf.get_seq(0).size() - suffix_offset  );
		m_max_chr_coord = maf.get_max_pos(ref_spid);
		mfi++;
	}
	m_base_aln_coor = m_base_global_aln_coord;
	create_maps();
}

void GAlignsMafsCache::append_seq_by_maf_entry(MafParse& maf , int max_spi , int active_sp, int offset , int suffix_offset){
	for(int spi = 0; spi < max_spi; spi++) {
		int sp = maf.get_species_id(m_phylo, spi);
		if(sp == -1) {
			//ASSERT(sp >= 0, "Cannot match maf species " << maf.get_desc(spi) << " in phylo table");
			continue;
		}
		m_frag_to_genome[sp].push_back(AlnFrag(maf.get_desc(spi), maf.get_pos(spi), maf.get_strand(spi)));

		if(m_max_missing_align < m_seq.size() - active_sp) {//if we have too few active speceis, put N's instead of the original alignment
			//int len = maf.get_seq(0).size() - offset - suffix_offset ;
			//m_seq[sp].append(len , 'N');
			append_seq_of_N(m_seq[sp] , maf.get_seq(spi).substr(offset , maf.get_seq(spi).size() - offset - suffix_offset) );
		} else {
			if(offset > 0 || suffix_offset > 0) {
				m_seq[sp].append(maf.get_seq(spi), offset, maf.get_seq(spi).size() - offset - suffix_offset);
			} else {
				m_seq[sp] += maf.get_seq(spi);
			}
		}
	}

// appending N's to species that are not in the current focus interval
	unsigned int max_length = m_seq[0].length();
	unsigned int max_sp_length=0;
	for(unsigned int sp = 1; sp < m_seq.size(); sp++) {
		if (m_seq[sp].length() > max_length ){
			max_length = m_seq[sp].length();
			max_sp_length = sp;
		}
	}
	ASSERT (max_sp_length == 0 , "ref species is not in interval");
	for(unsigned int sp = 0; sp < m_seq.size(); sp++) {
		if(m_seq[max_sp_length].length() != m_seq[sp].length()) {
			m_frag_to_genome[sp].push_back(AlnFrag("null", -1, 1));
			m_seq[sp].append(m_seq[max_sp_length].length() - m_seq[sp].length(), 'N');
		}
	}
}

void GAlignsMafsCache::create_maps(){
	int delta = (m_base_ref_coord == -1 ? 1 :0); // for binning the base_ref_coord should not start at -1
	m_base_ref_coord += delta;
	int coord = m_base_ref_coord-1;
	string::const_iterator s = m_seq[0].begin();
	string::const_iterator max_s = m_seq[0].end();
	int aln_coord = 0;
	m_ref_genome_bins.resize(int((max_s-s)/m_spat_bin)+1);
	while(s != max_s) {
		if(*s != '-') {
			coord++;
		}
		if(aln_coord % m_spat_bin == 0) {
			m_ref_genome_bins[int(aln_coord/m_spat_bin)] = coord;
		}
		if((coord - m_base_ref_coord) % m_spat_bin == 0 && *s != '-') {
			int rbin = int((coord-m_base_ref_coord)/m_spat_bin);
			m_alnpos_to_ref_genome_bins.resize(rbin + 1);
			m_alnpos_to_ref_genome_bins[rbin] = aln_coord;
		}
		s++;
		aln_coord++;
	}
	m_base_ref_coord -= delta;
	
/*	int frag_i = 0;
	for(int aln_coord = 0; aln_coord < m_seq[0].size(); aln_coord += m_spat_bin) {
		while((frag_i + 1) < m_frag_start.size() && m_frag_start[frag_i+1] < aln_coord) {
			frag_i++;
		}
    	m_alnpos_frag[int(aln_coord/m_spat_bin)] = frag_i;
    }*/
}

void GAlignsMafsCache::append_seq_of_N (string& str, const string& to_append  ){
	string to_append_cpy = to_append;
	string::iterator s = to_append_cpy.begin();
	while (s != to_append_cpy.end()){
		if (*s != '-'){
			*s = 'N';
		}
		s++;
	}
	str.append(to_append_cpy);
}

int GAlignsMafsCache::get_refgenome_pos(int pos) const
{
	//Rcpp::Rcerr << "GAlignsMafsCache::get_refgenome_pos (" << pos << ")" << endl;
	//check to see that this if max align pos
	if (pos == m_seq[0].length()) {
		//Rcpp::Rcerr << "pos is last pos, returning " << m_max_ref_coord << endl;
		return(m_max_ref_coord);
	}

	//look for a seed
	int bin = int(pos/m_spat_bin);
	int base = bin*m_spat_bin;
	int refpos = m_ref_genome_bins[bin];

	string::const_iterator s = m_seq[0].begin() + base;
	while(base < pos) {
		s++;
		if(*s != '-') {
			refpos++;
		}
		base++;
	}
	return(refpos);
}

int GAlignsMafsCache::get_alnpos_by_refgenome(int refpos)
{
	//look for a seed
	int delta = (m_base_ref_coord == -1 ? 1 :0); // for binning the base ref coord should not start at -1
	m_base_ref_coord += delta;
	int bin = int((refpos-m_base_ref_coord)/m_spat_bin);
	int baseref = m_base_ref_coord + bin*m_spat_bin;
	int pos = m_alnpos_to_ref_genome_bins[bin];

	string::const_iterator s = m_seq[0].begin() + pos;
	ASSERT (*s != '-' , "mapping problem a reg_genome bin points to a gap " << bin * m_spat_bin << "\t" << pos );
	while(baseref < refpos) {
		pos++;
		s++;
		while(*s == '-') {
			pos++;
			s++;
		}
		baseref++;
	}
	m_base_ref_coord -= delta;
	return(pos);
}

int GAlignsMafsCache::step_align_coord(int fr_aln_coord, int fr_chr_coord, int to_chr_coord) const
{
	int chr_coord = fr_chr_coord;
	int aln_coord = fr_aln_coord;
	string::const_iterator seq = m_seq[0].begin() + fr_aln_coord;
	int max_aln_coord = m_seq[0].size();
	while(chr_coord < to_chr_coord && aln_coord < max_aln_coord) {
		if(*seq != '-') {
			chr_coord++;
		}
		seq++;
		aln_coord++;
	}
	return(aln_coord);
}
