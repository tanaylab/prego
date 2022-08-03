#ifndef GENOMEINTERVSUBSETS_H_
#define GENOMEINTERVSUBSETS_H_

#include "GenomeSeqIntervSet.h"
#include "BitVec.h"
#include <map>

class GenomeSeqIntervSubsets {
	
protected:

	GenomeSeqIntervSet &m_intervs;
	
	vector<ds_bitvec> m_subsets;
	
	vector<vector<pair<int, int> > > m_offsets;
	
	map<string, int> m_setname_to_id;
	
	vector<string> m_subset_name;
	
public:

	const GenomeSeqIntervSet &get_intervs() const {
		return(m_intervs);
	}
	
	GenomeSeqIntervSubsets(GenomeSeqIntervSet &intr) :
		m_intervs(intr)
	{}
	
	int max_subset_id() const {
		return(m_subsets.size());
	}
	
	const ds_bitvec &subset(int i) const {
		return(m_subsets[i]);
	}
	
	int subset_size(int i) const {
		return(m_subsets[i].on_bits());
	}
	
	void subset(int i, const ds_bitvec &set) {
		if(i >= m_subsets.size()) {
			m_subsets.resize(i + 1);
		}
		m_subsets[i] = set;
	}
	
	const string &subset_name(int i) {
		return(m_subset_name[i]);
	}

	int fr_offset(int setid, int intrvid) const {
		return(m_offsets[setid][intrvid].first);
	}
	int to_offset(int setid, int intrvid) const {
		return(m_offsets[setid][intrvid].second);
	}
	
	vector<pair<int, int> > *subset_offsets(int setid) {
		return(&(m_offsets[setid]));
	}
	
	void shift_intervs(int shift) { 
		m_intervs.shift_intervs(shift); 
	}
	
	void read(istream &sets, int generic_format = 0);
};

#endif /*GENOMEINTERVSUBSETS_H_*/
