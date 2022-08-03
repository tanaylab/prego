#ifndef GENOMEINTERVSET_H_
#define GENOMEINTERVSET_H_

#include <vector>
#include <string>
#include <map>
#include <list>
#include "GenomeSequence.h"

struct GenomeInterval {

public:
	string m_chrom;

	int m_from;
	int m_to;
	int m_strand;
};

class GenomeIntervSet {

public:

	void read(ifstream &tab, int nostrand=0);
	void write(ostream &tab);

	void init_from_key(ifstream &key);
    void read_bed(istream &bed);

	const GenomeInterval &interv(int id) const {
		return(m_Intervs[id]);
	}

	const vector<GenomeInterval> &intervs() const {
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

	// returns id of the added interval
    int add_interv(const string &chrom, int from, int to, int strand);

    int add_interv(const GenomeInterval &val) { return add_interv(val.m_chrom, val.m_from, val.m_to, val.m_strand); }

 //initialize by merging overlaps and soring (by chrom,start coord) an input set of intervals
    void merge_and_sort(const GenomeIntervSet &other);

    void gen_coverage_stat(const GenomeIntervSet &ref_probes, int directed_margins, vector<float> &probe_coverage);

    void reset() {
    	m_Intervs.resize(0);
    	m_interv_map.clear();
    }

private:

	vector<GenomeInterval> m_Intervs;

	map<pair<string, int>, int> m_interv_map;
};

#endif /*GENOMEINTERVSET_H_*/
