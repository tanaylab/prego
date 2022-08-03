#ifndef GENOMESEQUENCE_H_
#define GENOMESEQUENCE_H_

#include <vector>
#include <string>
#include <map>

class GenomeSequence {
	
public:

	GenomeSequence() :
		m_capitalize(0)
	{}
	
	void capitalize(int cap = 1) {
		m_capitalize = cap;
	}
	
	void read_key_fasta(istream &key_file);
	
	const string &chrom_seq(const string &chrom) const {
		int id = chrom_code(chrom);
		if(id == -1) {
			return(m_null);
		} else {
    		return(m_Seq[id]);
		}
	}
	int chrom_size(int id) const {
		return(m_Seq[id].size());
	}
	
	int chrom_size(const string &chrom) const {
		int id = chrom_code(chrom);
		if(id == -1) {
			return(-1);
		} else {
    		return(m_Seq[id].size());
		}
	}

	int chrom_code(const string &chrom) const {
		map<string, int>::const_iterator i= m_ChromId.find(chrom);
		if(i == m_ChromId.end()) {
			return(-1);
		} else {
    		return(i->second);
		}
	}

	string::const_iterator locus(const string &chrom, int coord) const {
		int id = chrom_code(chrom);
		if(id == -1) {
			return(m_null.begin());
		} else {
			return(locus(id, coord));
		}
	}
	string::const_iterator locus(int id, int coord) const {
		return(m_Seq[id].begin() + coord);
	}
	
	const string &chrom_name(int cid) const {
		return(m_ChromName[cid]);
	}
	
	int max_chrom_id() const {
		return(m_ChromName.size());
	}
	
private:

	void read_chrom_fasta(int id);

	string m_null;
	
	int m_capitalize;
	
	vector<string> m_ChromFn;
	vector<string> m_ChromName;
	
	map<string, int> m_ChromId;
	
	vector<string> m_Seq;
	
};
#endif /*GENOMESEQUENCE_H_*/
