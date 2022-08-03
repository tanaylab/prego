#ifndef biodata_MafParse_h
#define biodata_MafParse_h 1

/*=================================================
 *
=================================================*/

#include <vector>
#include <string>

class PhyloTree;

class MafParse {

protected:

	istream &m_file;

	vector<string> m_seq;
	vector<int> m_pos;
	vector<int> m_max_pos;
	vector<int> m_src_size;
	vector<char> m_strand;
	vector<string> m_desc;

public:

	int get_num_of_seq() {
		return(m_seq.size());
	}
	const string &get_desc(int i) {
		return(m_desc[i]);
	}

	int get_pos(int i) {
		return(m_pos[i]);
	}

	int get_max_pos(int i) {
		return(m_max_pos[i]);
	}

	int get_src_size(int i) {
		return(m_src_size[i]);
	}

	char get_strand(int i) {
		return(m_strand[i]);
	}

	const string &get_seq(int i) {
		return(m_seq[i]);
	}

	MafParse(ifstream &in) :
		m_file(in)
	{}

	bool next();

	int get_species_id(const PhyloTree &phylo, int i);
};

#endif // MafParse 
