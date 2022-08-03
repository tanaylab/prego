#ifndef biodata_PhyloTree_h
#define biodata_PhyloTree_h 1

/*=================================================
=================================================*/

#include <vector>
#include <string>
#include <map>

class PhyloTree {

protected:

	int m_root_id;

	vector<int> m_pred;
	vector<int> m_right;
	vector<int> m_left;

	vector<float> m_left_time;
	vector<float> m_right_time;
	vector<float> m_pred_time;

	vector<bool> m_is_inner;

	vector<string> m_node_name;

	int m_max_sp_node;
	vector<vector<int> > m_up_partit;
	vector<int> m_up_part_size;

//	vector<vector<int> > m_uppart_multid;

	map<string, int> m_spname_map;

public:

	int get_root_id() const {
		return(m_root_id);
	}

	int get_max_node_id() const {
		return(m_is_inner.size());
	}
	int get_max_sp_node_id() const {
		return(m_max_sp_node);
	}

	int get_pred(int nid) const {
		return(m_pred[nid]);
	}
	int get_right(int nid) const {
		return(m_right[nid]);
	}
	int get_left(int nid) const {
		return(m_left[nid]);
	}

	float get_right_time(int nid) const {
		return(m_right_time[nid]);
	}
	float get_left_time(int nid) const {
		return(m_left_time[nid]);
	}
	float get_pred_time(int nid) const {
		return(m_pred_time[nid]);
	}

	bool is_inner(int nid) const {
		return(m_is_inner[nid]);
	}
	bool is_leaf(int nid) const {
		return(!m_is_inner[nid]);
	}

//	int comp_uppart_multid(int lin, int id) const;
//	int get_uppart_multid(int lin, int id) const {
//		return(m_uppart_multid[lin][id]);
//	}
	int get_uppart_max_multid(int lin) const {
		return(m_up_part_size[lin]);
	}
	int is_in_uppart(int lin, int node) const {
		return(m_up_partit[lin][node]);
	}

	void read(istream &tab);

	const string &get_name(int sp) const {
		return(m_node_name[sp]);
	}

	int get_sp_by_name(const string &nm) const {
		map<string, int>::const_iterator i = m_spname_map.find(nm);
		if(i == m_spname_map.end()) {
			return(-1);
		} else {
			return(i->second);
		}
	}
};

#endif // PhyloTree
