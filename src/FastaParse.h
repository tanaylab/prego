#ifndef biodata_FastaParse_h
#define biodata_FastaParse_h 1

/*=================================================
 *
=================================================*/

#include <vector>
#include <string>

class FastaParse {

protected:

	ifstream m_file;

	string *m_seq;
	
	bool m_ext_storage;
	string m_desc;

	bool m_is_valid;

	bool m_capitalize;

private:

	vector<string> m_aux_fields;
public:

	const string &get_desc() {
		return(m_desc);
	}

	const string &get_seq() {
		return(*m_seq);
	}

	FastaParse(const string &fn, bool capit = false) :
		m_file(fn.c_str()),
		m_capitalize(capit)
	{
		m_seq = new string;
		m_ext_storage = false;
		first();
	}
	FastaParse(const string &fn, string *store, bool capit = false) :
		m_file(fn.c_str()),
		m_capitalize(capit)
	{
		m_seq = store;
		m_ext_storage = true;
		first();
	}
	
	~FastaParse();
	
	void use_external_storage(string *store) {
		if(m_ext_storage == false) {
			delete m_seq;
		} else {
    		m_ext_storage = true;
		}
		m_seq = store;
	}

	void first();
	bool next();

	bool is_valid() {
		return(m_is_valid);
	}
};

#endif // FastaParse 
