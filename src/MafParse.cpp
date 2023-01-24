#include "port.h"

#include "MafParse.h"
#include "PhyloTree.h"

bool MafParse::next()
{
	if(!m_file) { 
		return(false); 
	}
	string slack;
	char code;
	m_file >> code;
	while(m_file && code == '#') {
		getline(m_file, slack, '\n');//eat comment
		m_file >> code;
	}
	if(!m_file) {
		return(false);
	}
	getline(m_file, slack, '\n');//eat a line
	

	bool not_blank = 1;
	int si = 0;
	string src;
	int size;
	while(not_blank) {
		code = m_file.get();
		if(code == '\n') {
			not_blank = 0;
			continue;
		}
		if(si >= m_desc.size()) {
			m_pos.resize(si + 1);
			m_max_pos.resize(si + 1);
			m_desc.resize(si + 1);
			m_seq.resize(si + 1);
			m_src_size.resize(si + 1);
			m_strand.resize(si + 1);
		}
		if(code != 's') {
			getline(m_file, slack, '\n');//eat a line
			
		} else {
    		m_file >> m_desc[si] >> m_pos[si] >> size >> m_strand[si] >> m_src_size[si] >> m_seq[si];
    		m_max_pos[si] = m_pos[si] + size;
    		code = m_file.get(); //eat trailing newline
    		si++;
		}
	}
	m_seq.resize(si);
	m_desc.resize(si);
	return(true);
}

int MafParse::get_species_id(const PhyloTree &phylo, int i)
{
	//which desc?
	string spname;
	string::const_iterator dsc = m_desc[i].begin();
	string::const_iterator max_dsc = m_desc[i].end();
	while(dsc < max_dsc && *dsc != '.') {
		spname.push_back(*dsc);
		dsc++;
	}
	return(phylo.get_sp_by_name(spname));
}
