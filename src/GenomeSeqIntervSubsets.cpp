#include "port.h"
BASE_CC_FILE
#include "GenomeSeqIntervSubsets.h"
#include "strutil.h"

void GenomeSeqIntervSubsets::read(istream &tab, int generic_format)
{
	string chrom;
	int coord,tocoord;
	string setname;
	int fr, to;

	int max_interv_id = m_intervs.max_interv_id();
	int missed = 0;
	int lcount = 0;
	vector<string> fields;
	tab >> chrom;
	while(tab) {
		if(generic_format) {
			split_line(tab, fields, '\t');
			fr = 0;
			to = 0;
			ASSERT(fields.size() > 3, "Bad interv line " << lcount << " when reading sets\n");
			coord = atoi(fields[1].c_str());
			tocoord = atoi(fields[2].c_str());
			setname = fields[3];
		} else {
    		tab >> coord >> tocoord >> setname >> fr >> to;
		}
		lcount++;
		int setid;
		if(m_setname_to_id.find(setname) == m_setname_to_id.end()) {
			m_setname_to_id[setname] = m_subsets.size();
			m_subset_name.push_back(setname);
			setid = m_subsets.size();
			m_subsets.resize(setid + 1, ds_bitvec(max_interv_id));
			m_subsets[setid].set(false);
			m_offsets.resize(setid + 1, vector<pair<int, int> >(max_interv_id, pair<int, int>(0, 0)));
		} else {
			setid = m_setname_to_id[setname];
		}
		int id = m_intervs.interv_id(chrom, coord);
		if(id != -1) {
			m_subsets[setid][id] = 1;
			m_offsets[setid][id] = pair<int, int>(fr,to);
		} else {
			missed++;
		}
		tab >> chrom;
	}
	cerr << "lcount " << lcount << " missed " << missed << " probes" << endl;
}
