#include "port.h"
BASE_CC_FILE
#include "DnaPWML.h"
#include "Random.h"


void DnaPWML::write_tab(ostream &pssmd, ostream &spat, int id) const
{
	m_pssm.write_tab(pssmd, id);

	for(size_t i = 0; i < m_spat_factors.size(); i++) {
		spat << id << "\t" << i * m_spat_bin_size << "\t" << m_spat_factors[i] << endl;
	}
}

void DnaPWML::randomize_pos(int pos) {
	vector<float> p(4);
	p[0] = Random::fraction();
	p[1] = Random::fraction();
	p[2] = Random::fraction();
	p[3] = Random::fraction();
	float tot = p[0] + p[1] + p[2] + p[3];
	p[0] /= tot;
	p[1] /= tot;
	p[2] /= tot;
	p[3] /= tot;

	m_pssm[pos].reset(p);
}

ostream &operator<<(ostream &out, const DnaPWML &pwml)
{
	out << pwml.get_pssm() << endl;
	out << "Spat: ";
	
	for(int i = 0; i < pwml.get_num_spat_bins(); i++) {
		out << "  " << pwml.get_spat_bin_factor(i);
	}
	out << endl;
	return(out);
}

