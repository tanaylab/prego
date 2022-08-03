#ifndef seqpack_DnaPWML_h
#define seqpack_DnaPWML_h 1

/*=================================================
=================================================*/

#include <string>
#include <list>
#include "DnaPSSM.h"

class DnaPWML  {

protected:

	DnaPSSM m_pssm;

	vector<float> m_spat_factors;
	int m_spat_bin_size;

public:

	const DnaPSSM &get_pssm() const {
		return(m_pssm);
	}
	DnaPSSM &get_pssm() {
		return(m_pssm);
	}
	void set_pssm(const DnaPSSM &pssm) {
		m_pssm = pssm;
	}

	int get_max_pos() const {
		return(m_pssm.size());
	}
	int size() const {
		return(m_pssm.size());
	}

	float get_p(int pos, int nuc_i) {
		return(m_pssm[pos].get_direct_prob(nuc_i));
	}
	void  set_p(int pos, int nuc_i, float p) {
		return(m_pssm[pos].set_direct_prob(nuc_i, p));
	}

	int get_spat_bin_size() const {
		return(m_spat_bin_size);
	}
	int set_spat_bin_size(int sz) {
		m_spat_bin_size = sz;
	}
	int set_num_spat_bins(int nm) {
		m_spat_factors.resize(nm);
	}

	float get_spat_factor(int pos) const {
		return(m_spat_factors[int(pos/m_spat_bin_size)]);
	}

	float get_spat_bin_factor(int bin) const{
		return(m_spat_factors[bin]);
	}
	int get_num_spat_bins() const {
		return(m_spat_factors.size());
	}

	void set_spat_factor(int bin, float f) {
		if(bin >= m_spat_factors.size()) {
			m_spat_factors.resize(bin + 1, 0);
		}
		m_spat_factors[bin] = f;
	}
	const vector<float> &get_spat_dist() {
		return(m_spat_factors);
	}
	void set_spat_dist(const vector<float> &dst) {
		m_spat_factors = dst;
	}

	DnaPWML() {}

	DnaPWML(const DnaPSSM &pssm, vector<float> spat_fac, int bin_size) :
		m_pssm(pssm),
		m_spat_factors(spat_fac),
		m_spat_bin_size(bin_size)
	{}

	void init_from_seed(const string &seed, float prior) {
		m_pssm.init_from_seed(seed, prior);
		fill(m_spat_factors.begin(), m_spat_factors.end(), 1);
	}

	void integrate_energy(const string &target, float &energy) {
		m_pssm.integrate_energy(target, energy, m_spat_factors, m_spat_bin_size);
	}

	void randomize_pos(int pos);
		
	void write_tab(ostream &pssmd, ostream &spat, int id) const;
};

ostream &operator<<(ostream &out, const DnaPWML &pat);

#endif // DnaPWML
