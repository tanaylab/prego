#ifndef seqpack_PssmRegression_h
#define seqpack_PssmRegression_h 1

/*=================================================
=================================================*/

#include <vector>
#include "BitVec.h"

class PssmRegression  {

protected:

//m_chars[pos]['A'] = weight for A in position pos
	vector<vector<float> > m_chars;
	int m_bidirect;

	vector<vector<float> > m_coefs;

//Data to fit to
	vector<float> m_interv_data;

//Relevant vars
	const vector<string> &m_sequences;
	
	const ds_bitvec &m_relevant_vars;

	float m_cur_rms;

	vector<float> m_pos_coefs;

	float m_base_coef;
	float m_epsilon;

	int m_no_neg_mode;

	vector<bool> m_is_wildcard;

	float m_min_rms_for_star;

public:

	PssmRegression(
		const vector<string> &loci,
		const string &init_mot,
		int isbid,
		const ds_bitvec &relev,
		vector<float> &chip_stat,
		float eps,
		int no_reg_mode, float min_rms_for_star);
	
	float go();
	void compute_wgts(int pos);
	void reopti_pos(int pos);

	float test(const ds_bitvec &test_set, 
				const vector<float> &chip_data, ostream &out);
	float output_pssm(ostream &out, int psid);

	void update_yvals(const ds_bitvec &relevant_set, 
					vector<float> &chip_data);
private:
	void update_coefs(vector<float> &coefs, int vid, int pos);
};

#endif //PssmRegression_h
