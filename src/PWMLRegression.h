#ifndef seqpack_PWMLRegression_h
#define seqpack_PWMLRegression_h 1

/*=================================================
=================================================*/

#include <vector>
#include "DnaPWML.h"
#include "BitVec.h"

class PWMLRegression  {

private:
	struct UpdAux {
		vector<float>::iterator p;
		float factor;
	};

	struct NeighStep {
		char nuc;
		float diff;
		NeighStep(char n, float d) :
			nuc(n),
			diff(d)
		{}
	};

protected:

//Relevant chromosomes and probes (if in tiling mode)
	const vector<string> &m_sequences;
	const vector<int> &m_train_mask;
	int m_train_n;

//response dimension
	int m_rdim;
//epigenomic readouts to fit to (blocks of m_rdim in the vector)
	vector<float> m_interv_stat;
	vector<double> m_data_avg;
	vector<double> m_data_var;

	int m_min_range;
	int m_max_range;
	
	vector<bool> m_is_wildcard;
	
//m_bidirect[i] is mixture pssm i bidirectional?
	bool m_bidirect;

//m_chars[mix_id][pos]['A'] = weight for A in position pos
	vector<vector<float> > m_nuc_factors;
	vector<vector<vector<float> > > m_derivs; //deriv[vid][pos][char]

	float m_min_prob;

//step sizes at each phase
	vector<float> m_resolutions;
	vector<float> m_spat_resolutions;

	vector<float> m_spat_factors;
	vector<vector<float> > m_spat_derivs; //deriv[vid][pos][char]
	int m_spat_bin_size;
	float m_spat_factor_step;

	float m_cur_r2;

//cur best fit line m_a + m_b*pred(seq)
	vector<float> m_a;
	vector<float> m_b;

	float m_min_rms_for_star;

	float m_unif_prior;

	vector<UpdAux> m_aux_upds;

//defining the local search policty
	vector<vector<NeighStep> > m_cur_neigh;
public:

	bool m_logit;

	float get_cur_r2() const {
		return(m_cur_r2);
	}

//init in tiling mode
	PWMLRegression(
		const vector<string> &loci,
		const vector<int> &train_mask,
		int min_range,
		int max_range,
		float min_prob,
		int spat_bin_size,
		const vector<float> &resolutions,
		const vector<float> &s_resolutions,
		float eps,
		float min_improv_for_star,
		float unif_prior);

	void add_responses(const vector<vector<float> > &stats);

	void init_seed(const string &init_mot, int isbid);
	void init_pwm(DnaPSSM &pwm);

	void optimize();

	void output_pssm(ostream &out, ostream &spat, int psid);
	void report_cur_lpwm(ostream &out);

	void init_energies();
	void fill_predictions(vector<float> &pred);
	void get_diff_vector(vector<float> &to_diff, float &avg, float &var);
	float get_r2();
	void get_model(DnaPWML &pwml);
private:
	void init_neighborhood(float resolution);

	void update_seq_interval(int prb_id, string::const_iterator min_i, 
			string::const_iterator max_i, int sign, int pos);

	void take_best_step();

	float compute_cur_r2(int pos, vector<float> &probs);
	float compute_cur_r2_spat();
};

#endif //PWMLRegression_h
