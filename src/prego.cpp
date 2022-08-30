/*=================================================

  
=================================================*/

#include "port.h"
BASE_CC_FILE
#include "options.h"
#include "strutil.h"
#include "dnastrutil.h"
#include "Random.h"
#include "LeastSquare.h"
#include "SpecialFunc.h"
#include "BitVecIter.h"

#include "PssmRegression.h"
#include "PWMLRegression.h"
#include "KMerMultiStat.h"

const int OPT_DEFS=54;
const char *c_opt_defaults[OPT_DEFS] = {
"mode", "screen",
"OnlyScreen", "0",
"motif", "******",
"rseed", "19",
"xmin", "0",
"xmax", "-1",
"gmin", "0",
"gmax", "0",
"min_n", "50",
"min_cor", "0.08",
"L", "6",
"kmer_tab", "prego_kmers.txt",
"verbose", "0",
"pwml", "1",
"bidirect", "1",
"eps", "0.001",
"no_reg", "0",
"min_rms_for_star", "0.001",
"mod_id", "0",
"spat_min", "0",
"spat_max", "-1",
"min_nuc_prob", "0.001",
"spat_bin", "50",
"test_fn", "prego_test.txt",
"pssm_fn", "prego_pssm.txt",
"spat_fn", "prego_spat.txt",
"preds_fn", "prego_preds.txt"
};


void read_fasta(char *seq_fn, vector<string> &sequences)
{
	ifstream seq_file(seq_fn);
	string line;
	while(seq_file) {
		getline(seq_file, line, '\n');
		if(line.length() != 0) {
			sequences.push_back(line);
		}
	}
}

void read_responses(char *response_fn, vector<vector<float> > &response_stat, 
								vector<int> &is_train, int &n_in_train)
{
	ifstream stat(response_fn);
	n_in_train = 0;
	char *test, *test2;
	vector<string> fields;
	while(stat) {
		split_line(stat, fields, '\t');
		if(fields.size() == 0) {
			continue;
		}
		if((fields.size()-1) > response_stat.size()) {
			response_stat.resize(fields.size()-1);
		}
		int train = atoi(fields[0].c_str());
    	try {
			is_train.push_back(stoi(fields[0]));
    	}
    	catch (std::invalid_argument const &e) {
        std::cout << "Bad input: std::invalid_argument thrown" << std::endl;
    	}
    	catch (std::out_of_range const &e) {
        std::cout << "Integer overflow: std::out_of_range thrown" << std::endl;
    	}
		if(train) {
			n_in_train++;
		}

		for(int i = 1; i < fields.size(); i++) {
			float response = strtof(fields[i].c_str(), &test);
			if(test != fields[i].c_str()) {
				response_stat[i-1].push_back(response);
			}
		}
	}
}

void screen_kmers(vector<string> &sequences, vector<vector<float> > &response_stat,
						vector<int> &is_train, int n_in_train, 
						options &opt)
{
	int resp_dim = response_stat.size();
	int L = opt.get_int("", "L");
	int from_range = opt.get_int("", "xmin");
	int to_range = opt.get_int("", "xmax");
	if(to_range == -1) {
		to_range = sequences[0].length();
	}
	float corr_thresh = opt.get_float("", "min_cor");
	float r2_thresh = corr_thresh*corr_thresh;
	int min_n = opt.get_int("", "min_n");
	ofstream kmer_tab(opt.get_str("", "kmer_tab"));
	int bin_num = 20;
	int norm = 0;
	float norm_factor = 1;

	KMerMultiStat multi(L, 0,
			opt.get_int("", "gmin"),
			opt.get_int("", "gmax"),
			&sequences,
			&is_train,
			bin_num, norm, norm_factor,
			response_stat,
			from_range, to_range);

	string best_mot = "";
	float best_r2 = 0;
	vector<string> foc_mots;
	vector<float> foc_scores;
	vector<int> foc_ids;
	vector<float> response_avg(resp_dim, 0);
	vector<float> response_var(resp_dim, 0);

	for(int ri = 0; ri < resp_dim; ri++) {
		vector<int>::iterator train = is_train.begin();
		for(vector<float>::iterator i = response_stat[ri].begin(); 
			 i != response_stat[ri].end(); 
			 i++) {
			if(*train) {
				response_avg[ri] += *i;
				response_var[ri] += (*i)*(*i);
			}
			train++;
		}
		response_avg[ri] /= n_in_train;
		response_var[ri] /= n_in_train;
		response_var[ri] -= response_avg[ri]*response_avg[ri];
	}
	cerr << "done normalizing response " << endl;

	for(map<const string, vector<pair<int, vector<float> > > >::const_iterator k = multi.get_pat_begin();
   	 k != multi.get_pat_end();
    	k++) {
		vector<float> cov(resp_dim, 0);
		vector<float> corr(resp_dim, 0);
		float avg_multi = 0;
		float tot_multi2 = 0;
//ODO - compute cov over the rdim
		const vector<pair<int, vector<float> > > &multi = k->second;
		for(int m = 1; m < multi.size(); m++) {
			avg_multi += multi[m].first*m;
			tot_multi2 += multi[m].first*m*m;
			for(int ri = 0; ri < resp_dim; ri++) {
				cov[ri] += m * multi[m].second[ri];
			}
		}
		avg_multi /= n_in_train;
		float multi_var = tot_multi2/n_in_train - avg_multi*avg_multi;
		float max_r2 = 0;
		for(int ri = 0; ri < resp_dim; ri++) {
			cov[ri] /= n_in_train;
			cov[ri] -= avg_multi * response_avg[ri];
			corr[ri] = cov[ri]/sqrt(multi_var * response_var[ri]);
			if(max_r2  < corr[ri]*corr[ri]) {
				max_r2 = corr[ri]*corr[ri];
			}
		}
		if(max_r2 > r2_thresh) {
			kmer_tab << k->first << "\t" << max_r2 << "\t" << avg_multi << "\t" << multi_var;
			for(int ri = 0; ri < resp_dim; ri++) {
				kmer_tab << "\t" << corr[ri];
			}
			kmer_tab << "\n";

			foc_mots.push_back(k->first);
			foc_scores.push_back(max_r2);
			foc_ids.push_back(foc_ids.size());
			if(max_r2 > best_r2 && (avg_multi*n_in_train) > min_n) {
				best_r2 = max_r2;
				best_mot = k->first;
				cerr << "new best " << best_mot << "  " << best_r2 << endl;
			}
		}
	}
	//seedmot = "**"+best_mot+"**";
//		if(k->first == "TCCCAG" | k->first == "AAAAAA") {
	//		for(int m = 1; m < multi.size(); m++) {
	//			cerr << k->first << " " << multi[m].first << " " << multi[m].second[0] << " " << multi[m].second[0]/(multi[m].first) << endl;
	//		}
	//		cerr << "response " << response_avg[0] << " " << response_var[0] << endl;
	//		cerr << "avg multi " << avg_multi << " " << multi_var << endl;
	//		cerr << "cov " << cov[0] << " cor " << corr[0] << endl;
	//	}
}

int main(int argc, char *argv[])
{
	options opt;
	opt.load_defaults(c_opt_defaults, OPT_DEFS);
	opt.parse_argv(argc, argv);

	const string &mode = opt.get_str("", "mode", "screen");
	if(mode != "regress" && mode != "screen" && mode != "full") {
		cerr << "mode must be regress, screen or full" << endl;
		return(1);
	}
	if(argc < 3) {
		cerr << "usage prego fasta response options [options]" << endl;
		return(1);
	}
	char *seq_fn = argv[1];
	char *response_fn = argv[2];

	Random::reset(opt.get_int("", "rseed"));

//read sequences
	vector<string> sequences;
	read_fasta(seq_fn, sequences);

//read responses
	vector<vector<float> > response_stat;
	vector<int> is_train;
	int n_in_train;
	int resp_dim = response_stat.size();

	read_responses(response_fn, response_stat, is_train, n_in_train);


	cerr << "done reading " << sequences.size() << " sequences" << endl;
	cerr << "response dim " << resp_dim << endl;
	ASSERT(sequences.size() == response_stat[0].size(), 
				"incompatible sequence and response table size " 
							<< sequences.size() << " " << response_stat[0].size());

	cerr << "mode is " << mode << endl;
	if(mode == "screen" || mode == "full") {
		screen_kmers(sequences, response_stat, is_train, n_in_train, opt);
	}

	if(mode == "screen") {
		return(0);
	}
	const string &motif = opt.get_str("", "motif");
	string seedmot(motif);

	if(seedmot == "") {
		cerr << "no motif found by screen, reverting to stars" << endl;
		seedmot = "******";
	} else {
		cerr << "using seed " << seedmot << endl;
	}
	int isbid = opt.get_int("", "bidirect");

	float epsilon = opt.get_float("", "eps");
	cerr << "epislon is " << epsilon << endl;
	int no_reg = opt.get_float("", "no_reg");
	float min_rms_for_star = opt.get_float("", "min_rms_for_star");

	int psid = opt.get_int("", "mod_id");
	ofstream test_out(opt.get_str("", "test_fn").c_str());
	ofstream pssm_out(opt.get_str("", "pssm_fn").c_str());
	const string &spat_fn = opt.get_str("", "spat_fn");
	const string &preds_fn = opt.get_str("", "preds_fn");
	ofstream spat_out(spat_fn.c_str());
	ASSERT(spat_out, "Cannot open spat dist file for writing");

	int with_spat = opt.get_int("", "pwml");
	vector<float> res(4); 
	res[0] = 0.05; res[1]=0.02; res[2] = 0.01; res[3]= 0.005; //res[4]= 0.005;
	vector<float> spres(4); 
	spres[0] = 0.01; spres[1]=0.01; spres[2] = 0.01; spres[3]= 0.005; //spres[4]= 0.005;

	int smin = opt.get_int("", "spat_min");
	int smax = opt.get_int("", "spat_max");
	if(smax == -1) {
		smax = sequences[0].length();
	}
	if(with_spat) {
		cerr << "into pwmlreg" << endl;
		PWMLRegression pwmlreg(
			sequences,
			is_train,	
			smin, smax,
			opt.get_float("", "min_nuc_prob"),
			opt.get_int("", "spat_bin"),
			res, spres,
			0.001, 0.001, 0.01);

		pwmlreg.add_responses(response_stat);

		pwmlreg.m_logit = opt.get_int("", "verbose");
		pwmlreg.init_seed(seedmot, isbid);
		cerr << "done init seed " << seedmot << endl;
		pwmlreg.optimize();
		pwmlreg.output_pssm(pssm_out, spat_out, psid);
		pssm_out.close();
		DnaPWML pwml;
		pwmlreg.get_model(pwml);
		ofstream preds_tab(preds_fn);

		for(int i = 0; i < sequences.size(); i++) {
			float energy;
			pwml.integrate_energy(sequences[i], energy);
			preds_tab << i << "\t" << energy 
						<< "\t" << is_train[i] << endl;
		}
	}
}

/*
 
		int tst_n = 0;
		int tr_n = 0;
		float tr_xy = 0;
		float tr_x = 0;
		float tr_x2 = 0;
		float tr_y = 0;
		float tr_y2 = 0;
		float tst_xy = 0;
		float tst_x = 0;
		float tst_x2 = 0;
		float tst_y = 0;
		float tst_y2 = 0;

		for(int i = 0; i < sequences.size(); i++) {
			float energy;
			pwml.integrate_energy(sequences[i], energy);
			if(is_train[i]) {
				tr_xy += energy*response_stat[i];
				tr_x += energy;
				tr_x2 += energy*energy;
				tr_y += response_stat[i];
				tr_y2 += response_stat[i]*response_stat[i];
				tr_n++;
			} else {
				tst_xy += energy*response_stat[i];
				tst_x += energy;
				tst_x2 += energy*energy;
				tst_y += response_stat[i];
				tst_y2 += response_stat[i]*response_stat[i];
				tst_n++;
			}
			preds_tab << i << "\t" << energy 
						<< "\t" << response_stat[i] 
						<< "\t" << is_train[i] << endl;
		}
		float cur_corr = pwmlreg.get_cur_corr();
		tr_xy = tr_xy/tr_n;
		tr_x = tr_x/tr_n;
		tr_y = tr_y/tr_n;
		tr_x2 = tr_x2/tr_n;
		tr_y2 = tr_y2/tr_n;
		tst_xy = tst_xy/tst_n;
		tst_x = tst_x/tst_n;
		tst_y = tst_y/tst_n;
		tst_x2 = tst_x2/tst_n;
		tst_y2 = tst_y2/tst_n;
		float tr_cor =  (tr_xy- tr_x*tr_y) / 
								sqrt((tr_x2-tr_x*tr_x)*(tr_y2-tr_y*tr_y));
		float tst_cor =  (tst_xy- tst_x*tst_y) / 
								sqrt((tst_x2-tst_x*tst_x)*(tst_y2-tst_y*tst_y));
		cerr << "pwml tr cor " << cur_corr << " recomp " << tr_cor << " tst cor " << tst_cor << endl;
		preds_tab.close();
*/
