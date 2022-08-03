#ifndef stdalg_funcs_SpecialFunc_h
#define stdalg_funcs_SpecialFunc_h 1

//NRC implementation for usefull combinatorial functions. See NRC chapter 6
//for full details. Note that we are wrapping the factorial related functions
//in classes to provide storage for lookup table without using static
//initialization or singletons. The user of these function can always introduce
//a singleton himself if he likes to.
//

#include <vector>

//Using Lanczos power series decomposition (based on gamma complex completions)

float gamma_ln(float xx);
double dbl_gamma_ln(float xx);
double betai(double a, double b, double x);
double betacf(double a, double b, double x);
double cummulative_binom(float p, int n, int k);

float gammp(float a, float x);
float gammq(float a, float x);
void gcf(float *gammcf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);
float chip(float chi2, float neu);
float chiq(float chi2, float neu);

class FactorialFuncTable {

private:
	std::vector<float> m_table;

public:

	FactorialFuncTable(int nmax = 32);

	//Factorial (using a lookup table for small readings)
	float factorial(int n);

	//Factorial (using a lookup table for small readings)
	float factorial_ln(int n);
};

class BinomFuncTable {

private:

	std::vector<float> m_factor_table;
	std::vector<float> m_bin_table;
	std::vector<float> m_lnbin_table;
	std::vector<double> m_dbl_lnbin_table;

	int m_max_cache_n;
	int m_max_cache_k;

public:

//Using gamma
	BinomFuncTable(int max_n = 32, int max_k = 32);

	float binom_ln(int n, int k);
	double dbl_binom_ln(int n, int k);
	float log_hg(int n, int n_i, int m, int k);
	float log_neg_hg(int n, int n_i, int m, int k);
	double dbl_log_hg(int n, int n_i, int m, int k);
	float binom(int n, int k);
};

#endif // stdalg_funcs_SpecialFunc_h
