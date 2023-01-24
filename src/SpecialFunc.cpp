#include <cmath>
#include "port.h"
BASE_CC_FILE

#include "SpecialFunc.h"

float gamma_ln(float xx)
{
	double x, y, tmp, ser;

	static double cof[6] = {
		76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5};
	
	y = x = xx;

	tmp = x+5.5;
	tmp -= (x+0.5) * log(tmp);
	ser = 1.000000000190015;
	for(int j = 0; j <= 5; j++) {
		ser += cof[j]/++y;
	}
	return(-tmp+log(2.5066282746310005 * ser/x));
}

double dbl_gamma_ln(float xx)
{
	double x, y, tmp, ser;

	static double cof[6] = {
		76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5};
	
	y = x = xx;

	tmp = x+5.5;
	tmp -= (x+0.5) * log(tmp);
	ser = 1.000000000190015;
	for(int j = 0; j <= 5; j++) {
		ser += cof[j]/++y;
	}
	return(-tmp+log(2.5066282746310005 * ser/x));
}

FactorialFuncTable::FactorialFuncTable(int max_n) :
	m_table(max_n)
{


}

float FactorialFuncTable::factorial(int n) { return 0; }

float FactorialFuncTable::factorial_ln(int n) { return 0; }

BinomFuncTable::BinomFuncTable(int max_n, int max_k) :
	m_factor_table(max_n, -1),
	m_bin_table(max_n * max_k, -1),
	m_lnbin_table(max_n * max_k, -1),
	m_dbl_lnbin_table(max_n * max_k, -1),
	m_max_cache_n(max_n),
	m_max_cache_k(max_k)
{

	int i = 2;
	m_factor_table[0] = 1;
	m_factor_table[1] = 1;
	m_factor_table[2] = 2;

	while(i < max_n-1) {
		int k = i++;
		m_factor_table[i] = m_factor_table[k] * i;
	}
}

double BinomFuncTable::dbl_log_hg(int n, int n_i, int m, int k) 
{
	if(float(m)*n_i/float(n) > k) {
		return(0.0);
	}
	return(dbl_binom_ln(m, k) + dbl_binom_ln(n-m, n_i-k) - 
						dbl_binom_ln(n, n_i));
}

float BinomFuncTable::log_hg(int n, int n_i, int m, int k) 
{
	if(float(m)*n_i/float(n) > k) {
		return(0.0);
	}
	return(binom_ln(m, k) + binom_ln(n-m, n_i-k) - binom_ln(n, n_i));
}
float BinomFuncTable::log_neg_hg(int n, int n_i, int m, int k) 
{
	if(float(m)*n_i/float(n) < k) {
		return(0.0);
	}
	return(binom_ln(m, k) + binom_ln(n-m, n_i-k) - binom_ln(n, n_i));
}

float BinomFuncTable::binom_ln(int n, int k)
{
	if(n >= m_max_cache_n || k >= m_max_cache_k) {
		return((n == k || k==0) ? 0 : gamma_ln(n+1) - gamma_ln(k+1) - gamma_ln(n-k+1));
	}
	int ndx = n + k * m_max_cache_n;
	float &bico = m_lnbin_table[ndx];
	if(bico == -1) {
		bico = ((n == k || k==0) ? 0 : gamma_ln(n+1) - gamma_ln(k+1) - gamma_ln(n-k+1));
	}
	return(bico);
}

double BinomFuncTable::dbl_binom_ln(int n, int k)
{
	if(n >= m_max_cache_n || k >= m_max_cache_k) {
		return((n==k || k==0) ? 0 : dbl_gamma_ln(n+1) - dbl_gamma_ln(k+1) - dbl_gamma_ln(n-k+1));
	}
	int ndx = n + k * m_max_cache_n;
	double &bico = m_dbl_lnbin_table[ndx];
	if(bico == -1) {
		bico = ((n == k || k==0) ? 0 : dbl_gamma_ln(n+1) - 
					dbl_gamma_ln(k+1) - dbl_gamma_ln(n-k+1));
	}
	return(bico);
}

float BinomFuncTable::binom(int n, int k)
{
	int ndx = n + k * m_max_cache_n;
	if((size_t)ndx >= m_bin_table.size()) {
		return(exp(gamma_ln(n+1) - gamma_ln(k+1) - gamma_ln(n-k+1)));
	}
	float &bico = m_bin_table[n + k * m_max_cache_n];
	if(bico == -1) {
		m_bin_table[ndx] = m_factor_table[n]/
				(m_factor_table[k] * m_factor_table[n - k]);
	}
	return(bico);
}


//Returns the incomplete beta function I x (a; b).

double betai(double a, double b, double x)
{
	double betacf(double a, double b, double x);
	double dbl_gamma_ln(float xx);
	void nrerror(char error_text[]);
	double bt;
	if(x < 0.0 || x > 1.0) {
		Rcpp::Rcerr << "Bad x " << x<< " in routine betai";
		return(-1);
	}
	if(x == 0.0 || x == 1.0) {
		bt=0.0;
	} else {	//Factors in front of the continued fraction.
		bt=exp(dbl_gamma_ln(a+b)-dbl_gamma_ln(a)-dbl_gamma_ln(b)+a*log(x)+b*log(1.0-x));
	}
	if(x < (a+1.0)/(a+b+2.0)) { //Use continued fraction directly.
		return bt*betacf(a,b,x)/a;
	} else {	//Use continued fraction after making the 
			//symmetry transformation. 
		return 1.0-bt*betacf(b,a,1.0-x)/b;
	}
}

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

//Used by betai: Evaluates continued fraction for incomplete beta function by 
//modified Lentz's method ( x 5.2).

double betacf(double a, double b, double x)
{
	void nrerror(char error_text[]);

	int m = 0;
	int m2 = 0;
	double aa,c,d,del,h,qab,qam,qap;
	qab=a+b; 		//These q's will be used in factors that occur
				//in the coecients (6.4.6). 
	qap=a+1.0;
	qam=a-1.0;
	c=1.0; 			//First step of Lentz's method.
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) 
		d=FPMIN;

	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d; //One step (the even one) of the recurrence.
		if (fabs(d) < FPMIN) 
			d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) 
			c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d; //Next step of the recurrence (the odd one).
		if (fabs(d) < FPMIN) 
			d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) 
			c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS)
			break; 	//Are we done?
	}
	if (m > MAXIT)  {
		Rcpp::Rcerr << "a " << a << " or b " << b << " too big, or MAXIT too small in betacf, x = " << x << endl;
	}
	return h;
}

double cummulative_binom(float p, int n, int k) {
	double pv = betai(k, n-k+1, p);
	return(std::isnan(pv) ? 1 : pv);
}

//Returns the incomplete gamma function P(a; x). 
float gammp(float a, float x) 
{ 
	float gamser,gammcf,gln; 
	
	if (x < 0.0 || a <= 0.0) 
		Rcpp::Rcerr << "Invalid arguments in routine gammp" << endl;; 
	if (x < (a+1.0)) { 
		//Use the series representation. 
		gser(&gamser,a,x,&gln); 
		return gamser; 
	} else { 
		//Use the continued fraction representation 
		gcf(&gammcf,a,x,&gln); 
		return 1.0-gammcf; //and take its complement. 
	} 
} 

float gammq(float a, float x) // Returns the incomplete gamma function Q(a; x)   1   P(a; x). 
{ 
	float gamser,gammcf,gln; 

	if (x < 0.0 || a <= 0.0) 
		Rcpp::Rcerr << "Invalid arguments in routine gammq" << endl; 
	if (x < (a+1.0)) { 
		//Use the series representation 
		gser(&gamser,a,x,&gln); 
		return 1.0-gamser; //and take its complement. 
	} else { 
		//Use the continued fraction representation. 
		gcf(&gammcf,a,x,&gln); 
		return gammcf; 
	}
}

#define ITMAX 100     //Maximum allowed number of iterations. 
#define EPS 3.0e-7    //Relative accuracy. 
#define FPMIN 1.0e-30 //Number near the smallest representable floating-point number.

//Returns the incomplete gamma function P(a; x) evaluated by its series 
//representation as gamser. Also returns ln 
void gser(float *gamser, float a, float x, float *gln) 
{
	int n; 
	float sum,del,ap; 
	*gln=gamma_ln(a); 
	if (x <= 0.0) { 
		if (x < 0.0) 
			Rcpp::Rcerr << "x less than 0 in routine gser" << endl;; 
		*gamser=0.0; 
		return; 
	} else { 
		ap=a; 
		del=sum=1.0/a; 
		for (n=1;n<=ITMAX;n++) { 
			++ap; 
			del *= x/ap; 
			sum += del; 
			if (fabs(del) < fabs(sum)*EPS) { 
				*gamser=sum*exp(-x+a*log(x)-(*gln)); 
				return; 
			} 
		} 
		Rcpp::Rcerr << "a too large, ITMAX too small in routine gser" << endl; 
		return; 
	} 
}

//Returns the incomplete gamma function Q(a; x) evaluated by its continued fraction representation as gammcf. Also returns ln 
void gcf(float *gammcf, float a, float x, float *gln) 
{
	int i; 
	float an,b,c,d,del,h; 
	
	*gln=gamma_ln(a); 
	b=x+1.0-a; //Set up for evaluating continued fraction by modi ed Lentz's method (x5.2) with b0 = 0. 
	c=1.0/FPMIN; 
	d=1.0/b; 
	h=d; 
	for (i=1;i<=ITMAX;i++) { //Iterate to convergence. 
		an = -i*(i-a); 
		b += 2.0; d=an*d+b; 
		if (fabs(d) < FPMIN) 
			d=FPMIN; 
		c=b+an/c; 
		if (fabs(c) < FPMIN) 
			c=FPMIN; 
		d=1.0/d; 
		del=d*c; 
		h *= del; 
		if (fabs(del-1.0) < EPS) 
			break; 
	} 
	if (i > ITMAX) 
		Rcpp::Rcerr << "a too large, ITMAX too small in gcf" << endl;; 
	*gammcf=exp(-x+a*log(x)-(*gln))*h; //Put factors in front. 
}

float chip(float chi2, float neu) 
{
	return(gammp(neu/2, chi2/2));
}
float chiq(float chi2, float neu) 
{
	return(gammq(neu/2, chi2/2));
}
