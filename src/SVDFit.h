#ifndef stdalg_alg_SVDFit_h
#define stdalg_alg_SVDFit_h 1

#include <vector>

void svbksb(vector<vector<double> > &u, vector<double> &w, 
		vector<vector<double> > &v, 
		int m, int n, vector<double> &b, 
		vector<double> &x);
void svdcmp(vector<vector<double> > &a, int m, int n, vector<double> &w, 
		vector<vector<double> > &v);
double pythag(double a, double b);
void svdfit(vector<vector<double> > &x, vector<double> &y, vector<double> &sig, 
		int ndata, vector<double> &a, int ma, 
		vector<vector<double> > &u, 
		vector<vector<double> > &v, 
		vector<double> &w, 
		double *chisq);
void svdvar(vector<vector<double> > &v, int ma, vector<double> &w, vector<vector<double> > &cvm);
#endif 
