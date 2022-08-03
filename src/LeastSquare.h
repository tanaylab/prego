#ifndef stdalg_alg_LeastSquare_h
#define stdalg_alg_LeastSquare_h 1

#include <vector>

void least_square(vector<vector<double> > &x, 
	    vector<double> &y, 
	    vector<double> &a, int dim, double &chisq);

void least_square(vector<vector<float> > &x, 
	    vector<float> &y, 
	    vector<float> &a, int dim, float &chisq);
#endif 

