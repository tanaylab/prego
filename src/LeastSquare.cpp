#include "port.h"

#include "LeastSquare.h"
#include "SVDFit.h"

void least_square(vector<vector<double> > &x, 
	    vector<double> &y, 
	    vector<double> &a, int dim, double &chisq)
{
	vector<double> sig(y.size(), 1);

	vector<vector<double> > u; 
	vector<vector<double> > v(dim + 2, vector<double>(dim + 2)); 
	vector<double> w(dim + 2);

	svdfit(x, y, sig, x.size() - 1, a, dim, u, v, w, &chisq);

	vector<vector<double> > covar(dim + 1, vector<double>(dim + 1));
	svdvar(v, dim, w, covar);
}

void least_square(vector<vector<float> > &x, 
	    vector<float> &y, 
	    vector<float> &a, int dim, float &chisq)
{
	vector<vector<double> > dx(x.size());
	vector<vector<float> >::iterator ix = x.begin();
	for(vector<vector<double> >::iterator idx = dx.begin(); idx != dx.end(); idx++) {
		idx->resize(ix->size());
		vector<float>::iterator ixx = (*ix).begin();
		for(vector<double>::iterator idxx = (*idx).begin(); 
		    idxx != (*idx).end(); 
		    idxx++) {
			*idxx = double(*ixx);
			ixx++;
		}
		ix++;
	}
	vector<double> dy(y.size());
	vector<float>::iterator iy = y.begin();
	for(vector<double>::iterator idy = dy.begin(); idy != dy.end(); idy++) {
		*idy = double(*iy);
		iy++;
	}
	vector<double> da(a.size());
	vector<float>::iterator ia = a.begin();
	for(vector<double>::iterator ida = da.begin(); ida != da.end(); ida++) {
		*ida = double(*ia);
		ia++;
	}
	double dchi;
	least_square(dx, dy, da, dim, dchi);
	ia = a.begin();
	for(vector<double>::iterator ida = da.begin(); ida != da.end(); ida++) {
		*ia = float(*ida);
		ia++;
	}
	chisq = dchi;
}
