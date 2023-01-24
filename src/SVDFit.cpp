#include "port.h"
BASE_CC_FILE
#include "SVDFit.h"

//Default value for single precision and variables scaled to order unity. 

#define TOL 1.0e-5 

// some macros are commented in order to avoid compilation warnings.

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//Solves AX = B for a vector X, where A is speci ed by the arrays u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector. No input quantities are destroyed, so the routine may be called sequentially with di
void svbksb(vector<vector<double> > &u, vector<double> &w, 
		vector<vector<double> > &v, 
		int m, int n, vector<double> &b, 
		vector<double> &x)

{ 
	int jj,j,i;
	double s;
	vector<double> tmp(n+1);
	for (j=1;j<=n;j++) { //Calculate UTB. 
		s=0.0;
		if (w[j]) { //Nonzero result only if wj is nonzero. 
			//Rcpp::Rcerr << "bksb j " << j << " w " << w[j] << endl;
			for (i=1;i<=m;i++) {
				s += u[i][j]*b[i];
			}
			s /= w[j]; //This is the divide by wj . 
			//Rcpp::Rcerr << "bksb j " << j << " s " << s << endl;
		} 
		tmp[j]=s;
	} 
	for (j=1;j<=n;j++) { //Matrix multiply by V to get answer. 
		s=0.0;
		for (jj=1;jj<=n;jj++) {
			s += v[j][jj]*tmp[jj];
		}
		x[j]=s;
	} 
}

void svdcmp(vector<vector<double> > &a, int m, int n, vector<double> &w, 
		vector<vector<double> > &v)
//Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A = UWV T. Thematrix U replaces a on output. The diagonal matrix of singular values W is output as a vector w[1..n]. Thematrix V (not the transpose V T ) is output as v[1..n][1..n]. 
{ 
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;

//	ASSERT(a.size() == m + 1, "a size is " << a.size() << " but m " << m);
//	ASSERT(a[1].size() == n + 1, "a[1] size is " << a[0].size() << " but n " << n);
//	ASSERT(w.size() == n + 1, "w size is " << w.size() << " but n " << n);
//	ASSERT(v.size() == n + 1, "v size is " << v.size() << " but n " << n);
//	ASSERT(v[1].size() == n + 1, "v[1] size is " << v[1].size() << " but n " << n);

	vector<double> rv1(n+1);
	g=scale=anorm=0.0;
	//Householder reduction to bidiagonal form. 
	for (i=1;i<=n;i++) { 
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) { 
			for(k=i;k<=m;k++) {
				scale += fabs(a[k][i]);
			}
			if (scale) { 
				for (k=i;k<=m;k++) { 
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				} 
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) { 
					s = 0.0;
					for (k=i;k<=m;k++) {
						s += a[k][i]*a[k][j];
					}
					f=s/h;
					for (k=i;k<=m;k++) {
						a[k][j] += f*a[k][i];
					}
				} 
				for(k=i;k<=m;k++) {
					 a[k][i] *= scale;
				}
			} 
		} 
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) { 
			for (k=l;k<=n;k++) {
				scale += fabs(a[i][k]);
			}
			if(scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				} 
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) {
					rv1[k]=a[i][k]/h;
				}
				for(j=l;j<=m;j++) { 
					s=0.0;
					for(k=l;k<=n;k++) {
						s += a[j][k]*a[i][k];
					}
					for(k=l;k<=n;k++) {
						a[j][k] += s*rv1[k];
					}
				} 
				for(k=l;k<=n;k++) {
					a[i][k] *= scale;
				}
			} 
		} 
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	} 
	for(i=n;i>=1;i--) { //Accumulation of right-hand transformations. 
		if(i < n) { 
			if (g) { 
				for (j=l;j<=n;j++) { //Double division to avoid possible under ow. 
					v[j][i]=(a[i][j]/a[i][l])/g;
				}
				for(j=l;j<=n;j++) { 
					s=0.0;
					for (k=l;k<=n;k++) {
						s += a[i][k]*v[k][j];
					}
					for(k=l;k<=n;k++) { 
						v[k][j] += s*v[k][i];
					}
				} 
			} 
			for(j=l;j<=n;j++) {
				v[i][j]=v[j][i]=0.0;
			}
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	} 
	for(i=IMIN(m,n);i>=1;i--) { //Accumulation of left-hand transformations. 
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) {
			a[i][j]=0.0;
		}
		if(g) { 
			g=1.0/g;
			for(j=l;j<=n;j++) { 
				s=0.0;
				for(k=l;k<=m;k++) {
					s += a[k][i]*a[k][j];
				}
				f=(s/a[i][i])*g;
				for(k=i;k<=m;k++) {
					a[k][j] += f*a[k][i];
				}
			} 
			for(j=i;j<=m;j++) {
				a[j][i] *= g;
			}
		} else {
			for(j=i;j<=m;j++) {
				 a[j][i]=0.0;
			}
		}
		++a[i][i];
	}
	for(k=n;k>=1;k--) { //Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations. 
		for(its=1;its<=30;its++) { 
			flag=1;
			for(l=k;l>=1;l--) { //Test for splitting. 
				nm=l-1; //Note that rv1[1] is always zero. 
				if ((double)(fabs(rv1[l])+anorm) == anorm) { 
					flag=0;
					break;
				}
				if((double)(fabs(w[nm])+anorm) == anorm) {
					 break;
				}
			} 
			if(flag) { 
				c=0.0; //Cancellation of rv1[l], if l > 1. 
				s=1.0;
				for(i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if((double)(fabs(f)+anorm) == anorm) {
						 break;
					}
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for(j=1;j<=m;j++) { 
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					} 
				} 
			} 
			z=w[k];
			if(l == k) { //Convergence. 
				if(z < 0.0) {// Singular value is made nonnegative. 
					w[k] = -z;
					for(j=1;j<=n;j++) {
						 v[j][k] = -v[j][k];
					}
				} 
				break;
			} 
			if(its == 30) {
				 ASSERT(false, "no convergence in 30 svdcmp iterations");
			}
			x=w[l];
			//Shift from bottom 2-by-2 minor. 
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			//Next QR transformation: 
			for(j=l;j<=nm;j++) { 
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for(jj=1;jj<=n;jj++) { 
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				} 
				z=pythag(f,h);
				w[j]=z;
				//Rotation can be arbitrary if z = 0. 
				if(z) { 
					z=1.0/z;
					c=f*z;
					s=h*z;
				} 
				f=c*g+s*y;
				x=c*y-s*g;
				for(jj=1;jj<=m;jj++) { 
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				} 
			} 
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	} 
} 

//Computes (a2 + b2)1/2 without destructive under ow or over ow. 
double pythag(double a, double b) 
{ 
	double absa,absb; 
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) 
		return absa*sqrt(1.0+SQR(absb/absa));
	else 
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

//Code from NRC
//Inputs: x - the data points, 
//y - the value to fit to.
//ndata - number of points,
//a - the coefficients
//ma number of dimensions
//u - auxilary ndata on ma matrix
//v auxilary ma on ma matrix
//w - auxilary ma vector

void svdfit(vector<vector<double> > &x, vector<double> &y, vector<double> &sig, 
		int ndata, vector<double> &a, int ma, 
		vector<vector<double> > &u, 
		vector<vector<double> > &v, 
		vector<double> &w, 
		double *chisq)
{
	int j,i; 
	double wmax,tmp,thresh,sum;

	vector<double> b(ndata + 1);

	u = x;

	for(i=1;i<=ndata;i++) { //Accumulate coe cients of the  tting matrix. 
		tmp=1.0/sig[i]; 
		for(vector<double>::iterator j = u[i].begin() + 1; 
		    j != u[i].end();
		    j++) {
			*j *= tmp;
		}
		b[i]=y[i]*tmp; 
	}

	svdcmp(u, ndata, ma, w, v); //Singular value decomposition. 
/*dump matrices
	Rcpp::Rcerr << "V mat:" << endl;
	for(i = 1; i <= ma ; i++) {
		Rcpp::Rcerr << i;
		for(j = 1; j <= ma; j++) {
			Rcpp::Rcerr << "\t" << v[i][j];
		}
		Rcpp::Rcerr << endl;
	}
	Rcpp::Rcerr << "U mat:" << endl;
	for(i = 1; i <= ndata; i++) {
		Rcpp::Rcerr << i;
		for(j = 1; j <= ma; j++) {
			Rcpp::Rcerr << "\t" << u[i][j];
		}
		Rcpp::Rcerr << endl;
	}
*/
	wmax=0.0; //Edit the singular values, given TOL 
	for(j=1;j<=ma;j++) {
	//	Rcpp::Rcerr << "w val " << j << " was " << w[j] << endl;
		if (w[j] > wmax) {
			wmax=w[j]; 
		}
	}
//	Rcpp::Rcerr << "Max w val was " << wmax << endl;
	thresh=TOL*wmax; 
	for (j=1;j<=ma;j++) {
		if (w[j] < thresh) {
			Rcpp::Rcerr << "nullify singular val " << w[j] << " at " << j << endl;
			w[j]=0.0; 
		}
	}
	svbksb(u,w,v,ndata,ma,b,a); 
/*	for(j = 1; j <= ma; j++) {
		Rcpp::Rcerr << "a[" << j << "] = " << a[j] << endl;
	}
*/
	*chisq=0.0; //Evaluate chi-square. 
	for (i=1;i<=ndata;i++) { 
		for (sum=0.0,j=1;j<=ma;j++) {
			 sum += a[j]*x[i][j];
		} 
		*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp); 
	} 
}

//To evaluate the covariance matrix cvm[1..ma][1..ma] of the  t for ma parameters obtained by svdfit, call this routine with matrices v[1..ma][1..ma], w[1..ma] as returned from svdfit. 
void svdvar(vector<vector<double> > &v, int ma, vector<double> &w, vector<vector<double> > &cvm) 
{ 
	int k,j,i; 
	double sum;

	vector<double> wti(ma + 1);

	for (i=1;i<=ma;i++) { 
		wti[i]=0.0; 
		if (w[i]) 
			wti[i]=1.0/(w[i]*w[i]); 
	} 
	for (i=1;i<=ma;i++) { //Sum contributions to covariance matrix (15.4.20). 
		for (j=1;j<=i;j++) { 
			for (sum=0.0,k=1;k<=ma;k++) 
				sum += v[i][k]*v[j][k]*wti[k]; 
			cvm[j][i]=cvm[i][j]=sum; 
		} 
	} 
}

