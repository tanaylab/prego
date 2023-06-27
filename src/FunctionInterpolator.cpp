#include "FunctionInterpolator.h"

FunctionInterpolator::FunctionInterpolator(Rcpp::Function func_, float xmin_, float xmax_,
                                           int npts_) {
    init(func_, xmin_, xmax_, npts_);
}

void FunctionInterpolator::init(Rcpp::Function func_, float xmin_, float xmax_, int npts_) {    
    xmin = xmin_;
    xmax = xmax_;
    npts = npts_;
    float step = (xmax - xmin) / (npts - 1);
    xvals.resize(npts);
    yvals.resize(npts);
    for (int i = 0; i < npts; i++) {
        xvals[i] = xmin + i * step;
        yvals[i] = Rcpp::as<float>(func_(xvals[i]));
    }
}

float FunctionInterpolator::operator()(float x) const {
    if (x < xmin){
        Rcpp::warning("x is out of range. Please set xmin and xmax (the domain of the energy function) correctly. The current value of x is: " + std::to_string(x) + ". while the domain is [" + std::to_string(xmin) + ", " + std::to_string(xmax) + "].");
        return yvals[0];
    }
    if (x > xmax){
        Rcpp::warning("x is out of range. Please set xmin and xmax (the domain of the energy function) correctly. The current value of x is: " + std::to_string(x) + ". while the domain is [" + std::to_string(xmin) + ", " + std::to_string(xmax) + "].");
        return yvals[npts - 1];
    }
    
    int i = std::lower_bound(xvals.begin(), xvals.end(), x) - xvals.begin();
    if (i == 0) {
        return yvals[0];
    } else if (i == npts) {
        return yvals[npts - 1];
    } else {
        float x0 = xvals[i - 1];
        float x1 = xvals[i];
        float y0 = yvals[i - 1];
        float y1 = yvals[i];
        return y0 + (y1 - y0) / (x1 - x0) * (x - x0);
    }
}

std::vector<float> FunctionInterpolator::interpolate(const std::vector<float> &xvals) const {
    std::vector<float> y_interp(xvals.size());    
    for (size_t i = 0; i < xvals.size(); i++) {
        y_interp[i] = operator()(xvals[i]);        
    }
    return y_interp;
}
