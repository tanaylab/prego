#ifndef FUNCTIONINTERPOLATOR_H
#define FUNCTIONINTERPOLATOR_H

#include <Rcpp.h>
#include <vector>

class FunctionInterpolator {
  private:    
    float xmin, xmax;
    int npts = 0;
    std::vector<float> xvals;
    std::vector<float> yvals;

  public:
    FunctionInterpolator() {}
    FunctionInterpolator(Rcpp::Function func_, float xmin_, float xmax_, int npts_);

    void init(Rcpp::Function func_, float xmin_, float xmax_, int npts_);

    float operator()(float x) const;

    std::vector<float> interpolate(const std::vector<float> &xvals) const;

    float get_xmin() const { return xmin; }
    float get_xmax() const { return xmax; }
    int get_npts() const { return npts; }

    std::vector<float> get_xvals() const {         
        std::vector<float> xvals_copy(xvals);
        return xvals_copy;
    }
    std::vector<float> get_yvals() const {
        std::vector<float> yvals_copy(yvals);
        return yvals_copy;
    }

    bool is_initialized() const { return npts > 0; }
};

#endif
