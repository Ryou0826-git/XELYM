// spline.hpp
#ifndef SPLINE_HPP
#define SPLINE_HPP

#include <vector>

using namespace std;

struct s_coeff {
    vector<double> a;
    vector<double> b;
    vector<double> c;
    vector<double> d;
};

double Func3rd(const double& a, const double& b, const double& c, const double& d, 
               const double& xi, const double& x);

s_coeff Spline3rd(const vector<double>& x, 
                  const vector<double>& y);

vector<double> PerformSpline3rd(const s_coeff& coeff, 
                                const vector<double>& x, 
                                const vector<double>& xnew);

vector<vector<double>> 
PerformSpline3rd2d(const vector<double>& x, 
                   const vector<double>& y, 
                   const vector<vector<double>>& z, 
                   const vector<double>& xnew, 
                   const vector<double>& ynew);

#endif
