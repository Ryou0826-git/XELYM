// akima.hpp
#ifndef AKIMA_HPP
#define AKIMA_HPP

#include <vector>

using namespace std;

vector<double> PerformAkima1d(const vector<double>& x, 
                              const vector<double>& y, 
                              const vector<double>& xnew);

vector<vector<double>> 
PerformAkima2d(const vector<double>&         x, 
               const vector<double>&         y, 
               const vector<vector<double>>& z, 
               const vector<double>&         xnew, 
               const vector<double>&         ynew);

#endif
