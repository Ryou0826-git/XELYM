// extract.hpp
#ifndef EXTRACT_HPP
#define EXTRACT_HPP

#include <vector>

using namespace std;

vector<vector<double>> 
ExtractCoord(const vector<vector<double>>& coord, 
             const vector<vector<int>>& molinfo);

#endif
