// xtc_tool.hpp

#ifndef XTC_TOOLS_HPP
#define XTC_TOOLS_HPP

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <stdexcept>
#include <xdrfile_xtc.h>

using namespace std;

struct s_xtc {
    int natoms = 0;
    int nframes = 0; 
    vector<double> time;
    vector<vector<double>> box;
    vector<vector<vector<double>>> coord;
};

s_xtc ReadXTC(const string &filename);

#endif
