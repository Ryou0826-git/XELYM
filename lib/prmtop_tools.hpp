// prmtop_tools.hpp
#ifndef PRMTOP_TOOLS_HPP
#define PRMTOP_TOOLS_HPP

#include <vector>
#include <string>

using namespace std;

struct s_prmtop {
    vector<double> charge;
    vector<double> mass;
    vector<double> Alj;
    vector<double> Blj;
};

class PRMTOPOPR {
  public:
    //
    // Read file section...
    s_prmtop Loadfromprmtop(const string& prmtopfile);
};

#endif
