// inpcrd_tools.hpp
#ifndef INPCRD_TOOLS_HPP
#define INPCRD_TOOLS_HPP

#include <boost/algorithm/string.hpp>
#include <pdb_tools.hpp>

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iomanip>

using namespace std;

struct s_inpcrd {
    string                 title;
    int                    natm;
    vector<vector<double>> coord;
    vector<double>         box;
    vector<double>         angle;
};

class INPCRDOPR {
  public:
    s_inpcrd ReadInfo(const string& filename);

    s_inpcrd PdbToInpcrd(const s_pdb& pdb);

    void OutInpcrd(const vector<s_inpcrd>& inpcrd, 
                   const string& filename);
    void PlungeBoxAngle(const s_inpcrd& inpcrd_origin, 
                        s_inpcrd& inpcrd);
};

#endif
