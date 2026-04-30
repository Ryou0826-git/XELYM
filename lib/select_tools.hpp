// selection_tools.hpp

#ifndef SELECTION_TOOLS_HPP
#define SELECTION_TOOLS_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include <pdb_tools.hpp>

using namespace std;

struct s_sel {
    optional<int>    resid;
    optional<string> atmnm;
    optional<string> resnm;
};

class SELOPR {
  public:
  //
  // Research selection
  s_sel ParseSel(const string& expr);

  vector<int> ResearchPdb(const s_pdb& pdb, 
                          const s_sel& sel);
};

vector<s_sel> ExtractSel(const string& sel);

#endif
