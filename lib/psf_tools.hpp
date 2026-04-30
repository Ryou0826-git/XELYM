// psf_tools.hpp
#ifndef PSF_TOOLS_HPP
#define PSF_TOOLS_HPP

#include <cstdlib>
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <utility>
#include <vector>
#include <string>

using namespace std;

struct s_psf {
     vector<string>      segname;
     vector<string>      resname;
     vector<string>      atmname;
     vector<string>      katm;
     vector<int>         molid;
     vector<double>      mass;
     vector<double>      charge;
     vector<vector<int>> bond;
};

struct s_branch {
    vector<vector<int>> molinfo;
    vector<vector<string>> atmnm;
    vector<vector<string>> resnm;
};

struct s_reschain {
    vector<vector<int>> resinfo;
    vector<int>         residinfo;
};

class PSFOPR {
  public:
    //
    // Read file section
    s_psf loadFile(const string& filename);

    //
    // getmolinfo. section
    vector<vector<int>> 
    getmolinfo(const s_psf& psf, 
               const string& resnm);
    vector<vector<int>> 
    getmolinfo_segnm(const s_psf& psf, 
                     const string& segnm);
    vector<vector<int>> 
    getmolinfo_endtoend(const s_psf& psf, 
                        const string& segnm, 
                        const vector<string>& selatm);
    vector<vector<int>> 
    getmolinfo_segres(const s_psf& psf, 
                      const string& segnm, 
                      const vector<int>& reslist);
    vector<vector<int>> 
    getmolinfo_seg_all(const s_psf& psf, 
                       const string& segnm);
    vector<vector<int>> 
    getmolinfo_segatm_res(const s_psf& psf, 
                          const string& segnm, 
                          const string& atmnm, 
                          const vector<int>& reslist);
    vector<vector<int>> 
    getmolinfo_segatm_all(const s_psf& psf, 
                          const string& pronm, 
                          const string& atmnm);
    s_branch getmolinfo_segnm_heavy(const s_psf& psf, 
                                    const string& segnm);
    s_reschain SplitListResid(const s_psf& psf, const vector<int>& list);
    //
    // bond section
    void LoadpsffileBond(s_psf& psf, const string& psffile);
};

#endif
