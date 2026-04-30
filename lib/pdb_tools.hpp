// pdb_tools.hpp
#ifndef PDB_TOOLS_HPP
#define PDB_TOOLS_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <boost/algorithm/string.hpp>

using namespace std;

struct s_pdb {
    vector<string>         atmnm;
    vector<string>         resnm;
    vector<string>         chainid;
    vector<int>            resid;
    vector<vector<double>> coord;
    vector<double>         occupancy;
    vector<double>         tempFactor;
    vector<double>         box;
    vector<string>         segnm;
    vector<string>         katm;
};

inline void trim(string &s) {
    s.erase(s.begin(),
            find_if(s.begin(), s.end(),
                         [](unsigned char ch) { return !isspace(ch); }));

    s.erase(find_if(s.rbegin(), s.rend(),
                         [](unsigned char ch) { return !isspace(ch); }).base(),
            s.end());
}

class PDBOPR {
  public:
    //
    // Read file section
    vector<vector<string>> 
    Readpdball(const string& pdb_file);
    s_pdb Loadfrompdb(const string& pdbfile, 
                      const int& natm);
    s_pdb LoadfrompdbNeo(const string& pdbfile);
    s_pdb LoadfrompdbAmber(const string& pdbfile);
    s_pdb Removepdb(const s_pdb& pdb, const vector<int> nonbond_list);

    //
    // getmolinfo. section
    s_pdb ExtractSegnm(const s_pdb& refspdb, 
                       const string& resnm);
    vector<vector<int>> 
    getmolinfo(const s_pdb& pdb, 
               const vector<int>& reslist);
    vector<vector<int>> 
    getmolinfo_atmnm(const s_pdb& pdb, 
                     const vector<int>& reslist, 
                     const string& atmnm);

    vector<vector<int>> 
    getmolinfo_atmnm_all(const s_pdb& pdb, 
                         const string& atmnm);
    
    vector<vector<int>> 
    getmolinfo_atmnm_all_amber(const s_pdb& pdb, 
                               const string& atmnm);
    
    vector<vector<int>> SearchProtein(const s_pdb& pdb);
    s_pdb Extractfrominfo(const s_pdb& refspdb, 
                          const vector<vector<int>>& molinfo);
    s_pdb ExtractfrominfoNOT(const s_pdb& refspdb, 
                             const vector<vector<int>>& molinfo);

    //
    // Operate coord section
    void ScaleMOL(s_pdb& pdb, const double& scale);
    void PlungeCoord(s_pdb& pdb, const vector<vector<double>>& coord);
    
    //
    // Out section
    void Outpdb(const vector<vector<string>>& line, 
                const vector<vector<double>>& coord, 
                const string& filename);
    void CatandOutpdb(const s_pdb& pdb1, 
                      const s_pdb& pdb2, 
                      const vector<vector<string>>& line1, 
                      const vector<vector<string>>& line2, 
                      const string& filename);
    void OutInsertpdb(const vector<vector<double>>& coord, 
                      const vector<vector<string>>& sysline, 
                      const vector<vector<string>>& ligline, 
                      const string& filename);
    void Outxyz(const s_pdb& pdb, 
                const vector<vector<double>>& coord, 
                const string& filename);
    void OutPdbNeo(const vector<s_pdb>& pdblist, const string& filename);
    void OutPdbAmber(const vector<s_pdb>& pdblist, const string& filename);
    void OutPdbAmber2(const vector<s_pdb>& pdblist, const string& filename);
    
};

#endif
