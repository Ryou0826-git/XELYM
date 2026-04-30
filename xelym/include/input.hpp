// input.hpp

#ifndef INPUT_HPP
#define INPUT_HPP

#include "json.hpp"

#include <string>
#include <vector>
#include <optional>
#include <iostream>
#include <fstream>

using namespace std;

inline void PrintExample() {
    cerr << "Usage: ./xelym.x <input (json format)>"                                      << endl;
    cerr << endl;
    cerr << "========================= Example of inputfile ========================"     << endl;
    cerr << "{"                                                                           << endl;
    cerr << "\"mode\": \"crosslinking-radical\","                                         << endl;
    cerr << "\"trajfile\": \"../02_nvt/nvt.xtc\","                                        << endl;
    cerr << "\"itpfile\": [\"../sys/nipam_bis.itp\", \"../data_set/itp/imit.itp\"],"      << endl;
    cerr << "\"outitp\": \"nipam_bis_temp.itp\","                                         << endl;
    cerr << "\"selpol\": [\"C1_Ni\", \"C2_Ni\"],"                                         << endl;
    cerr << "\"selcross\": [\"C2_BIS\", \"C1_BIS\", \"C6_BIS\", \"C7_BIS\"],"             << endl;
    cerr << "\"Np\": 100,"                                                                << endl;
    cerr << "\"Nc\": 15,"                                                                 << endl;
    cerr << "\"rc\": 4.0"                                                                 << endl;
    cerr << "}"                                                                           << endl;
    cerr << "======================================================================="     << endl;
    cerr << endl;
    cerr << "If you want to use other options, check source code";
    cerr << " or contanct to developer (R. Okabe e-mail:okabe-ry@cheng.es.osaka-u.ac.jp)" << endl;
}


class Config {
  public:
    //
    string                   mode;
    optional<vector<string>> pdbfile;
    optional<string>         outpdb;
    optional<vector<string>> chargefile;
    optional<string>         outcharge;
    optional<vector<string>> itpfile;
    optional<string>         outitp;
    optional<string>         topfile;
    optional<string>         itpimitfile;
    optional<string>         trajfile;
    optional<vector<string>> trajlist;
    optional<string>         captype;
    optional<int>            npol;
    optional<double>         joint_dist;
    optional<vector<int>>    joint_ref;
    optional<vector<int>>    joint_mov;
    optional<vector<string>> joint_string;
    optional<vector<string>> joint_string_polymer;
    optional<vector<string>> joint_string_cross;
    optional<vector<string>> joint_string_remove;
    optional<int>            Np;
    optional<int>            Nc;
    optional<int>            Next;
    optional<int>            Ni;
    optional<double>         dt;
    optional<vector<string>> selpol;
    optional<vector<string>> selcross;
    optional<vector<string>> selion;
    optional<vector<string>> selmemb;
    optional<vector<string>> sel;
    optional<vector<string>> sel1;
    optional<vector<string>> sel2;
    optional<vector<string>> drypol;
    optional<vector<string>> drycross;
    optional<double>         rc;
    optional<double>         bond_c0;
    optional<double>         bond_c1;
    optional<double>         netcharge;
    optional<vector<string>> atomlist;
    optional<string>         resid;
    optional<string>         rule;
    optional<int>            nbond;
    optional<double>         distance;
    optional<int>            ncycle;
    optional<double>         prob;
    optional<vector<string>> ions_string;
    optional<int>            rst_sign;
    optional<vector<double>> Min;
    optional<vector<double>> Max;
    optional<int>            grid;
    optional<double>         aft_charge;
    optional<int>            nbin;
    //
    void print() const;
};

Config load_config(const string& filename);

#endif
