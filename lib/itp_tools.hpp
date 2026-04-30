// itp_tools.hpp

#ifndef ITP_TOOLS_HPP
#define ITP_TOOLS_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cctype>
#include <iomanip>
#include <algorithm>
#include <unordered_set>

using namespace std;

//-----------------------------------------------------------------------------
// Structure of itp
//-----------------------------------------------------------------------------
struct s_atoms {
    vector<int>    nr;
    vector<string> type;
    vector<int>    resnr;
    vector<string> resid;
    vector<string> atom;
    vector<int>    cgnr;
    vector<double> charge;
    vector<double> mass; 
    int            n;
};

struct s_bonds {
    vector<int>    ai;
    vector<int>    aj;
    vector<int>    funct;
    vector<double> c0;
    vector<double> c1;
};

struct s_pairs {
    vector<int> ai;
    vector<int> aj;
    vector<int> funct;
};

struct s_angles {
    vector<int>    ai;
    vector<int>    aj;
    vector<int>    ak;
    vector<int>    funct;
    vector<double> angle;
    vector<double> fc;
};

struct s_dihedrals {
    vector<int>    ai;
    vector<int>    aj;
    vector<int>    ak;
    vector<int>    al;
    vector<int>    funct;
    vector<double> ph0;
    vector<double> cp;
    vector<int>    mult;
    vector<double> angle;
    vector<double> fc;
    //
    // AMBER FIELD
    vector<double> phase;
    vector<double> kd;
    vector<int>    pn;
    // OPLS, funct = 3
    vector<double> c0;
    vector<double> c1;
    vector<double> c2;
    vector<double> c3;
    vector<double> c4;
    vector<double> c5;
};

struct s_exclusion {
    vector<int> ai;
    vector<int> aj;
};

struct s_itp {
   s_atoms     atoms;
   s_bonds     bonds;
   s_pairs     pairs;
   s_angles    angles;
   s_dihedrals dihedrals;
   s_exclusion exclusion;
};
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Structure of itp_imit
//-----------------------------------------------------------------------------
struct s_bonds_imit {
    vector<string> ai;
    vector<string> aj;
    vector<int>    funct;
    vector<double> c0;
    vector<double> c1;
};

struct s_angles_imit {
    vector<string> ai;
    vector<string> aj;
    vector<string> ak;
    vector<int>    funct;
    vector<double> angle;
    vector<double> fc;
};

struct s_dihedrals_imit {
    vector<string> ai;
    vector<string> aj;
    vector<string> ak;
    vector<string> al;
    vector<int>    funct;
    vector<double> ph0;
    vector<double> cp;
    vector<int>    mult;
    vector<double> angle;
    vector<double> fc;
    vector<double> phase;
    vector<double> kd;
    vector<int>    pn;
    // OPLS, funct = 3
    vector<double> c0;
    vector<double> c1;
    vector<double> c2;
    vector<double> c3;
    vector<double> c4;
    vector<double> c5;
};

struct s_itp_imit {
   s_bonds_imit     bonds;
   s_angles_imit    angles;
   s_dihedrals_imit dihedrals;
};

//-----------------------------------------------------------------------------
struct s_dihed_parts {
    double angle;
    double cp;
    double fc;
    double ph0;
    double mult;
    double phase;
    double kd;
    int    pn;
    double c0, c1, c2, c3, c4, c5;
};
//-----------------------------------------------------------------------------

static inline string trim(const string &s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
};

enum Section {
    NONE, ATOMS, BONDS, PAIRS, ANGLES, DIHEDRALS, EXCLUSIONS
};

///////////////////////////////////////////////////////////////////////////////
void Readitpfile(const string &filename, s_itp &itp);

void RemoveIndex(s_itp &itp, int target);

void RemoveAtoms(s_itp& itp,
                 const std::vector<int>& atoms_to_remove);

void RemoveAtomsKeepIndex(s_itp& itp, const std::vector<int>& atoms_to_remove);

void RemoveConnections(s_itp& itp, int target);

void CompactAtomIndices(s_itp &itp);

void Writeitp(const s_itp &itp, const string &filename);

s_itp CombineItp(const s_itp& A, const s_itp& B);

void AppendBond(s_itp& itp,
                const int idx1,
                const int idx2,
                const int funct, 
                const double c0,
                const double c1);

void AppendAngle(s_itp& itp,
                 const int idx1,
                 const int idx2,
                 const int idx3,
                 const int funct,
                 const double angle,
                 const double fc);

void AppendDihedral(s_itp& itp,
                    const int idx1,
                    const int idx2,
                    const int idx3,
                    const int idx4,
                    const int funct,
                    const s_dihed_parts dihed_parts);

void AppendPair(s_itp& itp,
                const int idx1,
                const int idx2,
                const int funct);

int CountBonds(const s_itp &itp, int atom_index);

vector<string> GetBondedResidues(const s_itp& itp, int atom_index);

vector<string> GetbondResnmandResid(const s_itp& itp,
                                    int atom_index,
                                    string polnm);

pair<vector<vector<string>>, vector<int>> GetBondedAtomandResid(const s_itp itp,
                                                                int atom_index);

vector<int> GetBondedIndex(const s_itp& itp,
                           int atom_index);

bool is_bonded(const s_itp& itp, int atom1, int atom2);

void Readitpimitfile(const string &filename,
                     s_itp_imit &itp_imit);

#endif
