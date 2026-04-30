// analyze.hpp
#ifndef ANALYZE_HPP
#define ANALYZE_HPP

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <regex>
#include <random>

#include <pdb_tools.hpp>
#include <psf_tools.hpp>
#include <prmtop_tools.hpp>
#include <pos_tools.hpp>
#include <traj_tools.hpp>
#include <dcd_r.hpp>
#include <output.hpp>
#include <read_file.hpp>
#include <protein_tools.hpp>
#include <nc_tools.hpp>
#include <select_tools.hpp>
#include <math_tools.hpp>
#include <extension.hpp>
#include <itp_tools.hpp>
#include <xtc_tools.hpp>
#include <make_hist.hpp>

using namespace std;

// Structure --->
struct s_bond {
    double c0, c1;
};

struct s_charge {
   int nlist = 0;
    vector<string> name;
    vector<string> atmnm;
    vector<string> resnm;
    vector<double> value;
};

struct s_joint {
    vector<string> resnm;
    vector<string> atmnm;
};

struct s_connect {
    vector<int> joint1, joint2;
    vector<string> joint_resid1, joint_resid2;
};

struct s_pair {
    vector<int> list1, list2;
};


struct s_bond_set {
  int                    id1, id2;
  vector<string>         id1_set, id2_set;
  vector<vector<string>> id1_bond, id2_bond;
  vector<int>            id1_index, id2_index;
};
// <---

s_pdb ChangeLength(s_pdb monopdb,
                   const double joint_dist,
                   const vector<int>& joint_ref,
                   const vector<int>& joint_mov);

void SetupFitConnect(const s_pdb& polpdb,
                     const s_pdb& monopdb,
                     const vector<int>& joint_ref,
                     const vector<int>& joint_mov,
                     s_fit& fit);

s_pdb OperateTrrotConnect(s_fit& fit,
                          s_pdb& polpdb,
                          s_pdb& monopdb, 
                          const vector<int>& joint_ref, 
                          const vector<int>& joint_mov, 
                          const int resid);

void MAKE_POLYMER(const string& pdbfile,
                  const string& outfile,
                  const string& chargefile,
                  const string& outcharge,
                  double& joint_dist,
                  vector<int>& joint_ref,
                  vector<int>& joint_mov,
                  int& npol, 
                  const string& captype, 
                  const double netcharge);

void MAKE_TOPOLOGY(const string& itpfile,
                   const string& outitp,
                   vector<int>& joint_mov,
                   vector<int>& joint_ref,
                   const int npol, 
                   const s_bond& bond);


void MAKE_SYS_TOP(const vector<string>& itpfile,
                  const string& outitp,
                  vector<int>& joint_mov,
                  vector<int>& joint_ref,
                  const int npol,
                  const int Np,
                  const int Nc, 
                  const s_bond& bond);

void MAKE_CROSSLINKING(const vector<string>& itpfile,
                       const string& trajfile,
                       const string& outitp,
                       const vector<string>& chargefile,
                       const vector<string>& selpol,
                       const vector<string>& selcross,
                       const vector<string>& joint_string,
                       const int Np,
                       const int Nc,
                       const double rc, 
                       const string rule);

void COMBINE_POL_CROSS(const vector<string>& itpfile,
                       const string& outitp,
                       const int Np,
                       const int Nc);

void PolymerConvertItp(const string& itpfile,
                       const string& outitp,
                       const string& resid,
                       const vector<string>& atomlist,
                       const vector<int>& joint_ref,
                       const vector<int>& joint_mov,
                       const int npol);

void PolymerConvertPdb(const string& pdbfile,
                       const string& outpdb,
                       const string& resid,
                       const vector<string>& atomlist,
                       const vector<int>& joint_ref,
                       const vector<int>& joint_mov,
                       const int npol);

void REMOVE_ANALYSIS(const string& itpfile,
                     const string& outitp,
                     const string& pdbfile, 
                     const string& outpdb, 
                     const vector<string>& selpol,
                     const vector<string>& selcross,
                     const int Nc);

void REMOVE_RESID_ANALYSIS(const string& itpfile,
                           const string& outitp,
                           const string& pdbfile,
                           const string& outpdb,
                           const vector<string>& resid_list);

void CROSS_LINKING_CHECK(const string& itpfile,
                         const vector<string>& selpol,
                         const vector<string>& selcross,
                         const int Nc);

void ANALYZE_BOND_AND_REMOVE(const vector<string>& pdbfile,
                             const string& outpdb,
                             const vector<string>& itpfile, 
                             const string& outitp,
                             const string& chargefile,
                             const vector<string>& joint_string, 
                             const vector<string>& joint_string_polymer, 
                             const double distance, 
                             optional<int> ncycle=1000000, 
                             optional<double> prob=1.0);

void ANALYZE_BOND_AND_REMOVE_AC(const vector<string>& pdbfile,
                                const string& outpdb,
                                const vector<string>& itpfile, 
                                const string& outitp,
                                const string& chargefile, 
                                const vector<string>& joint_string, 
                                const vector<string>& joint_string_polymer, 
                                const vector<string>& joint_string_cross, 
                                const vector<string>& joint_string_remove, 
                                const int nbond, 
                                const double distance, 
                                optional<int> ncycle=1000000);

void MAKE_CROSSLINKING_DRY(const vector<string>& itpfile,
                           const string& trajfile,
                           const string& outitp,
                           const string& pdbfile,
                           const string& outpdb,
                           const vector<string>& chargefile,
                           const vector<string>& selpol,
                           const vector<string>& selcross,
                           const vector<string>& drypol,
                           const vector<string>& drycross,
                           const vector<string>& joint_string,
                           const int Np,
                           const int Nc,
                           const double rc);

void MAKE_CROSSLINKING_RADICAL(const vector<string>& itpfile,
                               const string& trajfile,
                               const string& outitp,
                               const vector<string>& chargefile,
                               const vector<string>& selpol,
                               const vector<string>& selcross,
                               const int Np,
                               const int Nc,
                               const double rc);

void ANALYZE_REMOVE_IONS(const string& itpfile,
                         const string& pdbfile,
                         const string& outitp,
                         const string& outpdb,
                         const vector<string>& ions_string);

void ANALYZE_CL_POSITION(const string& itpfile,
                         const vector<string>& trajlist,
                         const vector<string>& selpol,
                         const vector<string>& selcross,
                         const int Nc, 
                         const vector<double>& Min,
                         const vector<double>& Max,
                         const int& grid);

void ANALYZE_HOPPING(const string& topfile,
                     const vector<string>& trajlist,
                     const vector<string>& selion,
                     const vector<string>& selmemb, 
                     const double dt, 
                     const double rc,
                     const int rst_sign);

void ANALYZE_SD(const string& topfile,
                const string& trajfile,
                const vector<string>& sel,
                const int Next, 
                const double dt);

void ANALYZE_COORD_NUM(const string& topfile,
                       const vector<string>& trajlist,
                       const vector<string>& sel1,
                       const vector<string>& sel2,
                       const double& rc);

void ANALYZE_END_TO_END(const string& topfile,
                        const vector<string>& trajlist,
                        const vector<string>& selpol, 
                        const int Np,
                        const int npol, 
                        const int nbin);

void ANALYZE_MASS(const string& itpfile);

void CHANGE_CHARGE(const string& itpfile,
                   const string& outitp,
                   const vector<string>& sel, 
                   const double& aft_charge);

// inline function --->
inline vector<string> return_atom_set(const s_itp itp, const int idx) {
    //
    string resid       = itp.atoms.resid[idx-1];
    string resid_trim  = remove_digits(resid);
    string atom        = itp.atoms.atom[idx-1];
    vector<string> set = {atom, resid_trim};
    //
    return set;
}

inline s_charge prepare_charge(const string& chargefile) {
    //
    auto c_tmp = ReadChargefile(chargefile);
    s_charge charge;
    charge.name = move(c_tmp.first);
    charge.value = move(c_tmp.second);
    charge.nlist = charge.name.size();
    //
    for (int i = 0; i < charge.nlist; ++i) {
        vector<string> name_split = split(charge.name[i]);
        charge.atmnm.push_back(name_split[0]);
        charge.resnm.push_back(name_split[1]);
    }
    return charge; 
}
// <----

#endif
