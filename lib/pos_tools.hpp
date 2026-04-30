// pos_tools.hpp
#ifndef POS_TOOLS_HPP
#define POS_TOOLS_HPP

#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <string>

#include <dcd_r.hpp>
#include <traj_tools.hpp>
#include <psf_tools.hpp>
#include <pdb_tools.hpp>

using namespace std;

struct s_fit {
    int natm;
    vector<double>         mass;
    vector<vector<double>> refcoord;
    vector<vector<double>> movcoord;
    vector<double>         refcom;
    vector<double>         movcom;
    Eigen::MatrixXd        rot_mat;
    s_fit() 
        : refcom(3, 0.0f), 
          movcom(3, 0.0f) {
    }
};

//
// inline function --->
inline double pbc_diff(double dx, double box) {
    return dx - box * round(dx / box);
}

inline double pbc_wrap(double x, double box) {
    return x - box * floor(x / box);
}

inline double dist2_pbc(const vector<double>& r1,
                        const vector<double>& r2,
                        const vector<double>& box){
    //
    double dx = pbc_diff(r2[0] - r1[0], box[0]);
    double dy = pbc_diff(r2[1] - r1[1], box[1]);
    double dz = pbc_diff(r2[2] - r1[2], box[2]);
    //
    return dx*dx + dy*dy + dz*dz;
}

inline double pbc_midpoint(double x1, double x2, double box) {
    double dx = pbc_diff(x2 - x1, box);
    return x1 + 0.5 * dx;
}

inline vector<double> pbc_midpoint_3d(const vector<double>& r1,
                                      const vector<double>& r2,
                                      const vector<double>& box) {
    vector<double> mid(3);
    //
    for (int d = 0; d < 3; ++d) {
        double dr = pbc_diff(r2[d] - r1[d], box[d]);
        mid[d] = r1[d] + 0.5 * dr;
    }
    return mid;
}
// <-----

vector<vector<double>>
CalcAveStructure(const vector<vector<vector<double>>>& coord);

vector<vector<double>>
CalcCoM(const vector<double>& mass,
        const vector<vector<vector<double>>>& coord);

vector<vector<vector<double>>>
CalcCoMres(const vector<vector<int>>& molinfo, 
           const s_psf& psf, 
           const vector<vector<vector<double>>>& coord);

vector<vector<double>>
CalcCoMresOneframe(const vector<vector<int>>& molinfo, 
                    const s_psf& psf, 
                    const vector<vector<double>>& coord);

void RefSetupFit(const s_psf& psf, 
                 const vector<vector<int>>& psfinfo, 
                 const vector<vector<double>>& refcoord, 
                 const vector<vector<int>>& refinfo, 
                 s_fit& fit);

void SetupFitConnect(const s_pdb& polpdb,
                     const s_pdb& monopdb,
                     const vector<int>& joint_st,
                     const vector<int>& joint_ed,
                     s_fit& fit);

void RefSetupFitNeo(const vector<double>& mass, 
                    const vector<int>& massinfo, 
                    const vector<vector<double>>& refcoord, 
                    const vector<int>& refinfo, 
                    s_fit& fit);

void SetupFit(const vector<vector<int>>& psfinfo, 
              const vector<vector<double>>& coord, 
              s_fit& fit);

void GetTrrot(s_fit& fit);

void ChangeMovcom(s_fit& fit, 
                  const s_psf& psf, 
                  const vector<vector<int>>& psfinfo_fit, 
                  const vector<vector<double>>& coord);

void InverseFit(s_fit& fit);

void OperateTrrot(s_fit& fit, 
                  vector<vector<double>>& coord);

s_pdb OperateTrrotConnect(s_fit& fit, 
                          s_pdb& polpdb, 
                          s_pdb& monopdb);

void ShiftOrigin(const s_fit& fit, 
                 vector<vector<double>>& coord);

void ShiftOriginAll(const s_psf& psf, 
                    vector<vector<double>>& coord);

void WrapCoordinatesToBoxCenter(vector<vector<vector<double>>>& coord,
                                const vector<vector<double>>& box);

vector<vector<vector<double>>>
GetFitCoord(vector<vector<vector<double>>>& coord, s_fit& fit);

vector<vector<vector<double>>>
StoreFitTraj(s_fit& fit,
             const int& frames,
             const s_dcd& dcd,
             const vector<vector<int>>& psfinfo_fit,
             const vector<vector<int>>& psfinfo_get);

void StoreFitTrajAllAtm(s_fit& fit,
                        const int& frames,
                        s_dcd& dcd,
                        const vector<vector<int>>& psfinfo);

void GetFitAllAtm(const string& dcdfile,
                  const int& sumatm,
                  const vector<vector<int>>& psfinfo_fit,
                  s_fit& fit,
                  s_dcd& dcd);

vector<vector<vector<double>>>
GetFitcoordall(const vector<string>& dcdfile,
               const int& sumatm,
               const vector<vector<int>>& psfinfo_fit,
               const vector<vector<int>>& psfinfo_get,
               s_fit& fit);

vector<vector<vector<double>>>
Extractcoordall(const vector<string>& dcdfile,
                const int& sumatm,
                const vector<vector<int>>& molinfo);

vector<vector<int>> GenerateNewID(vector<vector<int>>& molinfo);

vector<vector<vector<double>>>
ExtractAtmnmCoord(const vector<vector<vector<double>>>& coord,
                  const s_branch& branch,
                  const string& atmnm);

pair<vector<vector<vector<double>>>, s_branch>
OperateGeneralFit(const string& reffile,
                  const string& psffile,
                  const vector<string>& dcdfile,
                  string& segnm,
                  const string& fitatmnm);

vector<double> CalcThickness(const vector<vector<vector<double>>>& coord);

#endif
