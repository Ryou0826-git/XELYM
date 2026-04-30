// pos_tools.cpp

#include "psf_tools.hpp"
#include "prmtop_tools.hpp"
#include "math_tools.hpp"
#include "pos_tools.hpp"
#include "omp.h"

#include <cstdlib>
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <cmath>
#include <numeric>

using namespace std;

vector<vector<double>> 
CalcAveStructure(const vector<vector<vector<double>>>& coord) {
    int frames = coord.size();
    int natms  = coord[0].size();
    vector<vector<double>> coordave(natms, vector<double>(3, 0.0));

    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < natms; ++j) {
            for (int k = 0; k < 3; ++k) {
                coordave[j][k] += coord[i][j][k];
            }
        }
    }
   
    for (int i = 0; i < natms; ++i) {
        for (int j = 0; j < 3; ++j) {
            coordave[i][j] /= frames;
        }
    }

    return coordave;
}

vector<vector<double>> 
CalcCoM(const vector<double>& mass, 
        const vector<vector<vector<double>>>& coord) {
    //
    int frames = coord.size();
    int natms = coord[0].size();
    double summass = 0.0;
    for (int i = 0; i < natms; ++i) {
        summass += mass[i];
    }
    //
    vector<vector<double>> CoMcoord;
    for (int f = 0; f < frames; ++f) {
        vector<double> cellcoord(3, 0.0);
        for (int n = 0; n < natms; ++n) {
            for (int k = 0; k < 3; ++k) {
                cellcoord[k] += mass[n] * coord[f][n][k];
            }
        }
        //
        for (int k = 0; k < 3; ++k) {
            cellcoord[k] /= summass;
        }
        CoMcoord.push_back(cellcoord);
    }

    return CoMcoord;
}

vector<vector<vector<double>>>
CalcCoMres(const vector<vector<int>>& molinfo, 
           const s_psf& psf, 
           const vector<vector<vector<double>>>& coord) {
    int frames = coord.size();
    //
    int numres = molinfo.size();
    int nummol = molinfo[0].size();
    double summass = 0;
    //
    for (int i = 0; i < nummol; ++i) {
        int j = molinfo[0][i];
        summass = summass + psf.mass[j-1];
    }

    //
    // Allocate Memory of CoM.
    int dim = 3;
    vector<vector<vector<double>>>\
    CoMres(frames, vector<vector<double>>(numres, vector<double>(3, 0.0)));

    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < numres; ++j) {
            for (int k = 0; k < nummol; ++k) {
                int l = molinfo[j][k];
                for (int m = 0; m < dim; ++m) {
                  CoMres[i][j][m] = CoMres[i][j][m] + psf.mass[l-1] * coord[i][l-1][m];
                }
            }
            //
            for (int m = 0; m < dim; ++m) {
                CoMres[i][j][m] = CoMres[i][j][m] / summass;
            }
        }
    }
    return CoMres;
}

vector<vector<double>>
CalcCoMresOneframe(const vector<vector<int>>& molinfo, 
                   const s_psf& psf, 
                   const vector<vector<double>>& coord) {
    //
    int numres = molinfo.size();
    int nummol = molinfo[0].size();
    double summass = 0;
    //
    for (int i = 0; i < nummol; ++i) {
        int j = molinfo[0][i];
        summass = summass + psf.mass[j-1];
    }

    //
    // Allocate Memory of CoM.
    int dim = 3;
    vector<vector<double>> 
    CoMres(numres, vector<double>(3, 0.0));

    for (int j = 0; j < numres; ++j) {
        for (int k = 0; k < nummol; ++k) {
            int l = molinfo[j][k];
            for (int m = 0; m < dim; ++m) {
              CoMres[j][m] = CoMres[j][m] + psf.mass[l-1] * coord[l-1][m];
            }
        }
        //
        for (int m = 0; m < dim; ++m) {
            CoMres[j][m] = CoMres[j][m] / summass;
        }
    }
    return CoMres;
}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Fitting section ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void RefSetupFit(const s_psf& psf, 
                 const vector<vector<int>>& psfinfo, 
                 const vector<vector<double>>& refcoord, 
                 const vector<vector<int>>& refinfo, 
                 s_fit& fit) {
    int nres = psfinfo.size();
    fit.refcoord.resize(nres, vector<double>(3, 0.0));
    fit.mass.resize(nres, 0.0);
    //
    for (int i = 0; i < nres; ++i) {
        int resatm = psfinfo[i].size();
        for (int j = 0; j < resatm; ++j) {
            int label = psfinfo[i][j];
            fit.mass[i] = psf.mass[label-1];
            //fit.mass[i] = 1.0;
        }
    }
    //
    int l = 0;
    for (int i = 0; i < nres; ++i) {
        int resatm = refinfo[i].size();
        for (int j = 0; j < resatm; ++j) {
            int label = refinfo[i][j];
            cout << "Resid (reference): " << label << endl;
            for (int k = 0; k < 3; ++k) {
                fit.refcoord[l][k] = refcoord[label-1][k];
            }
            l += 1;
        }
    }

    fit.natm = nres;

    double summass = 0.0;
    for (int i = 0; i < fit.natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            fit.refcom[j] += fit.mass[i] * fit.refcoord[i][j];
        }
        summass += fit.mass[i];
    }
    //
    cout << "Reference CoM: ";
    for (int i = 0; i < 3; ++i) {
        fit.refcom[i] /= summass;
        cout << fit.refcom[i] << " ";
    }
    cout << endl;
    //
    return;
}

void RefSetupFitNeo(const vector<double>& mass, 
                    const vector<int>& massinfo, 
                    const vector<vector<double>>& refcoord, 
                    const vector<int>& refinfo, 
                    s_fit& fit) {
    //
    for (const int label : massinfo) {
        fit.mass.push_back(mass[label-1]);
    }

    for (const int label : refinfo) {
        cout << "Resid (reference): " << label << endl;
        double x = refcoord[label-1][0];
        double y = refcoord[label-1][1];
        double z = refcoord[label-1][2];
        fit.refcoord.push_back({x, y, z});
    }

    fit.natm = refinfo.size();

    double summass = 0.0;
    for (int i = 0; i < fit.natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            fit.refcom[j] += fit.mass[i] * fit.refcoord[i][j];
        }
        summass += fit.mass[i];
    }
    //
    cout << "Reference CoM: ";
    for (int i = 0; i < 3; ++i) {
        fit.refcom[i] /= summass;
        cout << fit.refcom[i] << " ";
    }
    cout << endl;
    //
    return;
}

//void SetupFitConnect(const s_pdb& polpdb,
//                     const s_pdb& monopdb,
//                     const vector<int>& joint_st, 
//                     const vector<int>& joint_ed, 
//                     s_fit& fit) {
//    //
//    if (joint_st.size() != joint_ed.size()) {
//        cout << endl;
//        cout << "Error: fitting selection is not correct..." << endl;
//    } else {
//        fit.natm = joint_st.size();
//    }
//    
//    for (int i = 0; i < fit.natm; ++i) {
//        fit.mass.push_back(1.0);
//    }
//
//    for (int j = joint_st.size(); j > 0; --j) {
//        int i = joint_st[j-1];
//        double x = monopdb.coord[i-1][0];
//        double y = monopdb.coord[i-1][1];
//        double z = monopdb.coord[i-1][2];
//        fit.movcoord.push_back({x, y, z});
//        fit.movcom[0] += x;
//        fit.movcom[1] += y;
//        fit.movcom[2] += z;
//    }
//
//    for (int i = 0; i < 3; ++i) {
//        fit.movcom[i] /= joint_st.size();
//    }
//    //cout << fit.movcom[0] << " " << fit.movcom[1] << " " << fit.movcom[2] << endl;
//
//    // for ref
//    for (int j : joint_ed) {
//        double x = polpdb.coord[j-1][0];
//        double y = polpdb.coord[j-1][1];
//        double z = polpdb.coord[j-1][2];
//        fit.refcoord.push_back({x, y, z});
//        fit.refcom[0] += x;
//        fit.refcom[1] += y;
//        fit.refcom[2] += z;
//        //cout << x << " " << y << " " << z << endl;
//    }
//
//    for (int i = 0; i < 3; ++i) {
//        fit.refcom[i] /= joint_ed.size();
//    }
//    //cout << fit.refcom[0] << " " << fit.refcom[1] << " " << fit.refcom[2] << endl;
//}
//----------------------------------------------------------------------------

void SetupFit(const vector<vector<int>>& psfinfo, 
              const vector<vector<double>>& coord, 
              s_fit& fit) {
    int natm = fit.natm;
    fit.movcoord.resize(natm, vector<double>(3, 0.0));
    //
    for (int i = 0; i < 3; ++i) {
        fit.movcom[i] = 0.0;
    }
    //
    int nres = psfinfo.size();
    int k = 0;
    //
    for (int i = 0; i < nres; ++i) {
        int resatm = psfinfo[i].size();
        for (int j = 0; j < resatm; ++j) {
            int label = psfinfo[i][j];
            //cout << label << endl;
            for (int l = 0; l < 3; ++l) {
                fit.movcoord[k][l] = coord[label-1][l];
                //cout << fit.movcoord[k][l] << " ";
            }
            //cout << endl;
            k += 1;
        }
    }
    //
    double summass = 0.0;
    for (int i = 0; i < fit.natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            fit.movcom[j] += fit.mass[i] * fit.movcoord[i][j];
            //fit.movcom[j] += fit.movcoord[i][j];
        }
        summass += fit.mass[i];
    }
    //
    for (int i = 0; i < 3; ++i) {
        fit.movcom[i] /= summass;
    }
    //
    return;
}

void SetupFitNonsel(const vector<vector<double>>& coord, 
                    s_fit& fit) {
    //
    int natm = fit.natm;
    if (natm != static_cast<int>(coord.size())) {
        cout << "Error: number of atoms to fit is incorrect!!" << endl;
        return;
    }
    //
    fit.movcoord.resize(natm, vector<double>(3, 0.0));
    //
    for (int i = 0; i < 3; ++i) {
        fit.movcom[i] = 0.0;
    }
    //
    for (int i = 0; i < natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            fit.movcoord[i][j] = coord[i][j];
        }
    }
    //
    double summass = 0.0;
    for (int i = 0; i < fit.natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            fit.movcom[j] += fit.mass[i] * fit.movcoord[i][j];
            //fit.movcom[j] += fit.movcoord[i][j];
        }
        summass += fit.mass[i];
    }
    //
    for (int i = 0; i < 3; ++i) {
        fit.movcom[i] /= summass;
    }
    //
    return;
}

void GetTrrot(s_fit& fit) {

    Eigen::MatrixXd sym_mat(4, 4);
    sym_mat.setZero();

    Eigen::Vector3d dref, dmov, dsub, dadd;

    sym_mat = Eigen::MatrixXd::Zero(4, 4);

    for (int i = 0; i < fit.natm; ++i) {
        dref << fit.refcoord[i][0] - fit.refcom[0],
                fit.refcoord[i][1] - fit.refcom[1],
                fit.refcoord[i][2] - fit.refcom[2];
        dmov << fit.movcoord[i][0] - fit.movcom[0],
                fit.movcoord[i][1] - fit.movcom[1],
                fit.movcoord[i][2] - fit.movcom[2];

        dsub = fit.mass[i] * (dref - dmov);
        dadd = fit.mass[i] * (dref + dmov);

        sym_mat(0, 0) += dsub[0] * dsub[0] + dsub[1] * dsub[1] + dsub[2] * dsub[2];
        sym_mat(0, 1) += dadd[1] * dsub[2] - dsub[1] * dadd[2];
        sym_mat(0, 2) += dsub[0] * dadd[2] - dadd[0] * dsub[2];
        sym_mat(0, 3) += dadd[0] * dsub[1] - dsub[0] * dadd[1];

        sym_mat(1, 1) += dsub[0] * dsub[0] + dadd[1] * dadd[1] + dadd[2] * dadd[2];
        sym_mat(1, 2) += dsub[0] * dsub[1] - dadd[0] * dadd[1];
        sym_mat(1, 3) += dsub[0] * dsub[2] - dadd[0] * dadd[2];

        sym_mat(2, 2) += dadd[0] * dadd[0] + dsub[1] * dsub[1] + dadd[2] * dadd[2];
        sym_mat(2, 3) += dsub[1] * dsub[2] - dadd[1] * dadd[2];

        sym_mat(3, 3) += dadd[0] * dadd[0] + dadd[1] * dadd[1] + dsub[2] * dsub[2];
    }
    //
    sym_mat(1, 0) = sym_mat(0, 1);
    sym_mat(2, 0) = sym_mat(0, 2);
    sym_mat(2, 1) = sym_mat(1, 2);
    sym_mat(3, 0) = sym_mat(0, 3);
    sym_mat(3, 1) = sym_mat(1, 3);
    sym_mat(3, 2) = sym_mat(2, 3);
    //
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(sym_mat);

    if (solver.info() != Eigen::Success) {
        cerr << "Eigenvalue decomposition failed!" << endl;
    }
    //
    Eigen::VectorXd eigenvalues  = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();
    //
    Eigen::VectorXd evec = eigenvectors.col(0);
    //
    fit.rot_mat = Eigen::MatrixXd::Zero(3, 3);

    fit.rot_mat(0, 0) = evec(0)*evec(0) + evec(1)*evec(1) 
                        - evec(2)*evec(2) - evec(3)*evec(3);
    fit.rot_mat(0, 1) = 2.0 * (evec(1)*evec(2) + evec(0)*evec(3));
    fit.rot_mat(0, 2) = 2.0 * (evec(1)*evec(3) - evec(0)*evec(2));
    fit.rot_mat(1, 0) = 2.0 * (evec(1)*evec(2) - evec(0)*evec(3));
    fit.rot_mat(1, 1) = evec(0)*evec(0) - evec(1)*evec(1) 
                        + evec(2)*evec(2) - evec(3)*evec(3);
    fit.rot_mat(1, 2) = 2.0 * (evec(2)*evec(3) + evec(0)*evec(1));
    fit.rot_mat(2, 0) = 2.0 * (evec(1)*evec(3) + evec(0)*evec(2));
    fit.rot_mat(2, 1) = 2.0 * (evec(2)*evec(3) - evec(0)*evec(1));
    fit.rot_mat(2, 2) = evec(0)*evec(0) - evec(1)*evec(1) 
                        - evec(2)*evec(2) + evec(3)*evec(3);
    
    return;
}

void InverseFit(s_fit& fit) {
    //
    fit.rot_mat = fit.rot_mat.inverse();
    return;
}

void OperateTrrot(s_fit& fit, 
                  vector<vector<double>>& coord) {
    int natm = coord.size();
    //
    vector<double> dmov(3, 0.0);
    vector<double> rotmov(3, 0.0);
    
    for (int iatm = 0; iatm < natm; ++iatm) {
        //
        // Set coordinate to origin...
        for (int j = 0; j < 3; ++j) {
            dmov[j] = coord[iatm][j] - fit.movcom[j];
        }
        //
        // Rotate coordinate...
        for (int j = 0; j < 3; ++j) {
            rotmov[j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                rotmov[j] += fit.rot_mat(j, k) * dmov[k];
            }
        }
        //
        // Set coordinate to reference origin...
        for (int j = 0; j < 3; ++j) {
            coord[iatm][j] = rotmov[j] + fit.refcom[j];
        }
    }
    return;
}

//s_pdb OperateTrrotConnect(s_fit& fit, 
//                          s_pdb& polpdb, 
//                          s_pdb& monopdb) {
//    s_pdb conn_pdb;
//    int natm = monopdb.atmnm.size();
//
//    vector<double> dmov(3, 0.0);
//    vector<double> rotmov(3, 0.0);
//
//    for (int iatm = 0; iatm < natm; ++iatm) {
//        //
//        // Set coordinate to origin...
//        for (int j = 0; j < 3; ++j) {
//            dmov[j] = monopdb.coord[iatm][j] - fit.movcom[j];
//        }
//        //
//        // Rotate coordinate...
//        for (int j = 0; j < 3; ++j) {
//            rotmov[j] = 0.0;
//            for (int k = 0; k < 3; ++k) {
//                rotmov[j] += fit.rot_mat(j, k) * dmov[k];
//            }
//        }
//        //
//        // Set coordinate to reference origin...
//        for (int j = 0; j < 3; ++j) {
//            monopdb.coord[iatm][j] = rotmov[j] + fit.refcom[j];
//        }
//    }
//
//    return monopdb;
//}

////////////////////////////////////////////////////////////////////////////////

void ChangeMovcom(s_fit& fit, 
                  const s_psf& psf, 
                  const vector<vector<int>>& psfinfo_fit, 
                  const vector<vector<double>>& coord) {
    //
    int    nres    = psfinfo_fit.size();
    double summass = 0.0;

    for (int i = 0; i < 3; ++i) {
        fit.movcom[i] = 0.0;
    }

    for (int i = 0; i < nres; ++i) {
        int resatm = psfinfo_fit[i].size();
        for (int j = 0; j < resatm; ++j) {
            int label = psfinfo_fit[i][j];
            summass += psf.mass[label-1];
            for (int k = 0; k < 3; ++k) {
                fit.movcom[k] += psf.mass[label-1] * coord[label-1][k];
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        fit.movcom[i] /= summass;
    }

    return;
}

void ShiftOriginAll(const s_psf& psf, 
                    vector<vector<double>>& coord) {
    int natm = coord.size();
    vector<double> com(3, 0.0);
    double summass = 0.0f;

    for (int i = 0; i < natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            com[j] += psf.mass[i] * coord[i][j];
        }
        summass += psf.mass[i];
    }

    for (int i = 0; i < 3; ++i) {
        com[i] /= summass;
    }

    for (int i = 0; i < natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            coord[i][j] -= com[j];
        }
    }

    return;
}

void ShiftOrigin(const s_fit& fit, 
                 vector<vector<double>>& coord) {
    int natm = coord[0].size();
    //
    for (int i = 0; i < natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            coord[i][j] -= fit.refcom[j];
        }
    }
    return;
}

void WrapCoordinatesToBoxCenter(vector<vector<vector<double>>>& coord,
                                const vector<vector<double>>& box) {
    //
    int natm = coord[0].size();
    int frames = coord[0][0].size();
    //
    for (int frame = 0; frame < frames; ++frame) {
        vector<double> box_center(3, 0.0);
        //
        for (int dim = 0; dim < 3; ++dim) {
            box_center[dim] = box[frame][dim] / 2.0f;
        }
        //
        for (int atom = 0; atom < natm; ++atom) {
            for (int dim = 0; dim < 3; ++dim) {
                coord[frame][atom][dim] -= box_center[dim];
                coord[frame][atom][dim] = fmod(coord[frame][atom][dim], box[frame][dim]);
                if (coord[frame][atom][dim] < 0) {
                    coord[frame][atom][dim] += box[frame][dim];
                }
                coord[frame][atom][dim] -= box_center[dim];
            }
        }
    }
    return;
}

vector<vector<vector<double>>> 
GetFitCoord(vector<vector<vector<double>>>& coord, 
            s_fit& fit) {
    //
    //
    int frames = coord.size();
    for (int f = 0; f < frames; ++f) {
        SetupFitNonsel(coord[f], fit);
        GetTrrot(fit);
        OperateTrrot(fit, coord[f]);
    }

    return coord;
}

vector<vector<vector<double>>> 
StoreFitTraj(s_fit& fit, 
             const int& frames, 
             const s_dcd& dcd, 
             const vector<vector<int>>& psfinfo_fit, 
             const vector<vector<int>>& psfinfo_get) {
    //
    int sumatm = 0;
    int nget = psfinfo_get.size();
    for (int i = 0; i < nget; ++i) {
       int natm_get = psfinfo_get[i].size();
       for (int j = 0; j < natm_get; ++j) {
           sumatm += 1;
       }
    }

    vector<vector<vector<double>>>
    coord(frames, vector<vector<double>>(sumatm, vector<double>(3, 0.0)));
    //
    for (int iframe = 0; iframe < frames; ++iframe) {
        SetupFit(psfinfo_fit, dcd.coord[iframe], fit);
        GetTrrot(fit);
        //
        int l = 0;
        int nget = psfinfo_get.size();
        for (int i = 0; i < nget; ++i) {
            int natm_get = psfinfo_get[i].size();
            for (int j = 0; j < natm_get; ++j) {
                int label = psfinfo_get[i][j];
                for (int k = 0; k < 3; ++k) {
                    coord[iframe][l][k] = dcd.coord[iframe][label-1][k];
                }
                l += 1;
            }
        }
        OperateTrrot(fit, coord[iframe]);
    }
    return coord;
}

void StoreFitTrajAllAtm(s_fit& fit, 
                        const int& frames, 
                        s_dcd& dcd, 
                        const vector<vector<int>>& psfinfo) {
    //
    for (int iframe = 0; iframe < frames; ++iframe) {
        SetupFit(psfinfo, dcd.coord[iframe], fit);
        GetTrrot(fit);
        OperateTrrot(fit, dcd.coord[iframe]);
    }
}


void GetFitAllAtm(const string& dcdfile, 
                  const int& sumatm, 
                  const vector<vector<int>>& psfinfo_fit, 
                  s_fit& fit, 
                  s_dcd& dcd) {
    //
    // Gain info. of traj.
    DCD_R dcdf(dcdfile.c_str());
    //
    dcdf.read_header();
    cout << "Read "\
         << dcdfile << endl;
    int natm   = sumatm;
    int frames = dcdf.getNFILE();
    //
    AllocateDCD(dcd, natm, frames);
    readAndStoreCoordBox(dcdf, dcd, natm, frames);
    StoreFitTrajAllAtm(fit, frames, dcd, psfinfo_fit);
}

vector<vector<vector<double>>>
GetFitcoordall(const vector<string>& dcdfile,
               const int& sumatm,
               const vector<vector<int>>& psfinfo_fit,
               const vector<vector<int>>& psfinfo_get,
               s_fit& fit) {
    //
    int sum = 0;
    int nget = psfinfo_get.size();
    for (int i = 0; i < nget; ++i) {
        int natm_get = psfinfo_get[i].size();
        for (int j = 0; j < natm_get; ++j) {
            sum += 1;
        }
    }
    //
    int nfile = dcdfile.size();
    int totfrms = 0;
    //
    for (int ifile = 0; ifile < nfile; ++ifile) {
        DCD_R dcdf(dcdfile[ifile].c_str());
        dcdf.read_header();
        int frames = dcdf.getNFILE();
        totfrms += frames;
    }
    cout << endl;
    cout << "TOTAL FRAMES: " << totfrms << endl;
    cout << endl;

    vector<vector<vector<double>>>
    coordall(totfrms, vector<vector<double>>(sum, vector<double>(3, 0.0)));
    //
    totfrms = 0;
    for (int ifile = 0; ifile < nfile; ++ifile) {
        //
        // Gain info. of traj.
        DCD_R dcdf(dcdfile[ifile].c_str());
        //
        dcdf.read_header();
        cout << "Read "
             << dcdfile[ifile] << endl;
        int natm   = sumatm;
        int frames = dcdf.getNFILE();
        //
        s_dcd dcd;
        AllocateDCD(dcd, natm, frames);
        readAndStoreCoordBox(dcdf, dcd, natm, frames);
        //
        vector<vector<vector<double>>> 
        coord = StoreFitTraj(fit, frames, dcd, psfinfo_fit, psfinfo_get);
        //
        for (int j = 0; j < frames; ++j) {
            for (int l = 0; l < sum; ++l) {
                for (int k = 0; k < 3; ++k) {
                    coordall[totfrms][l][k] = coord[j][l][k];
                }
            }
            totfrms += 1;
        }
    }
    //
    return coordall;
}

vector<vector<vector<double>>> 
Extractcoordall(const vector<string>& dcdfile, 
                const int& sumatm, 
                const vector<vector<int>>& molinfo) {
    //
    vector<vector<vector<double>>> coordall;

    for (string fdcd : dcdfile) {
        //
        // Gain info. of traj.
        DCD_R dcdf(fdcd.c_str());
        //
        dcdf.read_header();
        cout << "Read "
             << fdcd << endl;
        int natm   = sumatm;
        int frames = dcdf.getNFILE();
        //
        s_dcd dcd;
        AllocateDCD(dcd, natm, frames);
        readAndStoreCoordBox(dcdf, dcd, natm, frames);
        //
        for (int iframe = 0; iframe < frames; ++iframe) {
            //
            vector<vector<double>> cell;
            for (vector<int> cellinfo : molinfo) {
                for (int iatm : cellinfo) {
                    double x = dcd.coord[iframe][iatm-1][0];
                    double y = dcd.coord[iframe][iatm-1][1];
                    double z = dcd.coord[iframe][iatm-1][2];
                    cell.push_back({x, y, z});
                }
            }
            //
            coordall.push_back(cell);
        }
    }

    return coordall;
}

vector<vector<int>> GenerateNewID(vector<vector<int>>& molinfo) {
    //
    vector<vector<int>> molinfo_new;
    vector<int> cellinfo;

    int l = 1;
    int nres = molinfo.size();
    for (int i = 0; i < nres; ++i) {
        int natm = molinfo[i].size();
        for (int j = 0; j < natm; ++j) {
            cellinfo.push_back(l);
            l += 1;
        }
        molinfo_new.push_back(cellinfo);
        cellinfo.clear();
    }
    return molinfo_new;
}

vector<vector<vector<double>>>
ExtractAtmnmCoord(const vector<vector<vector<double>>>& coord,
                  const s_branch& branch, 
                  const string& atmnm) {
    int frames = coord.size();
    vector<vector<vector<double>>> coord_ex;
    //
    for (int i = 0; i < frames; ++i) {
        vector<vector<double>> coord_cell;
        int nres = branch.atmnm.size();
        int l = 0;
        for (int j = 0; j < nres; ++j) {
            int resatm = branch.atmnm[j].size();
            for (int k = 0; k < resatm; ++k) {
                string atmtrim = branch.atmnm[j][k];
                atmtrim.erase(remove_if(atmtrim.begin(), atmtrim.end(), ::isspace),
                              atmtrim.end());
                if (atmtrim == atmnm) {
                    double x = coord[i][l][0];
                    double y = coord[i][l][1];
                    double z = coord[i][l][2];
                    coord_cell.push_back({x, y, z});
                }
                l += 1;
            }
        }
        coord_ex.push_back(coord_cell);
        for (auto& row : coord_cell) {
            row.clear();
        }
        coord_cell.clear();
    }
    cout << endl;
    cout << "Extract number of atmnm : " << coord_ex[0].size() << endl;
    return coord_ex;
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////// combine the above module /////////////////////////////
///////////////////////////////////////////////////////////////////////////////
pair<vector<vector<vector<double>>>, s_branch>
OperateGeneralFit(const string& reffile, 
                  const string& psffile, 
                  const vector<string>& dcdfile, 
                  string& segnm, 
                  const string& fitatmnm) {
    PSFOPR PsfOpr;
    PDBOPR PdbOpr;
    //
    cout << endl;
    cout << "STEP1: Read reffile (file format is pdbfile)" << endl;
    s_pdb refpdb = PdbOpr.LoadfrompdbNeo(reffile);
    //
    cout << endl;
    cout << "STEP2: Read psffile, and gain info. of fitting atoms..." << endl;
    s_psf psf = PsfOpr.loadFile(psffile);
    int natm = psf.segname.size();
    //
    // Gain molinfo of protein.
    cout << endl;
    cout << "Select fitting mol..." << endl;
    //
    vector<vector<int>> 
    refinfo_fit = PdbOpr.getmolinfo_atmnm_all(refpdb, fitatmnm);
    vector<vector<int>> 
    psfinfo_fit = PsfOpr.getmolinfo_segatm_all(psf, segnm, fitatmnm);
    //
    cout << endl;
    cout << "STEP3: fitting coordinate" << endl;
    s_fit fit;
    RefSetupFit(psf, psfinfo_fit, refpdb.coord, refinfo_fit, fit);
    
    //if (size_t pos = segnm.find("not H"); pos != string::npos) {
    //    segnm.erase(pos, 5);
        //psfinfo_get = PsfOpr.getmolinfo_segnm_heavy(psf, segnm);
    s_branch branch = PsfOpr.getmolinfo_segnm_heavy(psf, segnm);
    //} else {
    //    psfinfo_get = PsfOpr.getmolinfo_seg_all(psf, segnm);
    //}
    vector<vector<vector<double>>>
    coordall = GetFitcoordall(dcdfile, natm, psfinfo_fit, branch.molinfo, fit);
    //
    return {coordall, branch};
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////// analyze thickness ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
vector<double> CalcThickness(const vector<vector<vector<double>>>& coord) {
    int frames = coord.size();
    int natm   = coord[0].size();
    vector<double> thickness;
    
    for (int i = 0; i < frames; ++i) {
        double length;
        double len_upp = 0.0;
        double len_low = 0.0;
        int    n_upp = 0;
        int    n_low = 0;
        for (int j = 0; j < natm; ++j) {
            if (coord[i][j][2] > 0) {
                len_upp += coord[i][j][2];
                n_upp   += 1;
            } else if (coord[i][j][2] < 0) {
                len_low += coord[i][j][2];
                n_low   += 1;
            }
        }
        len_upp /= n_upp;
        len_low /= n_low;
        length = len_upp - len_low;
        thickness.push_back(length);
    }
    return thickness;
}
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
/////////////////// variance covariance matrix (length) /////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//vector<vector<double>> 
//CalcVarCovarMat(const vector<vector<vector<double>>>& coord, 
//                const ) {
//    vector<vector<double>> matrix;
//
//    return matrix;
//}
/////////////////////////////////////////////////////////////////////////////////
