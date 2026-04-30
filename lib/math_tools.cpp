//math_tools.cpp

#include "math_tools.hpp"

using namespace std;

vector<double> crossProduct(const vector<double>& a, 
                            const vector<double>& b) {
    if (a.size() != 3 || b.size() != 3) {
        throw invalid_argument("Both vectors must be of size 3.");
    }

    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

double dotProduct(const vector<double>& a, 
                  const vector<double>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double norm(const vector<double>& v) { return sqrt(dotProduct(v, v)); }

void PerformRodrigues(vector<vector<double>>& coord, 
                      const vector<double>& nvec, 
                      const double& theta) {
    //
    int num = coord.size();
    vector<vector<double>> 
    Rotcoord(num, vector<double>(3, 0.0f));
    //
    for (int i = 0; i < num; ++i) {
        vector<double> rvec(3, 0.0f);
        for (int j = 0; j < 3; ++j) {
            rvec[j] = coord[i][j];
        }
        double dot = inner_product(nvec.begin(), nvec.end(), rvec.begin(), 0.0f);
        vector<double> crossvec = crossProduct(nvec, rvec);
        //
        for (int k = 0; k < 3; ++k) {
            Rotcoord[i][k] = rvec[k] * cos(theta) 
                             + nvec[k] * dot * (1 - cos(theta)) 
                             + crossvec[k] * sin(theta);
        }
    }
    //
    for (int i = 0; i < num; ++i) {
        for (int j = 0; j < 3; ++j) {
            coord[i][j] = Rotcoord[i][j];
        }
    }
    return;
}

vector<double> ReturnRandomvec(const int& n) {
     vector<double> nvec(n);
     //
     random_device rd;
     mt19937 gen(rd());
     uniform_real_distribution<> dis(-1.0, 1.0);
     //
     for (int i = 0; i < n; ++i) {
         nvec[i] = dis(gen);
     }
     double absn = 0.0f;
     for (int i = 0; i < n; ++i) {
         absn += nvec[i] * nvec[i];
     }
     absn = sqrt(absn);
     for (int i = 0; i < n; ++i) {
         nvec[i] /= absn;
     }
     //
     return nvec;
}

double ReturnRandomtheta() {
    mt19937 gen(random_device{}());
    uniform_real_distribution<> dis(0.0, M_PI);

    return dis(gen);
}

double ReturnRandomphi() {
    mt19937 gen(random_device{}());
    uniform_real_distribution<> dis(0.0, 2*M_PI);

    return dis(gen);
}

pair<mt19937, uniform_real_distribution<>> 
InitRandUni(const double& zmin, 
            const double& zmax) {
    //
    mt19937 gen(random_device{}());
    uniform_real_distribution<> dist(zmin, zmax);
    //
    return {gen, dist};
}

pair<mt19937, normal_distribution<>> 
InitRand(const double& mean, const double& stddev) {
    //
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dist(mean, stddev);

    return {gen, dist};
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////// For eigen  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
pair<Eigen::VectorXcd, Eigen::MatrixXcd> PerformDiag(const Eigen::MatrixXd& Mat) {
    //
    Eigen::EigenSolver<Eigen::MatrixXd> solver(Mat);
    Eigen::VectorXcd eigenvalues  = solver.eigenvalues();
    Eigen::MatrixXcd eigenvectors = solver.eigenvectors();

    vector<int> indices(eigenvalues.size());

    for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;

    sort(indices.begin(), indices.end(), [&](int i, int j) {
        return eigenvalues[i].real() > eigenvalues[j].real();
    });

    Eigen::VectorXcd sorted_eigenvalues(eigenvalues.size());
    Eigen::MatrixXcd sorted_eigenvectors(eigenvectors.rows(), eigenvectors.cols());

    for (size_t i = 0; i < indices.size(); ++i) {
        sorted_eigenvalues[i]      = eigenvalues[indices[i]];
        sorted_eigenvectors.col(i) = eigenvectors.col(indices[i]);
    }

    return {sorted_eigenvalues, sorted_eigenvectors};
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////// For covariance : length //////////////////////////////
///////////////////////////////////////////////////////////////////////////////
vector<vector<double>> NormVarCovL(const vector<vector<vector<double>>>& coord) {

    int frames = coord.size();
    int natm   = coord[0].size();
    //
    Eigen::MatrixXd VarCovMat(natm, natm);
    VarCovMat.setZero();
    //
    // Calc. average coordinate...
    vector<vector<double>> avecoord(natm, vector<double>(3, 0.0));

    //
    // Average coordinate:
    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < natm; ++j) {
            for (int k = 0; k < 3; ++k) {
                avecoord[j][k] += coord[i][j][k];
            }
        }
    }
    //
    for (int j = 0; j < natm; ++j) {
        for (int k = 0; k < 3; ++k) {
            avecoord[j][k] /= frames;
        }
    }

    //
    // Calc. variance-convariance matrix
    vector<vector<double>> delr_jl(natm, vector<double>(natm, 0.0));
    vector<double> delr2_ave(natm, 0.0);

    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < natm; ++j) {
            vector<double> delr_j(3, 0.0);
            for (int k = 0; k < 3; ++k) {
                delr_j[k] = coord[i][j][k] - avecoord[j][k];
            }
            //
            delr2_ave[j] += delr_j[0]*delr_j[0] + delr_j[1]*delr_j[1] + delr_j[2]*delr_j[2];

            for (int l = 0; l < natm; ++l) {
                vector<double> delr_l(3, 0.0);
                for (int m = 0; m < 3; ++m) {
                    delr_l[m] = coord[i][l][m] - avecoord[l][m];
                }
                delr_jl[j][l] += delr_j[0]*delr_l[0] + delr_j[1]*delr_l[1] + delr_j[2]*delr_l[2];
            }
        }
    }

    for (int i = 0; i < natm; ++i) {
        delr2_ave[i] /= frames;
    }

    for (int j = 0; j < natm; ++j) {
        for (int l = 0; l < natm; ++l) {
            double norm = delr2_ave[j] * delr2_ave[l];
            norm = sqrt(norm);
            delr_jl[j][l] /= frames;
            delr_jl[j][l] /= norm;
        }
    }

    return delr_jl;
}
///////////////////////////////////////////////////////////////////////////////

