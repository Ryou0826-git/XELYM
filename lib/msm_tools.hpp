// msm_tools.hpp
#ifndef MSM_TOOLS_HPP
#define MSM_TOOLS_HPP

#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>

#include "math_tools.hpp"
#include <utility>
#include <boost/algorithm/string.hpp>
#include <complex>

using namespace std;

struct s_msm {
    vector<double> pop;
    vector<double> fe;
    Eigen::MatrixXd TrMat;
};

Eigen::MatrixXd GainTrMat(const vector<int>& data,
                          const int& tau);

s_msm PerformMSM(const Eigen::MatrixXd& TrMat_norm,
                 const double& kT);

void OutMSMresult(const s_msm& msm,
                  const string& filename);

#endif
