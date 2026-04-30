// math_tools.hpp
#ifndef MATH_TOOLS_HPP
#define MATH_TOOLS_HPP

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <numeric>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <random>

using namespace std;

vector<double> crossProduct(const vector<double>& a, 
                            const vector<double>& b);

double dotProduct(const vector<double>& a, 
                  const vector<double>& b);

double norm(const vector<double>& v);

void PerformRodrigues(vector<vector<double>>& coord, 
                      const vector<double>& nvec, 
                      const double& theta);

vector<double> ReturnRandomvec(const int& n);

double ReturnRandomtheta();

double ReturnRandomphi();

pair<mt19937, normal_distribution<>> 
InitRand(const double& mean, const double& stddev);

pair<mt19937, uniform_real_distribution<>>
InitRandUni(const double& zmin,
            const double& zmax);

pair<Eigen::VectorXcd, Eigen::MatrixXcd> PerformDiag(const Eigen::MatrixXd& Mat);

vector<vector<double>> NormVarCovL(const vector<vector<vector<double>>>& coord);

///////////////////////////////////////////////////////////////////////////////
class RandomGenerator {
public:
    RandomGenerator() : gen(rd()), dist_minus1_1(-1.0, 1.0), dist_0_2pi(0.0, 2.0 * M_PI) {}

    vector<double> RandomUnitVector(int n) {
        vector<double> vec(n);
        for (int i = 0; i < n; ++i) {
            vec[i] = dist_minus1_1(gen);
        }
        double norm = sqrt(inner_product(vec.begin(), vec.end(), vec.begin(), 0.0));
        for (auto& v : vec) v /= norm;
        return vec;
    }

    vector<double> RandomVector(int n, double len) {
        vector<double> vec(n);
        for (int i = 0; i < n; ++i) {
            vec[i] = dist_minus1_1(gen);
        }
        double norm = sqrt(inner_product(vec.begin(), vec.end(), vec.begin(), 0.0));
        norm /= len;
        for (auto& v : vec) v /= norm;
        return vec ;
    }

    double RandomAngle() {
        return dist_0_2pi(gen);
    }

private:
    random_device rd;
    mt19937 gen;
    uniform_real_distribution<> dist_minus1_1;
    uniform_real_distribution<> dist_0_2pi;
};
///////////////////////////////////////////////////////////////////////////////

#endif
