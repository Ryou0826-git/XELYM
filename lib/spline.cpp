// spline.cpp

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include <spline.hpp>

using namespace std;

double Func3rd(const double& a, const double& b, const double& c, const double& d,\
               const double& xi, const double& x) {
    double delx = x - xi;
    double y = a + b * delx + c * delx * delx + d * delx * delx * delx;
    return y;
}

s_coeff Spline3rd(const vector<double>& x, 
                  const vector<double>& y) {
    int dimx = x.size();
    int dimy = y.size();
    s_coeff coeff;

    if (dimy != dimx) {
        cout << "Error: dim of x != dim of y" << endl;
    }
    coeff.a.resize(dimx, 0.0);
    coeff.b.resize(dimx-1, 0.0);
    coeff.c.resize(dimx, 0.0);
    coeff.d.resize(dimx-1, 0.0);

    for (int i = 0; i < dimx; ++i) {
        coeff.a[i] = y[i];
    }

    vector<double> h(dimx-1, 0);
    //
    for (int i = 0; i < dimx-1; ++i) {
        h[i] = x[i+1] - x[i];
    }

    Eigen::MatrixXd Amat = Eigen::MatrixXd::Zero(dimx, dimx);

    //
    // Create matrix of A
    for (int i = 0; i < dimx; ++i) {
        for (int j = 0; j < dimx; ++j) {
            if (i == j) {
                if ((i == 0) || (i == dimx-1)) {
                    Amat(i, j) = 1.0f;
                } else {
                    Amat(i, j)   = 2.0f * (h[i-1] + h[i]);
                    Amat(i, j-1) = h[i-1];
                    Amat(i, j+1) = h[i];
                }
            }
        }
    }

    //
    // Calc. inv. of matrix of A(i, j)
    Eigen::MatrixXd Ainv = Amat.inverse();

    vector<double> bmat(vector<double>(dimx, 0));
    for (int i = 0; i < dimx; ++i) {
        if ((i == 0) || (i == dimx-1)) {
            bmat[i] = 0.0;
        } else {
            bmat[i] = 3.0f / h[i] * (coeff.a[i+1] - coeff.a[i]) 
                      - 3.0f / h[i-1] * (coeff.a[i] - coeff.a[i-1]);
        }
    }

    for (int i = 0; i < dimx; ++i) {
        for (int j = 0; j < dimx; ++j) {
            coeff.c[i] += Ainv(i, j) * bmat[j];
        }
    }

    for (int i = 0; i < dimx-1; ++i){
        coeff.b[i] = (coeff.a[i+1] - coeff.a[i]) / h[i] 
                     - h[i] * (coeff.c[i+1] + 2.0f * coeff.c[i]) / 3.0f;
        coeff.d[i] = (coeff.c[i+1] - coeff.c[i]) / 3.0f / h[i];
    }

    return coeff;
}

vector<double> 
PerformSpline3rd(const s_coeff& coeff, 
                 const vector<double>& x, 
                 const vector<double>& xnew) {
  //
  int bin    = x.size();
  int binnew = xnew.size();
  vector<double> ynew(binnew, 0.0);
  //
  for (int i = 0; i < binnew; ++i) {
      for (int j = 0; j < bin-1; ++j) {
          if ((x[j] <= xnew[i]) && (xnew[i] < x[j+1])) {
              ynew[i] = Func3rd(coeff.a[j], coeff.b[j], coeff.c[j], coeff.d[j], 
                                x[j], xnew[i]);
              break;
          }
      }
  }

  return ynew;
}

vector<s_coeff> 
Createcoefflist(const vector<double>& datax, 
                const vector<vector<double>>& dataz){
    //
    int dimx = dataz[0].size();
    int dimy = dataz.size();
    vector<double> splitz(dimx, 0.0);
    vector<s_coeff> coefflist(dimy);
    //
    for (int i = 0; i < dimy; ++i) {
        for (int j = 0; j < dimx; ++j) {
            splitz[j] = dataz[i][j];
        }
        coefflist[i] = Spline3rd(datax, splitz);
    }

    return coefflist;
}

s_coeff Spline3rd2d(const vector<s_coeff> coefflist, 
                    const vector<double> datax, 
                    const vector<double> datay, 
                    const double& x) {
    //
    int nx   = datax.size();
    int dimy = coefflist.size();
    vector<double> zy(dimy, 0.0);
    s_coeff coeff;
    //

    for (int i = 0; i < nx; ++i) {
        if ((datax[i] <= x) && (x < datax[i+1])) {
            for (int j = 0; j < dimy; ++j) {
                zy[j] = Func3rd(coefflist[j].a[i], 
                                coefflist[j].b[i], 
                                coefflist[j].c[i], 
                                coefflist[j].d[i], 
                                datax[i], x);
            }
            coeff = Spline3rd(datay, zy);
            break;
        }
    }
    return coeff;
}

vector<s_coeff> 
GainCoeffxSpline3rd2d(const vector<double>& x, 
                      const vector<double>& y, 
                      const vector<vector<double>>& z, 
                      const vector<double>& xnew, 
                      const vector<double>& ynew) {
    //
    int binx = xnew.size();
    int biny = ynew.size();
    //
    vector<vector<double>> 
    znew(binx, vector<double>(biny, 0.0));

    vector<s_coeff> 
    coefflistx = Createcoefflist(x, z);

    vector<s_coeff> coefflist(binx);

    for (int i = 0; i < binx; ++i) {
        coefflist[i] = Spline3rd2d(coefflistx, x, y, xnew[i]);
    }

    return coefflist;
}

vector<vector<double>> 
PerformSpline3rd2d(const vector<double>&         x, 
                   const vector<double>&         y, 
                   const vector<vector<double>>& z, 
                   const vector<double>&         xnew, 
                   const vector<double>&         ynew) {
    //
    vector<s_coeff> 
    coefflist = GainCoeffxSpline3rd2d(x, y, z, xnew, ynew);
    //
    int ny   = y.size();
    int binx = xnew.size();
    int biny = ynew.size();
    //
    vector<vector<double>> 
    znew(binx, vector<double>(biny, 0.0));

    for (int i = 0; i < binx; ++i) {
        for (int j = 0; j < biny; ++j) {
            for (int k = 0; k < ny; ++ k) {
                if (k != ny-1) {
                    if ((y[k] <= ynew[j]) && (ynew[j] < y[k+1])) {
                        znew[i][j] = Func3rd(coefflist[i].a[k], 
                                             coefflist[i].b[k], 
                                             coefflist[i].c[k], 
                                             coefflist[i].d[k], 
                                             y[k], ynew[j]);
                        break;
                    }
                } else {
                    znew[i][j] = Func3rd(coefflist[i].a[k], 
                                         coefflist[i].b[k], 
                                         coefflist[i].c[k], 
                                         coefflist[i].d[k], 
                                         y[k], ynew[j]);
                }
            }
        }
    }

    return znew;
}
