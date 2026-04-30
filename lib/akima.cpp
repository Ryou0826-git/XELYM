// akima.cpp

#include <vector>
#include <cmath>
#include <iostream>

#include <akima.hpp>
#include <spline.hpp>

using namespace std;

double Func3rdakima(const double& a, const double& b, const double& c, const double& d,\
                    const double& xi, const double& x) {
    double delx = x - xi;
    double y = a + b * delx + c * delx * delx + d * delx * delx * delx;
    return y;
}

void AkimaRegion1(const double& x3, 
                  const double& y3, 
                  const double& x4, 
                  const double& y4, 
                  const double& x5, 
                  const double& y5, 
                  double& x1, double& y1, 
                  double& x2, double& y2) {
    x1 = -x5 + 2.0f * x3;
    x2 = x3 + x4 -x5;
    //
    double m3 = (y4 - y3) / (x4 - x3);
    double m4 = (y5 - y4) / (x5 - x4);
    y2 = y3 - (x3 - x2) / (-m4 + 2.0f*m3);
    y1 = y2 - (x2 - x1) / (-2.0f*m4 + 3.0f*m3);
    return;
}

void AkimaRegion3(const double& x1, 
                  const double& y1, 
                  const double& x2, 
                  const double& y2, 
                  const double& x3, 
                  const double& y3, 
                  double& x4, double& y4, 
                  double& x5, double& y5) {
    x4 = - x1 + x2 + x3;
    x5 = - x1 + 2.0f*x3;
    //
    double m1 = (y2 - y1) / (x2 - x1);
    double m2 = (y3 - y2) / (x3 - x2);
    y4 = y3 + (x4 - x3) / (2.0f * m2 - m1);
    y5 = y4 + (x5 - x4) / (3.0f * m2 - 2.0f * m1);
    return;
}

s_coeff AkimaRegion2(const vector<double>& datax, 
                     const vector<double>& datay) {
    //
    double m0, m1, m2, m3, m4;
    double w1, w2;
    double q1, q2;
    int dimx = datax.size();
    int dimy = datay.size();
    //
    s_coeff coeff;
    coeff.a.resize(dimx-4, 0.0);
    coeff.b.resize(dimx-4, 0.0);
    coeff.c.resize(dimx-4, 0.0);
    coeff.d.resize(dimx-4, 0.0);
    //
    if (dimx != dimy) {
        cout << "Error: dimx != dimy..." << endl;
    }
    //
    for (int i = 0; i < dimx; ++i) {
        if ((i > 1)  && (i < dimx-2)) {
            m0 = (datay[i-1] - datay[i-2]) / (datax[i-1] - datax[i-2]);
            m1 = (datay[i] - datay[i-1]) / (datax[i] - datax[i-1]);
            m2 = (datay[i+1] - datay[i]) / (datax[i+1] - datax[i]);
            m3 = (datay[i+2] - datay[i+1]) / (datax[i+2] - datax[i+1]);
            m4 = (datay[i+3] - datay[i+2]) / (datax[i+3] - datax[i+2]);

            //
            w1 = abs(m1 - m0) / abs(m3 - m2);
            w2 = abs(m2 - m1) / abs(m4 - m3);
            //
            q1 = (m1 - m2) / (1.0f + w1);
            q2 = (m3 - m2) / (1.0f + w2);
            //
            coeff.a[i-2] = datay[i];
            coeff.b[i-2] = q1 + m2;
            coeff.c[i-2] = - (2.0f*q1 + q2) / (datax[i+1] - datax[i]);
            coeff.d[i-2] = (q1 + q2) / (datax[i+1] - datax[i]) / (datax[i+1] - datax[i]);
        }
    }
    return coeff;
}

vector<double> 
PerformAkima1d(const vector<double>& x, 
               const vector<double>& y, 
               const vector<double>& xnew) {
    //
    int dimx = x.size();
    int nx = xnew.size();
    s_coeff coeff;
    //
    vector<double> datax(dimx+4, 0.0);
    vector<double> datay(dimx+4, 0.0);
    //
    for (int i = 0; i < dimx; ++i) {
        datax[i+2] = x[i];
        datay[i+2] = y[i];
    }
    //
    AkimaRegion1(datax[2], datay[2], 
                 datax[3], datay[3], 
                 datax[4], datay[4], 
                 datax[0], datay[0],
                 datax[1], datay[1]);
    //
    AkimaRegion3(datax[dimx-1], datay[dimx-1], 
                 datax[dimx], datay[dimx], 
                 datax[dimx+1], datay[dimx+1], 
                 datax[dimx+2], datay[dimx+2], 
                 datax[dimx+3], datay[dimx+3]);

    coeff = AkimaRegion2(datax, datay);
    //
    vector<double> ynew(nx, 0.0);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < dimx+4; ++j){
            if ((datax[j] <= xnew[i]) && (xnew[i] < datax[j+1])) {
                ynew[i] = Func3rdakima(coeff.a[j-2], coeff.b[j-2], coeff.c[j-2], coeff.d[j-2], 
                                       datax[j], xnew[i]);
                break;
            }
        }
    }

    return ynew;
}

vector<s_coeff> 
CreateCoefflistAkima(const vector<double>& datax, 
                     const vector<vector<double>>& dataz) {
    //
    int dimx = dataz[0].size();
    int dimy = dataz.size();
    vector<double> splitz(dimx+4, 0.0);
    vector<double> xnew(dimx+4, 0.0);
    vector<s_coeff> coefflist(dimy);
    //
    for (int i = 2; i < dimx+2; ++i) {
       xnew[i] = datax[i-2];
    }
    //
    for (int i = 0; i < dimy; ++i) {
        for (int j = 2; j < dimx+2; ++j) {
            splitz[j] = dataz[i][j];
        }
        AkimaRegion1(xnew[2], splitz[2], 
                     xnew[3], splitz[3], 
                     xnew[4], splitz[4], 
                     xnew[0], splitz[0], 
                     xnew[1], splitz[1]);
        //
        AkimaRegion3(xnew[dimx-1], splitz[dimx-1], 
                     xnew[dimx],   splitz[dimx],   
                     xnew[dimx+1], splitz[dimx+1], 
                     xnew[dimx+2], splitz[dimx+2], 
                     xnew[dimx+3], splitz[dimx+3]);
        //
        coefflist[i] =  AkimaRegion2(xnew, splitz);
    }

    return coefflist;
}

s_coeff Akima2d(const vector<s_coeff>& coefflist, 
                const vector<double>& datax, 
                const vector<double>& datay, 
                const double& x) {
    //
    int nx   = datax.size();
    int dimy = coefflist.size();
    vector<double> zy(dimy+4, 0.0);
    s_coeff coeff;
    //
    vector<double> ynew(dimy+4, 0.0);
    for (int i = 2; i < dimy+2; ++i) {
        ynew[i] = datay[i-2];
    }

    for (int i = 0; i < nx; ++i) {
        if ((datax[i] <= x) && (x < datax[i+1])) {
            for (int j = 2; j < dimy+2; ++j) {
                zy[j] = Func3rdakima(coefflist[j-2].a[i], 
                                     coefflist[j-2].b[i], 
                                     coefflist[j-2].c[i], 
                                     coefflist[j-2].d[i], 
                                     datax[i], x);
                //cout << zy[j] << endl;
            }
            AkimaRegion1(ynew[2], zy[2], 
                         ynew[3], zy[3], 
                         ynew[4], zy[4], 
                         ynew[0], zy[0], 
                         ynew[1], zy[1]);
            //
            AkimaRegion3(ynew[dimy-1], zy[dimy-1], 
                         ynew[dimy],   zy[dimy],   
                         ynew[dimy+1], zy[dimy+1], 
                         ynew[dimy+2], zy[dimy+2], 
                         ynew[dimy+3], zy[dimy+3]);
            //
            coeff = AkimaRegion2(ynew, zy);
        }
    }

    return coeff;
}

vector<s_coeff> 
GainCoeffxAkima2d(const vector<double>&         x, 
                  const vector<double>&         y, 
                  const vector<vector<double>>& z, 
                  const vector<double>&         xnew, 
                  const vector<double>&         ynew) {
    //
    int binx = xnew.size();
    int biny = ynew.size();
    //
    vector<vector<double>> 
    znew(binx, vector<double>(biny, 0.0));

    vector<s_coeff> 
    coefflistx = CreateCoefflistAkima(x, z);

    vector<s_coeff> coefflist(binx);

    for (int i = 0; i < binx; ++i) {
        coefflist[i] = Akima2d(coefflistx, x, y, xnew[i]);
    }

    return coefflist;
}

vector<vector<double>> 
PerformAkima2d(const vector<double>&         x, 
               const vector<double>&         y, 
               const vector<vector<double>>& z, 
               const vector<double>&         xnew, 
               const vector<double>&         ynew) {
    //
    vector<s_coeff> 
    coefflist = GainCoeffxAkima2d(x, y, z, xnew, ynew);
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
                        znew[i][j] = Func3rdakima(coefflist[i].a[k], 
                                                  coefflist[i].b[k], 
                                                  coefflist[i].c[k], 
                                                  coefflist[i].d[k], 
                                                  y[k], ynew[j]);
                        break;
                    }
                } else {
                    znew[i][j] = Func3rdakima(coefflist[i].a[k], 
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
