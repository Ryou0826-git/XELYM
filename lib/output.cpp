// output.cpp
#include "output.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

void CoordToFile(const string& filename, 
                 const vector<vector<vector<double>>>& data, 
                 int selid) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    int frames = data.size();

    for (int i = 0; i < frames; ++i) {
        outfile << setw(10) << i+1 
                << "   "
                << setw(15) << setprecision(10) << data[i][selid-1][0] 
                << "   "
                << setw(15) << setprecision(10) << data[i][selid-1][1] 
                << "   "
                << setw(15) << setprecision(10) << data[i][selid-1][2] << endl;
    }

    outfile.close();
}

void data3dToFile(const string& filename, 
                  const vector<vector<double>>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    int frames = data[0].size();

    for (int i = 0; i < frames; ++i) {
        outfile << setw(15) << setprecision(10) << data[0][i] << "   "
                << setw(15) << setprecision(10) << data[1][i] << "   "
                << setw(15) << setprecision(10) << data[2][i] << endl;
    }
    //
    outfile.close();
}

void data2dToFile(const string& filename, 
                  const vector<vector<double>>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int frames = data[0].size();
    //
    for (int i = 0; i < frames; ++i) {
        outfile << setw(15) << setprecision(10) << data[0][i] 
                << "   "
                << setw(15) << setprecision(10) << data[1][i] << endl;
    }
    outfile.close();
}

void data2dToFileLines(const string& filename, 
                       const vector<vector<double>>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int frames = data[0].size();
    //
    for (int i = 0; i < frames; ++i) {
        outfile << setw(10) << i + 1 
                << "   "
                << setw(15) << setprecision(10) << data[0][i] 
                << "   "
                << setw(15) << setprecision(10) << data[1][i] << endl;
    }
    outfile.close();
}

void AxisToFile(const string& filename, 
                const vector<double>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int frames = data.size();

    for (int i = 0; i < frames; ++i) {
        outfile << setw(20) << setprecision(10) << data[i] << flush;
    }
    outfile.close();
}

void Only1dToFile(const string& filename, 
                  const vector<double>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int frames = data.size();

    for (int i = 0; i < frames; ++i) {
        outfile << setw(20) << setprecision(10) << data[i] << endl;;
    }
    outfile.close();
}

void Mapdata2dToFile(const string& filename, 
                     const vector<vector<double>>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int binx = data.size();
    int biny = data[0].size();
    //
    for (int i = 0; i < biny; ++i) {
        for (int j = 0; j < binx; ++j) {
            outfile << setw(20) << setprecision(10) << data[j][i] 
                    << flush;
        }
        outfile << endl;
    }
    outfile.close();
}

void SplitdataToFile(const string& filename, 
                     const vector<double>& datax, 
                     const vector<double>& datay) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int frames = datax.size();
    //
    for (int i = 0; i < frames; ++i) {
        outfile << setw(10) << setprecision(5) << datax[i] 
                << "   "
                << setw(15) << setprecision(10) << datay[i] << endl;
    }
    outfile.close();
}

void SplitdataToFileINTdouble(const string& filename, 
                             const vector<double>& datax, 
                             const vector<double>& datay) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int frames = datax.size();
    //
    for (int i = 0; i < frames; ++i) {
        outfile << setw(10) << datax[i] 
                << "   "
                << setw(15) << setprecision(10) << datay[i] << endl;
    }
    outfile.close();
}

void SplitdataToFileNUMdouble(const string& filename, 
                              const vector<double>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int frames = data.size();
    //
    for (int i = 0; i < frames; ++i) {
        outfile << setw(10) << i + 1 
                << "   "
                << setw(15) << setprecision(10) << data[i] << endl;
    }
    outfile.close();
}

void SplitINTdataToFileNUMdouble(const string& filename, 
                                 const vector<int>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int frames = data.size();
    //
    for (int i = 0; i < frames; ++i) {
        outfile << setw(10) << i + 1 
                << "   "
                << data[i] << endl;
    }
    outfile.close();
}

void OutdxFile(const vector<int>& ng3, 
               const vector<double>& del3, 
               const vector<double>& origin, 
               const vector<vector<vector<double>>>& data, 
               const string& filename) {
    //
    double zero = 0.0;
    int sgrid = ng3[0] * ng3[1] * ng3[2];
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    outfile << scientific
            << uppercase
            << setprecision(6);

    outfile << "object 1 class gridpositions counts " << ng3[0] << "  " << ng3[1] << "  " << ng3[2] << endl;
    outfile << "origin " << origin[0] << "  " << origin[1] << "  " << origin[2] << endl;
    outfile << "delta " << del3[0] << "  " <<  zero    << "  " << zero    << endl;
    outfile << "delta " << zero    << "  " <<  del3[1] << "  " << zero    << endl;
    outfile << "delta " << zero    << "  " <<  zero    << "  " << del3[2] << endl;
    outfile << "object 2 class gridpositions counts " << ng3[0] << "  " <<  ng3[1] << "  " << ng3[2] << endl;
    outfile << "object 3 class array type double rank 0 items " << sgrid << " data follows" << endl;
    
    int count = 0;
    for (int i = 0; i < ng3[2]; ++i) {
        for (int j = 0; j < ng3[1]; ++j) {
            for (int k = 0; k < ng3[0]; ++k) {
                outfile << data[k][j][i] << " ";
                count++;
                if (count % 3 == 0) {
                    outfile << "\n";
                }
            }
        }
    }

    outfile.close();
}

void dataXdToFile(const string& filename, 
                  const vector<vector<double>>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int ncomp  = data.size();
    int frames = data[0].size();
    //
    for (int i = 0; i < frames; ++i) {
        outfile << setw(10) << i+1  << "   ";
        for (int j = 0; j < ncomp; ++j) {      
            outfile << setw(15) << setprecision(10) << data[j][i] << "   ";
        }
        outfile << endl; 
    }
    outfile.close();
}

void INTdataXdToFile(const string& filename,
                     const vector<int>& times, 
                     const vector<vector<double>>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int tot = times.size();
    int cvdim = data[0].size();
    //
    for (int i = 0; i < tot; ++i) {
        outfile << setw(10) << times[i]  << "   ";
        for (int j = 0; j < cvdim; ++j) {
            outfile << setw(15) << setprecision(10) << data[i][j] << "   ";
        }
        outfile << endl;
    }
    outfile.close();
}

void TimedataXdToFile(const string& filename,
                     const vector<double>& times, 
                     const vector<vector<double>>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int tot = times.size();
    int cvdim = data[0].size();
    //
    for (int i = 0; i < tot; ++i) {
        outfile << setw(15) << setprecision(10) << times[i]  << "   ";
        for (int j = 0; j < cvdim; ++j) {
            outfile << setw(15) << setprecision(10) << data[i][j] << "   ";
        }
        outfile << endl;
    }
    outfile.close();
}

void TimeINTdataXdToFile(const string& filename,
                         const vector<double>& times, 
                         const vector<vector<int>>& data) {
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int tot = times.size();
    int cvdim = data[0].size();
    //
    for (int i = 0; i < tot; ++i) {
        outfile << setw(15) << setprecision(10) << times[i]  << "   ";
        for (int j = 0; j < cvdim; ++j) {
            outfile << setw(15) << setprecision(10) << data[i][j] << "   ";
        }
        outfile << endl;
    }
    outfile.close();
}

void EigenMapToFile(const string& filename, 
                    const Eigen::MatrixXd data) {
    ofstream file(filename);
    if (file.is_open()) {
        for (int i = 0; i < data.rows(); ++i) {
            for (int j = 0; j < data.cols(); ++j) {
                file << setw(15) << setprecision(10) << data(i, j) << "   ";
                //file << data(i, j);
            }
            file << "\n";
        }
        file.close();
    } else {
        cout << "Error opennig file; " << filename << endl;
    }
}

//
//  --> Output 3d info. to file
void write_3d_hist_gnuplot(const string& filename,
                           const vector<double>& hist,
                           const vector<double>& Min,
                           const vector<double>& Max,
                           int grid) {
    //
    ofstream ofs(filename);
    //
    double dx[3];
    for (int d = 0; d < 3; ++d)
        dx[d] = (Max[d] - Min[d]) / grid;

    for (int iz = 0; iz < grid; ++iz) {
        for (int iy = 0; iy < grid; ++iy) {
            for (int ix = 0; ix < grid; ++ix) {
                double x = Min[0] + (ix + 0.5) * dx[0];
                double y = Min[1] + (iy + 0.5) * dx[1];
                double z = Min[2] + (iz + 0.5) * dx[2];

                ofs << x << " "
                    << y << " "
                    << z << " "
                    << hist[idx({ix, iy, iz}, grid)] 
                    << "\n";
            }
        }
    }
}

