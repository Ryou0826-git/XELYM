// make_hist.hpp
#ifndef MAKE_HIST_HPP
#define MAKE_HIST_HPP

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

using namespace std;

enum class NormType {
    COUNT,
    PROB,
    PDF
};

struct Histogram {
    vector<double> bin_center;
    vector<double> value;
    double bin_width;
};

vector<double> Create1dbin(const double& xmin, 
                           const double& xmax, 
                           const double& delx);

Histogram make_histgram(const vector<double>& data,
                        int   nbins,
                        NormType norm = NormType::COUNT);

vector<vector<double>> make_3d_grid(const vector<double>& Min,
                                    const vector<double>& Max,
                                    int grid);

//--> iline function
inline int grid_index_1d(double x,
                         double xmin,
                         double xmax,
                         int grid) {
    double t = (x - xmin) / (xmax - xmin);
    int i = static_cast<int>(floor(t * grid));
    //
    if (i < 0) i = 0;
    if (i >= grid) i = grid - 1;
    //
    return i;
}

inline vector<int> grid_index_3d(const vector<double>& r,
                                 const vector<double>& Min,
                                 const vector<double>& Max,
                                 int grid) {
    //
    int ix = grid_index_1d(r[0], Min[0], Max[0], grid);
    int iy = grid_index_1d(r[1], Min[1], Max[1], grid);
    int iz = grid_index_1d(r[2], Min[2], Max[2], grid);
    //
    return {ix, iy, iz};
}

inline int idx(vector<int> index_vec, int grid) {
    //
    int ix = index_vec[0];
    int iy = index_vec[1];
    int iz = index_vec[2];
    //
    return ix + grid * (iy + grid * iz);
}
//<--

#endif
