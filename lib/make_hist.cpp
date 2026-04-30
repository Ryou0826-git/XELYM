//make_hist.cpp
#include <make_hist.hpp>

#include <vector>
#include <iostream>

using namespace std;

vector<double> Create1dbin(const double& xmin, 
                           const double& xmax, 
                           const double& delx) {
    //
    int binx = static_cast<int>((xmax - xmin) / delx);

    //
    // Allocate data of hist
    vector<double> hist(binx, 0.0);
    //
    for (int i = 0; i < binx; ++i) {
        hist[i] = xmin + (i + 0.5f) * delx;
    }
    //
    return hist;
}

//
// {{x, y, z}, {x, y, z}, ...}
vector<vector<double>> make_3d_grid(const vector<double>& Min,
                                    const vector<double>& Max,
                                    int grid) {
    if (Min.size() != 3 || Max.size() != 3)
        throw invalid_argument("Min and Max must be size 3");
    if (grid <= 0)
        throw invalid_argument("grid must be positive");

    vector<vector<double>> points;
    points.reserve(grid * grid * grid);

    double dx[3];
    for (int d = 0; d < 3; ++d)
        dx[d] = (Max[d] - Min[d]) / grid;

    for (int i = 0; i < grid; ++i) {
        for (int j = 0; j < grid; ++j) {
            for (int k = 0; k < grid; ++k) {
                vector<double> r(3);
                r[0] = Min[0] + (i + 0.5) * dx[0];
                r[1] = Min[1] + (j + 0.5) * dx[1];
                r[2] = Min[2] + (k + 0.5) * dx[2];
                points.push_back(r);
            }
        }
    }
    //
    return points;
}

Histogram make_histgram(const vector<double>& data,
                        int nbins,
                        NormType norm) {
    if (data.empty()) {
        throw runtime_error("data is empty");
    }

    auto [min_it, max_it] = minmax_element(data.begin(), data.end());
    double xmin = *min_it;
    double xmax = *max_it;

    xmax += 1e-12;

    double bin_width = (xmax - xmin) / nbins;

    vector<int> count(nbins, 0);
    vector<double> bin_center(nbins);

    for (int i = 0; i < nbins; ++i) {
        bin_center[i] = xmin + (i + 0.5) * bin_width;
    }

    for (double x : data) {
        int ibin = static_cast<int>((x - xmin) / bin_width);
        if (ibin < nbins) {
            count[ibin]++;
        }
    }

    int N = data.size();
    vector<double> value(nbins, 0.0);

    for (int i = 0; i < nbins; ++i) {
        if (norm == NormType::COUNT) {
            value[i] = static_cast<double>(count[i]);
        }
        else if (norm == NormType::PROB) {
            value[i] = static_cast<double>(count[i]) / N;
        }
        else if (norm == NormType::PDF) {
            value[i] = static_cast<double>(count[i]) / (N * bin_width);
        }
    }
    //
    return Histogram{bin_center, value, bin_width};
}

