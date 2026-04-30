// extract.cpp
#include <extract.hpp>

using namespace std;

vector<vector<double>> 
ExtractCoord(const vector<vector<double>>& coord, 
             const vector<vector<int>>& molinfo) {
    int nres = molinfo.size();
    vector<vector<double>> ExCoord(nres, vector<double>(3, 0.0));
    //
    int l = 0;
    for (int i = 0; i < nres; ++i) {
        int natm = molinfo[i].size();
        for (int j = 0; j < natm; ++j) {
            int label = molinfo[i][j];
            for (int k = 0; k < 3; ++k) {
                ExCoord[l][k] = coord[label-1][k];
            }
            l += 1;
        }
    }

    return ExCoord;
}
