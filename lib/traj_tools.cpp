// traj_tools.cpp

#include <dcd_r.hpp>
#include <psf_tools.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;

void AllocateDCD(s_dcd& dcd, int natms, int frames) {
    int dim = 3;
    //
    dcd.coord = vector<vector<vector<double>>>(
        frames, vector<vector<double>>(natms, vector<double>(dim, 0)
        )
    );
    //
    dcd.box   = vector<vector<double>>( 
        frames, vector<double>(dim, 0)
    );
}

void readAndStoreCoordBox(DCD_R& dcdf, 
                          s_dcd& dcd, 
                          int nAtoms, 
                          int frames) {
    const float *x_dumm, *y_dumm, *z_dumm;
    const double *pbc;

    for (int i = 0; i < frames; ++i) {
        dcdf.read_oneFrame();
        x_dumm = dcdf.getX();
        y_dumm = dcdf.getY();
        z_dumm = dcdf.getZ();
        //
        pbc    = dcdf.getPbc();
        dcd.box[i][0] = pbc[0];
        dcd.box[i][1] = pbc[2];
        dcd.box[i][2] = pbc[5];

        for (int j = 0; j < nAtoms; ++j) {
            dcd.coord[i][j][0] = x_dumm[j];
            dcd.coord[i][j][1] = y_dumm[j];
            dcd.coord[i][j][2] = z_dumm[j];
        }
    }
}

void CatCoordinate(const s_dcd& dcd, 
                   const vector<vector<double>>& ligcoord, 
                   const int& frame, 
                   const int& insf, 
                   s_dcd& dcdins) {
    //
    int natm = dcd.coord[0].size();
    int nlig = ligcoord.size();

    for (int i = 0; i < 3; ++i) {
        dcdins.box[insf][i] = dcd.box[insf][i];
    }

    for (int i = 0; i < natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            dcdins.coord[insf][i][j] = dcd.coord[frame-1][i][j];
        }
    }

    for (int i = 0; i < nlig; ++i) {
        for (int j = 0; j < 3; ++j) {
            dcdins.coord[insf][i+natm][j] = ligcoord[i][j];
        }
    }
}

s_dcd ReaddcdFile(const string& dcdfile, 
                  const s_psf& psf) {
    DCD_R dcdf(dcdfile.c_str());
    dcdf.read_header();

    cout << endl;
    cout << "Read " << dcdfile << endl;
    int frames = dcdf.getNFILE();
    int natm   = psf.resname.size();
    //
    s_dcd dcd;
    AllocateDCD(dcd, natm, frames);
    readAndStoreCoordBox(dcdf, dcd, natm, frames);
    //
    return dcd;
}

s_dcd ReaddcdFileAmber(const string& dcdfile, 
                       const int& natm) {
    DCD_R dcdf(dcdfile.c_str());
    dcdf.read_header();

    cout << endl;
    cout << "Read " << dcdfile << endl;
    int frames = dcdf.getNFILE();
    //
    s_dcd dcd;
    AllocateDCD(dcd, natm, frames);
    readAndStoreCoordBox(dcdf, dcd, natm, frames);
    //
    return dcd;
}


