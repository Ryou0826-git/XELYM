// nc_tools.cpp

#include <nc_tools.hpp>

NcVar FindVar(const NcFile &dataFile, 
              const string& keyword) {
    //
    NcVar keyvar = dataFile.getVar(keyword);
    //
    return keyvar;
}

size_t GetDim(const NcFile &dataFile, 
              const string &keyword) {
    //
    NcDim KeyDim = dataFile.getDim(keyword);
    size_t KeySize = KeyDim.getSize();
    //
    return KeySize;
}

vector<NcVar*> getAllVariables(const NcFile &dataFile) {
    //
    multimap<string, NcVar> varMap = dataFile.getVars();
    vector<NcVar*> varList;
    //
    for (auto &entry : varMap) {
        varList.push_back(&entry.second);
    }
    //
    return varList;
}

NcVar* FindVar(const vector<NcVar*>& varList, 
               const string& keyword) {
    //
    NcVar* keyvar = nullptr;
    //
    for (const auto &var : varList) {
        if (var->getName() == keyword) {
            keyvar = var;
            break;
        }
    }
    //
    return keyvar;
}

vector<vector<double>> SplitCoordOneframe(const NcVar& coordvar, 
                                          const int& natm, 
                                          const int& f) {
    if (f <= 0) {
        throw invalid_argument("Frame index must be >= 1");
    }

    vector<size_t> startp = {static_cast<size_t>(f - 1), 0, 0};
    vector<size_t> countp = {1, static_cast<size_t>(natm), 3};

    vector<double> coord_buffer(natm * 3, 0.0);
    vector<vector<double>> coord(natm, vector<double>(3, 0.0));

    coordvar.getVar(startp, countp, coord_buffer.data());

    for (int i = 0; i < natm; ++i) {
        coord[i][0] = coord_buffer[i * 3];
        coord[i][1] = coord_buffer[i * 3 + 1];
        coord[i][2] = coord_buffer[i * 3 + 2];
    }

    return coord;
}

////vector<float> CoordOneframe(const NcVar& coordvar, 
////                            const int& natm, 
////                            const int& f) {
////    if (f <= 0) {
////        throw invalid_argument("Frame index must be >= 1");
////    }
////
////    vector<size_t> startp = {static_cast<size_t>(f - 1), 0, 0};
////    vector<size_t> countp = {1, static_cast<size_t>(natm), 3};
////
////    vector<float> coord(natm * 3, 0.0);
////
////    coordvar.getVar(startp, countp, coord.data());
////
////    return coord;
////}

//-----------------------------------------------------------------------------
//                      Read netcdf format trajectory
//-----------------------------------------------------------------------------
vector<vector<vector<double>>> GetCoordfromNC(const string& trajfile) {
    vector<vector<vector<double>>> coord;

    NcFile dataFile(trajfile, NcFile::read);
    NcVar coordvar = FindVar(dataFile, "coordinates");
    size_t frames  = GetDim(dataFile, "frame");
    size_t natm    = GetDim(dataFile, "atom");

    for (size_t f = 0; f < frames; ++f) {
        vector<vector<double>> cell_coord;

        for (size_t n = 0; n < natm; ++n) {
            vector<double> xyz(3);
            vector<size_t> startp = {f, n, 0};
            vector<size_t> countp = {1, 1, 3};

            coordvar.getVar(startp, countp, xyz.data());
            cell_coord.push_back(xyz);
        }

        coord.push_back(cell_coord);
    }

    return coord;
}

//vector<float> GetOneframefromNC_vec1_float(const string& trajfile, 
//                                           const int f) {
//    //
//    NcFile dataFile(trajfile, NcFile::read);
//    NcVar coordvar = FindVar(dataFile, "coordinates");
//    size_t frames  = GetDim(dataFile, "frame");
//    size_t natm    = GetDim(dataFile, "atom");
//    
//    if (f <= 0) {
//        throw invalid_argument("Frame index must be >= 1");
//    }
//
//    vector<size_t> startp = {static_cast<size_t>(f - 1), 0, 0};
//    vector<size_t> countp = {1, static_cast<size_t>(natm), 3};
//
//    vector<float> coord(natm * 3, 0.0);
//
//    coordvar.getVar(startp, countp, coord.data());
//
//    return coord;
//}

vector<vector<vector<double>>> 
GetSelCoordfromNC(const vector<string>& trajlist, 
                  const vector<int>& sel_list) {
    //
    vector<vector<vector<double>>> coord;

    for (const string& trajfile : trajlist) {
        //
        if (!fs::exists(trajfile)) throw runtime_error("File not found: " + trajfile);
        //
        NcFile dataFile(trajfile, NcFile::read);
        NcVar coordvar = FindVar(dataFile, "coordinates");
        size_t frames  = GetDim(dataFile, "frame");

        for (size_t f = 0; f < frames; ++f) {
            vector<vector<double>> cell_coord;
            for (int n : sel_list) {
                vector<double> xyz(3);
                vector<size_t> startp = {f, static_cast<size_t>(n-1), 0};
                vector<size_t> countp = {1, 1, 3};

                coordvar.getVar(startp, countp, xyz.data());
                cell_coord.push_back(xyz);
            }

            coord.push_back(cell_coord);
        }
    }

    return coord;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Read restart file (Amber)
//-----------------------------------------------------------------------------
////s_ncrst ReadNcrstAmber(const string& rstfile) {
////    s_ncrst ncrst;
////    //
////    NcFile dataFile(rstfile, NcFile::read);
////    NcVar coordvar = FindVar(dataFile, "coordinates");
////    NcVar velocityvar = FindVar(dataFile, "velocities");
////    size_t natm    = GetDim(dataFile, "atom");
////    ncrst.natm = natm;
////    //
////    cout << endl;
////    cout << "Read coordinates from amber netcdf rstfile..." << endl;
////    for (size_t n = 0; n < natm; ++n) {
////        vector<double> xyz(3);
////        vector<size_t> startp = {n, 0};
////        vector<size_t> countp = {1, 3};
////        //
////        coordvar.getVar(startp, countp, xyz.data());
////        ncrst.coord.push_back(xyz);
////    }
////    //
////    cout << endl;
////    cout << "Read velocities from amber netcdf rstfile..." << endl;
////    for (size_t n = 0; n < natm; ++n) {
////        vector<double> xyz(3);
////        vector<size_t> startp = {n, 0};
////        vector<size_t> countp = {1, 3};
////        //
////        velocityvar.getVar(startp, countp, xyz.data());
////        ncrst.velocity.push_back(xyz);
////    }
////
////    return ncrst;
////}

s_ncrst ReadNcrstAmber(const string& rstfile) {
    s_ncrst ncrst;
    //
    NcFile dataFile(rstfile, NcFile::read);
    NcVar coordvar    = FindVar(dataFile, "coordinates");
    NcVar velocityvar = FindVar(dataFile, "velocities");
    NcVar boxVar      = FindVar(dataFile, "cell_lengths");
    
    size_t natm       = GetDim(dataFile, "atom");
    ncrst.natm = natm;
    
    ncrst.coord_x.resize(natm);
    ncrst.coord_y.resize(natm);
    ncrst.coord_z.resize(natm);
    ncrst.velocity_x.resize(natm);
    ncrst.velocity_y.resize(natm);
    ncrst.velocity_z.resize(natm);

    //
    // coord section
    cout << endl;
    cout << "Read coordinates and velocities from amber netcdf rstfile..." << endl;
    
    vector<size_t> startp = {0, 0};
    vector<size_t> countp = {static_cast<size_t>(natm), 1};
    coordvar.getVar(startp, countp, ncrst.coord_x.data());
    velocityvar.getVar(startp, countp, ncrst.velocity_x.data());
    
    startp = {0, 1};
    coordvar.getVar(startp, countp, ncrst.coord_y.data());
    velocityvar.getVar(startp, countp, ncrst.velocity_y.data());
    
    startp = {0, 2};
    coordvar.getVar(startp, countp, ncrst.coord_z.data());
    velocityvar.getVar(startp, countp, ncrst.velocity_z.data());

    vector<double> box(3);
    boxVar.getVar(box.data());
    ncrst.box_size_x = box[0];
    ncrst.box_size_y = box[1];
    ncrst.box_size_z = box[2];

    return ncrst;
}
//-----------------------------------------------------------------------------
