// nc_tools.hpp
#ifndef NC_TOOLS_HPP
#define NC_TOOLS_HPP

#include <vector>
#include <string>
//#include <netcdf.h>
#include <netcdf>
#include <iostream>
#include <fstream>
#include <map>
#include <filesystem>

using namespace netCDF;
using namespace std;
namespace fs = std::filesystem;

struct s_ncrst {
    size_t         natm;
    vector<char>   bnatm;
    vector<double> coord_x;
    vector<char>   bcoord_x;
    vector<double> coord_y;
    vector<char>   bcoord_y;
    vector<double> coord_z;
    vector<char>   bcoord_z;
    vector<double> velocity_x;
    vector<char>   bvelocity_x;
    vector<double> velocity_y;
    vector<char>   bvelocity_y;
    vector<double> velocity_z;
    vector<char>   bvelocity_z;
    //
    double       box_size_x;
    vector<char> bbox_size_x;
    double       box_size_y;
    vector<char> bbox_size_y;
    double       box_size_z;
    vector<char> bbox_size_z;
};

NcVar FindVar(const NcFile &dataFile,
              const string& keyword);

size_t GetDim(const NcFile &dataFile,
              const string &keyword);

vector<NcVar*> getAllVariables(const NcFile &dataFile);

NcVar* FindVar(const vector<NcVar*>& varList,
               const string& keyword);

vector<vector<double>> SplitCoordOneframe(const NcVar& coordvar,
                                         const int& natm,
                                         const int& f);

//GetOneframefromNC_vec1_float(const string& trajfile);
vector<vector<vector<double>>> GetCoordfromNC(const string& trajfile);
vector<vector<vector<double>>> GetSelCoordfromNC(const vector<string>& trajfile,
                                                 const vector<int>& sel_list);

s_ncrst ReadNcrstAmber(const string& rstfile);

#endif
