// read_file.hpp
#ifndef READ_FILE_HPP
#define READ_FILE_HPP

#include <vector>
#include <string>

using namespace std;

struct s_region {
    vector<double> x;
    vector<double> y;
    vector<int>    reg;
};

struct s_col1others {
    vector<int> col1;
    vector<vector<double>> others;
};

s_col1others ReadCol1others(const string filename);

int linenum(const string& filename);

vector<double> Read1dfile(const string& filename);

vector<string> ReadList(const string& filename);

vector<vector<double>> Read2dfile(const string& filename);

vector<double> Readaxisfile(const string& filename);

vector<vector<double>> 
Read2dmapfile(const string& filename);

vector<vector<double>> Readxyfile(const string& filename);

vector<vector<vector<double>>> 
ReadMultixyfile(const vector<string>& filelist);

vector<vector<vector<int>>> 
ReadMultixyfileINT(const vector<string>& filelist);

s_region Readxyreg(const string& filename);

struct s_dx {
    int    ng3[3];
    double origin[3];
    double delta[3][3];
    vector<vector<vector<double>>> data;
};

s_dx ReaddxFile(const string& filename);

pair<vector<string>, vector<double>>
ReadChargefile(const string& filename);

#endif
