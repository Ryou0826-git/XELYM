// read_file.cpp

#include <read_file.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;

int linenum(const string& filename) {
     ifstream file(filename);
     string line;
     //
     int i = 0;
     //
     if (file.is_open()) {
         while (getline(file, line)) {
             i += 1;
         }
     } else {
         cerr << "Unable to open file: "\
                   << filename << endl;
     }
     return i;
}

vector<double> Read1dfile(const string& filename) {
     ifstream file(filename);
     string line;
     //
     int num = linenum(filename);
     vector<double> data(num, 0.0);
     //
     if (file.is_open()) {
         for (int i = 0; i < num && getline(file, line); i++) {
             double dumm;
             istringstream ss(line);
             ss >> dumm >> data[i];
         }
     } else {
         cerr << "Unable to open file: " << filename << endl;
     }
     return data;
}

vector<string> ReadList(const string& filename) {
    ifstream file(filename);
    string line;
    //
    int num = linenum(filename);
    vector<string> data(num);
    //
    if (file.is_open()) {
        for (int i = 0; i < num && getline(file, line); i++) {
            istringstream ss(line);
            ss >> data[i];
        }
    } else {
        cerr << "Unable to open file: " << filename << endl;
    }
    return data;
}

vector<vector<double>> Read2dfile(const string& filename) {
     ifstream file(filename);
     string line;
     //
     int num = linenum(filename);
     vector<vector<double>> 
     data(2, vector<double>(num, 0.0));
     //
     if (file.is_open()) {
         for (int i = 0; i < num && getline(file, line); i++) {
             double dumm;
             istringstream ss(line);
             ss >> dumm >> data[0][i] >> data[1][i];
         }
     } else {
         cerr << "Unable to open file: " << filename << endl;
     }
     return data;
}

vector<vector<double>> Readxyfile(const string& filename) {
     ifstream file(filename);
     string line;
     //
     int num = linenum(filename);
     vector<vector<double>> data(2, vector<double>(num, 0.0));
     //
     if (file.is_open()) {
         for (int i = 0; i < num && getline(file, line); i++) {
             istringstream ss(line);
             ss >> data[0][i] >> data[1][i];
         }
     } else {
         cerr << "Unable to open file: " << filename << endl;
     }
     return data;
}

vector<vector<vector<double>>> 
ReadMultixyfile(const vector<string>& filelist) {
     int nfile = filelist.size();
     vector<vector<vector<double>>> data(nfile);
     for (int i = 0; i < nfile; ++i) {
         string filename = filelist[i];
         ifstream file(filename);
         string line;
         //
         int num = linenum(filename);
         data[i] = vector<vector<double>>(2, vector<double>(num, 0.0));
         //
         if (file.is_open()) {
             for (int j = 0; j < num && getline(file, line); j++) {
                 istringstream ss(line);
                 ss >> data[i][0][j] >> data[i][1][j];
             }
         } else {
             cerr << "Unable to open file: " << filename << endl;
         }
     }
     return data;
}

vector<vector<vector<int>>> 
ReadMultixyfileINT(const vector<string>& filelist) {
     int nfile = filelist.size();
     vector<vector<vector<int>>> data(nfile);
     for (int i = 0; i < nfile; ++i) {
         string filename = filelist[i];
         ifstream file(filename);
         string line;
         //
         int num = linenum(filename);
         data[i] = vector<vector<int>>(2, vector<int>(num));
         //
         if (file.is_open()) {
             for (int j = 0; j < num && getline(file, line); j++) {
                 istringstream ss(line);
                 ss >> data[i][0][j] >> data[i][1][j];
             }
         } else {
             cerr << "Unable to open file: " << filename << endl;
         }
     }
     return data;
}

s_region Readxyreg(const string& filename) {
     ifstream file(filename);
     string line;
     //
     int num = linenum(filename);
     double dumm;
     s_region data;
     data.x.resize(num, 0.0);
     data.y.resize(num, 0.0);
     data.reg.resize(num, 0);

     //
     if (file.is_open()) {
         for (int i = 0; i < num && getline(file, line); i++) {
             istringstream ss(line);
             ss >> dumm >> data.x[i] >> data.y[i] >> data.reg[i];
         }
     } else {
         cerr << "Unable to open file: " << filename << endl;
     }

     return data;
}

vector<double> Readaxisfile(const string& filename) {
    //
    vector<double> data;
    double num;
    ifstream infile(filename);
    //
    if (!infile) {
         cerr << "Cannot read axis file..." << endl;
    }
    //
    while (infile >> num) {
        data.push_back(num);
    }
    //
    return data;
}

s_col1others ReadCol1others(const string filename) {
    //
    s_col1others col1others;
    //
    ifstream infile(filename);
    if (!infile) {
        cerr << "Unable to open " << filename << endl;
    }

    string line;
    while (getline(infile, line)) {
        if (line.empty()) continue;

        istringstream iss(line);

        double first;
        iss >> first;

        col1others.col1.push_back(static_cast<int>(first));

        vector<double> row;
        double val;
        while (iss >> val) {
            row.push_back(val);
        }
        col1others.others.push_back(row);
    }
    infile.close();
    //
    return col1others;
}

vector<vector<double>> 
Read2dmapfile(const string& filename) {
    ifstream infile(filename);
    //
    if (!infile) {
         cerr << "Cannot read 2d-map file..." << endl;
    }
    vector<vector<double>> mapdata;
    string line;

    while (getline(infile, line)) {
        istringstream ss(line);
        vector<double> data;
        double z;
        //
        while (ss >> z) {
            data.push_back(z);
        }
        mapdata.push_back(data);
    }
    return mapdata;
}


s_dx ReaddxFile(const string& filename) {
    //
    s_dx dx;
    ifstream infile(filename);
    string line;
    //
    if (!infile) {
         cerr << "Cannot read dxfile..." << endl;
    }
    //
    int i = 0;
    bool dataSection = false;
    //
    while (getline(infile, line)) {
        istringstream ss(line);
        if (line.find("object 1 class gridpositions counts") != string::npos) {
            ss.ignore(35);
            ss >> dx.ng3[0] >> dx.ng3[1] >> dx.ng3[2];
        } else if (line.find("origin") != string::npos) {
            ss.ignore(6);
            ss >> dx.origin[0] >> dx.origin[1] >> dx.origin[2];
        } else if (line.find("delta") != string::npos) {
            ss.ignore(5);
            ss >> dx.delta[i][0] >> dx.delta[i][1] >> dx.delta[i][2];
            i += 1;
        } else if (line.find("data follows") != string::npos) {
            dataSection = true;
            break;
        }
    }
    //
    dx.data.resize(dx.ng3[0], vector<vector<double>>(dx.ng3[1], vector<double>(dx.ng3[2], 0.0)));
    if (dataSection) {
        int nx = 0, ny = 0, nz = 0;
        double value;
        while (infile >> value) {
            dx.data[nx][ny][nz] = value;
            nz += 1;
            if (nz == dx.ng3[2]) {
                nz = 0;
                ny += 1;
                if (ny == dx.ng3[1]) {
                    ny = 0;
                    nx += 1;
                }
            }
        }
    }
    infile.close();
    cout << "Finish reading dxfile..." << endl;

    return dx;
}

pair<vector<string>, vector<double>> 
ReadChargefile(const string& filename) {
    ifstream fin(filename);
   if (!fin) {
        throw runtime_error("Error: cannot open" + filename);
    }

    vector<string> names;
    vector<double> charges;

    string line;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        istringstream iss(line);

        string name;
        double value;

        if (iss >> name >> value) {
            names.push_back(name);
            charges.push_back(value);
        }
    }

    return make_pair(names, charges);
}
