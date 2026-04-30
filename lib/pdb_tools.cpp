// pdb_tools.cpp

#include <pdb_tools.hpp>

//---------------------------------------------------------------------------//
// tool of contains
//---------------------------------------------------------------------------//
bool contains(const vector<vector<int>>& matrix, 
              const int& target) {
    //
    for (const auto& row : matrix) {
        for (int value : row) {
            if (value == target) return true;
        }
    }
    //
    return false;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<string>> 
PDBOPR::Readpdball(const string& pdbfile) {
    int upper = 1000000;
    ifstream file(pdbfile);
    vector<vector<string>> 
    line(2, vector<string>(upper));
    string dumm;
    int natm = -1;
    //
    if (file.is_open()) {
        while (getline(file, dumm)) {
            if (dumm.find("ATOM") != string::npos) {
                natm += 1;
                line[0][natm] = dumm.substr(0, 29);
                line[1][natm] = dumm.substr(54, 77);
            }
        }
    } else {
        cerr << "Unable to open file: " << pdbfile << endl;
    }
    //
    line[0].resize(natm+1);
    line[1].resize(natm+1);
    //
    return line;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_pdb PDBOPR::LoadfrompdbNeo(const string& pdbfile) {
    ifstream file(pdbfile);
    string line;
    s_pdb  pdb;
    bool   sign = true;
    //
    if (file.is_open()) {
        while (getline(file, line)) {
            if (line.find("ATOM") != string::npos) {
                string dumm1 = line.substr(0, 3);
                string dumm2 = line.substr(4, 7);
                string atmnm = line.substr(11, 6);
                string resnm = line.substr(17, 4);
                string chainid  = line.substr(21, 1);
                int    resid = stoi(line.substr(22, 4));
                //
                double x = stod(line.substr(31, 8));
                double y = stod(line.substr(38, 8));
                double z = stod(line.substr(46, 8));
                //
                double occupancy = stod(line.substr(54, 5));
                double tempFactor = stod(line.substr(60, 5));
                //
                string segnm = line.substr(72, 4);
                //
                // Push back info. to list...
                pdb.atmnm.push_back(atmnm);
                pdb.resnm.push_back(resnm);
                pdb.chainid.push_back(chainid);
                pdb.resid.push_back(resid);
                pdb.coord.push_back({x, y, z});
                pdb.occupancy.push_back(occupancy);
                pdb.tempFactor.push_back(tempFactor);
                pdb.segnm.push_back(segnm);
                if ((line.size() > 77) && (sign == true)) {
                    string katm  = line.substr(77, 1);
                    pdb.katm.push_back(katm);
                } else if (line.size() <= 76) {
                    sign = false;
                }
            }
        }
    } else {
        cerr << "Unable to open file: " << pdbfile << endl;
    }
    return pdb;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_pdb PDBOPR::LoadfrompdbAmber(const string& pdbfile) {
    ifstream file(pdbfile);
    string line;
    s_pdb  pdb;
    //
    if (file.is_open()) {
        while (getline(file, line)) {
            if (line.find("ATOM") != string::npos) {
                string dumm1 = line.substr(0, 3);
                string dumm2 = line.substr(4, 7);
                string atmnm = line.substr(11, 6);
                string resnm = line.substr(17, 4);
                int    resid = stoi(line.substr(22, 4));
                //
                double x = stod(line.substr(30, 8));
                double y = stod(line.substr(38, 8));
                double z = stod(line.substr(46, 8));
                //
                double occupancy = stod(line.substr(54, 5));
                double tempFactor = stod(line.substr(60, 5));
                //
                // Push back info. to list...
                pdb.atmnm.push_back(atmnm);
                pdb.resnm.push_back(resnm);
                pdb.resid.push_back(resid);
                pdb.coord.push_back({x, y, z});
                pdb.occupancy.push_back(occupancy);
                pdb.tempFactor.push_back(tempFactor);
            }
            //
            if (line.find("CRYST1") != string::npos) {
                stringstream ss(line.substr(6));
                double a, b, c, alpha, beta, gamma;
                ss >> a >> b >> c >> alpha >> beta >> gamma;
                //
                pdb.box.push_back(a);
                pdb.box.push_back(b);
                pdb.box.push_back(c);
                pdb.box.push_back(alpha);
                pdb.box.push_back(beta);
                pdb.box.push_back(gamma);
            }
        }
    } else {
        cerr << "Unable to open file: " << pdbfile << endl;
    }
    return pdb;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_pdb PDBOPR::Removepdb(const s_pdb& pdb, 
                        const vector<int> nonbond_list) {
    s_pdb remove_pdb;
    //
    int natm = pdb.atmnm.size();
    for (int i = 0; i < natm; ++i) {
        if (std::find(nonbond_list.begin(), nonbond_list.end(), i+1) == nonbond_list.end()) {
            remove_pdb.atmnm.push_back(pdb.atmnm[i]);
            remove_pdb.resnm.push_back(pdb.resnm[i]);
            remove_pdb.resid.push_back(pdb.resid[i]);
            remove_pdb.coord.push_back(pdb.coord[i]);
            remove_pdb.occupancy.push_back(pdb.occupancy[i]);
            remove_pdb.tempFactor.push_back(pdb.tempFactor[i]);
        }
    }
    return remove_pdb;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_pdb PDBOPR::Loadfrompdb(const string& pdbfile, 
                          const int& natm) {
    ifstream file(pdbfile);
    string line;
    s_pdb  pdb;
    string dumm1;
    int    dumm2;
    string ast = "*****";
    //
    // Allocate memory...
    pdb.atmnm.resize(natm);
    pdb.resnm.resize(natm);
    pdb.chainid.resize(natm);
    pdb.resid.resize(natm);
    pdb.coord.resize(natm, vector<double>(3, 0.0));
    pdb.segnm.resize(natm);
    //
    int i = 0;
    if (file.is_open()) {
        while (getline(file, line)) {
            if (line.find("ATOM") != string::npos) {
                //
                if (line.find(ast) != std::string::npos) {
                    size_t pos;
                    while ((pos = line.find(ast)) != std::string::npos) {
                        line.erase(pos, ast.length());  // erase "****"
                    }
                    istringstream ss(line);
                    ss >> dumm1           >>
                          pdb.atmnm[i]    >>
                          pdb.resnm[i]    >>
                          pdb.chainid[i]  >>
                          pdb.resid[i]    >>
                          pdb.coord[i][0] >>
                          pdb.coord[i][1] >>
                          pdb.coord[i][2];
                    i += 1;
                } else {
                    istringstream ss(line);
                    ss >> dumm1           >>
                          dumm2           >>
                          pdb.atmnm[i]    >>
                          pdb.resnm[i]    >>
                          pdb.chainid[i]  >>
                          pdb.resid[i]    >>
                          pdb.coord[i][0] >>
                          pdb.coord[i][1] >>
                          pdb.coord[i][2];
                    i += 1;
                }
            }
        }
    } else {
        cerr << "Unable to open file: " << pdbfile << endl;
    }
    return pdb;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_pdb PDBOPR::ExtractSegnm(const s_pdb& refspdb, 
                           const string& segnm) {
    //
    s_pdb segnmpdb;
    int n = refspdb.resnm.size();
    //
    for (int i = 0; i < n; ++i) {
        string name = refspdb.segnm[i];
        name.erase(name.begin(), std::find_if_not(name.begin(), name.end(), ::isspace));
        name.erase(std::find_if_not(name.rbegin(), name.rend(), ::isspace).base(), name.end());
        if (boost::iequals(name, segnm)){
            segnmpdb.atmnm.push_back(refspdb.atmnm[i]);
            segnmpdb.resnm.push_back(refspdb.resnm[i]);
            segnmpdb.chainid.push_back(refspdb.chainid[i]);
            segnmpdb.resid.push_back(refspdb.resid[i]);
            segnmpdb.coord.push_back(refspdb.coord[i]);
            segnmpdb.occupancy.push_back(refspdb.occupancy[i]);
            segnmpdb.tempFactor.push_back(refspdb.tempFactor[i]);
            segnmpdb.segnm.push_back(refspdb.segnm[i]);
            //
            if (refspdb.katm.size() > 0) {
                segnmpdb.katm.push_back(refspdb.katm[i]);
            }
        }
    }

    return segnmpdb;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>> 
PDBOPR::getmolinfo(const s_pdb& pdb, 
                   const vector<int>& reslist) {
    //
    int num = pdb.atmnm.size();
    int nlist = reslist.size();
    int dumm = 100;
    int nres = 0;
    //
    vector<vector<int>> molinfo(nlist, vector<int>(dumm, 0));
    //
    int j = 0;
    for (int i = 0; i < num; ++i) {
       if (pdb.resid[i] == reslist[j]) {
           molinfo[j][nres] = i + 1;
           nres += 1;
           if (pdb.resid[i+1] != reslist[j]) {
               molinfo[j].resize(nres);
               j += 1;
               nres = 0;
           }
       }
       if (j == nlist) {
           break;
       }
    }
    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>> 
PDBOPR::getmolinfo_atmnm(const s_pdb& pdb, 
                         const vector<int>& reslist, 
                         const string& atmnm) {
    //
    int num = pdb.atmnm.size();
    int nlist = reslist.size();
    //
    vector<vector<int>> molinfo(nlist, vector<int>(1, 0));
    //
    int j = 0;
    for (int i = 0; i < num; ++i) {
       if ((pdb.atmnm[i] == atmnm) && (pdb.resid[i] == reslist[j])) {
           molinfo[j][0] = i + 1;
           j += 1;
       }
       //
       if (j == nlist) {
           break;
       }
    }
    //
    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>>
PDBOPR::getmolinfo_atmnm_all(const s_pdb& pdb, 
                             const string& atmnm) {
    //
    int num = pdb.atmnm.size();
    vector<vector<int>> molinfo;
    vector<int> cellinfo;

    for (int i = 0; i < num; ++i) {
        string atmtrim = pdb.atmnm[i];
        atmtrim.erase(std::remove_if(atmtrim.begin(), atmtrim.end(), ::isspace), atmtrim.end());
        if (boost::iequals(atmtrim, atmnm)) {
            cellinfo.push_back(i+1);
        }
        if (i < num-1) {
            if (pdb.resid[i] != pdb.resid[i+1]) {
                molinfo.push_back(cellinfo);
                cellinfo.clear();
            }
        } else if (i == num-1) {
            molinfo.push_back(cellinfo);
            cellinfo.clear();
        }
    }

    return molinfo;
}

vector<vector<int>> 
PDBOPR::getmolinfo_atmnm_all_amber(const s_pdb& pdb, 
                                   const string& atmnm) {
    //
    int num = pdb.atmnm.size();
    vector<vector<int>> molinfo;
    vector<int> cellinfo;

    for (int i = 0; i < num; ++i) {
        //
        string atmtrim = pdb.atmnm[i];
        atmtrim.erase(std::remove_if(atmtrim.begin(), atmtrim.end(), ::isspace), atmtrim.end());
        //
        if (boost::iequals(atmtrim, atmnm)) {
            cellinfo.push_back(i+1);
        }
    }

    molinfo.push_back(cellinfo);

    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void PDBOPR::ScaleMOL(s_pdb& pdb, 
                      const double& scale) {
    //
    int natm = pdb.coord.size();
    vector<vector<double>> ScaleCoord(natm, vector<double>(3, 0.0));
    //
    for (int i = 0; i < natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            pdb.coord[i][j] *= scale;
        }
    }
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void PDBOPR::PlungeCoord(s_pdb& pdb, 
                         const vector<vector<double>>& coord){
    //
    int natm = pdb.coord.size();
    //
    for (int i = 0; i < natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            pdb.coord[i][j] = coord[i][j];
        }
    }
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void PDBOPR::Outpdb(const vector<vector<string>>& line, 
                    const vector<vector<double>>& coord, 
                    const string& filename) {
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    int natm = line[0].size();
    //
    for (int i = 0; i < natm; ++i) {
        outfile << line[0][i] << " "
                << fixed << setprecision(3) << setw(8) << coord[i][0]
                << fixed << setprecision(3) << setw(8) << coord[i][1]
                << fixed << setprecision(3) << setw(8) << coord[i][2]
                << line[1][i] << endl;
    }
    outfile.close();
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_pdb CatPdb(const vector<s_pdb>& pdblist) {
    s_pdb pdb;
    for (s_pdb cellpdb : pdblist) {
        int n = cellpdb.resnm.size();
        for (int i = 0; i < n; ++i) {
            pdb.atmnm.push_back(cellpdb.atmnm[i]);
            pdb.resnm.push_back(cellpdb.resnm[i]);
            pdb.chainid.push_back(cellpdb.chainid[i]);
            pdb.resid.push_back(cellpdb.resid[i]);
            pdb.coord.push_back(cellpdb.coord[i]);
        }
    }
    return pdb;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void PDBOPR::CatandOutpdb(const s_pdb& pdb1, 
                          const s_pdb& pdb2, 
                          const vector<vector<string>>& line1, 
                          const vector<vector<string>>& line2, 
                          const string& filename) {
    int natm1 = pdb1.atmnm.size();
    int natm2 = pdb2.atmnm.size();
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    for (int i = 0; i < natm1; ++i) {
        outfile << line1[0][i] << " "
                << fixed  << setprecision(3) << setw(8) << pdb1.coord[i][0]
                << fixed  << setprecision(3) << setw(8) << pdb1.coord[i][1]
                << fixed  << setprecision(3) << setw(8) << pdb1.coord[i][2]
                << line1[1][i] << endl;
    }
    for (int i = 0; i < natm2; ++i) {
        outfile << line2[0][i] << " "
                << fixed  << setprecision(3) << setw(8) << pdb2.coord[i][0]
                << fixed  << setprecision(3) << setw(8) << pdb2.coord[i][1]
                << fixed  << setprecision(3) << setw(8) << pdb2.coord[i][2]
                << line2[1][i] << endl;
    }
    outfile.close();
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void PDBOPR::OutInsertpdb(const vector<vector<double>>& coord, 
                          const vector<vector<string>>& sysline, 
                          const vector<vector<string>>& ligline, 
                          const string& filename) {
    //
    int natm = coord.size();
    int nsys = sysline[0].size();
    int nlig = ligline[0].size();
    //
    if (natm != nsys + nlig) {
        cerr << "Error natm != natm of system + natm of ligand !!" << endl;
        return;
    }

    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    cout << "Output data of traj. to pdb file..." << endl;
    int j = 0;
    for (int i = 0; i < natm; ++i) {
        if (j < nsys) {
            outfile << sysline[0][i] << " "
                    << fixed << setprecision(3) << setw(8) << coord[i][0]
                    << fixed << setprecision(3) << setw(8) << coord[i][1]
                    << fixed << setprecision(3) << setw(8) << coord[i][2]
                    << sysline[1][i] << endl;
            j += 1;
        } else {
            outfile << ligline[0][i-nsys] << " "
                    << fixed << setprecision(3) << setw(8) << coord[i][0]
                    << fixed << setprecision(3) << setw(8) << coord[i][1]
                    << fixed << setprecision(3) << setw(8) << coord[i][2]
                    << ligline[1][i-nsys] << endl;
        }
    }
    outfile.close();
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void PDBOPR::Outxyz(const s_pdb& pdb, 
                    const vector<vector<double>>& coord, 
                    const string& filename) {
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    cout << "Output data of traj. to xyz file..." << endl;
    int natm = pdb.atmnm.size();
    outfile << natm << endl;
    outfile << "system" << endl;
    for (int i = 0; i < natm; ++i) {
        outfile << pdb.atmnm[i][0] << " "
                << fixed << setprecision(3) << setw(8) << coord[i][0]
                << fixed << setprecision(3) << setw(8) << coord[i][1]
                << fixed << setprecision(3) << setw(8) << coord[i][2]
                << endl;
    }
    outfile.close();
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void PDBOPR::OutPdbNeo(const vector<s_pdb>& pdblist, 
                       const string& filename) {
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    cout << "output pdb data to pdbfile..." << endl;
    int natm = 0;
    for (const s_pdb& pdb : pdblist) {
        int patm = pdb.resnm.size();
        for (int i = 0; i < patm; ++i) {
            if (natm < 99999) {
                outfile << "ATOM" 
                        << setw(7) << right << natm + 1;
            } else {
                outfile << "ATOM  *****";
            }
            //
            outfile << setw(6) << left  << pdb.atmnm[i] 
                    << setw(4) << right << pdb.resnm[i] 
                    << setw(1) << pdb.chainid[i] 
                    << setw(4) << right << pdb.resid[i] 
                    << fixed << setprecision(3) << setw(12) << right << pdb.coord[i][0]
                    << fixed << setprecision(3) << setw(8)  << right << pdb.coord[i][1]
                    << fixed << setprecision(3) << setw(8)  << right << pdb.coord[i][2]
                    << fixed << setprecision(2) << setw(6)  << pdb.occupancy[i] 
                    << fixed << setprecision(2) << setw(6)  << pdb.tempFactor[i] 
                    << "      "
                    << setw(4) << right << pdb.segnm[i];
            if (pdb.katm.size() > 0) {
                outfile << setw(2) << right << pdb.katm[i] << endl;
            } else {
                outfile <<  "  " << endl;
            }
            //
            natm += 1;
        }
    }

    outfile.close();
}

void PDBOPR::OutPdbAmber(const vector<s_pdb>& pdblist, 
                         const string& filename) {
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    cout << "output pdb data to pdbfile..." << endl;
    if (pdblist[0].box.size() == 6) {
        outfile << "CRYST1 "  
                << pdblist[0].box[0] << " " 
                << pdblist[0].box[1] << " "
                << pdblist[0].box[2] << " "
                << pdblist[0].box[3] << " "
                << pdblist[0].box[4] << " "
                << pdblist[0].box[5] << " "
                << "P 1 1" << endl;
    }
    int natm = 0;
    for (const s_pdb& pdb : pdblist) {
        int patm = pdb.resnm.size();
        for (int i = 0; i < patm; ++i) {
            //if (natm < 100000) {
            if (natm < 99999) {
                outfile << "ATOM" 
                        << setw(7) << right << natm + 1;
            } else {
                outfile << "ATOM  *****";
            }
            //
            //cout << pdb.atmnm[i] << endl;
            //cout << pdb.resnm[i] << endl;
            //cout << pdb.resid[i] << endl;
            //cout << pdb.coord[i][0] << endl;
            //cout << pdb.coord[i][1] << endl;
            //cout << pdb.coord[i][2] << endl;
            //cout << pdb.occupancy[i] << endl;
            //cout << pdb.tempFactor[i] << endl;
            //
            outfile << setw(6) << left  << pdb.atmnm[i] 
                    << setw(4) << right << pdb.resnm[i] 
                    << setw(4) << right << pdb.resid[i] 
                    << fixed << setprecision(3) << setw(12) << right << pdb.coord[i][0]
                    << fixed << setprecision(3) << setw(8)  << right << pdb.coord[i][1]
                    << fixed << setprecision(3) << setw(8)  << right << pdb.coord[i][2]
                    << fixed << setprecision(2) << setw(6)  << pdb.occupancy[i] 
                    << fixed << setprecision(2) << setw(6)  << pdb.tempFactor[i];
            if (pdb.katm.size() > 0) {
                outfile << setw(2) << right << pdb.katm[i] << endl;
            } else {
                outfile <<  "  " << endl;
            }
            //
            natm += 1;
        }
    }

    outfile.close();
}

void PDBOPR::OutPdbAmber2(const vector<s_pdb>& pdblist, 
                          const string& filename) {
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    //
    cout << "output pdb data to pdbfile..." << endl;
    int natm = 0;
    for (const s_pdb& pdb : pdblist) {
        int patm = pdb.resnm.size();
        for (int i = 0; i < patm; ++i) {
            //if (natm < 100000) {
            if (natm < 99999) {
                outfile << "ATOM" 
                        << setw(7) << right << natm + 1;
            } else {
                outfile << "ATOM  *****";
            }
            //
            //cout << pdb.atmnm[i] << endl;
            //cout << pdb.resnm[i] << endl;
            //cout << pdb.resid[i] << endl;
            //cout << pdb.coord[i][0] << endl;
            //cout << pdb.coord[i][1] << endl;
            //cout << pdb.coord[i][2] << endl;
            //cout << pdb.occupancy[i] << endl;
            //cout << pdb.tempFactor[i] << endl;
            //
            outfile << " " << setw(4) << left  << pdb.atmnm[i] 
                    << setw(4) << right << pdb.resnm[i] 
                    << "  " << setw(4) << right << pdb.resid[i] 
                    << fixed << setprecision(3) << setw(12) << right << pdb.coord[i][0]
                    << fixed << setprecision(3) << setw(8)  << right << pdb.coord[i][1]
                    << fixed << setprecision(3) << setw(8)  << right << pdb.coord[i][2]
                    << fixed << setprecision(2) << setw(6)  << pdb.occupancy[i] 
                    << fixed << setprecision(2) << setw(6)  << pdb.tempFactor[i];
            if (pdb.katm.size() > 0) {
                outfile << setw(2) << right << pdb.katm[i] << endl;
            } else {
                outfile <<  "  " << endl;
            }
            //
            natm += 1;
        }
    }

    outfile.close();
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Search protein
//---------------------------------------------------------------------------//
string TrimLett(const string& str) {
    size_t first = str.find_first_not_of(" \t\n\r\f\v");
    if (first == string::npos)
        return ""; 

    size_t last = str.find_last_not_of(" \t\n\r\f\v");
    return str.substr(first, last - first + 1);
}

vector<vector<int>> PDBOPR::SearchProtein(const s_pdb& pdb) {
    //
    vector<vector<int>> proinfo;
    vector<string> residue_names = {
        "ALA", "GLY", "SER", "THR", "LEU",
        "ILE", "VAL", "ASP", "GLU", "ASN",
        "GLN", "LYS", "ARG", "HIS", "PHE",
        "TYR", "TRP", "CYS", "MET", "PRO", 
        "HIE"
    };

    int label = 1;
    vector<int> cell_proinfo;

    for (const string& resnm : pdb.resnm) {
       string resnm_trim = TrimLett(resnm);
       if (find(residue_names.begin(), residue_names.end(), resnm_trim) != residue_names.end()) {
           //cout << resnm_trim << " is in the list.\n";
           cell_proinfo.push_back(label);
       }
       //} else {
       //    cout << resnm << " is not in the list.\n";
       //}
       label += 1;
    }

    proinfo.push_back(cell_proinfo);

    return proinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_pdb PDBOPR::Extractfrominfo(const s_pdb& refspdb, 
                              const vector<vector<int>>& molinfo) {
    //
    s_pdb pdb;
    //
    for (const vector<int>& resinfo : molinfo) {
        for (const int& i : resinfo) {
           pdb.atmnm.push_back(refspdb.atmnm[i-1]);
           pdb.resnm.push_back(refspdb.resnm[i-1]);
           //pdb.chainid.push_back(refspdb.chainid[i-1]);
           pdb.resid.push_back(refspdb.resid[i-1]);
           pdb.coord.push_back(refspdb.coord[i-1]);
           pdb.occupancy.push_back(refspdb.occupancy[i-1]);
           pdb.tempFactor.push_back(refspdb.tempFactor[i-1]);
        }
    }

    return pdb;
}

s_pdb PDBOPR::ExtractfrominfoNOT(const s_pdb& refspdb, 
                                 const vector<vector<int>>& molinfo) {
    //
    s_pdb pdb;
    //
    int n = refspdb.atmnm.size();
    vector<int> nonlist;

    //cout << n << endl;
    for (int i = 0; i < n; ++i) {
        bool sign = contains(molinfo, i + 1);
        //cout << n << endl;
        if (sign == false) {
            nonlist.push_back(i + 1);
            //cout << i + 1 << endl;
        }
    }
    
    for (const int& i : nonlist) {
       pdb.atmnm.push_back(refspdb.atmnm[i-1]);
       pdb.resnm.push_back(refspdb.resnm[i-1]);
       //pdb.chainid.push_back(refspdb.chainid[i-1]);
       pdb.resid.push_back(refspdb.resid[i-1]);
       pdb.coord.push_back(refspdb.coord[i-1]);
       pdb.occupancy.push_back(refspdb.occupancy[i-1]);
       pdb.tempFactor.push_back(refspdb.tempFactor[i-1]);
    }
    
    return pdb;
}
//---------------------------------------------------------------------------//
