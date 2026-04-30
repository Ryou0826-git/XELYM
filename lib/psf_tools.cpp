// psf_tools.cpp

#include <psf_tools.hpp>

using namespace std;

//---------------------------------------------------------------------------//
int extractNumber(const string& line) {
    stringstream ss(line);
    string numberStr;
    ss >> numberStr;

    try {
        return stoi(numberStr);
    } catch (const invalid_argument& e) {
        cerr << "Invalid number format: " << numberStr << endl;
        return 0;
    } catch (const out_of_range& e) {
        cerr << "Number out of range: " << numberStr << endl;
        return 0;
    }
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_psf PSFOPR::loadFile(const string& filename) {
    ifstream file(filename);
    string line;
    s_psf psf;
  
    if (file.is_open()) {
        while (getline(file, line)) {
            int  natm = 0;
            bool found = false;
  
            if (line.find("!NATOM") != string::npos) {
                found = true;
                natm  = extractNumber(line);
                //
                psf.segname.resize(natm);
                psf.molid.resize(natm);
                psf.resname.resize(natm);
                psf.atmname.resize(natm);
                psf.katm.resize(natm);
                psf.charge.resize(natm);
                psf.mass.resize(natm);
                //
                cout << "Number of atoms: " << natm << endl;
                if (found) {
                    for (int i = 0; i < natm && getline(file, line); i++) {
                        istringstream ss(line);
                        double dumm1, dumm2;
                        double mass, charge;
                        int    molid;
                        string segname, resname, atmname, katm;
                        //
                        ss >> dumm1   >> 
                              segname >> 
                              molid   >> 
                              resname >> 
                              atmname >> 
                              katm    >> 
                              charge  >> 
                              mass    >> 
                              dumm2;
                        //
                        psf.segname[i] = segname;
                        psf.molid[i]   = molid;
                        psf.resname[i] = resname;
                        psf.atmname[i] = atmname;
                        psf.katm[i]    = katm;
                        psf.charge[i]  = charge;
                        psf.mass[i]    = mass;
                    }
                }
                break;
            }
        }
        file.close();
    } else {
        cerr << "Unable to open file: " << filename << endl;
    }

    return psf;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>> 
PSFOPR::getmolinfo(const s_psf& psf, 
                   const string& resnm) {
    int num = psf.resname.size();
    int sumres = 0;
    int numres = 0;

    for (int i = 0; i < num; i++){
        if (psf.resname[i] == resnm) {
            sumres = sumres + 1;
            if (psf.molid[i+1] != psf.molid[i]) {
                numres = numres + 1;
            } else {
                if (psf.resname[i] != psf.resname[i+1]) {
                  numres = numres + 1;
                }
            }
        }
    }
    int nummol = sumres / numres;

    vector<vector<int>> molinfo(numres, vector<int>(nummol, 0));

    int j = 0; // About numres
    int k = 0; // About nummol
    int l = 0;
    int del = 0;

    for (int i = 0; i < num; i++) {
        if (psf.resname[i] == resnm) {
            if (psf.molid[i] == j + 1 + del) {
                molinfo[j][k] = i + 1;
                k += 1;
            }
            if (k == nummol) {
                k = 0;
                j += 1;
                l += 1;
                if (l == 9999) {
                    del += -9999;
                    l = 0;
                }
            }
        }
    }

    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>> 
PSFOPR::getmolinfo_segnm(const s_psf& psf, 
                         const string& segnm) {
    int num    = psf.segname.size();
    int sumseg = 0;
    int numseg = 0;

    for (int i = 0; i < num; i++){
        string segtrim;
        segtrim = psf.segname[i];
        segtrim.erase(remove_if(segtrim.begin(), segtrim.end(), ::isdigit), 
                      segtrim.end());
        if (segtrim == segnm) {
            sumseg = sumseg + 1;
            if (psf.segname[i+1] != psf.segname[i]) {
                numseg = numseg + 1;
            }
        }
    }
    //
    int nummol = sumseg / numseg;

    //
    // Allocate memory
    vector<vector<int>> molinfo(numseg, vector<int>(nummol, 0));

    int j = 0; // About numseg
    int k = 0; // About nummol

    for (int i = 0; i < num; ++i) {
        string segtrim;
        segtrim = psf.segname[i];
        segtrim.erase(remove_if(segtrim.begin(), segtrim.end(), ::isdigit), segtrim.end());
        if (segtrim == segnm) {
            molinfo[j][k] = i+1;
            if (psf.segname[i+1] != psf.segname[i]) {
                j += 1;
                k = 0;
            } else {
                k += 1;
            }
        }
    }

    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>> 
PSFOPR::getmolinfo_endtoend(const s_psf& psf, 
                            const string& segnm, 
                            const vector<string>& selatm) {
    int num    = psf.segname.size();
    int sumseg = 0;
    int numseg = 0;

    for (int i = 0; i < num; i++){
        string segtrim;
        segtrim = psf.segname[i];
        segtrim.erase(remove_if(segtrim.begin(), segtrim.end(), ::isdigit), segtrim.end());
        if (segtrim == segnm) {
            sumseg = sumseg + 1;
            if (psf.segname[i+1] != psf.segname[i]) {
                numseg = numseg + 1;
            }
        }
    }

    //
    // Allocate memory
    vector<vector<int>> molinfo(numseg, vector<int>(2, 0));

    int j = 0; // About numseg

    for (int i = 0; i < num; ++i) {
        string segtrim;
        segtrim = psf.segname[i];
        segtrim.erase(remove_if(segtrim.begin(), segtrim.end(), ::isdigit), segtrim.end());
        if (segtrim == segnm) {
            if (psf.atmname[i] == selatm[0]) {
                molinfo[j][0] = i+1;
            } else if (psf.atmname[i] == selatm[1]) {
                molinfo[j][1] = i+1;
            }
            if (psf.segname[i+1] != psf.segname[i]) {
                j += 1;
            }
        }
    }

    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>>
PSFOPR::getmolinfo_segres(const s_psf& psf,
                          const string& segnm,
                          const vector<int>& reslist) {
    int num   = psf.resname.size();
    int nlist = reslist.size();
    int nres  = 0;
    int dumm  = 100;
    //
    vector<vector<int>> proinfo(nlist, vector<int>(dumm, 0));
    //
    int j = 0;
    for (int i = 0; i < num; ++i) {
        if (boost::iequals(psf.segname[i], segnm)) {
            if (psf.molid[i] == reslist[j]) {
                proinfo[j][nres] = i + 1;
                nres += 1;
                if (psf.molid[i+1] != reslist[j]) {
                    proinfo[j].resize(nres);
                    j += 1;
                    nres = 0;
                }
            }
        }
        if (j == nlist) {
            break;
        }
    }
    return proinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>>
PSFOPR::getmolinfo_seg_all(const s_psf& psf,
                           const string& segnm) {
    int num   = psf.resname.size();
    int nres  = 0;
    int dumm1 = 10000;
    int dumm2 = 100;
    //
    vector<vector<int>> molinfo(dumm1, vector<int>(dumm2, 0));
    //
    int j = 0;
    for (int i = 0; i < num; ++i) {
        if (boost::iequals(psf.segname[i], segnm)) {
            molinfo[j][nres] = i + 1;
            nres += 1;
            if (psf.molid[i+1] != psf.molid[i] - 1) {
                molinfo[j].resize(nres);
                j += 1;
                nres = 0;
            }
        }
    }

    molinfo.resize(j);

    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>> 
PSFOPR::getmolinfo_segatm_res(const s_psf& psf, 
                              const string& segnm, 
                              const string& atmnm, 
                              const vector<int>& reslist) {
    int num  = psf.resname.size();
    int nlist = reslist.size();
    //
    vector<vector<int>> molinfo(nlist, vector<int>(1, 0));
    //
    int j = 0;
    for (int i = 0; i < num; ++i) {
        if ((boost::iequals(psf.segname[i], segnm))
            && (boost::iequals(psf.atmname[i], atmnm))
            && (psf.molid[i] == reslist[j])) {
            molinfo[j][0] = i + 1;
            j += 1;
        }
    }

    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<int>> 
PSFOPR::getmolinfo_segatm_all(const s_psf& psf, 
                              const string& segnm, 
                              const string& atmnm) {
    int num  = psf.resname.size();
    int dumm = 100000;
    //
    vector<vector<int>> molinfo(dumm, vector<int>(1, 0));
    //
    int j = 0;
    for (int i = 0; i < num; ++i) {
        if (psf.segname[i] == segnm) {
            if (psf.atmname[i] == atmnm) {
                molinfo[j][0] = i + 1;
                j += 1;
            }
        }
    }
    //
    molinfo.resize(j);

    return molinfo;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_branch PSFOPR::getmolinfo_segnm_heavy(const s_psf& psf, 
                                        const string& segnm) {
    s_branch branch;
    int num = psf.segname.size();

    vector<int> current_mol;
    vector<string> current_atmnm;
    //
    for (int i = 0; i < num; ++i) {
        if (boost::iequals(psf.segname[i], segnm) 
            && !boost::iequals(psf.atmname[i].substr(0, 1), "H")) {
            current_mol.push_back(i + 1);
            current_atmnm.push_back(psf.atmname[i]);
        }

        if (i == num - 1 
            || !boost::iequals(psf.segname[i], psf.segname[i + 1]) 
            || abs(psf.molid[i] - psf.molid[i + 1]) > 1) {
            if (!current_mol.empty()) {
                branch.molinfo.push_back(current_mol);
                branch.atmnm.push_back(current_atmnm);
                current_mol.clear();
                current_atmnm.clear();
            }
        }
    }

    int nseg = branch.molinfo.size();
    for (int i = 0; i < nseg; ++i) {
        int natm = branch.molinfo[i].size();
        cout << "Segname (Heavy atom): " 
             << segnm << " " << i+1 << endl;
        for (int j = 0; j < natm; ++j) {
            cout << branch.atmnm[i][j] << " " 
                 << branch.molinfo[i][j] << " ";
        }
        cout << endl;
    }

    return branch;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_reschain PSFOPR::SplitListResid(const s_psf& psf, 
                                  const vector<int>& list) {
    //
    s_reschain reschain;
    //
    int nlist = list.size();
    vector<int> cellinfo;

    for (int i = 0; i < nlist; ++i) {
        int atm, atm_next;
        int resid, resid_next;
        //
        if (i != nlist-1) {
            atm        = list[i];
            atm_next   = list[i+1];
            resid      = psf.molid[atm-1];
            resid_next = psf.molid[atm_next-1];
            cellinfo.push_back(atm);
            if (resid != resid_next) {
                reschain.residinfo.push_back(resid);
                reschain.resinfo.push_back(cellinfo);
                cellinfo.clear();
            }
        } else {
            atm   = list[i];
            resid = psf.molid[atm-1];
            cellinfo.push_back(atm);
            reschain.resinfo.push_back(cellinfo);
            reschain.residinfo.push_back(resid);
        }
    }

    return reschain;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void PSFOPR::LoadpsffileBond(s_psf& psf, 
                             const string& psffile) {
    //
    ifstream file(psffile);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << psffile << endl;
        return;
    }

    string line;
    bool in_bond_section = false;

    while (getline(file, line)) {
        if (line.find("!NBOND") != string::npos) {
            in_bond_section = true;
            continue;
        }

        if (in_bond_section) {
            if (line.empty()) break;
            stringstream ss(line);
            int atom1, atom2;
            while (ss >> atom1 >> atom2) {
                psf.bond.push_back({atom1, atom2});
            }
        }
    }
    //
    file.close();
}
//---------------------------------------------------------------------------//
