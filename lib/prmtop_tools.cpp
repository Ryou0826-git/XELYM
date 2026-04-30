// prmtop_tools.cpp

#include <prmtop_tools.hpp>

#include <cstdlib>
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iomanip>

using namespace std;

s_prmtop LoadfromprmtopKey(const string& prmtopfile, 
                           const string& keyword, 
                           s_prmtop& prmtop) {
    //
    ifstream file(prmtopfile);
    string   line;
    //
    bool inMassSection       = false;
    bool inChargeSection = false;
    bool inAljSection    = false;
    bool inBljSection    = false;

    //
    if (file.is_open()) {
        while (getline(file, line)) {
            istringstream iss(line);
            string token1, token2;
            iss >> token1 >> token2;

            if (token1 == "%FLAG" && token2 == keyword) {
                if (keyword == "MASS") {
                  inMassSection = true;
                } else if (keyword == "CHARGE") {
                  inChargeSection = true;
                } else if (keyword == "LENNARD_JONES_ACOEF") {
                   inAljSection = true;
                } else if (keyword == "LENNARD_JONES_BCOEF") {
                   inBljSection = true;
                } 
                continue;
            }

            if (inMassSection) {
                if (line.find("%FORMAT") != std::string::npos) {
                    continue;
                }

                if (line.find("%FLAG") != std::string::npos) {
                    break;
                }

                istringstream massStream(line);
                double mass;
                while (massStream >> mass) {
                    prmtop.mass.push_back(mass);
                }
            }
            
            if (inChargeSection) {
                if (line.find("%FORMAT") != std::string::npos) {
                    continue;
                }

                if (line.find("%FLAG") != std::string::npos) {
                    break;
                }

                istringstream chargeStream(line);
                double charge;
                while (chargeStream >> charge) {
                    charge /= 18.2223;
                    prmtop.charge.push_back(charge);
                }
            }

            if (inAljSection) {
                if (line.find("%FORMAT") != std::string::npos) {
                    continue;
                }

                if (line.find("%FLAG") != std::string::npos) {
                    break;
                }

                istringstream AljStream(line);
                double Alj;
                while (AljStream >> Alj) {
                    prmtop.Alj.push_back(Alj);
                }
            }
            
            if (inBljSection) {
                if (line.find("%FORMAT") != std::string::npos) {
                    continue;
                }

                if (line.find("%FLAG") != std::string::npos) {
                    break;
                }

                istringstream BljStream(line);
                double Blj;
                while (BljStream >> Blj) {
                    prmtop.Blj.push_back(Blj);
                }
            }
        }
    } else {
        cout << "Error: connot open " 
             << prmtopfile << ", please try again..." 
             << endl;
    }

    return prmtop;
}

s_prmtop PRMTOPOPR::Loadfromprmtop(const string& prmtopfile) {
    //
    s_prmtop prmtop;
    prmtop = LoadfromprmtopKey(prmtopfile, "MASS", prmtop);
    prmtop = LoadfromprmtopKey(prmtopfile, "CHARGE", prmtop);
    prmtop = LoadfromprmtopKey(prmtopfile, "LENNARD_JONES_ACOEF", prmtop);
    prmtop = LoadfromprmtopKey(prmtopfile, "LENNARD_JONES_BCOEF", prmtop);
    //
    return prmtop;
}

//s_prmtop PRMTOPOPR::Loadfromprmtop(const string& prmtopfile) {
//    //
//    ifstream file(prmtopfile);
//    string   line;
//    s_prmtop prmtop;
//    //
//    bool inMassSection   = false;
//    bool inChargeSection = false;
//    bool inAljSection    = false;
//    bool inBljSection    = false;
//
//    //
//    if (file.is_open()) {
//        while (getline(file, line)) {
//            istringstream iss(line);
//            string token1, token2;
//            iss >> token1 >> token2;
//
//            if (token1 == "%FLAG" && token2 == "MASS") {
//                inMassSection = true;
//                continue;
//            }
//
//            if (inMassSection) {
//                if (line.find("%FORMAT") != std::string::npos) {
//                    continue;
//                }
//
//                if (line.find("%FLAG") != std::string::npos) {
//                    break;
//                }
//
//                istringstream massStream(line);
//                double mass;
//                while (massStream >> mass) {
//                    prmtop.mass.push_back(mass);
//                }
//            }
//        }
//    } else {
//        cout << "Error: connot open " 
//             << prmtopfile << ", please try again..." 
//             << endl;
//    }
//
//    return prmtop;
//}
