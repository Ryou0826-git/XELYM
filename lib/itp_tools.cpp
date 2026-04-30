// itp_tools.cpp

#include <itp_tools.hpp>

//-----------------------------------------------------------------------------
// Functions of itp
//-----------------------------------------------------------------------------

void Readitpfile(const string &filename, s_itp &itp) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Error: cannot open " << filename << endl;
        return;
    }

    string line;
    Section current = NONE;

    while (getline(fin, line)) {
        line = trim(line);
        if (line.empty() || line[0] == ';' || line[0] == '#') continue;

        if (line[0] == '[') {
            if      (line.find("atoms")      != string::npos) current = ATOMS;
            else if (line.find("bonds")      != string::npos) current = BONDS;
            else if (line.find("pairs")      != string::npos) current = PAIRS;
            else if (line.find("angles")     != string::npos) current = ANGLES;
            else if (line.find("dihedrals")  != string::npos) current = DIHEDRALS;
            else if (line.find("exclusions") != string::npos) current = EXCLUSIONS;
            else current = NONE;
            continue;
        }

        istringstream iss(line);
        switch (current) {
        case ATOMS: {
            string type, resid, atom;
            int nr, resnr, cgnr;
            double charge, mass;
            if (iss >> nr >> type >> resnr >> resid >> atom >> cgnr >> charge) {
                itp.atoms.nr.push_back(nr);
                itp.atoms.type.push_back(type);
                itp.atoms.resnr.push_back(resnr);
                itp.atoms.resid.push_back(resid);
                itp.atoms.atom.push_back(atom);
                itp.atoms.cgnr.push_back(cgnr);
                itp.atoms.charge.push_back(charge);
                
                if (iss >> mass) {
                    itp.atoms.mass.push_back(mass);
                } else {
                    itp.atoms.mass.push_back(999);
                }
            }
            break;
        }
        case BONDS: {
            int ai, aj, funct;
            double c0 = 0.0, c1 = 0.0;
            if (iss >> ai >> aj >> funct) {
                itp.bonds.ai.push_back(ai);
                itp.bonds.aj.push_back(aj);
                itp.bonds.funct.push_back(funct);
                if (iss >> c0 >> c1) {
                    itp.bonds.c0.push_back(c0);
                    itp.bonds.c1.push_back(c1);
                } else {
                    itp.bonds.c0.push_back(0.0);
                    itp.bonds.c1.push_back(0.0);
                }
            }
            break;
        }
        case PAIRS: {
            int ai, aj, funct;
            if (iss >> ai >> aj >> funct) {
                itp.pairs.ai.push_back(ai);
                itp.pairs.aj.push_back(aj);
                itp.pairs.funct.push_back(funct);
            }
            break;
        }
        case ANGLES: {
            int ai, aj, ak, funct;
            double angle = 0.0, fc = 0.0;
            if (iss >> ai >> aj >> ak >> funct) {
                itp.angles.ai.push_back(ai);
                itp.angles.aj.push_back(aj);
                itp.angles.ak.push_back(ak);
                itp.angles.funct.push_back(funct);
                if (iss >> angle >> fc) {
                    itp.angles.angle.push_back(angle);
                    itp.angles.fc.push_back(fc);
                } else {
                    itp.angles.angle.push_back(0.0);
                    itp.angles.fc.push_back(0.0);
                }
            }
            break;
        }
        case DIHEDRALS: {
            int ai, aj, ak, al, funct, mult;
            double ph0 = 0.0, cp = 0.0;
            double angle, fc;
            double phase, kd;
            int    pn;
            double c0, c1, c2, c3, c4, c5;
            if (iss >> ai >> aj >> ak >> al >> funct) {
                itp.dihedrals.ai.push_back(ai);
                itp.dihedrals.aj.push_back(aj);
                itp.dihedrals.ak.push_back(ak);
                itp.dihedrals.al.push_back(al);
                itp.dihedrals.funct.push_back(funct);
                if (funct == 1) {
                    iss >> ph0 >> cp >> mult;
                    itp.dihedrals.ph0.push_back(ph0);
                    itp.dihedrals.cp.push_back(cp);
                    itp.dihedrals.mult.push_back(mult);
                    itp.dihedrals.angle.push_back(0.0);
                    itp.dihedrals.fc.push_back(0.0);
                    itp.dihedrals.phase.push_back(0.0);
                    itp.dihedrals.kd.push_back(0.0);
                    itp.dihedrals.pn.push_back(0.0);
                    itp.dihedrals.c0.push_back(0.0);
                    itp.dihedrals.c1.push_back(0.0);
                    itp.dihedrals.c2.push_back(0.0);
                    itp.dihedrals.c3.push_back(0.0);
                    itp.dihedrals.c4.push_back(0.0);
                    itp.dihedrals.c5.push_back(0.0);
                } else if (funct == 2) {
                    iss >> angle >> fc;
                    itp.dihedrals.angle.push_back(angle);
                    itp.dihedrals.fc.push_back(fc);
                    itp.dihedrals.ph0.push_back(0.0);
                    itp.dihedrals.cp.push_back(0.0);
                    itp.dihedrals.mult.push_back(0.0);
                    itp.dihedrals.phase.push_back(0.0);
                    itp.dihedrals.kd.push_back(0.0);
                    itp.dihedrals.pn.push_back(0.0);
                    itp.dihedrals.c0.push_back(0.0);
                    itp.dihedrals.c1.push_back(0.0);
                    itp.dihedrals.c2.push_back(0.0);
                    itp.dihedrals.c3.push_back(0.0);
                    itp.dihedrals.c4.push_back(0.0);
                    itp.dihedrals.c5.push_back(0.0);
                } else if ((funct == 9) || (funct == 4)) {
                    iss >> phase >> kd >> pn;
                    itp.dihedrals.phase.push_back(phase);
                    itp.dihedrals.kd.push_back(kd);
                    itp.dihedrals.pn.push_back(pn);
                    itp.dihedrals.ph0.push_back(0.0);
                    itp.dihedrals.cp.push_back(0.0);
                    itp.dihedrals.mult.push_back(0.0);
                    itp.dihedrals.angle.push_back(0.0);
                    itp.dihedrals.fc.push_back(0.0);
                    itp.dihedrals.c0.push_back(0.0);
                    itp.dihedrals.c1.push_back(0.0);
                    itp.dihedrals.c2.push_back(0.0);
                    itp.dihedrals.c3.push_back(0.0);
                    itp.dihedrals.c4.push_back(0.0);
                    itp.dihedrals.c5.push_back(0.0);
                } else if (funct == 3) {
                    if (iss >> c0 >> c1 >> c2 >> c3 >> c4 >> c5) {
                        itp.dihedrals.c0.push_back(c0);
                        itp.dihedrals.c1.push_back(c1);
                        itp.dihedrals.c2.push_back(c2);
                        itp.dihedrals.c3.push_back(c3);
                        itp.dihedrals.c4.push_back(c4);
                        itp.dihedrals.c5.push_back(c5);
                    } else {
                        itp.dihedrals.c0.push_back(999);
                        itp.dihedrals.c1.push_back(999);
                        itp.dihedrals.c2.push_back(999);
                        itp.dihedrals.c3.push_back(999);
                        itp.dihedrals.c4.push_back(999);
                        itp.dihedrals.c5.push_back(999);
                    }
                    itp.dihedrals.phase.push_back(0.0);
                    itp.dihedrals.kd.push_back(0.0);
                    itp.dihedrals.pn.push_back(0.0);
                    itp.dihedrals.ph0.push_back(0.0);
                    itp.dihedrals.cp.push_back(0.0);
                    itp.dihedrals.mult.push_back(0.0);
                    itp.dihedrals.angle.push_back(0.0);
                    itp.dihedrals.fc.push_back(0.0);
                }
            }
            break;
        }
        case EXCLUSIONS: {
            int ai, aj;
            if (iss >> ai) {
                while (iss >> aj) {
                    itp.exclusion.ai.push_back(ai);
                    itp.exclusion.aj.push_back(aj);
                }
            }
            break;
        }
        default: break;
        }
    }
    itp.atoms.n = itp.atoms.atom.size();
}


void RemoveIndex(s_itp &itp, int target) {
    // --- [1] atoms ---
    auto erase_atom_entry = [&](auto &v) {
        if (target - 1 < (int)v.size())
            v.erase(v.begin() + (target - 1));
    };
    erase_atom_entry(itp.atoms.nr);
    erase_atom_entry(itp.atoms.type);
    erase_atom_entry(itp.atoms.resnr);
    erase_atom_entry(itp.atoms.resid);
    erase_atom_entry(itp.atoms.atom);
    erase_atom_entry(itp.atoms.cgnr);
    erase_atom_entry(itp.atoms.charge);
    erase_atom_entry(itp.atoms.mass);

    // --- [2] bonds ---
    auto erase_related_bonds = [&](int idx){
        for (int i = (int)itp.bonds.ai.size() - 1; i >= 0; --i) {
            if (itp.bonds.ai[i] == idx || itp.bonds.aj[i] == idx) {
                itp.bonds.ai.erase(itp.bonds.ai.begin() + i);
                itp.bonds.aj.erase(itp.bonds.aj.begin() + i);
                itp.bonds.funct.erase(itp.bonds.funct.begin() + i);
                itp.bonds.c0.erase(itp.bonds.c0.begin() + i);
                itp.bonds.c1.erase(itp.bonds.c1.begin() + i);
            }
        }
    };
    erase_related_bonds(target);

    // --- [3] pairs ---
    for (int i = (int)itp.pairs.ai.size() - 1; i >= 0; --i) {
        if (itp.pairs.ai[i] == target || itp.pairs.aj[i] == target) {
            itp.pairs.ai.erase(itp.pairs.ai.begin() + i);
            itp.pairs.aj.erase(itp.pairs.aj.begin() + i);
            itp.pairs.funct.erase(itp.pairs.funct.begin() + i);
        }
    }

    // --- [4] angles ---
    for (int i = (int)itp.angles.ai.size() - 1; i >= 0; --i) {
        if (itp.angles.ai[i] == target || itp.angles.aj[i] == target || itp.angles.ak[i] == target) {
            itp.angles.ai.erase(itp.angles.ai.begin() + i);
            itp.angles.aj.erase(itp.angles.aj.begin() + i);
            itp.angles.ak.erase(itp.angles.ak.begin() + i);
            itp.angles.funct.erase(itp.angles.funct.begin() + i);
            itp.angles.angle.erase(itp.angles.angle.begin() + i);
            itp.angles.fc.erase(itp.angles.fc.begin() + i);
        }
    }

    // --- [5] dihedrals ---
    for (int i = (int)itp.dihedrals.ai.size() - 1; i >= 0; --i) {
        if (itp.dihedrals.ai[i] == target || itp.dihedrals.aj[i] == target ||
            itp.dihedrals.ak[i] == target || itp.dihedrals.al[i] == target) {
            itp.dihedrals.ai.erase(itp.dihedrals.ai.begin() + i);
            itp.dihedrals.aj.erase(itp.dihedrals.aj.begin() + i);
            itp.dihedrals.ak.erase(itp.dihedrals.ak.begin() + i);
            itp.dihedrals.al.erase(itp.dihedrals.al.begin() + i);
            itp.dihedrals.funct.erase(itp.dihedrals.funct.begin() + i);
            itp.dihedrals.ph0.erase(itp.dihedrals.ph0.begin() + i);
            itp.dihedrals.cp.erase(itp.dihedrals.cp.begin() + i);
            itp.dihedrals.mult.erase(itp.dihedrals.mult.begin() + i);
            itp.dihedrals.angle.erase(itp.dihedrals.angle.begin() + i);
            itp.dihedrals.fc.erase(itp.dihedrals.fc.begin() + i);
            itp.dihedrals.phase.erase(itp.dihedrals.phase.begin() + i);
            itp.dihedrals.kd.erase(itp.dihedrals.kd.begin() + i);
            itp.dihedrals.pn.erase(itp.dihedrals.pn.begin() + i);
            itp.dihedrals.c0.erase(itp.dihedrals.c0.begin() + i);
            itp.dihedrals.c1.erase(itp.dihedrals.c1.begin() + i);
            itp.dihedrals.c2.erase(itp.dihedrals.c2.begin() + i);
            itp.dihedrals.c3.erase(itp.dihedrals.c3.begin() + i);
            itp.dihedrals.c4.erase(itp.dihedrals.c4.begin() + i);
            itp.dihedrals.c5.erase(itp.dihedrals.c5.begin() + i);
        }
    }

    // --- [6] exclusion ---
    for (int i = (int)itp.exclusion.ai.size() - 1; i >= 0; --i) {
        if (itp.exclusion.ai[i] == target || itp.exclusion.aj[i] == target) {
            itp.exclusion.ai.erase(itp.exclusion.ai.begin() + i);
            itp.exclusion.aj.erase(itp.exclusion.aj.begin() + i);
        }
    }

    auto renumber = [&](int &x) { if (x > target) x--; };

    for (auto &x : itp.atoms.cgnr) renumber(x);
    for (auto &x : itp.atoms.nr) renumber(x);
    for (auto &x : itp.bonds.ai) renumber(x);
    for (auto &x : itp.bonds.aj) renumber(x);
    for (auto &x : itp.pairs.ai) renumber(x);
    for (auto &x : itp.pairs.aj) renumber(x);
    for (auto &x : itp.angles.ai) renumber(x);
    for (auto &x : itp.angles.aj) renumber(x);
    for (auto &x : itp.angles.ak) renumber(x);
    for (auto &x : itp.dihedrals.ai) renumber(x);
    for (auto &x : itp.dihedrals.aj) renumber(x);
    for (auto &x : itp.dihedrals.ak) renumber(x);
    for (auto &x : itp.dihedrals.al) renumber(x);
    for (auto &x : itp.exclusion.ai) renumber(x);
    for (auto &x : itp.exclusion.aj) renumber(x);
}

void RemoveAtomsKeepIndex(s_itp& itp, const vector<int>& atoms_to_remove) {
    unordered_set<int> remove_set(atoms_to_remove.begin(), atoms_to_remove.end());

    // --- atoms ---
    s_atoms new_atoms;
    for (size_t i = 0; i < itp.atoms.nr.size(); ++i) {
        int nr = itp.atoms.nr[i];
        if (remove_set.count(nr) == 0) {
            new_atoms.nr.push_back(nr);
            new_atoms.type.push_back(itp.atoms.type[i]);
            new_atoms.resnr.push_back(itp.atoms.resnr[i]);
            new_atoms.resid.push_back(itp.atoms.resid[i]);
            new_atoms.atom.push_back(itp.atoms.atom[i]);
            new_atoms.cgnr.push_back(itp.atoms.cgnr[i]);
            new_atoms.charge.push_back(itp.atoms.charge[i]);
            new_atoms.mass.push_back(itp.atoms.mass[i]);
        }
    }
    itp.atoms = move(new_atoms);

    // --- bonds ---
    s_bonds new_b;
    for (size_t i = 0; i < itp.bonds.ai.size(); ++i) {
        int ai = itp.bonds.ai[i];
        int aj = itp.bonds.aj[i];
        if (remove_set.count(ai) == 0 && remove_set.count(aj) == 0) {
            new_b.ai.push_back(ai);
            new_b.aj.push_back(aj);
            new_b.funct.push_back(itp.bonds.funct[i]);
            if (i < itp.bonds.c0.size()) new_b.c0.push_back(itp.bonds.c0[i]);
            if (i < itp.bonds.c1.size()) new_b.c1.push_back(itp.bonds.c1[i]);
        }
    }
    itp.bonds = move(new_b);

    // --- pairs ---
    s_pairs new_p;
    for (size_t i = 0; i < itp.pairs.ai.size(); ++i) {
        int ai = itp.pairs.ai[i];
        int aj = itp.pairs.aj[i];
        if (remove_set.count(ai) == 0 && remove_set.count(aj) == 0) {
            new_p.ai.push_back(ai);
            new_p.aj.push_back(aj);
            new_p.funct.push_back(itp.pairs.funct[i]);
        }
    }
    itp.pairs = move(new_p);

    // --- angles ---
    s_angles new_a;
    for (size_t i = 0; i < itp.angles.ai.size(); ++i) {
        int ai = itp.angles.ai[i];
        int aj = itp.angles.aj[i];
        int ak = itp.angles.ak[i];
        if (remove_set.count(ai) == 0 && remove_set.count(aj) == 0 && remove_set.count(ak) == 0) {
            new_a.ai.push_back(ai);
            new_a.aj.push_back(aj);
            new_a.ak.push_back(ak);
            new_a.funct.push_back(itp.angles.funct[i]);
            if (i < itp.angles.angle.size()) new_a.angle.push_back(itp.angles.angle[i]);
            if (i < itp.angles.fc.size()) new_a.fc.push_back(itp.angles.fc[i]);
        }
    }
    itp.angles = move(new_a);

    // --- dihedrals ---
    s_dihedrals new_d;
    for (size_t i = 0; i < itp.dihedrals.ai.size(); ++i) {
        int ai = itp.dihedrals.ai[i];
        int aj = itp.dihedrals.aj[i];
        int ak = itp.dihedrals.ak[i];
        int al = itp.dihedrals.al[i];
        if (remove_set.count(ai) == 0 && remove_set.count(aj) == 0 &&
            remove_set.count(ak) == 0 && remove_set.count(al) == 0) {
            new_d.ai.push_back(ai);
            new_d.aj.push_back(aj);
            new_d.ak.push_back(ak);
            new_d.al.push_back(al);
            new_d.funct.push_back(itp.dihedrals.funct[i]);
            //
            new_d.ph0.push_back(itp.dihedrals.ph0[i]);
            new_d.cp.push_back(itp.dihedrals.cp[i]);
            new_d.mult.push_back(itp.dihedrals.mult[i]);
            new_d.angle.push_back(itp.dihedrals.angle[i]);
            new_d.fc.push_back(itp.dihedrals.fc[i]);
            
            new_d.phase.push_back(itp.dihedrals.phase[i]);
            new_d.kd.push_back(itp.dihedrals.kd[i]);
            new_d.pn.push_back(itp.dihedrals.pn[i]);
            //
            new_d.c0.push_back(itp.dihedrals.c0[i]);
            new_d.c1.push_back(itp.dihedrals.c1[i]);
            new_d.c2.push_back(itp.dihedrals.c2[i]);
            new_d.c3.push_back(itp.dihedrals.c3[i]);
            new_d.c4.push_back(itp.dihedrals.c4[i]);
            new_d.c5.push_back(itp.dihedrals.c5[i]);
            
            //if ( (itp.dihedrals.funct[i] == 4) 
            //      || (itp.dihedrals.funct[i] == 9) ) {
            //    new_d.phase.push_back(itp.dihedrals.phase[i]);
            //    new_d.kd.push_back(itp.dihedrals.kd[i]);
            //    new_d.pn.push_back(itp.dihedrals.pn[i]);
            //}
        }
    }
    itp.dihedrals = move(new_d);

    // --- exclusions ---
    s_exclusion new_e;
    for (size_t i = 0; i < itp.exclusion.ai.size(); ++i) {
        int ai = itp.exclusion.ai[i];
        int aj = itp.exclusion.aj[i];
        if (remove_set.count(ai) == 0 && remove_set.count(aj) == 0) {
            new_e.ai.push_back(ai);
            new_e.aj.push_back(aj);
        }
    }
    itp.exclusion = move(new_e);
}


void RemoveConnections(s_itp& itp, int target) {
    // --- [2] bonds ---
    auto erase_related_bonds = [&](int idx){
        for (int i = (int)itp.bonds.ai.size() - 1; i >= 0; --i) {
            if (itp.bonds.ai[i] == idx || itp.bonds.aj[i] == idx) {
                itp.bonds.ai.erase(itp.bonds.ai.begin() + i);
                itp.bonds.aj.erase(itp.bonds.aj.begin() + i);
                itp.bonds.funct.erase(itp.bonds.funct.begin() + i);
                itp.bonds.c0.erase(itp.bonds.c0.begin() + i);
                itp.bonds.c1.erase(itp.bonds.c1.begin() + i);
            }
        }
    };
    erase_related_bonds(target);

    // --- [3] pairs ---
    for (int i = (int)itp.pairs.ai.size() - 1; i >= 0; --i) {
        if (itp.pairs.ai[i] == target || itp.pairs.aj[i] == target) {
            itp.pairs.ai.erase(itp.pairs.ai.begin() + i);
            itp.pairs.aj.erase(itp.pairs.aj.begin() + i);
            itp.pairs.funct.erase(itp.pairs.funct.begin() + i);
        }
    }

    // --- [4] angles ---
    for (int i = (int)itp.angles.ai.size() - 1; i >= 0; --i) {
        if (itp.angles.ai[i] == target || itp.angles.aj[i] == target || itp.angles.ak[i] == target) {
            itp.angles.ai.erase(itp.angles.ai.begin() + i);
            itp.angles.aj.erase(itp.angles.aj.begin() + i);
            itp.angles.ak.erase(itp.angles.ak.begin() + i);
            itp.angles.funct.erase(itp.angles.funct.begin() + i);
            itp.angles.angle.erase(itp.angles.angle.begin() + i);
            itp.angles.fc.erase(itp.angles.fc.begin() + i);
        }
    }

    // --- [5] dihedrals ---
    for (int i = (int)itp.dihedrals.ai.size() - 1; i >= 0; --i) {
        if (itp.dihedrals.ai[i] == target || itp.dihedrals.aj[i] == target ||
             itp.dihedrals.ak[i] == target || itp.dihedrals.al[i] == target) {
            itp.dihedrals.ai.erase(itp.dihedrals.ai.begin() + i);
            itp.dihedrals.aj.erase(itp.dihedrals.aj.begin() + i);
            itp.dihedrals.ak.erase(itp.dihedrals.ak.begin() + i);
            itp.dihedrals.al.erase(itp.dihedrals.al.begin() + i);
            itp.dihedrals.funct.erase(itp.dihedrals.funct.begin() + i);
            itp.dihedrals.ph0.erase(itp.dihedrals.ph0.begin() + i);
            itp.dihedrals.cp.erase(itp.dihedrals.cp.begin() + i);
            itp.dihedrals.mult.erase(itp.dihedrals.mult.begin() + i);
            itp.dihedrals.angle.erase(itp.dihedrals.angle.begin() + i);
            itp.dihedrals.fc.erase(itp.dihedrals.fc.begin() + i);
            itp.dihedrals.phase.erase(itp.dihedrals.phase.begin() + i);
            itp.dihedrals.kd.erase(itp.dihedrals.kd.begin() + i);
            itp.dihedrals.pn.erase(itp.dihedrals.pn.begin() + i);
            itp.dihedrals.c0.erase(itp.dihedrals.c0.begin() + i);
            itp.dihedrals.c1.erase(itp.dihedrals.c1.begin() + i);
            itp.dihedrals.c2.erase(itp.dihedrals.c2.begin() + i);
            itp.dihedrals.c3.erase(itp.dihedrals.c3.begin() + i);
            itp.dihedrals.c4.erase(itp.dihedrals.c4.begin() + i);
            itp.dihedrals.c5.erase(itp.dihedrals.c5.begin() + i);
            //if ((itp.dihedrals.funct[i] == 1) || (itp.dihedrals.funct[i] == 2)) {
            //    itp.dihedrals.ai.erase(itp.dihedrals.ai.begin() + i);
            //    itp.dihedrals.aj.erase(itp.dihedrals.aj.begin() + i);
            //    itp.dihedrals.ak.erase(itp.dihedrals.ak.begin() + i);
            //    itp.dihedrals.al.erase(itp.dihedrals.al.begin() + i);
            //    itp.dihedrals.funct.erase(itp.dihedrals.funct.begin() + i);
            //    itp.dihedrals.ph0.erase(itp.dihedrals.ph0.begin() + i);
            //    itp.dihedrals.cp.erase(itp.dihedrals.cp.begin() + i);
            //    itp.dihedrals.mult.erase(itp.dihedrals.mult.begin() + i);
            //    itp.dihedrals.angle.erase(itp.dihedrals.angle.begin() + i);
            //    itp.dihedrals.fc.erase(itp.dihedrals.fc.begin() + i);
            //} else if ((itp.dihedrals.funct[i] == 4) || (itp.dihedrals.funct[i] == 9)) {
            //    itp.dihedrals.ai.erase(itp.dihedrals.ai.begin() + i);
            //    itp.dihedrals.aj.erase(itp.dihedrals.aj.begin() + i);
            //    itp.dihedrals.ak.erase(itp.dihedrals.ak.begin() + i);
            //    itp.dihedrals.al.erase(itp.dihedrals.al.begin() + i);
            //    itp.dihedrals.funct.erase(itp.dihedrals.funct.begin() + i);
            //    itp.dihedrals.phase.erase(itp.dihedrals.phase.begin() + i);
            //    itp.dihedrals.kd.erase(itp.dihedrals.kd.begin() + i);
            //    itp.dihedrals.pn.erase(itp.dihedrals.pn.begin() + i);
            //}
        }
    }

    // --- [6] exclusions ---
    for (int i = (int)itp.exclusion.ai.size() - 1; i >= 0; --i) {
        if (itp.exclusion.ai[i] == target || itp.exclusion.aj[i] == target) {
            itp.exclusion.ai.erase(itp.exclusion.ai.begin() + i);
            itp.exclusion.aj.erase(itp.exclusion.aj.begin() + i);
        }
    }
}

void CompactAtomIndices(s_itp &itp) {
    vector<int> old_indices;
    for (int nr : itp.atoms.nr) {
        old_indices.push_back(nr);
    }

    unordered_map<int, int> remap;
    for (size_t i = 0; i < old_indices.size(); ++i) {
        remap[old_indices[i]] = static_cast<int>(i + 1);
    }

    for (size_t i = 0; i < itp.atoms.nr.size(); ++i) {
        itp.atoms.nr[i] = static_cast<int>(i + 1);
    }

    for (auto &b : itp.bonds.ai) b = remap[b];
    for (auto &b : itp.bonds.aj) b = remap[b];

    for (auto &a : itp.angles.ai) a = remap[a];
    for (auto &a : itp.angles.aj) a = remap[a];
    for (auto &a : itp.angles.ak) a = remap[a];

    for (auto &d : itp.dihedrals.ai) d = remap[d];
    for (auto &d : itp.dihedrals.aj) d = remap[d];
    for (auto &d : itp.dihedrals.ak) d = remap[d];
    for (auto &d : itp.dihedrals.al) d = remap[d];

    for (auto &p : itp.pairs.ai) p = remap[p];
    for (auto &p : itp.pairs.aj) p = remap[p];

    for (auto &e : itp.exclusion.ai) e = remap[e];
    for (auto &e : itp.exclusion.aj) e = remap[e];

    cout << "CompactAtomIndices: atom indices renumbered (1-based)." << endl;
}

void Writeitp(const s_itp &itp, const string &filename) {
    ofstream ofs(filename);
    if (!ofs.is_open()) {
        cerr << "Error: cannot open file " << filename << endl;
        return;
    }

    ofs << "; generated by custom writer\n\n";

    //----------------------------
    // [ moleculetype ]
    //----------------------------
    ofs << "[ moleculetype ]\n";
    ofs << "MOL 3\n\n";
    
    // ---------------------------
    // [ atoms ]
    // ---------------------------
    ofs << "[ atoms ]\n";
    ofs << "; nr  type  resnr  resid  atom  cgnr  charge  mass\n";
    for (size_t i = 0; i < itp.atoms.nr.size(); ++i) {
        ofs << setw(7)  << itp.atoms.nr[i]
            << setw(10)  << itp.atoms.type[i]
            << setw(8)  << itp.atoms.resnr[i]
            << setw(8)  << itp.atoms.resid[i]
            << setw(8)  << itp.atoms.atom[i]
            << setw(8)  << itp.atoms.cgnr[i]
            << setw(12) << fixed << setprecision(7) << itp.atoms.charge[i];
            if (itp.atoms.mass[i] < 998) {
                ofs << setw(12) << fixed << setprecision(5) << itp.atoms.mass[i];
            }
            ofs << "\n";
    }
    ofs << "\n";

    // ---------------------------
    // [ bonds ]
    // ---------------------------
    ofs << "[ bonds ]\n";
    ofs << "; ai  aj  funct  c0  c1\n";
    for (size_t i = 0; i < itp.bonds.ai.size(); ++i) {
        ofs << setw(7) << itp.bonds.ai[i]
            << setw(7) << itp.bonds.aj[i]
            << setw(7) << itp.bonds.funct[i];
        if (!itp.bonds.c0.empty())
            ofs << setw(12) << fixed << setprecision(5) << itp.bonds.c0[i];
        if (!itp.bonds.c1.empty())
            ofs << setw(13) << scientific << setprecision(5) << itp.bonds.c1[i];
        ofs << "\n";
    }
    ofs << "\n";

    // ---------------------------
    // [ pairs ]
    // ---------------------------
    ofs << "[ pairs ]\n";
    ofs << "; ai  aj  funct\n";
    for (size_t i = 0; i < itp.pairs.ai.size(); ++i) {
        ofs << setw(7) << itp.pairs.ai[i]
            << setw(7) << itp.pairs.aj[i]
            << setw(7) << itp.pairs.funct[i]
            << "\n";
    }
    ofs << "\n";

    // ---------------------------
    // [ angles ]
    // ---------------------------
    ofs << "[ angles ]\n";
    ofs << "; ai  aj  ak  funct  angle  fc\n";
    for (size_t i = 0; i < itp.angles.ai.size(); ++i) {
        ofs << setw(7) << itp.angles.ai[i]
            << setw(7) << itp.angles.aj[i]
            << setw(7) << itp.angles.ak[i]
            << setw(7) << itp.angles.funct[i];
        if (!itp.angles.angle.empty())
            ofs << setw(12) << fixed << setprecision(3) << itp.angles.angle[i];
        if (!itp.angles.fc.empty())
            ofs << setw(12) << fixed << setprecision(3) << itp.angles.fc[i];
        ofs << "\n";
    }
    ofs << "\n";

    // ---------------------------
    // [ dihedrals ]
    // ---------------------------
    ofs << "[ dihedrals ]\n";
    ofs << "; ai  aj  ak  al  funct  ph0  cp  mult\n";
    for (size_t i = 0; i < itp.dihedrals.ai.size(); ++i) {
        ofs << setw(7) << itp.dihedrals.ai[i]
            << setw(7) << itp.dihedrals.aj[i]
            << setw(7) << itp.dihedrals.ak[i]
            << setw(7) << itp.dihedrals.al[i]
            << setw(7) << itp.dihedrals.funct[i];
        if (itp.dihedrals.funct[i] == 1) {
            ofs << setw(10) << fixed << setprecision(3) << itp.dihedrals.ph0[i];
            ofs << setw(10) << fixed << setprecision(3) << itp.dihedrals.cp[i];
            ofs << setw(5) << itp.dihedrals.mult[i];
        } else if (itp.dihedrals.funct[i] == 2) {
            ofs << setw(10) << fixed << setprecision(3) << itp.dihedrals.angle[i];
            ofs << setw(10) << fixed << setprecision(3) << itp.dihedrals.fc[i];
        } else if ( (itp.dihedrals.funct[i] == 4) 
                   || (itp.dihedrals.funct[i] == 9) ) {
            //cout << itp.dihedrals.phase[i] 
            //     << " " << itp.dihedrals.kd[i]
            //     << " " << itp.dihedrals.pn[i] << endl;
            ofs << setw(10) << fixed << setprecision(5) << itp.dihedrals.phase[i];
            ofs << setw(10) << fixed << setprecision(5) << itp.dihedrals.kd[i];
            ofs << setw(5)  << itp.dihedrals.pn[i];
        } else if (itp.dihedrals.funct[i] == 3) {
            if ( (itp.dihedrals.c0[i] < 998) 
                  && (itp.dihedrals.c1[i] < 998) 
                  && (itp.dihedrals.c2[i] < 998) 
                  && (itp.dihedrals.c3[i] < 998) 
                  && (itp.dihedrals.c4[i] < 998) 
                  && (itp.dihedrals.c5[i] < 998) ) {
                ofs << setw(10) << fixed << setprecision(5) << itp.dihedrals.c0[i];
                ofs << setw(10) << fixed << setprecision(5) << itp.dihedrals.c1[i];
                ofs << setw(10) << fixed << setprecision(5) << itp.dihedrals.c2[i];
                ofs << setw(10) << fixed << setprecision(5) << itp.dihedrals.c3[i];
                ofs << setw(10) << fixed << setprecision(5) << itp.dihedrals.c4[i];
                ofs << setw(10) << fixed << setprecision(5) << itp.dihedrals.c5[i];
            }
        }
        ofs << "\n";
    }
    ofs << "\n";

    // ---------------------------
    // [ exclusions ]
    // ---------------------------
    ofs << "[ exclusions ]\n";
    ofs << "; ai  aj\n";
    for (size_t i = 0; i < itp.exclusion.ai.size(); ++i) {
        ofs << setw(7) << itp.exclusion.ai[i]
            << setw(7) << itp.exclusion.aj[i]
            << "\n";
    }
    ofs << "\n";

    ofs.close();
}

s_itp CombineItp(const s_itp& A, const s_itp& B) {
    s_itp C = A;
    int offset = A.atoms.nr.size();

    // --- atoms ---
    for (size_t i = 0; i < B.atoms.nr.size(); ++i) {
        C.atoms.nr.push_back(B.atoms.nr[i] + offset);
        C.atoms.type.push_back(B.atoms.type[i]);
        C.atoms.resnr.push_back(B.atoms.resnr[i]);
        C.atoms.resid.push_back(B.atoms.resid[i]);
        C.atoms.atom.push_back(B.atoms.atom[i]);
        C.atoms.cgnr.push_back(B.atoms.cgnr[i] + offset);
        C.atoms.charge.push_back(B.atoms.charge[i]);
        C.atoms.mass.push_back(B.atoms.mass[i]);
    }

    // --- bonds ---
    for (size_t i = 0; i < B.bonds.ai.size(); ++i) {
        C.bonds.ai.push_back(B.bonds.ai[i] + offset);
        C.bonds.aj.push_back(B.bonds.aj[i] + offset);
        C.bonds.funct.push_back(B.bonds.funct[i]);
        C.bonds.c0.push_back(B.bonds.c0[i]);
        C.bonds.c1.push_back(B.bonds.c1[i]);
    }

    // --- pairs ---
    for (size_t i = 0; i < B.pairs.ai.size(); ++i) {
        C.pairs.ai.push_back(B.pairs.ai[i] + offset);
        C.pairs.aj.push_back(B.pairs.aj[i] + offset);
        C.pairs.funct.push_back(B.pairs.funct[i]);
    }

    // --- angles ---
    for (size_t i = 0; i < B.angles.ai.size(); ++i) {
        C.angles.ai.push_back(B.angles.ai[i] + offset);
        C.angles.aj.push_back(B.angles.aj[i] + offset);
        C.angles.ak.push_back(B.angles.ak[i] + offset);
        C.angles.funct.push_back(B.angles.funct[i]);
        C.angles.angle.push_back(B.angles.angle[i]);
        C.angles.fc.push_back(B.angles.fc[i]);
    }

    // --- dihedrals ---
    for (size_t i = 0; i < B.dihedrals.ai.size(); ++i) {
        C.dihedrals.ai.push_back(B.dihedrals.ai[i] + offset);
        C.dihedrals.aj.push_back(B.dihedrals.aj[i] + offset);
        C.dihedrals.ak.push_back(B.dihedrals.ak[i] + offset);
        C.dihedrals.al.push_back(B.dihedrals.al[i] + offset);
        C.dihedrals.funct.push_back(B.dihedrals.funct[i]);
        C.dihedrals.ph0.push_back(B.dihedrals.ph0[i]);
        C.dihedrals.cp.push_back(B.dihedrals.cp[i]);
        C.dihedrals.mult.push_back(B.dihedrals.mult[i]);
        C.dihedrals.angle.push_back(B.dihedrals.angle[i]);
        C.dihedrals.fc.push_back(B.dihedrals.fc[i]);
        C.dihedrals.phase.push_back(B.dihedrals.phase[i]);
        C.dihedrals.kd.push_back(B.dihedrals.kd[i]);
        C.dihedrals.pn.push_back(B.dihedrals.pn[i]);
        C.dihedrals.c0.push_back(B.dihedrals.c0[i]);
        C.dihedrals.c1.push_back(B.dihedrals.c1[i]);
        C.dihedrals.c2.push_back(B.dihedrals.c2[i]);
        C.dihedrals.c3.push_back(B.dihedrals.c3[i]);
        C.dihedrals.c4.push_back(B.dihedrals.c4[i]);
        C.dihedrals.c5.push_back(B.dihedrals.c5[i]);
        //if ((B.dihedrals.funct[i] == 1) || (B.dihedrals.funct[i] == 2)) {
        //    C.dihedrals.ph0.push_back(B.dihedrals.ph0[i]);
        //    C.dihedrals.cp.push_back(B.dihedrals.cp[i]);
        //    C.dihedrals.mult.push_back(B.dihedrals.mult[i]);
        //    C.dihedrals.angle.push_back(B.dihedrals.angle[i]);
        //    C.dihedrals.fc.push_back(B.dihedrals.fc[i]);
        //} else if ((B.dihedrals.funct[i] == 4) || (B.dihedrals.funct[i] == 9)) {
        //    C.dihedrals.phase.push_back(B.dihedrals.phase[i]);
        //    C.dihedrals.kd.push_back(B.dihedrals.kd[i]);
        //    C.dihedrals.pn.push_back(B.dihedrals.pn[i]);
        //} else if (B.dihedrals.funct[i] == 3) {
        //    C.dihedrals.c0.push_back(B.dihedrals.c0[i]);
        //    C.dihedrals.c1.push_back(B.dihedrals.c1[i]);
        //    C.dihedrals.c2.push_back(B.dihedrals.c2[i]);
        //    C.dihedrals.c3.push_back(B.dihedrals.c3[i]);
        //    C.dihedrals.c4.push_back(B.dihedrals.c4[i]);
        //    C.dihedrals.c5.push_back(B.dihedrals.c5[i]);
        //}
    }

    // --- exclusions ---
    for (size_t i = 0; i < B.exclusion.ai.size(); ++i) {
        C.exclusion.ai.push_back(B.exclusion.ai[i] + offset);
        C.exclusion.aj.push_back(B.exclusion.aj[i] + offset);
    }

    return C;
}

void AppendBond(s_itp& itp, 
                const int idx1, 
                const int idx2,
                const int funct,
                const double c0, 
                const double c1) {
    //
    itp.bonds.ai.push_back(idx1);
    itp.bonds.aj.push_back(idx2);
    itp.bonds.funct.push_back(funct);
    itp.bonds.c0.push_back(c0);
    itp.bonds.c1.push_back(c1);
}

void AppendAngle(s_itp& itp, 
                 const int idx1, 
                 const int idx2, 
                 const int idx3, 
                 const int funct, 
                 const double angle, 
                 const double fc) {
    //
    itp.angles.ai.push_back(idx1);
    itp.angles.aj.push_back(idx2);
    itp.angles.ak.push_back(idx3);
    itp.angles.funct.push_back(funct);
    itp.angles.angle.push_back(angle);
    itp.angles.fc.push_back(fc);
}

void AppendDihedral(s_itp& itp, 
                    const int idx1, 
                    const int idx2, 
                    const int idx3, 
                    const int idx4, 
                    const int funct,
                    const s_dihed_parts dihed_parts) {
    //
    itp.dihedrals.ai.push_back(idx1);
    itp.dihedrals.aj.push_back(idx2);
    itp.dihedrals.ak.push_back(idx3);
    itp.dihedrals.al.push_back(idx4);
    itp.dihedrals.funct.push_back(funct);
    if (funct == 1) {
        itp.dihedrals.ph0.push_back(dihed_parts.ph0);
        itp.dihedrals.cp.push_back(dihed_parts.cp);
        itp.dihedrals.mult.push_back(dihed_parts.mult);
        itp.dihedrals.angle.push_back(0.0);
        itp.dihedrals.fc.push_back(0.0);
        itp.dihedrals.phase.push_back(0.0);
        itp.dihedrals.kd.push_back(0.0);
        itp.dihedrals.pn.push_back(0.0);
        itp.dihedrals.c0.push_back(0.0);
        itp.dihedrals.c1.push_back(0.0);
        itp.dihedrals.c2.push_back(0.0);
        itp.dihedrals.c3.push_back(0.0);
        itp.dihedrals.c4.push_back(0.0);
        itp.dihedrals.c5.push_back(0.0);
    } else if (funct == 2) {
        itp.dihedrals.angle.push_back(dihed_parts.angle);
        itp.dihedrals.fc.push_back(dihed_parts.fc);
        itp.dihedrals.ph0.push_back(0.0);
        itp.dihedrals.cp.push_back(0.0);
        itp.dihedrals.mult.push_back(0.0);
        itp.dihedrals.phase.push_back(0.0);
        itp.dihedrals.kd.push_back(0.0);
        itp.dihedrals.pn.push_back(0.0);
        itp.dihedrals.c0.push_back(0.0);
        itp.dihedrals.c1.push_back(0.0);
        itp.dihedrals.c2.push_back(0.0);
        itp.dihedrals.c3.push_back(0.0);
        itp.dihedrals.c4.push_back(0.0);
        itp.dihedrals.c5.push_back(0.0);
    } else if ((funct == 4) || (funct == 9)) {
        itp.dihedrals.phase.push_back(dihed_parts.phase);
        itp.dihedrals.kd.push_back(dihed_parts.kd);
        itp.dihedrals.pn.push_back(dihed_parts.pn);
        itp.dihedrals.ph0.push_back(0.0);
        itp.dihedrals.cp.push_back(0.0);
        itp.dihedrals.mult.push_back(0.0);
        itp.dihedrals.angle.push_back(0.0);
        itp.dihedrals.fc.push_back(0.0);
        itp.dihedrals.c0.push_back(0.0);
        itp.dihedrals.c1.push_back(0.0);
        itp.dihedrals.c2.push_back(0.0);
        itp.dihedrals.c3.push_back(0.0);
        itp.dihedrals.c4.push_back(0.0);
        itp.dihedrals.c5.push_back(0.0);
    } else if (funct == 3) {
        itp.dihedrals.c0.push_back(dihed_parts.c0);
        itp.dihedrals.c1.push_back(dihed_parts.c1);
        itp.dihedrals.c2.push_back(dihed_parts.c2);
        itp.dihedrals.c3.push_back(dihed_parts.c3);
        itp.dihedrals.c4.push_back(dihed_parts.c4);
        itp.dihedrals.c5.push_back(dihed_parts.c5);
        itp.dihedrals.phase.push_back(0.0);
        itp.dihedrals.kd.push_back(0.0);
        itp.dihedrals.pn.push_back(0.0);
        itp.dihedrals.ph0.push_back(0.0);
        itp.dihedrals.cp.push_back(0.0);
        itp.dihedrals.mult.push_back(0.0);
        itp.dihedrals.angle.push_back(0.0);
        itp.dihedrals.fc.push_back(0.0);
    }
    //
}

void AppendPair(s_itp& itp, 
                const int idx1, 
                const int idx2, 
                const int funct) {
    //
    itp.pairs.ai.push_back(idx1);
    itp.pairs.aj.push_back(idx2);
    itp.pairs.funct.push_back(funct);
}

int CountBonds(const s_itp &itp, int atom_index) {
    int count = 0;
    for (size_t i = 0; i < itp.bonds.ai.size(); ++i) {
        if (itp.bonds.ai[i] == atom_index || itp.bonds.aj[i] == atom_index) {
            ++count;
        }
    }
    return count;
}

bool is_bonded(const s_itp& itp, int atom1, int atom2) {
    for (size_t i = 0; i < itp.bonds.ai.size(); ++i) {
        int ai = itp.bonds.ai[i];
        int aj = itp.bonds.aj[i];
        if ((ai == atom1 && aj == atom2) || (ai == atom2 && aj == atom1)) {
            return true;
        }
    }
    return false;
}

vector<string> GetBondedResidues(const s_itp& itp, 
                                 int atom_index) {
    vector<string> bonded_resids;
    const auto& ai = itp.bonds.ai;
    const auto& aj = itp.bonds.aj;

    for (size_t i = 0; i < ai.size(); ++i) {
        int neighbor = -1;

        if (ai[i] == atom_index) {
            neighbor = aj[i]-1;
        } else if (aj[i] == atom_index) {
            neighbor = ai[i]-1;
        }

        if (neighbor >= 0 && neighbor < (int)itp.atoms.resid.size()) {
            bonded_resids.push_back(itp.atoms.resid[neighbor]);
        }
    }

    return bonded_resids;
}

vector<int> GetBondedIndex(const s_itp& itp, 
                           int atom_index) {
     vector<int> bonded_list;
     const auto& ai = itp.bonds.ai;
     const auto& aj = itp.bonds.aj;
     //
     for (size_t i = 0; i < ai.size(); ++i) {
         int neighbor = -1;
         //
         if (ai[i] == atom_index) {
             neighbor = aj[i];
         } else if (aj[i] == atom_index) {
             neighbor = ai[i];
         }
         //
         if (neighbor >= 0 && neighbor < (int)itp.atoms.resid.size()) {
             bonded_list.push_back(neighbor);
         }
     }
     //
     return bonded_list;
}

vector<string> GetbondResnmandResid(const s_itp& itp, 
                                    int atom_index, 
                                    string polnm) {
   vector<string> bonded_resnm_resid;
   const auto& ai = itp.bonds.ai;
   const auto& aj = itp.bonds.aj;

   for (size_t i = 0; i < ai.size(); ++i) {
       int neighbor = -1;

       if (ai[i] == atom_index) {
           neighbor = aj[i]-1;
       } else if (aj[i] == atom_index) {
           neighbor = ai[i]-1;
       }
       if (neighbor >= 0 && neighbor < (int)itp.atoms.resid.size()) {
           if (itp.atoms.resid[neighbor].find(polnm) != string::npos) {
               bonded_resnm_resid.push_back(itp.atoms.resid[neighbor]);
               bonded_resnm_resid.push_back(to_string(itp.atoms.resnr[neighbor]));
           }
       }
   }

   return bonded_resnm_resid;
}

pair<vector<vector<string>>, vector<int>> 
GetBondedAtomandResid(const s_itp itp, 
                      int atom_index) {
    vector<vector<string>> bonded_atom_resids;
    vector<int>            bonded_index;
    const auto& ai = itp.bonds.ai;
    const auto& aj = itp.bonds.aj;

    for (size_t i = 0; i < ai.size(); ++i) {
        int neighbor = -1;

        if (ai[i] == atom_index) {
            neighbor = aj[i]-1;
        } else if (aj[i] == atom_index) {
            neighbor = ai[i]-1;
        }

        if (neighbor >= 0 && neighbor < (int)itp.atoms.resid.size()) {
            bonded_atom_resids.push_back( {itp.atoms.atom[neighbor],
                                           itp.atoms.resid[neighbor]} );
            bonded_index.push_back(neighbor+1);
        }
    }

    return {bonded_atom_resids, bonded_index};
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Functions of itp_imit
//-----------------------------------------------------------------------------
void Readitpimitfile(const string &filename, 
                     s_itp_imit &itp_imit) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Error: cannot open " << filename << endl;
        return;
    }

    string line;
    Section current = NONE;

    while (getline(fin, line)) {
        line = trim(line);
        if (line.empty() || line[0] == ';' || line[0] == '#') continue;

        if (line[0] == '[') {
            if (line.find("bonds")           != string::npos) current = BONDS;
            else if (line.find("angles")     != string::npos) current = ANGLES;
            else if (line.find("dihedrals")  != string::npos) current = DIHEDRALS;
            else current = NONE;
            continue;
        }

        istringstream iss(line);
        switch (current) {
        case BONDS: {
            string ai, aj;
            int    funct;
            double c0 = 0.0, c1 = 0.0;
            if (iss >> ai >> aj >> funct) {
                itp_imit.bonds.ai.push_back(ai);
                itp_imit.bonds.aj.push_back(aj);
                itp_imit.bonds.funct.push_back(funct);
                if (iss >> c0 >> c1) {
                    itp_imit.bonds.c0.push_back(c0);
                    itp_imit.bonds.c1.push_back(c1);
                } else {
                    itp_imit.bonds.c0.push_back(0.0);
                    itp_imit.bonds.c1.push_back(0.0);
                }
            }
            break;
        }
        case ANGLES: {
            string ai, aj, ak;
            int funct;
            double angle = 0.0, fc = 0.0;
            if (iss >> ai >> aj >> ak >> funct) {
                itp_imit.angles.ai.push_back(ai);
                itp_imit.angles.aj.push_back(aj);
                itp_imit.angles.ak.push_back(ak);
                itp_imit.angles.funct.push_back(funct);
                if (iss >> angle >> fc) {
                    itp_imit.angles.angle.push_back(angle);
                    itp_imit.angles.fc.push_back(fc);
                } else {
                    itp_imit.angles.angle.push_back(0.0);
                    itp_imit.angles.fc.push_back(0.0);
                }
            }
            break;
        }
        case DIHEDRALS: {
            string ai, aj, ak, al;
            int    funct, mult;
            double ph0 = 0.0, cp = 0.0;
            double angle, fc;
            double phase, kd;
            int    pn;
            double c0, c1, c2, c3, c4, c5;
            if (iss >> ai >> aj >> ak >> al >> funct) {
                itp_imit.dihedrals.ai.push_back(ai);
                itp_imit.dihedrals.aj.push_back(aj);
                itp_imit.dihedrals.ak.push_back(ak);
                itp_imit.dihedrals.al.push_back(al);
                itp_imit.dihedrals.funct.push_back(funct);
                if (funct == 1) {
                    iss >> ph0 >> cp >> mult;
                    itp_imit.dihedrals.ph0.push_back(ph0);
                    itp_imit.dihedrals.cp.push_back(cp);
                    itp_imit.dihedrals.mult.push_back(mult);
                    itp_imit.dihedrals.angle.push_back(0.0);
                    itp_imit.dihedrals.fc.push_back(0.0);
                    itp_imit.dihedrals.phase.push_back(0.0);
                    itp_imit.dihedrals.kd.push_back(0.0);
                    itp_imit.dihedrals.pn.push_back(0.0);
                    itp_imit.dihedrals.c0.push_back(0.0);
                    itp_imit.dihedrals.c1.push_back(0.0);
                    itp_imit.dihedrals.c2.push_back(0.0);
                    itp_imit.dihedrals.c3.push_back(0.0);
                    itp_imit.dihedrals.c4.push_back(0.0);
                    itp_imit.dihedrals.c5.push_back(0.0);
                } else if (funct == 2) {
                    iss >> angle >> fc;
                    itp_imit.dihedrals.angle.push_back(angle);
                    itp_imit.dihedrals.fc.push_back(fc);
                    itp_imit.dihedrals.ph0.push_back(0.0);
                    itp_imit.dihedrals.cp.push_back(0.0);
                    itp_imit.dihedrals.mult.push_back(0.0);
                    itp_imit.dihedrals.phase.push_back(0.0);
                    itp_imit.dihedrals.kd.push_back(0.0);
                    itp_imit.dihedrals.pn.push_back(0.0);
                    itp_imit.dihedrals.c0.push_back(0.0);
                    itp_imit.dihedrals.c1.push_back(0.0);
                    itp_imit.dihedrals.c2.push_back(0.0);
                    itp_imit.dihedrals.c3.push_back(0.0);
                    itp_imit.dihedrals.c4.push_back(0.0);
                    itp_imit.dihedrals.c5.push_back(0.0);
                } else if ((funct == 4) || (funct == 9)) {
                    iss >> phase >> kd >> pn;
                    //cout << phase << " " << kd << " " << pn << endl;
                    itp_imit.dihedrals.phase.push_back(phase);
                    itp_imit.dihedrals.kd.push_back(kd);
                    itp_imit.dihedrals.pn.push_back(pn);
                    itp_imit.dihedrals.ph0.push_back(0.0);
                    itp_imit.dihedrals.cp.push_back(0.0);
                    itp_imit.dihedrals.mult.push_back(0.0);
                    itp_imit.dihedrals.angle.push_back(0.0);
                    itp_imit.dihedrals.fc.push_back(0.0);
                    itp_imit.dihedrals.c0.push_back(0.0);
                    itp_imit.dihedrals.c1.push_back(0.0);
                    itp_imit.dihedrals.c2.push_back(0.0);
                    itp_imit.dihedrals.c3.push_back(0.0);
                    itp_imit.dihedrals.c4.push_back(0.0);
                    itp_imit.dihedrals.c5.push_back(0.0);
                } else if (funct == 3) {
                    iss >> c0 >> c1 >> c2 >> c3 >> c4 >> c5;
                    itp_imit.dihedrals.c0.push_back(c0);
                    itp_imit.dihedrals.c1.push_back(c1);
                    itp_imit.dihedrals.c2.push_back(c2);
                    itp_imit.dihedrals.c3.push_back(c3);
                    itp_imit.dihedrals.c4.push_back(c4);
                    itp_imit.dihedrals.c5.push_back(c5);
                    // 
                    itp_imit.dihedrals.phase.push_back(0.0);
                    itp_imit.dihedrals.kd.push_back(0.0);
                    itp_imit.dihedrals.pn.push_back(0.0);
                    itp_imit.dihedrals.ph0.push_back(0.0);
                    itp_imit.dihedrals.cp.push_back(0.0);
                    itp_imit.dihedrals.mult.push_back(0.0);
                    itp_imit.dihedrals.angle.push_back(0.0);
                    itp_imit.dihedrals.fc.push_back(0.0);
                }
            }
            break;
        }
        default: break;
        }
    }
}
//-----------------------------------------------------------------------------
