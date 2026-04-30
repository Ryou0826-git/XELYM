// analyze.cpp

#include <analyze.hpp>

///////////////////////////////////////////////////////////////////////////////
///////////////////////////// Fine tools //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
s_pdb ChangeLength(s_pdb monopdb,
                   const double joint_dist, 
                   const vector<int>& joint_ref, 
                   const vector<int>& joint_mov) {
    //
    vector<double> n_st(3), n_ed(3);
    double         len_st=0.0, len_ed=0.0;
    //
    for (int i = 0; i < 3; ++i) {
        n_st[i] = monopdb.coord[joint_ref[1]-1][i] - monopdb.coord[joint_ref[0]-1][i];
        n_ed[i] = monopdb.coord[joint_mov[1]-1][i] - monopdb.coord[joint_mov[0]-1][i];
        len_st += n_st[i] * n_st[i];
        len_ed += n_ed[i] * n_ed[i];
    }

    len_st = sqrt(len_st);
    len_ed = sqrt(len_ed);
    
    for (int i = 0; i < 3; ++i) {
        monopdb.coord[joint_ref[1]-1][i] = monopdb.coord[joint_ref[0]-1][i] 
                                          + joint_dist * n_st[i] / len_st;
        monopdb.coord[joint_mov[1]-1][i] = monopdb.coord[joint_mov[0]-1][i] 
                                          + joint_dist * n_ed[i] / len_ed;
    }
    //
    return monopdb;
}

s_pdb ChangeLengthSingle(s_pdb monopdb,
                         const double joint_dist, 
                         const vector<int>& joint_mov) {
    //
    vector<double> n_mov(3);
    double         len_mov=0.0;
    //
    for (int i = 0; i < 3; ++i) {
        n_mov[i] = monopdb.coord[joint_mov[1]-1][i] - monopdb.coord[joint_mov[0]-1][i];
        len_mov += n_mov[i] * n_mov[i];
    }
    
    len_mov = sqrt(len_mov);
    
    for (int i = 0; i < 3; ++i) {
        monopdb.coord[joint_mov[1]-1][i] = monopdb.coord[joint_mov[0]-1][i] 
                                          + joint_dist * n_mov[i] / len_mov;
    }
    //
    return monopdb;
}

void apply_pbc(vector<vector<double>>& coord,
               double Lx, 
               double Ly, 
               double Lz) {
    const int natm = coord.size();

    for (int i = 0; i < natm; i++) {
        // x
        if (coord[i][0] > Lx || coord[i][0] < 0) {
            cout << "X; Larger than box, pbc apply..." << endl;
            coord[i][0] -= floor(coord[i][0] / Lx) * Lx;
        }
        // y
        if (coord[i][1] > Ly || coord[i][1] < 0) {
            cout << "Y: Larger than box, pbc apply..." << endl;
            coord[i][1] -= floor(coord[i][1] / Ly) * Ly;
        }
        // z
        if (coord[i][2] > Lz || coord[i][2] < 0) {
            cout << "Z: Larger than box, pbc apply..." << endl;
            coord[i][2] -= floor(coord[i][2] / Lz) * Lz;
        }
    }
}

s_pdb PolplusMonoPdb(s_pdb polpdb, 
                     s_pdb monopdb, 
                     int erase_p, 
                     int erase_m, 
                     int resid) {
    s_pdb connpdb;
    //
    int m = monopdb.atmnm.size();
    int p = polpdb.atmnm.size();
    int ppm = m + p;
    //
    for (int i = 0; i < ppm; ++i) {
        if ((i < p) && i != erase_p-1) {
            connpdb.atmnm.push_back(polpdb.atmnm[i]);
            connpdb.resnm.push_back(polpdb.resnm[i]);
            connpdb.resid.push_back(polpdb.resid[i]);
            connpdb.coord.push_back(polpdb.coord[i]);
            connpdb.occupancy.push_back(polpdb.occupancy[i]);
            connpdb.tempFactor.push_back(polpdb.tempFactor[i]);
        } else if ((i >= p) && i-p != erase_m-1) {
            connpdb.atmnm.push_back(monopdb.atmnm[i-p]);
            connpdb.resnm.push_back(monopdb.resnm[i-p]);
            connpdb.resid.push_back(resid);
            connpdb.coord.push_back(monopdb.coord[i-p]);
            connpdb.occupancy.push_back(monopdb.occupancy[i-p]);
            connpdb.tempFactor.push_back(monopdb.tempFactor[i-p]);
        }
    }
    //
    return connpdb;
}

void SetupFitConnect(const s_pdb& polpdb,
                     const s_pdb& monopdb,
                     const vector<int>& joint_ref,
                     const vector<int>& joint_mov,
                     s_fit& fit) {
    //
    if (joint_ref.size() != joint_mov.size()) {
        cout << endl;
        cout << "Error: fitting selection is not correct..." << endl;
    } else {
        fit.natm = joint_ref.size();
    }

    for (int i = 0; i < fit.natm; ++i) {
        fit.mass.push_back(1.0);
    }

    for (int j = joint_mov.size(); j > 0; --j) {
        int i = joint_mov[j-1];
        double x = monopdb.coord[i-1][0];
        double y = monopdb.coord[i-1][1];
        double z = monopdb.coord[i-1][2];
        fit.movcoord.push_back({x, y, z});
        fit.movcom[0] += x;
        fit.movcom[1] += y;
        fit.movcom[2] += z;
    }

    for (int i = 0; i < 3; ++i) {
        fit.movcom[i] /= joint_mov.size();
    }
    // for ref
    for (int j : joint_ref) {
        double x = polpdb.coord[j-1][0];
        double y = polpdb.coord[j-1][1];
        double z = polpdb.coord[j-1][2];
        fit.refcoord.push_back({x, y, z});
        fit.refcom[0] += x;
        fit.refcom[1] += y;
        fit.refcom[2] += z;
    }

    for (int i = 0; i < 3; ++i) {
        fit.refcom[i] /= joint_mov.size();
    }
}

s_pdb OperateTrrotConnect(s_fit& fit,
                          s_pdb& polpdb,
                          s_pdb& monopdb, 
                          const vector<int>& joint_ref, 
                          const vector<int>& joint_mov, 
                          const int resid) {
    //
    int natm = monopdb.atmnm.size();
    //
    vector<double> dmov(3, 0.0);
    vector<double> rotmov(3, 0.0);

    for (int iatm = 0; iatm < natm; ++iatm) {
        //
        // Set coordinate to origin...
        for (int j = 0; j < 3; ++j) {
            dmov[j] = monopdb.coord[iatm][j] - fit.movcom[j];
        }
        //
        // Rotate coordinate...
        for (int j = 0; j < 3; ++j) {
            rotmov[j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                rotmov[j] += fit.rot_mat(j, k) * dmov[k];
            }
        }
        //
        // Set coordinate to reference origin...
        for (int j = 0; j < 3; ++j) {
            monopdb.coord[iatm][j] = rotmov[j] + fit.refcom[j];
        }
    }

    s_pdb connpdb = PolplusMonoPdb(polpdb, monopdb, joint_ref[1], joint_mov[1], resid);

    return connpdb;
}

s_pdb OperateTrrotCap(s_fit& fit,
                      s_pdb& monopdb) {
    //
    int natm = monopdb.atmnm.size();
    //
    vector<double> dmov(3, 0.0);
    vector<double> rotmov(3, 0.0);

    for (int iatm = 0; iatm < natm; ++iatm) {
        //
        // Set coordinate to origin...
        for (int j = 0; j < 3; ++j) {
            dmov[j] = monopdb.coord[iatm][j] - fit.movcom[j];
        }
        //
        // Rotate coordinate...
        for (int j = 0; j < 3; ++j) {
            rotmov[j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                rotmov[j] += fit.rot_mat(j, k) * dmov[k];
            }
        }
        //
        // Set coordinate to reference origin...
        for (int j = 0; j < 3; ++j) {
            monopdb.coord[iatm][j] = rotmov[j] + fit.refcom[j];
        }
    }

    return monopdb;
}

vector<s_pdb> CapPolymer(s_pdb& polpdb, 
                         const double joint_dist,
                         const vector<int>& joint_ref, 
                         const vector<int>& joint_mov, 
                         const string group) {
    //
    s_pdb cappdb, cappdb_a, cappdb_e;
    vector<s_pdb> connpdb;
    int resid = polpdb.resid.back();
    string resnm = polpdb.resnm[0];

    if (boost::iequals(group, "methyl")) {
        //---------------------------------------------------------------------
        // C
        cappdb.coord.push_back({0.0, 0.0, 0.0});
        cappdb.atmnm.push_back("  CME");
        cappdb.resnm.push_back(resnm);
        cappdb.resid.push_back(1);
        cappdb.occupancy.push_back(1.0);
        cappdb.tempFactor.push_back(0.0);
        
        //
        // H
        cappdb.coord.push_back({0.635186, 0.635186, 0.635186});
        cappdb.atmnm.push_back(" HME1");
        cappdb.resnm.push_back(resnm);
        cappdb.resid.push_back(1);
        cappdb.occupancy.push_back(1.0);
        cappdb.tempFactor.push_back(0.0);

        //
        // H
        cappdb.coord.push_back({-0.635186, -0.635186, 0.635186});
        cappdb.atmnm.push_back(" HME2");
        cappdb.resnm.push_back(resnm);
        cappdb.resid.push_back(1);
        cappdb.occupancy.push_back(1.0);
        cappdb.tempFactor.push_back(0.0);

        //
        // H
        cappdb.coord.push_back({-0.635186, 0.635186, -0.635186});
        cappdb.atmnm.push_back(" HME3");
        cappdb.resnm.push_back(resnm);
        cappdb.resid.push_back(1);
        cappdb.occupancy.push_back(1.0);
        cappdb.tempFactor.push_back(0.0);

        //
        // H
        cappdb.coord.push_back({0.635186, -0.635186, -0.635186});
        cappdb.atmnm.push_back(" HME4");
        cappdb.resnm.push_back(resnm);
        cappdb.resid.push_back(1);
        cappdb.occupancy.push_back(1.0);
        cappdb.tempFactor.push_back(0.0);
        
        //
        cappdb = ChangeLengthSingle(cappdb, joint_dist, {1, 5});
        //
        s_fit fit_e;
        SetupFitConnect(polpdb, cappdb, joint_mov, {1, 5}, fit_e);
        GetTrrot(fit_e);
        s_pdb cappdb_e = OperateTrrotCap(fit_e, cappdb);
        cappdb_e.coord.pop_back();
        cappdb_e.atmnm.pop_back();
        cappdb_e.resnm.pop_back();
        cappdb_e.resid.pop_back();
        cappdb_e.occupancy.pop_back();
        cappdb_e.tempFactor.pop_back();
        //
        connpdb.push_back(cappdb_e);
        //---------------------------------------------------------------------
        //
        cappdb.atmnm[0] = "  CMA";
        cappdb.atmnm[1] = " HMA1";
        cappdb.atmnm[2] = " HMA2";
        cappdb.atmnm[3] = " HMA3";
        cappdb.atmnm[4] = " HMA4";
        //
        cappdb.resid[0] = resid;
        cappdb.resid[1] = resid;
        cappdb.resid[2] = resid;
        cappdb.resid[3] = resid;
        cappdb.resid[4] = resid;
        //
        s_fit fit_a;
        SetupFitConnect(polpdb, cappdb, joint_ref, {1, 5}, fit_a);
        GetTrrot(fit_a);
        cappdb_a = OperateTrrotCap(fit_a, cappdb);
        cappdb_a.coord.pop_back();
        cappdb_a.atmnm.pop_back();
        cappdb_a.resnm.pop_back();
        cappdb_a.resid.pop_back();
        cappdb_a.occupancy.pop_back();
        cappdb_a.tempFactor.pop_back();
        
        int index = joint_mov[1] - 1;
        polpdb.coord.erase(polpdb.coord.begin() + index);
        polpdb.atmnm.erase(polpdb.atmnm.begin() + index);
        polpdb.resnm.erase(polpdb.resnm.begin() + index);
        polpdb.resid.erase(polpdb.resid.begin() + index);
        polpdb.occupancy.erase(polpdb.occupancy.begin() + index);
        polpdb.tempFactor.erase(polpdb.tempFactor.begin() + index);
        
        int index2 = joint_ref[1] -2;
        polpdb.coord.erase(polpdb.coord.begin() + index2);
        polpdb.atmnm.erase(polpdb.atmnm.begin() + index2);
        polpdb.resnm.erase(polpdb.resnm.begin() + index2);
        polpdb.resid.erase(polpdb.resid.begin() + index2);
        polpdb.occupancy.erase(polpdb.occupancy.begin() + index2);
        polpdb.tempFactor.erase(polpdb.tempFactor.begin() + index2);
        
        connpdb.push_back(polpdb);
        connpdb.push_back(cappdb_a);
    } else {
        connpdb.push_back(polpdb);
    }
    return connpdb;
}


//
// For charge
vector<double> UpdateCharge(vector<double>& polcharge, 
                            vector<double>& monocharge, 
                            vector<int>& joint_ref, 
                            vector<int>& joint_mov) {
    //
    int index = joint_ref[1] - 1;
    polcharge.erase(polcharge.begin() + index);
    
    int index2 = joint_mov[1] - 1;
    vector<double> erasecharge = monocharge;
    erasecharge.erase(erasecharge.begin() + index2);
    int e = erasecharge.size();
    
    for (int j = 0; j < e; ++j) {
        polcharge.push_back(erasecharge[j]);
    }
    return polcharge;
}

vector<double> UpdateChargeCap(vector<double>& polcharge, 
                               const vector<int>& joint_ref, 
                               const vector<int>& joint_mov, 
                               const string& captype) {
     //
     vector<double> capcharge;
     //
     if (boost::iequals(captype, "methyl")) {
         int index = joint_mov[1] - 1;
         polcharge.erase(polcharge.begin() + index);
         int index2 = joint_ref[1] - 2;
         polcharge.erase(polcharge.begin() + index2);
         //
         capcharge = {0.0, 0.0, 0.0, 0.0};
         for (double c : polcharge) {
             capcharge.push_back(c);
         }
         //
         for (int i = 0; i < 4; ++i) {
             capcharge.push_back(0.0);
         }
     } else {
         capcharge = polcharge;
     }

     return capcharge;
}

vector<double> ControlSide(const vector<double>& monocharge, 
                           const vector<int>& joint_ref, 
                           const vector<int>& joint_mov, 
                           const double netcharge) {
   //
   cout << endl;
   cout << "joint_charge -> 0" << endl;
   int nmono = monocharge.size();
   double charge_sum = -netcharge;
   for (int i = 0; i <  nmono; ++i) {
       if ( (i != joint_ref[1] - 1) && (i != joint_mov[1] - 1) ) {
           charge_sum += monocharge[i];
       }
   }
   
   vector<double> charge(nmono, 0.0);
   for (int i = 0; i < nmono; ++i) {
       if ( (i != joint_ref[1] - 1) && (i != joint_mov[1] - 1) ) {
           charge[i] += round((monocharge[i] - charge_sum / static_cast<double>(nmono - 2)) * 100000.0) / 100000.0;
       }
   }
   
   charge_sum = - netcharge;
   for (int i = 0; i < nmono; ++i) {
       charge_sum += charge[i];
   }
   
   for (int i = 0; i < nmono; ++i) {
       if ( (i != joint_ref[1] - 1) && (i != joint_mov[1] - 1) ) {
           charge[i] -= charge_sum;
       }
   }
   //
   return charge;
}

vector<vector<int>> ExtractList(const s_itp& itp,
                                const int Np, 
                                const int natoms, 
                                const vector<string>& selpol) {
    //
    vector<vector<int>> selpol_list;
    for (int i = 0; i < Np; ++i) {
        vector<int> sel_list;
        string resid_pol = selpol[0] + to_string(i+1);
        for (int j = 0; j < natoms; ++j) {
            if ( boost::equals(resid_pol, itp.atoms.resid[j]) ) {
                if ( (boost::iequals(selpol[1], itp.atoms.atom[j])) || 
                     (boost::iequals(selpol[2], itp.atoms.atom[j])) ) {
                    sel_list.push_back(itp.atoms.nr[j]);
                }
            }
        }
        selpol_list.push_back(sel_list);
    }
    return selpol_list;
}

vector<vector<int>> ExtractListSingle(const s_itp& itp,
                                      const int Np, 
                                      const int natoms, 
                                      const vector<string>& selpol) {
    //
    vector<vector<int>> selpol_list;
    for (int i = 0; i < Np; ++i) {
        vector<int> sel_list;
        string resid_pol = selpol[0] + to_string(i+1);
        for (int j = 0; j < natoms; ++j) {
            if ( boost::equals(resid_pol, itp.atoms.resid[j]) ) {
                if (boost::iequals(selpol[1], itp.atoms.atom[j])) {
                    sel_list.push_back(itp.atoms.nr[j]);
                }
            }
        }
        selpol_list.push_back(sel_list);
    }
    //
    return selpol_list;
}

s_joint MakeJoint(const vector<string>& joint_string) {
    //
    s_joint joint;
    for (string j : joint_string) {
        vector<string> name_split = split(j);
        joint.atmnm.push_back(name_split[0]);
        joint.resnm.push_back(name_split[1]);
    }
    return joint;
}

vector<int> ExtractfromITP(s_itp& itp, 
                           const s_joint& joint) {
    int nj = joint.resnm.size();
    vector<int> selection;
    for (int i = 0; i < itp.atoms.n; ++i) {
        string trim_resid = remove_digits(itp.atoms.resid[i]);
        string atmnm = itp.atoms.atom[i];
        for (int j = 0; j < nj; ++j) {
            if ((trim_resid == joint.resnm[j]) 
                 && (atmnm == joint.atmnm[j]) ) {
                selection.push_back(i+1);
            }
        }
    }
    return selection;
}

vector<int> ExtractfromITPnonTrim(s_itp& itp, 
                                  const s_joint& joint) {
    int nj = joint.resnm.size();
    vector<int> selection;
    for (int i = 0; i < itp.atoms.n; ++i) {
        string resid = itp.atoms.resid[i];
        string atmnm = itp.atoms.atom[i];
        for (int j = 0; j < nj; ++j) {
            if ((resid == joint.resnm[j]) 
                 && (atmnm == joint.atmnm[j]) ) {
                selection.push_back(i+1);
            }
        }
    }
    return selection;
}

//-----------------------------------------------------------------------------
// Charge tools
//-----------------------------------------------------------------------------
//
// Trim joint single side charge info.
s_charge TrimJointSingleSideCharge(const string& chargefile, 
                                   const s_joint& joint) {
    //
    s_charge charge_single = prepare_charge(chargefile);

    double net_charge_single = 0.0;
    double sum_charge_single = 0.0;
    for (int i = 0; i < charge_single.nlist; ++i) {
        if ( (charge_single.resnm[i] == joint.resnm[0]) 
              && (charge_single.atmnm[i] == joint.atmnm[0]) ) {
            continue;
        } else if ( (charge_single.resnm[i] == joint.resnm[1]) 
                 && (charge_single.atmnm[i] == joint.atmnm[1]) ) {
            continue;
        } else {
            sum_charge_single += charge_single.value[i];
        }
    }

    net_charge_single -= sum_charge_single;
    //
    for (int i = 0; i < charge_single.nlist; ++i) {
        if ( (charge_single.resnm[i] == joint.resnm[0]) 
              && (charge_single.atmnm[i] == joint.atmnm[0]) ) {
            charge_single.value[i] = 0.0;
        } else if ( (charge_single.resnm[i] == joint.resnm[1]) 
                 && (charge_single.atmnm[i] == joint.atmnm[1]) ) {
            charge_single.value[i] = 0.0;
        } else {
            charge_single.value[i] 
            = round( (charge_single.value[i] + net_charge_single / static_cast<double>(charge_single.nlist - 2) ) * 10000000.0 ) / 10000000.0;
        }
    }
    //
    return charge_single;
}

//
// Trim joint single side charge info.
s_charge TrimJointDoubleSideCharge(const string& chargefile, 
                                   const s_joint& joint) {
    //
    s_charge charge_double = prepare_charge(chargefile);
    //
    double net_charge_double = 0.0;
    double sum_charge_double = 0.0;
    int cross_count = 0;
    int pol_count = 0;
    for (int i = 0; i < charge_double.nlist; ++i) {
        if (charge_double.resnm[i] == joint.resnm[0]) {
            if (charge_double.atmnm[i] == joint.atmnm[0]) {
               continue; 
            } else if (charge_double.atmnm[i] == joint.atmnm[1]) {
               continue;
            } else {
               sum_charge_double += 2.0 * charge_double.value[i];
               pol_count += 1;
            }
        } else {
            sum_charge_double += charge_double.value[i];
            cross_count += 1;
        }
    }

    net_charge_double -= sum_charge_double;
    for (int i = 0; i < charge_double.nlist; ++i) {
        if (charge_double.resnm[i] == joint.resnm[0]) {
            if (charge_double.atmnm[i] == joint.atmnm[0]) {
                charge_double.value[i] = 0.0;
            } else if (charge_double.atmnm[i] == joint.atmnm[1]) {
                charge_double.value[i] = 0.0;
            } else {
                charge_double.value[i] 
                = round( (charge_double.value[i] + net_charge_double / static_cast<double>(2*charge_double.nlist - cross_count - 4)) * 10000000.0 ) / 10000000.0;
            }
        } else {
           charge_double.value[i] 
           = round( (charge_double.value[i] + net_charge_double / static_cast<double>(2*charge_double.nlist - cross_count - 4)) * 10000000.0 ) / 10000000.0;
        }
    }
    return charge_double;
}

void ChangeCharge(s_itp& itp, 
                  const vector<int>& list,
                  const s_charge& charge) {
    //
    int ncount = 0;
    double sum_charge = 0.0;
    for (auto l : list) {
       string atom = itp.atoms.atom[l-1];
       string resid = itp.atoms.resid[l-1];
       string resid_trim = remove_digits(resid);
       for (int i = 0; i < charge.nlist; ++i) {
           if ( (charge.resnm[i] == resid_trim) && (charge.atmnm[i] == atom) ) {
               itp.atoms.charge[l-1] = charge.value[i];
               ncount += 1;
               sum_charge += charge.value[i];
               break;
           }
       }
    }

    //cout << "# of atoms: " << ncount << endl;
    //cout << "Total charge: " << sum_charge << endl;

    //double total = 0.0;
    //for (auto cc : itp.atoms.charge) total += cc;
    //cout << "Tatal of itp: " <<  total << endl;
}

//-----------------------------------------------------------------------------
s_connect ReturnConnectionPairSame(const s_itp& itp, 
                                   const s_pair& pair,
                                   const vector<vector<double>>& coord, 
                                   const vector<double>& box, 
                                   const double& rc) {
    //
    int l1 = pair.list1.size();
    int l2 = pair.list2.size();
    if (l1 != l2) throw runtime_error("Size mismatch between list1 and list2");
    //
    s_connect connect_first;
    for (int i = 0; i < l2; ++i) {
        int selpol2 = pair.list2[i];
        int selpol1_pair = pair.list1[i];
        vector<double> selpol2_crd = coord[selpol2-1];
        double length_min = 100.0;
        int    selpol1_min = -1;
        for (int j = 0; j < l1; ++j) {
            if (j != i) {
                int selpol1 = pair.list1[j];
                vector<double> selpol1_crd = coord[selpol1-1];
                double length = dist2_pbc(selpol2_crd, selpol1_crd, box);
                length = sqrt(length);
                if (length < length_min) {
                    length_min = length;
                    selpol1_min = selpol1;
                }
            }
        }
        //
        if (selpol1_min == -1) continue;
        //
        string bonded_sp1m_resid = itp.atoms.resid[selpol1_min-1];
        vector<string> bonded_sp1p_resids = GetBondedResidues(itp, selpol1_pair);
        vector<string> bonded_sp1m_resids = GetBondedResidues(itp, selpol1_min);
        vector<string> bonded_sp2_resids = GetBondedResidues(itp, selpol2);
        //
        int bonded_sp1m_resids_size = bonded_sp1m_resids.size();
        int nbond_sp1m = return_nbond_from_itp(itp, selpol1_min);
        int bonded_sp2_resids_size = bonded_sp2_resids.size();
        int nbond_sp2 = return_nbond_from_itp(itp, selpol2);
        //
        if ( find(bonded_sp1p_resids.begin(), bonded_sp1p_resids.end(), bonded_sp1m_resid) 
            != bonded_sp1p_resids.end() ) continue;
        if (bonded_sp1m_resids_size >= nbond_sp1m) continue;
        if (bonded_sp2_resids_size >= nbond_sp2) continue;
        //
        bool b_sign = is_bonded(itp, selpol1_min, selpol2);
        if (!b_sign && length_min < rc) {
            connect_first.joint1.push_back(selpol1_min);
            connect_first.joint2.push_back(selpol2);
            connect_first.joint_resid1.push_back(itp.atoms.resid[selpol1_min-1]);
            connect_first.joint_resid2.push_back(itp.atoms.resid[selpol2-1]);
        }
    }

    //
    // second check (polymer-polymer)
    s_connect connect_second;
    int jp1_size = connect_first.joint1.size(); 
    vector<string> vec_dup, vec_alr;
    for (int i = 0; i < jp1_size; ++i) {
        int    sel2 = connect_first.joint2[i];
        string resid2 = connect_first.joint_resid2[i];
        vector<double> sel2_crd = coord[sel2-1];
        double length_min = 100.0;
        int    idx_min = -1;
        for (int j = 0; j < jp1_size; ++j) {
            int sel1 = connect_first.joint1[j];
            vector<double> sel1_crd = coord[sel1-1];
            string resid1 = connect_first.joint_resid1[j];
            //
            if (resid2 == resid1) continue;
            double length = dist2_pbc(sel1_crd, sel2_crd, box);
            length = sqrt(length);
            if (length < length_min) {
                length_min = length;
                idx_min = j;
            }
        }
        //
        if (idx_min== -1) throw runtime_error("index of minimum distance is not found!");
        //
        int sel_min = connect_first.joint1[idx_min];
        if ( find(vec_dup.begin(), vec_dup.end(), itp.atoms.resid[sel_min-1]) != vec_dup.end()
             || find(vec_alr.begin(), vec_alr.end(), itp.atoms.resid[sel_min-1]) != vec_alr.end())
          continue;
        vec_dup.push_back(itp.atoms.resid[sel_min-1]);
        vec_alr.push_back(itp.atoms.resid[sel2-1]);
        connect_second.joint1.push_back(sel_min);
        connect_second.joint2.push_back(sel2);
        connect_second.joint_resid1.push_back(itp.atoms.resid[sel2-1]);
        connect_second.joint_resid2.push_back(itp.atoms.resid[sel_min-1]);
    }
    //
    return connect_second;
}

s_connect ReturnConnectionPairDiff(const s_itp& itp, 
                                   const s_pair& pair1,
                                   const vector<int>& list2,
                                   const vector<vector<double>>& coord, 
                                   const vector<double>& box, 
                                   const double& rc) {
    //
    int l1 = pair1.list1.size();
    int pl1 = pair1.list2.size();
    int l2 = list2.size();
    bool pair_sign = false;
    //
    if (l1 == pl1) {
        pair_sign = true;
    } else if (pl1 == 0) {
        pair_sign = false;
    } else {
        throw runtime_error("(l1 != pl1) || (pl1 != 0)");
    }
    //
    s_connect connect_first;
    //
    if (pair_sign) {
        for (int i = 0; i < l2; ++i) {
            int idx2 = list2[i];
            string idx2_resid = itp.atoms.resid[idx2-1];
            vector<double> selpol2_crd = coord[idx2-1];
            double length_min = 100.0;
            int    idx1_min = -1, idx1_pair_min = -1;
            for (int j = 0; j < l1; ++j) {
                int idx1 = pair1.list1[j];
                int idx1_pair = pair1.list2[j];
                vector<double> selpol1_crd = coord[idx1-1];
                double length = dist2_pbc(selpol2_crd, selpol1_crd, box);
                length = sqrt(length);
                if (length < length_min) {
                    length_min = length;
                    idx1_min = idx1;
                    idx1_pair_min = idx1_pair;
                }
            }
            //
            //if (length_min < rc) {cout << length_min << endl;}
            //
            vector<string> bonded_i1pm_resids = GetBondedResidues(itp, idx1_pair_min);
            vector<string> bonded_i1m_resids = GetBondedResidues(itp, idx1_min);
            vector<string> bonded_i2_resids = GetBondedResidues(itp, idx2);
            int bonded_i2_resids_size = bonded_i2_resids.size();
            int nbond_idx2 = return_nbond_from_itp(itp, idx2);
            if ( find(bonded_i1pm_resids.begin(), bonded_i1pm_resids.end(), idx2_resid)
                != bonded_i1pm_resids.end() ) continue;
            if (bonded_i2_resids_size >= nbond_idx2) continue;
            //
            bool b_sign = is_bonded(itp, idx1_min, idx2); 
            if (!b_sign && length_min < rc) {
                connect_first.joint1.push_back(idx1_min);
                connect_first.joint2.push_back(idx2);
                connect_first.joint_resid1.push_back(itp.atoms.resid[idx1_min-1]);
                connect_first.joint_resid2.push_back(itp.atoms.resid[idx2-1]);
            }
        }
    } else {
        for (int i = 0; i < l2; ++i) {
            int idx2 = list2[i];
            string idx2_resid = itp.atoms.resid[idx2-1];
            vector<double> selpol2_crd = coord[idx2-1];
            double length_min = 100.0;
            int    idx1_min = -1;
            for (int j = 1; j < l1; ++j) {
                int idx1 = pair1.list1[j];
                vector<double> selpol1_crd = coord[idx1-1];
                double length = dist2_pbc(selpol2_crd, selpol1_crd, box);
                length = sqrt(length);
                if (length < length_min) {
                    length_min = length;
                    idx1_min = idx1;
                }
            }
            //
            //if (length_min < rc) {cout << length_min << endl;}
            //
            vector<string> bonded_i1m_resids = GetBondedResidues(itp, idx1_min);
            vector<string> bonded_i2_resids = GetBondedResidues(itp, idx2);
            int bonded_i2_resids_size = bonded_i2_resids.size();
            int nbond_idx2 = return_nbond_from_itp(itp, idx2);
            if (bonded_i2_resids_size >= nbond_idx2) continue;
            //
            bool b_sign = is_bonded(itp, idx1_min, idx2); 
            if (!b_sign && length_min < rc) {
                connect_first.joint1.push_back(idx1_min);
                connect_first.joint2.push_back(idx2);
                connect_first.joint_resid1.push_back(itp.atoms.resid[idx1_min-1]);
                connect_first.joint_resid2.push_back(itp.atoms.resid[idx2-1]);
            }
        }
    }

    //
    // second check
    s_connect connect_second;
    int jp1_size = connect_first.joint1.size(); 
    vector<string> vec_dup, vec_alr;
    for (int i = 0; i < jp1_size; ++i) {
        int    sel2 = connect_first.joint2[i];
        string resid2 = connect_first.joint_resid2[i];
        vector<double> sel2_crd = coord[sel2-1];
        double length_min = 100.0;
        int    idx_min = -1;
        for (int j = 0; j < jp1_size; ++j) {
            int sel1 = connect_first.joint1[j];
            vector<double> sel1_crd = coord[sel1-1];
            string resid1 = connect_first.joint_resid1[j];
            //
            if (resid2 == resid1) continue;
            double length = dist2_pbc(sel1_crd, sel2_crd, box);
            length = sqrt(length);
            if (length < length_min) {
                length_min = length;
                idx_min = j;
            }
        }
        //
        if (idx_min== -1) throw runtime_error("index of minimum distance is not found!");
        //
        int sel_min = connect_first.joint1[idx_min];
        if ( find(vec_dup.begin(), vec_dup.end(), itp.atoms.resid[sel_min-1]) != vec_dup.end()
             || find(vec_alr.begin(), vec_alr.end(), itp.atoms.resid[sel_min-1]) != vec_alr.end())
          continue;
        vec_dup.push_back(itp.atoms.resid[sel_min-1]);
        vec_alr.push_back(itp.atoms.resid[sel2-1]);
        connect_second.joint1.push_back(sel_min);
        connect_second.joint2.push_back(sel2);
        connect_second.joint_resid1.push_back(itp.atoms.resid[sel2-1]);
        connect_second.joint_resid2.push_back(itp.atoms.resid[sel_min-1]);
    }
    //
    return connect_second;
}

s_connect erase_s_connect_by_list(s_connect& connect, 
                                  const vector<int>& list) {
    //
    vector<int> list_tmp = list;
    erase_by_index_list(connect.joint1, list_tmp);
    erase_by_index_list(connect.joint_resid1, list_tmp);
    erase_by_index_list(connect.joint2, list_tmp);
    erase_by_index_list(connect.joint_resid2, list_tmp);
    //
    return connect;
}

s_itp return_s_itp_change_charge(s_itp& itp, 
                                 const vector<int>& list, 
                                 const double& aft_charge) {
    for (int l : list) {
        int idx = l - 1;
        itp.atoms.charge[idx] = aft_charge;
    }
    //
    double sum_charge = 0.0;
    for (auto c : itp.atoms.charge) {
        sum_charge += c;
    }
    cout << endl;
    cout << "Total charge: " << sum_charge << endl;
    //
    return itp;
}

//
// get data from topology file --->
vector<int> GetSelselection(const s_top& top, 
                            const s_joint& j_sel) {
    vector<int> sel_index;
    int sel_seq = -1;
    for (auto itp : top.itp_list) {
        sel_index = ExtractfromITP(itp, j_sel);
        int ns = sel_index.size();
        sel_seq += 1;
        if (ns > 0) break; 
    }
    if (sel_seq == -1) throw runtime_error("itp is not included in topfile...");
    // 
    int nbf_sel = 0;
    for (int i = 0; i < sel_seq; ++i) { nbf_sel += top.molecules.natom[i]; }
    //
    int Ns = top.molecules.nmol[sel_seq];

    vector<int> sel_selection;
    if (Ns > 1) {
        if (sel_index.size() == 1) {
            int i_temp = nbf_sel + sel_index[0];
            for (int i = 0; i < Ns; ++i) {
                sel_selection.push_back(i_temp);
                i_temp += top.itp_list[sel_seq].atoms.n;
            }
        } else {
            throw runtime_error("index of sel is not 1 (this analysis is not supported in this case)");
        }
    } else {
        for (auto sel_int : sel_index) {
            int sel_idx = sel_int + nbf_sel;
            sel_selection.push_back(sel_idx);
        }
    }

    return sel_selection;
}

vector<int> GetSelselectionNeo(const s_top& top, 
                               const s_joint& j_sel, 
                               bool trim_sign) {
    vector<int> sel_index;
    int sel_seq = -1;
    for (auto itp : top.itp_list) {
        if (trim_sign) {
            sel_index = ExtractfromITP(itp, j_sel);
        } else {
            sel_index = ExtractfromITPnonTrim(itp, j_sel);
        }
        int ns = sel_index.size();
        sel_seq += 1;
        if (ns > 0) break; 
    }
    if (sel_seq == -1) throw runtime_error("itp is not included in topfile...");
    // 
    int nbf_sel = 0;
    for (int i = 0; i < sel_seq; ++i) { nbf_sel += top.molecules.natom[i]; }
    //
    int Ns = top.molecules.nmol[sel_seq];

    vector<int> sel_selection;
    if (Ns > 1) {
        if (sel_index.size() == 1) {
            int i_temp = nbf_sel + sel_index[0];
            for (int i = 0; i < Ns; ++i) {
                sel_selection.push_back(i_temp);
                i_temp += top.itp_list[sel_seq].atoms.n;
            }
        } else {
            throw runtime_error("index of sel is not 1 (this analysis is not supported in this case)");
        }
    } else {
        for (auto sel_int : sel_index) {
            int sel_idx = sel_int + nbf_sel;
            sel_selection.push_back(sel_idx);
        }
    }

    return sel_selection;
}

// <----

s_charge ChangeChargeSumZero(const string& chargefile) {
    //
    s_charge charge_pol = prepare_charge(chargefile);
    double cp_sum = 0.0;
    for (auto cp : charge_pol.value) {
        cp_sum += cp;
    }
    //
    double cp_sum_new = 0.0;
    if (cp_sum != 0.0) {
        cout << endl;
        cout << "Total charge is not zero, correcting..." << endl;
        double cp_split = cp_sum / charge_pol.nlist;
        for (int i = 0; i < charge_pol.nlist; ++i) {
            cout << charge_pol.atmnm[i] << "_" << charge_pol.resnm[i] 
                 << ": " << charge_pol.value[i] << " --> ";
            //charge_pol.value[i] = round( (charge_pol.value[i] - cp_split) * 10000000000.0 ) / 10000000000.0;
            charge_pol.value[i] = charge_pol.value[i] - cp_split;
            cout << charge_pol.value[i] << endl;
            cp_sum_new += charge_pol.value[i];
        }
    } else {
        cout << endl;
        cout << "Total charge is zero, need not correcting..." << endl;
    }

    return charge_pol;
}

vector<s_charge> ChangeChargeListSumZero(const vector<string>& chargefile) {
    //
    vector<s_charge> charge_x_list;
    s_charge charge_void;
    double cx_sum = 0.0, cx_sum_new = 0.0;
    for (auto cf : chargefile) {
        //cout << cf << endl;
        if (cf == "") {
            charge_x_list.push_back(charge_void);
        } else {
            s_charge charge_x = prepare_charge(cf);
            charge_x_list.push_back(charge_x);
            for (auto cx : charge_x.value) cx_sum += cx;
        }
    }
    //
    //cout << "# of charge_x_list: " << charge_x_list.size() << endl;
    if (cx_sum != 0.0) {
        //cout << endl;
        //cout << "Total charge is not zero, correcting..." << endl;
        //
        int n = 0;
        int nc = charge_x_list.size();
        for (int i = 0; i < nc; ++i) {
            n += charge_x_list[i].nlist;
        }
        //
        //cout << cx_sum << endl;
        cx_sum /= n;
        //cout << cx_sum << endl;

        for (int i = 0; i < nc; ++i) {
            int natom =  charge_x_list[i].value.size();
            for (int j = 0; j < natom; ++j) {
               //cout << charge_x_list[i].atmnm[j] << "_" << charge_x_list[i].resnm[j] 
               //     << ": " << charge_x_list[i].value[j] << " --> ";
               //charge_x_list[i].value[j] = round( (charge_x_list[i].value[j] - cx_sum) * 10000000000.0 ) / 10000000000.0;
               charge_x_list[i].value[j] = charge_x_list[i].value[j] - cx_sum;
               cx_sum_new += charge_x_list[i].value[j];
               //cout << charge_x_list[i].value[j] << endl;
            }
        }
        //cout << "size of charge_x_list: " << charge_x_list.size() << endl;
        //cout << "TTTTTTTTTT: " << cx_sum_new << endl;
    } else {
        cout << endl;
        cout << "Total charge is zero, need not correcting..." << endl;
    }
    return charge_x_list;
}

// Bonded routine ---->
vector<vector<int>> ReturnBondedListFour(const s_itp& itp,
                                         const vector<int> selcross1_single, 
                                         const vector<int> selcross2_single,
                                         const vector<int> selcross3_single,
                                         const vector<int> selcross4_single,
                                         const string& resid) {
    int Nc = selcross1_single.size();
    vector<vector<int>> sel_bonded_list;
    for (int i = 0; i < Nc; ++i) {
        int sel1 = selcross1_single[i];
        int sel2 = selcross2_single[i];
        int sel3 = selcross3_single[i];
        int sel4 = selcross4_single[i];
        //
        vector<int> sel1_bonded_list = GetBondedIndex(itp, sel1);
        int sel1_bonded = GetIndexMatchResid(itp, sel1_bonded_list, resid);
        vector<int> sel2_bonded_list = GetBondedIndex(itp, sel2);
        int sel2_bonded = GetIndexMatchResid(itp, sel2_bonded_list, resid);
        vector<int> sel3_bonded_list = GetBondedIndex(itp, sel3);
        int sel3_bonded = GetIndexMatchResid(itp, sel3_bonded_list, resid);
        vector<int> sel4_bonded_list = GetBondedIndex(itp, sel4);
        int sel4_bonded = GetIndexMatchResid(itp, sel4_bonded_list, resid);
        //
        vector<int> cell_bonded_list = {sel1_bonded, sel2_bonded, 
                                        sel3_bonded, sel4_bonded};
        //cout << i+1 << ": " << sel1_bonded << " " << sel2_bonded 
        //     << " " << sel3_bonded << " " << sel4_bonded << endl;
        sel_bonded_list.push_back(cell_bonded_list);
    }

    return sel_bonded_list;
}
// <----

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////// FINE SUBROUTINE (TOPOLOGY) ///////////////////////////
///////////////////////////////////////////////////////////////////////////////
void UpdateBond(s_itp& itp,
                const s_itp_imit& imit_itp, 
                const s_bond_set& bond_set) {
    //
    int id1 = bond_set.id1;
    int id2 = bond_set.id2;
    vector<string> id1_set = bond_set.id1_set;
    vector<string> id2_set = bond_set.id2_set;
    //
    cout << "BOND SECTION: " << endl;
    int nb = imit_itp.bonds.ai.size();
    double bond_c0 = 0.0, bond_c1 = 0.0;
    int funct = 1;
    for (int j = 0; j < nb; ++j) {
        vector<string> imit_ai = split(imit_itp.bonds.ai[j], '_');
        vector<string> imit_aj = split(imit_itp.bonds.aj[j], '_');
        if ( ( imit_ai == id1_set && imit_aj == id2_set ) 
             || ( imit_ai == id2_set && imit_aj == id1_set ) ) {
            funct = imit_itp.bonds.funct[j];
            bond_c0 = imit_itp.bonds.c0[j];
            bond_c1 = imit_itp.bonds.c1[j];
        }
    }
    cout << "Create " << id1 << " and " << id2 << " bond" << endl;
    AppendBond(itp, id1, id2, funct, bond_c0, bond_c1);
}

void ChangeBond(s_itp& itp, 
                const s_itp_imit& imit_itp, 
                const vector<int>& list1, 
                const vector<int>& list2) {
    //
    int l1 = list1.size();
    int l2 = list2.size();
    if (l1 != l2) throw runtime_error("l1 != l2");

    cout << endl;
    cout << "CHANGE BOND SECTION: " << endl;
    int nb  = itp.bonds.ai.size();
    int nbi = imit_itp.change_bonds.ai.size();
    for (int i = 0; i < nb; ++i) {
        int ai = itp.bonds.ai[i];
        int aj = itp.bonds.aj[i];
        //
        for (int l = 0; l < l1; ++l) {
            int id1 = list1[l];
            int id2 = list2[l];
            //
            if ( (ai == id2  && aj == id1) || (ai == id1  && aj == id2) ) continue;
            //
            if ( (ai == id2 || aj == id2) || (ai == id1 || aj == id1) ) {
                vector<string> ai_set = return_atom_set(itp, ai);
                vector<string> aj_set = return_atom_set(itp, aj);
                for (int j = 0; j < nbi; ++j) {
                    vector<string> imit_ai = split(imit_itp.change_bonds.ai[j], '_');
                    vector<string> imit_aj = split(imit_itp.change_bonds.aj[j], '_');
                    
                    if ( ( imit_ai == ai_set && imit_aj == aj_set ) 
                         || ( imit_ai == aj_set && imit_aj == ai_set ) ) {
                        int    funct   = imit_itp.change_bonds.funct[j];
                        double bond_c0 = imit_itp.change_bonds.c0[j];
                        double bond_c1 = imit_itp.change_bonds.c1[j];
                        itp.bonds.funct[i] = funct;
                        itp.bonds.c0[i]    = bond_c0;
                        itp.bonds.c1[i]    = bond_c1;
                        cout << "Change "<< ai << " and " << aj << " bond ";
                        cout << "( " <<  imit_itp.change_bonds.ai[j] 
                             << " and " << imit_itp.change_bonds.aj[j] << " )" << endl;
                    }
                }
            }
        }
    }
}

void UpdateAngle(s_itp& itp, 
                 const s_itp_imit& imit_itp, 
                 s_bond_set& bond_set) {
    //
    cout << endl;
    cout << "ANGLE SECTION: " << endl;
    int na = imit_itp.angles.ai.size();
    //
    // Creating angles...
    // don. side
    int nda = bond_set.id1_bond.size();
    for (int j = 0; j < nda; ++j) {
        //
        vector<vector<string>> angle_set;
        vector<vector<string>> angle_imit_set;
        bond_set.id1_bond[j][1] = remove_digits(bond_set.id1_bond[j][1]);
        //
        angle_set.push_back(bond_set.id1_bond[j]);
        angle_set.push_back(bond_set.id1_set);
        angle_set.push_back(bond_set.id2_set);
        //
        for (int k = 0; k < na; ++k) {
            vector<string> imit_ai = split(imit_itp.angles.ai[k], '_');
            vector<string> imit_aj = split(imit_itp.angles.aj[k], '_');
            vector<string> imit_ak = split(imit_itp.angles.ak[k], '_');
            angle_imit_set = {imit_ai, imit_aj, imit_ak};
            //
            if ( (bond_set.id1 == bond_set.id1_index[j]) 
                 || (bond_set.id2 == bond_set.id1_index[j]) ){
                 continue;
            }
            vector<int> angle_index = {bond_set.id1_index[j], bond_set.id1, bond_set.id2};
            //auto [osign, angle_index_replace] = compare_and_reindex(angle_imit_set, angle_set, angle_index);
            auto [osign, angle_index_replace] = check_and_return_index(angle_imit_set, angle_set, angle_index);
            if (osign) {
                cout << "Create "
                     << angle_index_replace[0] << " and "
                     << angle_index_replace[1] << " and "
                     << angle_index_replace[2] << " angle " << endl;
                AppendAngle(itp, angle_index_replace[0], angle_index_replace[1], 
                            angle_index_replace[2], imit_itp.angles.funct[k], 
                            imit_itp.angles.angle[k], imit_itp.angles.fc[k]);
            }
        }
    }
    
    //
    // acc. side
    int naa = bond_set.id2_bond.size();
    for (int j = 0; j < naa; ++j) {
        //
        vector<vector<string>> angle_set;
        vector<vector<string>> angle_imit_set;
        bond_set.id2_bond[j][1] = remove_digits(bond_set.id2_bond[j][1]);
        //
        angle_set.push_back(bond_set.id2_bond[j]);
        angle_set.push_back(bond_set.id2_set);
        angle_set.push_back(bond_set.id1_set);
        //
        for (int k = 0; k < na; ++k) {
            //
            vector<string> imit_ai = split(imit_itp.angles.ai[k], '_');
            vector<string> imit_aj = split(imit_itp.angles.aj[k], '_');
            vector<string> imit_ak = split(imit_itp.angles.ak[k], '_');
            angle_imit_set = {imit_ai, imit_aj, imit_ak};
            //
            if ( (bond_set.id2 == bond_set.id2_index[j]) 
                || (bond_set.id1 == bond_set.id2_index[j]) ) {
                continue;
            }
            vector<int> angle_index = {bond_set.id2_index[j], bond_set.id2, bond_set.id1};
            //auto [osign, angle_index_replace] = compare_and_reindex(angle_imit_set, angle_set, angle_index);
            auto [osign, angle_index_replace] = check_and_return_index(angle_imit_set, angle_set, angle_index);
            if (osign) {
                cout << "Create " 
                     << angle_index_replace[0] << " and "
                     << angle_index_replace[1] << " and " 
                     << angle_index_replace[2] << " angle" << endl;
                AppendAngle(itp, angle_index_replace[0], angle_index_replace[1],
                            angle_index_replace[2], imit_itp.angles.funct[k],
                            imit_itp.angles.angle[k], imit_itp.angles.fc[k]);
            }
        }
    }
    //
}

void ChangeAngle(s_itp& itp, 
                 const s_itp_imit& imit_itp, 
                 const vector<int>& list1, 
                 const vector<int>& list2) {
    //
    cout << endl;
    cout << "CHANGE ANGLE SECTION: " << endl;
    int l1 = list1.size();
    int l2 = list2.size();
    if (l1 != l2) throw runtime_error("l1 != l2");
    //
    int na = itp.angles.ai.size();
    int nai = imit_itp.change_angles.ai.size();
    //
    for (int i = 0; i < na; ++i) {
        int ai = itp.angles.ai[i];
        vector<string> ai_set = return_atom_set(itp, ai);
        int aj = itp.angles.aj[i];
        vector<string> aj_set = return_atom_set(itp, aj);
        int ak = itp.angles.ak[i];
        vector<string> ak_set = return_atom_set(itp, ak);
        //
        vector<int> angle_index = {ai, aj, ak};
        vector<vector<string>> angle_set = {ai_set, aj_set, ak_set};
        
        for (int l = 0; l < l1; ++l) {
            int id1 = list1[l];
            int id2 = list2[l];
            //
            bool has_id1 = find(angle_index.begin(), angle_index.end(), id1) != angle_index.end();
            bool has_id2 = find(angle_index.begin(), angle_index.end(), id2) != angle_index.end();
            //
            if ( has_id1 && has_id2 ) continue;
            //
            if ( has_id1 || has_id2 ) {
                for (int j = 0; j < nai; ++j) {
                    vector<string> imit_ai = split(imit_itp.change_angles.ai[j], '_');
                    vector<string> imit_aj = split(imit_itp.change_angles.aj[j], '_');
                    vector<string> imit_ak = split(imit_itp.change_angles.ak[j], '_');
                    vector<vector<string>> angle_imit_set = {imit_ai, imit_aj, imit_ak};
                    auto [osign, angle_index_replace] = check_and_return_index(angle_imit_set, angle_set, angle_index);
                    if (osign) {
                        cout << "Change " 
                             << angle_index_replace[0] << " and "
                             << angle_index_replace[1] << " and " 
                             << angle_index_replace[2] << " angle ";
                        cout << "( " << imit_itp.change_angles.ai[j]
                             << " and " << imit_itp.change_angles.aj[j]
                             << " and " << imit_itp.change_angles.ak[j] << " )" << endl;
                        //
                        itp.angles.ai[i] = angle_index_replace[0];
                        itp.angles.aj[i] = angle_index_replace[1];
                        itp.angles.ak[i] = angle_index_replace[2];
                        itp.angles.funct[i] = imit_itp.change_angles.funct[j];
                        itp.angles.angle[i] = imit_itp.change_angles.angle[j];
                        itp.angles.fc[i]    = imit_itp.change_angles.fc[j];
                        // 
                        break;
                    }
                }
            }
        }
    }
}

void UpdateDihedral(s_itp& itp, 
                    const s_itp_imit& imit_itp, 
                    s_bond_set& bond_set) {
    cout << endl;
    cout << "DIHEDRAL SECTION:" << endl;
    int nd = imit_itp.dihedrals.ai.size();
    int nda = bond_set.id1_bond.size();
    int naa = bond_set.id2_bond.size();
    //
    // don. side
    for (int l = 0; l < nda; ++l) {
        //
        auto [id1_bond2, id1_index2] = GetBondedAtomandResid(itp, bond_set.id1_index[l]);
        int nb2 = id1_index2.size();
        
        if (nb2 == 1) continue;
    
        for (int m = 0; m < nb2; ++m) {
            id1_bond2[m][1] = remove_digits(id1_bond2[m][1]);
            vector<vector<string>> dihedral_set;
            dihedral_set.push_back(id1_bond2[m]);
            dihedral_set.push_back(bond_set.id1_bond[l]);
            dihedral_set.push_back(bond_set.id1_set);
            dihedral_set.push_back(bond_set.id2_set);
            //
            if ( (bond_set.id1 == id1_index2[m]) 
                || (bond_set.id1 == bond_set.id1_index[l]) 
                || (bond_set.id2 == id1_index2[m]) 
                || (bond_set.id2 == bond_set.id1_index[l]) ) {
                continue;
            }
            //
            vector<int> dihed_index = {id1_index2[m], bond_set.id1_index[l], bond_set.id1, bond_set.id2};
            for (int n = 0; n < nd; ++n) {
                //
                vector<string> imit_ai = split(imit_itp.dihedrals.ai[n], '_');
                vector<string> imit_aj = split(imit_itp.dihedrals.aj[n], '_');
                vector<string> imit_ak = split(imit_itp.dihedrals.ak[n], '_');
                vector<string> imit_al = split(imit_itp.dihedrals.al[n], '_');
                vector<vector<string>> 
                dihedral_imit_set = {imit_ai, imit_aj, imit_ak, imit_al};
                auto [osign, dihed_index_replace] = check_and_return_index(dihedral_imit_set, dihedral_set, dihed_index);
                if (osign) {
                    cout << "Create " 
                         << dihed_index_replace[0] << " and " 
                         << dihed_index_replace[1] << " and " 
                         << dihed_index_replace[2] << " and " 
                         << dihed_index_replace[3] << " dihedral" << endl;
                    s_dihed_parts dihed_parts;
                    if ( (imit_itp.dihedrals.funct[n] == 1) 
                        || (imit_itp.dihedrals.funct[n] == 2) ) {
                        //
                        dihed_parts.ph0   = imit_itp.dihedrals.ph0[n];
                        dihed_parts.cp    = imit_itp.dihedrals.cp[n];
                        dihed_parts.mult  = imit_itp.dihedrals.mult[n];
                        dihed_parts.angle = imit_itp.dihedrals.angle[n];
                        dihed_parts.fc    = imit_itp.dihedrals.fc[n];
                    } else if ( (imit_itp.dihedrals.funct[n] == 4) 
                        || (imit_itp.dihedrals.funct[n] == 9) ) {
                        //
                        dihed_parts.phase = imit_itp.dihedrals.phase[n];
                        dihed_parts.kd    = imit_itp.dihedrals.kd[n];
                        dihed_parts.pn    = imit_itp.dihedrals.pn[n];
                    } else if (imit_itp.dihedrals.funct[n] == 3) {
                        dihed_parts.c0    = imit_itp.dihedrals.c0[n];
                        dihed_parts.c1    = imit_itp.dihedrals.c1[n];
                        dihed_parts.c2    = imit_itp.dihedrals.c2[n];
                        dihed_parts.c3    = imit_itp.dihedrals.c3[n];
                        dihed_parts.c4    = imit_itp.dihedrals.c4[n];
                        dihed_parts.c5    = imit_itp.dihedrals.c5[n];
                    }
                    AppendDihedral(itp, 
                                   dihed_index_replace[0], 
                                   dihed_index_replace[1], 
                                   dihed_index_replace[2], 
                                   dihed_index_replace[3], 
                                   imit_itp.dihedrals.funct[n],
                                   dihed_parts);
                    AppendPair(itp, 
                               dihed_index_replace[0], 
                               dihed_index_replace[3], 1);
                }       
            }           
            //}              
        }              
    }
    //------------------------------------//
    //
    // acc. side
    cout << endl;
    for (int l = 0; l < naa; ++l) {
        auto [id2_bond2, id2_index2] = GetBondedAtomandResid(itp, bond_set.id2_index[l]);
        int nb2 = id2_index2.size();
        for (int m = 0; m < nb2; ++m) {
            id2_bond2[m][1] = remove_digits(id2_bond2[m][1]);
            vector<vector<string>> dihedral_set;
            dihedral_set.push_back(id2_bond2[m]);
            dihedral_set.push_back(bond_set.id2_bond[l]);
            dihedral_set.push_back(bond_set.id2_set);
            dihedral_set.push_back(bond_set.id1_set);
            //if ( (acc_bond2[m][0] != don_set[0]) && (acc_bond2[m][1] != don_set[1]) ) {
            if (id2_bond2[m][0] != bond_set.id1_set[0]) {
                for (int n = 0; n < nd; ++n) {
                    if ( (bond_set.id2 == id2_index2[m]) 
                        || (bond_set.id2 == bond_set.id2_index[l]) 
                        || (bond_set.id1 == id2_index2[m]) 
                        || (bond_set.id1 == bond_set.id2_index[l]) ) {
                        continue;
                    }
                    vector<string> imit_ai = split(imit_itp.dihedrals.ai[n], '_');
                    vector<string> imit_aj = split(imit_itp.dihedrals.aj[n], '_');
                    vector<string> imit_ak = split(imit_itp.dihedrals.ak[n], '_');
                    vector<string> imit_al = split(imit_itp.dihedrals.al[n], '_');
                    vector<vector<string>>
                    dihedral_imit_set = {imit_ai, imit_aj, imit_ak, imit_al};
                    vector<int> dihed_index = {id2_index2[m], bond_set.id2_index[l], bond_set.id2, bond_set.id1};
                    //auto [osign, dihed_index_replace] = compare_and_reindex(dihedral_imit_set, dihedral_set, dihed_index);
                    auto [osign, dihed_index_replace] = check_and_return_index(dihedral_imit_set, dihedral_set, dihed_index);
                    if (osign) {
                         cout << "Create " 
                              << dihed_index_replace[0] << " and "
                              << dihed_index_replace[1] << " and " 
                              << dihed_index_replace[2] << " and " 
                              << dihed_index_replace[3] << " dihedral" << endl;
                         s_dihed_parts dihed_parts;
                         if ( (imit_itp.dihedrals.funct[n] == 1) 
                             || (imit_itp.dihedrals.funct[n] == 2) ) {
                             //
                             dihed_parts.ph0   = imit_itp.dihedrals.ph0[n];
                             dihed_parts.cp    = imit_itp.dihedrals.cp[n];
                             dihed_parts.mult  = imit_itp.dihedrals.mult[n];
                             dihed_parts.angle = imit_itp.dihedrals.angle[n];
                             dihed_parts.fc    = imit_itp.dihedrals.fc[n];
                         } else if ( (imit_itp.dihedrals.funct[n] == 4) 
                             || (imit_itp.dihedrals.funct[n] == 9) ) {
                             //
                             dihed_parts.phase = imit_itp.dihedrals.phase[n];
                             dihed_parts.kd    = imit_itp.dihedrals.kd[n];
                             dihed_parts.pn    = imit_itp.dihedrals.pn[n];
                         } else if (imit_itp.dihedrals.funct[n] == 3) {
                             dihed_parts.c0    = imit_itp.dihedrals.c0[n];
                             dihed_parts.c1    = imit_itp.dihedrals.c1[n];
                             dihed_parts.c2    = imit_itp.dihedrals.c2[n];
                             dihed_parts.c3    = imit_itp.dihedrals.c3[n];
                             dihed_parts.c4    = imit_itp.dihedrals.c4[n];
                             dihed_parts.c5    = imit_itp.dihedrals.c5[n];
                         }
                         AppendDihedral(itp, 
                                       dihed_index_replace[0], 
                                       dihed_index_replace[1], 
                                       dihed_index_replace[2], 
                                       dihed_index_replace[3], 
                                       imit_itp.dihedrals.funct[n], dihed_parts);
                         AppendPair(itp, 
                                    dihed_index_replace[0], 
                                    dihed_index_replace[3], 1);
                     }
                }
            }
        }
    }
    //
    // Between don. and acc.
    cout << endl;
    for (int l = 0; l < naa; ++l) {
        for (int m = 0; m < nda; ++m) {
            vector<vector<string>> dihedral_set;
            dihedral_set.push_back(bond_set.id2_bond[l]);
            dihedral_set.push_back(bond_set.id2_set);
            dihedral_set.push_back(bond_set.id1_set);
            dihedral_set.push_back(bond_set.id1_bond[m]);
            for (int n = 0; n < nd; ++n) {
                 if ( (bond_set.id2 == bond_set.id2_index[l]) 
                     || (bond_set.id2 == bond_set.id1_index[m]) 
                     || (bond_set.id1 == bond_set.id2_index[l]) 
                     || (bond_set.id1 == bond_set.id1_index[m]) ) {
                     continue;
                 }
                 vector<string> imit_ai = split(imit_itp.dihedrals.ai[n], '_');
                 vector<string> imit_aj = split(imit_itp.dihedrals.aj[n], '_');
                 vector<string> imit_ak = split(imit_itp.dihedrals.ak[n], '_');
                 vector<string> imit_al = split(imit_itp.dihedrals.al[n], '_');
                 vector<vector<string>>
                 dihedral_imit_set = {imit_ai, imit_aj, imit_ak, imit_al};
                 //bool osign = unordered_equal(dihedral_imit_set, dihedral_set);
                 vector<int> dihed_index = {bond_set.id2_index[l], bond_set.id2, bond_set.id1, bond_set.id1_index[m]};
                 //auto [osign, dihed_index_replace] = compare_and_reindex(dihedral_imit_set, dihedral_set, dihed_index);
                 auto [osign, dihed_index_replace] = check_and_return_index(dihedral_imit_set, dihedral_set, dihed_index);
                 if (osign) {
                 //if (dihedral_imit_set == dihedral_set) {
                      cout << "Create " 
                           << dihed_index_replace[0] << " and " 
                           << dihed_index_replace[1] << " and "
                           << dihed_index_replace[2] << " and " 
                           << dihed_index_replace[3] << " dihedral" << endl;
                      s_dihed_parts dihed_parts;
                      if ( (imit_itp.dihedrals.funct[n] == 1) 
                          || (imit_itp.dihedrals.funct[n] == 2) ) {
                          //
                          dihed_parts.ph0   = imit_itp.dihedrals.ph0[n];
                          dihed_parts.cp    = imit_itp.dihedrals.cp[n];
                          dihed_parts.mult  = imit_itp.dihedrals.mult[n];
                          dihed_parts.angle = imit_itp.dihedrals.angle[n];
                          dihed_parts.fc    = imit_itp.dihedrals.fc[n];
                      } else if ( (imit_itp.dihedrals.funct[n] == 4) 
                          || (imit_itp.dihedrals.funct[n] == 9) ) {
                          //
                          dihed_parts.phase = imit_itp.dihedrals.phase[n];
                          dihed_parts.kd    = imit_itp.dihedrals.kd[n];
                          dihed_parts.pn    = imit_itp.dihedrals.pn[n];
                      } else if (imit_itp.dihedrals.funct[n] == 3) {
                          dihed_parts.c0    = imit_itp.dihedrals.c0[n];
                          dihed_parts.c1    = imit_itp.dihedrals.c1[n];
                          dihed_parts.c2    = imit_itp.dihedrals.c2[n];
                          dihed_parts.c3    = imit_itp.dihedrals.c3[n];
                          dihed_parts.c4    = imit_itp.dihedrals.c4[n];
                          dihed_parts.c5    = imit_itp.dihedrals.c5[n];
                      }
                      //
                      AppendDihedral(itp, 
                                     dihed_index_replace[0], 
                                     dihed_index_replace[1], 
                                     dihed_index_replace[2], 
                                     dihed_index_replace[3], 
                                     imit_itp.dihedrals.funct[n], dihed_parts);
                      AppendPair(itp, 
                                 dihed_index_replace[0], 
                                 dihed_index_replace[3], 1);
                 }
            }
        }
    }
}

void RemoveDihedrals(s_itp& itp, 
                     const s_itp_imit& imit_itp, 
                     s_bond_set& bond_set) {
    cout << endl;
    cout << "REMOVE DIHEDRAL SECTION: " << endl;
    vector<int> remove_list; 
    vector<int> remove_list_id1 = ReturnRemoveDihedList(itp, imit_itp, bond_set.id1, bond_set.id1_set);
    vector<int> remove_list_id2 = ReturnRemoveDihedList(itp, imit_itp, bond_set.id2, bond_set.id2_set);
    //
    for (auto r : remove_list_id2) {
        remove_list.push_back(r);
    }
    for (auto r : remove_list_id1) {
        remove_list.push_back(r);
    }
    //
    // Remove dihedral info.
    for (auto r : remove_list) {
      int ai = itp.dihedrals.ai[r];
      int aj = itp.dihedrals.aj[r];
      int ak = itp.dihedrals.ak[r];
      int al = itp.dihedrals.al[r];
      cout <<  "Remove " 
           << ai << " and " 
           << aj << " and "
           << ak << " and " 
           << al << " dihedral" << endl;
    }
    itp = RemoveDihedrals(itp, remove_list);
}


void ChangeDihedral(s_itp& itp, 
                    const s_itp_imit& imit_itp, 
                    const vector<int>& list1, 
                    const vector<int>& list2) {
    cout << endl;
    cout << "CHANGE DIHEDRAL SECTION: " << endl;
    //
    int l1 = list1.size();
    int l2 = list2.size();
    if (l1 != l2) throw runtime_error("l1 != l2");
    //
    int nd  = itp.dihedrals.ai.size();
    int ndi = imit_itp.change_dihedrals.ai.size();
    
    vector<vector<int>> remove_index;
    vector<int>         cell_index;
    //
    vector<vector<string>>         bf_dihed_set;
    vector<vector<vector<string>>> dihed_set_list;
    
    for (int i = 0; i <  nd; ++i) {
        int ai = itp.dihedrals.ai[i];
        vector<string> ai_set = return_atom_set(itp, ai);
        int aj = itp.dihedrals.aj[i];
        vector<string> aj_set = return_atom_set(itp, aj);
        int ak = itp.dihedrals.ak[i];
        vector<string> ak_set = return_atom_set(itp, ak);
        int al = itp.dihedrals.al[i];
        vector<string> al_set = return_atom_set(itp, al);
        //
        vector<int>  dihed_index = {ai, aj, ak, al};
        vector<vector<string>> dihed_set = {ai_set, aj_set, ak_set, al_set};
        //
        for (int l = 0; l < l1; ++l) {
            int id1 = list1[l];
            int id2 = list2[l];
            //
            bool has_id1 = find(dihed_index.begin(), dihed_index.end(), id1) != dihed_index.end();
            bool has_id2 = find(dihed_index.begin(), dihed_index.end(), id2) != dihed_index.end();
            //
            if ( has_id1 && has_id2 ) continue;
            //
            if ( has_id1 || has_id2 ) {
                //
                for (int j = 0; j < ndi; ++j) {
                    vector<string> imit_ai = split(imit_itp.change_dihedrals.ai[j], '_');
                    vector<string> imit_aj = split(imit_itp.change_dihedrals.aj[j], '_');
                    vector<string> imit_ak = split(imit_itp.change_dihedrals.ak[j], '_');
                    vector<string> imit_al = split(imit_itp.change_dihedrals.al[j], '_');
                    vector<vector<string>> dihed_imit_set = {imit_ai, imit_aj, imit_ak, imit_al};
                    auto [osign, dihed_index_replace] = check_and_return_index(dihed_imit_set, dihed_set, dihed_index);
                    if (osign) {
                        if (bf_dihed_set.empty()) {
                            cell_index = {i+1};
                        } else if (dihed_set != bf_dihed_set) {
                            remove_index.push_back(cell_index);
                            dihed_set_list.push_back(bf_dihed_set);
                            cell_index = {i+1};
                        } else {
                            cell_index.push_back(i+1);
                        }
                        bf_dihed_set = dihed_set;
                        
                        cout << "Change " 
                             << dihed_index_replace[0] << " and "
                             << dihed_index_replace[1] << " and " 
                             << dihed_index_replace[2] << " and " 
                             << dihed_index_replace[3] << " dihedral ";
                        cout << "( " << imit_itp.change_dihedrals.ai[j]  
                             << " and " << imit_itp.change_dihedrals.aj[j] 
                             << " and " << imit_itp.change_dihedrals.ak[j]
                             << " and " << imit_itp.change_dihedrals.al[j] << " )" << endl;
                        break;
                    }
                }
            }
        }
    }
    
    remove_index.push_back(cell_index);
    dihed_set_list.push_back(bf_dihed_set);

    int nr = remove_index.size();
    int nds = dihed_set_list.size();
    
    if (nr != nds) throw runtime_error("nr != nds");

    int nappend = 0;

    for (int r = 0; r < nr; ++r) {
        //
        vector<int> dihed_seq = remove_index[r];
        vector<int> dihed_index;
        vector<vector<string>> dihed_set = dihed_set_list[r];
        int nc = dihed_seq.size();
        //
        cout << endl;
        cout << "Dihedral type: ";
        for (int d = 0; d < 4; ++d) {
          cout << dihed_set[d][0] << "_" << dihed_set[d][1] << " ";
        }
        cout << endl;
        cout << "# of same label diherals: " << nc << endl;
        //
        int start_label = dihed_seq[0] + nappend;
        dihed_index = {itp.dihedrals.ai[start_label-1], 
                       itp.dihedrals.aj[start_label-1], 
                       itp.dihedrals.ak[start_label-1], 
                       itp.dihedrals.al[start_label-1]};
        for (int c = 0; c < nc; ++c) {
            int label = dihed_seq[c] + nappend;
            cout << "Remove label: " << dihed_seq[c] << endl;
            itp.dihedrals.ai.erase(itp.dihedrals.ai.begin() + label - 1);
            itp.dihedrals.aj.erase(itp.dihedrals.aj.begin() + label - 1);
            itp.dihedrals.ak.erase(itp.dihedrals.ak.begin() + label - 1);
            itp.dihedrals.al.erase(itp.dihedrals.al.begin() + label - 1);
            //
            itp.dihedrals.funct.erase(itp.dihedrals.funct.begin() + label - 1);
            itp.dihedrals.ph0.erase(itp.dihedrals.ph0.begin() + label - 1);
            itp.dihedrals.cp.erase(itp.dihedrals.cp.begin() + label - 1);
            itp.dihedrals.mult.erase(itp.dihedrals.mult.begin() + label - 1);
            itp.dihedrals.fc.erase(itp.dihedrals.fc.begin() + label - 1);
            //
            itp.dihedrals.phase.erase(itp.dihedrals.phase.begin() + label - 1);
            itp.dihedrals.kd.erase(itp.dihedrals.kd.begin() + label - 1);
            itp.dihedrals.pn.erase(itp.dihedrals.pn.begin() + label - 1);
            //
            itp.dihedrals.c0.erase(itp.dihedrals.c0.begin() + label - 1);
            itp.dihedrals.c1.erase(itp.dihedrals.c1.begin() + label - 1);
            itp.dihedrals.c2.erase(itp.dihedrals.c2.begin() + label - 1);
            itp.dihedrals.c3.erase(itp.dihedrals.c3.begin() + label - 1);
            itp.dihedrals.c4.erase(itp.dihedrals.c4.begin() + label - 1);
            itp.dihedrals.c5.erase(itp.dihedrals.c5.begin() + label - 1);
        }
        
        cout << "Erase " << nc << " dihedrals" << endl;
        //
        // insert new dihedrals
        vector<int> ai_list, aj_list, ak_list, al_list, funct_list, mult_list, pn_list;
        vector<double> ph0_list, cp_list, angle_list, fc_list, phase_list, kd_list, 
                       c0_list, c1_list, c2_list, c3_list, c4_list, c5_list;
        //
        int imit_count = 0;
        for (int j = 0; j < ndi; ++j) {
            vector<string> imit_ai = split(imit_itp.change_dihedrals.ai[j], '_');
            vector<string> imit_aj = split(imit_itp.change_dihedrals.aj[j], '_');
            vector<string> imit_ak = split(imit_itp.change_dihedrals.ak[j], '_');
            vector<string> imit_al = split(imit_itp.change_dihedrals.al[j], '_');
            vector<vector<string>> dihed_imit_set = {imit_ai, imit_aj, imit_ak, imit_al};
            auto [osign, dihed_index_replace] = check_and_return_index(dihed_imit_set, dihed_set, dihed_index);
            if (osign) {
                ai_list.push_back(dihed_index_replace[0]);
                aj_list.push_back(dihed_index_replace[1]);
                ak_list.push_back(dihed_index_replace[2]);
                al_list.push_back(dihed_index_replace[3]);
                //
                funct_list.push_back(imit_itp.change_dihedrals.funct[j]);
                ph0_list.push_back(imit_itp.change_dihedrals.ph0[j]);
                cp_list.push_back(imit_itp.change_dihedrals.cp[j]);
                mult_list.push_back(imit_itp.change_dihedrals.mult[j]);
                angle_list.push_back(imit_itp.change_dihedrals.angle[j]);
                fc_list.push_back(imit_itp.change_dihedrals.fc[j]);
                //
                phase_list.push_back(imit_itp.change_dihedrals.phase[j]);
                kd_list.push_back(imit_itp.change_dihedrals.kd[j]);
                pn_list.push_back(imit_itp.change_dihedrals.pn[j]);
                //
                c0_list.push_back(imit_itp.change_dihedrals.c0[j]);
                c1_list.push_back(imit_itp.change_dihedrals.c1[j]);
                c2_list.push_back(imit_itp.change_dihedrals.c2[j]);
                c3_list.push_back(imit_itp.change_dihedrals.c3[j]);
                c4_list.push_back(imit_itp.change_dihedrals.c4[j]);
                c5_list.push_back(imit_itp.change_dihedrals.c5[j]);
                //
                imit_count += 1;
            }
        }
        //
        itp.dihedrals.ai.insert(itp.dihedrals.ai.begin() + start_label-1, ai_list.begin(), ai_list.end());
        itp.dihedrals.aj.insert(itp.dihedrals.aj.begin() + start_label-1, aj_list.begin(), aj_list.end());
        itp.dihedrals.ak.insert(itp.dihedrals.ak.begin() + start_label-1, ak_list.begin(), ak_list.end());
        itp.dihedrals.al.insert(itp.dihedrals.al.begin() + start_label-1, al_list.begin(), al_list.end());
        //
        itp.dihedrals.funct.insert(itp.dihedrals.funct.begin() + start_label-1, funct_list.begin(), funct_list.end());
        itp.dihedrals.ph0.insert(itp.dihedrals.ph0.begin() + start_label-1, ph0_list.begin(), ph0_list.end());
        itp.dihedrals.cp.insert(itp.dihedrals.cp.begin() + start_label-1, cp_list.begin(), cp_list.end());
        itp.dihedrals.mult.insert(itp.dihedrals.mult.begin() + start_label-1, mult_list.begin(), mult_list.end());
        itp.dihedrals.angle.insert(itp.dihedrals.angle.begin() + start_label-1, angle_list.begin(), angle_list.end());
        itp.dihedrals.fc.insert(itp.dihedrals.fc.begin() + start_label-1, fc_list.begin(), fc_list.end());
        itp.dihedrals.phase.insert(itp.dihedrals.phase.begin() + start_label-1, phase_list.begin(), phase_list.end());
        itp.dihedrals.kd.insert(itp.dihedrals.kd.begin() + start_label-1, kd_list.begin(), kd_list.end());
        itp.dihedrals.pn.insert(itp.dihedrals.pn.begin() + start_label-1, pn_list.begin(), pn_list.end());
        itp.dihedrals.c0.insert(itp.dihedrals.c0.begin() + start_label-1, c0_list.begin(), c0_list.end());
        itp.dihedrals.c1.insert(itp.dihedrals.c1.begin() + start_label-1, c1_list.begin(), c1_list.end());
        itp.dihedrals.c2.insert(itp.dihedrals.c2.begin() + start_label-1, c2_list.begin(), c2_list.end());
        itp.dihedrals.c3.insert(itp.dihedrals.c3.begin() + start_label-1, c3_list.begin(), c3_list.end());
        itp.dihedrals.c4.insert(itp.dihedrals.c4.begin() + start_label-1, c4_list.begin(), c4_list.end());
        itp.dihedrals.c5.insert(itp.dihedrals.c5.begin() + start_label-1, c5_list.begin(), c5_list.end());
        //
        cout << "Append " << imit_count << " dihedrals" << endl;
        nappend += imit_count - nc; // increase - decrease
    }
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////// Update Charge ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// Update single side charge
void UpdateSingleSideCharge(s_itp& connect_itp,
                            const s_charge& charge_single,
                            const vector<vector<string>>& polid_list_single, 
                            const vector<string>& resid_list_single) {
    //
    if (polid_list_single.size() != resid_list_single.size()) {
        throw runtime_error("Size mismatch between polid_list and resid_list");
    }
    int natoms = connect_itp.atoms.atom.size();
    int ns = polid_list_single.size();
    for (int i = 0; i < ns; ++i) {
        string pol_resnm = polid_list_single[i][0];
        int    pol_resid = stoi(polid_list_single[i][1]);
        string cross_resnm = resid_list_single[i];
        cout << pol_resnm << " " << pol_resid 
             << " <== " << cross_resnm << endl;
        for (int j = 0; j < natoms; ++j) {
            string atmnm = connect_itp.atoms.atom[j];
            if ( (connect_itp.atoms.resid[j] == pol_resnm) && 
                 (connect_itp.atoms.resnr[j] == pol_resid) ){
                //string atmnm = connect_itp.atoms.atom[j];
                for (int k = 0; k < charge_single.nlist; ++k) {
                    if ( (charge_single.atmnm[k] == atmnm) && 
                         (charge_single.resnm[k] == remove_digits(pol_resnm)) ) {
                        connect_itp.atoms.charge[j] = charge_single.value[k];
                    }
                }
            } else if ( (connect_itp.atoms.resid[j] == cross_resnm) ) {
                for (int k = 0; k < charge_single.nlist; ++k) {
                    if ( (charge_single.atmnm[k] == atmnm) && 
                         (charge_single.resnm[k] == remove_digits(cross_resnm)) ) {
                        connect_itp.atoms.charge[j] = charge_single.value[k];
                    }
                }
            }
        }
    }
}

//
// update double side charge
void UpdateDoubleSideCharge(s_itp& connect_itp,
                            const s_charge& charge_double,
                            const vector<vector<string>>& polid_list_double, 
                            const vector<string>& resid_list_double) {
    int nd = polid_list_double.size();
    int natoms = connect_itp.atoms.atom.size();
    for (int i = 0; i < nd; ++i) {
        string pol_resnm1 = polid_list_double[i][0];
        int    pol_resid1 = stoi(polid_list_double[i][1]);
        string pol_resnm2 = polid_list_double[i][2];
        int    pol_resid2 = stoi(polid_list_double[i][3]);
        string cross_resnm = resid_list_double[i];
        cout << pol_resnm1 << " " << pol_resid1 << " " 
             << pol_resnm2 << " " << pol_resid2 << " " 
             << "<== " <<cross_resnm << endl;
        //
        for (int j = 0; j < natoms; ++j) {
            string atmnm = connect_itp.atoms.atom[j];
            if ( (connect_itp.atoms.resid[j] == pol_resnm1) &&
                 (connect_itp.atoms.resnr[j] == pol_resid1) ) {
                for (int k = 0; k < charge_double.nlist; ++k) {
                    if ( (charge_double.atmnm[k] == atmnm) && 
                         (charge_double.resnm[k] == remove_digits(pol_resnm1)) ) {
                        connect_itp.atoms.charge[j] = charge_double.value[k];
                    }
                }
            } else if ( (connect_itp.atoms.resid[j] == pol_resnm2) &&
                         (connect_itp.atoms.resnr[j] == pol_resid2) ) {
                for (int k = 0; k < charge_double.nlist; ++k) {
                    if ( (charge_double.atmnm[k] == atmnm) && 
                         (charge_double.resnm[k] == remove_digits(pol_resnm2)) ) {
                        connect_itp.atoms.charge[j] = charge_double.value[k];
                    }
                }
            } else if ( (connect_itp.atoms.resid[j] == cross_resnm) ) {
                for (int k = 0; k < charge_double.nlist; ++k) {
                    if ( (charge_double.atmnm[k] == atmnm) && 
                         (charge_double.resnm[k] == remove_digits(cross_resnm)) ) {
                        connect_itp.atoms.charge[j] = charge_double.value[k];
                    }
                }
            }
        }
    }
}

void UpdateChargeXlinkerMonomer(s_itp& itp,
                                const s_charge& charge_x,
                                const s_charge& charge_m, 
                                const vector<int>& list_x, 
                                const vector<int>& list_m) {
    //
    int lx = list_x.size();
    int lm = list_m.size();
    if (lx != lm) throw runtime_error("lx != lm");
    //vector<double> sum_list;
    //double c_sum = 0.0;

    for (int i = 0; i < lx; ++i) {
        int idxx = list_x[i];
        int idxm = list_m[i];
        bool b_sign = false;
        //
        b_sign = is_bonded(itp, idxx, idxm);
        if (b_sign) {
            //
            // for idxx and idxm
            string residx = itp.atoms.resid[idxx-1];
            string residm = itp.atoms.resid[idxm-1];
            vector<int> same_idxx_list, same_idxm_list;
            for (int j = 0; j < itp.atoms.n; ++j) {
                string ri = itp.atoms.resid[j];
                if (ri == residx) same_idxx_list.push_back(j+1);
                if (ri == residm) same_idxm_list.push_back(j+1);
            }
            s_itp itp_bf = itp;
            ChangeCharge(itp, same_idxx_list, charge_x);
           
            //double diff_tot = 0.0;
            //for (int i = 0; i < itp.atoms.n; ++i) {
            //    double a = itp.atoms.charge[i];
            //    double b = itp_bf.atoms.charge[i];
            //    double diff = 0.0;
            //    if (a != b) diff = a - b;
            //    //if (diff < -0.05) {
            //    if (fabs(diff) > 10e-5) {
            //      cout << itp.atoms.resid[i] << ", " << itp.atoms.atom[i] << endl;
            //      cout << "a = " << a << ", " << "b = " << b << endl;
            //      cout << "Diff.: " << diff << endl;
            //    }
            //    diff_tot += diff;
            //}
            //cout << "Diff._total: " << diff_tot << endl;
            //double c_tmp = 0.0;
            //for (auto c : itp.atoms.charge) c_tmp += c;
            //cout << "Total charge: " << c_tmp << endl;
            //cout << "======================" << endl;

            //for (auto x : same_idxx_list) {
            //    string atom = itp.atoms.atom[x-1];
            //    string resid = itp.atoms.resid[x-1];
            //    string resid_trim = remove_digits(resid);
            //    for (int j = 0; j < charge_x.nlist; ++j) {
            //        if ( (charge_x.resnm[j] == resid_trim) && (charge_x.atmnm[j] == atom) ) {
            //            itp.atoms.charge[x-1] = charge_x.value[j];
            //            break;
            //        }
            //    }
            //}
            ChangeCharge(itp, same_idxm_list, charge_m);
            //for (auto m : same_idxm_list) {
            //    string atom = itp.atoms.atom[m-1];
            //    string resid = itp.atoms.resid[m-1];
            //    string resid_trim = remove_digits(resid);
            //    for (int j = 0; j < charge_m.nlist; ++j) {
            //        if ( (charge_m.resnm[j] == resid_trim) && (charge_m.atmnm[j] == atom) ) {
            //            itp.atoms.charge[m-1] = charge_m.value[j];
            //            break;
            //        }
            //    }
            //}
            
            //double x_c = 0.0;
            //for (auto same : same_idxx_list) {
            //    sum_list.push_back(same);
            //    x_c += itp.atoms.charge[same-1];
            //}
            //cout << "Charge of Xlinker: " << x_c << endl;
            //for (auto same : same_idxm_list) sum_list.push_back(same);
        } else {
            throw runtime_error("program error");
        }
    }

    //for (auto s : sum_list) c_sum += itp.atoms.charge[s-1];
    //cout << "BBBBBB: " << c_sum << endl;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Update topology ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void UpdateTopology(const vector<int>& sys_joint_list, 
                    const vector<int>& bond_joint_list, 
                    s_itp& connect_itp, 
                    s_itp_imit& imit_itp) {
    //
    int npair = bond_joint_list.size();
    for (int i = 0; i < npair; ++i) {
        //
        s_bond_set bond_set;
        //
        // Prepare grobal variables --->
        bond_set.id1 = sys_joint_list[i];
        bond_set.id2 = bond_joint_list[i];
        //
        bond_set.id1_set = return_atom_set(connect_itp, bond_set.id1);
        bond_set.id2_set = return_atom_set(connect_itp, bond_set.id2);
        //
        tie(bond_set.id1_bond, bond_set.id1_index) = GetBondedAtomandResid(connect_itp, bond_set.id1);
        tie(bond_set.id2_bond, bond_set.id2_index) = GetBondedAtomandResid(connect_itp, bond_set.id2);
        // <---

        cout << endl;
        cout << "--------------------------------------------" << endl;
        UpdateBond(connect_itp, imit_itp, bond_set); 
        UpdateAngle(connect_itp, imit_itp, bond_set);
        UpdateDihedral(connect_itp, imit_itp, bond_set);

        if (imit_itp.remove_sign) {
            //
            RemoveDihedrals(connect_itp, imit_itp, bond_set);
        }
        cout << "--------------------------------------------" << endl;
    }
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// CHANGE_TOPOLOGY ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void CHANGE_TOPOLOGY(s_itp& itp, 
                     const s_itp_imit& imit_itp,
                     const vector<int>& list1, 
                     const vector<int>& list2) {
    cout << endl;
    cout << "CHANGE TPOLOGY SECTION: " << endl;
    // Bond section
    ChangeBond(itp, imit_itp, list1, list2);
    // Angle section
    ChangeAngle(itp, imit_itp, list1, list2);
    // Dihedral section
    ChangeDihedral(itp, imit_itp, list1, list2);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// MAKE_POLYMER /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void MAKE_POLYMER(const string& pdbfile,
                  const string& outfile,
                  const string& chargefile,
                  const string& outcharge,
                  double& joint_dist,
                  vector<int>& joint_ref, 
                  vector<int>& joint_mov, 
                  int& npol, 
                  const string& captype, 
                  const double netcharge) {
    //
    PDBOPR PdbOpr;
    //
    cout << endl;
    cout << "STEP1: Read pdbfile, chargefile..." << endl;
    s_pdb monopdb = PdbOpr.LoadfrompdbAmber(pdbfile);
    vector<double> monocharge = Readaxisfile(chargefile);

    cout << endl;
    cout << "joint_charge -> 0" << endl;
    int nmono = monocharge.size();
    double charge_sum = -netcharge;
    for (int i = 0; i <  nmono; ++i) {
        if ( (i != joint_ref[1] - 1) && (i != joint_mov[1] - 1) ) {
            charge_sum += monocharge[i];
        }
    }
   
    vector<double> charge(nmono, 0.0);
    for (int i = 0; i < nmono; ++i) {
        if ( (i != joint_ref[1] - 1) && (i != joint_mov[1] - 1) ) {
            charge[i] += round((monocharge[i] - charge_sum / static_cast<double>(nmono - 2)) * 100000.0) / 100000.0;
        }
    }
    
    charge_sum = - netcharge;
    for (int i = 0; i < nmono; ++i) {
        charge_sum += charge[i];
    }

    for (int i = 0; i < nmono; ++i) {
        if ( (i != joint_ref[1] - 1) && (i != joint_mov[1] - 1) ) {
            charge[i] -= charge_sum;
        }
        //cout << charge[i] << endl;
    }

    if (boost::iequals(captype, "none")) {
        vector<double> polcharge;
        for (int i = 0; i < npol; ++i) {
            if (i == 0) {
                for (int j = 0; j < nmono; ++j) {
                    if (j != joint_ref[1] - 1) {
                        polcharge.push_back(charge[j]);
                    }
                }
            } else if (i == npol - 1) {
                for (int j = 0; j < nmono; ++j) {
                    if (j != joint_mov[1] - 1) {
                        polcharge.push_back(charge[j]);
                    }
                }
            } else {
                for (int j = 0; j < nmono; ++j) {
                    if ( (j != joint_ref[1] - 1) && (j != joint_mov[1] - 1) ) {
                        polcharge.push_back(charge[j]);
                    }
                }
            }
        }
        Only1dToFile(outcharge, polcharge);
    }

    cout << endl;
    cout << "STEP2: Start connecting..." << endl;
    monopdb = ChangeLength(monopdb, joint_dist, joint_ref, joint_mov);
    s_pdb polpdb;
    s_pdb connpdb;
    s_pdb resultpdb;

    for (int i = 0; i < npol - 1; ++i) {
        //
        if (i == 0) {
            polpdb = monopdb;
        }
        //
        s_fit fit;
        SetupFitConnect(polpdb, monopdb, joint_ref, joint_mov, fit);
        GetTrrot(fit);
        s_pdb connpdb = OperateTrrotConnect(fit, polpdb, monopdb, joint_ref, joint_mov, i+2);
        //
        joint_ref[0] += monopdb.atmnm.size() - 2;
        joint_ref[1] += monopdb.atmnm.size() - 2;
        polpdb = connpdb;
        //
        if (i == npol-2) {
            resultpdb = connpdb;
        }
    }
    //
    vector<s_pdb> resultpdb_cap = CapPolymer(resultpdb, joint_dist, 
                                             joint_ref, joint_mov, captype); 

    cout << endl;
    cout << "STTEP3: Output data to file" << endl;
    PdbOpr.OutPdbAmber(resultpdb_cap, outfile);
}

void MAKE_TOPOLOGY(const string& itpfile,
                   const string& outitp,
                   vector<int>& joint_mov, 
                   vector<int>& joint_ref, 
                   const int npol, 
                   const s_bond& bond) {
    s_itp itp;
    Readitpfile(itpfile, itp);
    s_itp  politp;
    int deln = itp.atoms.nr.size() - 1;
    
    int H2 = joint_mov[1];
    s_itp itpH2 = itp;
    RemoveIndex(itpH2, H2);
    //
    int H1 = joint_ref[1];
    s_itp itpH1 = itp;
    //
    int C2 = joint_mov[0] + deln;
    int C1 = joint_ref[0];
    //
    for (int i = 0; i <  npol-1; ++i) {
        RemoveIndex(itpH1, H1); 
        politp = CombineItp(itpH1, itpH2);
        AppendBond(politp, C1, C2, 2, bond.c0, bond.c1);
        cout << "Combine " << C1 << " and " << C2 << " bond..." << endl;
        C1 += deln - 1;
        C2 += deln - 1;
        itpH1 = politp;
        H1 += deln - 1;
    }
    Writeitp(politp, outitp);
}

void MAKE_SYS_TOP(const vector<string>& itpfile,
                  const string& outitp,
                  vector<int>& joint_mov,
                  vector<int>& joint_ref,
                  const int npol, 
                  const int Np, 
                  const int Nc, 
                  const s_bond& bond) {
    //
    s_itp monoitp;
    Readitpfile(itpfile[0], monoitp);
    s_itp  politp;
    int deln = monoitp.atoms.nr.size() - 1;
    
    int H2 = joint_mov[1];
    s_itp itpH2 = monoitp;
    RemoveIndex(itpH2, H2);
    //
    int nH2 = itpH2.atoms.resnr.size();
    for (int j = 0; j < nH2; ++j) {
        itpH2.atoms.resnr[j] = 2;
    }
    //
    int H1 = joint_ref[1];
    s_itp itpH1 = monoitp;
    //
    int C2 = joint_mov[0] + deln;
    int C1 = joint_ref[0];
    //
    for (int i = 0; i <  npol-1; ++i) {
        RemoveIndex(itpH1, H1); 
        politp = CombineItp(itpH1, itpH2);
        int nH2 = itpH2.atoms.resnr.size();
        for (int j = 0; j < nH2; ++j) {
            itpH2.atoms.resnr[j] = i + 3;
        }
        AppendBond(politp, C1, C2, 2, bond.c0, bond.c1);
        cout << "Combine " << C1 << " and " << C2 << " bond..." << endl;
        C1 += deln - 1;
        C2 += deln - 1;
        itpH1 = politp;
        H1 += deln - 1;
    }
    //
    s_itp sysitp = politp;
    for (int i = 0; i < Np-1; ++i) {
        int polsize = politp.atoms.resnr.size();
        for (int j = 0; j < polsize; ++j) {
             politp.atoms.resid[j] = replaceNumber(politp.atoms.resid[j], i+2);
        }
        sysitp = CombineItp(sysitp, politp);
    }

    //
    // Read itp info. of cross.
    s_itp crossitp;
    Readitpfile(itpfile[1], crossitp);
    cout << itpfile[1] << endl;

    if (Nc > 0) {
        for (int i = 0; i < Nc; ++i) {
            int crosssize = crossitp.atoms.resnr.size();
            for (int j = 0; j < crosssize; ++j) {
                 crossitp.atoms.resid[j] = replaceNumber(crossitp.atoms.resid[j], i+1);
            }
            sysitp = CombineItp(sysitp, crossitp);
        }
    } else {
        cout << "Unset crosslinking material..." << endl;
    }

    Writeitp(sysitp, outitp);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////// COMBINE POLYMER and CROSSLINKING MATER. /////////////////////
///////////////////////////////////////////////////////////////////////////////
void COMBINE_POL_CROSS(const vector<string>& itpfile, 
                       const string& outitp,
                       const int Np, 
                       const int Nc) {
    //
    s_itp politp;
    s_itp crossitp;
    Readitpfile(itpfile[0], politp);
    Readitpfile(itpfile[1], crossitp);
    s_itp combineitp = politp;
    int npol = politp.atoms.resid.size();
    int ncross = crossitp.atoms.resid.size();
    //
    for (int i = 0; i < Np-1; ++i) {
        for (int j = 0; j < npol; ++j) {
            politp.atoms.resid[j] = replaceNumber(politp.atoms.resid[j], i+2);
        }
        combineitp = CombineItp(combineitp, politp);
    }
    for (int i = 0; i < Nc; ++i) {
        combineitp = CombineItp(combineitp, crossitp); 
        for (int j = 0; j < ncross; ++j) {
            crossitp.atoms.resid[j] = replaceNumber(crossitp.atoms.resid[j], i+2);
        }
    }
    Writeitp(combineitp, outitp);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////// Crosslinking analysis ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void MAKE_CROSSLINKING(const vector<string>& itpfile, 
                       const string& trajfile,
                       const string& outitp,
                       const vector<string>& chargefile,
                       const vector<string>& selpol,
                       const vector<string>& selcross,
                       const vector<string>& joint_string,
                       const int Np, 
                       const int Nc, 
                       const double rc, 
                       const string rule) {
    cout << endl;
    cout << "STEP1: Read itpfile and trajfile..." << endl;
    s_itp itp;
    Readitpfile(itpfile[0], itp);
    //
    cout << endl; 
    cout << "Use last frame..." << endl;
    vector<vector<double>> lastcoord; 
    vector<double> lastbox;
    if (boost::iequals(get_extension(trajfile), "xtc")) {
        //
        s_xtc xtc = ReadXTC(trajfile);
        int nlast = xtc.nframes - 1;
        lastcoord = xtc.coord[nlast];
        lastbox.push_back(xtc.box[nlast][0]);
        lastbox.push_back(xtc.box[nlast][4]);
        lastbox.push_back(xtc.box[nlast][8]);
    } else if (boost::iequals(get_extension(trajfile), "nc")) {
        //
        vector<vector<vector<double>>> coord = GetCoordfromNC(trajfile);
        int nlast = coord.size() - 1;
        lastcoord = coord[nlast];
        lastbox = {0.0, 0.0, 0.0};
    }
    //
    cout << endl;
    cout << "STEP2: Creating bonds of crosslinking mater. and polymers..." << endl;
    
    //--------------------------------------------------------------------------
    // For polymer part...
    int natoms = itp.atoms.resid.size();
    //
    vector<string> selpol_split = split_reverse(selpol[0]);
    vector<string> remove_split = split_reverse(selpol[1]);
    if (selpol_split[0] != remove_split[0]) 
        throw runtime_error("kind of selpol_split != kind of remove_split");
    
    vector<string> pr_tmp = {selpol_split[0], selpol_split[1], remove_split[1]};

    vector<vector<int>> selpol_list = ExtractList(itp, Np, natoms, pr_tmp);
   
    cout << endl;
    cout << "Polymer select: " << endl;
    int pol_label = 1;
    for (vector<int> sel_list : selpol_list) {
        cout << pr_tmp[0] << pol_label << ": ";
        for (int sel : sel_list) {
            cout << sel << " ";
        }
        cout << endl;
        pol_label += 1;
    }
    int npol = selpol_list[0].size() / 2;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // For crosslinking material part...
    vector<string> selcross_split1 = split_reverse(selcross[0]);
    vector<string> selcross_split2 = split_reverse(selcross[1]);
    if (selcross_split1[0] != selcross_split2[0]) 
        throw runtime_error("kind of selcross_split1 != kind of selcross_split2");
    vector<string> c_tmp = {selcross_split1[0], selcross_split1[1], selcross_split2[1]};
    vector<vector<int>> selcross_list = ExtractList(itp, Nc, natoms, c_tmp);
    cout << endl;
    cout << "Crosslinking material select: " << endl;
    int cross_label = 1;
    for (vector<int> sel_list : selcross_list) {
        cout << c_tmp[0] << cross_label << ": ";
        for (int sel : sel_list) {
            cout << sel << " ";
        }
        cout << endl;
        cross_label += 1;
    }
    //--------------------------------------------------------------------------

    // cout << endl; 
    // cout << "Use last frame: " << xtc.nframes << endl;
    // vector<array<double, 3>> lastcoord = xtc.coord[xtc.nframes-1];
    // vector<vector<double>> lastcoord; 
    
    vector<vector<vector<int>>> connect_list;
    for (int i = 0; i < Np; ++i) {
        vector<vector<int>> connect_pol;
        for (int j = 0; j < npol; ++j) {
            vector<int> connect_cell;
            int joint_pol = selpol_list[i][2*j+1] - 1;
            //
            vector<double> joint_crd = lastcoord[joint_pol];
            //
            for (int k = 0; k < Nc; ++k) {
                int c1 = selcross_list[k][0] - 1;
                int c2 = selcross_list[k][1] - 1;
                //
                vector<double> c1_crd = lastcoord[c1];
                vector<double> c2_crd = lastcoord[c2];
                double length_jc1 = dist2_pbc(c1_crd, joint_crd, lastbox);
                length_jc1 = sqrt(length_jc1);
                double length_jc2 = dist2_pbc(c2_crd, joint_crd, lastbox);
                length_jc2 = sqrt(length_jc2);
                if (length_jc1 < rc) {
                    connect_cell.push_back(c1+1);
                }
                //
                if (length_jc2 < rc) {
                    connect_cell.push_back(c2+1);
                }
            }
            connect_pol.push_back(connect_cell);
        }
        connect_list.push_back(connect_pol);
    }

    // 
    if (boost::iequals(rule, "random")) {
        cout << endl;
        for (int i = 0; i < Np; ++i) {
            for (int j = 0; j < npol; ++j) {
               int ns = connect_list[i][j].size();
               //cout << ns << endl;
               if (ns > 1) {
                   cout << pr_tmp[0] << Np << " (resid " << j+1 << ") is connected ";
                   for (int k = 0; k < ns; ++k) {
                       cout << connect_list[i][j][k] << " ";
                   }
                   cout << endl;
                   cout << "==>  Pick up 1-component ";
                   connect_list[i][j] = pick_one_random(connect_list[i][j]);
                   cout << connect_list[i][j][0] << endl;
               }
            }
        }
    } else if (boost::iequals(rule, "min")) {
        cout << endl;
        for (int i = 0; i < Np; ++i) {
            for (int j = 0; j < npol; ++j) {
               int joint_pol = selpol_list[i][2*j+1] - 1;
               int ns = connect_list[i][j].size();
               if (ns > 1) {
                   cout << pr_tmp[0] << i+1 << " (resid " << j+1 << ") is connected ";
                   for (int k = 0; k < ns; ++k) {
                       cout << connect_list[i][j][k] << " ";
                   }
                   double length_min = 100.0;
                   int idx_min = 0;
                   for (int k = 0; k < ns; ++k) {
                      int idx = connect_list[i][j][k];
                      double length = dist2_pbc(lastcoord[idx-1], lastcoord[joint_pol], lastbox);
                      length = sqrt(length);
                      if (length < length_min) {
                          length_min = length;
                          idx_min = idx;
                      }
                   }
                   connect_list[i][j] = filter_to_one(connect_list[i][j], idx_min); 
                   cout << endl;
                   cout << "==> Pick up 1-component ";
                   cout << connect_list[i][j][0] << endl;
               }
            }
        }
    } else if (boost::iequals(rule, "non-self")) {
        vector<vector<vector<int>>> connect_list_tmp;
        cout << endl;
        for (int i = 0; i < Np; ++i) {
            vector<vector<int>> cell_vv_tmp;
            string polname = combine_string_int(pr_tmp[0], i+1);
            for (int j = 0; j < npol; ++j) {
                vector<int> cell_v_tmp;
                for (auto c : connect_list[i][j]) {
                    vector<string> boneded_resids_c = GetBondedResidues(itp, c);
                    int c_pair = -1;
                    for (int k = 0; k < Nc; ++k) {
                        vector<int> x_tmp = selcross_list[k];
                        if (c == x_tmp[0]) {
                            c_pair = x_tmp[1];
                            break;
                        } else if (c == x_tmp[1]) {
                            c_pair = x_tmp[0];
                            break;
                        }
                    }
                    //
                    if (c_pair == -1) throw runtime_error("Xlinked atom does not have pair atom...");
                    //
                    vector<string> boneded_resids_cp = GetBondedResidues(itp, c_pair);
                    string r_c;
                    for (auto r_c_tmp : boneded_resids_c) {
                        string r_c_trim = remove_digits(r_c_tmp);
                        if ( r_c_trim == pr_tmp[0] ) {
                            r_c = r_c_tmp;
                            break;
                        }
                    }

                    if ( !r_c.empty() ) break;
                    //
                    bool self_sign = false;
                    for (auto r_cp_tmp : boneded_resids_cp) {
                        //cout << r_cp_tmp << " ";
                        if ( r_cp_tmp == polname ) {
                            self_sign = true;
                        }
                    }
                    //
                    cout << endl;
                    if (!self_sign) {
                        cell_v_tmp.push_back(c);
                    } else {
                        break;
                    }
                }
                cell_vv_tmp.push_back(cell_v_tmp);
            }
            connect_list_tmp.push_back(cell_vv_tmp);
        }
        connect_list = connect_list_tmp;
        //
        cout << endl;
        for (int i = 0; i < Np; ++i) {
            for (int j = 0; j < npol; ++j) {
               int joint_pol = selpol_list[i][2*j+1] - 1;
               int ns = connect_list[i][j].size();
               if (ns > 1) {
                   cout << pr_tmp[0] << i+1 << " (resid " << j+1 << ") is connected ";
                   for (int k = 0; k < ns; ++k) {
                       cout << connect_list[i][j][k] << " ";
                   }
                   double length_min = 100.0;
                   int idx_min = 0;
                   for (int k = 0; k < ns; ++k) {
                      int idx = connect_list[i][j][k];
                      double length = dist2_pbc(lastcoord[idx-1], lastcoord[joint_pol], lastbox);
                      length = sqrt(length);
                      if (length < length_min) {
                          length_min = length;
                          idx_min = idx;
                      }
                   }
                   connect_list[i][j] = filter_to_one(connect_list[i][j], idx_min); 
                   cout << endl;
                   cout << "==> Pick up 1-component ";
                   cout << connect_list[i][j][0] << endl;
               }
            }
        }
    }

    //
    //
    cout << endl;
    cout << "Read itpimitfile..." << endl;
    s_itp_imit itp_imit;
    Readitpimitfile(itpfile[1], itp_imit);

    s_itp connect_itp = itp;
    cout << endl;
    for (int i = 0; i < Np; ++i) {
        for (int j = 0; j < npol; ++j) {
             int ncross_connect = connect_list[i][j].size();
             //cout << ncross_connect << endl;
             if (ncross_connect > 0) {
                 for (int k = 0; k < ncross_connect; ++k) {
                     //
                     int cross_id = connect_list[i][j][k];
                     int remove_id = selpol_list[i][2*j];
                     int don_id = selpol_list[i][2*j+1];
                     vector<string> bonded_resids = GetBondedResidues(connect_itp, cross_id);
                     //for (auto b : bonded_resids) {
                     //    cout << b << " ";
                     //}
                     //cout << endl;
                     bool csign = has_bonded_resid(bonded_resids, pr_tmp[0]);
                     bool bsign = is_bonded(connect_itp, remove_id, don_id);
                     //cout << csign << " " << bsign <<  selpol[0] << endl; 
                     //
                     //
                     if (csign == true || bsign == false) {
                          cout << "--------------------------------------------" << endl;
                          cout << cross_id << " and " << pr_tmp[0] 
                               << " has already bond, skipping..." << endl;
                          cout << "--------------------------------------------" << endl;
                     } else {
                          //
                          string don_resid       = itp.atoms.resid[don_id-1];
                          string don_resid_trim  = remove_digits(don_resid);
                          string don_atom        = itp.atoms.atom[don_id-1];
                          vector<string> don_set = {don_atom, don_resid_trim};
                          //
                          string cross_resid       = itp.atoms.resid[cross_id - 1];
                          string cross_resid_trim  = remove_digits(cross_resid);
                          string cross_atom        = itp.atoms.atom[cross_id - 1];
                          vector<string> cross_set = {cross_atom, cross_resid_trim};
                          //
                          cout << endl;
                          cout << "--------------------------------------------" << endl;
                          cout << "BOND SECTION: " << endl;
                          cout << "Search bond info..." << endl;
                          int nb = itp_imit.bonds.ai.size();
                          double bond_c0 = 0.0, bond_c1 = 0.0;
                          int funct = 1;
                          for (int l = 0; l < nb; ++l) {
                              vector<string> imit_ai = split(itp_imit.bonds.ai[l], '_');
                              vector<string> imit_aj = split(itp_imit.bonds.aj[l], '_');
                              if ( ( imit_ai == don_set && imit_aj == cross_set ) 
                                   || ( imit_ai == cross_set && imit_aj == don_set ) ) {
                                  funct = itp_imit.bonds.funct[l];
                                  bond_c0 = itp_imit.bonds.c0[l];
                                  bond_c1 = itp_imit.bonds.c1[l];
                              }
                          }
                          cout << cross_id << " and " << don_id 
                               << " does not have bond, connecting..." << endl;
                          AppendBond(connect_itp, don_id, cross_id, funct, bond_c0, bond_c1);
                          cout << "Remove " << remove_id << " conections..." << endl;
                          RemoveConnections(connect_itp, remove_id);
                          //connect_itp.atoms.resid[remove_id-1] = "R" + to_string(remove_id);
                          //connect_itp.atoms.resid[remove_id-1] = "R";
                          //
                          cout << endl;
                          cout << "ANGLE SECTION:" << endl;
                          int na = itp_imit.angles.ai.size();
                          cout << "bond check" << endl;
                          auto [don_bond, don_index] = GetBondedAtomandResid(connect_itp, don_id);
                          cout << don_id << " "
                               << itp.atoms.atom[don_id-1] << ": " << endl;
                          for (auto vec : don_bond) {
                              for (auto dumm : vec) {
                                  cout << dumm << " ";
                              }
                              cout << endl;
                          }
                          
                          auto [cross_bond, cross_index] = GetBondedAtomandResid(connect_itp, cross_id);
                          cout << cross_id << " " << itp.atoms.atom[cross_id-1] << ": " << endl;
                          for (auto vec : cross_bond) {
                              for (auto dumm : vec) {
                                  cout << dumm << " ";
                              }
                              cout << endl;
                          }
                          cout << endl;
                          //
                          // Creating angles...
                          //
                          // don. side
                          int nda = don_bond.size() - 1;
                          for (int l = 0; l < nda; ++l) {
                               vector<vector<string>> angle_set;
                               vector<vector<string>> angle_imit_set;
                               don_bond[l][1] = remove_digits(don_bond[l][1]);
                               //
                               angle_set.push_back(don_bond[l]);
                               angle_set.push_back(don_set);
                               angle_set.push_back(cross_set);
                               //
                               for (int m = 0; m < na; ++m) {
                                   vector<string> imit_ai = split(itp_imit.angles.ai[m], '_');
                                   vector<string> imit_aj = split(itp_imit.angles.aj[m], '_');
                                   vector<string> imit_ak = split(itp_imit.angles.ak[m], '_');
                                   angle_imit_set = {imit_ai, imit_aj, imit_ak};
                                   //bool osign = unordered_equal(angle_imit_set, angle_set);
                                   vector<int> angle_index = {don_index[l], don_id, cross_id};
                                   //for (auto v : angle_index) {
                                   //    cout << v << endl;
                                   //}
                                   auto [osign, angle_index_replace] = compare_and_reindex(angle_imit_set, angle_set, angle_index);
                                   //if (angle_imit_set == angle_set) {
                                   if (osign) {
                                       cout << "Create angle: " << angle_index_replace[0] << " "
                                            << angle_index_replace[1] << " " << angle_index_replace[2] << endl;
                                       AppendAngle(connect_itp, angle_index_replace[0], angle_index_replace[1], 
                                                   angle_index_replace[2], itp_imit.angles.funct[m], 
                                                   itp_imit.angles.angle[m], itp_imit.angles.fc[m]);
                                   }
                               }
                          }
                          //
                          // cross. side
                          int naa = cross_bond.size() - 1;
                          for (int l = 0; l < naa; ++l) {
                               vector<vector<string>> angle_set;
                               vector<vector<string>> angle_imit_set;
                               cross_bond[l][1] = remove_digits(cross_bond[l][1]);
                               //
                               angle_set.push_back(cross_bond[l]);
                               angle_set.push_back(cross_set);
                               angle_set.push_back(don_set);
                               //
                               for (int m = 0; m < na; ++m) {
                                   vector<string> imit_ai = split(itp_imit.angles.ai[m], '_');
                                   vector<string> imit_aj = split(itp_imit.angles.aj[m], '_');
                                   vector<string> imit_ak = split(itp_imit.angles.ak[m], '_');
                                   angle_imit_set = {imit_ai, imit_aj, imit_ak};
                                   //bool osign = unordered_equal(angle_imit_set, angle_set);
                                   //if (angle_imit_set == angle_set) {
                                   vector<int> angle_index = {cross_index[l], cross_id, don_id};
                                   auto [osign, angle_index_replace] = compare_and_reindex(angle_imit_set, angle_set, angle_index);
                                   if (osign) {
                                       cout << "Create angle: " << angle_index_replace[0] << " "
                                            << angle_index_replace[1] << " " << angle_index_replace[2] << endl;
                                       AppendAngle(connect_itp, angle_index_replace[0], angle_index_replace[1], 
                                                   angle_index_replace[2], itp_imit.angles.funct[m], 
                                                   itp_imit.angles.angle[m], itp_imit.angles.fc[m]);
                                   }
                               }
                          }
                          //
                          cout << endl;
                          cout << "DIHEDRAL SECTION:" << endl;
                          int nd = itp_imit.dihedrals.ai.size();
                          //
                          // don. side
                          for (int l = 0; l < nda; ++l) {
                              //cout << don_index[l] << ": ";
                              auto [don_bond2, don_index2] = GetBondedAtomandResid(connect_itp, don_index[l]);
                              int nb2 = don_index2.size();
                              //for (auto v : don_index2) {
                              //    cout << v << " ";
                              //}
                              if (nb2 == 1) continue;

                              for (int m = 0; m < nb2; ++m) {
                                  don_bond2[m][1] = remove_digits(don_bond2[m][1]);
                                  vector<vector<string>> dihedral_set;
                                  dihedral_set.push_back(don_bond2[m]);
                                  dihedral_set.push_back(don_bond[l]);
                                  dihedral_set.push_back(don_set);
                                  dihedral_set.push_back(cross_set);
                                  if (don_bond2[m][0] != pr_tmp[1]) {
                                      //cout << "Create " << don_index2[m] << " and " 
                                      //     << don_index[l] << " and " << selpol_list[i][2*j+1] 
                                      //     << " and " << cross_id << endl;
                                     vector<int> dihed_index = {don_index2[m], don_index[l], don_id, cross_id};
                                     //for (auto v : dihed_index) {
                                     //    cout << v << " ";
                                     //}
                                     for (int n = 0; n < nd; ++n) {
                                         vector<string> imit_ai = split(itp_imit.dihedrals.ai[n], '_');
                                         vector<string> imit_aj = split(itp_imit.dihedrals.aj[n], '_');
                                         vector<string> imit_ak = split(itp_imit.dihedrals.ak[n], '_');
                                         vector<string> imit_al = split(itp_imit.dihedrals.al[n], '_');
                                         vector<vector<string>> 
                                         dihedral_imit_set = {imit_ai, imit_aj, imit_ak, imit_al};
                                         //cout << "---------" << endl;
                                         //cout << imit_ai[0] << " " << imit_ai[1] 
                                         //     <<  " " << imit_aj[0] << " " << imit_aj[1] 
                                         //     << " " << imit_ak[0] << " "  << imit_ak[1] 
                                         //     << " " << imit_al[0] << " "  << imit_al[1] << endl;
                                         //cout << dihedral_set[0][0] << " " << dihedral_set[0][1] 
                                         //     << " " << dihedral_set[1][0] << " " << dihedral_set[1][1] 
                                         //     << " " << dihedral_set[2][0] << " " << dihedral_set[2][1]
                                         //     << " " << dihedral_set[3][0] << " " << dihedral_set[3][1] << endl;
                                         //cout << "---------" << endl;
                                         //bool osign = unordered_equal(dihedral_imit_set, dihedral_set);
                                         //if (dihedral_imit_set == dihedral_set) {
                                         //vector<int> dihed_index = {don_index2[m], don_index[l], don_id, cross_id};
                                         //for (auto v : dihed_index) {
                                         //    cout << v << " ";
                                         //}
                                         //cout << endl;
                                         auto [osign, dihed_index_replace] = compare_and_reindex(dihedral_imit_set, dihedral_set, dihed_index);
                                         if (osign) {
                                              cout << "Create " 
                                                   << dihed_index_replace[0] << " and " 
                                                   << dihed_index_replace[1] << " and " 
                                                   << dihed_index_replace[2] << " and " 
                                                   << dihed_index_replace[3] << " dihedral" << endl;
                                              s_dihed_parts dihed_parts;
                                              if ( (itp_imit.dihedrals.funct[n] == 1) 
                                                  || (itp_imit.dihedrals.funct[n] == 2) ) {
                                                  //
                                                  dihed_parts.ph0   = itp_imit.dihedrals.ph0[n];
                                                  dihed_parts.cp    = itp_imit.dihedrals.cp[n];
                                                  dihed_parts.mult  = itp_imit.dihedrals.mult[n];
                                                  dihed_parts.angle = itp_imit.dihedrals.angle[n];
                                                  dihed_parts.fc    = itp_imit.dihedrals.fc[n];
                                              } else if ( (itp_imit.dihedrals.funct[n] == 4) 
                                                  || (itp_imit.dihedrals.funct[n] == 9) ) {
                                                  //
                                                  dihed_parts.phase = itp_imit.dihedrals.phase[n];
                                                  dihed_parts.kd    = itp_imit.dihedrals.kd[n];
                                                  dihed_parts.pn    = itp_imit.dihedrals.pn[n];
                                              }
                                              AppendDihedral(connect_itp, 
                                                             dihed_index_replace[0], 
                                                             dihed_index_replace[1], 
                                                             dihed_index_replace[2], 
                                                             dihed_index_replace[3], 
                                                             itp_imit.dihedrals.funct[n],
                                                             dihed_parts);
                                               AppendPair(connect_itp, 
                                                          dihed_index_replace[0], 
                                                          dihed_index_replace[3], 1);
                                         }       
                                     }           
                                  }              
                              }                  
                          }                      
                          //                     
                          // cross. side
                          cout << endl;
                          for (int l = 0; l < naa; ++l) {
                              auto [cross_bond2, cross_index2] = GetBondedAtomandResid(connect_itp, cross_index[l]);
                              int nb2 = cross_index2.size();
                              for (int m = 0; m < nb2; ++m) {
                                  cross_bond2[m][1] = remove_digits(cross_bond2[m][1]);
                                  //cout << cross_bond2[m][0] << " " << cross_bond2[m][1] << endl;
                                  vector<vector<string>> dihedral_set;
                                  dihedral_set.push_back(cross_bond2[m]);
                                  dihedral_set.push_back(cross_bond[l]);
                                  dihedral_set.push_back(cross_set);
                                  dihedral_set.push_back(don_set);
                                  if ( cross_bond2[m][0] != c_tmp[1]
                                       && cross_bond2[m][0] != c_tmp[2] ) {
                                      //cout << "Create " << cross_index2[m] << " and "
                                      //     << cross_index[l] << " and " << selpol_list[i][2*j+1]
                                      //     << " and " << cross_id << endl;
                                      for (int n = 0; n < nd; ++n) {
                                          vector<string> imit_ai = split(itp_imit.dihedrals.ai[n], '_');
                                          vector<string> imit_aj = split(itp_imit.dihedrals.aj[n], '_');
                                          vector<string> imit_ak = split(itp_imit.dihedrals.ak[n], '_');
                                          vector<string> imit_al = split(itp_imit.dihedrals.al[n], '_');
                                          vector<vector<string>>
                                          dihedral_imit_set = {imit_ai, imit_aj, imit_ak, imit_al};
                                          //cout << "---------" << endl;
                                          //cout << imit_ai[0] << " " << imit_ai[1] 
                                          //     <<  " " << imit_aj[0] << " " << imit_aj[1] 
                                          //     << " " << imit_ak[0] << " "  << imit_ak[1] 
                                          //     << " " << imit_al[0] << " "  << imit_al[1] << endl;
                                          //cout << dihedral_set[0][0] << " " << dihedral_set[0][1] 
                                          //     << " " << dihedral_set[1][0] << " " << dihedral_set[1][1] 
                                          //     << " " << dihedral_set[2][0] << " " << dihedral_set[2][1]
                                          //     << " " << dihedral_set[3][0] << " " << dihedral_set[3][1] << endl;
                                          //cout << "---------" << endl;
                                          //bool osign = unordered_equal(dihedral_imit_set, dihedral_set);
                                          //if (dihedral_imit_set == dihedral_set) {
                                          vector<int> dihed_index = {cross_index2[m], cross_index[l], cross_id, don_id};
                                          auto [osign, dihed_index_replace] = compare_and_reindex(dihedral_imit_set, dihedral_set, dihed_index);
                                          if (osign) {
                                               cout << "Create " 
                                                    << dihed_index_replace[0] << " and "
                                                    << dihed_index_replace[1] << " and " 
                                                    << dihed_index_replace[2] << " and " 
                                                    << dihed_index_replace[3] << " dihedral" << endl;
                                               s_dihed_parts dihed_parts;
                                               if ( (itp_imit.dihedrals.funct[n] == 1) 
                                                   || (itp_imit.dihedrals.funct[n] == 2) ) {
                                                   //
                                                   dihed_parts.ph0   = itp_imit.dihedrals.ph0[n];
                                                   dihed_parts.cp    = itp_imit.dihedrals.cp[n];
                                                   dihed_parts.mult  = itp_imit.dihedrals.mult[n];
                                                   dihed_parts.angle = itp_imit.dihedrals.angle[n];
                                                   dihed_parts.fc    = itp_imit.dihedrals.fc[n];
                                               } else if ( (itp_imit.dihedrals.funct[n] == 4) 
                                                   || (itp_imit.dihedrals.funct[n] == 9) ) {
                                                   //
                                                   dihed_parts.phase = itp_imit.dihedrals.phase[n];
                                                   dihed_parts.kd    = itp_imit.dihedrals.kd[n];
                                                   dihed_parts.pn    = itp_imit.dihedrals.pn[n];
                                               }
                                               AppendDihedral(connect_itp, 
                                                             dihed_index_replace[0], 
                                                             dihed_index_replace[1], 
                                                             dihed_index_replace[2], 
                                                             dihed_index_replace[3], 
                                                             itp_imit.dihedrals.funct[n], dihed_parts);
                                               AppendPair(connect_itp, 
                                                          dihed_index_replace[0], 
                                                          dihed_index_replace[3], 1);
                                           }
                                      }
                                  }
                              }
                          }
                          //
                          // Between don. and cross.
                          cout << endl;
                          for (int l = 0; l < naa; ++l) {
                              for (int m = 0; m < nda; ++m) {
                                  //cout << "aCreate " << cross_index[l] << " and " << cross_id << " and "
                                  //     << don_index[m] << " and " << selpol_list[i][2*j+1] << endl;
                                  vector<vector<string>> dihedral_set;
                                  dihedral_set.push_back(cross_bond[l]);
                                  dihedral_set.push_back(cross_set);
                                  dihedral_set.push_back(don_set);
                                  dihedral_set.push_back(don_bond[m]);
                                  for (int n = 0; n < nd; ++n) {
                                       vector<string> imit_ai = split(itp_imit.dihedrals.ai[n], '_');
                                       vector<string> imit_aj = split(itp_imit.dihedrals.aj[n], '_');
                                       vector<string> imit_ak = split(itp_imit.dihedrals.ak[n], '_');
                                       vector<string> imit_al = split(itp_imit.dihedrals.al[n], '_');
                                       vector<vector<string>>
                                       dihedral_imit_set = {imit_ai, imit_aj, imit_ak, imit_al};
                                       //bool osign = unordered_equal(dihedral_imit_set, dihedral_set);
                                       vector<int> dihed_index = {cross_index[l], cross_id, don_id, don_index[m]};
                                       auto [osign, dihed_index_replace] = compare_and_reindex(dihedral_imit_set, dihedral_set, dihed_index);
                                       if (osign) {
                                       //if (dihedral_imit_set == dihedral_set) {
                                            cout << "Create " 
                                                 << cross_index[l] << " and " 
                                                 << cross_id       << " and "
                                                 << don_id         << " and " 
                                                 << don_index[m]   << " dihedral" << endl;
                                            s_dihed_parts dihed_parts;
                                            if ( (itp_imit.dihedrals.funct[n] == 1) 
                                                || (itp_imit.dihedrals.funct[n] == 2) ) {
                                                //
                                                dihed_parts.ph0   = itp_imit.dihedrals.ph0[n];
                                                dihed_parts.cp    = itp_imit.dihedrals.cp[n];
                                                dihed_parts.mult  = itp_imit.dihedrals.mult[n];
                                                dihed_parts.angle = itp_imit.dihedrals.angle[n];
                                                dihed_parts.fc    = itp_imit.dihedrals.fc[n];
                                            } else if ( (itp_imit.dihedrals.funct[n] == 4) 
                                                || (itp_imit.dihedrals.funct[n] == 9) ) {
                                                //
                                                dihed_parts.phase = itp_imit.dihedrals.phase[n];
                                                dihed_parts.kd    = itp_imit.dihedrals.kd[n];
                                                dihed_parts.pn    = itp_imit.dihedrals.pn[n];
                                            }
                                            //
                                            AppendDihedral(connect_itp, 
                                                          dihed_index_replace[0], 
                                                          dihed_index_replace[1], 
                                                          dihed_index_replace[2], 
                                                          dihed_index_replace[3], 
                                                          itp_imit.dihedrals.funct[n], dihed_parts);
                                            AppendPair(connect_itp, 
                                                       dihed_index_replace[0], 
                                                       dihed_index_replace[3], 1);
                                       }
                                  }
                              }
                          }
                          cout << "--------------------------------------------" << endl;
                     }
                 }
             }
        }
    }

    cout << endl;
    cout << "Update charge..."<< endl;
    //
    s_joint joint;
    for (string j : joint_string) {
        vector<string> name_split = split(j);
        joint.atmnm.push_back(name_split[0]);
        joint.resnm.push_back(name_split[1]);
    }

    //
    // Case single cross-linking bond
     
    //
    // for 1-side
    s_charge charge_single1;
    auto c1_temp = ReadChargefile(chargefile[0]);
    charge_single1.name = move(c1_temp.first);
    charge_single1.value = move(c1_temp.second);
    charge_single1.nlist = charge_single1.name.size();
    for (int i = 0; i < charge_single1.nlist; ++i) {
        vector<string> name_split = split(charge_single1.name[i]);
        charge_single1.atmnm.push_back(name_split[0]);
        charge_single1.resnm.push_back(name_split[1]);
    }
    //
    double net_charge_single1 = 0.0;
    double sum_charge_single1 = 0.0;
    for (int i = 0; i < charge_single1.nlist; ++i) {
        if ( (charge_single1.resnm[i] == joint.resnm[0]) 
              && (charge_single1.atmnm[i] == joint.atmnm[0]) ) {
            continue;
        } else if ( (charge_single1.resnm[i] == joint.resnm[1]) 
                 && (charge_single1.atmnm[i] == joint.atmnm[1]) ) {
            continue;
        } else if ( (charge_single1.resnm[i] == pr_tmp[0])
                 && (charge_single1.atmnm[i] == pr_tmp[2]) ) {
            net_charge_single1 = - charge_single1.value[i];
        } else {
            sum_charge_single1 += charge_single1.value[i];
        }
    }

    net_charge_single1 -= sum_charge_single1;
    //
    for (int i = 0; i < charge_single1.nlist; ++i) {
        if ( (charge_single1.resnm[i] == joint.resnm[0]) 
              && (charge_single1.atmnm[i] == joint.atmnm[0]) ) {
            charge_single1.value[i] = 0.0;
        } else if ( (charge_single1.resnm[i] == joint.resnm[1]) 
                 && (charge_single1.atmnm[i] == joint.atmnm[1]) ) {
            charge_single1.value[i] = 0.0;
        } else if ( (charge_single1.resnm[i] == pr_tmp[0])
                 && (charge_single1.atmnm[i] == pr_tmp[2]) ) {
            continue;
        } else {
            //charge_single1.value[i] += net_charge_single1 / (charge_single1.nlist - 2);
            charge_single1.value[i] 
            = round( (charge_single1.value[i] + net_charge_single1 / static_cast<double>(charge_single1.nlist - 3) ) * 10000000.0 ) / 10000000.0;
        }
    }

    //double cs1 = 0.0;
    //for (auto c : charge_single1.value) {
    //    cs1 += c;
    //}
    //cout << "Total sytem charge: " << cs1 << endl;

    //
    // for 2-side
    s_charge charge_single2;
    auto c2_temp = ReadChargefile(chargefile[1]);
    charge_single2.name = move(c2_temp.first);
    charge_single2.value = move(c2_temp.second);
    charge_single2.nlist = charge_single2.name.size();
    for (int i = 0; i < charge_single2.nlist; ++i) {
        vector<string> name_split = split(charge_single2.name[i]);
        charge_single2.atmnm.push_back(name_split[0]);
        charge_single2.resnm.push_back(name_split[1]);
    }
    //
    double net_charge_single2 = 0.0;
    double sum_charge_single2 = 0.0;
    for (int i = 0; i < charge_single2.nlist; ++i) {
        if ( (charge_single2.resnm[i] == joint.resnm[0]) 
              && (charge_single2.atmnm[i] == joint.atmnm[0]) ) {
            continue;
        } else if ( (charge_single2.resnm[i] == joint.resnm[1]) 
                 && (charge_single2.atmnm[i] == joint.atmnm[1]) ) {
            continue;
        } else if ( (charge_single2.resnm[i] == pr_tmp[0])
                 && (charge_single2.atmnm[i] == pr_tmp[2]) ) {
            net_charge_single2 = - charge_single2.value[i];
        } else {
            sum_charge_single2 += charge_single2.value[i];
        }
    }

    net_charge_single2 -= sum_charge_single2;
    //
    for (int i = 0; i < charge_single2.nlist; ++i) {
        if ( (charge_single2.resnm[i] == joint.resnm[0]) 
              && (charge_single2.atmnm[i] == joint.atmnm[0]) ) {
            charge_single2.value[i] = 0.0;
        } else if ( (charge_single2.resnm[i] == joint.resnm[1]) 
                 && (charge_single2.atmnm[i] == joint.atmnm[1]) ) {
            charge_single2.value[i] = 0.0;
        } else if ( (charge_single2.resnm[i] == pr_tmp[0])
                 && (charge_single2.atmnm[i] == pr_tmp[2]) ) {
            continue;
        } else {
            //charge_single2.value[i] += net_charge_single2 / (charge_single2.nlist - 2);
            charge_single2.value[i] 
            = round( (charge_single2.value[i]  + net_charge_single2 / static_cast<double>(charge_single2.nlist - 3) ) * 10000000.0 ) / 10000000.0;
        }
    }

    //for (auto c : charge_single.value) {
    //    cout << c << endl;
    //}

    //
    //  Case double cross-linking bond
    s_charge charge_double;
    auto c3_temp = ReadChargefile(chargefile[2]);
    charge_double.name = move(c3_temp.first);
    charge_double.value = move(c3_temp.second);
    charge_double.nlist = charge_double.name.size();
    for (int i = 0; i < charge_double.nlist; ++i) {
        vector<string> name_split = split(charge_double.name[i]);
        charge_double.atmnm.push_back(name_split[0]);
        charge_double.resnm.push_back(name_split[1]);
    }
    //
    double net_charge_double = 0.0;
    double sum_charge_double = 0.0;
    int cross_count = 0;
    int pol_count = 0;
    for (int i = 0; i < charge_double.nlist; ++i) {
        if (charge_double.resnm[i] == joint.resnm[0]) {
            if (charge_double.atmnm[i] == joint.atmnm[0]) {
               continue; 
            } else if (charge_double.atmnm[i] == joint.atmnm[1]) {
               continue;
            } else if (charge_double.atmnm[i] == pr_tmp[2]) {
               net_charge_double = - 2.0 * charge_double.value[i];
            } else {
               sum_charge_double += 2.0 * charge_double.value[i];
               pol_count += 1;
            }
        } else {
            cross_count += 1;
            sum_charge_double += charge_double.value[i];
        }
    }

    net_charge_double -= sum_charge_double;
    //cout << net_charge_double / (2*charge_double.nlist - cross_count - 4) << endl;
    for (int i = 0; i < charge_double.nlist; ++i) {
        if (charge_double.resnm[i] == joint.resnm[0]) {
            if (charge_double.atmnm[i] == joint.atmnm[0]) {
                charge_double.value[i] = 1.0;
            } else if (charge_double.atmnm[i] == joint.atmnm[1]) {
                charge_double.value[i] = 0.0;
            } else if (charge_double.atmnm[i] == pr_tmp[2]) {
                continue;
            } else {
                //charge_double.value[i] += net_charge_double / (2*charge_double.nlist - cross_count - 4);
                charge_double.value[i] 
                = round( (charge_double.value[i] + net_charge_double / static_cast<double>(2*charge_double.nlist - cross_count - 6)) * 10000000.0 ) / 10000000.0;
               //charge[i] += round((monocharge[i] - charge_sum / static_cast<double>(nmono - 2)) * 100000.0) / 100000.0;
            }
        } else {
           // charge_double.value[i] += net_charge_double / (2*charge_double.nlist - cross_count - 4);
           charge_double.value[i] 
           = round( (charge_double.value[i] + net_charge_double / static_cast<double>(2*charge_double.nlist - cross_count - 6)) * 10000000.0 ) / 10000000.0;
        }
    }
    
    //for (auto c : charge_double.value) {
    //    cout << c << endl;
    //}
    //
    //double s_temp = 0.0;
    //for (int i = 0; i < charge_double.nlist; ++i) {
    //    if (i < 21) {
    //        s_temp += 2 * charge_double.value[i];
    //    } else {
    //        s_temp += charge_double.value[i];
    //    }
    //}
    //cout << endl;
    //cout << "Total: " << s_temp << endl;

    vector<string> resid_list_single1, resid_list_single2, resid_list_double;
    vector<vector<string>> polid_list_single1, polid_list_single2, polid_list_double;
    // --> format (single): [[PC1, 1], [PC8. 10], .... ]
    // --> format (double): [[PC1, 1, PC2, 8], [PC8, 10, PC7, 9], .... ]
    
    for (int i = 0; i < Nc; ++i) {
        if (selcross_list[i].size() > 0) {
            int ncount = 0;
            vector<string> polid_list;
            int c = selcross_list[i].size();
            string cross_resid;
            for (int j = 0; j < c; ++j) {
                int cross_id = selcross_list[i][j];
                cross_resid = connect_itp.atoms.resid[cross_id-1];
                vector<string> bonded_resids = GetBondedResidues(connect_itp, cross_id);
                vector<string> resid_trim = TrimResid(bonded_resids);
                vector<string> polid_list_temp = GetbondResnmandResid(connect_itp, cross_id, pr_tmp[0]);
                polid_list.insert(polid_list.end(), polid_list_temp.begin(), polid_list_temp.end());
                if (count(resid_trim.begin(), resid_trim.end(), pr_tmp[0]) == 1) {
                    ncount += j+1;
                }
            }
            //
            if (ncount == 1) {
                polid_list_single1.push_back(polid_list);
                resid_list_single1.push_back(cross_resid);
            } else if (ncount == 2) {
                polid_list_single2.push_back(polid_list);
                resid_list_single2.push_back(cross_resid);
            } else if (ncount == 3) {
                polid_list_double.push_back(polid_list);
                resid_list_double.push_back(cross_resid);
            }
        }
    }
    
    cout << endl;
    cout << "Single cross linking (1): " << endl;
    //for (auto p : polid_list_single) {
    //    for (auto q : p) {
    //       cout << q << " ";
    //    }
    //    cout << endl;
    //}
    int ns1 = polid_list_single1.size();
    for (int i = 0; i < ns1; ++i) {
        string pol_resnm = polid_list_single1[i][0];
        int    pol_resid = stoi(polid_list_single1[i][1]);
        string cross_resnm = resid_list_single1[i];
        cout << pol_resnm << " " << pol_resid << " " 
             << "<== " << cross_resnm << endl;
        for (int j = 0; j < natoms; ++j) {
            string atmnm = connect_itp.atoms.atom[j];
            if ( (connect_itp.atoms.resid[j] == pol_resnm) && 
                 (connect_itp.atoms.resnr[j] == pol_resid) ){
                //string atmnm = connect_itp.atoms.atom[j];
                for (int k = 0; k < charge_single1.nlist; ++k) {
                    if ( (charge_single1.atmnm[k] == atmnm) && 
                         (charge_single1.resnm[k] == remove_digits(pol_resnm)) ) {
                        connect_itp.atoms.charge[j] = charge_single1.value[k];
                    }
                }
            } else if ( (connect_itp.atoms.resid[j] == cross_resnm) ) {
                for (int k = 0; k < charge_single1.nlist; ++k) {
                    if ( (charge_single1.atmnm[k] == atmnm) && 
                         (charge_single1.resnm[k] == remove_digits(cross_resnm)) ) {
                        connect_itp.atoms.charge[j] = charge_single1.value[k];
                    }
                }
            }
        }
    }
    
    cout << endl;
    cout << "Single cross linking (2): " << endl;
    int ns2 = polid_list_single2.size();
    for (int i = 0; i < ns2; ++i) {
        string pol_resnm = polid_list_single2[i][0];
        int    pol_resid = stoi(polid_list_single2[i][1]);
        string cross_resnm = resid_list_single2[i];
        cout << pol_resnm << " " << pol_resid << " " 
             << "<== " << cross_resnm << endl;
        for (int j = 0; j < natoms; ++j) {
            string atmnm = connect_itp.atoms.atom[j];
            if ( (connect_itp.atoms.resid[j] == pol_resnm) && 
                 (connect_itp.atoms.resnr[j] == pol_resid) ){
                //string atmnm = connect_itp.atoms.atom[j];
                for (int k = 0; k < charge_single2.nlist; ++k) {
                    if ( (charge_single2.atmnm[k] == atmnm) && 
                         (charge_single2.resnm[k] == remove_digits(pol_resnm)) ) {
                        connect_itp.atoms.charge[j] = charge_single2.value[k];
                    }
                }
            } else if ( (connect_itp.atoms.resid[j] == cross_resnm) ) {
                for (int k = 0; k < charge_single2.nlist; ++k) {
                    if ( (charge_single2.atmnm[k] == atmnm) && 
                         (charge_single2.resnm[k] == remove_digits(cross_resnm)) ) {
                        connect_itp.atoms.charge[j] = charge_single2.value[k];
                    }
                }
            }
        }
    }
    
    cout << endl;
    cout << "Double cross linking: " << endl;
    int nd = polid_list_double.size();
    for (int i = 0; i < nd; ++i) {
        string pol_resnm1 = polid_list_double[i][0];
        int    pol_resid1 = stoi(polid_list_double[i][1]);
        string pol_resnm2 = polid_list_double[i][2];
        int    pol_resid2 = stoi(polid_list_double[i][3]);
        string cross_resnm = resid_list_double[i];
        cout << pol_resnm1 << " " << pol_resid1 << " " 
             << pol_resnm2 << " " << pol_resid2 << " " 
             << "<== " <<cross_resnm << endl;
        //
        for (int j = 0; j < natoms; ++j) {
            string atmnm = connect_itp.atoms.atom[j];
            if ( (connect_itp.atoms.resid[j] == pol_resnm1) &&
                 (connect_itp.atoms.resnr[j] == pol_resid1) ) {
                for (int k = 0; k < charge_double.nlist; ++k) {
                    if ( (charge_double.atmnm[k] == atmnm) && 
                         (charge_double.resnm[k] == remove_digits(pol_resnm1)) ) {
                        connect_itp.atoms.charge[j] = charge_double.value[k];
                    }
                }
            } else if ( (connect_itp.atoms.resid[j] == pol_resnm2) &&
                         (connect_itp.atoms.resnr[j] == pol_resid2) ) {
                for (int k = 0; k < charge_double.nlist; ++k) {
                    if ( (charge_double.atmnm[k] == atmnm) && 
                         (charge_double.resnm[k] == remove_digits(pol_resnm2)) ) {
                        connect_itp.atoms.charge[j] = charge_double.value[k];
                    }
                }
            } else if ( (connect_itp.atoms.resid[j] == cross_resnm) ) {
                for (int k = 0; k < charge_double.nlist; ++k) {
                    if ( (charge_double.atmnm[k] == atmnm) && 
                         (charge_double.resnm[k] == remove_digits(cross_resnm)) ) {
                        connect_itp.atoms.charge[j] = charge_double.value[k];
                    }
                }
            }
        }
    }

    cout << endl;
    cout << "Write results to " << outitp << "..." << endl;
    Writeitp(connect_itp, outitp);
}
///////////////////////////////////////////////////////////////////////////////

void PolymerConvertItp(const string& itpfile, 
                       const string& outitp,
                       const string& resid,
                       const vector<string>& atomlist,
                       const vector<int>& joint_ref, 
                       const vector<int>& joint_mov, 
                       const int npol) {
    //
    s_itp itp;
    Readitpfile(itpfile, itp);
    //
    int nres = atomlist.size();
    int ref_index = joint_ref[1];
    int mov_index = joint_mov[1];
    vector<string> atomlist_ref, atomlist_mov, atomlist_refmov; 
    for (int i = 0; i < nres; ++i) {
       if (i != ref_index-1) {
           atomlist_ref.push_back(atomlist[i]);
       }
       if (i != mov_index-1) {
           atomlist_mov.push_back(atomlist[i]);
       }
       if ((i != ref_index-1) && (i != mov_index-1)) {
           atomlist_refmov.push_back(atomlist[i]);
       }
    }

    int resnr = 1;
    int count = 0;
    for (int i = 0; i < nres-1; ++i) {
        itp.atoms.resnr[count] = resnr;
        itp.atoms.resid[count] = resid;
        itp.atoms.atom[count] = atomlist_ref[i];
        count += 1;
    }
    resnr += 1;
    //
    for (int i = 0; i < npol-2; ++i) {
        for (int j = 0; j < nres-2; ++j) {
            itp.atoms.resnr[count] = resnr;
            itp.atoms.resid[count] = resid;
            itp.atoms.atom[count] = atomlist_refmov[j];
            count += 1;
        }
        resnr += 1;
    }
    //
    for (int i = 0; i < nres-1; ++i) {
        itp.atoms.resnr[count] = resnr;
        itp.atoms.resid[count] = resid;
        itp.atoms.atom[count] = atomlist_mov[i];
        count += 1;
    }
    cout << endl;
    cout << "Write results to " << outitp << "..." << endl;
    Writeitp(itp, outitp);
}


void PolymerConvertPdb(const string& pdbfile, 
                       const string& outpdb, 
                       const string& resid,
                       const vector<string>& atomlist,
                       const vector<int>& joint_ref,
                       const vector<int>& joint_mov,
                       const int npol) {
    PDBOPR PdbOpr;
    //
    s_pdb pdb = PdbOpr.LoadfrompdbAmber(pdbfile);
    int nres = atomlist.size();
    int ref_index = joint_ref[1];
    int mov_index = joint_mov[1];
    vector<string> atomlist_ref, atomlist_mov, atomlist_refmov; 
    for (int i = 0; i < nres; ++i) {
       if (i != ref_index-1) {
           atomlist_ref.push_back(atomlist[i]);
       }
       if (i != mov_index-1) {
           atomlist_mov.push_back(atomlist[i]);
       }
       if ((i != ref_index-1) && (i != mov_index-1)) {
           atomlist_refmov.push_back(atomlist[i]);
       }
    }

    int resnr = 1;
    int count = 0;
    for (int i = 0; i < nres-1; ++i) {
        pdb.resid[count] = resnr;
        pdb.resnm[count] = resid;
        pdb.atmnm[count] = atomlist_ref[i];
        count += 1;
    }
    resnr += 1;
    //
    for (int i = 0; i < npol-2; ++i) {
        for (int j = 0; j < nres-2; ++j) {
            pdb.resid[count] = resnr;
            pdb.resnm[count] = resid;
            pdb.atmnm[count] = atomlist_refmov[j];
            count += 1;
        }
        resnr += 1;
    }
    //
    for (int i = 0; i < nres-1; ++i) {
        pdb.resid[count] = resnr;
        pdb.resnm[count] = resid;
        pdb.atmnm[count] = atomlist_mov[i];
        count += 1;
    }
    cout << endl;
    cout << "Write results to " << outpdb << "..." << endl;
    PdbOpr.OutPdbAmber2({pdb}, outpdb);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// REMOVE_ANALYSIS ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void REMOVE_ANALYSIS(const string& itpfile, 
                     const string& outitp, 
                     const string& pdbfile,
                     const string& outpdb,
                     const vector<string>& selpol, 
                     const vector<string>& selcross, 
                     const int Nc) {
    cout << endl;
    cout << "STEP1: Read itpfile and trajfile..." << endl;
    s_itp itp;
    Readitpfile(itpfile, itp);
    int natoms = itp.atoms.resnr.size();
    //
    cout << endl;
    cout << "STEP2: Search crosslinking mater. " << endl;
    //--------------------------------------------------------------------------
    // For crosslinking material part...
    int c = selcross.size() - 1;
    vector<vector<int>> selcross_list = ExtractList(itp, Nc, natoms, selcross);
    cout << endl;
    cout << "Crosslinking material select: " << endl;
    int cross_label = 1;
    for (vector<int> sel_list : selcross_list) {
        cout << selcross[0] << cross_label << ": ";
        for (int sel : sel_list) {
            cout << sel << " ";
        }
        cout << endl;
        cross_label += 1;
    }
    //--------------------------------------------------------------------------
    vector<int> nonbond_list;
    for (int i = 0; i < Nc; ++i) {
        vector<bool> sign_list;
        for (int j = 0; j < c; ++j) {
            int cross_id =  selcross_list[i][j];
            vector<string> bonded_resids = GetBondedResidues(itp, cross_id);
            //for (auto b : bonded_resids) {
            //    cout << b << " ";
            //}
            //cout << endl;
            bool csign = has_bonded_resid(bonded_resids, selpol[0]);
            sign_list.push_back(csign);
        }
        if (any_of(sign_list.begin(), sign_list.end(), [](bool b){ return b; })) {
            cout << "Cross linking mater. "<< selcross[0] 
                 << i+1 << " has polymer-bond..." << endl;;
        } else {
            cout << "Cross linking mater. " << selcross[0] 
                 << i+1 << " does not have polymer-bond..." << endl;;
            vector<int> cell_nonbond;
            for (int j = 0; j < natoms; ++j) {
                string resid = selcross[0] + to_string(i+1);
                if (boost::iequals(resid, itp.atoms.resid[j])) {
                    nonbond_list.push_back(j+1);
                }
            }
            //nonbond_list.push_back(cell_nonbond);
        }
    }
    int Nb = nonbond_list.size();
    cout << "-->" << Nb << " atoms will be erased..." << endl;
    RemoveAtomsKeepIndex(itp, nonbond_list);
    CompactAtomIndices(itp);
    //
    Writeitp(itp, outitp);
    //
    PDBOPR PdbOpr;
    s_pdb pdb = PdbOpr.LoadfrompdbAmber(pdbfile);
    s_pdb remove_pdb = PdbOpr.Removepdb(pdb, nonbond_list);
    PdbOpr.OutPdbAmber({remove_pdb}, outpdb);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////// REMOVE_RESID_ANALYSIS ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void REMOVE_RESID_ANALYSIS(const string& itpfile, 
                           const string& outitp, 
                           const string& pdbfile, 
                           const string& outpdb, 
                           const vector<string>& resid_list) {
    cout << endl;
    cout << "STEP1: Read itpfile..." << endl;
    s_itp itp;
    Readitpfile(itpfile, itp);
    //
    cout << endl;
    cout << "STEP2: Search resid to remove..." << endl;
    vector<int> remove_list;
    int natoms = itp.atoms.n;
    //
    for (auto ri : resid_list) {
        for (int j = 0; j < natoms; ++j) {
            if (itp.atoms.resid[j] == ri) {
                remove_list.push_back(j+1);
            }
        }
    }
    //
    cout << endl;
    cout << "Select remove atoms: " << endl;
    for (auto r : remove_list) cout << r << " ";
    cout << endl;
    cout << "==> # of atoms: " << remove_list.size() << endl;
    //
    RemoveAtomsKeepIndex(itp, remove_list);
    CompactAtomIndices(itp);
    //
    Writeitp(itp, outitp);
    //
    PDBOPR PdbOpr;
    s_pdb pdb = PdbOpr.LoadfrompdbAmber(pdbfile);
    s_pdb remove_pdb = PdbOpr.Removepdb(pdb, remove_list);
    PdbOpr.OutPdbAmber({remove_pdb}, outpdb);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////// CROSS_LINKING_CHECK /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void CROSS_LINKING_CHECK(const string& itpfile, 
                         const vector<string>& selpol, 
                         const vector<string>& selcross, 
                         const int Nc) {
    cout << endl;
    cout << "STEP1: Read itpfile and trajfile..." << endl;
    s_itp itp;
    Readitpfile(itpfile, itp);
    int natoms = itp.atoms.resnr.size();
    //
    cout << endl;
    cout << "STEP2: Search crosslinking mater. " << endl;
    //--------------------------------------------------------------------------
    // For crosslinking material part...
    int c = selcross.size() - 1;
    vector<vector<int>> selcross_list = ExtractList(itp, Nc, natoms, selcross);
    cout << endl;
    cout << "Crosslinking material select: " << endl;
    int cross_label = 1;
    for (vector<int> sel_list : selcross_list) {
        cout << selcross[0] << cross_label << ": ";
        for (int sel : sel_list) {
            cout << sel << " ";
        }
        cout << endl;
        cross_label += 1;
    }
    //--------------------------------------------------------------------------
    vector<int> resid_list;
    vector<vector<int>> bonded_list;
    for (int i = 0; i < Nc; ++i) {
        if (selcross_list[i].size() > 0) {
            int ncount = 0;
            for (int j = 0; j < c; ++j) {
                int cross_id = selcross_list[i][j];
                vector<string> bonded_resids = GetBondedResidues(itp, cross_id);
                vector<string> resid_trim = TrimResid(bonded_resids);
               // for (auto b : resid_trim) {
               //     cout << b  << " ";
               // }
               // cout << endl;
                if (count(resid_trim.begin(), resid_trim.end(), selpol[0]) == 1) {
                    ncount += 1;
                    vector<int> bonded_cell = GetBondedIndex(itp, cross_id);
                    bonded_cell.push_back(cross_id);
                    bonded_list.push_back(bonded_cell);
                }
            }
            if (ncount == 2) {
                resid_list.push_back(i+1);
            }
        }
    }
    //
    for (int r : resid_list) {
        string resid = selcross[0] + to_string(r);
        cout << resid << endl;
        cout << "index: ";
        for (int i = 0; i < natoms; ++i) {
           if (boost::iequals(resid, itp.atoms.resid[i])) {
               cout << i+1 << " ";
           }
        }
        cout << endl;
    }
    //
    for (auto bonded_cell : bonded_list) {
        cout << "Bond: ";
        for (auto b : bonded_cell) {
            cout << b << " ";
        }
        cout << endl;
    }
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////// ANALYZE_BOND_AND_REMOVE /////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_BOND_AND_REMOVE(const vector<string>& pdbfile, 
                             const string& outpdb,
                             const vector<string>& itpfile, 
                             const string& outitp,
                             const string& chargefile,
                             const vector<string>& joint_string, 
                             const vector<string>& joint_string_polymer,
                             const double distance, 
                             optional<int> ncycle, 
                             optional<double> prob) {
    //
    cout << endl;
    cout << "STEP1: Read pdbfile and itpfile..." << endl;
    PDBOPR PdbOpr;
    
    cout << "=> System Section: ";
    string sys_pdbfile = pdbfile[0];
    string sys_itpfile = itpfile[0];
    s_pdb  sys_pdb = PdbOpr.LoadfrompdbAmber(sys_pdbfile);
    s_itp  sys_itp; 
    Readitpfile(sys_itpfile, sys_itp);
    int sys_natms = sys_itp.atoms.nr.size();
    cout << "OK" << endl;
    //
    cout << "=> Bond Section: ";
    string bond_pdbfile = pdbfile[1];
    string bond_itpfile = itpfile[1];
    s_pdb  bond_pdb = PdbOpr.LoadfrompdbAmber(bond_pdbfile);
    s_itp  bond_itp; 
    Readitpfile(bond_itpfile, bond_itp);
    int bond_natms = bond_itp.atoms.nr.size();
    cout << "OK" << endl;
    //
    cout << "=> imitate itp section: ";
    string imit_itpfile = itpfile[2];
    s_itp_imit imit_itp;
    Readitpimitfile(imit_itpfile, imit_itp);
    cout << "OK" << endl;

    cout << endl;
    cout << "STEP2: Bond molecules insetion" << endl;
    //
    s_itp connect_itp = sys_itp;
    //
    s_joint joint;
    for (string j : joint_string) {
        vector<string> name_split = split(j);
        joint.atmnm.push_back(name_split[0]);
        joint.resnm.push_back(name_split[1]);
    }
    //
    // get index of polymer's joint
    // cout << sys_natms << " " << bond_natms << endl;
    int bond_index = -1;
    for (int i = 0; i < bond_natms; ++i) {
        string resnm_trim = remove_digits(bond_itp.atoms.resid[i]);
        if ( (boost::iequals(bond_itp.atoms.atom[i], joint.atmnm[1]))
             && (boost::iequals(resnm_trim, joint.resnm[1])) ) {
            bond_index = i + 1;
        }
    }
    if (bond_index > 0) {
        cout << "Index (BOND): " << bond_index << endl;
    } else {
        cout << "Index (BOND) does not match your input, check your input file again..." << endl;
        return;
    }

    vector<s_pdb> pdb_list;
    s_pdb sys_pdb_new = sys_pdb;
    //
    // random number generator
    RandomGenerator rng;
    
    vector<int> sys_joint_list;
    vector<int> bond_joint_list;
    int joint_label = sys_natms + bond_index;
    int acc_count = 1;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(0.0, 1.0);
    
    for (int i = 0; i < sys_natms; ++i) {
        string resnm_trim = remove_digits(sys_itp.atoms.resid[i]);
        if ( (boost::iequals(sys_itp.atoms.atom[i], joint.atmnm[0]))
             && (boost::iequals(resnm_trim, joint.resnm[0])) ) {
            int rem_index = i + 1;
            vector<int> con_index = GetBondedIndex(sys_itp, rem_index);
            vector<int> con_index2;
            if (con_index.size() == 0) {
                continue;
            } else {
                if (dist(gen) < 1 - *prob) {
                    continue;   // skipping 50% 
                }
                //
                for (auto c1 : con_index) {
                    con_index2.push_back(c1);
                    vector<int> con_temp = GetBondedIndex(sys_itp, c1);
                    for (auto c2 : con_temp) {
                        con_index2.push_back(c2);
                    }
                }
            }
            vector<double> rem_coord = sys_pdb.coord[rem_index-1];
            //
            bool isign = false;
            s_pdb bond_pdb_new;
            for (int iter = 0; iter < *ncycle; ++iter) {
                s_pdb bond_pdb_temp = bond_pdb;
                auto vec = rng.RandomUnitVector(3);
                double theta = rng.RandomAngle();
                PerformRodrigues(bond_pdb_temp.coord, vec, theta);
                vector<double> bond_coord = bond_pdb_temp.coord[bond_index-1];
                for (int j = 0; j < bond_natms; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        bond_pdb_temp.coord[j][k] += rem_coord[k] - bond_coord[k] ;
                    }
                }
                //
                double length_min = 100.0;
                for (auto c : con_index2) {
                    for (int j = 0; j < bond_natms; ++j) {
                        double length = 0.0;
                        for (int k = 0; k < 3; ++k) {
                            length += (bond_pdb_temp.coord[j][k] - sys_pdb.coord[c][k]) 
                                       * (bond_pdb_temp.coord[j][k] - sys_pdb.coord[c][k]);
                        }
                        length = sqrt(length);
                        if (length < length_min) {
                            length_min = length;
                        }
                    }
                }
                //
                if (length_min > distance) {
                    isign = true;
                    bond_pdb_new = bond_pdb_temp;
                    break;
                }
            }
            //
            bool rsign = false;
            //
            for (int iter = 0; iter < *ncycle; ++iter) {
                double length_rem_min = 100.0;
                auto vec = rng.RandomVector(3, distance);
                //auto vec = rng.RandomUnitVector(3);
                vector<double> rem_iter_coord 
                = { rem_coord[0] + vec[0], rem_coord[1] + vec[1], rem_coord[2] + vec[2] };
                for (auto c : con_index) {
                    double length = 0.0;
                    for (int k = 0; k < 3; ++k) {
                        length += (rem_iter_coord[k] - sys_pdb.coord[c-1][k]) 
                                   * (rem_iter_coord[k] - sys_pdb.coord[c-1][k]);
                    }
                    //
                    length = sqrt(length);
                    if (length < length_rem_min) {
                        length_rem_min = length;
                    }
                }
                for (int j = 0; j < bond_natms; ++j) {
                    double length = 0.0;
                    for (int k = 0; k < 3; ++k) {
                        length += (rem_iter_coord[k] - bond_pdb_new.coord[j][k]) 
                                   * (rem_iter_coord[k] - bond_pdb_new.coord[j][k]);
                    }
                    //
                    length = sqrt(length);
                    if (length < length_rem_min) {
                        length_rem_min = length;
                    }
                }
                //
                if (length_rem_min > distance) {
                    rsign = true;
                    for (int k = 0; k < 3; ++k) {
                       sys_pdb_new.coord[rem_index-1][k] += vec[k];
                    }
                    break;
                }
            }
            //
            if (isign && rsign) {
                for (int k = 0; k < bond_natms; ++k) {
                    bond_itp.atoms.resid[k] = joint.resnm[1] + to_string(acc_count);
                }
                RemoveConnections(connect_itp, rem_index);
                connect_itp = CombineItp(connect_itp, bond_itp);
                pdb_list.push_back(bond_pdb_new);
                //
                sys_joint_list.push_back(con_index[0]);
                bond_joint_list.push_back(joint_label);
                joint_label += bond_natms;
                acc_count += 1;
            } else {
                cout << "Error: iteration step stop!" << endl;
                return;
            }
        }
    }
    //
    vector<s_pdb> pdb_list_sort;
    pdb_list_sort.push_back(sys_pdb_new);
    //
    //cout << "BOX: " << sys_pdb_new.box[0] << " " << sys_pdb_new.box[1] << " " << sys_pdb_new.box[2] << endl;
    int pcount = 1;
    for (auto pdb : pdb_list) {
        cout << "Number "  << pcount << ": " << endl;
        apply_pbc(pdb.coord, sys_pdb_new.box[0], sys_pdb_new.box[1], sys_pdb_new.box[2]);
        pdb_list_sort.push_back(pdb);
        pcount += 1;
    }
    PdbOpr.OutPdbAmber(pdb_list_sort, outpdb);
    // 
    cout << endl;
    cout << "STEP3: Update topology file..." << endl;
    UpdateTopology(sys_joint_list, bond_joint_list, connect_itp, imit_itp);
    cout << "==> main topology section is complete!!" << endl;
    //
    cout <<  endl;
    cout << "Next; Update charge..." << endl;
    s_charge charge;
    auto charge_temp = ReadChargefile(chargefile);
    charge.name = move(charge_temp.first);
    charge.value = move(charge_temp.second);
    charge.nlist = charge.name.size();
    for (int i = 0; i < charge.nlist; ++i) {
        vector<string> name_split = split(charge.name[i]);
        charge.atmnm.push_back(name_split[0]);
        charge.resnm.push_back(name_split[1]);
    }
    //
    s_joint joint_polymer;
    for (string j : joint_string_polymer) {
        vector<string> name_split = split(j);
        joint_polymer.atmnm.push_back(name_split[0]);
        joint_polymer.resnm.push_back(name_split[1]);
    }
    //
    double net_charge = 0.0;
    double sum_charge = 0.0;
    for (int i = 0; i < charge.nlist; ++i) {
        if ( (charge.resnm[i] == joint_polymer.resnm[0]) 
              && (charge.atmnm[i] == joint_polymer.atmnm[0]) ) {
            continue;
        } else if ( (charge.resnm[i] == joint_polymer.resnm[1]) 
                  && (charge.atmnm[i] == joint_polymer.atmnm[1]) ) {
            continue;
        } else if ( (charge.resnm[i] == joint.resnm[0]) 
                  && (charge.atmnm[i] == joint.atmnm[0]) ) {
            net_charge -= charge.value[i];
        } else {
            sum_charge += charge.value[i];
        }
    }
   
    net_charge -= sum_charge;
    //
    for (int i = 0; i < charge.nlist; ++i) {
        if ( (charge.resnm[i] == joint_polymer.resnm[0]) 
              && (charge.atmnm[i] == joint_polymer.atmnm[0]) ) {
            charge.value[i] = 0.0;
        } else if ( (charge.resnm[i] == joint_polymer.resnm[1]) 
                  && (charge.atmnm[i] == joint_polymer.atmnm[1]) ) {
            charge.value[i] = 0.0;
        } else if ( (charge.resnm[i] == joint.resnm[0]) 
                  && (charge.atmnm[i] == joint.atmnm[0]) ) {
            continue;
        } else {
            charge.value[i] 
            = round( (charge.value[i] + net_charge / static_cast<double>(charge.nlist - 3) ) * 10000000.0 ) / 10000000.0;
        }
    }
    //
    // Using sys_joint_list, bond_joint_list, connect_itp
    int npair = sys_joint_list.size();
    int sum_atms = connect_itp.atoms.atom.size();
    for (int i = 0; i < npair; ++i) {
        int don_id = sys_joint_list[i];
        int acc_id = bond_joint_list[i];
        //
        string don_resid       = connect_itp.atoms.resid[don_id-1];
        string don_resid_trim  = remove_digits(don_resid);
        string don_atom        = connect_itp.atoms.atom[don_id-1];
        int    don_resnr       = connect_itp.atoms.resnr[don_id - 1];
        //
        string acc_resid       = connect_itp.atoms.resid[acc_id - 1];
        string acc_resid_trim  = remove_digits(acc_resid);
        string acc_atom        = connect_itp.atoms.atom[acc_id - 1];
        //
        cout << "----------------------------------------------------" << endl;
        for (int j = 0; j < sum_atms; ++j) {
            if (connect_itp.atoms.resid[j] == acc_resid) {
               string acc_resid_temp = remove_digits(connect_itp.atoms.resid[j]);
               for (int k = 0; k < charge.nlist; ++k) {
                   if ( (acc_resid_temp == charge.resnm[k]) 
                         && (connect_itp.atoms.atom[j] == charge.atmnm[k]) ) {
                       connect_itp.atoms.charge[j] = charge.value[k];
                       cout << acc_resid << " " << connect_itp.atoms.atom[j] << "->" << charge.value[k]<< endl;
                   }
               }
            } else if ( (connect_itp.atoms.resid[j] == don_resid) 
                         && (connect_itp.atoms.resnr[j] == don_resnr) ) {
               string don_resid_temp = remove_digits(connect_itp.atoms.resid[j]);
               for (int k = 0; k <  charge.nlist; ++k) {
                   if ( (don_resid_temp == charge.resnm[k]) 
                        && (connect_itp.atoms.atom[j] == charge.atmnm[k]) ) {
                       connect_itp.atoms.charge[j] = charge.value[k];
                       cout << don_resid << " (" << connect_itp.atoms.resnr[j] << ") "  
                            << connect_itp.atoms.atom[j] << "->" << charge.value[k]<< endl;
                   }
               }
            }
        }
        cout << "----------------------------------------------------" << endl;
    }
    //
    cout << "Write results to " << outitp << "..." << endl;
    Writeitp(connect_itp, outitp);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/////////////////////// ANALYZE_BOND_AND_REMOVE_AC ////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_BOND_AND_REMOVE_AC(const vector<string>& pdbfile,
                                const string& outpdb,
                                const vector<string>& itpfile,
                                const string& outitp,
                                const string& chargefile,
                                const vector<string>& joint_string, 
                                const vector<string>& joint_string_polymer,
                                const vector<string>& joint_string_cross,
                                const vector<string>& joint_string_remove,
                                const int nbond, 
                                const double distance, 
                                optional<int> ncycle) {
    //
    cout << endl;
    cout << "STEP1: Read pdbfile and itpfile..." << endl;
    PDBOPR PdbOpr;
    
    cout << "=> System Section: ";
    string sys_pdbfile = pdbfile[0];
    string sys_itpfile = itpfile[0];
    s_pdb  sys_pdb = PdbOpr.LoadfrompdbAmber(sys_pdbfile);
    s_itp  sys_itp; 
    Readitpfile(sys_itpfile, sys_itp);
    int sys_natms = sys_itp.atoms.nr.size();
    cout << "OK" << endl;
    //
    cout << "=> Bond Section: ";
    string bond_pdbfile = pdbfile[1];
    string bond_itpfile = itpfile[1];
    s_pdb  bond_pdb = PdbOpr.LoadfrompdbAmber(bond_pdbfile);
    s_itp  bond_itp; 
    Readitpfile(bond_itpfile, bond_itp);
    int bond_natms = bond_itp.atoms.nr.size();
    cout << "OK" << endl;
    //
    cout << "=> imitate itp section: ";
    string imit_itpfile = itpfile[2];
    s_itp_imit imit_itp;
    Readitpimitfile(imit_itpfile, imit_itp);
    cout << "OK" << endl;

    cout << endl;
    cout << "STEP2: Bond molecules insetion" << endl;
    //
    s_itp connect_itp = sys_itp;
    //
    s_joint joint;
    for (string j : joint_string) {
        vector<string> name_split = split(j);
        joint.atmnm.push_back(name_split[0]);
        joint.resnm.push_back(name_split[1]);
    }
    //
    // get index of polymer's joint
    // cout << sys_natms << " " << bond_natms << endl;
    int bond_index = -1;
    for (int i = 0; i < bond_natms; ++i) {
        string resnm_trim = remove_digits(bond_itp.atoms.resid[i]);
        if ( (boost::iequals(bond_itp.atoms.atom[i], joint.atmnm[1]))
             && (boost::iequals(resnm_trim, joint.resnm[1])) ) {
            bond_index = i + 1;
        }
    }
    if (bond_index > 0) {
        cout << "Index (BOND): " << bond_index << endl;
    } else {
        cout << "Index (BOND) does not match your input, check your input file again..." << endl;
        return;
    }
    //
    int rem_index = -1;
    for (int i = 0; i < bond_natms; ++i) {
        string resnm_trim = remove_digits(bond_itp.atoms.resid[i]);
        if ( (boost::iequals(bond_itp.atoms.atom[i], joint.atmnm[2]))
             && (boost::iequals(resnm_trim, joint.resnm[2])) ) {
            rem_index = i + 1;
        }
    }
    if (rem_index > 0) {
        cout << "Index (REMOVE): " << rem_index << endl;
    } else {
        cout << "Index (REMOVE) does not match your input, check your input file again..." << endl;
        return;
    }

    vector<s_pdb> pdb_list;
    pdb_list.push_back(sys_pdb);
    //
    // random number generator
    RandomGenerator rng;
    RandomGenerator rng2;
    
    vector<int> sys_joint_list;
    vector<int> bond_joint_list;
    int joint_label = sys_natms + bond_index;
    int acc_count = 1;
    //
    RemoveConnections(bond_itp, rem_index);
    //
    for (int i = 0; i <  sys_natms; ++i) {
        string resnm_trim = remove_digits(sys_itp.atoms.resid[i]);
        if ( (boost::iequals(sys_itp.atoms.atom[i], joint.atmnm[0]))
             && (boost::iequals(resnm_trim, joint.resnm[0])) ) {
            int con_id = i + 1;
            vector<double> sys_coord = sys_pdb.coord[con_id-1];
            vector<int> con_index = GetBondedIndex(sys_itp, con_id);
            vector<int> con_index2;
            if (con_index.size() > static_cast<size_t>(nbond)) {
                continue;
            } else {
                for (auto c1 : con_index) {
                    con_index2.push_back(c1);
                    vector<int> con_temp = GetBondedIndex(sys_itp, c1);
                    for (auto c2 : con_temp) {
                        con_index2.push_back(c2);
                    }
                }
            }
            //
            bool isign = false;
            s_pdb bond_pdb_new;
            for (int iter = 0; iter < *ncycle; ++iter) {
                s_pdb bond_pdb_temp = bond_pdb;
                auto vec_norm = rng.RandomUnitVector(3);
                double theta  = rng.RandomAngle();
                auto vec      = rng2.RandomVector(3, distance);
                PerformRodrigues(bond_pdb_temp.coord, vec_norm, theta);
                vector<double> rem_coord = bond_pdb_temp.coord[rem_index-1];
                vector<double> bond_coord = bond_pdb_temp.coord[bond_index-1];
                for (int j = 0; j < bond_natms; ++j) {
                     for (int k = 0; k < 3; ++k) {
                         bond_pdb_temp.coord[j][k] += sys_coord[k] - rem_coord[k];
                     }
                }
                //
                for (int k = 0; k < 3; ++k) {
                    bond_pdb_temp.coord[rem_index-1][k] += vec[k] + bond_coord[k] - rem_coord[k];
                }
                //
                double length_min = 100.0;
                for (auto c : con_index2) {
                    for (int j = 0; j < bond_natms; ++j) {
                        double length = 0.0;
                        for (int k = 0; k < 3; ++k) {
                            length += (bond_pdb_temp.coord[j][k] - sys_pdb.coord[c][k]) 
                                       * (bond_pdb_temp.coord[j][k] - sys_pdb.coord[c][k]);
                        }
                        length = sqrt(length);
                        if (length < length_min) {
                            length_min = length;
                        }
                    }
                }
                //
                if (length_min > distance) {
                    isign = true;
                    bond_pdb_new = bond_pdb_temp;
                    break;
                }
            }
            //
            if (isign) {
                for (int j = 0; j < bond_natms; ++j) {
                    bond_itp.atoms.resid[j] = joint.resnm[1] + to_string(acc_count);
                }
                connect_itp = CombineItp(connect_itp, bond_itp);
                sys_joint_list.push_back(con_id);
                bond_joint_list.push_back(joint_label);
                pdb_list.push_back(bond_pdb_new);
                joint_label += bond_natms;
                acc_count += 1;
            } else {
                cout << "Error: iteration step stop!" << endl;
                return;
            }
        }
    }
    PdbOpr.OutPdbAmber(pdb_list, outpdb);
    // 
    cout << endl;
    cout << "STEP3: Update topology file..." << endl;
    UpdateTopology(sys_joint_list, bond_joint_list, connect_itp, imit_itp);
    cout << "==> main topology section is complete!!" << endl;
    //
    cout <<  endl;
    cout << "Next: Update charge..." << endl;
    s_charge charge;
    auto charge_temp = ReadChargefile(chargefile);
    charge.name = move(charge_temp.first);
    charge.value = move(charge_temp.second);
    charge.nlist = charge.name.size();
    for (int i = 0; i < charge.nlist; ++i) {
        vector<string> name_split = split(charge.name[i]);
        charge.atmnm.push_back(name_split[0]);
        charge.resnm.push_back(name_split[1]);
    }
    
    //--------------------------------------------------//
    s_joint joint_polymer = MakeJoint(joint_string_polymer);
    s_joint joint_cross   = MakeJoint(joint_string_cross);
    s_joint joint_remove  = MakeJoint(joint_string_remove);
    //--------------------------------------------------//

    double net_charge = 0.0;
    double sum_charge = 0.0;
    for (int i = 0; i < charge.nlist; ++i) {
        if ( (charge.resnm[i] == joint_polymer.resnm[0]) 
              && (charge.atmnm[i] == joint_polymer.atmnm[0]) ) {
            continue;
        } else if ( (charge.resnm[i] == joint_polymer.resnm[1]) 
                  && (charge.atmnm[i] == joint_polymer.atmnm[1]) ) {
            continue;
        } else if ( (charge.resnm[i] == joint_remove.resnm[0]) 
                  && (charge.atmnm[i] == joint_remove.atmnm[0]) ) {
            net_charge -= charge.value[i];
        } else if ( (charge.resnm[i] == joint_remove.resnm[1]) 
                  && (charge.atmnm[i] == joint_remove.atmnm[1]) ) {
            net_charge -= charge.value[i];
        } else {
            sum_charge += charge.value[i];
        }
    }

    net_charge -= sum_charge;
    //
    for (int i = 0; i < charge.nlist; ++i) {
        if ( (charge.resnm[i] == joint_polymer.resnm[0]) 
              && (charge.atmnm[i] == joint_polymer.atmnm[0]) ) {
            charge.value[i] = 0.0;
        } else if ( (charge.resnm[i] == joint_polymer.resnm[1]) 
                  && (charge.atmnm[i] == joint_polymer.atmnm[1]) ) {
            charge.value[i] = 0.0;
        } else if ( (charge.resnm[i] == joint_remove.resnm[0]) 
                  && (charge.atmnm[i] == joint_remove.atmnm[0]) ) {
            continue;
        } else if ( (charge.resnm[i] == joint_remove.resnm[1]) 
                  && (charge.atmnm[i] == joint_remove.atmnm[1]) ) {
            continue;
        } else {
            charge.value[i] 
            = round( (charge.value[i] + net_charge / static_cast<double>(charge.nlist - 4) ) * 10000000.0 ) / 10000000.0;
        }
    }
   
    int npair = bond_joint_list.size();
    int sum_atms = connect_itp.atoms.atom.size();
    for (int i = 0; i < npair; ++i) {
        int don_id = sys_joint_list[i];
        int acc_id = bond_joint_list[i];
        //
        string don_resid       = connect_itp.atoms.resid[don_id-1];
        string don_resid_trim  = remove_digits(don_resid);
        string don_atom        = connect_itp.atoms.atom[don_id-1];
        //
        string acc_resid       = connect_itp.atoms.resid[acc_id - 1];
        string acc_resid_trim  = remove_digits(acc_resid);
        string acc_atom        = connect_itp.atoms.atom[acc_id - 1];
        //
        // Search cross_id
        int cross_id = -1;
        for (int j = 0; j < sum_atms; ++j) {
           if ( (connect_itp.atoms.resid[j] == don_resid) 
                && (connect_itp.atoms.atom[j] == joint_cross.atmnm[0]) ) {
               cross_id = j + 1;
               break;
           }
        }
        if (cross_id < 0) {
            cout << "Error: Cannot search cross_id" << endl;
            return;
        }
        //
        int pol_id = -1;
        vector<int> bonded_index = GetBondedIndex(connect_itp, cross_id);
        for (auto b : bonded_index) {
            string bonded_resid_trim = remove_digits(connect_itp.atoms.resid[b-1]);
            if ( ( bonded_resid_trim == joint_cross.resnm[1] ) 
                 && ( connect_itp.atoms.atom[b-1] == joint_cross.atmnm[1] ) ) {
                pol_id = b;
                break;
            }
        }
        if (pol_id < 0) {
            cout << "CROSS: " << cross_id << endl;
            cout << "Error: Cannot search pol_id" << endl;
            return;
        }
        cout << "----------------------------------------------------" << endl;
        cout << cross_id << " and " << pol_id << " is connected..." << endl;
        //
        string cross_resid       = connect_itp.atoms.resid[cross_id-1];
        string cross_resid_trim  = remove_digits(cross_resid);
        //
        string pol_resid       = connect_itp.atoms.resid[pol_id-1];
        string pol_resid_trim  = remove_digits(pol_resid);
        int pol_resnr          = connect_itp.atoms.resnr[pol_id - 1];
        //
        for (int j = 0; j < sum_atms; ++j) {
            //
            // For cross
            if (connect_itp.atoms.resid[j] == cross_resid) {
                for (int k = 0; k < charge.nlist; ++k) {
                    if ( (cross_resid_trim == charge.resnm[k]) 
                        && (connect_itp.atoms.atom[j] == charge.atmnm[k]) ) {
                        connect_itp.atoms.charge[j] = charge.value[k];
                        cout << cross_resid << " " << connect_itp.atoms.atom[j] << " -> " << charge.value[k] << endl;
                    }
                }
            }
            //
            // For acc
            if (connect_itp.atoms.resid[j] == acc_resid) {
                for (int k = 0; k < charge.nlist; ++k) {
                    if ( (acc_resid_trim == charge.resnm[k]) 
                        && (connect_itp.atoms.atom[j] == charge.atmnm[k]) ) {
                        connect_itp.atoms.charge[j] = charge.value[k];
                        cout << acc_resid << " " << connect_itp.atoms.atom[j] << " -> " << charge.value[k] << endl;
                    }
                }
            }
            //
            // For pol
            if ( (connect_itp.atoms.resid[j] == pol_resid) 
                 && (connect_itp.atoms.resnr[j] == pol_resnr) ) {
                for (int k = 0; k < charge.nlist; ++k) {
                    if ( (pol_resid_trim == charge.resnm[k]) 
                        && (connect_itp.atoms.atom[j] == charge.atmnm[k]) ) {
                        connect_itp.atoms.charge[j] = charge.value[k];
                        cout << pol_resid << "(" << pol_resnr << ") " 
                             << connect_itp.atoms.atom[j] << " -> " << charge.value[k] << endl;
                    }
                }
            }
        }
        //
        cout << "----------------------------------------------------" << endl;
    }

    cout << "Write results to " << outitp << "..." << endl;
    Writeitp(connect_itp, outitp);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////// MAKE_CROSSLINKING_DRY ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void MAKE_CROSSLINKING_DRY(const vector<string>& itpfile, 
                           const string& trajfile,
                           const string& outitp,
                           const string& pdbfile,
                           const string& outpdb,
                           const vector<string>& chargefile,
                           const vector<string>& selpol,
                           const vector<string>& selcross,
                           const vector<string>& drypol,
                           const vector<string>& drycross,
                           const vector<string>& joint_string,
                           const int Np,
                           const int Nc,
                           const double rc) {
    cout << endl;
    cout << "STEP1: Read itpfile and trajfile..." << endl;
    s_itp itp;
    Readitpfile(itpfile[0], itp);
    //
    cout << endl; 
    cout << "Use last frame..." << endl;
    vector<vector<double>> lastcoord; 
    vector<double> lastbox;
    if (boost::iequals(get_extension(trajfile), "xtc")) {
        //
        s_xtc xtc = ReadXTC(trajfile);
        int nlast = xtc.nframes - 1;
        lastcoord = xtc.coord[nlast];
        lastbox.push_back(xtc.box[nlast][0]);
        lastbox.push_back(xtc.box[nlast][4]);
        lastbox.push_back(xtc.box[nlast][8]);
    } else {
        cout << "Error: your traj. format is not supported..." << endl;
        return;
    }
    //
    cout << endl;
    cout << "STEP2: Creating bonds of crosslinking mater. and polymers..." << endl;
    
    //--------------------------------------------------------------------------
    // For polymer part...
    int natoms = itp.atoms.resid.size();
    vector<vector<int>> selpol_list = ExtractListSingle(itp, Np, natoms, selpol);
    //
    cout << endl;
    cout << "Polymer select: " << endl;
    int pol_label = 1;
    for (vector<int> sel_list : selpol_list) {
        cout << selpol[0] << pol_label << ": ";
        for (int sel : sel_list) {
            cout << sel << " ";
        }
        cout << endl;
        pol_label += 1;
    }
    int npol = selpol_list[0].size();
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // For crosslinking material part...
    vector<vector<int>> selcross_list = ExtractList(itp, Nc, natoms, selcross);
    cout << endl;
    cout << "Crosslinking material select: " << endl;
    int cross_label = 1;
    for (vector<int> sel_list : selcross_list) {
        cout << selcross[0] << cross_label << ": ";
        for (int sel : sel_list) {
            cout << sel << " ";
        }
        cout << endl;
        cross_label += 1;
    }
    //--------------------------------------------------------------------------

    vector<vector<vector<int>>> connect_list;
    for (int i = 0; i < Np; ++i) {
        vector<vector<int>> connect_pol;
        for (int j = 0; j < npol; ++j) {
            vector<int> connect_cell;
            int joint_pol = selpol_list[i][j] - 1;
            //
            vector<string> bonded_resids = GetBondedResidues(itp, selpol_list[i][j]);
            vector<string> resid_trim = TrimResid(bonded_resids);
            if (find(resid_trim.begin(), resid_trim.end(), selcross[0]) != resid_trim.end()) {
                //cout << selpol_list[i][j] << " is already bonded..." << endl;
                connect_cell = {-1};
                connect_pol.push_back(connect_cell);
            } else {
                //
                vector<double> joint_crd = lastcoord[joint_pol];
                //
                for (int k = 0; k < Nc; ++k) {
                    int c1 = selcross_list[k][0] - 1;
                    int c2 = selcross_list[k][1] - 1;
                    //
                    vector<double> c1_crd = lastcoord[c1];
                    vector<double> c2_crd = lastcoord[c2];
                    //
                    double length_jc1 = dist2_pbc(c1_crd, joint_crd, lastbox);
                    length_jc1 = sqrt(length_jc1);
                    double length_jc2 = dist2_pbc(c2_crd, joint_crd, lastbox);
                    length_jc2 = sqrt(length_jc2);
                    if (length_jc1 < rc) {
                        connect_cell.push_back(c1+1);
                    }
                    //
                    if (length_jc2 < rc) {
                        connect_cell.push_back(c2+1);
                    }
                }
                // 
                connect_pol.push_back(connect_cell);
            }
        }
        connect_list.push_back(connect_pol);
    }

    //
    vector<int> don_joint_list;
    vector<int> acc_joint_list;
    //
    cout << endl;
    for (int i = 0; i < Np; ++i) {
        for (int j = 0; j < npol; ++j) {
           int joint_pol = selpol_list[i][j] - 1;
           int ns = connect_list[i][j].size();
           //cout << i << " " << j << " " << ns << endl; 
           if (ns > 1) {
               cout << selpol[0] << i+1 << " (resid " << j+1 << ") is connected ";
               for (int k = 0; k < ns; ++k) {
                   cout << connect_list[i][j][k] << " ";
               }
               double length_min = 100.0;
               int idx_min = 0;
               for (int k = 0; k < ns; ++k) {
                  int idx = connect_list[i][j][k];
                  double length = dist2_pbc(lastcoord[idx-1], lastcoord[joint_pol], lastbox);
                  length = sqrt(length);
                  if (length < length_min) {
                      length_min = length;
                      idx_min = idx;
                  }
               }
               connect_list[i][j] = filter_to_one(connect_list[i][j], idx_min); 
               cout << endl;
               cout << "==> Pick up 1-component ";
               cout << connect_list[i][j][0] << endl;
           }
        }
    }
    //
    // Explain vector size: connect_list[Np][npol][0]
    cout << endl;
    cout << "Connection pair-->" << endl;
    vector<int> don_joint_list_temp;
    vector<int> acc_joint_list_temp;
    for (int i = 0; i < Np; ++i) {
        for (int j = 0; j < npol; ++j) {
            int pol_id = selpol_list[i][j];
            if (connect_list[i][j].size() > 0) {
                int cross_id = connect_list[i][j][0];
                if (cross_id > 0) {
                    cout << selpol[0] << i+1 << " (" << j+1 << ")" << " == " 
                         << pol_id << ": " << cross_id  << endl;
                    don_joint_list_temp.push_back(pol_id);
                    acc_joint_list_temp.push_back(cross_id);
                }
            }
        }
    }
    cout << endl;
    //
    vector<int> a_list = find_duplicates(acc_joint_list_temp);
    int acc_size = acc_joint_list_temp.size();
    for (int a : a_list) {
        vector<int> dup_list;
        for (int i = 0; i < acc_size; ++i) {
            if (acc_joint_list_temp[i] == a) {
                dup_list.push_back(i+1);
            }
        }
        //
        double length_min = 100.0;
        int    idx_min = 0;
        for (int d : dup_list) {
            int acc_id = acc_joint_list_temp[d-1];
            int don_id = don_joint_list_temp[d-1];
            double length = dist2_pbc(lastcoord[acc_id-1], lastcoord[don_id-1], lastbox);
            if (length < length_min) {
                length_min = length;
                idx_min = d;
            }
        }
        vector<int> dup_list_except = except_id(dup_list, idx_min);
        erase_by_1based_indices(acc_joint_list_temp, dup_list_except);
        erase_by_1based_indices(don_joint_list_temp, dup_list_except);
    }
    //
    int temp_index = acc_joint_list_temp.size();
    for (int idx = 0; idx < temp_index; ++idx) {
        vector<string> bonded_resids = GetBondedResidues(itp, acc_joint_list_temp[idx]);
        vector<string> resid_trim = TrimResid(bonded_resids);
        if (find(resid_trim.begin(), resid_trim.end(), selpol[0]) != resid_trim.end()) {
            //cout << selcross[0] << " and " << selpol[0] << " already bonded: " << endl;
            cout << don_joint_list_temp[idx] << " and " << acc_joint_list_temp[idx] << " is already bonded..." << endl;
            continue;
        }
        bool b_sign = is_bonded(itp, acc_joint_list_temp[idx], don_joint_list_temp[idx]);
        if (b_sign) {
            cout << don_joint_list_temp[idx] << " and " << acc_joint_list_temp[idx] << " is already bonded..." << endl;
        } else {
            cout << acc_joint_list_temp[idx] << " and " << don_joint_list_temp[idx] << " does not bonded, push back..." << endl;
            acc_joint_list.push_back(acc_joint_list_temp[idx]);
            don_joint_list.push_back(don_joint_list_temp[idx]);
        }
    }
   
    if (acc_joint_list.size() != don_joint_list.size()) {
        cout << "Error: acc_joint_list.size() != don_joint_list.size()" << endl;
        return;
    }

    //// Debug
    //cout << endl;
    //vector<int> a_index = find_duplicates(acc_joint_list);
    //for (auto a : a_index) {
    //    cout << "!!!" << a << endl;
    //}
    
    //
    cout << endl;
    cout << "STEP3: Update topology..." << endl;
    s_itp connect_itp = itp; 
    //
    cout << endl;
    cout << "Read itpimitfile..." << endl;
    s_itp_imit imit_itp;
    Readitpimitfile(itpfile[1], imit_itp);
    //
    UpdateTopology(don_joint_list, acc_joint_list, connect_itp, imit_itp);
    if (imit_itp.change_sign) {
        cout << endl;
        cout << "CHANGE TPOLOGY SECTION: " << endl;
        // Bond section
        ChangeBond(connect_itp, imit_itp, don_joint_list, acc_joint_list);
        // Angle section
        ChangeAngle(connect_itp, imit_itp, don_joint_list, acc_joint_list);
        // Dihedral section
        ChangeDihedral(connect_itp, imit_itp, don_joint_list, acc_joint_list);
    }
    //
    cout << "--------------------------------------------" << endl;
    cout << endl;
    cout << "REMOVE SECTION: " << endl;
    vector<int> remove_index; // remove index of acc. and don. side
    //
    // don. side
    cout << "==> don. side" << endl;
    int nrem_p = drypol.size();
    cout << "# of remove atoms in acc. side: " << nrem_p << endl;
    cout << endl;
    //
    for (auto don_id : don_joint_list) {
        cout << "REMOVE INFO: " << connect_itp.atoms.resid[don_id-1] << " " 
             << connect_itp.atoms.resnr[don_id-1] << endl;
        for (int r = 0; r < nrem_p; ++r) { 
          vector<string> name_split = split(drypol[r]);
          string trim_resid = remove_digits(connect_itp.atoms.resid[don_id-1]);
          if ( boost::iequals(trim_resid, name_split[1]) ) {
              for (int n = 0; n < natoms; ++n) {
                  if ( boost::iequals(connect_itp.atoms.resid[don_id-1], connect_itp.atoms.resid[n])
                       && connect_itp.atoms.resnr[don_id-1] == connect_itp.atoms.resnr[n]
                       && boost::iequals(name_split[0], connect_itp.atoms.atom[n]) ) {
                      cout << "--> " << n+1 << " index removed..." << endl;
                      remove_index.push_back(n+1);
                      break;
                  } else if (n == natoms-1) {
                      cout << "--> " << "not found remove index" << endl;
                  }
              }
          }
        }
    }
    // acc. side
    cout << endl;
    cout << "==> acc. side" << endl;
    int nrem = drycross.size();
    cout << "# of remove atoms in acc. side: " << nrem << endl;
    cout << endl;
    
    for (int acc_id : acc_joint_list) {
        //
        cout << "REMOVE INFO: " << connect_itp.atoms.resid[acc_id-1] << " ";
        string trim_resid = remove_digits(connect_itp.atoms.resid[acc_id-1]);
        //
        if ( boost::iequals(trim_resid, selcross[0])
             && boost::iequals(connect_itp.atoms.atom[acc_id-1], selcross[1]) ) {
            cout << connect_itp.atoms.atom[acc_id-1] << endl;
            for (int r = 0; r < nrem/2; ++r) { // double side cross-link
                vector<string> name_split = split(drycross[r]);
                if ( boost::iequals(trim_resid, name_split[1]) ) {
                    for (int n = 0; n < natoms; ++n) {
                        if ( boost::iequals(connect_itp.atoms.resid[acc_id-1], connect_itp.atoms.resid[n])
                             && connect_itp.atoms.resnr[acc_id-1] == connect_itp.atoms.resnr[n]
                             && boost::iequals(name_split[0], connect_itp.atoms.atom[n]) ) {
                            cout << "--> " << n+1 << " index removed..." << endl;
                            remove_index.push_back(n+1);
                            break;
                        } else if (n == natoms-1) {
                            cout << "--> " << "not found remove index" << endl;
                        }
                    }
                }
            }
        } else if ( boost::iequals(trim_resid, selcross[0])
                    && boost::iequals(connect_itp.atoms.atom[acc_id-1], selcross[2]) ) {
            cout << connect_itp.atoms.atom[acc_id-1] << endl;
            for (int r = 2; r < nrem; ++r) {
                vector<string> name_split = split(drycross[r]);
                if ( boost::iequals(trim_resid, name_split[1]) ) {
                    for (int n = 0; n < natoms; ++n) {
                        if ( boost::iequals(connect_itp.atoms.resid[acc_id-1], connect_itp.atoms.resid[n])
                             && connect_itp.atoms.resnr[acc_id-1] == connect_itp.atoms.resnr[n]
                             && boost::iequals(name_split[0], connect_itp.atoms.atom[n]) ) {
                            cout << "--> " << n+1 << " index removed..." << endl;
                            remove_index.push_back(n+1);
                            break;
                        } else if (n == natoms-1) {
                            cout << "--> " << "not found remove index" << endl;
                        }
                    }
                }
            }
        }
    }

    cout << endl;
    cout << "# of total remove atoms: " << remove_index.size() << endl;
    
    // Remove index
    RemoveAtomsKeepIndex(connect_itp, remove_index);
    // Compact
    CompactAtomIndices(connect_itp);
    natoms = connect_itp.atoms.atom.size();
    cout << "--------------------------------------------" << endl;
    
    cout << endl;
    cout << "--------------------------------------------" << endl;
    cout << "CHARGE SECTION: " << endl;
    vector<string> resid_list_single1, resid_list_single2, resid_list_double;
    vector<vector<string>> polid_list_single1, polid_list_single2, polid_list_double;
    // --> format (single): [[PC1, 1], [PC8. 10], .... ]
    // --> format (double): [[PC1, 1, PC2, 8], [PC8, 10, PC7, 9], .... ]
  
    // ---->
    s_joint joint = MakeJoint(joint_string); 
    
    // Case 1-side
    s_charge charge_single1 = TrimJointSingleSideCharge(chargefile[0], joint); 
    
    // Case 2-side
    s_charge charge_single2 = TrimJointSingleSideCharge(chargefile[1], joint);
    
    // Case double cross-linking bond
    s_charge charge_double = TrimJointDoubleSideCharge(chargefile[2], joint);
    // <-----

    vector<vector<int>> selcross_list_compact = ExtractList(connect_itp, Nc, natoms, selcross);

    for (int i = 0; i < Nc; ++i) {
        int ncount = 0;
        vector<string> polid_list;
        int c = selcross_list_compact[i].size();
        string cross_resid;
        for (int j = 0; j < c; ++j) {
            int cross_id = selcross_list_compact[i][j];
            cross_resid = connect_itp.atoms.resid[cross_id-1];
            vector<string> bonded_resids = GetBondedResidues(connect_itp, cross_id);
            vector<string> resid_trim = TrimResid(bonded_resids);
            vector<string> polid_list_temp = GetbondResnmandResid(connect_itp, cross_id, selpol[0]);
            polid_list.insert(polid_list.end(), polid_list_temp.begin(), polid_list_temp.end());
            //
            if (count(resid_trim.begin(), resid_trim.end(), selpol[0]) == 1) {
                ncount += j+1;
            }
        }
        //
        if (ncount == 1) {
            polid_list_single1.push_back(polid_list);
            resid_list_single1.push_back(cross_resid);
        } else if (ncount == 2) {
            polid_list_single2.push_back(polid_list);
            resid_list_single2.push_back(cross_resid);
        } else if (ncount == 3) {
            polid_list_double.push_back(polid_list);
            resid_list_double.push_back(cross_resid);
        }
    }

    cout << endl;
    cout << "Single cross linking (1): " << endl;
    UpdateSingleSideCharge(connect_itp, charge_single1, polid_list_single1, resid_list_single1);
    
    cout << endl;
    cout << "Single cross linking (2): " << endl;
    UpdateSingleSideCharge(connect_itp, charge_single2, polid_list_single2, resid_list_single2);

    cout << endl;
    cout << "Double cross linking: " << endl;
    UpdateDoubleSideCharge(connect_itp, charge_double, polid_list_double, resid_list_double);

    ////
    //// Debug
    //double ics = 0.0;
    //for (auto c : connect_itp.atoms.charge) {
    //    ics += c;
    //}
    //cout << "Toatal system charge: " << ics << endl;
    cout << "--------------------------------------------" << endl;

    cout << endl;
    cout << "Write results to " << outitp << "..." << endl;
    Writeitp(connect_itp, outitp);
    
    PDBOPR PdbOpr;
    s_pdb pdb = PdbOpr.LoadfrompdbAmber(pdbfile);
    for (int n = 0; n <  natoms; ++n) {
        for (int k = 0; k < 3; ++k) {
            pdb.coord[n][k] = lastcoord[n][k];
        }
    }
    s_pdb remove_pdb = PdbOpr.Removepdb(pdb, remove_index);
    remove_pdb.box = {lastbox[0], lastbox[1], lastbox[2], 90.0, 90.0, 90.0};
    PdbOpr.OutPdbAmber({remove_pdb}, outpdb);
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/////////////////////// MAKE_CROSSLINKING_RADICAL /////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void MAKE_CROSSLINKING_RADICAL(const vector<string>& itpfile, 
                               const string& trajfile,
                               const string& outitp,
                               const vector<string>& chargefile,
                               const vector<string>& selpol,
                               const vector<string>& selcross,
                               const int Np,
                               const int Nc,
                               const double rc) {
    cout << endl;
    cout << "STEP1: Read itpfile and trajfile..." << endl;
    s_itp itp;
    Readitpfile(itpfile[0], itp);

    s_itp_imit imit_itp;
    Readitpimitfile(itpfile[1], imit_itp);
    
    // <----
    cout << endl; 
    cout << "Use last frame..." << endl;
    vector<vector<double>> lastcoord; 
    vector<double> lastbox;
    if (boost::iequals(get_extension(trajfile), "xtc")) {
        //
        s_xtc xtc = ReadXTC(trajfile);
        int nlast = xtc.nframes - 1;
        lastcoord = xtc.coord[nlast];
        lastbox.push_back(xtc.box[nlast][0]);
        lastbox.push_back(xtc.box[nlast][4]);
        lastbox.push_back(xtc.box[nlast][8]);
    } else {
        cout << "Error: your traj. format is not supported..." << endl;
        return;
    }
    //
    cout << endl;
    cout << "STEP2: Creating bonds of crosslinking mater. and polymers..." << endl;
    
    //--------------------------------------------------------------------------
    // For polymer part...
    int natoms = itp.atoms.resid.size();
    if (selpol.size() != 2) {
      cerr << "Error: selpol.size() != 2, check your input file, again..." << endl;
      return;
    }
    vector<string> selpol_split1 = split_reverse(selpol[0]);
    vector<string> selpol_split2 = split_reverse(selpol[1]);
    vector<vector<int>> selpol1_list = ExtractListSingle(itp, Np, natoms, selpol_split1);
    vector<vector<int>> selpol2_list = ExtractListSingle(itp, Np, natoms, selpol_split2);
    vector<int> selpol1_single, selpol2_single;
    //
    if (selpol_split1[0] != selpol_split2[0]) {
        cout << "Error: pair of resid of selpol is different..." << endl;
        cout << selpol_split1[0] << "!=" << selpol_split2[0] << endl;
        return;
    }
    cout << endl;
    cout << "Polymer select: " << endl;
    //
    int selpol_size = selpol1_list.size();
    for (int i = 0; i < selpol_size; ++i) {
        int sel_size = selpol2_list[i].size();
        cout << selpol_split1[0] << i + 1 << ": ";
        if (sel_size != 1) {
            throw runtime_error("Your topology is not supported...");
        }
        for (int j = 0; j < sel_size; ++j) {
            cout << selpol1_list[i][j] << " "
                 << selpol2_list[i][j] << endl;
            selpol1_single.push_back(selpol1_list[i][j]);
            selpol2_single.push_back(selpol2_list[i][j]);
        }
    }
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // For crosslinking material part...
    if (selcross.size() != 4) {
      cerr << "Error: selcross.size() != 4, check your input file, again..." << endl;
      return;
    }

    vector<string> selcross_split1 = split_reverse(selcross[0]);
    vector<string> selcross_split2 = split_reverse(selcross[1]);
    vector<string> selcross_split3 = split_reverse(selcross[2]);
    vector<string> selcross_split4 = split_reverse(selcross[3]);
    //
    vector<vector<int>> selcross1_list = ExtractListSingle(itp, Nc, natoms, selcross_split1);
    vector<vector<int>> selcross2_list = ExtractListSingle(itp, Nc, natoms, selcross_split2);
    vector<vector<int>> selcross3_list = ExtractListSingle(itp, Nc, natoms, selcross_split3);
    vector<vector<int>> selcross4_list = ExtractListSingle(itp, Nc, natoms, selcross_split4);
    vector<int> selcross1_single, selcross2_single, selcross3_single, selcross4_single;
    //
    if ( (selcross_split1[0] != selcross_split2[0])
         || (selcross_split2[0] != selcross_split3[0])
         || (selcross_split3[0] != selcross_split4[0]) ) {
        cout << "Error: pair of resid of selcross is different..." << endl;
        cout << selcross_split1[0] << "!=" << selcross_split2[0] << endl;
        cout << selcross_split2[0] << "!=" << selcross_split3[0] << endl;
        cout << selcross_split3[0] << "!=" << selcross_split4[0] << endl;
        return;
    
    }
    cout << endl;
    cout << "Crosslinking material select: " << endl;
    int selcross_size = selcross1_list.size();
    for (int i = 0; i < selcross_size; ++i) {
        int sel_size = selcross1_list[i].size();
        cout << selcross_split1[0] << i + 1 << ": ";
        if (sel_size != 1) {
            throw runtime_error("Your topology is not supported...");
        }
        for (int j = 0; j < sel_size; ++j) {
            cout << selcross1_list[i][j] << " "
                 << selcross2_list[i][j] << " "
                 << selcross3_list[i][j] << " "
                 << selcross4_list[i][j] << endl;
            selcross1_single.push_back(selcross1_list[i][j]);
            selcross2_single.push_back(selcross2_list[i][j]);
            selcross3_single.push_back(selcross3_list[i][j]);
            selcross4_single.push_back(selcross4_list[i][j]);
        }
    }
    //
    vector<int> selcross13_single, selcross24_single;
    for (auto sel1 : selcross1_single) {selcross13_single.push_back(sel1);}
    for (auto sel3 : selcross3_single) {selcross13_single.push_back(sel3);}
    for (auto sel2 : selcross2_single) {selcross24_single.push_back(sel2);}
    for (auto sel4 : selcross4_single) {selcross24_single.push_back(sel4);}
    //-------------------------------------------------------------------------

    cout << endl;
    cout << "STEP3: Check connect parts..." << endl;
    
    vector<int> change_list1, change_list2;
    //
    // polymer(1) and polymer(2) --------->
    s_pair polpair;
    for (int i = 0; i < Np; ++i) {
        int pol1_id = selpol1_single[i];
        int pol2_id = selpol2_single[i];
        vector<string> bonded_pol1_resids = GetBondedResidues(itp, pol1_id);
        int bonded_pol1_resids_size = bonded_pol1_resids.size();
        vector<string> bonded_pol2_resids = GetBondedResidues(itp, pol2_id);
        int bonded_pol2_resids_size = bonded_pol2_resids.size();
        //
        int nbond_pol1 = return_nbond_from_itp(itp, pol1_id);
        int nbond_pol2 = return_nbond_from_itp(itp, pol2_id);
        if ( (bonded_pol1_resids_size < nbond_pol1) 
             || (bonded_pol2_resids_size < nbond_pol2) ) {
            polpair.list1.push_back(pol1_id);
            polpair.list2.push_back(pol2_id);
        }
    }
     
    s_connect polpol = ReturnConnectionPairSame(itp, polpair, lastcoord, lastbox, rc); 
    
    int pp = polpol.joint1.size();
    if (pp != 0) {
        cout << endl;
        cout << "Case polymer(1) and polymer(2): " << endl;
        for (int i = 0; i < pp; ++i) {
            cout << polpol.joint1[i] << ", "
                 << polpol.joint2[i] << endl;
        }
        //
        cout << endl;
        cout << "Update topology: " << endl;
        UpdateTopology(polpol.joint1, polpol.joint2, itp, imit_itp);
         
        if (imit_itp.change_sign) {
            int nj1 = polpol.joint1.size();
            int nj2 = polpol.joint2.size();
            if (nj1 != nj2) throw runtime_error("nj1 != nj2");
            for (int n = 0; n < nj1; ++n) {
                change_list1.push_back(polpol.joint1[n]);
                change_list2.push_back(polpol.joint2[n]);
            }
        }
    } else {
        cout << endl;
        cout << "Case polymer(1) and polymer(2): Need not update..." << endl;
    }
    // <------------

    //////////////////////////////////////////////////////////////////
    ////// This method is slightly rough..., need debug(2026/02/27) //
    ////// Start debug: 2026/03/04~2026/03/23                       //
    //////////////////////////////////////////////////////////////////
   
    //
    //  crosslinked side, cross24 and pol2 ------->
    s_pair selcross24_pair;
    for (int i = 0; i < Nc; ++i) {
        //
        int cross1_id = selcross1_single[i];
        int cross2_id = selcross2_single[i];
        int cross3_id = selcross3_single[i];
        int cross4_id = selcross4_single[i];
        //
        vector<string> bonded_cross2_resids = GetBondedResidues(itp, cross2_id);
        int bonded_cross2_resids_size = bonded_cross2_resids.size();
        vector<string> bonded_cross4_resids = GetBondedResidues(itp, cross4_id);
        int bonded_cross4_resids_size = bonded_cross4_resids.size();
        //
        int nbond_cross2 = return_nbond_from_itp(itp, cross2_id);
        int nbond_cross4 = return_nbond_from_itp(itp, cross4_id);
        //
        if (bonded_cross2_resids_size < nbond_cross2) {
            selcross24_pair.list1.push_back(cross2_id);
            selcross24_pair.list2.push_back(cross1_id);
        }
        if (bonded_cross4_resids_size < nbond_cross4) {
            selcross24_pair.list1.push_back(cross4_id);
            selcross24_pair.list2.push_back(cross3_id);
        }
    }
    //
    s_connect cross24pol2 = ReturnConnectionPairDiff(itp, selcross24_pair, polpair.list2, 
                                                     lastcoord, lastbox, rc); 
    //
    int c24p2 = cross24pol2.joint1.size();
    if (c24p2 != 0) { 
        cout << endl;
        cout << "Case polymer(2), crosslinking mater. (2, 4): " << endl;
        vector<int> erase_c24p2_list;
        for (int i = 0; i < c24p2; ++i) {
            int pol24_id = cross24pol2.joint2[i];
            vector<string> bonded_resids = GetBondedResidues(itp, pol24_id);
            int bonded_resids_size = bonded_resids.size();
            int nbond_pol24 = return_nbond_from_itp(itp, pol24_id);
            int c24_id = cross24pol2.joint1[i];
            int p2_id = cross24pol2.joint2[i];
            bool neighbor_sign = false;
            
            int p2_pos = -1;
            for (int j = 0; j < Np; ++j) {
                int p2_id_tmp = selpol2_single[j];
                if (p2_id == p2_id_tmp) {
                    p2_pos = j;
                    break;
                }
            }

            int p1_id = selpol1_single[p2_pos];
            //cout << p1_id << " " << p2_id << endl;
            vector<string> p1_bonded_resids = GetBondedResidues(itp, p1_id);
            for (auto p1ri : p1_bonded_resids) {
                string p1ri_trim = remove_digits(p1ri);
                //cout << p1ri_trim << " ";
                if (p1ri_trim == selcross_split1[0]) {
                    neighbor_sign = true;
                    break;
                }
            }
            cout << endl;

            if (bonded_resids_size >= nbond_pol24) {
                erase_c24p2_list.push_back(i);
            } else if (neighbor_sign) {
                erase_c24p2_list.push_back(i);
            } else {
                cout << c24_id << ", " << p2_id << endl;
            }
        }
        cross24pol2 = erase_s_connect_by_list(cross24pol2, erase_c24p2_list);
        //
        cout << endl;
        cout << "Update topology: "<< endl;
        UpdateTopology(cross24pol2.joint1, cross24pol2.joint2, itp, imit_itp);
        //
        if (imit_itp.change_sign) {
            int nj1 = cross24pol2.joint1.size();
            int nj2 = cross24pol2.joint2.size();
            if (nj1 != nj2) throw runtime_error("nj1 != nj2");
            for (int n = 0; n < nj1; ++n) {
                change_list1.push_back(cross24pol2.joint1[n]);
                change_list2.push_back(cross24pol2.joint2[n]);
            }
        }
    } else {
        cout << endl;
        cout << "Case polymer(2), crosslinking mater. (2, 4): Need not update..." << endl;
    }
    // <-------
    
    // cross13 and pol1 ------>
    s_pair selcross13_pair;
    for (int i = 0; i < Nc; ++i) {
        int cross1_id = selcross1_single[i];
        int cross2_id = selcross2_single[i];
        int cross3_id = selcross3_single[i];
        int cross4_id = selcross4_single[i];
        //
        vector<string> bonded_cross1_resids = GetBondedResidues(itp, cross1_id);
        int bonded_cross1_resids_size = bonded_cross1_resids.size();
        vector<string> bonded_cross2_resids = GetBondedResidues(itp, cross2_id);
        int bonded_cross2_resids_size = bonded_cross2_resids.size();
        vector<string> bonded_cross3_resids = GetBondedResidues(itp, cross3_id);
        int bonded_cross3_resids_size = bonded_cross3_resids.size();
        vector<string> bonded_cross4_resids = GetBondedResidues(itp, cross4_id);
        int bonded_cross4_resids_size = bonded_cross4_resids.size();
        //
        int nbond_cross1 = return_nbond_from_itp(itp, cross1_id);
        int nbond_cross2 = return_nbond_from_itp(itp, cross2_id);
        int nbond_cross3 = return_nbond_from_itp(itp, cross3_id);
        int nbond_cross4 = return_nbond_from_itp(itp, cross4_id);
        //
        if ( (bonded_cross2_resids_size == nbond_cross2) 
             && (bonded_cross1_resids_size < nbond_cross1) ){
            selcross13_pair.list1.push_back(cross1_id);
            selcross13_pair.list2.push_back(cross2_id);
        }
        if ( (bonded_cross4_resids_size == nbond_cross4) 
             && (bonded_cross3_resids_size < nbond_cross3) ){
            selcross13_pair.list1.push_back(cross3_id);
            selcross13_pair.list2.push_back(cross4_id);
        }
    }
    //
    s_connect cross13pol1 = ReturnConnectionPairDiff(itp, selcross13_pair, polpair.list1, 
                                                     lastcoord, lastbox, rc); 
    //
    int c13p1 = cross13pol1.joint1.size();
    if (c13p1 != 0) {
        cout << endl;
        cout << "Case polymer(1), crosslinking mater. (1, 3): " << endl;
        vector<int> erase_c13p1_list;
        for (int i = 0; i < c13p1; ++i) {
            int c13_id = cross13pol1.joint1[i];
            int p1_id   = cross13pol1.joint2[i];
            //
            vector<string> bonded_resids = GetBondedResidues(itp, p1_id);
            int bonded_resids_size = bonded_resids.size();
            int nbond_pol13 = return_nbond_from_itp(itp, p1_id);
            //
            // neghbor check --->
            bool neighbor_sign = false;
            int p1_pos = -1;
            for (int j = 0; j < Np; ++j) {
                int p1_id_tmp = selpol1_single[j];
                if (p1_id == p1_id_tmp) {
                    p1_pos = j;
                    break;
                }
            }
            int p2_id = selpol2_single[p1_pos];
            //cout << "Neighbor selection: " << p1_id << " and " << p2_id << endl;
            vector<string> p2_bonded_resids = GetBondedResidues(itp, p2_id);
            //cout << itp.atoms.resid[p2_id-1] << ", " << itp.atoms.atom[p2_id-1] << ": ";
            //for (auto p2_br : p2_bonded_resids) cout << p2_br << " ";
            //cout << endl;
            for (auto p2ri : p2_bonded_resids) {
                string p2ri_trim = remove_digits(p2ri);
                if (p2ri_trim == selcross_split1[0]) {
                    neighbor_sign = true;
                    break;
                }
            }
            // <----

            if (bonded_resids_size >= nbond_pol13) {
                erase_c13p1_list.push_back(i);
            } else if (neighbor_sign) {
                erase_c13p1_list.push_back(i);
            } else {
                // final check ------>
                bool n_sign = false;
                int idx_pol13   = find(polpair.list1.begin(), polpair.list1.end(), p1_id) - polpair.list1.begin();
                int idx_cross13 = find(selcross13_single.begin(), selcross13_single.end(), c13_id) - selcross13_single.begin();  
                int p1_id_connect  = polpair.list2[idx_pol13];
                int c13_id_connect = selcross24_single[idx_cross13];
                //
                vector<string> bonded_p13ic_resids = GetBondedResidues(itp, p1_id_connect);
                vector<string> bonded_c13ic_resids = GetBondedResidues(itp, c13_id_connect);
                string bonded_c13ic_resid_back = bonded_c13ic_resids.back();
                //
                for (auto b : bonded_p13ic_resids) {
                    if (b == bonded_c13ic_resid_back) {
                        n_sign = true;
                        break;
                    }
                }
                // <--------
                if (n_sign) {
                    erase_c13p1_list.push_back(i);
                } else {
                    cout << c13_id << ", " << p1_id << endl;
                }
            }
        }
        cross13pol1 = erase_s_connect_by_list(cross13pol1, erase_c13p1_list);
        int c13p1_2 = cross13pol1.joint1.size();
        if (c13p1_2 != 0) {
            cout << endl;
            cout << "Update topology: "<< endl;
            UpdateTopology(cross13pol1.joint1, cross13pol1.joint2, itp, imit_itp);
            if (imit_itp.change_sign) {
                int nj1 = cross13pol1.joint1.size();
                int nj2 = cross13pol1.joint2.size();
                if (nj1 != nj2) throw runtime_error("nj1 != nj2");
                for (int n = 0; n < nj1; ++n) {
                    change_list1.push_back(cross13pol1.joint1[n]);
                    change_list2.push_back(cross13pol1.joint2[n]);
                }
            }
        } else {
            cout << endl;
            cout << "Case polymer(1), crosslinking mater. (1, 3): Need not update..." << endl;
        }
    } else {
        cout << endl;
        cout << "Case polymer(1), crosslinking mater. (1, 3): Need not update..." << endl;
    }
    // <-------
    if (imit_itp.change_sign) {
        cout << "==============================================================" << endl; 
        int cl1 = change_list1.size();
        if (cl1 > 0) {
            cout << endl;
            cout << "# of change topology pairs: " << cl1 << endl;
            CHANGE_TOPOLOGY(itp, imit_itp, change_list1, change_list2);
        } else {
            cout << "Not change topology..." << endl;
        }
        cout << "==============================================================" << endl; 
    }
   
    cout << endl;
    cout << "CHARGE SECTION: " << endl;
    //
    int nc = chargefile.size();
    if (nc != 4) throw runtime_error("# of chargefile != 4");
    //
    cout << endl;
    cout << "Monomer side ==>" << endl;
    s_charge charge_pol = ChangeChargeSumZero(chargefile[0]);
    UpdateChargeXlinkerMonomer(itp, charge_pol, charge_pol, 
                               polpol.joint1, polpol.joint2);
    //
    cout << endl;
    cout << "Xlinked side ==>" << endl; 
    //cout << "Bond check-->" << endl;
    vector<vector<int>> 
    sel_bonded_list = ReturnBondedListFour(itp, selcross1_single, selcross2_single, 
                                           selcross3_single, selcross4_single, 
                                           selpol_split1[0]);
    ////
    //double charge_sum_old = 0.0;
    //for (auto c : itp.atoms.charge) charge_sum_old += c;
    //cout << "Total charge (old): " << charge_sum_old << endl;
    ////
    for (int i = 0; i < Nc; ++i) {
        string x_resid = combine_string_int(selcross_split1[0], i+1);
        vector<int> x_list = GetIndexResid(itp, x_resid);
        vector<int> sel_bonded = sel_bonded_list[i];
        //for (auto x : x_list) cout << x << " ";
        //cout << endl;
        //
        cout << endl;
        cout << x_resid << ": ";
        vector<string> chargefile_list = {chargefile[1]};
        //vector<string> chargefile_list;
        //int fcount = 0;
        for (int b = 0; b < 4; ++b) {
           string chargefile_void = "";
           int sb = sel_bonded[b];
           cout << sb << " ";
           //if (sb != -1) fcount += 1;
           if (sb != -1 && b%2 == 0) chargefile_list.push_back(chargefile[2]);
           if (sb != -1 && b%2 != 0) chargefile_list.push_back(chargefile[3]);
           if (sb == -1) chargefile_list.push_back(chargefile_void);
        }
        //if (fcount == 4) {
        //    chargefile_list = {chargefile[1], chargefile[2], chargefile[3], 
        //                       chargefile[2], chargefile[3]};
        //} else {
        //    continue;
        //}
        //cout << "Size of chargefile_list: " << chargefile_list.size() << endl;
        vector<s_charge> charge_list = ChangeChargeListSumZero(chargefile_list);
        if (charge_list.size() != 5) throw runtime_error("size of charge_list != 5");
        //
        //bool change_sign = false;
        for (int i = 4; i > -1; --i) {
            //cout << i << ",    " << endl;
            s_charge charge = charge_list[i];
            if (charge.nlist == 0) continue;
            if (i == 0) {
                //if (!change_sign) break;
                //cout << "Update xlinked side charge..." << endl;
                ChangeCharge(itp, x_list, charge);
            } else {
                //if (charge.nlist == 0) continue;
                //change_sign = true;
                //cout << "Update polymer side charge (" << i << ")" << endl;
                int sb = sel_bonded[i-1];
                vector<int> m_list = ReturnSameResidIndex(itp, sb);
                ChangeCharge(itp, m_list, charge);
            }
        }
        //double charge_sum_tmp = 0.0;
        //for (auto c : itp.atoms.charge) charge_sum_tmp += c;
        //cout << "==> Total charge (tmp): " << charge_sum_tmp << endl;
    }
    //
    double sum_charge = 0.0;
    for (auto c : itp.atoms.charge) sum_charge += c;
    cout << endl;
    cout << endl;
    cout << "Total charge: " << sum_charge << endl;
    //
    cout << endl;
    cout << "STEP4: Output itpfile..." << endl;
    //
    Writeitp(itp, outitp);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
////////////////////////// ANALYZE_REMOVE_IONS ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_REMOVE_IONS(const string& itpfile, 
                         const string& pdbfile, 
                         const string& outitp, 
                         const string& outpdb, 
                         const vector<string>& ions_string) {
    s_itp itp;
    Readitpfile(itpfile, itp);
    //
    int natoms = itp.atoms.resid.size();
    vector<int> ions_index;
    for (int i = 0; i < natoms; ++i) {
        string name = itp.atoms.atom[i];
        for (auto name_ions : ions_string) {
            if (name == name_ions) {
                ions_index.push_back(i+1);
            }
        }
    }
    //
    cout << endl;
    cout << "REMOVE INDEX: " << endl;
    for (auto idx : ions_index) {
        cout << idx << " ";
    }
    cout << endl;
    cout << "# of total ions: " << ions_index.size() << endl;
    //
    RemoveAtomsKeepIndex(itp, ions_index);
    CompactAtomIndices(itp);
    //
    Writeitp(itp, outitp);
    //
    PDBOPR PdbOpr;
    s_pdb pdb = PdbOpr.LoadfrompdbAmber(pdbfile);
    s_pdb remove_pdb = PdbOpr.Removepdb(pdb, ions_index);
    PdbOpr.OutPdbAmber({remove_pdb}, outpdb);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// ANALYZE_CL_POSITION ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_CL_POSITION(const string& itpfile, 
                         const vector<string>& trajlist, 
                         const vector<string>& selpol, 
                         const vector<string>& selcross, 
                         const int Nc, 
                         const vector<double>& Min,
                         const vector<double>& Max,
                         const int& grid) {
    cout << endl;
    cout << "STEP1: Read itpfile..." << endl;
    s_itp itp;
    Readitpfile(itpfile, itp);
    int natoms = itp.atoms.resnr.size();
    ////
    //s_xtc xtc;
    //if (boost::iequals(get_extension(trajfile), "xtc")) {
    //    xtc = ReadXTC(trajfile);
    //}
    ////
    cout << endl;
    cout << "STEP2: Search crosslinking mater. " << endl;
    //--------------------------------------------------------------------------
    // For crosslinking material part...
    int c = selcross.size() - 1;
    vector<vector<int>> selcross_list = ExtractList(itp, Nc, natoms, selcross);
    //--------------------------------------------------------------------------
    vector<int> resid_list;
    vector<vector<int>> bonded_list;
    vector<vector<int>> cross_id_list;
    for (int i = 0; i < Nc; ++i) {
        if (selcross_list[i].size() > 0) {
            int ncount = 0;
            vector<string> polname;
            vector<int> cross_id_cell;
            for (int j = 0; j < c; ++j) {
                int cross_id = selcross_list[i][j];
                vector<string> bonded_resids = GetBondedResidues(itp, cross_id);
                vector<string> resid_trim = TrimResid(bonded_resids);
                //
                for (auto s : bonded_resids) {
                    if (s.find(selpol[0]) != string::npos) {
                        polname.push_back(s);
                        cross_id_cell.push_back(cross_id);
                    }
                }
                if (count(resid_trim.begin(), resid_trim.end(), selpol[0]) == 1) {
                    ncount += 1;
                    vector<int> bonded_cell = GetBondedIndex(itp, cross_id);
                    bonded_cell.push_back(cross_id);
                    bonded_list.push_back(bonded_cell);
                }
            }
            if (ncount == 2) {
                if (polname[0] != polname[1]) {
                    resid_list.push_back(i+1);
                    cross_id_list.push_back(cross_id_cell);
                }
            }
        }
    }
    //
    cout << endl;
    cout << "Double cross-linked selection: " << endl;
    int ncl = resid_list.size();
    for (int i = 0; i < ncl; ++i) {
        cout << selcross[0] << resid_list[i] << ": ";
        cout << cross_id_list[i][0] << " " << cross_id_list[i][1] << endl;
    }
    //
    cout << endl;
    cout << "STEP3: Calc. end-to-end dist." << endl;
    double r_avg = 0.0;
    double r2_avg = 0.0;
    double v_avg = 0.0;
    int    tot_frames = 0;
    vector<double> r_list;
    vector<vector<double>> cl_cell = make_3d_grid(Min, Max, grid);
    vector<double> cl_avg(grid * grid * grid, 0);
    //
    for (auto trajfile : trajlist) {
        s_xtc xtc;
        if (boost::iequals(get_extension(trajfile), "xtc")) {
            cout << "Read " << trajfile << "..." << endl;
            xtc = ReadXTC(trajfile);
        } else {
            cout << trajfile << " format is not supported.." << endl;
        }
        for (int f = 0; f < xtc.nframes; ++f) {
            vector<double> box = {xtc.box[f][0], xtc.box[f][4], xtc.box[f][8]};
            for (int i = 0; i < ncl; ++i) {
                int cross_id1 = cross_id_list[i][0];
                int cross_id2 = cross_id_list[i][1];
                vector<double> vec1 = xtc.coord[f][cross_id1-1];
                vector<double> vec2 = xtc.coord[f][cross_id2-1];
                double r2 = dist2_pbc(vec1, vec2, box);
                double r  = sqrt(r2);
                r2_avg += r2; 
                r_avg  += r;
                r_list.push_back(r);
                //
                // Calc. midpoint
                vector<double> mid = pbc_midpoint_3d(vec1, vec2, box);
                vector<int> index_vec = grid_index_3d(mid, Min, Max, grid);
                cl_avg[idx(index_vec, grid)]++;
            }
            v_avg += box[0] * box[1] * box[2];
        }
        tot_frames += xtc.nframes;
    }
    //
    r2_avg /= tot_frames * ncl;
    r_avg  /= tot_frames * ncl;
    double del_r2 = r2_avg - r_avg * r_avg;
    v_avg /= tot_frames;
    //
    int cl_size = cl_avg.size();
    for (int c = 0; c < cl_size; ++c) {
        cl_avg[c] /= tot_frames;
    }
    //
    cout << endl;
    cout << "STEP4: Output result to results.dat and hist.dat and grid.dat..." << endl;
    cout << endl;
    ofstream outfile("results.dat");
    outfile << "# of cross-link points: " << ncl    << endl;
    outfile << "<r>:                    " << r_avg  << endl;
    outfile << "<r**2>:                 " << r2_avg << endl;
    outfile << "<(delta r)**2>:         " << del_r2 << endl;
    outfile << "<V>:                    " << v_avg  << endl;
    outfile.close();
    //
    auto h = make_histgram(r_list, 100, NormType::PDF);
    SplitdataToFile("hist.dat", h.bin_center, h.value);
    //
    write_3d_hist_gnuplot("grid.dat", cl_avg, Min, Max, grid);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////// ANALYZE_HOPPING /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_HOPPING(const string& topfile,
                     const vector<string>& trajlist, 
                     const vector<string>& selion, 
                     const vector<string>& selmemb, 
                     const double dt, 
                     const double rc, 
                     const int rst_sign) {
    //
    // Resname: selion & selmemb
    cout << endl;
    cout << "STEP1: select ions and polymer-membrane" << endl;
    s_joint ion  = MakeJoint(selion);
    s_joint memb = MakeJoint(selmemb);
    //
    cout << endl;
    cout << "STEP2: Read topology file..." << endl;
    s_top top;
    Readtopfile(topfile, top);

    vector<int> ion_index, memb_index;
    int ion_seq=-1, memb_seq=-1;
    
    for (auto itp : top.itp_list) {
        ion_index = ExtractfromITP(itp, ion);
        int ni = ion_index.size();
        ion_seq += 1;
        if (ni > 0) break; 
    }
    for (auto itp : top.itp_list) {
        memb_index = ExtractfromITP(itp, memb);
        int nm = memb_index.size();
        memb_seq += 1;
        if (nm > 0) break; 
    }
    if (ion_seq == -1 || memb_seq == -1) throw runtime_error("itp is not included in topfile...");

    cout << endl;
    cout << "MEMBRANE selection: ";
    vector<int> memb_selection;
    int nbf_memb = 0;
    for (int i = 0; i < memb_seq; ++i) { nbf_memb += top.molecules.natom[i]; }
    for (auto memb_int : memb_index) {
        int memb_idx = memb_int + nbf_memb;
        memb_selection.push_back(memb_idx);
        cout << memb_idx << " ";
    }
    cout << endl;
    //
    int nbf_ion = 0;
    for (int i = 0; i < ion_seq; ++i) { nbf_ion += top.molecules.natom[i]; }
    //
    int Ni = top.molecules.nmol[ion_seq];

    vector<int> ion_selection;
    if (ion_index.size() == 1) {
        cout << endl;
        cout << "ION selection: ";
        int i_temp = nbf_ion + ion_index[0];
        for (int i = 0; i < Ni; ++i) {
            cout << i_temp << " ";
            ion_selection.push_back(i_temp);
            i_temp += top.itp_list[ion_seq].atoms.n;
        }
        cout << endl;
    } else {
        cout << endl;
        cout << "index of ions is not 1 (this analysis is not supported in this case)" << endl;
        return;
    }
    //
    cout << endl;
    cout << "STEP2: Read traj. and analyze..." << endl;
    int Nm = memb_selection.size();
    vector<vector<double>> length_list;
    vector<vector<int>>    memb_list;
    int traj_count = 0;
    vector<double> time;
    double t_temp = 0.0;
    int tot_step = 0;
    //
    for (auto trajfile : trajlist) {
        cout << "Read " << trajfile << "..." << endl;
        s_xtc xtc = ReadXTC(trajfile);
        int restart = 0; 
        if (rst_sign == 1) {
           restart = (traj_count == 0) ? 0 : 1;
        }
        traj_count += 1;
        for (int f = restart; f < xtc.nframes; ++f) {
            vector<double> box = {xtc.box[f][0], xtc.box[f][4], xtc.box[f][8]};
            tot_step += 1;
            time.push_back(t_temp);
            t_temp += dt;
            if (fmod(tot_step, 100) == 0) {
                cout << "STEP: " << tot_step << endl;
            }
            // 
            vector<double> length_cell;
            vector<int>    memb_cell;
            for (int i = 0; i < Ni; ++i) {
                int ion_id = ion_selection[i]-1;
                int memb_id_min = -1;
                vector<double> i_coord = xtc.coord[f][ion_id];
                //
                double r_min = 1000.0;
                for (int m = 0; m < Nm; ++m) {
                    int memb_id = memb_selection[m]-1;
                    vector<double> m_coord = xtc.coord[f][memb_id];
                    double r2 = dist2_pbc(i_coord, m_coord, box);
                    double r = sqrt(r2);
                    if (r < r_min) {
                        r_min = r;
                        memb_id_min = memb_id + 1;
                    }
                }
                length_cell.push_back(r_min);
                memb_cell.push_back(memb_id_min);
            }
            length_list.push_back(length_cell);
            memb_list.push_back(memb_cell);
        }
    }
    // 
    //
    cout << endl;
    cout << "Output data to length.dat, index.dat..." << endl;
    TimedataXdToFile("length.dat", time, length_list);
    TimeINTdataXdToFile("index.dat", time, memb_list);
    //

    //int tot_step = time.size();
    cout << "Total step: " << tot_step << endl;
    //vector<vector<double>> tcf_den(Ni, vector<double>(Nm, 0.0));
    //vector<vector<vector<double>>> 
    //tcf_num(tot_step, vector<vector<double>>(Ni, vector<double>(Nm, 0.0)));
    double tcf_den = 0.0;
    vector<double> tcf_num(tot_step, 0.0);
    
    cout << endl;
    cout << "==> Calculate TCF..." << endl;
    //
    for (int f = 0; f < tot_step; ++f) {
        if ( fmod(f+1, 100) == 0) {
            cout << "STEP: " << f + 1 << endl;
        }
        for (int i = 0; i < Ni; ++i) {
            int memb_id_f =  memb_list[f][i];
            //int n = index_of(memb_selection, memb_id_f);
            //cout << "NUMBER: " << n << endl;
            if (length_list[f][i] < rc) {
                tcf_den += 1.0;
                // 
                for (int d = f; d < tot_step; ++d) {
                    int t = d - f;
                    int memb_id_d = memb_list[d][i];
                    if ( (memb_id_f == memb_id_d) &&
                         (length_list[d][i] < rc) ) {
                        tcf_num[t] += 1.0;
                    }
                }
            }
        }
    }

    // set avg.
    tcf_den /= tot_step;
    
    for (int f = 0; f < tot_step; ++f) {
        tcf_num[f] /= tot_step - f;
    }

    vector<double> tcf, tcf_time;
    for (int f = 0; f < tot_step; ++f) {
        tcf_time.push_back(dt*f);
        double tcf_temp = tcf_num[f] / tcf_den;
        tcf.push_back(tcf_temp);
    }
    //
    cout << endl;
    cout << "Output data to tcf.dat..." << endl;
    SplitdataToFile("tcf.dat", tcf_time, tcf);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////// ANALYZE_SD //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_SD(const string& topfile, 
                const string& trajfile, 
                const vector<string>& sel, 
                const int Next,
                const double dt) {
    cout << endl;
    cout << "STEP1: Selection part..." << endl;
    s_joint j_sel = MakeJoint(sel);
    //
    cout << endl;
    cout << "STEP2: Read topology file..." << endl;
    s_top top;
    Readtopfile(topfile, top);
    //
    vector<int> sel_selection = GetSelselection(top, j_sel);

    vector<int> random_list = random_sample(sel_selection, Next);
    cout << "Selection: ";
    for (auto r : random_list) { cout << r << " ";}
    cout << endl;

    cout << endl;
    cout << "STEP3: Calculate squere deviation..." << endl;
    s_xtc xtc;
    if (boost::iequals(get_extension(trajfile), "xtc")) {
        cout << "Read " << trajfile << "..." << endl;
        //xtc = ReadXTC(trajfile);
        xtc = ReadXTC_choose(trajfile, random_list);
        cout << xtc.coord[0].size() << endl;
    } else {
        cout << "Your trajectory format is not supported, sorry..." << endl;
        return;
    }
   
    vector<vector<double>> sd;
    vector<double> times;
    
    for (int f = 0; f < xtc.nframes; ++f) {
        times.push_back(f*dt);
    }

    for (int i = 0; i < Next; ++i) {
        vector<double> coord_origin;
        vector<double> sd_cell;
        for (int f = 0; f < xtc.nframes; ++f) {
            //
            vector<double> coord = xtc.coord[f][i];
            if (f == 0) {
                coord_origin = coord;
            }
            double r2 = dist2(coord, coord_origin);
            sd_cell.push_back(r2);
        }
        sd.push_back(sd_cell);
    }

    cout << endl;
    cout << "STEP4: Output result to files..." << endl;
    vector<string> outfiles = make_zero_padding_filenames(Next);
    for (int i = 0; i < Next; ++i) {
        string outfile = outfiles[i];
        SplitdataToFile(outfile, times, sd[i]);
    }
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// ANALZE_COORD_NUM /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_COORD_NUM(const string& topfile, 
                       const vector<string>& trajlist, 
                       const vector<string>& sel1, 
                       const vector<string>& sel2, 
                       const double& rc) {
    cout << endl;
    cout << "STEP1: Read topology file..." << endl;
    s_joint j_sel1 = MakeJoint(sel1);
    s_joint j_sel2 = MakeJoint(sel2);
    //
    s_top top;
    Readtopfile(topfile, top);
    //
    vector<int> sel1_selection = GetSelselection(top, j_sel1);
    cout << endl;
    cout << "Selection (sel1): ";
    for (auto s1 : sel1_selection) {cout << s1 << " ";}
    cout << endl;
    //
    vector<int> sel2_selection = GetSelselection(top, j_sel2);
    cout << endl;
    cout << "Selection (sel2): ";
    for (auto s2 : sel2_selection) {cout << s2 << " ";}
    cout << endl;
    
    cout << endl;
    cout << "STEP2: Start analysis..." << endl;
    //
    vector<double> coord_num(100, 0.0);
    for (auto trajfile : trajlist) {
        s_xtc xtc;
        if (boost::iequals(get_extension(trajfile), "xtc")) {
            cout << "Read " << trajfile << "..." << endl;
            xtc = ReadXTC(trajfile);
        } else {
           cout << trajfile << " format is not supported.." << endl;
        }
        //
        for (int f = 0; f  < xtc.nframes; ++f) {
            vector<double> box = {xtc.box[f][0], xtc.box[f][4], xtc.box[f][8]};
            for (auto s1 : sel1_selection) {
                int ncount = 0;
                vector<double> sel1_coord = xtc.coord[f][s1-1];
                for (auto s2 : sel2_selection) {
                    vector<double> sel2_coord = xtc.coord[f][s2-1];
                    double length = dist2_pbc(sel1_coord, sel2_coord, box);
                    length = sqrt(length);
                    if (length < rc) ncount += 1;
                }
                coord_num[ncount] += 1; 
            }
        }
    }
    //
    normalize_list(coord_num); 
    SplitdataToFileNUMdoubleZero("cn.dat", coord_num);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////// ANALYZE_END_TO_END ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_END_TO_END(const string& topfile, 
                        const vector<string>& trajlist,
                        const vector<string>& selpol, 
                        const int Np, 
                        const int npol, 
                        const int nbin) {
    cout << endl;
    cout << "STEP1: Read topology file..." << endl;
    s_joint j_end1 = MakeJoint({selpol[0]});
    s_joint j_end2 = MakeJoint({selpol[1]});
    //
    s_top top;
    Readtopfile(topfile, top);
    //
    vector<int> end1_selection = GetSelselection(top, j_end1);
    vector<int> end2_selection = GetSelselection(top, j_end2);
    //
    int end1_size = end1_selection.size();
    int end2_size = end2_selection.size();
    int ntot = Np * npol;
    if (end1_size != end2_size) throw runtime_error("end1_size != end2_size");
    if (ntot != end1_size) throw runtime_error("ntot != end1_size");
    //
    vector<vector<int>> end1_split, end2_split;
    //
    int ncount = 0;
    for (int i = 0; i < npol; ++i) {
        vector<int> end1_cell, end2_cell;
        for (int j = 0; j < Np; ++j) {
            end1_cell.push_back(end1_selection[ncount]);
            end2_cell.push_back(end2_selection[ncount]);
            ncount += 1;
        }
        end1_split.push_back(end1_cell);
        end2_split.push_back(end2_cell);
    }
    //
    cout << endl;
    cout << "Selection==> " << endl;
    //
    for (int i = 0; i < npol; ++i) {
        int e1 = end1_split[i][0];
        int e2 = end2_split[i][Np-1];
        cout << e1 << " and " << e2 << endl;
    }

    cout << endl;
    cout << "STEP2: Start Calc." << endl;
    
    vector<double> length_list;
    for (auto trajfile : trajlist) {
        s_xtc xtc;
        if (boost::iequals(get_extension(trajfile), "xtc")) {
            cout << "Read " << trajfile << "..." << endl;
            xtc = ReadXTC(trajfile);
        } else {
           cout << trajfile << " format is not supported.." << endl;
        }
        //
        for (int f = 0; f  < xtc.nframes; ++f) {
            vector<double> box = {xtc.box[f][0], xtc.box[f][4], xtc.box[f][8]};
            for (int i = 0; i < npol; ++i) {
                int e1 = end1_split[i][0];
                int e2 = end2_split[i][Np-1];
                vector<double> e1_coord = xtc.coord[f][e1-1];
                vector<double> e2_coord = xtc.coord[f][e2-1];
                double length = dist2_pbc(e1_coord, e2_coord, box);
                length = sqrt(length);
                length_list.push_back(length);
                //cout << e1 << " and " << e2 << " == " << length << endl;
            }
        }
    }
    //
    cout << endl;
    cout << "STEP3: Output data to file..." << endl;
    auto h = make_histgram(length_list, nbin, NormType::PDF);
    SplitdataToFile("hist.dat", h.bin_center, h.value);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
////////////////////////////// ANALYZE MASS ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ANALYZE_MASS(const string& itpfile) {
    cout << endl;
    cout << "STEP1: Read itp file..." << endl;
    s_itp itp;
    Readitpfile(itpfile, itp);
    double sum_mass = 0.0;
    for (double m : itp.atoms.mass) {
        sum_mass += m;
    }
    cout << endl;
    cout << "STEP2: Output result to mass.dat" << endl;
    cout << endl;
    ofstream outfile("mass.dat");
    outfile << "mass: " << fixed << setprecision(3) << sum_mass << endl;
    outfile.close();
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
////////////////////////////// CHANGE_CHARGE //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void CHANGE_CHARGE(const string& itpfile, 
                   const string& outitp, 
                   const vector<string>& selmemb, 
                   const double& aft_charge) {
    cout << endl;
    cout << "STEP1: Read itp file..." << endl;
    s_itp itp;
    Readitpfile(itpfile, itp);
    s_joint sel_joint = MakeJoint(selmemb);
    vector<int> sellist = ExtractfromITP(itp, sel_joint);
    //
    cout << endl;
    cout << "STEP2: Change charge..." << endl;
    itp = return_s_itp_change_charge(itp, sellist, aft_charge);
    //
    cout << endl;
    cout << "STEP3: Output itpfile..." << endl;
    Writeitp(itp, outitp);
}
/////////////////////////;/////////////////////////////////////////////////////
