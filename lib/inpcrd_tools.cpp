// inpcrd_tools.cpp

#include <inpcrd_tools.hpp>

s_inpcrd INPCRDOPR::ReadInfo(const string& filename) {
    //
    s_inpcrd inpcrd;
    //
    ifstream file(filename);
    if (file.is_open()) {
        string firstline;
        string secondline;
        //double x1, y1, z1, x2, y2, z2;
        //
        getline(file, firstline);
        getline(file, secondline);
        //
        inpcrd.title = firstline;
        inpcrd.natm  = stoi(secondline);
        //
        double x, y, z;
        int i = 0;
        while (file >> x >> y >> z) {
            if (i < inpcrd.natm) {
                inpcrd.coord.push_back({x, y, z});
                i += 1;
            } else if (i == inpcrd.natm) {
                inpcrd.box = {x, y, z};
                i += 1;
            } else if (i == inpcrd.natm + 1) {
                inpcrd.angle = {x, y, z};
            } else {
                cout << "This program does not support this data: " 
                     << x << " " << y << " " << z << endl;
            }
        }
        //
    } else {
        cerr << "Unable to open file: " << filename << endl;
    }
    file.close();

    return inpcrd;
}

s_inpcrd INPCRDOPR::PdbToInpcrd(const s_pdb& pdb) {
    //
    s_inpcrd inpcrd;
    int natm = pdb.coord.size();

    inpcrd.title = "default_name";
    inpcrd.natm  = natm;

    for (int i = 0; i < natm; ++i) {
        double x = pdb.coord[i][0];
        double y = pdb.coord[i][1];
        double z = pdb.coord[i][2];
        inpcrd.coord.push_back({x, y, z});
    }

    return inpcrd;
}

void INPCRDOPR::OutInpcrd(const vector<s_inpcrd>& sysinpcrd, 
                          const string& filename) {
    //
    s_inpcrd inpcrd;
    int natm = 0;
    //
    inpcrd.box = {sysinpcrd[0].box[0], sysinpcrd[0].box[1], sysinpcrd[0].box[2]};
    inpcrd.angle = {sysinpcrd[0].angle[0], sysinpcrd[0].angle[1], sysinpcrd[0].angle[2]};
    //
    for (s_inpcrd cellinpcrd : sysinpcrd) {
         natm += cellinpcrd.natm;
         for (int i = 0; i < cellinpcrd.natm; ++i) {
             double x = cellinpcrd.coord[i][0];
             double y = cellinpcrd.coord[i][1];
             double z = cellinpcrd.coord[i][2];
             inpcrd.coord.push_back({x, y, z});
         }
    }
    inpcrd.natm  = natm;
    inpcrd.title = "default_name";
    //
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    cout << "output inpcrd data to inpcrdfile..." << endl;

    outfile << inpcrd.title << endl;
    outfile << inpcrd.natm  << endl;
    //
    bool cont = true;
    for (int i = 0; i < natm; ++i) {
        if (cont) {
            outfile << fixed << setprecision(7) << setw(12) << inpcrd.coord[i][0] 
                    << fixed << setprecision(7) << setw(12) << inpcrd.coord[i][1] 
                    << fixed << setprecision(7) << setw(12) << inpcrd.coord[i][2];
            cont = false;
        } else {
            outfile << fixed << setprecision(7) << setw(12) << inpcrd.coord[i][0] 
                    << fixed << setprecision(7) << setw(12) << inpcrd.coord[i][1] 
                    << fixed << setprecision(7) << setw(12) << inpcrd.coord[i][2] << endl;
            cont = true;
        }
    }
    
    outfile << fixed << setprecision(7) << setw(12) << inpcrd.box[0] 
            << fixed << setprecision(7) << setw(12) << inpcrd.box[1] 
            << fixed << setprecision(7) << setw(12) << inpcrd.box[2] 
            << fixed << setprecision(7) << setw(12) << inpcrd.angle[0]
            << fixed << setprecision(7) << setw(12) << inpcrd.angle[1]
            << fixed << setprecision(7) << setw(12) << inpcrd.angle[2];

    return;
}

void INPCRDOPR::PlungeBoxAngle(const s_inpcrd& inpcrd_origin, 
                               s_inpcrd& inpcrd_out) {
    //
    inpcrd_out.box   = {inpcrd_origin.box[0], inpcrd_origin.box[1], inpcrd_origin.box[2]};
    inpcrd_out.angle = {inpcrd_origin.angle[0], inpcrd_origin.angle[1], inpcrd_origin.angle[2]};
    //
    return;
}
