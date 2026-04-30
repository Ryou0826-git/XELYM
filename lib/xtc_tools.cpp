// xtc_tool.cpp

#include <xtc_tools.hpp>

s_xtc ReadXTC(const string &filename) {
    s_xtc xtc_data;

    int natoms;
    if (read_xtc_natoms(const_cast<char*>(filename.c_str()), &natoms) != exdrOK) {
        cerr << "Error: cannot read number of atoms from " << filename << endl;
        return xtc_data;
    }
    xtc_data.natoms = natoms;
    
    XDRFILE *xd = xdrfile_open(const_cast<char*>(filename.c_str()), "r");
    if (!xd) {
        cerr << "Error: cannot open " << filename << endl;
        return xtc_data;
    }

    int step;
    float time;
    matrix box;
    rvec *x = new rvec[natoms];
    float prec;

    while (true) {
        int status = read_xtc(xd, natoms, &step, &time, box, x, &prec);
        if (status != exdrOK) break;

        xtc_data.nframes++;
        xtc_data.time.push_back(static_cast<double>(time));

        vector<double> box_flat(9);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                box_flat[i * 3 + j] = static_cast<double>(box[i][j] * 10.0);
        xtc_data.box.push_back(box_flat);

        vector<vector<double>> frame(natoms, vector<double>(3));
        for (int i = 0; i < natoms; ++i) {
            frame[i][0] = x[i][0] * 10.0; // nm --> \AA;
            frame[i][1] = x[i][1] * 10.0; // nm --> \AA;
            frame[i][2] = x[i][2] * 10.0; // nm --> \AA;
        }
        xtc_data.coord.push_back(frame);
    }

    delete[] x;
    xdrfile_close(xd);

    return xtc_data;
}
