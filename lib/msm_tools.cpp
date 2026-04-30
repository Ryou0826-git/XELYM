// msm_tools.hpp

#include "msm_tools.hpp"

//---------------------------------------------------------------------------//
Eigen::MatrixXd GainTrMat(const vector<int>& data,
                          const int& tau) {

    //
    // Search number of regions...
    int dim = -1;
    for (int nreg : data) {
        if (nreg > dim) {
            dim = nreg;
        }
    }

    int totstep = data.size();
    Eigen::MatrixXd TrMat(dim, dim);
    //
    for (int st = 0; st < totstep - tau; ++st) {
        int it_i = data[st] - 1;
        int it_j = data[st+tau] - 1;
        //
        TrMat(it_i, it_j) += 1.0;
    }

    Eigen::VectorXd SumTr = TrMat.rowwise().sum();
    for (int i = 0; i < dim; ++i) {
        if (SumTr(i) > 0) {
            TrMat.row(i) /= SumTr(i);
        }
    }
    Eigen::MatrixXd TrMat_tr = TrMat.transpose();

    return TrMat_tr;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_msm PerformMSM(const Eigen::MatrixXd& TrMat_norm,
                 const double& kT) {
    //
    bool fe_sign = false;
    s_msm msm;
    msm.TrMat = TrMat_norm;
    // 
    int dim = TrMat_norm.rows();
    vector<double> fe(dim, 0.0);

    auto [values, Pss] = PerformDiag(TrMat_norm);

    for (int i = 0; i < dim; ++i) {
        double pop = Pss(i).real() / Pss.real().sum();
        msm.pop.push_back(pop);
        if (Pss(i).real() > 0) {
            double fe = - kT * log(pop);
            msm.fe.push_back(fe);
            fe_sign = true;
        } else {
            cout << "Warning: " << Pss(i).real() << " < 0" << endl;
        }
    }

    double fe_min = *min_element(begin(msm.fe), end(msm.fe));
    if (fe_sign == true) {
        for (int i = 0; i < dim; ++i) {
            msm.fe[i] -= fe_min;
        }
    } else {
        cout << "Error: Free energy is not defined!!" << endl;
    }

    return msm;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void OutMSMresult(const s_msm& msm, 
                  const string& filename) {
    ofstream outfile(filename);
    //
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    } else {
        cout << "Output msm result to " << filename << endl; 
    }
    //
    outfile << "Trans Matrix: " << endl;
    outfile << msm.TrMat << endl;
    //
    outfile << endl;
    outfile << "Population: " << endl;
    int i = 0;
    for (double pop : msm.pop) {
        outfile << "state" << i+1 << " -> " << pop << endl;
        i += 1;
    }
    
    outfile << endl;
    outfile << "Free energy (kcal/mol): " << endl;
    int j = 0;
    for (double fe : msm.fe) {
        outfile << "state" << j+1 << " -> " << fe << endl;
        j += 1; 
    }
}
//---------------------------------------------------------------------------//
