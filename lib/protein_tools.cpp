// protein_tools.cpp

#include <protein_tools.hpp>

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// General Tools /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
vector<vector<vector<double>>>
ExtractCoordRes(const vector<vector<vector<double>>>& coord,
                const vector<int>& reslist) {
    //
    int frames = coord.size();
    int nres   = reslist.size();
    //
    vector<vector<vector<double>>>
    coordres(frames, vector<vector<double>>(nres, vector<double>(3, 0.0)));
    //
    for (int i = 0; i < frames; ++i) {
        int l = 0;
        for (int j = 0; j < nres; ++j) {
            int label = reslist[j];
            for (int k = 0; k < 3; ++k) {
                coordres[i][l][k] = coord[i][label-1][k];
            }
            l += 1;
        }
    }
    return coordres;
}

vector<vector<double>>
ReduceCoordSize(vector<vector<vector<double>>>& coord) {
    //
    int frames = coord.size();
    int natm   = coord[0].size();
    int dim    = natm * 3;

    vector<vector<double>> q(frames, vector<double>(dim, 0.0));
    for (int i = 0; i < frames; ++i) {
        int l = 0;
        for (int j = 0; j < natm; ++j) {
            for (int k = 0; k < 3; ++k) {
                q[i][l] = coord[i][j][k];
                l += 1;
            }
        }
    }

    return q;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// In order rmsd ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
double CalcRMSD(const vector<vector<double>>& coord,
               const vector<vector<double>>& refcoord) {
    int natm = coord.size();
    double rmsd = 0.0;

    for (int i = 0; i < natm; ++i) {
        double dx = coord[i][0] - refcoord[i][0];
        double dy = coord[i][1] - refcoord[i][1];
        double dz = coord[i][2] - refcoord[i][2];
        rmsd += dx*dx + dy*dy + dz*dz;
    }
    rmsd = sqrt(rmsd / natm);
    //
    return rmsd;
}

vector<vector<double>>
CalcRMSDall(const s_fit& fit,
            const double& dt,
            const vector<vector<vector<double>>>& coord) {
    //
    int frames = coord.size();
    vector<vector<double>> rmsd(2, vector<double>(frames, 0.0));
    //
    for (int i = 0; i < frames; ++i) {
        rmsd[0][i] = (i+1) * dt;
        for (int j = 0; j < fit.natm; ++j) {
            for (int k = 0; k < 3; ++k) {
                rmsd[1][i] += (coord[i][j][k] - fit.refcoord[j][k])
                              * (coord[i][j][k] - fit.refcoord[j][k]);
            }
        }
        rmsd[1][i] = sqrt(rmsd[1][i] / fit.natm);
    }
    return rmsd;
}

vector<vector<double>>
CalcRMSDres(const s_fit& fit,
            const vector<int>& reslist,
            const double& dt,
            const vector<vector<vector<double>>>& coord) {
    //
    int frames = coord.size();
    int nres   = reslist.size();

    vector<vector<double>> rmsd(2, vector<double>(frames, 0.0));
    //
    for (int i = 0; i < frames; ++i) {
        rmsd[0][i] = (i+1) * dt;
        for (int j = 0; j < nres; ++j) {
            int label = reslist[j];
            //cout << label << endl;
            for (int k = 0; k < 3; ++k) {
                rmsd[1][i] += (coord[i][j][k] - fit.refcoord[label-1][k])
                              * (coord[i][j][k] - fit.refcoord[label-1][k]);
            }
        }
        rmsd[1][i] = sqrt(rmsd[1][i] / fit.natm);
    }
    return rmsd;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// In order to RMSF /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
vector<double> CalcRMSF(const s_fit& fit,
                        const vector<vector<vector<double>>>& coord) {
    int frames = coord.size();
    vector<double> rmsf(fit.natm, 0.0);
    vector<vector<double>> avecoord(fit.natm, vector<double>(3, 0.0));

    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < fit.natm; ++j) {
            for (int k = 0; k < 3; ++k) {
                avecoord[j][k] += coord[i][j][k];
            }
        }
    }

    for (int j = 0; j < fit.natm; ++j) {
        for (int k = 0; k < 3; ++k) {
            avecoord[j][k] /= frames;
        }
    }

    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < fit.natm; ++j) {
            for (int k = 0; k < 3; ++k) {
                double diff = coord[i][j][k] - avecoord[j][k];
                rmsf[j] += diff * diff;
            }
        }
    }

    for (int i = 0; i < fit.natm; ++i) {
        rmsf[i] = sqrt(rmsf[i] / frames);
    }

    return rmsf;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////// In order to SDF /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
vector<vector<vector<double>>>
CalcSDF(const vector<string>& dcdfile,
        const vector<vector<int>>& psfinfo_fit,
        const vector<vector<int>>& psfinfo_sdf,
        const s_psf& psf,
        s_fit& fit,
        const vector<int>& ng3,
        const vector<double>& del3,
        const vector<double>& origin) {
    //
    int    natm    = psf.resname.size();
    int    nfile   = dcdfile.size();
    double vave    = 0.0;
    double rho     = 0.0;
    int    totfrms = 0;
    int    nres    = 0;

    cout << endl;
    cout << "Calcurate SDF..." << endl;
    vector<vector<vector<double>>>
    sdf(ng3[0], vector<vector<double>>(ng3[1], vector<double>(ng3[2], 0.0)));
    //
    for (int i = 0; i < nfile; ++i) {
        cout << endl;
        cout << "Read "
                  << dcdfile[i]
                  << ", and fitting..." << endl;
        s_dcd dcd;
        GetFitAllAtm(dcdfile[i], natm, psfinfo_fit, fit, dcd);
        totfrms += dcd.box.size();
        //
        vector<vector<vector<double>>>
        CoMRes = CalcCoMres(psfinfo_sdf, psf, dcd.coord);
        //
        int frames = CoMRes.size();
        nres   = CoMRes[0].size();
        for (int frame = 0; frame < frames; ++frame) {
            vave += dcd.box[frame][0] * dcd.box[frame][1] * dcd.box[frame][2];
            for (int j = 0; j < nres; ++j) {
                int binx = static_cast<int>((CoMRes[frame][j][0] - origin[0]) / del3[0]);
                int biny = static_cast<int>((CoMRes[frame][j][1] - origin[1]) / del3[1]);
                int binz = static_cast<int>((CoMRes[frame][j][2] - origin[2]) / del3[2]);
                //
                if (((binx >= 0) && (binx < ng3[0]))
                    && ((biny >= 0) && (biny < ng3[1]))
                    && ((binz >= 0) && (binz < ng3[2]))) {
                    sdf[binx][biny][binz] += 1.0;
                }
            }
        }
    }
    //
    vave /= totfrms;
    rho   = nres / vave;
    for (int i = 0; i < ng3[0]; ++i) {
        for (int j = 0; j < ng3[1]; ++j) {
            for (int k = 0; k < ng3[2]; ++k) {
                sdf[i][j][k] = sdf[i][j][k] / (rho * totfrms * del3[0] * del3[1] * del3[2]);
            }
        }
    }
    //
    return sdf;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// In order to PCA //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
s_pca PerformPCA(vector<vector<vector<double>>>& coord, 
                 const int& ncomp) {
    //
    s_pca pca;
    int frames = coord.size();
    int natm   = coord[0].size();
    int dim    = 3 * natm;
    //
    pca.VarCovMat.resize(dim, dim);
    pca.VarCovMat.setZero();
    //
    // Calc. average coordinate...
    vector<vector<double>> avecoord(natm, vector<double>(3, 0.0));

    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < natm; ++j) {
            for (int k = 0; k < 3; ++k) {
                avecoord[j][k] += coord[i][j][k];
            }
        }
    }
    //
    for (int j = 0; j < natm; ++j) {
        for (int k = 0; k < 3; ++k) {
            avecoord[j][k] /= frames;
        }
    }
    //
    //　summarize vector
    vector<double> qave(dim, 0.0);
    int l = 0;
    for (int i = 0; i < natm; ++i) {
        for (int j = 0; j < 3; ++j) {
            qave[l] = avecoord[i][j];
            l += 1;
        }
    }

    vector<vector<double>> q = ReduceCoordSize(coord);

    //
    // Calc. variance-convariance matrix
    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
                double delq_j = q[i][j] - qave[j];
                double delq_k = q[i][k] - qave[k];
                pca.VarCovMat(j, k) += delq_j * delq_k;
            }
        }
    }
    //
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            pca.VarCovMat(i, j) /= frames;
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(pca.VarCovMat);
    Eigen::VectorXd eigenvalues  = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();
    //
    Eigen::VectorXd sortedEigenvalues = eigenvalues.reverse();
    Eigen::MatrixXd sortedEigenvectors = eigenvectors.rowwise().reverse();
    //
    //s_pca pca;
    for (int i = 0; i < ncomp; ++i) {
        Eigen::VectorXd pc = sortedEigenvectors.col(i);
        pca.pcx.push_back(pc);
    }
    //
    pca.zeta.resize(ncomp, vector<double>(frames, 0.0));
    for (int i = 0; i < frames; ++i) {
        for (int j = 0; j < dim; ++j) {
            double delq = q[i][j] - qave[j];
            for (int k = 0; k < ncomp; ++k) {
                pca.zeta[k][i] += pca.pcx[k](j) * delq;
            }
        } 
    }
    
    return pca;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// In order to UMAP /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//vector<vector<float>>
//PerformUMAP(const vector<vector<vector<double>>>& coord,
//            const int& n_neighbors,
//            const int& out_dim,
//            const double& min_dist) {
//    //
//    vector<vector<double>>
//    q = ReduceCoordSize(vector<vector<vector<double>>>& coord);
//    //
//    vector<vector<float>> embedding(q.size(), vector<float>(out_dim));
//    umap.run(q, embedding, out_dim, n_neighbors, min_dist);
//
//    return embedding;
//}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// In order to dPCA /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
s_pca PerformdPCA(vector<vector<vector<double>>>& coord, 
                       const int& ncomp) {
    //
    s_pca pca;
    int frames = coord.size();
    int natm   = coord[0].size();
    int dim    = natm;
    //
    pca.VarCovMat.resize(dim, dim);
    pca.VarCovMat.setZero();

    vector<vector<double>> dave(dim, vector<double>(dim, 0));
    for (int f = 0; f < frames; ++f) {
        for (int i = 0; i < dim; ++ i) {
            for (int j = 0; j < dim; ++ j) {
                double d_ij = 0.0;
                for (int k = 0; k < 3; ++k) {
                    d_ij += pow(coord[f][i][k] - coord[f][j][k], 2);
                }
                d_ij = sqrt(d_ij);
                dave[i][j] += d_ij;
            }
        }
    }
    //
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            dave[i][j] /= frames;
        }
    }
    //
    // Calc. variance-convariance matrix
    for (int f = 0; f < frames; ++f) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double d_ij = 0.0;
                for (int k = 0; k < 3; ++k) {
                    d_ij += pow(coord[f][i][k] - coord[f][j][k], 2);
                }
                d_ij = sqrt(d_ij);
                pca.VarCovMat(i, j) += pow(d_ij - dave[i][j], 2);
            }
        }
    }
    
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            pca.VarCovMat(i, j) /= frames;
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(pca.VarCovMat);
    Eigen::VectorXd eigenvalues  = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();
    //
    Eigen::VectorXd sortedEigenvalues = eigenvalues.reverse();
    Eigen::MatrixXd sortedEigenvectors = eigenvectors.rowwise().reverse();
    //
    //s_pca pca;
    for (int i = 0; i < ncomp; ++i) {
        Eigen::VectorXd pc = sortedEigenvectors.col(i);
        pca.pcx.push_back(pc);
    }
    //
    pca.zeta.resize(ncomp, vector<double>(frames, 0.0));
    for (int f = 0; f < frames; ++f) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double deld = pow(coord[f][i][0] - coord[f][j][0], 2)
                              + pow(coord[f][i][1] - coord[f][j][1], 2)
                              + pow(coord[f][i][2] - coord[f][j][2], 2);
                deld = sqrt(deld);
                deld -= dave[i][j];
                for (int k = 0; k < ncomp; ++k) {
                    pca.zeta[k][i] += pca.pcx[k](j) * deld;
                }
            } 
        }
    }
    
    return pca;
}
///////////////////////////////////////////////////////////////////////////////
