// traj_tools.hpp
#ifndef TRAJ_TOOLS_HPP
#define TRAJ_TOOLS_HPP

#include "dcd_r.hpp"
#include "psf_tools.hpp"

using namespace std;

void allocateMemory(vector<vector<vector<double>>>& coord, 
                    int natms, 
                    int frames);

void AllocateDCD(s_dcd& dcd, 
                 int natms, 
                 int frames);

void readAndStoreCoordinates(DCD_R& dcdf, 
                             vector<vector<vector<double>>>& coord, 
                             int nAtoms, 
                             int frames);

void readAndStoreCoordBox(DCD_R& dcdf, 
                          s_dcd& dcd, 
                          int nAtoms, 
                          int frames);

void CatCoordinate(const s_dcd& dcd, 
                   const vector<vector<double>>& ligcoord, 
                   const int& frame, 
                   const int& insf, 
                   s_dcd& dcdins);

s_dcd ReaddcdFile(const string& dcdfile, 
                  const s_psf& psf);

s_dcd ReaddcdFileAmber(const string& dcdfile, 
                       const int& natm);

#endif
