// output.hpp
#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <make_hist.hpp>

using namespace std;

void CoordToFile(const string& inputfile, 
                 const vector<vector<vector<double>>>& data, 
                 int selid);

void data3dToFile(const string& filename, 
                  const vector<vector<double>>& data);

void data2dToFile(const string& filename, 
                  const vector<vector<double>>& data);

void data2dToFileLines(const string& filename, 
                       const vector<vector<double>>& data);

void AxisToFile(const string& filename, 
                const vector<double>& data);

void Only1dToFile(const string& filename,
                  const vector<double>& data);

void Mapdata2dToFile(const string& filename, 
                     const vector<vector<double>>& data);

void SplitdataToFile(const string& filename, 
                     const vector<double>& datax, 
                     const vector<double>& datay);

void SplitdataToFileINTdouble(const string& filename, 
                              const vector<double>& datax, 
                              const vector<double>& datay);

void SplitdataToFileNUMdouble(const string& filename, 
                              const vector<double>& data);

void SplitINTdataToFileNUMdouble(const string& filename, 
                                 const vector<int>& data);

void OutdxFile(const vector<int>& ng3, 
               const vector<double>& del3, 
               const vector<double>& origin, 
               const vector<vector<vector<double>>>& data, 
               const string& filename);

void dataXdToFile(const string& filename,
                  const vector<vector<double>>& data);

void INTdataXdToFile(const string& filename,
                     const vector<int>& times,
                     const vector<vector<double>>& data);

void TimedataXdToFile(const string& filename,
                      const vector<double>& times,
                      const vector<vector<double>>& data);

void TimeINTdataXdToFile(const string& filename,
                         const vector<double>& times,
                         const vector<vector<int>>& data);

void EigenMapToFile(const string& filename,
                    const Eigen::MatrixXd data);

void write_3d_hist_gnuplot(const string& filename,
                           const vector<double>& hist,
                           const vector<double>& Min,
                           const vector<double>& Max,
                           int grid);

#endif
