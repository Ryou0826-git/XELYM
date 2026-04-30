// def_func.hpp
#ifndef DEF_FUNC_HPP
#define DEF_FUNC_HPP

#include <exprtk.hpp>
#include <iostream>

using namespace std;

exprtk::expression<double> CreateExpression(const string& expression_str, 
                                            double& x);

#endif
