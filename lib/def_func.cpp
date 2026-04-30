// def_func.cpp

#include <def_func.hpp>

//using namespace std;

exprtk::expression<double> CreateExpression(const string& expression_str, 
                                            double& x) {
    //
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>   expression_t;
    typedef exprtk::parser<double>       parser_t;
    //
    symbol_table_t symbol_table;
    symbol_table.add_variable("x", x);

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    if (!parser.compile(expression_str, expression)) {
        throw runtime_error("Failed: " + expression_str);
    }

    return expression;
}

