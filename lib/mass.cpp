// mass.cpp

#include "mass.hpp"
#include <ostream>
#include <iostream>

double ReturnMass(const char& atm) {
    double mass;
    //
    if (atm == 'C') {
        mass = 12.0110;
    } else if (atm == 'H') {
        mass = 1.0080;
    } else if (atm == 'N') {
        mass = 14.0070;
    } else if (atm == 'O') {
        mass = 15.9994;
    } else if (atm == 'S') {
        mass = 32.0600;
    } else {
        std::cout << "Invalid atm: " << atm << std::endl;
        mass = 0.0f;
    }
    return mass;
}
