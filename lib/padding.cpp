// padding.cpp

#include <padding.hpp>

string ZeroPad(const string& head,
               const string& extension,
               const int frame) {
    ostringstream forigin;
    forigin << head
            << setfill('0') << setw(4) << frame
            << "." << extension;
    return forigin.str();
}
