// extension.cpp

#include <extension.hpp>

string get_extension(const string& filename) {
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos == string::npos || dot_pos == filename.length() - 1) {
        return "";
    }
    return filename.substr(dot_pos + 1);
}
