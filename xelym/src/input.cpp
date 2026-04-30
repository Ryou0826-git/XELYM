//input.cpp

#include "input.hpp"

using json = nlohmann::json;
using namespace std;

template <typename T>
static void get_optional(const json& j,
                         const string& key,
                         optional<T>& value)
{
    if (j.contains(key) && !j.at(key).is_null()) {
        value = j.at(key).get<T>();
    }
}

void from_json(const json& j, Config& c)
{
    // required
    j.at("mode").get_to(c.mode);

    // optional<string / vector<string>>
    get_optional(j, "pdbfile", c.pdbfile);
    get_optional(j, "outpdb", c.outpdb);
    get_optional(j, "chargefile", c.chargefile);
    get_optional(j, "outcharge", c.outcharge);
    get_optional(j, "itpfile", c.itpfile);
    get_optional(j, "outitp", c.outitp);
    get_optional(j, "topfile", c.topfile);
    get_optional(j, "itpimitfile", c.itpimitfile);
    get_optional(j, "trajfile", c.trajfile);
    get_optional(j, "trajlist", c.trajlist);
    get_optional(j, "captype", c.captype);

    // optional<int>
    get_optional(j, "npol", c.npol);
    get_optional(j, "Np", c.Np);
    get_optional(j, "Nc", c.Nc);
    get_optional(j, "Next", c.Next);
    get_optional(j, "Ni", c.Ni);
    get_optional(j, "nbond", c.nbond);
    get_optional(j, "ncycle", c.ncycle);
    get_optional(j, "rst_sign", c.rst_sign);
    get_optional(j, "grid", c.grid);
    get_optional(j, "nbin", c.nbin);

    // optional<double>
    get_optional(j, "joint_dist", c.joint_dist);
    get_optional(j, "dt", c.dt);
    get_optional(j, "rc", c.rc);
    get_optional(j, "bond_c0", c.bond_c0);
    get_optional(j, "bond_c1", c.bond_c1);
    get_optional(j, "netcharge", c.netcharge);
    get_optional(j, "distance", c.distance);
    get_optional(j, "prob", c.prob);
    get_optional(j, "aft_charge", c.aft_charge);

    // optional<vector<int>>
    get_optional(j, "joint_ref", c.joint_ref);
    get_optional(j, "joint_mov", c.joint_mov);

    // optional<vector<double>>
    get_optional(j, "Min", c.Min);
    get_optional(j, "Max", c.Max);

    // optional<vector<string>>
    get_optional(j, "joint_string", c.joint_string);
    get_optional(j, "joint_string_polymer", c.joint_string_polymer);
    get_optional(j, "joint_string_cross", c.joint_string_cross);
    get_optional(j, "joint_string_remove", c.joint_string_remove);
    get_optional(j, "selpol", c.selpol);
    get_optional(j, "selcross", c.selcross);
    get_optional(j, "selion", c.selion);
    get_optional(j, "selmemb", c.selmemb);
    get_optional(j, "sel", c.sel);
    get_optional(j, "sel1", c.sel1);
    get_optional(j, "sel2", c.sel2);
    get_optional(j, "drypol", c.drypol);
    get_optional(j, "drycross", c.drycross);
    get_optional(j, "atomlist", c.atomlist);
    get_optional(j, "ions_string", c.ions_string);

    // optional<string>
    get_optional(j, "resid", c.resid);
    get_optional(j, "rule", c.rule);
}

Config load_config(const string& filename)
{
    ifstream ifs(filename);
    if (!ifs) {
        throw runtime_error("Cannot open file: " + filename);
    }

    json j;
    ifs >> j;

    Config cfg = j.get<Config>(); // from_jsonが勝手に呼び出される

    return cfg;
}

template <typename T>
static void print_optional(const string& name,
                           const optional<T>& v)
{
    if (v) cout << name << " : " << *v << endl;
}

template <typename T>
static void print_optional_vector(const string& name,
                                  const optional<vector<T>>& v)
{
    if (!v) return;

    cout << name << " : [ ";
    for (const auto& x : *v) cout << x << " ";
    cout << "]" << endl;
}

//
// ★ Config::print()
void Config::print() const
{
    cout << "==== Your input ====" << endl;

    cout << "mode : " << mode << endl;
    print_optional_vector("pdbfile", pdbfile);
    print_optional("outpdb", outpdb);
    print_optional_vector("chargefile", chargefile);
    print_optional("outcharge", outcharge);
    print_optional_vector("itpfile", itpfile);
    print_optional("outitp", outitp);
    print_optional("topfile", topfile);
    print_optional("itpimitfile", itpimitfile);
    print_optional("trajfile", trajfile);
    print_optional_vector("trajlist", trajlist);
    print_optional("npol", npol);
    print_optional("joint_dist", joint_dist);
    print_optional_vector("joint_ref", joint_ref);
    print_optional_vector("joint_mov", joint_mov);
    print_optional_vector("joint_string", joint_string);
    print_optional_vector("joint_string_polymer", joint_string_polymer);
    print_optional_vector("joint_string_cross", joint_string_cross);
    print_optional_vector("joint_string_remove", joint_string_remove);
    print_optional("Np", Np);
    print_optional("Nc", Nc);
    print_optional("Next", Next);
    print_optional("Ni", Ni);
    print_optional("dt", dt);
    print_optional("rc", rc);
    print_optional("bond_c0", bond_c0);
    print_optional("bond_c1", bond_c1);
    print_optional("netcharge", netcharge);
    print_optional_vector("selpol", selpol);
    print_optional_vector("selcross", selcross);
    print_optional_vector("selion", selion);
    print_optional_vector("selmemb", selmemb);
    print_optional_vector("sel", sel);
    print_optional_vector("sel1", sel1);
    print_optional_vector("sel2", sel2);
    print_optional_vector("drypol", drypol);
    print_optional_vector("drycross", drycross);
    print_optional_vector("atomlist", atomlist);
    print_optional("resid", resid);
    print_optional("rule", rule);
    print_optional("nbond", nbond);
    print_optional("distance", distance);
    print_optional("ncycle", ncycle);
    print_optional("prob", prob);
    print_optional_vector("ions_string", ions_string);
    print_optional("rst_sign", rst_sign);
    print_optional_vector("Min", Min);
    print_optional_vector("Max", Max);
    print_optional("grid", grid);
    print_optional("aft_charge", aft_charge);
    print_optional("nbin", nbin);
    
    cout << "====================" << endl;
}
