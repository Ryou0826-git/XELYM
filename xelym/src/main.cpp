// main.cpp

#include <input.hpp>
#include <analyze.hpp>

using namespace std;

int main(int argc, char *args[]) {
    //
    if (argc != 2) {
        PrintExample();
        return 1;
    }
    //
    string inputfile = args[1];
    Config config = load_config(inputfile);
    config.print();

    if (boost::iequals(config.mode, "connect_only")) {
        //
        MAKE_POLYMER((*config.pdbfile)[0], *config.outpdb, (*config.chargefile)[0], *config.outcharge, 
                     *config.joint_dist, *config.joint_ref, *config.joint_mov, *config.npol, 
                     *config.captype, *config.netcharge);
    } else if (boost::iequals(config.mode, "connect")) {
        //
        MAKE_POLYMER((*config.pdbfile)[0], *config.outpdb, (*config.chargefile)[0], *config.outcharge, 
                     *config.joint_dist, *config.joint_ref, *config.joint_mov, *config.npol, 
                     *config.captype, *config.netcharge);
    } else if (boost::iequals(config.mode, "combine-pol-cross")) {
        //
        COMBINE_POL_CROSS(*config.itpfile, *config.outitp, *config.Np, *config.Nc);
    } else if (boost::iequals(config.mode, "systop")) {
        //
        MAKE_POLYMER((*config.pdbfile)[0], *config.outpdb, (*config.chargefile)[0], *config.outcharge, 
                     *config.joint_dist, *config.joint_ref, *config.joint_mov, *config.npol, 
                     *config.captype, *config.netcharge);
    } else if (boost::iequals(config.mode, "cross")) {
        //
        MAKE_CROSSLINKING(*config.itpfile, *config.trajfile, *config.outitp, 
                          *config.chargefile, *config.selpol, *config.selcross, 
                          *config.joint_string, *config.Np,
                          *config.Nc, *config.rc, *config.rule);
    } else if (boost::iequals(config.mode, "convert")) {
        //
        PolymerConvertItp((*config.itpfile)[0], *config.outitp, *config.resid,
                           *config.atomlist, *config.joint_ref,
                           *config.joint_mov, *config.npol);
        if (config.pdbfile) {
            PolymerConvertPdb((*config.pdbfile)[0], *config.outpdb, *config.resid,
                              *config.atomlist, *config.joint_ref,
                              *config.joint_mov, *config.npol);
        }
    } else if (boost::iequals(config.mode, "remove")) {
        //
        REMOVE_ANALYSIS((*config.itpfile)[0], *config.outitp, (*config.pdbfile)[0], 
                         *config.outpdb, *config.selpol, *config.selcross, *config.Nc);
    } else if (boost::iequals(config.mode, "remove-resid")) {
        //
        REMOVE_RESID_ANALYSIS((*config.itpfile)[0], *config.outitp, (*config.pdbfile)[0], 
                               *config.outpdb, *config.sel1);
    } else if (boost::iequals(config.mode, "check")) {
        //
        CROSS_LINKING_CHECK((*config.itpfile)[0], *config.selpol, *config.selcross, *config.Nc);
    } else if (boost::iequals(config.mode, "bond-and-remove")) {
        //
        if (config.ncycle && config.prob) {
            ANALYZE_BOND_AND_REMOVE(*config.pdbfile, *config.outpdb, 
                                    *config.itpfile, *config.outitp, (*config.chargefile)[0],
                                    *config.joint_string, 
                                    *config.joint_string_polymer, 
                                    *config.distance, config.ncycle, config.prob);
        } else {
            ANALYZE_BOND_AND_REMOVE(*config.pdbfile, *config.outpdb, 
                                 *config.itpfile, *config.outitp, (*config.chargefile)[0],
                                 *config.joint_string, 
                                 *config.joint_string_polymer, 
                                 *config.distance);
        }
    } else if (boost::iequals(config.mode, "bond-and-remove-ac")) {
        //
        if (config.ncycle) {
            ANALYZE_BOND_AND_REMOVE_AC(*config.pdbfile, *config.outpdb, 
                                       *config.itpfile, *config.outitp, 
                                       (*config.chargefile)[0], 
                                       *config.joint_string, 
                                       *config.joint_string_polymer,
                                       *config.joint_string_cross, 
                                       *config.joint_string_remove, 
                                       *config.nbond, 
                                       *config.distance, 
                                       config.ncycle);
        } else {
            ANALYZE_BOND_AND_REMOVE_AC(*config.pdbfile, *config.outpdb, 
                                       *config.itpfile, *config.outitp, 
                                       (*config.chargefile)[0], 
                                       *config.joint_string, 
                                       *config.joint_string_polymer,
                                       *config.joint_string_cross, 
                                       *config.joint_string_remove, 
                                       *config.nbond, 
                                       *config.distance);
        }
    } else if (boost::iequals(config.mode, "crosslinking-dry")) {
        MAKE_CROSSLINKING_DRY(*config.itpfile, *config.trajfile, *config.outitp, 
                              (*config.pdbfile)[0], *config.outpdb, *config.chargefile, 
                              *config.selpol, *config.selcross, 
                              *config.drypol, *config.drycross,
                              *config.joint_string, *config.Np,
                              *config.Nc, *config.rc);
    } else if (boost::iequals(config.mode, "crosslinking-radical")) {
        MAKE_CROSSLINKING_RADICAL(*config.itpfile, *config.trajfile, *config.outitp, 
                                  *config.chargefile, *config.selpol, *config.selcross, 
                                  *config.Np, *config.Nc, *config.rc);
    } else if (boost::iequals(config.mode, "remove-ions")) {
        ANALYZE_REMOVE_IONS((*config.itpfile)[0], (*config.pdbfile)[0], 
                             *config.outitp, *config.outpdb, *config.ions_string);
    } else if (boost::iequals(config.mode, "cl-position")) {
        //
        ANALYZE_CL_POSITION((*config.itpfile)[0], *config.trajlist, 
                            *config.selpol, *config.selcross, *config.Nc, 
                            *config.Min, *config.Max, *config.grid);
    } else if (boost::iequals(config.mode, "hopping")) {
        ANALYZE_HOPPING(*config.topfile, *config.trajlist, 
                         *config.selion, *config.selmemb, 
                         *config.dt, *config.rc, *config.rst_sign);
    } else if (boost::iequals(config.mode, "sd")) {
        ANALYZE_SD(*config.topfile, *config.trajfile, *config.sel, *config.Next, *config.dt);
    } else if (boost::iequals(config.mode, "coord_num")) {
        ANALYZE_COORD_NUM(*config.topfile, *config.trajlist, *config.sel1, *config.sel2, *config.rc);
    } else if (boost::iequals(config.mode, "end-to-end")) {
        ANALYZE_END_TO_END(*config.topfile, *config.trajlist, *config.selpol, *config.Np, *config.npol, *config.nbin);
    } else if (boost::iequals(config.mode, "mass")) {
        ANALYZE_MASS((*config.itpfile)[0]);
    } else if (boost::iequals(config.mode, "change-charge")) {
        CHANGE_CHARGE((*config.itpfile)[0], *config.outitp, *config.selmemb, *config.aft_charge);
    } else {
        cout << config.mode << " is not supported, try again..." << endl;
    } 

    return 0;
}
