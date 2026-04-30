// select_tools.cpp

#include <select_tools.hpp>

auto trim = [](const string& s) {
    return boost::algorithm::trim_copy(s);
};

vector<int> SELOPR::ResearchPdb(const s_pdb& pdb, 
                                const s_sel& sel) {
    vector<int> sellist;
    //

    int natm = pdb.resid.size();

    for (int i = 0; i < natm; ++i) {
        bool atm_match = !sel.atmnm ||
                         boost::algorithm::iequals(
                             trim(*sel.atmnm),
                             trim(pdb.atmnm[i])
                         );

        bool res_match = !sel.resnm ||
                         boost::algorithm::iequals(
                             trim(*sel.resnm),
                             trim(pdb.resnm[i])
                         );

        bool id_match  = !sel.resid || *sel.resid == pdb.resid[i];

        if (id_match && atm_match && res_match) {
            sellist.push_back(i+1);
        }
    }
    return sellist;
}

s_sel SELOPR::ParseSel(const string& expr) {
    istringstream in(expr);
    string tok;
    s_sel sel;

    auto lower = [](string& s) {
        transform(s.begin(), s.end(), s.begin(),
                  [](unsigned char c){ return tolower(c); });
    };

    while (in >> tok) {
        lower(tok);

        if (tok == "resid") {
            int v;  in >> v;
            sel.resid = v;
        } else if (tok == "atmnm") {
            string v;  in >> v;
            sel.atmnm = v;
        } else if (tok == "resnm") {
            string v;  in >> v;
            sel.resnm = v;
        } else if (tok == "and") {
            continue;
        } else {
            throw runtime_error("Unknown token: " + tok);
        }
    }
    return sel;
}

string trim_paren(const string &s) {
    string res = s;
    size_t start = res.find_first_not_of(" \t");
    size_t end   = res.find_last_not_of(" \t");
    if (start == string::npos) return "";
    res = res.substr(start, end - start + 1);

    if (!res.empty() && res.front() == '(' && res.back() == ')') {
        res = res.substr(1, res.size() - 2);
        start = res.find_first_not_of(" \t");
        end   = res.find_last_not_of(" \t");
        if (start != string::npos)
            res = res.substr(start, end - start + 1);
        else
            res.clear();
    }
    return res;
}

vector<s_sel> ExtractSel(const string& sel) {
    SELOPR SelOpr;
    vector<s_sel> sel_list;
    //
    vector<string> groups;

    string token;
    istringstream iss(sel);
    while (getline(iss, token, 'o')) {
    }
    
    size_t pos = 0, found;
    while ((found = sel.find(" or ", pos)) != string::npos) {
        groups.push_back(sel.substr(pos, found - pos));
        pos = found + 4;
    }
    groups.push_back(sel.substr(pos));

    for (string &g : groups) {
        sel_list.push_back(SelOpr.ParseSel(trim_paren(g)));
    }
    //
    return sel_list;
}
