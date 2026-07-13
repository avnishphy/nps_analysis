#include "utils.h"
#include "TSystem.h"
#include "TString.h"
#include <TSystemDirectory.h>


void logmsg(LogLevel level, const std::string &msg) {
    const char* levelstr = (level==INFO?"INFO":(level==WARN?"WARN":"ERROR"));
    std::cout << "[" << levelstr << "] " << msg << std::endl;
}

std::map<std::string,std::string> readKeyValueConfig(const std::string &filename) {
    std::map<std::string,std::string> kv;
    std::ifstream in(filename.c_str());
    if(!in.is_open()) {
        logmsg(ERROR, "Cannot open config file: " + filename);
        return kv;
    }
    std::string line;
    while(std::getline(in, line)) {
        // remove comments
        auto pos = line.find('#');
        if(pos != std::string::npos) line = line.substr(0,pos);
        // trim
        auto start = line.find_first_not_of(" \t\r\n");
        if(start==std::string::npos) continue;
        auto end = line.find_last_not_of(" \t\r\n");
        std::string t = line.substr(start, end-start+1);
        auto eq = t.find('=');
        if(eq == std::string::npos) continue;
        std::string key = trim(t.substr(0, eq));
        std::string val = trim(t.substr(eq+1));
        kv[key] = val;
    }
    in.close();
    return kv;
}

std::vector<std::string> readBranchList(const std::string &filename) {
    std::vector<std::string> branches;
    std::ifstream in(filename.c_str());
    if(!in.is_open()) {
        logmsg(WARN, "Branch list not found: " + filename + ". Will require manual enabling.");
        return branches;
    }
    std::string line;
    while(std::getline(in, line)) {
        // strip comments
        auto pos = line.find('#');
        if(pos != std::string::npos) line = line.substr(0,pos);
        std::string s = trim(line);
        if(s.empty()) continue;
        branches.push_back(s);
    }
    in.close();
    return branches;
}

void enableBranches(TTree* tree, const std::vector<std::string>& branches) {
    if(!tree) {
        logmsg(ERROR, "enableBranches: null tree.");
        return;
    }
    tree->SetBranchStatus("*", 0); // disable all
    for(const auto &b : branches) {
        tree->SetBranchStatus(b.c_str(), 1);
    }
}

TChain* buildChain(const std::string& treeName, const std::vector<int>& runs, const std::string& filePattern) {
    TChain* chain = new TChain(treeName.c_str());
    for(auto run : runs) {
        // Replace %d with run
        char fname[1024];
        snprintf(fname, sizeof(fname), filePattern.c_str(), run);
        // Use TSystem to expand a wildcarded pattern
        TSystemDirectory dir(".", ".");
        // We try to add name directly - TChain supports glob on some systems; to be safe, attempt Add with pattern
        int added = chain->Add(fname);
        if(added==0) {
            std::string msg = Form("No files matched pattern %s (run %d)", fname, run);
            logmsg(WARN, msg);
        } else {
            std::string msg = Form("Added %d files for pattern %s (run %d)", added, fname, run);
            logmsg(INFO, msg);
        }
    }
    return chain;
}

std::string makeOutputFilename(const std::string& outdir, const std::string& tag, int run, const std::string& suffix) {
    std::ostringstream ss;
    ss << outdir;
    if(outdir.back()!='/') ss << "/";
    ss << tag << "_run" << run << suffix;
    return ss.str();
}

std::string trim(const std::string &s) {
    if(s.empty()) return s;
    auto a = s.find_first_not_of(" \t\r\n"); if(a==std::string::npos) return "";
    auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b-a+1);
}

// ------------------------------------------------------------------
// Read run list (supports comments beginning '#')
// ------------------------------------------------------------------
vector<int> readRunList(const string &fname) {
    vector<int> runs;
    ifstream in(fname);
    if (!in.is_open()) {
        logmsg(WARN, "Cannot open runlist: " + fname);
        return runs;
    }
    string line;
    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        try { runs.push_back(stoi(line)); }
        catch (...) { logmsg(WARN, "Skipping invalid runlist line: " + line); }
    }
    return runs;
}
