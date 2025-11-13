#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TH1.h"

// Small logging helper
enum LogLevel { INFO, WARN, ERROR };
void logmsg(LogLevel level, const std::string &msg);

// Read lines from simple key=value config file
std::map<std::string,std::string> readKeyValueConfig(const std::string &filename);

// Read branch list file (ignores comments '#' and blank lines)
std::vector<std::string> readBranchList(const std::string &filename);

// Set branch status on TTree/TChain: disable all and enable those in vector
void enableBranches(TTree* tree, const std::vector<std::string>& branches);

// Build a TChain for runs with a filename pattern (e.g. nps_{run}_{segment}.root)
// pattern should include %d for run and optionally %d for segment if needed.
// Simpler usage: pattern = "/data/nps_%d_*.root"
TChain* buildChain(const std::string& treeName, const std::vector<int>& runs, const std::string& filePattern);

// Utility to make consistent filenames
std::string makeOutputFilename(const std::string& outdir, const std::string& tag, int run, const std::string& suffix);

// Simple string trim helpers
std::string trim(const std::string &s);

// ------------------------------------------------------------------
// Read run list (supports comments beginning '#')
// ------------------------------------------------------------------
vector<int> readRunList(const string &fname);

#endif
