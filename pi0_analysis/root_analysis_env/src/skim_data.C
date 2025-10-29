// ======================================================================
// skim_data.C
// ======================================================================
// Purpose  : Skim large HMS+NPS ROOT files to smaller ones with only 
//            selected branches enabled.
// Usage    : root -l -b -q 'src/skim_data.C+("data/","config/runlist.txt","output/skimmed/")'
// ======================================================================

#include "utils.C"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

// ----------------------------------------------------------------------
// Helper: Add files matching wildcard manually (for portability)
// ----------------------------------------------------------------------
void AddFilesWithWildcard(TChain* chain, const TString& pattern) {
    TString path = gSystem->DirName(pattern);
    TString base = gSystem->BaseName(pattern);

    void* dirp = gSystem->OpenDirectory(path);
    if(!dirp){
        logmsg(ERROR, Form("Cannot open directory: %s", path.Data()));
        return;
    }

    const char* fname;
    int added = 0;

    while((fname = gSystem->GetDirEntry(dirp))) {
        TString file = fname;
        if(file.Contains(base(0, base.First('*')))) { // match prefix before *
            TString fullpath = Form("%s/%s", path.Data(), file.Data());
            chain->Add(fullpath);
            ++added;
        }
    }

    gSystem->FreeDirectory(dirp);
    logmsg(INFO, Form("Added %d files matching %s", added, pattern.Data()));
}



// ----------------------------------------------------------------------
// Main skim function
// ----------------------------------------------------------------------
void skim_data(
               const TString &inputDir="/mnt/d/avnish/nps_data_rootfiles_pass2_updated/LD2",
            //    const TString &inputDir="/mnt/c/Users/as8oc/OneDrive/Documents/data_files_ifarm/LD2",
               const TString &runlistPath="config/runlist_x60_4b.txt",
               const TString &outDir="output/skimmed/") {

    TStopwatch sw; sw.Start();
    logmsg(INFO, "========================================================");
    logmsg(INFO, "Starting skim_data");
    logmsg(INFO, "InputDir  : " + string(inputDir.Data()));
    logmsg(INFO, "Runlist   : " + string(runlistPath.Data()));
    logmsg(INFO, "OutputDir : " + string(outDir.Data()));
    logmsg(INFO, "========================================================");

    gSystem->mkdir(outDir, true);

    TString indir = inputDir;
    if(!indir.EndsWith("/")) indir += "/";
    TString outdir = outDir;
    if(!outdir.EndsWith("/")) outdir += "/";

    // --- Read desired branches
    vector<string> branches = readBranchList("/home/ubuntu/nps_analysis/pi0_analysis/root_analysis_env/config/branches_to_read.txt");
    if(branches.empty()) logmsg(WARN, "Branch list empty. Will include all branches.");

    // --- Read runlist
    ifstream rl(runlistPath.Data());
    if(!rl.is_open()) {
        logmsg(ERROR, "Cannot open runlist file: " + string(runlistPath.Data()));
        return;
    }

    string line;
    while(getline(rl, line)) {
        line = trim(line);
        if(line.empty() || line[0]=='#') continue;
        TString token(line.c_str());

        vector<TString> filesToChain;

        // --- Handle absolute path or run:segment syntax
        if(token.BeginsWith("/")) {
            filesToChain.push_back(token);
        } 
        else if(token.Contains(":")) {
            Ssiz_t colon = token.Index(":");
            TString run = token(0, colon);
            TString segs = token(colon+1, token.Length()-colon-1);
            Ssiz_t dash = segs.Index("-");
            if(dash == kNPOS) {
                // single segment
                filesToChain.push_back(Form("%snps_hms_coin_%s_%s_1_-1.root", indir.Data(), run.Data(), segs.Data()));
            } else {
                int a = atoi(segs(0,dash).Data());
                int b = atoi(segs(dash+1, segs.Length()-dash-1).Data());
                for(int s=a; s<=b; ++s)
                    filesToChain.push_back(Form("%snps_hms_coin_%s_%d_1_-1.root", indir.Data(), run.Data(), s));
            }
        } 
        // --- Just run number
        else {
            filesToChain.push_back(Form("%snps_hms_coin_%s_*.root", indir.Data(), token.Data()));
        }


        // --- Build TChain
        TChain chain("T"); // Replace "T" with actual tree name (e.g. "T", "evtTree")
        for(auto &f : filesToChain) {
            if(f.Contains("*")) AddFilesWithWildcard(&chain, f);
            else chain.Add(f);
        }

        if(chain.GetNtrees() == 0) {
            logmsg(WARN, "No files found for token " + string(token.Data()));
            continue;
        }

        logmsg(INFO, TString::Format("Processing %lld entries across %d trees for run %s", 
               chain.GetEntries(), chain.GetNtrees(), token.Data()).Data());

        // --- Apply branch filtering
        if(!branches.empty()) enableBranches(&chain, branches);

        // --- Create output file
        TString outname = Form("%sskim_run%s.root", outdir.Data(), token.Data());
        TFile fout(outname, "RECREATE");
        TTree *tout = chain.CloneTree(0);

        // --- Loop and fill
        Long64_t nEntries = chain.GetEntries();
        for(Long64_t i=0; i<nEntries; ++i) {
            chain.GetEntry(i);
            tout->Fill();
            if(i>0 && i%100000==0)
                logmsg(INFO, TString::Format("Run %s: processed %lld / %lld", token.Data(), i, nEntries).Data());
        }

        tout->Write();
        fout.Close();
        logmsg(INFO, "Skim complete for run " + string(token.Data()) + " -> " + string(outname.Data()));
    }

    rl.close();
    sw.Stop();
    logmsg(INFO, TString::Format("skim_data completed in %.2f s", sw.RealTime()).Data());
}
