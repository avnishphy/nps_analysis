// compilation tcsh: g++ skim_elastic_data.cpp `root-config --cflags --libs` -o skim_elastic_data    
// usage single segment: ./skim_pi0_data 1234 1 false
// usage all segment: ./skim_pi0_data 1234


#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <filesystem>

void skim_elastic_data(int run, int maxSegment = 18, bool mergeAllSegments = true) {
    TString inputDir = "/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/hms_elastics/";
    TString outputDir = "/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/hms_elastics_skimmed/";
    TString outputFileName = Form("nps_hms_coin_skimmed_%d.root", run);
    TString fullOutputPath = outputDir + outputFileName;

    TChain chain("T");

    int filesAdded = 0;
    for (int seg = 0; seg < maxSegment; ++seg) {
        TString fileName = Form("nps_hms_elastics_%d_%d_1_-1.root", run, seg);
        TString fullPath = inputDir + fileName;

        if (!std::filesystem::exists(fullPath.Data())) continue;  // No more segments
        chain.Add(fullPath);
        ++filesAdded;

        if (!mergeAllSegments) break;  // Only add first segment if requested
    }

    if (filesAdded == 0) {
        std::cerr << "No input files found for run " << run << std::endl;
        return;
    }

    // Disable all branches by default
    chain.SetBranchStatus("*", 0);

    // Enable only needed branches
    chain.SetBranchStatus("H.kin.W", 1);
    chain.SetBranchStatus("H.kin.x_bj", 1);
    chain.SetBranchStatus("H.cal.etottracknorm", 1);
    chain.SetBranchStatus("H.hod.goodscinhit", 1);
    chain.SetBranchStatus("H.cer.npeSum", 1);
    chain.SetBranchStatus("H.gtr.dp", 1);
    chain.SetBranchStatus("H.gtr.th", 1);
    chain.SetBranchStatus("H.gtr.ph", 1);
    chain.SetBranchStatus("H.gtr.y", 1);
    chain.SetBranchStatus("H.react.z", 1);
    chain.SetBranchStatus("H.kin.scat_ang_deg", 1);
    chain.SetBranchStatus("H.kin.Q2", 1);
    chain.SetBranchStatus("H.dc.x_fp", 1);
    chain.SetBranchStatus("H.dc.y_fp", 1);
    chain.SetBranchStatus("H.dc.xp_fp", 1);
    chain.SetBranchStatus("H.dc.yp_fp", 1);

    // Create skimmed output file
    TFile *fout = TFile::Open(fullOutputPath, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error creating output file: " << fullOutputPath << std::endl;
        return;
    }

    // Clone skimmed tree
    TTree *skimmedTree = chain.CloneTree(0); // Empty tree with same structure

    Long64_t nentries = chain.GetEntries();
    std::cout << "Processing " << nentries << " entries..." << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        chain.GetEntry(i);
        skimmedTree->Fill();
    }

    skimmedTree->Write();
    fout->Close();

    std::cout << "Skimmed file written to: " << fullOutputPath << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run_number> [max_segments] [mergeAllSegments=1]" << std::endl;
        return 1;
    }

    int run = std::stoi(argv[1]);
    int maxSeg = (argc >= 3) ? std::stoi(argv[2]) : 100;
    bool mergeAll = (argc >= 4) ? std::stoi(argv[3]) != 0 : true;

    skim_elastic_data(run, maxSeg, mergeAll);
    return 0;
}
