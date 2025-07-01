// compilation: g++ skim_pi0_data.cpp $(root-config --cflags --libs) -o skim_pi0_data
// usage single segment: ./skim_pi0_data 1234 1 false
// usage all segment: ./skim_pi0_data 1234


#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <filesystem>

void skim_pi0_data(int run, int maxSegment = 100, bool mergeAllSegments = true) {
    TString inputDir = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated/";
    TString outputDir = "/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/pi0_skimmed/";
    TString outputFileName = Form("nps_hms_coin_skimmed_%d.root", run);
    TString fullOutputPath = outputDir + outputFileName;

    TChain chain("T");

    int filesAdded = 0;
    for (int seg = 0; seg < maxSegment; ++seg) {
        TString fileName = Form("nps_hms_coin_%d_%d_1_-1.root", run, seg);
        TString fullPath = inputDir + fileName;

        if (!std::filesystem::exists(fullPath.Data())) break;  // No more segments
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
    chain.SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    chain.SetBranchStatus("H.gtr.dp", 1);
    chain.SetBranchStatus("H.cal.etotnorm", 1);
    chain.SetBranchStatus("H.cer.npeSum", 1);
    chain.SetBranchStatus("H.gtr.th", 1);
    chain.SetBranchStatus("H.gtr.ph", 1);
    chain.SetBranchStatus("H.gtr.y", 1);
    chain.SetBranchStatus("NPS.cal.nclust", 1);
    chain.SetBranchStatus("NPS.cal.clusE", 1);
    chain.SetBranchStatus("NPS.cal.clusX", 1);
    chain.SetBranchStatus("NPS.cal.clusY", 1);
    chain.SetBranchStatus("NPS.cal.clusT", 1);

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

    skim_pi0_data(run, maxSeg, mergeAll);
    return 0;
}

