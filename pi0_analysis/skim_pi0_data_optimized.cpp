// compilation tcsh: g++ -std=c++17 skim_pi0_data.cpp `root-config --cflags --libs` -o skim_pi0_data
// usage single segment: ./skim_pi0_data 1234 1 0
// usage all segments: ./skim_pi0_data 1234

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <filesystem>
#include <cstdlib> // for std::stoi

void skim_pi0_data(int run, int maxSegment = 100, bool mergeAllSegments = true) {
    TString inputDir = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated/";
    TString outputDir = "/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/pi0_skimmed/";
    TString outputFileName = Form("nps_hms_coin_skimmed_%d.root", run);
    TString fullOutputPath = outputDir + outputFileName;

    TChain chain("T");
    TChain chain1("TSH");

    int filesAdded = 0;
    for (int seg = 0; seg < maxSegment; ++seg) {
        TString fileName = Form("nps_hms_coin_%d_%d_1_-1.root", run, seg);
        TString fullPath = inputDir + fileName;

        if (!std::filesystem::exists(std::string(fullPath.Data()))) break;  // No more segments
        chain.Add(fullPath);
        chain1.Add(fullPath); // add to chain1 as well
        ++filesAdded;

        if (!mergeAllSegments) break;  // Only add first segment if requested
    }

    if (filesAdded == 0) {
        std::cerr << "No input files found for run " << run << std::endl;
        return;
    }

    std::cout << "Added " << filesAdded << " file(s) to chains" << std::endl;

    // Disable all branches by default
    chain.SetBranchStatus("*", 0);

    // Enable only needed branches
    chain.SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    chain.SetBranchStatus("T.hms.hTRIG1_tdcTimeRaw", 1);
    chain.SetBranchStatus("T.hms.hTRIG2_tdcTimeRaw", 1);
    chain.SetBranchStatus("T.hms.hTRIG3_tdcTimeRaw", 1);
    chain.SetBranchStatus("T.hms.hTRIG4_tdcTimeRaw", 1);
    chain.SetBranchStatus("T.hms.hTRIG5_tdcTimeRaw", 1);
    chain.SetBranchStatus("T.hms.hTRIG6_tdcTimeRaw", 1);
    chain.SetBranchStatus("H.cal.etot", 1);
    chain.SetBranchStatus("H.cal.etotnorm", 1);
    chain.SetBranchStatus("H.cal.etottracknorm", 1);
    chain.SetBranchStatus("H.cal.etrack", 1);
    chain.SetBranchStatus("H.cal.etracknorm", 1);
    chain.SetBranchStatus("H.cer.npeSum", 1);
    chain.SetBranchStatus("H.gtr.dp", 1);
    chain.SetBranchStatus("H.gtr.th", 1);
    chain.SetBranchStatus("H.gtr.ph", 1);
    chain.SetBranchStatus("H.gtr.p", 1);
    chain.SetBranchStatus("H.gtr.px", 1);
    chain.SetBranchStatus("H.gtr.py", 1);
    chain.SetBranchStatus("H.gtr.pz", 1);
    chain.SetBranchStatus("H.gtr.x", 1);
    chain.SetBranchStatus("H.gtr.y", 1);
    chain.SetBranchStatus("NPS.cal.nclust", 1);
    chain.SetBranchStatus("NPS.cal.clusE", 1);
    chain.SetBranchStatus("NPS.cal.clusX", 1);
    chain.SetBranchStatus("NPS.cal.clusY", 1);
    chain.SetBranchStatus("NPS.cal.clusT", 1);

    // Disable all branches by default for chain1
    chain1.SetBranchStatus("*", 0);

    // Enable only needed branches for chain1
    chain1.SetBranchStatus("H.BCM4A.scalerCharge", 1);
    chain1.SetBranchStatus("H.BCM4A.scalerCurrent", 1);
    chain1.SetBranchStatus("H.BCM4A.scaler", 1);
    chain1.SetBranchStatus("H.BCM2.scalerCharge", 1);
    chain1.SetBranchStatus("H.BCM2.scalerCurrent", 1);
    chain1.SetBranchStatus("H.BCM2.scaler", 1);
    chain1.SetBranchStatus("H.EDTM.scaler", 1);
    chain1.SetBranchStatus("H.hTRIG4.scaler", 1);
    chain1.SetBranchStatus("H.hL1ACCP.scaler", 1);
    chain1.SetBranchStatus("H.hEL_REAL.scaler", 1);
    chain1.SetBranchStatus("H.1MHz.scalerTime", 1);

    // Create skimmed output file with optimizations
    TFile *fout = TFile::Open(fullOutputPath, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error creating output file: " << fullOutputPath << std::endl;
        return;
    }
    fout->SetCacheSize(10 * 1024 * 1024);  // 10 MB cache for faster I/O

    // Clone skimmed tree (empty)
    TTree *skimmedTree = chain.CloneTree(0);   // empty tree with selected branches
    if (!skimmedTree) {
        std::cerr << "Error cloning T tree" << std::endl;
        fout->Close();
        delete fout;
        return;
    }
    
    TTree *skimmedTree1 = chain1.CloneTree(0); // empty tree with selected branches
    if (!skimmedTree1) {
        std::cerr << "Error cloning TSH tree" << std::endl;
        delete skimmedTree;
        fout->Close();
        delete fout;
        return;
    }
    
    // Disable AutoSave to speed up writing
    skimmedTree->SetAutoSave(0);
    skimmedTree1->SetAutoSave(0);

    Long64_t nentries = chain.GetEntries();
    std::cout << "Processing T tree " << nentries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nentries; ++i) {
        if (chain.LoadTree(i) < 0) break;  // Load tree segment
        chain.GetEntry(i);
        skimmedTree->Fill();
        if ((i + 1) % 100000 == 0) {
            std::cout << "  Processed " << (i + 1) << " entries" << std::endl;
        }
    }

    Long64_t nentries1 = chain1.GetEntries();
    std::cout << "Processing TSH tree " << nentries1 << " entries..." << std::endl;
    for (Long64_t i = 0; i < nentries1; ++i) {
        if (chain1.LoadTree(i) < 0) break;  // Load tree segment
        chain1.GetEntry(i);
        skimmedTree1->Fill();
        if ((i + 1) % 100000 == 0) {
            std::cout << "  Processed " << (i + 1) << " entries" << std::endl;
        }
    }

    std::cout << "Writing trees to disk..." << std::endl;
    skimmedTree->Write(nullptr, TObject::kOverwrite);
    skimmedTree1->Write(nullptr, TObject::kOverwrite);

    // Proper resource cleanup
    delete skimmedTree;
    delete skimmedTree1;
    fout->Close();
    delete fout;

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
