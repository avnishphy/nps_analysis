// // compilation tcsh: g++ skim_pi0_data.cpp `root-config --cflags --libs` -o skim_pi0_data    
// // usage single segment: ./skim_pi0_data 1234 1 false
// // usage all segment: ./skim_pi0_data 1234


// #include <TFile.h>
// #include <TTree.h>
// #include <TChain.h>
// #include <TString.h>

// #include <iostream>
// #include <vector>
// #include <filesystem>

// void skim_pi0_data(int run, int maxSegment = 100, bool mergeAllSegments = true) {
//     TString inputDir = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated/";
//     TString outputDir = "/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/pi0_skimmed/";
//     TString outputFileName = Form("nps_hms_coin_skimmed_%d.root", run);
//     TString fullOutputPath = outputDir + outputFileName;

//     TChain chain("T");
//     TChain chain1("TSH");

//     int filesAdded = 0;
//     for (int seg = 0; seg < maxSegment; ++seg) {
//         TString fileName = Form("nps_hms_coin_%d_%d_1_-1.root", run, seg);
//         TString fullPath = inputDir + fileName;

//         if (!std::filesystem::exists(fullPath.Data())) break;  // No more segments
//         chain.Add(fullPath);
//         ++filesAdded;

//         if (!mergeAllSegments) break;  // Only add first segment if requested
//     }

//     if (filesAdded == 0) {
//         std::cerr << "No input files found for run " << run << std::endl;
//         return;
//     }

//     // Disable all branches by default
//     chain.SetBranchStatus("*", 0);

//     // Enable only needed branches
//     chain.SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
//     chain.SetBranchStatus("T.hms.hTRIG1_tdcTimeRaw", 1);
//     chain.SetBranchStatus("T.hms.hTRIG2_tdcTimeRaw", 1);
//     chain.SetBranchStatus("T.hms.hTRIG3_tdcTimeRaw", 1);
//     chain.SetBranchStatus("T.hms.hTRIG4_tdcTimeRaw", 1);
//     chain.SetBranchStatus("T.hms.hTRIG5_tdcTimeRaw", 1);
//     chain.SetBranchStatus("T.hms.hTRIG6_tdcTimeRaw", 1);
//     chain.SetBranchStatus("H.cal.etot", 1);
//     chain.SetBranchStatus("H.cal.etotnorm", 1);
//     chain.SetBranchStatus("H.cal.etottracknorm", 1);
//     chain.SetBranchStatus("H.cal.etrack", 1);
//     chain.SetBranchStatus("H.cal.etracknorm", 1);
//     chain.SetBranchStatus("H.cer.npeSum", 1);
//     chain.SetBranchStatus("H.gtr.dp", 1);
//     chain.SetBranchStatus("H.gtr.th", 1);
//     chain.SetBranchStatus("H.gtr.ph", 1);
//     chain.SetBranchStatus("H.gtr.p", 1);
//     chain.SetBranchStatus("H.gtr.px", 1);
//     chain.SetBranchStatus("H.gtr.py", 1);
//     chain.SetBranchStatus("H.gtr.pz", 1);
//     chain.SetBranchStatus("H.gtr.x", 1);
//     chain.SetBranchStatus("H.gtr.y", 1);
//     chain.SetBranchStatus("NPS.cal.nclust", 1);
//     chain.SetBranchStatus("NPS.cal.clusE", 1);
//     chain.SetBranchStatus("NPS.cal.clusX", 1);
//     chain.SetBranchStatus("NPS.cal.clusY", 1);
//     chain.SetBranchStatus("NPS.cal.clusT", 1);

//     // Disable all branches by default
//     chain1.SetBranchStatus("*", 0);

//     // Enable only needed branches
//     chain1.SetBranchStatus("H.BCM4A.scalerCharge", 1)
//     chain1.SetBranchStatus("H.BCM4A.scalerCurrent", 1)
//     chain1.SetBranchStatus("H.BCM4A.scaler", 1)
//     chain1.SetBranchStatus("H.EDTM.scaler", 1)
//     chain1.SetBranchStatus("H.hTRIG4.scaler", 1)
//     chain1.SetBranchStatus("H.hL1ACCP.scaler", 1)
//     chain1.SetBranchStatus("H.hEL_REAL.scaler", 1)
//     chain1.SetBranchStatus("H.1MHz.scalerTime", 1)

//     chain1.SetBranchAddress("H.BCM4A.scalerCharge", H_BCM4A_scalerCharge)
//     chain1.SetBranchAddress("H.BCM4A.scalerCurrent", H_BCM4A_scalerCurrent)
//     chain1.SetBranchAddress("H.BCM4A.scaler", H_BCM4A_scaler)

//     chain1.SetBranchAddress("H.BCM4A.scalerCurrent", H_BCM4A_scalerCurrent)
//     chain1.SetBranchAddress("H.EDTM.scaler", H_EDTM_scaler)
//     chain1.SetBranchAddress("H.hTRIG4.scaler", H_hTRIG4_scaler)
//     chain1.SetBranchAddress("H.hL1ACCP.scaler", H_hL1ACCP_scaler)
//     chain1.SetBranchAddress("H.hEL_REAL.scaler", H_hEL_REAL_scaler)
//     chain1.SetBranchAddress("H.1MHz.scalerTime", H_1MHz_scalerTime)
   

//     // Create skimmed output file
//     TFile *fout = TFile::Open(fullOutputPath, "RECREATE");
//     if (!fout || fout->IsZombie()) {
//         std::cerr << "Error creating output file: " << fullOutputPath << std::endl;
//         return;
//     }

//     // Clone skimmed tree
//     TTree *skimmedTree = chain.CloneTree(0); // Empty tree with same structure
//     TTree *skimmedTree1 = chain1.CloneTree(0); // Empty tree with same structure

//     Long64_t nentries = chain.GetEntries();
//     std::cout << "Processing T tree " << nentries << " entries..." << std::endl;

//     for (Long64_t i = 0; i < nentries; ++i) {
//         chain.GetEntry(i);
//         skimmedTree->Fill();
//     }

//     Long64_t nentries1 = chain1.GetEntries();
//     std::cout << "Processing TSH tree " << nentries1 << " entries..." << std::endl;

//     for (Long64_t i = 0; i < nentries1; ++i) {
//         chain1.GetEntry(i);
//         skimmedTree1->Fill();
//     }

//     skimmedTree->Write();
//     skimmedTree1->Write();

//     fout->Close();

//     std::cout << "Skimmed file written to: " << fullOutputPath << std::endl;
// }

// int main(int argc, char* argv[]) {
//     if (argc < 2) {
//         std::cerr << "Usage: " << argv[0] << " <run_number> [max_segments] [mergeAllSegments=1]" << std::endl;
//         return 1;
//     }

//     int run = std::stoi(argv[1]);
//     int maxSeg = (argc >= 3) ? std::stoi(argv[2]) : 100;
//     bool mergeAll = (argc >= 4) ? std::stoi(argv[3]) != 0 : true;

//     skim_pi0_data(run, maxSeg, mergeAll);
//     return 0;
// }


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

        if (!std::filesystem::exists(fullPath.Data())) break;  // No more segments
        chain.Add(fullPath);
        chain1.Add(fullPath); // <<-- add to chain1 as well
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

    // Create skimmed output file
    TFile *fout = TFile::Open(fullOutputPath, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error creating output file: " << fullOutputPath << std::endl;
        return;
    }

    // Clone skimmed tree (empty)
    TTree *skimmedTree = chain.CloneTree(0);   // empty tree with selected branches
    TTree *skimmedTree1 = chain1.CloneTree(0); // empty tree with selected branches

    Long64_t nentries = chain.GetEntries();
    std::cout << "Processing T tree " << nentries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nentries; ++i) {
        chain.GetEntry(i);
        skimmedTree->Fill();
    }

    Long64_t nentries1 = chain1.GetEntries();
    std::cout << "Processing TSH tree " << nentries1 << " entries..." << std::endl;
    for (Long64_t i = 0; i < nentries1; ++i) {
        chain1.GetEntry(i);
        skimmedTree1->Fill();
    }

    skimmedTree->Write();
    skimmedTree1->Write();

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
