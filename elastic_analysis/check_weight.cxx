#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

void plot_weights() {
    // Open the ROOT file
    TFile *file = TFile::Open("/u/group/nps/singhav/simc_gfortran/worksim/eep_hydrogen_6828_6841.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the tree
    TTree *tree = (TTree*)file->Get("h10");
    if (!tree) {
        std::cerr << "Error: Tree 'h10' not found!" << std::endl;
        file->Close();
        return;
    }

    // Check branches
    tree->Print();

    // Define histograms
    TH1D *hWeight = new TH1D("hWeight", "Weight Distribution;Weight;Counts", 50, 0, 10);
    TH1D *hFullWeight = new TH1D("hFullWeight", "Full Weight Distribution;Full Weight;Counts", 50, 0, 10);

    // Define variables and set branch
    double Weight;
    tree->SetBranchAddress("Weight", &Weight);

    // Compute normalization factor
    double normfac = 0.355463E+07;
    Long64_t nevents = tree->GetEntries();
    std::cout << "Number of events: " << nevents << std::endl;
    
    if (nevents == 0) {
        std::cerr << "No events in the tree!" << std::endl;
        file->Close();
        return;
    }

    double weight_factor = normfac / nevents;

    // Loop over tree entries and fill histograms
    for (Long64_t i = 0; i < nevents; i++) {
        tree->GetEntry(i);
        double full_weight = Weight * weight_factor;
        hWeight->Fill(Weight);
        hFullWeight->Fill(full_weight);
    }

    // Create canvas and draw histograms
    TCanvas *c1 = new TCanvas("c1", "Weight Histograms", 800, 600);
    hWeight->SetLineColor(kBlue);
    hWeight->Draw();
    hFullWeight->SetLineColor(kRed);
    hFullWeight->Draw("SAME");

    // Add legend
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hWeight, "Weight", "l");
    legend->AddEntry(hFullWeight, "Full Weight (Scaled)", "l");
    legend->Draw();

    c1->Draw();
    c1->Update();
    
    // Close file
    file->Close();
}
