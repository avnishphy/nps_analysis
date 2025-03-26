#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

// Open the two ROOT files
TFile *file1 = TFile::Open("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_6834_0_1_-1.root");
TFile *file2 = TFile::Open("/u/group/nps/singhav/simc_gfortran/worksim/eep_hydrogen_6828_6841_HMS4042.root");

// Get the trees from the files
TTree *tree1 = (TTree*)file1->Get("T");
TTree *tree2 = (TTree*)file2->Get("h10");

// Define variables to hold the values of the branches
Double_t H_dc_xpfp;
Float_t fhsxpfp;

// Set branches for tree1
tree1->SetBranchAddress("H.dc.xp_fp", &H_dc_xpfp);

// Set branches for tree2
tree2->SetBranchAddress("hsxpfp", &fhsxpfp);

// Create a canvas to plot the branches
TCanvas *canvas = new TCanvas("canvas", "Canvas for Branches", 800, 600);

// Create a 2D histogram to store the values of the branches
TH2F *hist = new TH2F("hist", "H.dc.xp_fp vs hsxpfp", 100, -10, 10, 100, -10, 10); 

// Loop over the entries in both trees and fill the 2D histogram
Int_t nEntries1 = tree1->GetEntries();
Int_t nEntries2 = tree2->GetEntries();

// Optionally, ensure the loop runs over the smallest number of entries in both trees
Int_t nEntries = (nEntries1 < nEntries2) ? nEntries1 : nEntries2;

for (Int_t i = 0; i < nEntries; i++) {
    tree1->GetEntry(i); // Get entry from tree1
    tree2->GetEntry(i); // Get entry from tree2
    
    // Fill the histogram with the values from both branches
    hist->Fill(H_dc_xpfp, fhsxpfp);
}

// Set axis labels
hist->GetXaxis()->SetTitle("H.dc.xp_fp");
hist->GetYaxis()->SetTitle("hsxpfp");

// Draw the 2D histogram
hist->Draw("COLZ");  // Use COLZ for color-coded 2D plot
