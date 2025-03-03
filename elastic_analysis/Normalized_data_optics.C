#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include "Variables_declaration.h" // For std::istringstream

void CollimatorStudy(){

 //Method to study various collimator cuts on the H(e,e'p) and D(e,e'p)n  Yield across Ytar, Y'tar, X'tar and delta

  cout << "Calling CollimatorStudy() . . . " << endl;
  
  //Scaling the HMS/SHMS Collimator Cuts
  hms_hsize = hms_scale*hms_hsize;  //The scale factor is read from set_heep_cuts.inp
  hms_vsize = hms_scale*hms_vsize;
  
  shms_hsize = shms_scale*shms_hsize;
  shms_vsize = shms_scale*shms_vsize;  

  //Define HMS Collimator Shape
  hms_Coll_gCut = new TCutG("hmsCollCut", 8 );
  hms_Coll_gCut->SetVarX("X");
  hms_Coll_gCut->SetVarY("Y");
 
  hms_Coll_gCut->SetPoint(0,  hms_hsize,     hms_vsize/2.);
  hms_Coll_gCut->SetPoint(1,  hms_hsize/2.,  hms_vsize   );
  hms_Coll_gCut->SetPoint(2, -hms_hsize/2.,  hms_vsize   );
  hms_Coll_gCut->SetPoint(3, -hms_hsize,     hms_vsize/2.);
  hms_Coll_gCut->SetPoint(4, -hms_hsize,    -hms_vsize/2.);
  hms_Coll_gCut->SetPoint(5, -hms_hsize/2., -hms_vsize   );
  hms_Coll_gCut->SetPoint(6,  hms_hsize/2., -hms_vsize   );
  hms_Coll_gCut->SetPoint(7,  hms_hsize,    -hms_vsize/2.);
  hms_Coll_gCut->SetPoint(8,  hms_hsize,     hms_vsize/2.);

  //Define SHMS Collimator Shape
  shms_Coll_gCut = new TCutG("shmsCollCut", 8 );
  shms_Coll_gCut->SetVarX("X");
  shms_Coll_gCut->SetVarY("Y");
 
  shms_Coll_gCut->SetPoint(0,  shms_hsize,     shms_vsize/2.);
  shms_Coll_gCut->SetPoint(1,  shms_hsize/2.,  shms_vsize   );
  shms_Coll_gCut->SetPoint(2, -shms_hsize/2.,  shms_vsize   );
  shms_Coll_gCut->SetPoint(3, -shms_hsize,     shms_vsize/2.);
  shms_Coll_gCut->SetPoint(4, -shms_hsize,    -shms_vsize/2.);
  shms_Coll_gCut->SetPoint(5, -shms_hsize/2., -shms_vsize   );
  shms_Coll_gCut->SetPoint(6,  shms_hsize/2., -shms_vsize   );
  shms_Coll_gCut->SetPoint(7,  shms_hsize,    -shms_vsize/2.);
  shms_Coll_gCut->SetPoint(8,  shms_hsize,     shms_vsize/2.);

  cout << "Ending CollimatorStudy() . . . " << endl;



}

 //  === initialize boolean for analysis cuts ===
     // Bool_t hmsColl_Cut = 0;
      Bool_t Em_cut = 0;
      Bool_t hdc_ntrk_cut ;
      Bool_t hScinGood_cut ;
      Bool_t hcer_NPE_Sum_cut ;
      Bool_t hetotnorm_cut ;
      Bool_t hBeta_notrk_cut ;

      Bool_t pdc_ntrk_cut ;
      Bool_t pScinGood_cut ;
      Bool_t phgcer_NPE_Sum_cut;
      Bool_t pngcer_NPE_Sum_cut ;
      Bool_t petotnorm_cut ;
      Bool_t pBeta_notrk_cut ;
     Bool_t phgcer_cut_flag = 0;
     
      


void analyze_tree(TString ana_type="") {

    // Declare common variables
    TFile *inFile;
   

    // Declare histograms
    TH1F *H_W, *H_Q2, *H_theta_e, *H_theta_p, *H_Pmx, *H_Pmy, *H_Pmz, *H_h_delta, *H_e_delta,  *H_scat_ang_rad_data, *H_h_theta_data, *H_xangle_data;
    TH1F *H_h_xptar, *H_e_xptar, *H_h_yptar, *H_e_yptar, *H_h_ytar, *H_e_ytar, *H_Em, *H_x_bj, *H_nu, *H_hXColl, *H_hYColl, *H_eXColl, *H_eYColl;
    TH1F *H_ztarDiff, *H_Pm, *H_W_sf, *H_Q2_sf, *H_Em_sf;
    const Double_t Mp = 0.938;
    const Double_t Eb = 10.542;





    //TH1F *H_W, *H_Q2, *H_Pmx, *H_Pmy, *H_Pmz, *H_e_delta, *H_h_delta, *H_h_xptar, *H_e_xptar, *H_h_yptar, *H_e_yptar, *H_Em;
    // Create a list to manage histograms
    TList *kin_HList = new TList();

    // Check analysis type and initialize accordingly
    if (ana_type == "simc") {
        std::cout << "Analyzing SIMC data..." << std::endl;

        // SIMC-specific input file
        TString simc_fname = "d2_heep_scan_rad_-8.root";
        inFile = new TFile(simc_fname.Data(), "READ");
        TTree *tree = (TTree*)inFile->Get("SNT");
        Long64_t nentries = tree->GetEntries();

        // Define SIMC-specific variables
        Double_t W, Q2, theta_e, theta_p, Pmx, Pmy, Pmz, h_delta, e_delta, hYColl, hXColl, pf;
        Double_t h_xptar, e_xptar, h_yptar, e_yptar, h_ytar, e_ytar, nu, h_xfp, h_yfp, h_ypfp, h_xpfp;
        Double_t Normfac, Weight, FullWeight, Em, x_bj, tar_x, htar_y, etar_y, etar_z, htar_z, Pm;
        Double_t Pcalc, Pmea, h_pf, deltaPsimc_frac, deltaPsimc_vs_yfp, deltaPsimc_vs_ypfp;
        Double_t  Em_vs_Pm, deltaPsimc_vs_xfp, deltaPsimc_vs_xpfp, deltaPsimc_frac_vs_yfp, deltaPsimc_frac_vs_ypfp, deltaPsimc_frac_vs_xfp, deltaPsimc_frac_vs_xpfp, Pmea_GeV;
        
        Double_t ztarDiff = htar_z - etar_z;
        //Double_t nu;
        // Set branch addresses
        
        tree->SetBranchAddress("W", &W);
        tree->SetBranchAddress("Q2", &Q2);
        tree->SetBranchAddress("theta_e", &theta_e);
        tree->SetBranchAddress("theta_p", &theta_p);
        tree->SetBranchAddress("Pm", &Pm);
        tree->SetBranchAddress("Pmx", &Pmx);
        tree->SetBranchAddress("Pmy", &Pmy);
        tree->SetBranchAddress("Pmz", &Pmz);
        tree->SetBranchAddress("h_delta", &h_delta);
        tree->SetBranchAddress("e_delta", &e_delta);
        tree->SetBranchAddress("h_xptar", &h_xptar);
        tree->SetBranchAddress("e_xptar", &e_xptar);
        tree->SetBranchAddress("h_ypfp", &h_ypfp);
        tree->SetBranchAddress("h_yfp", &h_yfp);
        tree->SetBranchAddress("h_xpfp", &h_xpfp);
        tree->SetBranchAddress("h_xfp", &h_xfp);
        tree->SetBranchAddress("h_yptar", &h_yptar);
        tree->SetBranchAddress("e_yptar", &e_yptar);
        tree->SetBranchAddress("h_ytar", &h_ytar);
        tree->SetBranchAddress("e_ytar", &e_ytar);
        tree->SetBranchAddress("Normfac", &Normfac);
        tree->SetBranchAddress("Weight", &Weight);
        tree->SetBranchAddress("Em", &Em);
        tree->SetBranchAddress("nu", &nu);
        tree->SetBranchAddress("h_pf", &Pmea);
          //(tarx, tary, tarz) in Hall Coord. System      
      tree->SetBranchAddress("tar_x", &tar_x);
      tree->SetBranchAddress("h_yv",  &htar_y);
      tree->SetBranchAddress("h_zv",  &htar_z);
      tree->SetBranchAddress("e_yv",  &etar_y);
      tree->SetBranchAddress("e_zv",  &etar_z);
     

        // Define histograms
        H_W = new TH1F("W", "W", 100, 0.5, 1.7);
        H_Q2 = new TH1F("Q2", "Q2", 100, 1.3, 6);
        H_theta_e = new TH1F("theta_e", "theta_e", 100, 9.15, 11.5);
        H_theta_p = new TH1F("theta_p", "theta_p", 100, 39, 45);
        H_Pm = new TH1F("Pm", "Pm", 100, -0.1, 0.1);
        H_Pmx = new TH1F("Pmx", "Pmx", 100, -0.1, 0.1);
        H_Pmy = new TH1F("Pmy", "Pmy", 100, -0.1, 0.1);
        H_Pmz = new TH1F("Pmz", "Pmz", 100, -0.1, 0.2);
        H_h_delta = new TH1F("h_delta", "h_delta", 100, -15, 15);
        H_e_delta = new TH1F("e_delta", "e_delta", 100, -12, 12);
        H_h_xptar = new TH1F("h_xptar", "h_xptar", 100, -0.2, 0.2);
        H_e_xptar = new TH1F("e_xptar", "e_xptar", 100, -0.03, 0.03);
        H_h_yptar = new TH1F("h_yptar", "h_yptar", 100, -0.06, 0.06);
        H_e_yptar = new TH1F("e_yptar", "e_yptar", 100, -0.02, 0.02);
        H_h_ytar = new TH1F("h_ytar", "h_ytar", 100, -15, 15);
        H_e_ytar = new TH1F("e_ytar", "e_ytar", 100, -2, 2);
        H_Em = new TH1F("Em","Em", 100, -0.1, 0.5);
        H_x_bj = new TH1F("x_bj", "x_bj", 100, 0.85, 1.05);
        H_nu = new TH1F("nu", "nu", 100, -1, 5);
        H_ztarDiff = new TH1F("ztarDiff", "ztar Difference(htar_z - etar_z); ztarDiff; Events", 100, -2.5,3);
        TH1F* H_Pmea = new TH1F("Pmea", "Pmea", 100, 1, 4);
        TH1F* H_Pcalc = new TH1F("Pcalc", "Pcalc", 100, 1, 4);
        TH1F* H_h_xfp = new TH1F("h_xfp", "h_xfp", 100, -40,40);
        TH1F* H_h_yfp = new TH1F("h_yfp", "h_yfp", 100, -25,25);
        TH1F* H_h_xpfp = new TH1F("h_xpfp", "h_xpfp", 100, -0.07,0.07);
        TH1F* H_h_ypfp = new TH1F("h_ypfp", "h_ypfp", 100, -0.02,0.02);
        TH2F* H_Em_vs_Pm = new TH2F("Em_vs_Pm", "Em_vs_Pm;  Pm;  Em", 100, -0.05, 0.05,  100, -0.01, 0.55);
        TH1F* H_deltaPsimc = new TH1F("deltaPsimc", "Pcalc-Pmea; deltaPsimc", 100, -0.1,0.1);
        TH1F* H_deltaPsimc_frac = new TH1F("deltaPsimc_frac", "(Pcalc-Pmea)/Pmea; deltaPsimc", 100, -0.1,0.1);
       
       TH2F* H_deltaPsimc_vs_xfp = new TH2F("deltaPsimc_frac_vs_xfp", "deltaPsimc_frac vs xfp; xfp; (Pcalc-Pmea)/Pmea", 100, -30, 30,  100, -0.05, 0.05);
       TH2F* H_deltaPsimc_vs_xpfp = new TH2F("deltaPsimc_frac_vs_xpfp", "deltaPsimc_frac vs xpfp; xpfp; (Pcalc-Pmea)/Pmea", 100, -0.04, 0.04,  100, -0.05, 0.05);
       TH2F* H_deltaPsimc_vs_yfp = new TH2F("deltaPsimc_frac_vs_yfp", "deltaPsimc_frac vs yfp; yfp; (Pcalc-Pmea)/Pmea", 100, -5, 17,  100, -0.05, 0.05);
       TH2F* H_deltaPsimc_vs_ypfp = new TH2F("deltaPsimc_frac_vs_ypfp", "deltaPsimc_frac vs ypfp; ypfp; (Pcalc-Pmea)/Pmea", 100, -0.01, 0.02,  100, -0.05, 0.05);
       
       // H_theta_p = new TH1F("theta_p","theta_p", 100, 0.06, 0.23);
        // Add histograms to the list
        kin_HList->Add(H_W);
        kin_HList->Add(H_Q2);
        kin_HList->Add(H_theta_e);
        kin_HList->Add(H_theta_p);
        kin_HList->Add(H_Pmx);
        kin_HList->Add(H_Pmy);
        kin_HList->Add(H_Pmz);
        kin_HList->Add(H_h_xfp);
        kin_HList->Add(H_h_xpfp);
        kin_HList->Add(H_h_yfp);
        kin_HList->Add(H_h_ypfp);
        kin_HList->Add(H_h_delta);
        kin_HList->Add(H_e_delta);
        kin_HList->Add(H_h_xptar);
        kin_HList->Add(H_e_xptar);
        kin_HList->Add(H_h_yptar);
        kin_HList->Add(H_e_yptar);
        kin_HList->Add(H_h_ytar);
        kin_HList->Add(H_e_ytar);
        kin_HList->Add(H_Em);
        kin_HList->Add(H_Pm);
        kin_HList->Add(H_x_bj);
        kin_HList->Add(H_nu);
        kin_HList->Add(H_ztarDiff);
        kin_HList->Add(H_Pmea);
        kin_HList->Add(H_Pcalc);
        kin_HList->Add(H_Em_vs_Pm);
        kin_HList->Add(H_deltaPsimc);
        kin_HList->Add(H_deltaPsimc_frac);
        kin_HList->Add(H_deltaPsimc_vs_xpfp);
        kin_HList->Add(H_deltaPsimc_vs_xfp);
        kin_HList->Add(H_deltaPsimc_vs_yfp);
        kin_HList->Add(H_deltaPsimc_vs_ypfp);
        //kin_HList->Add(H_theta_e);
        //kin_HList->Add(H_theta_p);
        // Event loop for SIMC

         // define cuts
         CollimatorStudy();
     
      

        for (Long64_t ientry = 0; ientry < nentries; ientry++) {
            tree->GetEntry(ientry);
            FullWeight = Normfac * Weight / nentries;
            Em_cut = Em>=-0.05 && Em<=0.05;
            hmsColl_Cut = hms_Coll_gCut->IsInside(hYColl, hXColl);
 
            // Convert radian angles to degrees
            Double_t theta_e_deg = theta_e * TMath::RadToDeg();
            Double_t theta_p_deg = theta_p * TMath::RadToDeg();
            Double_t x_bj = (Q2 / (2*Mp*nu));
            Double_t Pmea_GeV = Pmea/1000; //Changed into GeV/c
            Double_t ztarDiff = htar_z - etar_z;
            Double_t Pcalc = (2 * Mp * Eb * (Mp + Eb) * TMath::Cos(theta_p)) / 
                 (Mp * Mp + 2 * Mp * Eb + Eb * Eb * TMath::Sin(theta_p) * TMath::Sin(theta_p));
            Double_t deltaPsimc = (Pcalc - Pmea_GeV);
            Double_t deltaPsimc_frac = (Pcalc - Pmea_GeV)/Pmea_GeV;
            
         if(hmsColl_Cut && Em_cut && abs(h_delta)<= 10 &&e_delta>= -10 && e_delta <= 22 )
		    {
           H_W->Fill(W, FullWeight);
            H_Q2->Fill(Q2, FullWeight);
            H_theta_e->Fill(theta_e_deg, FullWeight);
            H_theta_p->Fill(theta_p_deg, FullWeight);
            H_Pm->Fill(Pm, FullWeight);
            H_Pmx->Fill(Pmx, FullWeight);
            H_Pmy->Fill(Pmy, FullWeight);
            H_Pmz->Fill(Pmz, FullWeight);
            H_h_xfp->Fill(h_xfp, FullWeight);
            H_h_xpfp->Fill(h_xpfp, FullWeight);
            H_h_yfp->Fill(h_yfp, FullWeight);
            H_h_ypfp->Fill(h_ypfp, FullWeight);
            H_h_delta->Fill(h_delta, FullWeight);
            H_e_delta->Fill(e_delta, FullWeight);
            H_h_xptar->Fill(h_xptar, FullWeight);
            H_e_xptar->Fill(e_xptar, FullWeight);
            H_h_yptar->Fill(h_yptar, FullWeight);
            H_e_yptar->Fill(e_yptar, FullWeight);
            H_h_ytar->Fill(h_ytar, FullWeight);
            H_e_ytar->Fill(e_ytar, FullWeight);
            H_Em->Fill(Em, FullWeight);
            H_Pmea->Fill(Pmea_GeV, FullWeight);
            H_Pcalc->Fill(Pcalc, FullWeight);
            H_nu->Fill(nu, FullWeight);
            H_x_bj->Fill(x_bj, FullWeight);
            H_ztarDiff->Fill(ztarDiff, FullWeight);
            H_Em_vs_Pm->Fill(Em, Pm, FullWeight);
            H_deltaPsimc->Fill(deltaPsimc, FullWeight);
            H_deltaPsimc_frac->Fill(deltaPsimc_frac, FullWeight);
           
            H_deltaPsimc_vs_xfp->Fill( h_xfp, deltaPsimc_frac);
            H_deltaPsimc_vs_xpfp->Fill(h_xpfp, deltaPsimc_frac);
            H_deltaPsimc_vs_yfp->Fill(h_yfp, deltaPsimc_frac);
            H_deltaPsimc_vs_ypfp->Fill(h_ypfp, deltaPsimc_frac);
           // H_theta_e->Fill(theta_p, FullWeight);
            //H_theta_p->Fill(theta_e, FullWeight);
        }
       
            std::cout << "SIMC Events Completed: " << std::setprecision(2)
                      << double(ientry) / nentries * 100. << " % " << std::flush << "\r";
        }
 // Write Histograms
       // H_Q2->Write();
       // H_W->Write();
        TString simc_OutputFileName = "d2_heep_-8_output_normalized.root";
        TFile *outROOT = new TFile(simc_OutputFileName, "RECREATE");
        outROOT->cd();
        kin_HList->Write();
        outROOT->Close();

    }
    
    //=====================
    //  DATA ANALYSIS 
    //======================
    else if (ana_type == "data") 
    {
      
      std::cout << "Analyzing DATA..." << std::endl;

      // DATA-specific input file
      TString data_fname = "deut_replay_prod_20846_-1.root";
      inFile = new TFile(data_fname.Data(), "READ");
        
      // DEFINE SCALER TREE
      //TTree *scaler_tree = (TTree*)inFile->Get("TSH");
      TTree *scaler_tree = (TTree*)inFile->Get("TSH");
      Long64_t scal_entries = scaler_tree->GetEntries();
   
      // DEFINE DATA TREE
      TTree *tree = (TTree*)inFile->Get("T");
      Long64_t nentries = tree->GetEntries();


      //=================
      // SCALER ANALYSIS
      //=================
    
      // Dynamically allocate arrays for event flags and event numbers
      Int_t *evt_flag_bcm = new Int_t[scal_entries]; // store 0 or 1, to determine which scaler read passed cut
      Int_t *scal_evt_num = new Int_t[scal_entries]; // store event associated with scaler read
    
    // Define DATA -specific variables
        Double_t W, Q2, e_delta, h_delta, Pm, Pmx, Pmy, Pmz, h_xptar, h_yptar, e_xptar, e_yptar;
        Double_t Em , x_bj, scat_ang_rad_data, xangle_data, theta_e,theta_e_rad,theta_p, gevnum, gevtyp; 
        Double_t etar_x, etar_y, etar_z, htar_x, htar_y, htar_z;
        Double_t eYColl, hXColl, hYColl, eXColl, ztarDiff, h_xfp, h_yfp, h_ypfp, h_xpfp;
        Double_t Pcalc, Pmea,nu;
        Double_t deltaPdata_frac, deltaPdata_frac_vs_xfp,deltaPdata_frac_vs_xpfp,deltaPdata_frac_vs_yfp,deltaPdata_frac_vs_ypfp;
  
  // ---- Set branch addresses to fetch data from the source tree -----
    scaler_tree->SetBranchAddress("evNumber", &evNum);
    //scaler_tree->SetBranchAddress("g.evnum", &gevnum);
    //scaler_tree->SetBranchAddress("g.evtyp", &gevtyp);
    scaler_tree->SetBranchAddress("H.BCM1.scalerCharge", &Scal_BCM1_charge);
    scaler_tree->SetBranchAddress("H.BCM1.scalerCurrent", &Scal_BCM1_current);
    scaler_tree->SetBranchAddress("H.BCM2.scalerCharge", &Scal_BCM2_charge);
    scaler_tree->SetBranchAddress("H.BCM2.scalerCurrent", &Scal_BCM2_current);
    scaler_tree->SetBranchAddress("H.BCM4A.scalerCharge", &Scal_BCM4A_charge);
    scaler_tree->SetBranchAddress("H.BCM4A.scalerCurrent", &Scal_BCM4A_current);
    scaler_tree->SetBranchAddress("H.BCM4B.scalerCharge", &Scal_BCM4B_charge);
    scaler_tree->SetBranchAddress("H.BCM4B.scalerCurrent", &Scal_BCM4B_current);
    scaler_tree->SetBranchAddress("H.BCM4C.scalerCharge", &Scal_BCM4C_charge);
    scaler_tree->SetBranchAddress("H.BCM4C.scalerCurrent", &Scal_BCM4C_current);
    scaler_tree->SetBranchAddress("H.1MHz.scalerTime", &Scal_time);
    scaler_tree->SetBranchAddress("H.S1X.scaler", &S1X_scaler);
    scaler_tree->SetBranchAddress("H.S2X.scaler",&S2X_scaler);  
    scaler_tree->SetBranchAddress("H.S2Y.scaler",&S2Y_scaler);  
    scaler_tree->SetBranchAddress("H.pTRIG1.scaler",&TRIG1_scaler);
    scaler_tree->SetBranchAddress("H.pTRIG2.scaler",&TRIG2_scaler);
    scaler_tree->SetBranchAddress("H.pTRIG3.scaler",&TRIG3_scaler);
    scaler_tree->SetBranchAddress("H.pTRIG4.scaler",&TRIG4_scaler);
    scaler_tree->SetBranchAddress("H.pTRIG5.scaler",&TRIG5_scaler);
    scaler_tree->SetBranchAddress("H.pTRIG6.scaler",&TRIG6_scaler);
    scaler_tree->SetBranchAddress("H.EDTM.scaler",  &EDTM_scaler);
   

    TString bcm_type="BCM4A";

    //-------------------
    // SCALER EVENT LOOP
    //--------------------

    for (int i = 0; i < scal_entries; i++) {
        scaler_tree->GetEntry(i);
        
        evt_flag_bcm[i] = 0;  // Default flag
        
        // Store the event number
        scal_evt_num[i] = evNum;
        
        // Your logic to determine if the flag should be set, etc.
         // Determine which bcm current to cut on (based on user input)
      if(bcm_type=="BCM1"){
	Scal_BCM_current = Scal_BCM1_current;
      }
      else if(bcm_type=="BCM2"){
	Scal_BCM_current = Scal_BCM2_current;
      }
      else if(bcm_type=="BCM4A"){
	Scal_BCM_current = Scal_BCM4A_current;
      }
      else if(bcm_type=="BCM4B"){
	Scal_BCM_current = Scal_BCM4B_current;
      }
      else if(bcm_type=="BCM4C"){
	Scal_BCM_current = Scal_BCM4C_current;
      }

      //Store Cumulative Quantities
      total_charge_bcm1 = Scal_BCM1_charge;
      total_charge_bcm2 = Scal_BCM2_charge;
      total_charge_bcm4a = Scal_BCM4A_charge;
      total_charge_bcm4b = Scal_BCM4B_charge;
      total_charge_bcm4c = Scal_BCM4C_charge;
      total_time = Scal_time;
      total_s1x_scaler = S1X_scaler;
      total_s1y_scaler = S1Y_scaler;
      total_s2x_scaler = S2X_scaler;
      total_s2y_scaler = S2Y_scaler;
      total_trig1_scaler = TRIG1_scaler;
      total_trig2_scaler = TRIG2_scaler;
      total_trig3_scaler = TRIG3_scaler;
      total_trig4_scaler = TRIG4_scaler;
      total_trig5_scaler = TRIG5_scaler;
      total_trig6_scaler = TRIG6_scaler;
      total_edtm_scaler = EDTM_scaler;
//std::cout << "Scal_time = " << Scal_time << std::endl;
// Print the running total of the charge
//std::cout << "Total BCM4A Charge so far = " << Scal_BCM4A_charge << std::endl;
//std::cout << "Total BCM4C Charge so far = " << Scal_BCM4C_charge << std::endl;
//std::cout << "Total BCM4A Current is = " << Scal_BCM4A_current << std::endl;
//std::cout << "Total BCM Current so far = " << Scal_BCM_current << std::endl;
   
    //TString bcm_type="BCM4A";
      if(Scal_BCM_current > bcm_thrs)
	{
    //std::cout << "Scal_time = " << (Scal_time - prev_time) << std::endl;
	  //std::cout << "Total BCM Current so far = " << Scal_BCM_current << std::endl;
    
	  //Turn Event Flag ON, if beam current is within threshold
	  evt_flag_bcm[i] = 1;
	  
	  //Store Quantities that Passed the Current Threshold
	  total_time_bcm_cut = total_time_bcm_cut + (Scal_time - prev_time);
	  total_charge_bcm1_cut = total_charge_bcm1_cut + (Scal_BCM1_charge - prev_charge_bcm1); 
    total_charge_bcm2_cut = total_charge_bcm2_cut + (Scal_BCM2_charge - prev_charge_bcm2);  
	  total_charge_bcm4a_cut = total_charge_bcm4a_cut + (Scal_BCM4A_charge - prev_charge_bcm4a);  
	  total_charge_bcm4b_cut = total_charge_bcm4b_cut + (Scal_BCM4B_charge - prev_charge_bcm4b);  
	  total_charge_bcm4c_cut = total_charge_bcm4c_cut + (Scal_BCM4C_charge - prev_charge_bcm4c);  
	  total_s1x_scaler_bcm_cut = total_s1x_scaler_bcm_cut + (S1X_scaler-prev_s1x_scaler);
	  total_s1y_scaler_bcm_cut = total_s1y_scaler_bcm_cut + (S1Y_scaler-prev_s1y_scaler);
	  total_s2x_scaler_bcm_cut = total_s2x_scaler_bcm_cut + (S2X_scaler-prev_s2x_scaler);
	  total_s2y_scaler_bcm_cut = total_s2y_scaler_bcm_cut + (S2Y_scaler-prev_s2y_scaler);
	  total_trig1_scaler_bcm_cut = total_trig1_scaler_bcm_cut + (TRIG1_scaler-prev_trig1_scaler);
	  total_trig2_scaler_bcm_cut = total_trig2_scaler_bcm_cut + (TRIG2_scaler-prev_trig2_scaler);
	  total_trig3_scaler_bcm_cut = total_trig3_scaler_bcm_cut + (TRIG3_scaler-prev_trig3_scaler);
	  total_trig4_scaler_bcm_cut = total_trig4_scaler_bcm_cut + (TRIG4_scaler-prev_trig4_scaler);
	  total_trig5_scaler_bcm_cut = total_trig5_scaler_bcm_cut + (TRIG5_scaler-prev_trig5_scaler);
	  total_trig6_scaler_bcm_cut = total_trig6_scaler_bcm_cut + (TRIG6_scaler-prev_trig6_scaler);
	  total_edtm_scaler_bcm_cut = total_edtm_scaler_bcm_cut + (EDTM_scaler - prev_edtm_scaler);

	} //End BCM Current Cut
    /*
std::cout << "Total BCM4A Charge cut so far = " << total_charge_bcm4a_cut << std::endl;
std::cout << "Total trigger 1 Scaler BCM cut so far = " << total_trig1_scaler_bcm_cut << std::endl;
    std::cout << "Total EDTM Scaler BCM cut so far = " << total_edtm_scaler_bcm_cut << std::endl;
    /*/
     //Previous Scaler Reads (Necessary to Take Average between S-1 and S scaler reads, to get values in between)
      prev_time = Scal_time;
      prev_charge_bcm1 = Scal_BCM1_charge;
      prev_charge_bcm2 = Scal_BCM2_charge;
      prev_charge_bcm4a = Scal_BCM4A_charge;
      prev_charge_bcm4b = Scal_BCM4B_charge;
      prev_charge_bcm4c = Scal_BCM4C_charge;
      prev_s1x_scaler = S1X_scaler;
      prev_s1y_scaler = S1Y_scaler;
      prev_s2x_scaler = S2X_scaler;
      prev_s2y_scaler = S2Y_scaler;
      prev_trig1_scaler = TRIG1_scaler;
      prev_trig2_scaler = TRIG2_scaler;
      prev_trig3_scaler = TRIG3_scaler;
      prev_trig4_scaler = TRIG4_scaler;
      prev_trig5_scaler = TRIG5_scaler;
      prev_trig6_scaler = TRIG6_scaler;
      prev_edtm_scaler = EDTM_scaler;
      
    } // End scaler event loop
    
      // Set generic bcm info to be used in charge normalization based on user input
  if(bcm_type=="BCM1"){
    total_charge_bcm     = total_charge_bcm1;
    total_charge_bcm_cut = total_charge_bcm1_cut;
  }
  else if(bcm_type=="BCM2"){
    total_charge_bcm     = total_charge_bcm2;
    total_charge_bcm_cut = total_charge_bcm2_cut;
  }
  else if(bcm_type=="BCM4A"){
    total_charge_bcm     = total_charge_bcm4a;
    total_charge_bcm_cut = total_charge_bcm4a_cut;
  }
  else if(bcm_type=="BCM4B"){
    total_charge_bcm     = total_charge_bcm4b;
    total_charge_bcm_cut = total_charge_bcm4b_cut;
  }
  else if(bcm_type=="BCM4C"){
    total_charge_bcm     = total_charge_bcm4c;
    total_charge_bcm_cut = total_charge_bcm4c_cut;
  }
  //std::cout << "Total BCM4C Current cut = " << total_charge_bcm4a_cut << std::endl;
  
   //Subtract EDTM counts from trigger scalers
  total_s1x_scaler_bcm_cut = total_s1x_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_s1y_scaler_bcm_cut = total_s1y_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_s2x_scaler_bcm_cut = total_s2x_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_s2y_scaler_bcm_cut = total_s2y_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig1_scaler_bcm_cut = total_trig1_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig2_scaler_bcm_cut = total_trig2_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig3_scaler_bcm_cut = total_trig3_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig4_scaler_bcm_cut = total_trig4_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig5_scaler_bcm_cut = total_trig5_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig6_scaler_bcm_cut = total_trig6_scaler_bcm_cut - total_edtm_scaler_bcm_cut;


   //Calculate Scaler Trigger Rates (EDTM subtracted already)
  S1XscalerRate_bcm_cut = total_s1x_scaler_bcm_cut / total_time_bcm_cut;
  S1YscalerRate_bcm_cut = total_s1y_scaler_bcm_cut / total_time_bcm_cut;
  S2XscalerRate_bcm_cut = total_s2x_scaler_bcm_cut / total_time_bcm_cut;
  S2YscalerRate_bcm_cut = total_s2y_scaler_bcm_cut / total_time_bcm_cut;
  TRIG1scalerRate_bcm_cut = total_trig1_scaler_bcm_cut / total_time_bcm_cut;
  TRIG2scalerRate_bcm_cut = total_trig2_scaler_bcm_cut / total_time_bcm_cut;
  TRIG3scalerRate_bcm_cut = total_trig3_scaler_bcm_cut / total_time_bcm_cut;
  TRIG4scalerRate_bcm_cut = total_trig4_scaler_bcm_cut / total_time_bcm_cut;
  TRIG5scalerRate_bcm_cut = total_trig5_scaler_bcm_cut / total_time_bcm_cut;
  TRIG6scalerRate_bcm_cut = total_trig6_scaler_bcm_cut / total_time_bcm_cut;
  EDTMscalerRate_bcm_cut =  total_edtm_scaler_bcm_cut / total_time_bcm_cut;

  // include changing the microAmperes into milliamperes

    // Scal_BCM4A_charge = Scal_BCM4A_charge/ 1000.;
     //Calculate Average BCM Current                                                  
  avg_current_bcm_cut = total_charge_bcm_cut / total_time_bcm_cut; //uA                              
  
  //Convert charge from uC to mC                                   
  total_charge_bcm_cut = total_charge_bcm_cut / 1000.; 
  total_charge_bcm1_cut  = total_charge_bcm1_cut / 1000.; 
  total_charge_bcm2_cut  = total_charge_bcm2_cut / 1000.; 
  total_charge_bcm4a_cut = total_charge_bcm4a_cut / 1000.; 
  total_charge_bcm4b_cut = total_charge_bcm4b_cut / 1000.; 
  total_charge_bcm4c_cut = total_charge_bcm4c_cut / 1000.; 

     
std::cout << "Total Charge micro coulomb = " << Scal_BCM4A_charge << std::endl;
std::cout << Form("Total Charge micro coulomb into milli coulomb = %.4f ", total_charge_bcm4a_cut) << std::endl; 
std::cout << "Total current = " << avg_current_bcm_cut << std::endl; 
std::cout <<"Total edtm scaler rate = " << EDTMscalerRate_bcm_cut << std::endl;

  //Convert Scaler Trigger/EDTM Rates from Hz to kHz 
  S1XscalerRate_bcm_cut   = S1XscalerRate_bcm_cut   / 1000.;
  S1YscalerRate_bcm_cut   = S1YscalerRate_bcm_cut   / 1000.;
  S2XscalerRate_bcm_cut   = S2XscalerRate_bcm_cut   / 1000.;
  S2YscalerRate_bcm_cut   = S2YscalerRate_bcm_cut   / 1000.;
  TRIG1scalerRate_bcm_cut = TRIG1scalerRate_bcm_cut / 1000.;
  TRIG2scalerRate_bcm_cut = TRIG2scalerRate_bcm_cut / 1000.;
  TRIG3scalerRate_bcm_cut = TRIG3scalerRate_bcm_cut / 1000.;
  TRIG4scalerRate_bcm_cut = TRIG4scalerRate_bcm_cut / 1000.;
  TRIG5scalerRate_bcm_cut = TRIG5scalerRate_bcm_cut / 1000.;
  TRIG6scalerRate_bcm_cut = TRIG6scalerRate_bcm_cut / 1000.;
  EDTMscalerRate_bcm_cut  = EDTMscalerRate_bcm_cut  / 1000.;


std::cout <<Form("Total edtm scaler rate per 1000 = %.4f ",  EDTMscalerRate_bcm_cut) << std::endl;
std::cout <<Form("Total trigger 6 scaler rate per 1000 = %.4f ", TRIG6scalerRate_bcm_cut) << std::endl;
std::cout <<Form("Total Trigger 5 scaler rate per 1000 = %.4f " , TRIG5scalerRate_bcm_cut) << std::endl;
std::cout <<Form("Total Trigger 3 scaler rate per 1000 = %.4f " , TRIG3scalerRate_bcm_cut) << std::endl;
std::cout <<Form("Total Trigger 2 scaler rate per 1000 = %.4f ", TRIG2scalerRate_bcm_cut) << std::endl;
std::cout <<Form("Total Trigger 1 scaler rate per 1000 = %.4f ", TRIG1scalerRate_bcm_cut) << std::endl;

  // ===============================================================
  
  //================
  // DATA ANALYSIS
  //=================


        // Define DATA-specific variables
       // Double_t W, Q2, e_delta, h_delta, Pm, Pmx, Pmy, Pmz, h_xptar, h_yptar, e_xptar, e_yptar;
       // Double_t Em , x_bj, scat_ang_rad_data, xangle_data, h_theta_data, gevnum, gevtyp; 
       //Double_t etar_x, etar_y, etar_z, htar_x, htar_y, htar_z;
        // Set branch addresses
       tree->SetBranchAddress("P.kin.primary.W", &W);
       tree->SetBranchAddress("g.evnum", &gevnum);
       tree->SetBranchAddress("g.evtyp", &gevtyp);
        tree->SetBranchAddress("P.kin.primary.Q2", &Q2);
        tree->SetBranchAddress("H.gtr.dp", &h_delta);
        tree->SetBranchAddress("P.gtr.dp", &e_delta);
        tree->SetBranchAddress("H.kin.secondary.Prec_x", &Pmx);
        tree->SetBranchAddress("H.kin.secondary.Prec_y", &Pmy);
        tree->SetBranchAddress("H.kin.secondary.Prec_z", &Pmz);
        tree->SetBranchAddress("H.gtr.th", &h_xptar);
        tree->SetBranchAddress("H.gtr.ph", &h_yptar);
        tree->SetBranchAddress("P.gtr.th", &e_xptar);
        tree->SetBranchAddress("P.gtr.ph", &e_yptar);
        tree->SetBranchAddress("H.gtr.p", &Pmea);
        tree->SetBranchAddress("H.kin.secondary.emiss", &Em);
        tree->SetBranchAddress("H.kin.secondary.pmiss", &Pm);
        tree->SetBranchAddress("H.kin.secondary.xangle", &xangle_data);
        tree->SetBranchAddress("P.kin.primary.scat_ang_rad", &theta_e_rad);
        tree->SetBranchAddress("P.kin.primary.x_bj", &x_bj);
        tree->SetBranchAddress("H.extcor.xsieve", &hXColl);
        tree->SetBranchAddress("H.extcor.ysieve", &hYColl);
        tree->SetBranchAddress("P.extcor.xsieve", &eXColl);
        tree->SetBranchAddress("P.extcor.ysieve", &eYColl);
        tree->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&EDTM_tdcTimeRaw);
         
         //Coin Time
  tree->SetBranchAddress("CTime.epCoinTime_ROC2", &epCoinTime);
  tree->SetBranchAddress("CTime.CoinTime_RAW_ROC2_NoTrack",  &epCoinTime_notrk); // used to cut on E/P (etotnorm) to account for multi-track events (without bias)

  //tree->SetBranchAddress("CTime.epCoinTime_ROC2_center", &epCoinTime_center);
  //tree->SetBranchAddress("CTime.CoinTime_RAW_ROC2_NoTrack_center", &epCoinTime_center_notrk);

         // Trigger Detector 
  tree->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw",&TRIG1_tdcTimeRaw);
  tree->SetBranchAddress("T.coin.pTRIG2_ROC2_tdcTimeRaw",&TRIG2_tdcTimeRaw);
  tree->SetBranchAddress("T.coin.pTRIG3_ROC2_tdcTimeRaw",&TRIG3_tdcTimeRaw);
  tree->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw",&TRIG4_tdcTimeRaw);
  tree->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcTimeRaw",&TRIG5_tdcTimeRaw);
  tree->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTimeRaw",&TRIG6_tdcTimeRaw);
  //SHMS DETECTORS
  tree->SetBranchAddress("P.ngcer.npeSum",       &pngcer_npesum);
  tree->SetBranchAddress("P.cal.etotnorm",       &pcal_etotnorm);
  tree->SetBranchAddress("P.cal.etottracknorm",  &pcal_etottracknorm);
  tree->SetBranchAddress("P.hod.betanotrack",    &phod_beta_ntrk);
  tree->SetBranchAddress("P.hod.beta",           &phod_beta);
  tree->SetBranchAddress("P.hod.goodscinhit",    &phod_GoodScinHit);    
  tree->SetBranchAddress("P.dc.ntrack",          &pdc_ntrack);
 // tree->SetBranchAddress("P.dc.TheRealGolden",   &pdc_TheRealGolden);
  //HMS DETECTORS
  tree->SetBranchAddress("H.cer.npeSum",         &hcer_npesum);
  tree->SetBranchAddress("H.cal.etot",           &hcal_etot);
  tree->SetBranchAddress("H.cal.etotnorm",       &hcal_etotnorm);
  tree->SetBranchAddress("H.cal.etottracknorm",  &hcal_etottracknorm);
  tree->SetBranchAddress("H.hod.betanotrack",    &hhod_beta_ntrk);
  tree->SetBranchAddress("H.hod.beta",           &hhod_beta);
  tree->SetBranchAddress("H.hod.goodscinhit",    &hhod_GoodScinHit);    
  tree->SetBranchAddress("H.dc.ntrack",          &hdc_ntrack);
  //tree->SetBranchAddress("H.dc.TheRealGolden",   &hdc_TheRealGolden);

        //(tarx, tary, tarz) in Hall Coord. System
	tree->SetBranchAddress("P.react.x", &etar_x);
	tree->SetBranchAddress("P.react.y", &etar_y);
	tree->SetBranchAddress("P.react.z", &etar_z);
  tree->SetBranchAddress("H.react.x", &htar_x);
	tree->SetBranchAddress("H.react.y", &htar_y);
	tree->SetBranchAddress("H.react.z", &htar_z);
  tree->SetBranchAddress("H.dc.yp_fp", &h_ypfp);
    tree->SetBranchAddress("H.dc.y_fp", &h_yfp);
    tree->SetBranchAddress("H.dc.xp_fp", &h_xpfp);
    tree->SetBranchAddress("H.dc.x_fp", &h_xfp);
        //tree->SetBranchAddress("xangle_data - scat_ang_rad_data", &h_theta_data);
       
        // Define histograms
         // Define histograms
        H_W = new TH1F("W", "W", 100, 0.5, 1.7);
        H_Q2 = new TH1F("Q2", "Q2", 100, 1.3, 6);
        H_h_delta = new TH1F("h_delta", "h_delta", 100, -15, 15);
        H_e_delta = new TH1F("e_delta", "e_delta", 100, -12, 12);
        H_Pm = new TH1F("Pm", "Pm", 100, -0.1, 0.1);
        H_Pmx = new TH1F("Pmx", "Pmx", 100, -0.1, 0.1);
        H_Pmy = new TH1F("Pmy", "Pmy", 100, -0.1, 0.1);
        H_Pmz = new TH1F("Pmz", "Pmz", 100, -0.1, 0.2);
        H_h_xptar = new TH1F("h_xptar", "h_xptar", 100, -0.2, 0.2);
        H_h_yptar = new TH1F("h_yptar", "h_yptar", 100, -0.06, 0.06);
        H_e_xptar = new TH1F("e_xptar", "e_xptar", 100, -0.03, 0.03);
        H_e_yptar = new TH1F("e_yptar", "e_yptar", 100, -0.02, 0.02);
        H_Em = new TH1F("Em", "Em", 100, -0.1, 0.5);
        H_x_bj = new TH1F("x_bj", "x_bj", 100, 0.85,1.05);
        H_theta_e= new TH1F("theta_e_rad", "theta_e_rad", 100, 5.5,20);
        H_xangle_data = new TH1F("xangle_data", "xangle_data",100, 25,50);
        H_theta_p = new TH1F("theta_p", "theta_p", 100, 30, 50);
        H_ztarDiff = new TH1F("ztarDiff", "ztar Difference(htar_z - etar_z); ztarDiff; Events", 100, -2.5,3);
        TH1F* H_Pmea = new TH1F("Pmea", "Pmea", 100, 1,4);
        TH1F* H_Pcalc = new TH1F("Pcalc", "Pcalc", 100, 1,4);
        TH1F* H_deltaPdata = new TH1F("deltaPdata", "(Pcalc - Pmea); deltaPdata", 100, -0.1,0.1);
        
        //TH2F* H_hXColl_vs_hYColl = new TH2F("H_hXColl_vs_hYColl", "HMS Collimator;  Y-Collimator [cm];  X-Collimator [cm]", 100, -15, 15,  100, -15, 15);
       // TH2F* H_Em_vs_Pm = new TH2F("H_Em_vs_Pm", "Em_vs_Pm;  Pm;  Em", 100, -0.05, 0.05,  100, -0.01, 0.55);
        TH1F* H_deltaPdata_frac = new TH1F("deltaPdata_frac", "(Pcalc-Pmea)/Pmea; deltaPdata", 100, -0.1,0.1);
       
       TH2F* H_deltaPdata_vs_xfp = new TH2F("deltaPdata_frac_vs_xfp", "deltaPdata_frac vs xfp; xfp; (Pcalc-Pmea)/Pmea", 100, -30, 30,  100, -0.05, 0.05);
       TH2F* H_deltaPdata_vs_xpfp = new TH2F("deltaPdata_frac_vs_xpfp", "deltaPdata_frac vs xpfp; xpfp; (Pcalc-Pmea)/Pmea", 100, -0.04, 0.04,  100, -0.05, 0.05);
       TH2F* H_deltaPdata_vs_yfp = new TH2F("deltaPdata_frac_vs_yfp", "deltaPdata_frac vs yfp; yfp; (Pcalc-Pmea)/Pmea", 100, -5, 17,  100, -0.05, 0.05);
       TH2F* H_deltaPdata_vs_ypfp = new TH2F("deltaPdata_frac_vs_ypfp", "deltaPdata_frac vs ypfp; ypfp; (Pcalc-Pmea)/Pmea", 100, -0.01, 0.02,  100, -0.05, 0.05);
       
        TH2F* H_hXColl_vs_hYColl = new TH2F("H_hXColl_vs_hYColl", "HMS Collimator;  Y-Collimator [cm];  X-Collimator [cm]", 100, -15, 15,  100, -15, 15);
        TH2F* H_Em_vs_Pm = new TH2F("H_Em_vs_Pm", "Em_vs_Pm;  Pm;  Em", 100, -0.05, 0.05,  100, -0.01, 0.55);
        
         //  Histogram to fill sample coin. time 
  TH1F *H_ctime_peak = new TH1F("H_ctime_peak", "Coin. Time Peak ", 200,-100,100);

  TH1F *H_ctime_peak_notrk = new TH1F("H_ctime_peak_notrk", "Coin. Time Peak (no track) ", 200,-100,100);

  
  
//Histograms without cuts
 TH1F *H_ctime_peak_centre = new TH1F("H_ctime_peak_centre", "ep Coincidence Time; ep Coincidence Time [ns]; Counts ", 200, -25, 25);
//  TH1F *H_the_noCUT      = new TH1F("H_the_noCUT", "Electron Scattering Angle, #theta_{e}", 70, 6, 18);
  TH1F *H_W_noCUT        = new TH1F("H_W_noCUT", "Invariant Mass, W", 100, 0.5, 1.7); 
  TH1F *H_Q2_noCUT       = new TH1F("H_Q2_noCUT","4-Momentum Transfer, Q^{2}", 100, 1, 6);
  TH1F *H_xbj_noCUT      = new TH1F("H_xbj_noCUT", "x-Bjorken", 90, 0.3, 1.7);
  TH1F *H_nu_noCUT       = new TH1F("H_nu_noCUT","Energy Transfer, #nu", 80, 0, 6); 
TH1F *H_Em_noCUT       = new TH1F("H_Em_noCUT","Missing Energy, Em", 70, -0.1, 0.250); 
TH1F *H_Pm_noCUT       = new TH1F("H_Pm_noCUT","Missing Momentum, Pm", 50, 0, 0.2); 
//Detectors
 // detector
 TH1F *H_pCalEtotTrkNorm_noCUT = new TH1F("H_pCalEtotTrkNorm_noCUT", "SHMS Calorimeter Total Normalized Track Energy; E_{tot} / P_{trk}; Counts ", 100, 0.01, 3.5);
  TH1F *H_pHodBetaTrk_noCUT     = new TH1F("H_pBetaTrk_noCUT", "SHMS Hodo #beta (golden track); #beta (golden track); Counts ", 100, 0.01, 1.5);
  TH1F *H_pNGCerNpeSum_noCUT    = new TH1F("H_pNGCerNpeSum_noCUT", "SHMS Noble Gas Cherenkov NPE Sum; Cherenkov NPE Sum; Counts  ", 100, 0.001, 25);
  TH1F *H_hCalEtotTrkNorm_noCUT = new TH1F("H_hCalEtotTrkNorm_noCUT", "HMS Calorimeter Total Normalized Track Energy; E_{tot} / P_{trk}; Counts ", 100, 0.001, 3.5);
  TH1F *H_hHodBetaTrk_noCUT     = new TH1F("H_hHodBetaTrk_noCUT", "HMS Hodo #beta (golden track); #beta (golden track); Counts ", 100, 0.01, 1.5);
  TH1F *H_hCerNpeSum_noCUT      = new TH1F("H_hCerNpeSum_noCUT", "HMS Cherenkov NPE Sum; Cherenkov NPE Sum; Counts ", 100, 0.001, 25);

        // ------ Add histograms to the list -----

        // kinematics
        kin_HList->Add(H_W);
        kin_HList->Add(H_Q2);
        kin_HList->Add(H_Em);
        kin_HList->Add(H_h_delta);
        kin_HList->Add(H_e_delta);
        kin_HList->Add(H_Pmx);
        kin_HList->Add(H_Pmy);
        kin_HList->Add(H_Pmz);
        kin_HList->Add(H_h_xptar);
        kin_HList->Add(H_h_yptar);
        kin_HList->Add(H_e_xptar);
        kin_HList->Add(H_e_yptar);
        kin_HList->Add(H_Pm);
        kin_HList->Add(H_Pmea);
        kin_HList->Add(H_Pcalc);
        kin_HList->Add(H_deltaPdata);
        kin_HList->Add(H_xangle_data);
        kin_HList->Add(H_theta_e);
        kin_HList->Add(H_x_bj);
        kin_HList->Add(H_theta_p);
        kin_HList->Add(H_ztarDiff);
        kin_HList->Add(H_ctime_peak);
        kin_HList->Add(H_ctime_peak_notrk);
        // acceptance (need to create accp_HList)
        kin_HList->Add(H_hXColl_vs_hYColl);
        kin_HList->Add(H_Em_vs_Pm);
        kin_HList->Add(H_deltaPdata_frac);
        kin_HList->Add(H_deltaPdata_vs_xpfp);
        kin_HList->Add(H_deltaPdata_vs_xfp);
        kin_HList->Add(H_deltaPdata_vs_yfp);
        kin_HList->Add(H_deltaPdata_vs_ypfp);
       // kin_HList->Add(hms_beta_peak);
       // kin_HList->Add(shms_beta_peak);
       // kin_HList->Add(shms_ecal_peak);

        //quality check
        kin_HList->Add(H_hCerNpeSum_noCUT);
        kin_HList->Add(H_hHodBetaTrk_noCUT);
        kin_HList->Add(H_hCalEtotTrkNorm_noCUT);
        kin_HList->Add(H_pNGCerNpeSum_noCUT);
        kin_HList->Add(H_pHodBetaTrk_noCUT);
        kin_HList->Add( H_pCalEtotTrkNorm_noCUT);

        //Quality check for Detectors
     
     //kin_HList->Add(H_the_noCUT);
      kin_HList->Add(H_W_noCUT);
      kin_HList->Add(H_xbj_noCUT);
      kin_HList->Add(H_nu_noCUT);
       kin_HList->Add(H_Em_noCUT);
       kin_HList->Add(H_Pm_noCUT);
       kin_HList->Add(H_Q2_noCUT);
      kin_HList->Add(H_ctime_peak_centre);
        // acceptance (need to create accp_HList)
        kin_HList->Add(H_hXColl_vs_hYColl);
        kin_HList->Add(H_Em_vs_Pm);


        //=====================
        // Event loop for DATA
        //======================
      
            CollimatorStudy();
             
 
  Long64_t sample_entries = 200000;
  for(int ientry=0; ientry<sample_entries; ientry++)
    {	  
      tree->GetEntry(ientry);
      
      //cout << "theRealGolden = " << pdc_TheRealGolden << endl;
      // Fill sample histo to find peak
      H_ctime_peak->Fill(epCoinTime);
      H_ctime_peak_notrk->Fill(epCoinTime_notrk);
     // shms_beta_peak->Fill(phod_beta);
     // hms_beta_peak->Fill(hhod_beta);
     // shms_ecal_peak->Fill(pcal_etottracknorm);
      H_pCalEtotTrkNorm_noCUT->Fill(pcal_etottracknorm);
      H_pHodBetaTrk_noCUT->Fill(phod_beta);
      H_pNGCerNpeSum_noCUT->Fill(pngcer_npesum);
      H_hCalEtotTrkNorm_noCUT->Fill(hcal_etottracknorm);
      H_hHodBetaTrk_noCUT->Fill(hhod_beta);
      H_hCerNpeSum_noCUT->Fill(hcer_npesum);
     // H_the_noCUT->Fill(scat_ang_data_deg);
      H_W_noCUT->Fill(W);
      H_Q2_noCUT->Fill(Q2);
      H_xbj_noCUT->Fill(x_bj);
     H_nu_noCUT->Fill(nu);
      H_Pm_noCUT->Fill( Pm);
      H_Em_noCUT->Fill( Em);
    H_ctime_peak_centre->Fill( epCoinTime-ctime_offset_peak_val );
      cout << "SampleEventLoop: " << std::setprecision(2) << double(ientry) / sample_entries * 100. << "  % " << std::flush << "\r";
    
  
  
  // bin number corresponding to maximum bin content
  int binmax_ctime = H_ctime_peak->GetMaximumBin();
  int binmax_ctime_notrk = H_ctime_peak_notrk->GetMaximumBin();
  //int binmax_pbeta =  shms_beta_peak->GetMaximumBin();
 // int binmax_hbeta =  hms_beta_peak->GetMaximumBin();
 // int binmax_pecal = shms_ecal_peak->GetMaximumBin();
  
  // x-value corresponding to bin number with max content (i.e., peak)
  double xmax_ctime = H_ctime_peak->GetXaxis()->GetBinCenter(binmax_ctime);
  double xmax_ctime_notrk = H_ctime_peak_notrk->GetXaxis()->GetBinCenter(binmax_ctime_notrk);
//  double xmax_pbeta = shms_beta_peak->GetXaxis()->GetBinCenter(binmax_pbeta);
//  double xmax_hbeta = hms_beta_peak->GetXaxis()->GetBinCenter(binmax_hbeta);
 // double xmax_pecal = shms_ecal_peak->GetXaxis()->GetBinCenter(binmax_pecal);

  ctime_offset_peak_val        = xmax_ctime;
  ctime_offset_peak_notrk_val  = xmax_ctime_notrk;      
 // hms_beta_peak_val            = xmax_hbeta;
 // shms_beta_peak_val           = xmax_pbeta;
//  shms_ecal_peak_val           = xmax_pecal;
  
}

     
 for (Long64_t ientry = 0; ientry < nentries; ientry++) {
            tree->GetEntry(ientry);


           // Convert radian angles to degrees
            Double_t xangle_data_deg = xangle_data * TMath::RadToDeg();
            Double_t theta_e_deg = theta_e_rad * TMath::RadToDeg();
            Double_t theta_p = (xangle_data - theta_e_rad);
            Double_t theta_p_deg = theta_p * TMath::RadToDeg();
            Double_t ztarDiff = htar_z - etar_z;
            Double_t Pcalc = (2 * Mp * Eb * (Mp + Eb) * TMath::Cos(theta_p)) / 
                 (Mp * Mp + 2 * Mp * Eb + Eb * Eb * TMath::Sin(theta_p) * TMath::Sin(theta_p));
            Double_t deltaPdata = (Pcalc - Pmea);
            Double_t deltaPdata_frac = (Pcalc - Pmea)/Pmea;


    
      // define cuts
      hmsColl_Cut = hms_Coll_gCut->IsInside(hYColl, hXColl);
      Em_cut = Em>=-0.05 && Em<=0.05;


          //CUTS USED IN EDTM LIVE TIME CALCULATION
	  c_noedtm = EDTM_tdcTimeRaw == 0.;
	  c_edtm   = EDTM_tdcTimeRaw  > 0.;
	  c_trig1  = TRIG1_tdcTimeRaw > 0.;
	  c_trig2  = TRIG2_tdcTimeRaw > 0.;
	  c_trig3  = TRIG3_tdcTimeRaw > 0.;
	  c_trig4  = TRIG4_tdcTimeRaw > 0.;
	  c_trig5  = TRIG5_tdcTimeRaw > 0.;
	  c_trig6  = TRIG6_tdcTimeRaw > 0.;

	  c_notrig1  = TRIG1_tdcTimeRaw == 0.;
	  c_notrig2  = TRIG2_tdcTimeRaw == 0.;
	  c_notrig3  = TRIG3_tdcTimeRaw == 0.;
	  c_notrig4  = TRIG4_tdcTimeRaw == 0.;
	  c_notrig5  = TRIG5_tdcTimeRaw == 0.;
	  c_notrig6  = TRIG6_tdcTimeRaw == 0.;



//-------------------------------------------------------------------------
    // ----------CUTS USED IN TRACKING EFFICIENCY CALCULATION -------
//-------------------------------------------------------------------------

//CUTS: HMS TRACKING EFFICIENCY
hdc_ntrk_cut =hdc_ntrack >= 1;
hScinGood_cut = hhod_GoodScinHit==1 ;
hcer_NPE_Sum_cut = hcer_npesum >= 0. && hcer_npesum <= 0.5;
hetotnorm_cut = hcal_etotnorm >= 0. && hcal_etotnorm <= 0.6;
hBeta_notrk_cut = hhod_beta_ntrk >= 0.5 && hhod_beta_ntrk <= 1.5;
//electrons (or hadrons) that 'SHOULD' have passed the cuts to form a track

good_hms_should = hScinGood_cut && hcer_NPE_Sum_cut && hetotnorm_cut && hBeta_notrk_cut;
//electrons (or hadrons) that 'DID' passed the cuts to form a track
good_hms_did = hdc_ntrk_cut && good_hms_should;

//CUTS: SHMS TRACKING EFFICIENCY

pdc_ntrk_cut =pdc_ntrack >= 1;
pScinGood_cut = phod_GoodScinHit==1 ;
pngcer_NPE_Sum_cut = pngcer_npesum >= 1.0 && pngcer_npesum <= 100.;
if(phgcer_cut_flag){phgcer_NPE_Sum_cut = phgcer_npesum >= 0. && phgcer_npesum <= 1.5;}
else{phgcer_NPE_Sum_cut = 1;}
petotnorm_cut = pcal_etotnorm >= 0.8 && pcal_etotnorm <= 1.3 ;
pBeta_notrk_cut = phod_beta_ntrk >= 0.5 && phod_beta_ntrk <= 1.5;
//electrons (or hadrons) that 'SHOULD' have passed the cuts to form a track

good_shms_should = pScinGood_cut && pngcer_NPE_Sum_cut && phgcer_NPE_Sum_cut && petotnorm_cut && pBeta_notrk_cut;
//electrons (or hadrons) that 'DID' passed the cuts to form a track
good_shms_did = pdc_ntrk_cut && good_shms_should;

// CUTS (SPECIFIC TO DATA)
	 
	  // -- proton coincidence time cut ----
	  

	    // main coincidence time window cut
	    eP_ctime_cut = (epCoinTime-ctime_offset_peak_val) >= -2.5 && (epCoinTime-ctime_offset_peak_val) <= 2.5;
// PID CUTS 
 //SHMS calorimeter total normalized track energy
	  pid_petot_trkNorm_cut = pcal_etottracknorm>=0.8 && pcal_etottracknorm<=1.3;
	 
     // End Defined Cuts

	  //Count Accepted EDTM events (no bcm current cut : this is just to compare with counts that have bcm cuts)
	  if(c_edtm){total_edtm_accp++;}
	  
	  //Count Accepted TRIG 1-6 events (no bcm current cut : this is just to compare with counts that have bcm cuts)
	  if(c_trig1){total_trig1_accp++;}
	  if(c_trig2){total_trig2_accp++;}
	  if(c_trig3){total_trig3_accp++;}
	  if(c_trig4){total_trig4_accp++;}
	  if(c_trig5){total_trig5_accp++;}
	  if(c_trig6){total_trig6_accp++;}
	  
	  //----------------------Check If BCM Current is within limits---------------------

    // Apply BCM Current Cut (check if BCM current passed threshold)
	      if(evt_flag_bcm[scal_read]==1) {

        if(c_edtm){ total_edtm_accp_bcm_cut++;}
	      
	      //Count Accepted TRIG1-6 events (without EDTM and with bcm current cut: to be used in the computer live time calculation)
	      if(c_trig1 && c_noedtm) { total_trig1_accp_bcm_cut++; }
	      if(c_trig2 && c_noedtm) { total_trig2_accp_bcm_cut++; }
	      if(c_trig3 && c_noedtm) { total_trig3_accp_bcm_cut++; }
	      if(c_trig4 && c_noedtm) { total_trig4_accp_bcm_cut++; }
	      if(c_trig5 && c_noedtm) { total_trig5_accp_bcm_cut++; }
	      if(c_trig6 && c_noedtm) { total_trig6_accp_bcm_cut++; }
	      
 //REQUIRE "NO EDTM" CUT TO FILL DATA HISTOGRAMS
	      if(c_noedtm)
		{
		  //cout << "passed NO EDTM Cut !" << endl;   
		  //Calculate HMS Tracking Efficiency Components
		  if(good_hms_did){ h_did++;}
		  if(good_hms_should){ h_should++; }
		  
		  //Calculate SHMS Tracking Efficiency Components
		  if(good_shms_did){ p_did++;}
		  if(good_shms_should){ p_should++; }
    
	    
	     if(hmsColl_Cut && Em_cut && pid_petot_trkNorm_cut && eP_ctime_cut && abs(h_delta)<= 10 &&e_delta>= -10 && e_delta <= 22)
		    {
           H_W->Fill(W);
            H_Q2->Fill(Q2);  
            H_h_delta->Fill(h_delta);
            H_e_delta->Fill(e_delta);
            H_Pmx->Fill(Pmx);
            H_Pmy->Fill(Pmy);
            H_Pmz->Fill(Pmz);
            H_e_xptar->Fill(e_xptar);
            H_e_yptar->Fill(e_yptar);
            H_h_xptar->Fill(h_xptar);
            H_h_yptar->Fill(h_yptar);
            H_Em->Fill(Em);
            H_x_bj->Fill(x_bj);
            H_xangle_data->Fill(xangle_data_deg);
            H_theta_e->Fill(theta_e_deg);
            H_theta_p->Fill(theta_p_deg);
            H_ztarDiff->Fill(ztarDiff);
            H_Pmea->Fill(Pmea);
            H_Pcalc->Fill(Pcalc);
            H_deltaPdata->Fill(deltaPdata);
            H_hXColl_vs_hYColl->Fill(hYColl, hXColl);
            H_Em_vs_Pm->Fill(Em,Pm);
            H_deltaPdata_frac->Fill(deltaPdata_frac);
            H_deltaPdata_vs_xfp->Fill( h_xfp, deltaPdata_frac);
            H_deltaPdata_vs_xpfp->Fill(h_xpfp, deltaPdata_frac);
            H_deltaPdata_vs_yfp->Fill(h_yfp, deltaPdata_frac);
            H_deltaPdata_vs_ypfp->Fill(h_ypfp, deltaPdata_frac);
        } // Ends data cuts

     } // end NOT edtm cut
            std::cout << "DATA Events Completed: " << std::setprecision(2)
                      << double(ientry) / nentries * 100. << " % " << std::flush << "\r";


    } // End BCM Current Cut

      

	  //gevnum = ientry + 1;

	  if(gevnum==scal_evt_num[scal_read]){ scal_read++; }
	  
//cout << "gevnum = " << gevnum << endl; 

  
    } // End DATA EVENT LOOP

    // Define the cut and count the number of events
   // const char * pdc_ntrk_cut ="pdc_ntrack >= 1";
   // const char * pScinGood_cut ="phod_GoodScinHit==1";
   // const char * pngcer_NPE_Sum_cut ="pngcer_npesum >= 1.0 && pngcer_npesum <= 100.";
   // const char * phgcer_NPE_Sum_cut ="phgcer_npesum >= 0. && phgcer_npesum <= 1.5";
  // // const char * petotnorm_cut ="pcal_etotnorm >= 0.8 && pcal_etotnorm <= 1.3 ";
   //const char *pBeta_notrk_cut = "phod_beta_ntrk >= 0.5 && phod_beta_ntrk <= 1.5";
   //const char*good_shms_should = "pScinGood_cut && pngcer_NPE_Sum_cut && phgcer_NPE_Sum_cut && petotnorm_cut && pBeta_notrk_cut";
   // Long64_t nEvents_pdc = tree->GetEntries(pdc_ntrk_cut);
   // Long64_t nEvents_pScin = tree->GetEntries(pScinGood_cut);
   // Long64_t nEvents_pngcer = tree->GetEntries(pngcer_NPE_Sum_cut);
   // Long64_t nEvents_phgcer = tree->GetEntries(phgcer_NPE_Sum_cut);
   // Long64_t nEvents_petot = tree->GetEntries(petotnorm_cut);
   // Long64_t nEvents_pBeta = tree->GetEntries(pBeta_notrk_cut);
   // Long64_t nEvents_should = tree->GetEntries(good_shms_should );


    // Print the result
   // std::cout << "Number of events passing the pdc ntrack cut: " << nEvents_pdc << std::endl;
// std::cout << "Number of events passing the ScinGood cut: " << nEvents_pScin << std::endl;
 //std::cout << "Number of events passing the Inert Gas Cerenkov cut: " << nEvents_pngcer << std::endl;
// std::cout << "Number of events passing the Heavy Gas Cerenkov cut: " << nEvents_phgcer << std::endl;
// std::cout << "Number of events passing the Calorimeter cut: " << nEvents_petot << std::endl;
 //std::cout << "Number of events passing the Beta Track cut: " << nEvents_pBeta << std::endl;
// std::cout << "Number of events passing the SHMS Should cut: " << nEvents_should << std::endl;

 //Count the EFFECT of CUTS

   // Initialize counters
double_t total_should_no_cut = 0;
double_t total_did_pdc_ntrk = 0;
double_t total_should_pScin =0;
double_t total_did_pScin_pdc =0;
double_t total_should_pScin_ngcer = 0;
double_t total_did_pScin_ngcer_pdc = 0;
double_t total_should_pScin_ngcer_petot = 0;
double_t total_did_pScin_ngcer_petot_pdc = 0;
double_t total_should_all_cuts = 0;
double_t total_did_all_cuts = 0;

// Loop through the events
Long64_t nEntries = tree->GetEntries();
for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
// Initialize conditions
    bool good_shms_should = true; // Placeholder for your initial condition
pdc_ntrk_cut =pdc_ntrack >= 1;
pScinGood_cut = phod_GoodScinHit==1 ;
pngcer_NPE_Sum_cut = pngcer_npesum >= 1.0 && pngcer_npesum <= 100.;
petotnorm_cut = pcal_etotnorm >= 0.8 && pcal_etotnorm <= 1.3 ;
pBeta_notrk_cut = phod_beta_ntrk >= 0.5 && phod_beta_ntrk <= 1.5;


    // Count events for each condition
    if (good_shms_should) {
        total_should_no_cut++;
    }

    if (good_shms_should && pdc_ntrk_cut) {
        total_did_pdc_ntrk++;
    }

 if ( pScinGood_cut ) {
        total_should_pScin++;
    }

    if (pScinGood_cut && pdc_ntrk_cut) {
        total_did_pScin_pdc++;
    }  
    if ( pScinGood_cut && pngcer_NPE_Sum_cut ) {
        total_should_pScin_ngcer++;
    }

    if (pScinGood_cut && pngcer_NPE_Sum_cut && pdc_ntrk_cut) {
        total_did_pScin_ngcer_pdc++;
    }

    if (pScinGood_cut && pngcer_NPE_Sum_cut && petotnorm_cut) {
        total_should_pScin_ngcer_petot++;
    }

    if (pScinGood_cut && pngcer_NPE_Sum_cut && petotnorm_cut && pdc_ntrk_cut) {
        total_did_pScin_ngcer_petot_pdc++;
    }

    if (pScinGood_cut && pngcer_NPE_Sum_cut && petotnorm_cut && pBeta_notrk_cut) {
        total_should_all_cuts++;
    }

    if (pScinGood_cut && pngcer_NPE_Sum_cut && petotnorm_cut && pBeta_notrk_cut && pdc_ntrk_cut) {
        total_did_all_cuts++;
    }
}

// Print the results
std::cout << Form("Total 'should' with no cut: = %.4f " , total_should_no_cut) << std::endl;
std::cout << Form("Total 'did' with pdc_ntrack cut: = %.4f ", total_did_pdc_ntrk) << std::endl;
std::cout << Form("Total 'should' with Scintillator hit cut:= %.4f ", total_should_pScin) << std::endl;
std::cout << Form("Total 'did' with pdc_ntrack and Good Scintillator hit cut: = %.4f" , total_did_pScin_pdc) << std::endl;
std::cout << Form("Total 'should' with good scintillator hit and noble gas cut: = %.4f " , total_should_pScin_ngcer) << std::endl;
std::cout << Form("Total 'did' with good scintillator hit, noble gas cut, and pdc_ntrack cut:=  %.4f " , total_did_pScin_ngcer_pdc) << std::endl;
std::cout << Form("Total 'should' with good scintillator hit, noble gas cut, and calorimeter cut: =%.4f " , total_should_pScin_ngcer_petot) << std::endl;
std::cout << Form("Total 'did' with good scintillator hit, noble gas cut, calorimeter cut, and pdc_ntrack cut: = %.4f ", total_did_pScin_ngcer_petot_pdc) << std::endl;
std::cout << Form("Total 'should' with all cuts: = %.4f " ,total_should_all_cuts) << std::endl;
std::cout << Form("Total 'did' with all cuts: = %.4f " , total_did_all_cuts) << std::endl;




   
 //good_shms_should = phod_GoodScinHit==1;
 

 //good_shms_should = "pScinGood_cut && pngcer_NPE_Sum_cut && phgcer_NPE_Sum_cut && petotnorm_cut && pBeta_notrk_cut";
  
//Calculate Ratesh

    //Calculate Accepted Trigger/EDTM Rates in kHz
  TRIG1accpRate_bcm_cut = (total_trig1_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG2accpRate_bcm_cut = (total_trig2_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG3accpRate_bcm_cut = (total_trig3_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG4accpRate_bcm_cut = (total_trig4_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG5accpRate_bcm_cut = (total_trig5_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG6accpRate_bcm_cut = (total_trig6_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  EDTMaccpRate_bcm_cut  = (total_edtm_accp_bcm_cut  / total_time_bcm_cut ) / 1000.;
  
 // std::cout << "Total time BCM cut = " << total_time_bcm_cut << std::endl;
 // std::cout << "total trigger 1 acceptance bcm cut = " << total_trig1_accp_bcm_cut << std::endl;
  
 
  //Calculate Pure Computer Live Time (numerator->accepted tdc trig requires NO EDTM :: denominator -> EDTM has already been subtracted from scaler counts)
  //Pre-Scale factor has been accounted 
  cpuLT_trig1 = total_trig1_accp_bcm_cut * Ps1_factor / (total_trig1_scaler_bcm_cut);
  cpuLT_trig2 = total_trig2_accp_bcm_cut * Ps2_factor / (total_trig2_scaler_bcm_cut);
  cpuLT_trig3 = total_trig3_accp_bcm_cut * Ps3_factor / (total_trig3_scaler_bcm_cut);
  cpuLT_trig4 = total_trig4_accp_bcm_cut * Ps4_factor / (total_trig4_scaler_bcm_cut);
  cpuLT_trig5 = total_trig5_accp_bcm_cut * Ps5_factor / (total_trig5_scaler_bcm_cut);
  cpuLT_trig6 = total_trig6_accp_bcm_cut * Ps6_factor / (total_trig6_scaler_bcm_cut);

  //Calculate Computer Live Time Error (Use Binomial Statistics Error formula: err^2 = N * P * (1-P), where S->total counts (or trials), and P->probability of success : accepted_triggers  / scalers 
  cpuLT_trig1_err_Bi = sqrt( total_trig1_accp_bcm_cut * (1. - (total_trig1_accp_bcm_cut / total_trig1_scaler_bcm_cut) ) ) * Ps1_factor / total_trig1_scaler_bcm_cut;
  cpuLT_trig2_err_Bi = sqrt( total_trig2_accp_bcm_cut * (1. - (total_trig2_accp_bcm_cut / total_trig2_scaler_bcm_cut) ) ) * Ps2_factor / total_trig2_scaler_bcm_cut;
  cpuLT_trig3_err_Bi = sqrt( total_trig3_accp_bcm_cut * (1. - (total_trig3_accp_bcm_cut / total_trig3_scaler_bcm_cut) ) ) * Ps3_factor / total_trig3_scaler_bcm_cut;
  cpuLT_trig4_err_Bi = sqrt( total_trig4_accp_bcm_cut * (1. - (total_trig4_accp_bcm_cut / total_trig4_scaler_bcm_cut) ) ) * Ps4_factor / total_trig4_scaler_bcm_cut;
  cpuLT_trig5_err_Bi = sqrt( total_trig5_accp_bcm_cut * (1. - (total_trig5_accp_bcm_cut / total_trig5_scaler_bcm_cut) ) ) * Ps5_factor / total_trig5_scaler_bcm_cut;
  cpuLT_trig6_err_Bi = sqrt( total_trig6_accp_bcm_cut * (1. - (total_trig6_accp_bcm_cut / total_trig6_scaler_bcm_cut) ) ) * Ps6_factor / total_trig6_scaler_bcm_cut;
 
 // std::cout << "cpuLT_trig6 = " << total_trig1_accp_bcm_cut << std::endl; 

  // why error is nan? --> total_trig5_accp_bcm_cut ~ total_trig5_scaler_bcm_cut 
  //cout << "total_trig5_accp_bcm_cut * (1. - (total_trig5_accp_bcm_cut )/total_trig5_scaler_bcm_cut ) = " << total_trig5_accp_bcm_cut * (1. - (total_trig5_accp_bcm_cut )/total_trig5_scaler_bcm_cut ) << endl;
  //cout << Form("total_trig5_accp_bcm_cut  = %.6f, total_trig5_scaler_bcm_cut  = %.6f ", total_trig5_accp_bcm_cut, total_trig5_scaler_bcm_cut) << endl; 
  
  //Calculate Computer Live Time Error (Use Bayesian Statistics Error formula)
  cpuLT_trig1_err_Bay = sqrt( (total_trig1_accp_bcm_cut * Ps1_factor + 1)*(total_trig1_accp_bcm_cut * Ps1_factor + 2)/((total_trig1_scaler_bcm_cut + 2)*(total_trig1_scaler_bcm_cut + 3)) - pow((total_trig1_accp_bcm_cut * Ps1_factor + 1),2)/pow((total_trig1_scaler_bcm_cut + 2),2) );
  cpuLT_trig2_err_Bay = sqrt( (total_trig2_accp_bcm_cut * Ps2_factor + 1)*(total_trig2_accp_bcm_cut * Ps2_factor + 2)/((total_trig2_scaler_bcm_cut + 2)*(total_trig2_scaler_bcm_cut + 3)) - pow((total_trig2_accp_bcm_cut * Ps2_factor + 1),2)/pow((total_trig2_scaler_bcm_cut + 2),2) );
  cpuLT_trig3_err_Bay = sqrt( (total_trig3_accp_bcm_cut * Ps3_factor + 1)*(total_trig3_accp_bcm_cut * Ps3_factor + 2)/((total_trig3_scaler_bcm_cut + 2)*(total_trig3_scaler_bcm_cut + 3)) - pow((total_trig3_accp_bcm_cut * Ps3_factor + 1),2)/pow((total_trig3_scaler_bcm_cut + 2),2) );
  cpuLT_trig4_err_Bay = sqrt( (total_trig4_accp_bcm_cut * Ps4_factor + 1)*(total_trig4_accp_bcm_cut * Ps4_factor + 2)/((total_trig4_scaler_bcm_cut + 2)*(total_trig4_scaler_bcm_cut + 3)) - pow((total_trig4_accp_bcm_cut * Ps4_factor + 1),2)/pow((total_trig4_scaler_bcm_cut + 2),2) );
  cpuLT_trig5_err_Bay = sqrt( (total_trig5_accp_bcm_cut * Ps5_factor + 1)*(total_trig5_accp_bcm_cut * Ps5_factor + 2)/((total_trig5_scaler_bcm_cut + 2)*(total_trig5_scaler_bcm_cut + 3)) - pow((total_trig5_accp_bcm_cut * Ps5_factor + 1),2)/pow((total_trig5_scaler_bcm_cut + 2),2) );
  cpuLT_trig6_err_Bay = sqrt( (total_trig6_accp_bcm_cut * Ps6_factor + 1)*(total_trig6_accp_bcm_cut * Ps6_factor + 2)/((total_trig6_scaler_bcm_cut + 2)*(total_trig6_scaler_bcm_cut + 3)) - pow((total_trig6_accp_bcm_cut * Ps6_factor + 1),2)/pow((total_trig6_scaler_bcm_cut + 2),2) );
  //cout << "cpuLT_trig5_err_Bay = " << cpuLT_trig5_err_Bay  << endl;
/**
// Print CPU Live Time and Errors
std::cout << "CPU Live Time for Trigger 1: " << cpuLT_trig1 << " (Error Binomial: " << cpuLT_trig1_err_Bi << ", Error Bayesian: " << cpuLT_trig1_err_Bay << ")" << std::endl;
std::cout << "CPU Live Time for Trigger 2: " << cpuLT_trig2 << " (Error Binomial: " << cpuLT_trig2_err_Bi << ", Error Bayesian: " << cpuLT_trig2_err_Bay << ")" << std::endl;
std::cout << "CPU Live Time for Trigger 3: " << cpuLT_trig3 << " (Error Binomial: " << cpuLT_trig3_err_Bi << ", Error Bayesian: " << cpuLT_trig3_err_Bay << ")" << std::endl;
std::cout << "CPU Live Time for Trigger 4: " << cpuLT_trig4 << " (Error Binomial: " << cpuLT_trig4_err_Bi << ", Error Bayesian: " << cpuLT_trig4_err_Bay << ")" << std::endl;
std::cout << "CPU Live Time for Trigger 5: " << cpuLT_trig5 << " (Error Binomial: " << cpuLT_trig5_err_Bi << ", Error Bayesian: " << cpuLT_trig5_err_Bay << ")" << std::endl;
std::cout << "CPU Live Time for Trigger 6: " << cpuLT_trig6 << " (Error Binomial: " << cpuLT_trig6_err_Bi << ", Error Bayesian: " << cpuLT_trig6_err_Bay << ")" << std::endl;
/*/
  //Calculated total EDTM Live Time
  tLT_trig1 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps1_factor);  
  tLT_trig2 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps2_factor);  
  tLT_trig3 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps3_factor);  
  tLT_trig4 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps4_factor);  
  tLT_trig5 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps5_factor);  
  tLT_trig6 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps6_factor);
//std::cout << "total livetime for trigger 6 = " << tLT_trig6 << std::endl;
std::cout << Form("total  edtm scalar  = %.4f ", total_edtm_scaler_bcm_cut) << std::endl;
std::cout << Form("total edtm accp  = %.4f", total_edtm_accp_bcm_cut) << std::endl;
std::cout << Form("total livetime for trigger 6 = %.4f ",  tLT_trig6) << std::endl;
  
  //Calculate EDTM Live Time Error (Use Binomial Error)
  tLT_trig1_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps1_factor / total_edtm_scaler_bcm_cut;
  tLT_trig2_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps2_factor / total_edtm_scaler_bcm_cut;
  tLT_trig3_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps3_factor / total_edtm_scaler_bcm_cut;
  tLT_trig4_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps4_factor / total_edtm_scaler_bcm_cut;
  tLT_trig5_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps5_factor / total_edtm_scaler_bcm_cut;
  tLT_trig6_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps6_factor / total_edtm_scaler_bcm_cut;

  //Calculate EDTM Live Time Error (Use Bayesian Error)
  tLT_trig1_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps1_factor + 1)*(total_edtm_accp_bcm_cut * Ps1_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps1_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig2_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps2_factor + 1)*(total_edtm_accp_bcm_cut * Ps2_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps2_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig3_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps3_factor + 1)*(total_edtm_accp_bcm_cut * Ps3_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps3_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig4_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps4_factor + 1)*(total_edtm_accp_bcm_cut * Ps4_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps4_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig5_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps5_factor + 1)*(total_edtm_accp_bcm_cut * Ps5_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps5_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig6_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps6_factor + 1)*(total_edtm_accp_bcm_cut * Ps6_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps6_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  
  //Ensure that if Ps_factor = -1 (trigger input OFF), then live times default to -1.0
  if(Ps1_factor==-1) { cpuLT_trig1 = -1.0, tLT_trig1 = -1.0; }
  if(Ps2_factor==-1) { cpuLT_trig2 = -1.0, tLT_trig2 = -1.0; }
  if(Ps3_factor==-1) { cpuLT_trig3 = -1.0, tLT_trig3 = -1.0; }
  if(Ps4_factor==-1) { cpuLT_trig4 = -1.0, tLT_trig4 = -1.0; }
  if(Ps5_factor==-1) { cpuLT_trig5 = -1.0, tLT_trig5 = -1.0; }
  if(Ps6_factor==-1) { cpuLT_trig6 = -1.0, tLT_trig6 = -1.0; }
  

  // select which pre-scale factors to apply in the weight (depend if looking at singles or coin)
  //if(trig_type_single=="T1") { trig_rate_single = TRIG1scalerRate_bcm_cut, cpuLT_trig_single = cpuLT_trig1, cpuLT_trig_err_Bi_single = cpuLT_trig1_err_Bi, cpuLT_trig_err_Bay_single = cpuLT_trig1_err_Bay, tLT_trig_single = tLT_trig1, tLT_trig_err_Bi_single = tLT_trig1_err_Bi, tLT_trig_err_Bay_single = tLT_trig1_err_Bay, Ps_factor_single = Ps1_factor, total_trig_scaler_bcm_cut_single = total_trig1_scaler_bcm_cut, total_trig_accp_bcm_cut_single = total_trig1_accp_bcm_cut; }
 // if(trig_type_single=="T2") { trig_rate_single = TRIG2scalerRate_bcm_cut, cpuLT_trig_single = cpuLT_trig2, cpuLT_trig_err_Bi_single = cpuLT_trig2_err_Bi, cpuLT_trig_err_Bay_single = cpuLT_trig2_err_Bay, tLT_trig_single = tLT_trig2, tLT_trig_err_Bi_single = tLT_trig2_err_Bi, tLT_trig_err_Bay_single = tLT_trig2_err_Bay, Ps_factor_single = Ps2_factor, total_trig_scaler_bcm_cut_single = total_trig2_scaler_bcm_cut, total_trig_accp_bcm_cut_single = total_trig2_accp_bcm_cut; }
 // if(trig_type_single=="T3") { trig_rate_single = TRIG3scalerRate_bcm_cut, cpuLT_trig_single = cpuLT_trig3, cpuLT_trig_err_Bi_single = cpuLT_trig3_err_Bi, cpuLT_trig_err_Bay_single = cpuLT_trig3_err_Bay, tLT_trig_single = tLT_trig3, tLT_trig_err_Bi_single = tLT_trig3_err_Bi, tLT_trig_err_Bay_single = tLT_trig3_err_Bay, Ps_factor_single = Ps3_factor, total_trig_scaler_bcm_cut_single = total_trig3_scaler_bcm_cut, total_trig_accp_bcm_cut_single = total_trig3_accp_bcm_cut; }
 // if(trig_type_single=="T4") { trig_rate_single = TRIG4scalerRate_bcm_cut, cpuLT_trig_single = cpuLT_trig4, cpuLT_trig_err_Bi_single = cpuLT_trig4_err_Bi, cpuLT_trig_err_Bay_single = cpuLT_trig4_err_Bay, tLT_trig_single = tLT_trig4, tLT_trig_err_Bi_single = tLT_trig4_err_Bi, tLT_trig_err_Bay_single = tLT_trig4_err_Bay, Ps_factor_single = Ps4_factor, total_trig_scaler_bcm_cut_single = total_trig4_scaler_bcm_cut, total_trig_accp_bcm_cut_single = total_trig4_accp_bcm_cut; }


  if(trig_type_coin=="T5") { trig_rate_coin = TRIG5scalerRate_bcm_cut, cpuLT_trig_coin = cpuLT_trig5, cpuLT_trig_err_Bi_coin = cpuLT_trig5_err_Bi, cpuLT_trig_err_Bay_coin = cpuLT_trig5_err_Bay, tLT_trig_coin = tLT_trig5, tLT_trig_err_Bi_coin = tLT_trig5_err_Bi, tLT_trig_err_Bay_coin = tLT_trig5_err_Bay, Ps_factor_coin = Ps5_factor, total_trig_scaler_bcm_cut_coin = total_trig5_scaler_bcm_cut, total_trig_accp_bcm_cut_coin = total_trig5_accp_bcm_cut; }
  if(trig_type_coin=="T6") { trig_rate_coin = TRIG6scalerRate_bcm_cut, cpuLT_trig_coin = cpuLT_trig6, cpuLT_trig_err_Bi_coin = cpuLT_trig6_err_Bi, cpuLT_trig_err_Bay_coin = cpuLT_trig6_err_Bay, tLT_trig_coin = tLT_trig6, tLT_trig_err_Bi_coin = tLT_trig6_err_Bi, tLT_trig_err_Bay_coin = tLT_trig6_err_Bay, Ps_factor_coin = Ps6_factor, total_trig_scaler_bcm_cut_coin = total_trig6_scaler_bcm_cut, total_trig_accp_bcm_cut_coin = total_trig6_accp_bcm_cut; }
  
 

 //Calculate HMS Tracking Efficiency                                                                                                                 
  hTrkEff = h_did / h_should;                                                                                                                  
  hTrkEff_err = sqrt(h_should-h_did) / h_should;
  
  //Calculate SHMS Tracking Efficiency                                                                                                
  pTrkEff = p_did / p_should; 
  pTrkEff_err = sqrt(p_should-p_did) / p_should; 

  std::cout << Form("HMS Should have tracked events = %.4f", h_should)<< std::endl;  
std::cout << Form("HMS did tracking of events = %.4f", h_did )<< std::endl;  
                                                           

std::cout << Form("HMS Tracking efficiency = %.4f", hTrkEff )<< std::endl;  

std::cout << Form("SHMS Should have tracked events = %.4f", p_should)<< std::endl;  
std::cout << Form("SHMS did tracking of events = %.4f", p_did )<< std::endl;  
                                                           

std::cout << Form("SHMS Tracking efficiency = %.4f", pTrkEff )<< std::endl;  

// Scale Histograms by Charge
Double_t scale_factor = (total_charge_bcm4a_cut*hTrkEff*pTrkEff*tLT_trig6);
  H_W->Scale(1/scale_factor);
  H_Q2->Scale(1/scale_factor);
  H_Em->Scale(1/scale_factor);
  H_h_delta->Scale(1/scale_factor);
  H_e_delta->Scale(1/scale_factor);
  H_e_xptar->Scale(1/scale_factor);
  H_e_yptar->Scale(1/scale_factor);
  H_h_xptar->Scale(1/scale_factor);
  H_h_yptar->Scale(1/scale_factor);
  H_Pmz->Scale(1/scale_factor);
  H_Pmy->Scale(1/scale_factor);
  H_Pmx->Scale(1/scale_factor);
  H_Pm->Scale(1/scale_factor);
  H_x_bj->Scale(1/scale_factor);
  H_xangle_data->Scale(1/scale_factor);
  H_theta_e->Scale(1/scale_factor);
  H_theta_p->Scale(1/scale_factor);
  H_ztarDiff->Scale(1/scale_factor);
  H_Em_vs_Pm->Scale(1/scale_factor);
  H_Pmea->Scale(1/scale_factor);
  H_Pcalc->Scale(1/scale_factor);
  H_deltaPdata->Scale(1/scale_factor);
  H_deltaPdata_vs_xfp->Scale(1/scale_factor);
  H_deltaPdata_vs_xpfp->Scale(1/scale_factor);
  H_deltaPdata_vs_yfp->Scale(1/scale_factor);
  H_deltaPdata_vs_ypfp->Scale(1/scale_factor);

std::cout << "All histograms scaled = "  << std::endl;
   // Write histograms to the output file

        TString data_OutputFileName = "deut_20846_-4_output_normalized.root";
        TFile *outROOT = new TFile(data_OutputFileName, "RECREATE");
        outROOT->cd();
        kin_HList->Write();
        outROOT->Close();
      

   {
        std::cerr << "Invalid analysis type. Use 'data' or 'simc'." << std::endl;
    }

    // Clean up
    delete inFile;
    delete kin_HList;


 }



}

