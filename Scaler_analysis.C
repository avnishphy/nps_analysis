void CreateScalerHMSData() {
    // Open the input ROOT file
    TFile *inFile = new TFile("/lustre24/expphy/cache/hallc/c-nps/analysis/online/replays/nps_hms_coin_4726_3_1_-1.root", "READ");

    // Access the relevant tree(s) from the file
    TTree *scaler_tree = (TTree*)inFile->Get("TSH");

    // Get the number of entries in the scaler tree
    Long64_t scal_entries = scaler_tree->GetEntries();
    
    // Dynamically allocate arrays for event flags and event numbers
    Int_t *evt_flag_bcm = new Int_t[scal_entries]; // store 0 or 1, to determine which scaler read passed cut
    Int_t *scal_evt_num = new Int_t[scal_entries]; // store event associated with scaler read

    // (The rest of your tree setup and histogram creation code goes here...)
     // Variables to hold data
    Double_t Scal_BCM1_charge, Scal_BCM1_current, S2X_scaler, TRIG5_scaler;
    Double_t Scal_BCM2_charge, Scal_BCM2_current, S2Y_scaler, TRIG6_scaler;
    Double_t Scal_BCM4A_charge, Scal_BCM4A_current, TRIG1_scaler, EDTM_scaler;
    Double_t Scal_BCM4B_charge, Scal_BCM4B_current, TRIG2_scaler, evNum;
    Double_t Scal_BCM4C_charge, Scal_BCM4C_current, TRIG4_scaler, TRIG3_scaler;
    Double_t Scal_time, S1X_scaler, S1Y_scaler, total_charge;
    TString bcm_type = "BCM4A"; 
//Double_t evNum;
  //Define Counter Quantities To Store Previous Reads
  Double_t prev_time = 0.;
  Double_t prev_charge_bcm1 = 0.;
  Double_t prev_charge_bcm2 = 0.;
  Double_t prev_charge_bcm4a = 0.;
  Double_t prev_charge_bcm4b = 0.;
  Double_t prev_charge_bcm4c = 0.;
  Double_t prev_s1x_scaler = 0;
  Double_t prev_s1y_scaler = 0;
  Double_t prev_s2x_scaler = 0;
  Double_t prev_s2y_scaler = 0;
  Double_t prev_trig1_scaler = 0;
  Double_t prev_trig2_scaler = 0;
  Double_t prev_trig3_scaler = 0;
  Double_t prev_trig4_scaler = 0;
  Double_t prev_trig5_scaler = 0;
  Double_t prev_trig6_scaler = 0;
  Double_t prev_edtm_scaler = 0;
  Double_t bcm_thrs = 5;
  //Double_t Scal_time;
  //Define Counter Quantities To Store Accumulated Reads
  Double_t total_time = 0.;
  Double_t total_charge_bcm = 0.; // placeholder for arbitrary BCM type (determined by user)
  Double_t total_charge_bcm1 = 0.;
  Double_t total_charge_bcm2 = 0.;
  Double_t total_charge_bcm4a = 0.;
  Double_t total_charge_bcm4b = 0.;
  Double_t total_charge_bcm4c = 0.;
  Double_t total_s1x_scaler = 0;
  Double_t total_s1y_scaler = 0;
  Double_t total_s2x_scaler = 0;
  Double_t total_s2y_scaler = 0;
  Double_t total_trig1_scaler = 0;
  Double_t total_trig2_scaler = 0;
  Double_t total_trig3_scaler = 0;
  Double_t total_trig4_scaler = 0;
  Double_t total_trig5_scaler = 0;
  Double_t total_trig6_scaler = 0;
  Double_t total_edtm_scaler = 0;
  //Store Accumulated Reads if they passed BCM Current Cut
  Double_t total_time_bcm_cut = 0.;
  Double_t total_charge_bcm_cut = 0.;  // placeholder for arbitrary BCM type (determined by user)
  Double_t total_charge_bcm1_cut = 0.;
  Double_t total_charge_bcm2_cut = 0.;
  Double_t total_charge_bcm4a_cut = 0.;
  Double_t total_charge_bcm4b_cut = 0.;
  Double_t total_charge_bcm4c_cut = 0.;
  Double_t total_s1x_scaler_bcm_cut = 0;
  Double_t total_s1y_scaler_bcm_cut = 0;
  Double_t total_s2x_scaler_bcm_cut = 0;
  Double_t total_s2y_scaler_bcm_cut = 0;
  Double_t total_trig1_scaler_bcm_cut = 0;
  Double_t total_trig2_scaler_bcm_cut = 0;
  Double_t total_trig3_scaler_bcm_cut = 0;
  Double_t total_trig4_scaler_bcm_cut = 0;
  Double_t total_trig5_scaler_bcm_cut = 0;
  Double_t total_trig6_scaler_bcm_cut = 0;
  Double_t total_edtm_scaler_bcm_cut = 0;
  Double_t Scal_BCM_charge;  // generic placeholder for BCM charge (depend on user input)
  Double_t Scal_BCM_current;
  Double_t avg_current_bcm_cut;
//total_charge = NULL;
  //Store Scaler Rates if current cut passed
  Double_t S1XscalerRate_bcm_cut;
  Double_t S1YscalerRate_bcm_cut;
  Double_t S2XscalerRate_bcm_cut;
  Double_t S2YscalerRate_bcm_cut;

  Double_t TRIG1scalerRate_bcm_cut;
  Double_t TRIG2scalerRate_bcm_cut;
  Double_t TRIG3scalerRate_bcm_cut;
  Double_t TRIG4scalerRate_bcm_cut;
  Double_t TRIG5scalerRate_bcm_cut;
  Double_t TRIG6scalerRate_bcm_cut;
  Double_t EDTMscalerRate_bcm_cut;
//delete h_total_charge; h_total_charge = NULL;
 //Int_t *evt_flag_bcm;   //flag (0 or 1) to determine whether the scaler read passed the cut
  //Int_t *scal_evt_num;   //store data event number associated with scaler read 
  Double_t set_current;  //Set current for each run. For now, take the current ->  maximum bin content of bcm current histogram 
  Int_t scal_read = 0; 

    // Set branch addresses to fetch data from the source tree
    scaler_tree->SetBranchAddress("evNumber", &evNum);
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

    // Create histograms
     TH1F *h_evNum = new TH1F("h_evNum", "evNum", 100, 410, 455);
    TH1F *h_Scal_BCM1_charge = new TH1F("h_Scal_BCM1_charge", "BCM1 Charge", 100, 0, 1000);
    TH1F *h_Scal_BCM1_current = new TH1F("h_Scal_BCM1_current", "BCM1 Current", 100, 26, 32);
    TH1F *h_Scal_BCM2_charge = new TH1F("h_Scal_BCM2_charge", "BCM2 Charge", 100, 0, 1000);
    TH1F *h_Scal_BCM2_current = new TH1F("h_Scal_BCM2_current", "BCM2 Current", 100, 26, 32);
    TH1F *h_Scal_BCM4A_charge = new TH1F("h_Scal_BCM4A_charge", "BCM4A Charge", 100, 0, 1000);
    TH1F *h_Scal_BCM4A_current = new TH1F("h_Scal_BCM4A_current", "BCM4A Current", 100, 26, 32);
    TH1F *h_Scal_BCM4B_charge = new TH1F("h_Scal_BCM4B_charge", "BCM4B Charge", 100, 0, 1000);
    TH1F *h_Scal_BCM4B_current = new TH1F("h_Scal_BCM4B_current", "BCM4B Current", 100, 24.5, 26);
    TH1F *h_Scal_BCM4C_charge = new TH1F("h_Scal_BCM4C_charge", "BCM4C Charge", 100, 0, 1000);
    TH1F *h_Scal_BCM4C_current = new TH1F("h_Scal_BCM4C_current", "BCM4C Current", 100, 26, 32);
    TH1F *h_Scal_time = new TH1F("h_Scal_time", "Scaler Time", 100, -125, 2000);
    TH1F *h_S1X_scaler = new TH1F("h_S1X_scaler", "S1X Scaler", 100, 410, 455);
    TH1F *h_S1Y_scaler = new TH1F("h_S1Y_scaler", "S1Y Scaler", 100, 410, 455);
    TH1F *h_S2X_scaler = new TH1F("h_S2X_scaler", "S2X Scaler", 100, 410, 455);
    TH1F *h_S2Y_scaler = new TH1F("h_S2Y_scaler", "S2Y Scaler", 100, 410, 455);
    TH1F *h_TRIG1_scaler = new TH1F("h_TRIG1_scaler", "TRIG1 Scaler", 100, 0, 20);
    TH1F *h_TRIG2_scaler = new TH1F("h_TRIG2_scaler", "TRIG2 Scaler", 100, 0, 20);
    TH1F *h_TRIG3_scaler = new TH1F("h_TRIG3_scaler", "TRIG3 Scaler", 100, 410, 455);
    TH1F *h_TRIG4_scaler = new TH1F("h_TRIG4_scaler", "TRIG4 Scaler", 100, 410, 455);
    TH1F *h_TRIG5_scaler = new TH1F("h_TRIG5_scaler", "TRIG5 Scaler", 100, 410, 455);
    TH1F *h_TRIG6_scaler = new TH1F("h_TRIG6_scaler", "TRIG6 Scaler", 100, 410, 455);
    TH1F *h_EDTM_scaler = new TH1F("h_EDTM_scaler", "EDTM Scaler", 100, 410, 455);
TH1F *h_total_charge = new TH1F ("total charge", "total Charge", 100, 47, 47.5);
TH1F *h_average_current = new TH1F("average current", "average current", 100, 27, 29);

//Create TLists to store categorical histograms (this is helpful to for organizing histograms later on, when you have many)
  TList *kin_HList  = new TList(); // kinematics histograms list
    
   kin_HList->Add(h_Scal_BCM1_charge);
   kin_HList->Add(h_Scal_BCM1_current);
   kin_HList->Add(h_Scal_BCM2_charge);
   kin_HList->Add(h_Scal_BCM2_current);
   kin_HList->Add(h_Scal_BCM4A_charge);
   kin_HList->Add(h_Scal_BCM4A_current);
   kin_HList->Add(h_Scal_BCM4B_charge);
   kin_HList->Add(h_Scal_BCM4B_current);
   kin_HList->Add(h_Scal_BCM4C_charge);
   kin_HList->Add(h_Scal_BCM4C_current);
   kin_HList->Add(h_Scal_time);
   kin_HList->Add(h_S1X_scaler);
   kin_HList->Add(h_S1Y_scaler);
    kin_HList->Add(h_S2X_scaler);
   kin_HList->Add(h_S2Y_scaler);
    kin_HList->Add(h_TRIG1_scaler);
   kin_HList->Add(h_TRIG1_scaler);
   kin_HList->Add(h_TRIG1_scaler);
    kin_HList->Add(h_TRIG2_scaler);
    kin_HList->Add(h_TRIG3_scaler);
   kin_HList->Add(h_TRIG4_scaler);
   kin_HList->Add(h_TRIG5_scaler);
   kin_HList->Add(h_TRIG6_scaler);
    kin_HList->Add(h_EDTM_scaler);
 kin_HList->Add(h_evNum);
kin_HList->Add(h_total_charge);
kin_HList->Add(h_average_current);
    
    
    
    for (int i = 0; i < scal_entries; i++) {
        scaler_tree->GetEntry(i);
        evt_flag_bcm[i] = 0;  // Default flag
        
        // Store the event number
        scal_evt_num[i] = evNum;
        
        // Your logic to determine if the flag should be set, etc.
         // Determine which bcm current to cut on (based on user input)
    //TString bcm_type = "BCM4A";  // Set this to the BCM you want to analyze
  
      if(bcm_type=="BCM4A"){
	Scal_BCM_current = Scal_BCM4A_current;
      }
      

      //Store Cumulative Quantities
      total_charge_bcm4a = Scal_BCM4A_charge;
      total_charge_bcm4b = Scal_BCM4B_charge;
      total_charge_bcm4c = Scal_BCM4C_charge;
      total_time = Scal_time;
      
std::cout << "Scal_time = " << Scal_time << std::endl;
// Print the running total of the charge
//std::cout << "Total BCM4A Charge so far = " << Scal_BCM4A_charge << std::endl;
//std::cout << "Total BCM4C Charge so far = " << Scal_BCM4C_charge << std::endl;
std::cout << "Total BCM4A Current is = " << Scal_BCM4A_current << std::endl;
std::cout << "Total BCM Current so far = " << Scal_BCM_current << std::endl;
   
      if(Scal_BCM_current > bcm_thrs)
	{
   
	  std::cout << "Total BCM Current so far = " << Scal_BCM_current << std::endl;
    
	  //Turn Event Flag ON, if beam current is within threshold
	  evt_flag_bcm[i] = 1;
	  
	  //Store Quantities that Passed the Current Threshold
	  total_time_bcm_cut = total_time_bcm_cut + (Scal_time - prev_time);
	  total_charge_bcm1_cut = total_charge_bcm1_cut + (Scal_BCM1_charge - prev_charge_bcm1); 
    total_charge_bcm2_cut = total_charge_bcm2_cut + (Scal_BCM2_charge - prev_charge_bcm2);  
	  total_charge_bcm4a_cut = total_charge_bcm4a_cut + (Scal_BCM4A_charge - prev_charge_bcm4a);  
	 

	} //End BCM Current Cut
std::cout << "Total BCM time cut so far = " << total_charge_bcm1_cut << std::endl;
    
     //Previous Scaler Reads (Necessary to Take Average between S-1 and S scaler reads, to get values in between)
      prev_time = Scal_time;
      prev_charge_bcm1 = Scal_BCM1_charge;
      prev_charge_bcm2 = Scal_BCM2_charge;
      prev_charge_bcm4a = Scal_BCM4A_charge;
      prev_charge_bcm4b = Scal_BCM4B_charge;
      
      
    }
      // Set generic bcm info to be used in charge normalization based on user input
    if(bcm_type=="BCM4A"){
    total_charge_bcm     = total_charge_bcm4a;
    total_charge_bcm_cut = total_charge_bcm4a_cut;
  }
  
  std::cout << "Total BCM4C Current cut = " << total_charge_bcm4a_cut << std::endl;
  
   //Calculate Average BCM Current                                                  
  avg_current_bcm_cut = total_charge_bcm_cut / total_time_bcm_cut; //uA                              
  
  //Convert charge from uC to mC                                   
  total_charge_bcm_cut = total_charge_bcm_cut / 1000.; 
  total_charge_bcm1_cut  = total_charge_bcm1_cut / 1000.; 
  total_charge_bcm2_cut  = total_charge_bcm2_cut / 1000.; 
  total_charge_bcm4a_cut = total_charge_bcm4a_cut / 1000.; 
  total_charge_bcm4b_cut = total_charge_bcm4b_cut / 1000.; 
  total_charge_bcm4c_cut = total_charge_bcm4c_cut / 1000.; 

  //h_total_charge->SetBinContent(46 ,total_charge_bcm4a_cut);
  

  
std::cout << "Total Charge micro coulomb into milli coulomb = " << total_charge_bcm_cut << std::endl; 
std::cout << "average current for bcm4a is = " << avg_current_bcm_cut << std::endl;
  // Loop over the entries in the source tree and fill the histograms
    Long64_t nentries = scaler_tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        scaler_tree->GetEntry(i);
        h_Scal_BCM1_charge->Fill(Scal_BCM1_charge);
        h_Scal_BCM1_current->Fill(Scal_BCM1_current);
        h_Scal_BCM2_charge->Fill(Scal_BCM2_charge);
        h_Scal_BCM2_current->Fill(Scal_BCM2_current);
        h_Scal_BCM4A_charge->Fill(Scal_BCM4A_charge);
        h_Scal_BCM4A_current->Fill(Scal_BCM4A_current);
        h_Scal_BCM4B_charge->Fill(Scal_BCM4B_charge);
        h_Scal_BCM4B_current->Fill(Scal_BCM4B_current);
        h_Scal_BCM4C_charge->Fill(Scal_BCM4C_charge);
        h_Scal_BCM4C_current->Fill(Scal_BCM4C_current);
        h_Scal_time->Fill(Scal_time);
        h_S1X_scaler->Fill(S1X_scaler);
        h_S1Y_scaler->Fill(S1Y_scaler);
        h_S2X_scaler->Fill(S2X_scaler);
        h_S2Y_scaler->Fill(S2Y_scaler);
        h_TRIG1_scaler->Fill(TRIG1_scaler);
        h_TRIG2_scaler->Fill(TRIG2_scaler);
        h_TRIG3_scaler->Fill(TRIG3_scaler);
        h_TRIG4_scaler->Fill(TRIG4_scaler);
        h_TRIG5_scaler->Fill(TRIG5_scaler);
        h_TRIG6_scaler->Fill(TRIG6_scaler);
        h_EDTM_scaler->Fill(EDTM_scaler);
         h_evNum->Fill(evNum);
         h_total_charge ->Fill(total_charge_bcm_cut);
         h_average_current-> Fill(avg_current_bcm_cut);
    }
 
    // (The rest of your code, including saving histograms, goes here...)
     // Create a new output file
    TFile *outFile = new TFile("/u/group/nps/singhav/nps_analysis/scaler_TSH_BCM4A.root", "RECREATE");

   // Write histograms to the output file
   kin_HList->Write();
    // Close the output file
    outFile->Close();

    // Clean up dynamically allocated memory
    delete[] evt_flag_bcm;
    delete[] scal_evt_num;

    inFile->Close();
}
