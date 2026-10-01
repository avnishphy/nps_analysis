#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <TChain.h>
#include <TH1F.h>

void Analyzer_new(Int_t Output_Flag, int RunNo, Double_t* result);
Double_t ReadTILT(Int_t RunNo, Int_t PS_Flag);
Double_t ReadTriggerRate(Int_t RunNo, Int_t PS_Flag);
void DT_Analyzer_new(Int_t batch_job = 1);

void DT_Analyzer_new(Int_t batch_job = 1){
    Double_t result[22];
    Int_t RunNo, Which_Trigger;
    Double_t TI_LT;
    // Analyzer(1, 3993, result);
    // Analyzer(0, RunNo, result);
    // Int_t Which_Trigger = int(result[11]);
    // Double_t TI_LT = ReadTILT(RunNo,Which_Trigger);

    // cout<<endl;
    // cout<<"EDTM Live Time: "<<result[0]<<" = "<<result[1]<<" / "<<result[2]<<"\t"<<endl;
    // cout<<"EDTM event count after low beam current cut: ";
    // cout<<result[3]<<" - "<<result[5]<<" = "<<result[1]<<endl;
    // cout<<"EDTM scaler count after low beam current cut: ";
    // cout<<result[4]<<" - "<<result[6]<<" + "<<result[7]<<" = "<<result[2]<<endl;
    // cout<<"Total Diff: "<<result[4]-result[3]<<endl;
    // cout<<"EDTM Scaler correction: "<<result[7]<<endl;
    // cout<<"Difference after 2uA cut: "<<result[2] - result[1]<<endl;
    // cout<<"Time Low Curr * EDTM_Rate: "<<result[9]<<" * "<<result[8]<<" = "<< result[9]*result[8]<<endl;
    // cout<<"Which trigger we're using: "<<result[11]<<endl;
    // cout<<"Trigger rate: "<<result[10]<<endl;
    // cout<<endl;
    // return;

    // Components in the array "result"
    // 0. EDTM LT after 2uA cut
    // 1. EDTM TDC events after 2uA cut
    // 2. EDTM Scaler count after 2uA cut
    // 3. Total EDTM TDC events
    // 4. Total EDTM Scaler count
    // 5. EDTM TDC events   (I<2uA)
    // 6. EDTM Scaler count (I<2uA)
    // 7. Missing EDTM Scaler count
    // 8. EDTM Rate
    // 9. Low beam current time(s)
    // 10. Trigger Rate (Hz)
    // 11. Which trigger we're using
    // 12. Minimum time difference between two events (ns)
    // 13. lambda[us^{-1}]
    // 14. lambda_err[us^{-1}]
    // 15. LT_Event_Tdiff_fit
    // 16. LT_Event_Tdiff_fit_err
    // 17. Total events
    // 18. Total time in second
    // 19. LT -- Good time percentage

    // cout<<result[0]<<"\t"<<result[1]/result[2]<<endl;
    // cout<<result[1]<<"\t"<<result[3]-result[5]<<endl;
    // cout<<result[2]<<"\t"<<result[4]-result[6]+result[7]<<endl;

    TFile *DT_output = new TFile(Form("./output/DeadTime_%d_new.root",batch_job), "RECREATE");
    TTree *tree = new TTree("DeadTime", "Dead_Time_Information");
    tree->Branch("RunNo", &RunNo, "RunNo/I");
    tree->Branch("Success_flag",&result[21],"Success_flag/D");
    tree->Branch("EDTM_LT", &result[0], "EDTM_LT/D");
    tree->Branch("EDTM_TDC", &result[1], "EDTM_TDC/D");
    tree->Branch("EDTM_Scaler", &result[2], "EDTM_Scaler/D");
    tree->Branch("EDTM_TDC_Total", &result[3], "EDTM_TDC_Total/D");
    tree->Branch("EDTM_Scaler_Total", &result[4], "EDTM_Scaler_Total/D");
    tree->Branch("EDTM_TDC_LowCurr", &result[5], "EDTM_TDC_LowCurr/D");
    tree->Branch("EDTM_Scaler_LowCurr", &result[6], "EDTM_Scaler_LowCurr/D");
    tree->Branch("EDTM_Scaler_Correction", &result[7], "EDTM_Scaler_Correction/D");
    tree->Branch("EDTM_Rate", &result[8], "EDTM_Rate/D");
    tree->Branch("Time_LowCurr", &result[9], "Time_LowCurr/D");
    tree->Branch("Trigger_Rate", &result[10], "Trigger_Rate/D");
    tree->Branch("Which_Trigger", &Which_Trigger, "Which_Trigger/I");
    tree->Branch("PS", &result[20], "PS/D");
    tree->Branch("TI_LT", &TI_LT, "TI_LT/D");
    tree->Branch("Event_Tdiff_Min",&result[12],"Event_Tdiff_Min/D");
    tree->Branch("Lambda",&result[13],"Lambda/D");
    tree->Branch("Lambda_err",&result[14],"Lambda_err/D");
    tree->Branch("LT_Event_Tdiff_fit",&result[15],"LT_Event_Tdiff_fit/D");
    tree->Branch("LT_Event_Tdiff_fit_err",&result[16],"LT_Event_Tdiff_fit_err/D");
    tree->Branch("Total_Events",&result[17],"Total_Events/D");
    tree->Branch("Total_Time",&result[18],"Total_Time/D");
    tree->Branch("LT_Time_Percent",&result[19],"LT_Time_Percent/D");

    std::string file_path = "run_list_sorted.txt";
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cerr << "Can't open the run_list file:" << file_path << std::endl;
    }

    Int_t NN = 0;

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int run_number;
        std::string type, cls;

        if (iss >> run_number >> type >> cls) {

            if(type == "production"){
                NN++;
                if(NN>100*(batch_job-1) && NN<=100*batch_job){ // 3571 production runs
                    RunNo = run_number;

                    // RunNo = 1845;
                    // RunNo = 1729;

                    // RunNo = 3448;
                    // RunNo = 3566;
                    // RunNo = 5803;
                    // RunNo = 5889;

                    Analyzer_new(0, RunNo, result);
                    Which_Trigger = int(result[11]);
                    TI_LT = ReadTILT(RunNo,Which_Trigger);
                    
                    tree->Fill();

                    memset(result,-1.,sizeof(result));
                    Which_Trigger = -1;
                    TI_LT = -1;

                    // break;
                }
            }
        } else {
            std::cerr << "Wrong Type: " << line << std::endl;
        }
    }

    infile.close();

    DT_output->cd();
    tree->Write();
    DT_output->Close();
    delete DT_output;

    // cout<<endl;
    // cout<<"EDTM Live Time: "<<result[0]<<" = "<<result[1]<<" / "<<result[2]<<"\t"<<endl;
    // cout<<"EDTM event count after low beam current cut: ";
    // cout<<result[3]<<" - "<<result[5]<<" = "<<result[1]<<endl;
    // cout<<"EDTM scaler count after low beam current cut: ";
    // cout<<result[4]<<" - "<<result[6]<<" + "<<result[7]<<" = "<<result[2]<<endl;
    // cout<<"Total Diff: "<<result[4]-result[3]<<endl;
    // cout<<"EDTM Scaler correction: "<<result[7]<<endl;
    // cout<<"Difference after 2uA cut: "<<result[2] - result[1]<<endl;
    // cout<<"Time Low Curr * EDTM_Rate: "<<result[9]<<" * "<<result[8]<<" = "<< result[9]*result[8]<<endl;
    // cout<<"Which trigger we're using: "<<result[11]<<endl;
    // cout<<"Trigger rate: "<<result[10]<<endl;
    // cout<<endl;

}

void Analyzer_new(Int_t Output_Flag, int RunNo, Double_t* result){
    // TString FilePath_PS = "/volatile/hallc/nps/nps-ana/REPORT_OUTPUT/REPORT_OUTPUT2/COIN/SKIM/";
    TString FilePath_PS = "./ReportFiles/"; // /work/hallc/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/SKIM
    TString Input_File_PS = FilePath_PS + Form("skim_NPS_HMS_report_%d_-1.report",RunNo);
    std::ifstream inputFile(Input_File_PS);
    if(!inputFile.is_open()){
        result[21] = -1.;
        return;
    } 
    else result[21] = 1;

    int Ps1_factor, Ps2_factor, Ps3_factor, Ps4_factor, Ps5_factor, Ps6_factor;
    Double_t NofCharge, Target_Type, BeamCurrent, EDTM_RATE;
    EDTM_RATE = 0.;
    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.find("Ps1_factor") != std::string::npos) {
            std::sscanf(line.c_str(), "Ps1_factor = %d", &Ps1_factor);
        } else if (line.find("Ps2_factor") != std::string::npos) {
            std::sscanf(line.c_str(), "Ps2_factor = %d", &Ps2_factor);
        } else if (line.find("Ps3_factor") != std::string::npos) {
            std::sscanf(line.c_str(), "Ps3_factor = %d", &Ps3_factor);
        } else if (line.find("Ps4_factor") != std::string::npos) {
            std::sscanf(line.c_str(), "Ps4_factor = %d", &Ps4_factor);
        } else if (line.find("Ps5_factor") != std::string::npos) {
            std::sscanf(line.c_str(), "Ps5_factor = %d", &Ps5_factor);
        } else if (line.find("Ps6_factor") != std::string::npos) {
            std::sscanf(line.c_str(), "Ps6_factor = %d", &Ps6_factor);
        }

        if (line.find("BCM4A Beam Cut Charge:") != std::string::npos) {
            std::sscanf(line.c_str(), "BCM4A Beam Cut Charge: %lf uC", &NofCharge);
        }

        if (line.find("BCM4A Beam Cut Current:") != std::string::npos)
        {
            std::sscanf(line.c_str(), "BCM4A Beam Cut Current: %lf uA", &BeamCurrent);
        }

        if (line.find("Target AMU") != std::string::npos)
        {
            std::sscanf(line.c_str(), "Target AMU  : %lf", &Target_Type);
        }

        if (line.find("EDTM Trigger Rate") != std::string::npos)
        {
            std::sscanf(line.c_str(), "EDTM Trigger Rate       : %lf kHz", &EDTM_RATE);
            break;
        }
 
    }
    EDTM_RATE = EDTM_RATE*1000;

    int ps[6];
    Double_t Ps[6];
    ps[0] = Ps1_factor;
    ps[1] = Ps2_factor;
    ps[2] = Ps3_factor;
    ps[3] = Ps4_factor;
    ps[4] = Ps5_factor;
    ps[5] = Ps6_factor;

    Double_t PSFACTOR;
    Double_t PS;
    for(int i=0;i<6;i++){
        if(ps[i]!=-1){
            PSFACTOR = (i+1)*1.;
            // cout<<ps[i]<<endl;
            PS = ps[i]*1.;
        }
        if(ps[i]==-1) Ps[i] = -1.;
        else if(ps[i]==1) Ps[i] = 0.;
        else Ps[i] = log2((ps[i]-1)*1.)+1.;
    }
    
    inputFile.close();
    
    Double_t TRIG_RATE = ReadTriggerRate(RunNo,PSFACTOR)*1000.;
    Double_t Trigger_Rate = 0;

    TChain	*ch1 = new TChain("T");
    TChain	*ch2 = new TChain("TSH");
    TString File_Path;
    // File_Path = "/cache/hallc/c-nps/analysis/online/replays/production/"; // pass0
    // File_Path = "/cache/hallc/c-nps/analysis/pass1/replays/skim/"; // pass1
    // TString FileName = File_Path + Form("nps_hms_skim_%d_1_-1.root",RunNo);
    // File_Path = "/cache/hallc/c-nps/analysis/online/replays/";
    // TString FileName = File_Path + Form("nps_hms_coin_%d_0_1_-1.root",RunNo);

    File_Path = "./";
    TString FileName = File_Path + Form("nps_hms_skim_%d_1_-1.root",RunNo);

    TFile* inputFile_root = new TFile(FileName, "READ");
    if (!inputFile_root || inputFile_root->IsZombie()) {
        result[21] = -1.;
        return;
    }
    else result[21] = 1;

    Double_t EDTM_TDC_Low, EDTM_TDC_High;
    TFile* file = new TFile(FileName, "READ");
    TTree* tree = (TTree*)file->Get("T");
    TH1F* histogram = new TH1F("histogram", "Branch to Histogram", 300,0,3000);
    tree->Draw("T.hms.hEDTM_tdcTimeRaw >> histogram");
    // TH1F* histogram = new TH1F("histogram", "EDTM TDC Time Raw",300,0,3000);
    // ch1->Draw("T.hms.hEDTM_tdcTimeRaw >> histogram");
    Int_t Main_Peak_Flag = 0;
    Int_t Main_Peak_Flag_last = 0;
    for (int i = 1; i <= 300; ++i) {
        double binContent = histogram->GetBinContent(i);
        double binLowEdge = histogram->GetBinLowEdge(i);
        double binWidth = histogram->GetBinWidth(i);
        double binUpEdge = binLowEdge + binWidth;
        if(binContent>10.){
            if(Main_Peak_Flag_last==0) EDTM_TDC_Low = binLowEdge;
            Main_Peak_Flag = 1;
        }
        else{
            if(Main_Peak_Flag_last==1) EDTM_TDC_High = binUpEdge;
            Main_Peak_Flag = 0;
        }
        Main_Peak_Flag_last = Main_Peak_Flag;
    }
    // cout<<EDTM_TDC_Low<<"\t"<<EDTM_TDC_High<<endl;
    file->Close();

    ch1->Add(FileName);
    ch2->Add(FileName);

    Long64_t nentries = ch1->GetEntries();
    Long64_t nentries_TSH = ch2->GetEntries();
    if(nentries_TSH<10 || nentries<100){
        result[21] = -1.;
        return;
    }

    ch1->SetBranchStatus("*",false);
    ch2->SetBranchStatus("*",false);

    ch1->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", true);
    ch1->SetBranchStatus("fEvtHdr.fEvtTime", true);
    // ch1->SetBranchStatus("H.EDTM_CP.scaler", true);
    // ch1->SetBranchStatus("H.1MHz.scaler", true);
    // ch1->SetBranchStatus("H.hTRIG3.scalerRate", true);
    // ch1->SetBranchStatus("H.hL1ACCP.scalerRate", true);

    ch2->SetBranchStatus("H.EDTM_CP.scaler", true);
    ch2->SetBranchStatus("H.1MHz.scaler", true);
    ch2->SetBranchStatus("H.BCM4A.scalerCurrent", true);
    ch2->SetBranchStatus("H.hL1ACCP.scalerRate", true);
    ch2->SetBranchStatus("evNumber", true);
    // ch2->SetBranchStatus("H.BCM4A.scaler", true);
    // ch2->SetBranchStatus("H.BCM4A.scalerRate", true);
    // ch2->SetBranchStatus("evcount", true);

    Double_t H_EDTM_CP_scaler, H_1MHz_scaler, H_BCM4A_scalerCurrent, H_hL1ACCP_scalerRate_TSH, evNumber;
    ch2->SetBranchAddress("H.EDTM_CP.scaler", &H_EDTM_CP_scaler);
    ch2->SetBranchAddress("H.1MHz.scaler", &H_1MHz_scaler);
    ch2->SetBranchAddress("H.BCM4A.scalerCurrent", &H_BCM4A_scalerCurrent);
    ch2->SetBranchAddress("H.hL1ACCP.scalerRate", &H_hL1ACCP_scalerRate_TSH);
    ch2->SetBranchAddress("evNumber", &evNumber);

    Double_t T_hms_hEDTM_tdcTimeRaw;
    ULong64_t fEvtHdr_fEvtTime;
    ch1->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &T_hms_hEDTM_tdcTimeRaw);
    ch1->SetBranchAddress("fEvtHdr.fEvtTime", &fEvtHdr_fEvtTime);

    // Double_t H_BCM4A_scaler, H_BCM4A_scalerRate, evcount;
    // ch2->SetBranchAddress("H.BCM4A.scaler", &H_BCM4A_scaler);
    // ch2->SetBranchAddress("H.BCM4A.scalerRate", &H_BCM4A_scalerRate);
    // ch2->SetBranchAddress("evcount", &evcount);

    // Double_t H_EDTM_CP_scaler_Ttree, H_1MHz_scaler_T, H_hTRIG3_scalerRate, H_hL1ACCP_scalerRate;
    // ch1->SetBranchAddress("H.EDTM_CP.scaler", &H_EDTM_CP_scaler_Ttree);     // Full replay beanches
    // ch1->SetBranchAddress("H.1MHz.scaler", &H_1MHz_scaler_T);               // Full replay beanches
    // ch1->SetBranchAddress("H.hTRIG3.scalerRate", &H_hTRIG3_scalerRate);     // Full replay beanches
    // ch1->SetBranchAddress("H.hL1ACCP.scalerRate", &H_hL1ACCP_scalerRate);   // Full replay beanches

    // ULong64_t Event_T_00;
    // ch2->GetEntry(0);
    // H_1MHz_scaler_0 = H_1MHz_scaler;
    // cout<<evNumber<<"\t"<<H_EDTM_CP_scaler<<endl;
    // ch1->GetEntry(0);
    // Event_T_00 = fEvtHdr_fEvtTime;
    // cout<<0<<"\t"<<H_EDTM_CP_scaler_Ttree<<endl;

    // ch2->GetEntry(1);
    // cout<<evNumber<<"\t"<<H_EDTM_CP_scaler<<"\t"<<(H_1MHz_scaler-H_1MHz_scaler_0)*1./1e6<<endl;
    // ch1->GetEntry(int(evNumber-1));
    // cout<<evNumber-1<<"\t"<<H_EDTM_CP_scaler_Ttree<<"\t"<<(fEvtHdr_fEvtTime-Event_T_00)*1./250e6<<endl;

    // return;

    Double_t EDTM_LT, EDTM_last, EDTM_now, EDTM_Count, EDTM_Event, EDTM_Event_LowCurr;
    Int_t LowCurr_Flag = 0;
    Int_t NofPeriod_LowCurr = 0;
    Double_t last_event = 1.;
    std::vector<double> start_event, end_event;
    std::vector<int> start_EDTM_scaler, end_EDTM_scaler;
    std::vector<double> start_event_time_TSH, end_event_time_TSH;
    std::vector<double> TSH_EventTime;
    std::vector<int> TSH_EventNo, TSH_EDTM_Scaler, TSH_hL1ACCP_ScalerRate;
    EDTM_Count = 0.;

    ch2->GetEntry(0);
    EDTM_last = H_EDTM_CP_scaler;
    Double_t H_1MHz_scaler_0 = H_1MHz_scaler;
    Double_t H_1MHz_scaler_last = H_1MHz_scaler;
    Double_t H_EDTM_CP_scaler_0 = H_EDTM_CP_scaler;
    Int_t NofTRIGRATE_Points = 0;
    Long64_t nentries2 = ch2->GetEntries();
    for(Long64_t i=0; i<nentries2; i++){
        ch2->GetEntry(i);
        if((H_1MHz_scaler-H_1MHz_scaler_last)*1./1e6>1.8||i<10){
            TSH_EventNo.push_back(int(evNumber-1));
            TSH_EventTime.push_back((H_1MHz_scaler-H_1MHz_scaler_0)*1./1e6);
            TSH_EDTM_Scaler.push_back(int(H_EDTM_CP_scaler-H_EDTM_CP_scaler_0));
            TSH_hL1ACCP_ScalerRate.push_back(H_hL1ACCP_scalerRate_TSH);
            if(fabs(H_hL1ACCP_scalerRate_TSH-TRIG_RATE)<TRIG_RATE*0.25){
                Trigger_Rate += H_hL1ACCP_scalerRate_TSH;
                NofTRIGRATE_Points++;
            }
        }
        
        EDTM_now = H_EDTM_CP_scaler;
        if(EDTM_now<10.) continue;
        if(H_BCM4A_scalerCurrent>2.){
            if(LowCurr_Flag == 1){
                end_event.push_back(last_event);
                end_event_time_TSH.push_back((H_1MHz_scaler_last-H_1MHz_scaler_0)*1./1e6);
                end_EDTM_scaler.push_back(int(EDTM_last-H_EDTM_CP_scaler_0));
            }
            LowCurr_Flag = 0;
            last_event = evNumber-1;
            H_1MHz_scaler_last = H_1MHz_scaler;
        }
        else{
            if(LowCurr_Flag == 0){
                start_event.push_back(last_event);
                start_event_time_TSH.push_back((H_1MHz_scaler-H_1MHz_scaler_0)*1./1e6);
                start_EDTM_scaler.push_back(int(H_EDTM_CP_scaler-H_EDTM_CP_scaler_0));
                NofPeriod_LowCurr++;
            }
            LowCurr_Flag = 1;
            last_event = evNumber-1;
            H_1MHz_scaler_last = H_1MHz_scaler;
        }
        EDTM_last = EDTM_now;
        H_1MHz_scaler_last = H_1MHz_scaler;
    }
    if(LowCurr_Flag == 1) end_event.push_back(evNumber-1);

    Trigger_Rate = Trigger_Rate/(NofTRIGRATE_Points*1.);
    Int_t TSH_Length = TSH_EventNo.size();

    for(int i=0;i<NofPeriod_LowCurr;i++){
        if(start_event[i]>TSH_EventNo[TSH_Length-1]) break;
        if(end_event[i]>TSH_EventNo[TSH_Length-1]){
            for(int k=0;k<TSH_Length;k++){
                if(TSH_EventNo[TSH_Length-1-k]<start_event[i]){
                    TSH_Length = TSH_Length - k;
                    break;
                }
            }
            break;
        }
    }

    // cout<<ch1->GetEntries()<<endl;
    // cout<<evNumber<<endl;
    // cout<<end_event[3]<<"\t"<<start_event[3]<<endl;
    // cout<<end_EDTM_scaler[3]<<"\t"<<start_EDTM_scaler[3]<<endl;
    // cout<<endl;
    // cout<<end_event[3]<<"\t"<<TSH_EventNo[TSH_Length-1]<<endl;
    // return;

    // cout<<TSH_EventNo[TSH_Length-4]<<"\t"<<TSH_EventNo[TSH_Length-3]<<"\t"<<TSH_EventNo[TSH_Length-2]<<"\t"<<TSH_EventNo[TSH_Length-1]<<endl;
    // cout<<TSH_EventTime[TSH_Length-4]<<"\t"<<TSH_EventTime[TSH_Length-3]<<"\t"<<TSH_EventTime[TSH_Length-2]<<"\t"<<TSH_EventTime[TSH_Length-1]<<endl;
    // cout<<TSH_EDTM_Scaler[TSH_Length-4]<<"\t"<<TSH_EDTM_Scaler[TSH_Length-3]<<"\t"<<TSH_EDTM_Scaler[TSH_Length-2]<<"\t"<<TSH_EDTM_Scaler[TSH_Length-1]<<endl;
    // return;

    Double_t EDTM_Scaler_correction = 0.;

    ch1->GetEntry(0);

    Int_t N_Period_LowCurr = 0;
    EDTM_Event = 0.;
    EDTM_Event_LowCurr = 0.;
    Int_t EDTM_Event_Nocut = 0;
    Int_t EDTM_Event_Nocut_last = 0;
    ULong64_t Event_T_0, Event_T_last, Event_T_last_event;
    Int_t N_TSH_Event = 0;
    Int_t Total_Events = 0;
    Double_t Total_Time = 0.;

    TFile* f1;
    if(Output_Flag==1) f1 = new TFile("./test.root","recreate");
    if(Output_Flag==2) f1 = new TFile("./test.root","update");
    TH1F* h_EDTM_tdcTimeRaw = new TH1F(Form("h_EDTM_tdcTimeRaw_%d",RunNo),Form("h_EDTM_tdcTimeRaw_%d;Raw tdc Time;Counts",RunNo),300,0,3000);
    TH1F* h_event_time_diff = new TH1F(Form("h_event_time_diff_%d",RunNo),Form("Event time difference (Run %d);T (#mus);Counts",RunNo),100,0,12000);
    TH1F* h_event_time_diff_ns = new TH1F(Form("h_event_time_diff_ns_%d",RunNo),Form("Event time difference (Run %d);T (ns);Counts",RunNo),100,0,1000);

    TMultiGraph *mg = new TMultiGraph();
    TGraph *graph1 = new TGraph();
    graph1->SetLineColor(1);
    graph1->SetLineWidth(2);
    TGraph *graph2 = new TGraph();
    graph2->SetLineColor(2);
    graph2->SetLineWidth(2);
    TGraph *graph3 = new TGraph();
    graph3->SetLineColor(4);
    graph3->SetLineWidth(2);

    TMultiGraph *mg_rate = new TMultiGraph();
    TGraph *graph1_rate = new TGraph();
    graph1_rate->SetLineColor(1);
    graph1_rate->SetLineWidth(2);
    TGraph *graph2_rate = new TGraph();
    graph2_rate->SetLineColor(2);
    graph2_rate->SetLineWidth(2);
    TGraph *graph3_rate = new TGraph();
    graph3_rate->SetLineColor(4);
    graph3_rate->SetLineWidth(2);

    TMultiGraph *mg_trig_rate = new TMultiGraph();
    TGraph *graph_trig_rate = new TGraph();
    graph_trig_rate->SetLineColor(4);
    graph_trig_rate->SetLineWidth(2);

    Event_T_last = 0;
    Event_T_last_event = 0;
    for(Long64_t evt=0; evt<nentries; evt++){
        ch1->GetEntry(evt);
        if(evt==0) Event_T_0 = fEvtHdr_fEvtTime;
        if(evt>TSH_EventNo[TSH_Length-1]) break;

        Total_Events++;

        if(T_hms_hEDTM_tdcTimeRaw<10.&&evt>0){
            h_event_time_diff->Fill((fEvtHdr_fEvtTime-Event_T_last_event)*4/1000.); //6257202 ~ 250 MHz clock EDTM_Rate=40Hz
            h_event_time_diff_ns->Fill((fEvtHdr_fEvtTime-Event_T_last_event)*4);
        }
        Event_T_last_event = fEvtHdr_fEvtTime;

        if(T_hms_hEDTM_tdcTimeRaw>=EDTM_TDC_Low&&T_hms_hEDTM_tdcTimeRaw<=EDTM_TDC_High) EDTM_Event_Nocut++;

        //  TSH_Length, TSH_EventNo, TSH_EventTime, TSH_EDTM_Scaler
        if(evt==TSH_EventNo[N_TSH_Event]){
            Double_t DeltaT = (fEvtHdr_fEvtTime-Event_T_last)*1./250e6;
            Double_t EDTM_Scaler_Rate = (TSH_EDTM_Scaler[N_TSH_Event]-TSH_EDTM_Scaler[N_TSH_Event-1])/DeltaT;
            Double_t EDTM_Event_Rate = (EDTM_Event_Nocut-EDTM_Event_Nocut_last)/DeltaT;
            if(EDTM_Scaler_Rate<EDTM_RATE-2){
                if(EDTM_Event_Rate>0.5*EDTM_RATE+0.5*EDTM_Scaler_Rate){
                    Double_t Diff_diff = EDTM_Event_Nocut-(TSH_EDTM_Scaler[N_TSH_Event]-TSH_EDTM_Scaler[0]) - (EDTM_Event_Nocut_last-(TSH_EDTM_Scaler[N_TSH_Event-1]-TSH_EDTM_Scaler[0]));
                    if(Diff_diff>EDTM_RATE*0.5) EDTM_Scaler_correction += Diff_diff;
                }
            }

            Total_Time = (fEvtHdr_fEvtTime-Event_T_0)*1./250e6;
            graph1->SetPoint(N_TSH_Event,(fEvtHdr_fEvtTime-Event_T_0)*1./250e6,0);
            graph2->SetPoint(N_TSH_Event,(fEvtHdr_fEvtTime-Event_T_0)*1./250e6,(fEvtHdr_fEvtTime-Event_T_0)*1./250e6*EDTM_RATE-(TSH_EDTM_Scaler[N_TSH_Event]-TSH_EDTM_Scaler[0]+EDTM_Scaler_correction));
            graph3->SetPoint(N_TSH_Event,(fEvtHdr_fEvtTime-Event_T_0)*1./250e6,EDTM_Event_Nocut-(TSH_EDTM_Scaler[N_TSH_Event]-TSH_EDTM_Scaler[0]+EDTM_Scaler_correction));

            if(N_TSH_Event>0.5){
                graph1_rate->SetPoint(N_TSH_Event-1,(fEvtHdr_fEvtTime-Event_T_0)*1./250e6,EDTM_Scaler_Rate);
                graph2_rate->SetPoint(N_TSH_Event-1,(fEvtHdr_fEvtTime-Event_T_0)*1./250e6,EDTM_RATE);
                graph3_rate->SetPoint(N_TSH_Event-1,(fEvtHdr_fEvtTime-Event_T_0)*1./250e6,EDTM_Event_Rate);
                graph_trig_rate->SetPoint(N_TSH_Event-1,(fEvtHdr_fEvtTime-Event_T_0)*1./250e6,TSH_hL1ACCP_ScalerRate[N_TSH_Event]);
            }

            Event_T_last = fEvtHdr_fEvtTime;
            EDTM_Event_Nocut_last = EDTM_Event_Nocut;
            N_TSH_Event++;
        }

        if(NofPeriod_LowCurr==0){
            if(T_hms_hEDTM_tdcTimeRaw>=EDTM_TDC_Low&&T_hms_hEDTM_tdcTimeRaw<=EDTM_TDC_High) EDTM_Event += 1.;
            h_EDTM_tdcTimeRaw->Fill(T_hms_hEDTM_tdcTimeRaw);
        }
        else{
            LowCurr_Flag = 0;
            for(int nn=N_Period_LowCurr;nn<NofPeriod_LowCurr;nn++){
                if(evt>start_event[nn] && evt<end_event[nn]) LowCurr_Flag = 1;
                if(evt<end_event[nn]){
                    N_Period_LowCurr = nn;
                    break;
                }
            }
            if(LowCurr_Flag==0){
                if(T_hms_hEDTM_tdcTimeRaw>=EDTM_TDC_Low&&T_hms_hEDTM_tdcTimeRaw<=EDTM_TDC_High){
                    EDTM_Event += 1.;
                    h_EDTM_tdcTimeRaw->Fill(T_hms_hEDTM_tdcTimeRaw);
                }
            }
            else if(T_hms_hEDTM_tdcTimeRaw>=EDTM_TDC_Low&&T_hms_hEDTM_tdcTimeRaw<=EDTM_TDC_High) EDTM_Event_LowCurr++;
        }
        
    }

    Double_t EDTM_Scaler_LowCurr = 0.;
    Double_t Time_LowCurr = 0.;
    std::vector<double> start_event_time, end_event_time;
    for(int i=0;i<NofPeriod_LowCurr;i++){
        if(end_event[i]>TSH_EventNo[TSH_Length-1]){
            NofPeriod_LowCurr = i;
            break;
        }
        ch1->GetEntry(int(start_event[i]));
        start_event_time.push_back((fEvtHdr_fEvtTime-Event_T_0)*1./250e6);
        ch1->GetEntry(int(end_event[i]));
        end_event_time.push_back((fEvtHdr_fEvtTime-Event_T_0)*1./250e6);
        EDTM_Scaler_LowCurr += end_EDTM_scaler[i] - start_EDTM_scaler[i];
        Time_LowCurr += end_event_time[i] - start_event_time[i];
    }
    EDTM_Count = TSH_EDTM_Scaler[TSH_Length-1]-TSH_EDTM_Scaler[0]-EDTM_Scaler_LowCurr;

    // cout<<"EDTM scaler count after low beam current cut:"<<endl;
    // cout<<TSH_EDTM_Scaler[TSH_Length-1]-TSH_EDTM_Scaler[0]<<" - "<<EDTM_Scaler_LowCurr<<" = "<<EDTM_Count<<endl;
    // cout<<"EDTM event count after low beam current cut:"<<endl;
    // cout<<EDTM_Event_LowCurr+EDTM_Event<<" - "<<EDTM_Event_LowCurr<<" = "<<EDTM_Event<<endl;
    // cout<<"Total Diff: "<<(TSH_EDTM_Scaler[TSH_Length-1]-TSH_EDTM_Scaler[0])-EDTM_Event_Nocut<<endl;
    // cout<<"Difference after cut: "<<EDTM_Count - EDTM_Event<<endl;
    // cout<<"Time Low Curr * 40: "<<Time_LowCurr<<" * "<<40<<" = "<< Time_LowCurr*40<<endl;
    // cout<<"EDTM Scaler correction: "<<EDTM_Scaler_correction<<endl;

    EDTM_Count += EDTM_Scaler_correction;

    if(EDTM_Count>0.1) EDTM_LT = EDTM_Event/EDTM_Count;
    else EDTM_LT = -1.;
    result[0] = EDTM_LT;                        // 0. EDTM LT after 2uA cut
    result[1] = EDTM_Event;                     // 1. EDTM TDC events after 2uA cut
    result[2] = EDTM_Count;                     // 2. EDTM Scaler count after 2uA cut
    result[3] = EDTM_Event_LowCurr+EDTM_Event;  // 3. Total EDTM TDC events
    result[4] = TSH_EDTM_Scaler[TSH_Length-1]-TSH_EDTM_Scaler[0];  
                                                // 4. Total EDTM Scaler count
    result[5] = EDTM_Event_LowCurr;             // 5. EDTM TDC events   (I<2uA)
    result[6] = EDTM_Scaler_LowCurr;            // 6. EDTM Scaler count (I<2uA)
    result[7] = EDTM_Scaler_correction;         // 7. Missing EDTM Scaler count
    result[8] = EDTM_RATE;                      // 8. EDTM Rate
    result[9] = Time_LowCurr;                   // 9. Low beam current time(s)
    result[10] = Trigger_Rate;                  // 10. Trigger Rate (Hz)
    result[11] = PSFACTOR;                      // 11. Which trigger we're using

    //--------------------------------------Fit Method----------------------------------------
    Double_t EventTimeDiff_min = 0.;
    for(int i=1;i<=100;i++){
        double binContent = h_event_time_diff_ns->GetBinContent(i);
        double binLowEdge = h_event_time_diff_ns->GetBinLowEdge(i);
        if(binContent>1){
            EventTimeDiff_min = binLowEdge;
            break;
        }
    }
    cout<<EventTimeDiff_min<<endl;

    if(EventTimeDiff_min<1.){
        result[12] = -1;          // 12. Minimum time difference between two events (ns)
        result[13] = -1;          // 13. lambda[us^{-1}]
        result[14] = -1;          // 14. lambda_err[us^{-1}]
        result[15] = -1;          // 15. LT_Event_Tdiff_fit
        result[16] = -1;          // 16. LT_Event_Tdiff_fit_err
        result[17] = -1;          // 17. Total events
        result[18] = -1;          // 18. Total time in second
        result[19] = -1;          // 19. LT -- Good time percentage
    }

    double BinContent = h_event_time_diff->GetBinContent(1);
    TF1 *fitFunc_event = new TF1("fitFunc_event", "[0]*exp(-[1]*x) + [2]", 0, 12000);
    fitFunc_event->SetParameter(0, BinContent);
    fitFunc_event->SetParameter(1, 0.0006);
    fitFunc_event->SetParameter(2, 0);
    h_event_time_diff->Fit(fitFunc_event, "R", "", 0, 12000);

    Double_t fitParameters_event[3];
    fitFunc_event->GetParameters(fitParameters_event);
    Double_t lambda = fitParameters_event[1];
    Double_t lambda_err = fitFunc_event->GetParError(1);

    // 0--lambda[us^{-1}]; 7--lambda_err[us^{-1}]
    Double_t LT_fit = 1./(lambda/1000.*EventTimeDiff_min+1.);
    Double_t LT_fit_err = (lambda_err/1000.*EventTimeDiff_min)/(lambda/1000.*EventTimeDiff_min+1.); // relative err for denominator
    LT_fit_err = LT_fit_err*LT_fit; // absolute err for LT_fit


    Double_t LT_Good_Time_percent = 1-(Total_Events*EventTimeDiff_min)*1./1e9/Total_Time;

    result[12] = EventTimeDiff_min;             // 12. Minimum time difference between two events (ns)
    result[13] = lambda;                        // 13. lambda[us^{-1}]
    result[14] = lambda_err;                    // 14. lambda_err[us^{-1}]
    result[15] = LT_fit;                        // 15. LT_Event_Tdiff_fit
    result[16] = LT_fit_err;                    // 16. LT_Event_Tdiff_fit_err
    result[17] = Total_Events*1.;               // 17. Total events
    result[18] = Total_Time;                    // 18. Total time in second
    result[19] = LT_Good_Time_percent;          // 19. LT -- Good time percentage

    // cout<<"Number of Low Beam Current Period: "<<NofPeriod_LowCurr<<endl;

    // for(int i=0;i<NofPeriod_LowCurr;i++){
    //     cout<<start_event_time_TSH[i]<<"\t"<<end_event_time_TSH[i]<<endl;
    //     cout<<start_event_time[i]<<"\t"<<end_event_time[i]<<endl;
    //     cout<<endl;
    // }

    result[20] = PS;

    cout<<result[0]<<"\t"<<result[1]<<"\t"<<result[2]<<endl;
    // 0.997406        47685   47809

    //-----------------------------------Make Plots-------------------------------
    mg->Add(graph1);
    // mg->Add(graph2);
    mg->Add(graph3);
    TCanvas *canvas = new TCanvas("canvas", "Multiple Lines", 800, 1300);
    canvas->Divide(1,4);
    canvas->cd(1);
    canvas->SetLogy(0);
    mg->Draw("AL");
    mg->SetTitle(Form("EDTM Counts Difference (Run %d)",RunNo));
    mg->GetXaxis()->SetTitle("Time (250MHz Clock) in second");
    mg->GetYaxis()->SetTitle("Difference");
    mg->GetXaxis()->SetTitleSize(0.04);
    mg->GetYaxis()->SetTitleSize(0.04);

    TLegend *legend;
    legend = new TLegend(0.11, 0.11, 0.35, 0.3);
    legend->AddEntry(graph3, "Events - Scaler", "l");
    // legend->AddEntry(graph2, Form("%dHz counts - Scaler",int(EDTM_RATE)), "l");
    legend->AddEntry(graph1, "Scaler - Scaler", "l");
    legend->SetTextSize(0.04); 
    legend->Draw();

    TAxis *yAxis = mg->GetYaxis();
    Double_t ymin = yAxis->GetXmin();
    Double_t ymax = yAxis->GetXmax();

    for(int i=0;i<NofPeriod_LowCurr;i++){
        TBox *box = new TBox(start_event_time[i], ymin, end_event_time[i], ymax);
        box->SetFillColorAlpha(kGreen, 0.5);
        box->Draw("SAME");
    }

    canvas->cd(2);
    mg_rate->Add(graph1_rate);
    mg_rate->Add(graph2_rate);
    mg_rate->Add(graph3_rate);
    mg_rate->Draw("AL");
    mg_rate->SetTitle(Form("EDTM Counts Rate (Run %d)",RunNo));
    mg_rate->GetXaxis()->SetTitle("Time (250MHz Clock) in second");
    mg_rate->GetYaxis()->SetTitle("Rate (Hz)");
    mg_rate->GetXaxis()->SetTitleSize(0.04);
    mg_rate->GetYaxis()->SetTitleSize(0.04);

    TAxis *yAxis_rate = mg_rate->GetYaxis();
    Double_t ymin_rate = yAxis_rate->GetXmin();
    Double_t ymax_rate = yAxis_rate->GetXmax();

    TLegend *legend_rate;
    if(ymax_rate<50) legend_rate = new TLegend(0.7, 0.11, 0.89, 0.35);
    else legend_rate = new TLegend(0.11, 0.7, 0.35, 0.89);
    legend_rate->AddEntry(graph3_rate, "Event Rate", "l");
    legend_rate->AddEntry(graph2_rate, Form("%dHz",int(EDTM_RATE)), "l");
    legend_rate->AddEntry(graph1_rate, "Scaler Rate", "l");
    legend_rate->SetTextSize(0.04); 
    legend_rate->Draw();

    for(int i=0;i<NofPeriod_LowCurr;i++){
        TBox *box = new TBox(start_event_time[i], ymin_rate, end_event_time[i], ymax_rate);
        box->SetFillColorAlpha(kGreen, 0.5);
        box->Draw("SAME");
    }

    canvas->cd(3);
    gPad->SetLogy(0);
    mg_trig_rate->Add(graph_trig_rate);
    mg_trig_rate->Draw("AL");
    mg_trig_rate->SetTitle(Form("hL1ACCP Rate (Run %d)",RunNo));
    mg_trig_rate->GetXaxis()->SetTitle("Time (250MHz Clock) in second");
    mg_trig_rate->GetYaxis()->SetTitle("Rate (Hz)");
    mg_trig_rate->GetXaxis()->SetTitleSize(0.04);
    mg_trig_rate->GetYaxis()->SetTitleSize(0.04);

    canvas->cd(4);
    gPad->SetLogy();
    h_EDTM_tdcTimeRaw->Draw();

    canvas->SaveAs(Form("pdf/EDTM_Check_%d.pdf",RunNo));

    if(Output_Flag==1||Output_Flag==2){
        f1->cd();
        h_EDTM_tdcTimeRaw->Write();
        h_event_time_diff->Write();
        h_event_time_diff_ns->Write();
        f1->Close();
    }

    delete ch1;
    delete ch2;
    delete file;
    // delete histogram;
    delete inputFile_root;
    delete h_EDTM_tdcTimeRaw;
    // delete legend;
    // delete legend_rate;
    delete h_event_time_diff;
    delete h_event_time_diff_ns;
    // delete graph1;
    // delete graph2;
    // delete graph3;
    // delete mg;
    // delete mg_rate;
    // delete mg_trig_rate;
    // delete graph1_rate, graph2_rate, graph3_rate;
    // delete graph_trig_rate;
    // delete fitFunc_event;
    delete canvas;
}

Double_t ReadTILT(Int_t RunNo, Int_t PS_Flag){
    // TString FilePath_PS = "/volatile/hallc/nps/nps-ana/REPORT_OUTPUT/REPORT_OUTPUT2/COIN/SKIM/";
    TString FilePath_PS = "./ReportFiles1/";
    TString Input_File_PS = FilePath_PS + Form("skim_NPS_HMS_report_%d_-1.report",RunNo);
    std::ifstream inputFile(Input_File_PS);
    if(!inputFile.is_open()){
        // FilePath_PS = "/volatile/hallc/nps/nps-ana/REPORT_OUTPUT/COIN/SKIM/";
        TString FilePath_PS = "./ReportFiles2/";
        Input_File_PS = FilePath_PS + Form("skim_NPS_HMS_report_%d_-1.report",RunNo);
        inputFile.open(Input_File_PS);
    }
    if(!inputFile.is_open()){
        return -10.;
    } 

    Double_t TILT[6];
    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.find("HMS TRIG1 Computer Live Time :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "HMS TRIG1 Computer Live Time : %lf ", &TILT[0]);
        }
        if (line.find("HMS TRIG2 Computer Live Time :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "HMS TRIG2 Computer Live Time : %lf ", &TILT[1]);
        }
        if (line.find("HMS TRIG3 Computer Live Time :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "HMS TRIG3 Computer Live Time : %lf ", &TILT[2]);
        }
        if (line.find("HMS TRIG4 Computer Live Time :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "HMS TRIG4 Computer Live Time : %lf ", &TILT[3]);
        }
        if (line.find("HMS TRIG5 Computer Live Time :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "HMS TRIG5 Computer Live Time : %lf ", &TILT[4]);
        }
        if (line.find("HMS TRIG6 Computer Live Time :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "HMS TRIG6 Computer Live Time : %lf ", &TILT[5]);
            break;
        }
    }
    return TILT[PS_Flag-1]/100.;
}

Double_t ReadTriggerRate(Int_t RunNo, Int_t PS_Flag){
    // TString FilePath_PS = "/volatile/hallc/nps/nps-ana/REPORT_OUTPUT/REPORT_OUTPUT2/COIN/SKIM/";
    TString FilePath_PS = "./ReportFiles1/";
    TString Input_File_PS = FilePath_PS + Form("skim_NPS_HMS_report_%d_-1.report",RunNo);
    std::ifstream inputFile(Input_File_PS);
    if(!inputFile.is_open()){
        // FilePath_PS = "/volatile/hallc/nps/nps-ana/REPORT_OUTPUT/COIN/SKIM/";
        TString FilePath_PS = "./ReportFiles2/";
        Input_File_PS = FilePath_PS + Form("skim_NPS_HMS_report_%d_-1.report",RunNo);
        inputFile.open(Input_File_PS);
    }
    if(!inputFile.is_open()){
        return -10.;
    } 

    Double_t TRIG1, TRIG2, TRIG3, TRIG4, TRIG5, TRIG6;
    Double_t Rate1, Rate2, Rate3, Rate4, Rate5, Rate6;
    std::string line;
    while (std::getline(inputFile, line)) {

        if (line.find("hTRIG1 :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "hTRIG1 :        %lf         [ %lf kHz ]", &TRIG1, &Rate1);
        }

        if (line.find("hTRIG2 :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "hTRIG2 :        %lf         [ %lf kHz ]", &TRIG2, &Rate2);
        }

        if (line.find("hTRIG3 :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "hTRIG3 :        %lf         [ %lf kHz ]", &TRIG3, &Rate3);
        }

        if (line.find("hTRIG4 :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "hTRIG4 :        %lf         [ %lf kHz ]", &TRIG4, &Rate4);
        }

        if (line.find("hTRIG5 :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "hTRIG5 :        %lf         [ %lf kHz ]", &TRIG5, &Rate5);
        }

        if (line.find("hTRIG6 :") != std::string::npos)
        {
            std::sscanf(line.c_str(), "hTRIG6 :        %lf         [ %lf kHz ]", &TRIG6, &Rate6);
            break;
        }
 
    }

    Double_t Rate;
    if(PS_Flag==1) Rate = Rate1;
    else if(PS_Flag==2) Rate = Rate2;
    else if(PS_Flag==3) Rate = Rate3;
    else if(PS_Flag==4) Rate = Rate4;
    else if(PS_Flag==5) Rate = Rate5;
    else Rate = Rate6;

    return Rate;
}