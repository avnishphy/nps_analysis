#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void reftime_quantify_1(TString infile, int run_number) {
    TFile *f = new TFile(infile);

    if (f->IsZombie() || !f->IsOpen()) {
        cout << "Error opening file: " << infile << endl;
        delete f;
        return;
    }

    // Open CSV file for writing
    std::ofstream csv_file(Form("results_%d.csv", run_number));
    csv_file << "RunNumber,Reference,TotalSignal,TotalBackground,BackgroundSignalRatio,SignalLost,SignalRemaining,SignalLossRatio,BackgroundSignalRatioCut\n";

    // Open TXT file for writing
    std::ofstream txt_file(Form("results_%d.txt", run_number));

    std::vector<int> run_numbers;
    std::vector<double> total_background_signal_ratios;
    std::vector<double> signal_loss_ratios;
    std::vector<double> background_signal_ratios_cut;

    // Access the tree from the file
    TTree *T = (TTree*) f->Get("T");
    if (!T) {
        std::cout << "Error: Tree 'T' not found in file: " << infile << std::endl;
        f->Close();
        delete f;
        return;
    }

    for (int ref_i = 1; ref_i <= 4; ++ref_i) {
        cout << "hDCREF" << ref_i << endl;

        int signal_tot = T->GetEntries(Form("T.hms.hDCREF%d_tdcTimeRaw && T.hms.hDCREF%d_tdcMultiplicity==1", ref_i, ref_i));

        cout << "Total signal (highest multiplicity) entries: " << signal_tot << endl;

        int background_tot = (T->GetEntries(Form("T.hms.hDCREF%d_tdcTimeRaw", ref_i))) - signal_tot;

        cout << "Total background entries: " << background_tot << endl;

        double back_sig_ratio_tot = static_cast<double>(background_tot) / signal_tot;

        cout << "Total background/signal ratio (in percentage): " << back_sig_ratio_tot * 100 << endl;

        int sig_lost_cut = T->GetEntries(Form("T.hms.hDCREF%d_tdcTimeRaw<14800 && T.hms.hDCREF%d_tdcMultiplicity==1", ref_i, ref_i));
        int sig_remain_cut = T->GetEntries(Form("T.hms.hDCREF%d_tdcTimeRaw>=14800 && T.hms.hDCREF%d_tdcMultiplicity==1", ref_i, ref_i));
    
        double signal_ratio_cut = static_cast<double>(sig_lost_cut) / sig_remain_cut;

        int background_remain_cut = (T->GetEntries(Form("T.hms.hDCREF%d_tdcTimeRaw>=14800", ref_i))) - sig_remain_cut;

        double back_sig_ratio_cut = static_cast<double>(background_remain_cut) / sig_remain_cut;
    
        cout << "Signal entries lost after the cut: " << sig_lost_cut << endl;
        cout << "Signal entries remaining after the cut: " << sig_remain_cut << endl;
        cout << "Signal loss ratio (in percentage): " << signal_ratio_cut * 100 << endl;
        cout << "Background/Signal ratio after cut applied (in percentage): " << back_sig_ratio_cut * 100 << endl;
        cout << "-------------------------------" << endl;

        // Write data to CSV and TXT files
        csv_file << run_number << "," << ref_i << "," << signal_tot << "," << background_tot << "," << back_sig_ratio_tot * 100 << ","
                 << sig_lost_cut << "," << sig_remain_cut << "," << signal_ratio_cut * 100 << "," << back_sig_ratio_cut * 100 << "\n";
        txt_file << "RunNumber: " << run_number << "\nReference: " << ref_i << "\nTotal Signal: " << signal_tot << "\nTotal Background: " << background_tot
                 << "\nBackground/Signal Ratio: " << back_sig_ratio_tot * 100 << "\nSignal Lost: " << sig_lost_cut << "\nSignal Remaining: " << sig_remain_cut
                 << "\nSignal Loss Ratio: " << signal_ratio_cut * 100 << "\nBackground/Signal Ratio Cut: " << back_sig_ratio_cut * 100 << "\n-------------------------------\n";

        // Save data for plotting
        run_numbers.push_back(run_number);
        total_background_signal_ratios.push_back(back_sig_ratio_tot * 100);
        signal_loss_ratios.push_back(signal_ratio_cut * 100);
        background_signal_ratios_cut.push_back(back_sig_ratio_cut * 100);
    }

    // Plot the data
    TCanvas *c1 = new TCanvas("c1", "Total Background/Signal Ratio", 800, 600);
    TGraph *gr1 = new TGraph(run_numbers.size(), &run_numbers[0], &total_background_signal_ratios[0]);
    gr1->SetTitle("Total Background/Signal Ratio;Run Number;Total Background/Signal Ratio (%)");
    gr1->Draw("ALP");
    c1->SaveAs(Form("Total_Background_Signal_Ratio_%d.png", run_number));

    TCanvas *c2 = new TCanvas("c2", "Signal Loss Ratio", 800, 600);
    TGraph *gr2 = new TGraph(run_numbers.size(), &run_numbers[0], &signal_loss_ratios[0]);
    gr2->SetTitle("Signal Loss Ratio;Run Number;Signal Loss Ratio (%)");
    gr2->Draw("ALP");
    c2->SaveAs(Form("Signal_Loss_Ratio_%d.png", run_number));

    TCanvas *c3 = new TCanvas("c3", "Background/Signal Ratio After Cut", 800, 600);
    TGraph *gr3 = new TGraph(run_numbers.size(), &run_numbers[0], &background_signal_ratios_cut[0]);
    gr3->SetTitle("Background/Signal Ratio After Cut;Run Number;Background/Signal Ratio After Cut (%)");
    gr3->Draw("ALP");
    c3->SaveAs(Form("Background_Signal_Ratio_After_Cut_%d.png", run_number));

    // Close files
    csv_file.close();
    txt_file.close();

    f->Close();
    delete f;
}

void process_all_runs() {
    std::string directory = "/lustre19/expphy/volatile/hallc/nps/singhav/ROOTfiles/NPS/TIMING/";
    std::vector<int> run_numbers = {3013, 3032};

    for (const auto &run_number : run_numbers) {
        TString infile = Form("%snps_hms_noReferenceTime_%d_1000000.root", directory.c_str(), run_number);
        reftime_quantify_1(infile, run_number);
    }
}
