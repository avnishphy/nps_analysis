/// @brief this is a script to analyze the reference times.
/// @param infile 
void reftime_quantify_1(TString infile) {
    TFile *f = new TFile(infile);

    if (f->IsZombie() || !f->IsOpen()) {
        cout << "Error opening file: " << infile << endl;
        delete f;
        return;
    }


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

        double back_sig_ratio_tot = static_cast<double>(background_tot)/signal_tot;

        cout << "Total background/signal ratio (in percentage): " << back_sig_ratio_tot * 100 << endl;

        int sig_lost_cut = T->GetEntries(Form("T.hms.hDCREF%d_tdcTimeRaw<14400 && T.hms.hDCREF%d_tdcMultiplicity==1", ref_i, ref_i));
        int sig_remain_cut = T->GetEntries(Form("T.hms.hDCREF%d_tdcTimeRaw>=14400 && T.hms.hDCREF%d_tdcMultiplicity==1", ref_i, ref_i));
    
        //double signal_ratio_cut = static_cast<double>(sig_lost_cut) / sig_remain_cut; this definition is flawed
        double signal_ratio_cut = static_cast<double>(sig_lost_cut) / signal_tot;

        int background_remain_cut = (T->GetEntries(Form("T.hms.hDCREF%d_tdcTimeRaw>=14400", ref_i))) - sig_remain_cut;

        double back_sig_ratio_cut = static_cast<double>(background_remain_cut) / sig_remain_cut;
    
        cout << "Signal entries lost after the cut: " << sig_lost_cut << endl;
        cout << "Signal entries remaining after the cut: " << sig_remain_cut << endl;
        cout << "Signal loss ratio (in percentage): " << signal_ratio_cut*100 << endl;
        cout << "Background/Signal ratio after cut applied (in percentage): " << back_sig_ratio_cut * 100 << endl;
        cout << "-------------------------------" << endl;
    }

    for (int ref_i = 1; ref_i <= 2; ++ref_i) {
        cout << "hT" << ref_i << endl;

        const char *branch_name = Form("T.hms.hT%d_tdcMultiplicity", ref_i);  // Replace with the name of your branch
        // Create a histogram from the branch data
        TH1F *histogram = new TH1F("histogram", "Title;X-axis;Y-axis", 10, 0, 10);  // Adjust bins, xmin, xmax as needed
        TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
        T->Draw(Form("%s>>histogram", branch_name));
        delete canvas;
        // Find the bin with the maximum content
        int max_bin = histogram->GetMaximumBin();

        int x_max_mult = max_bin - 1;

        int signal_tot = T->GetEntries(Form("T.hms.hT%d_tdcTimeRaw && T.hms.hT%d_tdcMultiplicity==%d", ref_i, ref_i, x_max_mult));

        cout << "Total signal (highest multiplicity) entries: " << signal_tot << endl;

        int background_tot = (T->GetEntries(Form("T.hms.hT%d_tdcTimeRaw", ref_i))) - signal_tot;
	
	cout << "Total background entries: " << background_tot << endl;

        double back_sig_ratio_tot = static_cast<double>(background_tot)/signal_tot;

        cout << "Total background/signal ratio (in percentage): " << back_sig_ratio_tot * 100<< endl;

        int sig_lost_cut = T->GetEntries(Form("T.hms.hT%d_tdcTimeRaw<1700 && T.hms.hT%d_tdcMultiplicity==%d", ref_i, ref_i, x_max_mult));
        int sig_remain_cut = T->GetEntries(Form("T.hms.hT%d_tdcTimeRaw>=1700 && T.hms.hT%d_tdcMultiplicity==%d", ref_i, ref_i, x_max_mult));
    
        //double signal_ratio_cut = static_cast<double>(sig_lost_cut) / sig_remain_cut; this definition is flawed
        double signal_ratio_cut = static_cast<double>(sig_lost_cut) / signal_tot;

        int background_remain_cut = (T->GetEntries(Form("T.hms.hT%d_tdcTimeRaw>=1700", ref_i))) - sig_remain_cut;

        double back_sig_ratio_cut = static_cast<double>(background_remain_cut) / sig_remain_cut;
    
        cout << "Signal entries lost after the cut: " << sig_lost_cut << endl;
        cout << "Signal entries remaining after the cut: " << sig_remain_cut << endl;
        cout << "Signal loss ratio (in percentage): " << signal_ratio_cut*100 << endl;
        cout << "Background/Signal ratio after cut applied (in percentage): " << back_sig_ratio_cut * 100 << endl;
        cout << "-------------------------------" << endl;

        delete histogram;
    }


    f->Close();
    delete f;
}




// /lustre19/expphy/volatile/hallc/nps/singhav/ROOTfiles/NPS/TIMING/nps_hms_noReferenceTime_3013_1000000.root
