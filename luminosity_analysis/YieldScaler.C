void extract_yield_scaler(Int_t run_start = 1524, Int_t run_end = 1528) {
  std::string base_dir = "/w/hallc-scshelf2102/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/SKIM/";
  std::vector<Double_t> yield_scalers, beam_currents;
  std::vector<Int_t> runs;
  std::vector<double> relative_yield_scalers;

  std::ofstream debug("yield_scaler_debug.txt");

  for (int run = run_start; run <= run_end; ++run) {
      std::string filename = base_dir + "skim_NPS_HMS_report_" + std::to_string(run) + "_-1.report";

      if (gSystem->AccessPathName(filename.c_str())) {
          std::cerr << "WARNING: File not found: " << filename << std::endl;
          debug << "WARNING: File not found: " << filename << std::endl;
          continue;
      }

      std::ifstream infile(filename);
      std::string line;
      int line_number = 0;

      double bcm_charge = -1;
      double bcm_current = -1;
      double edtm_trig = -1;
      double htrig_val = -1;
      int ps_index = -1;
      std::map<int, int> ps_factors;
      std::map<int, double> htrigs;

      while (std::getline(infile, line)) {
          ++line_number;

          if (line_number == 40 && line.find("BCM4A Beam Cut Current") != std::string::npos)
              sscanf(line.c_str(), "BCM4A Beam Cut Current: %lf uA", &bcm_current);

          if (line_number == 47 && line.find("BCM4A Beam Cut Charge") != std::string::npos)
              sscanf(line.c_str(), "BCM4A Beam Cut Charge: %lf uC", &bcm_charge);

          if (line_number >= 52 && line_number <= 61) {
              for (int i = 1; i <= 6; ++i) {
                  std::string key = "Ps" + std::to_string(i) + "_factor";
                  if (line.find(key) != std::string::npos) {
                      int val;
                      sscanf(line.c_str(), (key + " = %d").c_str(), &val);
                      ps_factors[i] = val;
                  }
              }
          }

          if (line_number >= 109 && line_number <= 114) {
              for (int i = 1; i <= 6; ++i) {
                  std::string key = "hTRIG" + std::to_string(i);
                  if (line.find(key) != std::string::npos) {
                      double val;
                      sscanf(line.c_str(), (key + " : %lf").c_str(), &val);
                      htrigs[i] = val;
                  }
              }
          }

          if (line_number == 120 && line.find("EDTM Triggers") != std::string::npos)
              sscanf(line.c_str(), "EDTM Triggers           : %lf", &edtm_trig);
      }

      infile.close();

      for (const auto& kv : ps_factors) {
          if (kv.second > 0) {
              ps_index = kv.first;
              break;
          }
      }

      if (ps_index == -1 || htrigs.find(ps_index) == htrigs.end()) {
          debug << "Run " << run << ": Missing PS or hTRIG info.\n";
          continue;
      }

      htrig_val = htrigs[ps_index];
      double number_scaler = htrig_val - edtm_trig;
      double yield = (bcm_charge > 0) ? number_scaler / bcm_charge : -1;

      debug << "Run: " << run << "\n";
      debug << "Beam Current: " << bcm_current << " uA\n";
      debug << "Charge: " << bcm_charge << " uC\n";
      debug << "hTRIG" << ps_index << ": " << htrig_val << "\n";
      debug << "EDTM Triggers: " << edtm_trig << "\n";
      debug << "Number Scaler: " << number_scaler << "\n";
      debug << "Yield Scaler: " << yield << "\n\n";

      runs.push_back(run);
      beam_currents.push_back(bcm_current);
      yield_scalers.push_back(yield);
  }

  double min_yield_scaler = *std::min_element(yield_scalers.begin(), yield_scalers.end());

  for (const auto& y : yield_scalers){
    relative_yield_scalers.push_back(y / min_yield_scaler);
  }

  // Plot
  TCanvas* c = new TCanvas("c", "Relative Yield Scaler vs Beam Current", 800, 600);
  TGraph* g = new TGraph(beam_currents.size(), &beam_currents[0], &relative_yield_scalers[0]);
  g->SetTitle("Relative Yield Scaler vs. Beam Current ;Beam Current (uA); Relative Yield Scaler (counts/uC)");
  g->SetMarkerStyle(20);
  g->Draw("AP");

  c->SaveAs("yield_scaler_vs_current.png");
}
