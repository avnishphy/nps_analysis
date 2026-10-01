#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

// Small valid, incompatible-type, or incomplete-branch input. No production data.
int main(int argc, char** argv) {
    if (argc != 3) return 64;
    const std::string dir = argv[1], mode = argv[2];
    if (mode == "check-output") {
        TFile f(dir.c_str(), "READ");
        auto* t = dynamic_cast<TTree*>(f.Get("simulation"));
        if (!t || t->GetEntries() != 120) return 1;
        Float_t weight = 0;
        if (t->SetBranchAddress("full_weight", &weight) < 0) return 1;
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            const Float_t input = 2.0f + .1f*((i%40)%7);
            if (t->GetEntry(i) <= 0 || std::abs(weight - input/40.0) > 1e-8) return 1;
        }
        return 0;
    }
    if (mode != "valid" && mode != "wrong-type" && mode != "short-branch") return 64;
    std::filesystem::create_directories(dir);
    for (const std::string channel : {"exclusive", "sidis", "delta"}) {
        TFile f((dir + "/" + channel + ".root").c_str(), "RECREATE");
        TTree t("nerd", "input validation fixture");
        std::vector<double> clust_E, clust_X, clust_Y;
        Float_t weight = 1, hsdelta = 0, hsyptar = 0, hsxptar = 0, hsytar = 0;
        Double_t wide_weight = 1, phot1_vx = 0;
        Float_t sc_e_E = 3200, sc_e_Px = 0, sc_e_Py = 1000, sc_e_Pz = 3000;
        Float_t vtx_z = 0, siglab = 1, sigcm = 1;
        t.Branch("clust_E", &clust_E); t.Branch("clust_X", &clust_X); t.Branch("clust_Y", &clust_Y);
        if (mode == "valid") t.Branch("Weight", &weight);
        if (mode == "wrong-type") t.Branch("Weight", &wide_weight);
        t.Branch("hsdelta", &hsdelta); t.Branch("hsyptar", &hsyptar);
        t.Branch("hsxptar", &hsxptar); t.Branch("hsytar", &hsytar);
        t.Branch("sc_e_E", &sc_e_E); t.Branch("sc_e_Px", &sc_e_Px);
        t.Branch("sc_e_Py", &sc_e_Py); t.Branch("sc_e_Pz", &sc_e_Pz);
        t.Branch("phot1_vx", &phot1_vx); t.Branch("vtx_z", &vtx_z);
        t.Branch("siglab", &siglab); t.Branch("sigcm", &sigcm);
        for (int i = 0; i < 40; ++i) {
            clust_E = {1.75 + .002*i, 1.55 + .001*i, .3};
            clust_X = {-15.0 + .02*i, 15.0 - .01*i, 0.};
            clust_Y = {-2. + .05*(i%20), 2. - .04*(i%20), 0.};
            weight = 2.0f + .1f*(i%7); wide_weight = weight;
            t.Fill();
        }
        if (mode == "short-branch") {
            auto* branch = t.Branch("Weight", &weight);
            for (int i = 0; i < 39; ++i) branch->Fill();
        }
        if (t.Write() <= 0) return 1;
    }
    return 0;
}
