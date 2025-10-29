// ======================================================================
// pid_analysis1.C
// ======================================================================
// Purpose  : PID diagnostics on skimmed HMS+NPS data using available branches
//            Produce 1D and 2D histograms for theta, phi, delta, y.
//            Supports per-run processing (runlist.txt) and optional
//            graphical 2D cuts (TCutG) overlay & saving, plus R-function.
// Usage    : root -l -b -q "src/pid_analysis1.C+()"
//           or: root -l -b -q "src/pid_analysis1.C+(\"/path/skim/\",\"/path/out/\",\"config/runlist.txt\",\"config/pid_cuts.root\")"
// ======================================================================

#include "utils.C"

#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

// ----------------------
// Helper structs & funcs
// ----------------------
struct Point2D { double x, y; };

// Minimum distance from (x,y) to segment p1-p2
double distanceToSegment(double x, double y, const Point2D &p1, const Point2D &p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double L2 = dx*dx + dy*dy;
    if (L2 == 0) return hypot(x - p1.x, y - p1.y);
    double t = ((x - p1.x) * dx + (y - p1.y) * dy) / L2;
    t = max(0.0, min(1.0, t));
    double projx = p1.x + t * dx;
    double projy = p1.y + t * dy;
    return hypot(x - projx, y - projy);
}

// Point-in-polygon test (ray-casting)
bool pointInsidePolygon(double x, double y, const vector<Point2D> &poly) {
    bool inside = false;
    int n = poly.size();
    for (int i = 0, j = n - 1; i < n; j = i++) {
        bool cond = ((poly[i].y > y) != (poly[j].y > y)) &&
                    (x < (poly[j].x - poly[i].x) * (y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x);
        if (cond) inside = !inside;
    }
    return inside;
}

// R-function for one polygon
double RFunction(double x, double y, const vector<Point2D> &polygon) {
    if (polygon.empty()) return -1e9;
    double min_dist = 1e9;
    int n = polygon.size();
    for (int i = 0; i < n; i++) {
        const Point2D &p1 = polygon[i];
        const Point2D &p2 = polygon[(i + 1) % n];
        double d = distanceToSegment(x, y, p1, p2);
        if (d < min_dist) min_dist = d;
    }
    return pointInsidePolygon(x, y, polygon) ? min_dist : -min_dist;
}

// Combine 6 R-functions
double R_total(double ytg, double dtg, double thtg, double phtg,
               const vector<Point2D> &poly_phi_dt,
               const vector<Point2D> &poly_phi_th,
               const vector<Point2D> &poly_phi_y,
               const vector<Point2D> &poly_th_dt,
               const vector<Point2D> &poly_th_y,
               const vector<Point2D> &poly_y_dt) {

    vector<double> Rs = {
        RFunction(phtg, dtg, poly_phi_dt),
        RFunction(phtg, thtg, poly_phi_th),
        RFunction(phtg, ytg, poly_phi_y),
        RFunction(thtg, dtg, poly_th_dt),
        RFunction(thtg, ytg, poly_th_y),
        RFunction(ytg, dtg, poly_y_dt)
    };
    return *min_element(Rs.begin(), Rs.end());
}

// Convert TCutG to polygon points
vector<Point2D> cutToPolygon(TCutG* cut) {
    vector<Point2D> poly;
    if (!cut) return poly;
    int n = cut->GetN();
    for (int i = 0; i < n; i++) poly.push_back({cut->GetX()[i], cut->GetY()[i]});
    return poly;
}

// ----------------------
// Runlist reader
// ----------------------
vector<int> readRunList(const string &fname) {
    vector<int> runs;
    ifstream in(fname);
    if (!in.is_open()) { logmsg(WARN, "Cannot open runlist: " + fname); return runs; }
    string line;
    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        try { runs.push_back(stoi(line)); } catch (...) { logmsg(WARN, "Skipping invalid runlist line: " + line); }
    }
    return runs;
}

// ======================
// Main PID analysis
// ======================
void pid_analysis(const TString &skimDir_in="output/skimmed/",
                   const TString &outPlotDir_in="output/plots/pid/",
                   const TString &runlistFile="config/runlist_x60_4a.txt",
                   const TString &cutFile="config/pid_cuts_x60_4a.root") {

    TStopwatch sw_total; sw_total.Start();
    logmsg(INFO, "================ PID Analysis ================");

    TString skimDir = skimDir_in.EndsWith("/") ? skimDir_in : skimDir_in + "/";
    TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
    gSystem->mkdir(outPlotDir, true);

    vector<int> runs = readRunList(runlistFile.Data());
    if (runs.empty()) { logmsg(ERROR, "No runs found in runlist"); return; }

    // Load cuts
    TFile *fcut = nullptr;
    if (!gSystem->AccessPathName(cutFile.Data())) {
        fcut = TFile::Open(cutFile.Data(), "READ");
        if (fcut && !fcut->IsZombie()) logmsg(INFO, "Loaded cuts from " + string(cutFile.Data()));
        else { logmsg(WARN, "Could not open cuts file"); if (fcut) { fcut->Close(); delete fcut; fcut = nullptr; } }
    } else logmsg(INFO, "No cuts file found — proceeding without cuts");

    for (auto run : runs) {
        logmsg(INFO, Form("---- Processing run %d ----", run));
        TStopwatch sw; sw.Start();
        TString runOutDir = Form("%s/run%d", outPlotDir.Data(), run);
        gSystem->mkdir(runOutDir, true);

        // Build TChain
        TString pattern = Form("%s/skim_run%d.root", skimDir.Data(), run);
        TChain chain("T");
        int added = chain.Add(pattern);
        if (added == 0) { logmsg(WARN, Form("No files for run %d", run)); continue; }

        Long64_t nentries = chain.GetEntries();
        if (nentries == 0) continue;

        // Enable branches
        vector<string> pid_branches = {"H.gtr.th", "H.gtr.ph", "H.gtr.dp", "H.gtr.y"};
        enableBranches(&chain, pid_branches);

        Double_t H_th=0, H_ph=0, H_dp=0, H_y=0;
        chain.SetBranchAddress("H.gtr.th", &H_th);
        chain.SetBranchAddress("H.gtr.ph", &H_ph);
        chain.SetBranchAddress("H.gtr.dp", &H_dp);
        chain.SetBranchAddress("H.gtr.y",  &H_y);

        // Histograms
        TH1D *h_th = new TH1D("h_th", Form("theta (run %d);#theta [rad];counts",run),100,-0.1,0.1);
        TH1D *h_ph = new TH1D("h_ph", Form("phi (run %d);#phi [rad];counts",run),100,-0.1,0.1);
        TH1D *h_dp = new TH1D("h_dp", Form("delta (run %d);#delta [%%];counts",run),100,-16,16);
        TH1D *h_y  = new TH1D("h_y",  Form("y (run %d);y [cm];counts",run),100,-6,6);

        TH2D *h_th_vs_dp = new TH2D("h_th_vs_dp", Form("#theta vs #delta (run %d);#delta [%%];#theta [rad]",run),100,-16,16,100,-0.1,0.1);
        TH2D *h_phi_vs_dp = new TH2D("h_phi_vs_dp", Form("#phi vs #delta (run %d);#delta [%%];#phi [rad]",run),100,-16,16,100,-0.1,0.1);
        TH2D *h_phi_vs_y  = new TH2D("h_phi_vs_y",  Form("#phi vs y (run %d);y [cm];#phi [rad]",run),100,-6,6,100,-0.1,0.1);
        TH2D *h_th_vs_ph  = new TH2D("h_th_vs_ph",  Form("#theta vs #phi (run %d);#phi [rad];#theta [rad]",run),100,-0.1,0.1,100,-0.1,0.1);
        TH2D *h_y_vs_dp   = new TH2D("h_y_vs_dp",   Form("y vs #delta (run %d);#delta [%%];y [cm]",run),100,-16,16,100,-6,6);
        TH2D *h_th_vs_y   = new TH2D("h_th_vs_y",   Form("#theta vs y (run %d);y [cm];#theta [rad]",run),100,-6,6,100,-0.1,0.1);
        TH1D *h_R = new TH1D("h_R", Form("R-function (run %d);R-value;counts",run),100,-1,1);

        // Load cuts
        TCutG *cut_th_dp=nullptr, *cut_phi_dp=nullptr, *cut_y_dp=nullptr;
        TCutG *cut_th_y=nullptr, *cut_phi_y=nullptr, *cut_th_ph=nullptr;

        if (fcut && !fcut->IsZombie()) {
            const char* names[] = {"cut_th_dp","cut_phi_dp","cut_y_dp","cut_th_y","cut_phi_y","cut_th_ph"};
            TCutG** ptrs[] = {&cut_th_dp,&cut_phi_dp,&cut_y_dp,&cut_th_y,&cut_phi_y,&cut_th_ph};
            for (int i=0;i<6;i++){
                *ptrs[i] = (TCutG*)fcut->Get(names[i]);
                if (*ptrs[i]) logmsg(INFO,Form("Loaded %s",(*ptrs[i])->GetName()));
            }
        }

        // Convert to polygons
        auto poly_th_dp  = cutToPolygon(cut_th_dp);
        auto poly_phi_dp = cutToPolygon(cut_phi_dp);
        auto poly_y_dp   = cutToPolygon(cut_y_dp);
        auto poly_th_y   = cutToPolygon(cut_th_y);
        auto poly_phi_y  = cutToPolygon(cut_phi_y);
        auto poly_th_ph  = cutToPolygon(cut_th_ph);

        // Event loop
        Long64_t count_in=0, count_th_dp=0, count_phi_dp=0, count_phi_y=0, count_y_dp=0, count_th_y=0, count_all=0, count_Rpass=0;
        double Rcut = 0.001;

        for (Long64_t i=0;i<nentries;i++){
            chain.GetEntry(i);
            count_in++;

            h_th->Fill(H_th); h_ph->Fill(H_ph); h_dp->Fill(H_dp); h_y->Fill(H_y);
            h_th_vs_dp->Fill(H_dp,H_th);
            h_phi_vs_dp->Fill(H_dp,H_ph);
            h_phi_vs_y->Fill(H_y,H_ph);
            h_th_vs_ph->Fill(H_ph,H_th);
            h_y_vs_dp->Fill(H_dp,H_y);
            h_th_vs_y->Fill(H_y,H_th);

            bool pass_th_dp  = cut_th_dp  ? cut_th_dp->IsInside(H_dp,H_th)  : true;
            bool pass_phi_dp = cut_phi_dp ? cut_phi_dp->IsInside(H_dp,H_ph) : true;
            bool pass_y_dp   = cut_y_dp   ? cut_y_dp->IsInside(H_dp,H_y)    : true;
            bool pass_th_y   = cut_th_y   ? cut_th_y->IsInside(H_y,H_th)    : true;
            bool pass_phi_y  = cut_phi_y  ? cut_phi_y->IsInside(H_y,H_ph)   : true;
            bool pass_th_ph  = cut_th_ph  ? cut_th_ph->IsInside(H_ph,H_th)  : true;

            if (pass_th_dp) count_th_dp++;
            if (pass_phi_dp) count_phi_dp++;
            if (pass_phi_y) count_phi_y++;
            if (pass_y_dp) count_y_dp++;
            if (pass_th_y) count_th_y++;

            if (pass_th_dp && pass_phi_dp && pass_phi_y && pass_y_dp && pass_th_y && pass_th_ph)
                count_all++;

            double Rval = R_total(H_y,H_dp,H_th,H_ph,poly_phi_dp,poly_th_ph,poly_phi_y,poly_th_dp,poly_th_y,poly_y_dp);
            h_R->Fill(Rval);
            if (Rval >= Rcut) count_Rpass++;
        }

        // Summary
        logmsg(INFO, Form("Run %d: total entries = %lld",run,count_in));
        logmsg(INFO, Form("Cuts: th_dp=%.2f%%, phi_dp=%.2f%%, phi_y=%.2f%%, y_dp=%.2f%%, th_y=%.2f%%",
            100.*count_th_dp/count_in, 100.*count_phi_dp/count_in,
            100.*count_phi_y/count_in, 100.*count_y_dp/count_in, 100.*count_th_y/count_in));
        logmsg(INFO, Form("Events passing all cuts = %lld (%.2f%%)",count_all,100.*count_all/count_in));
        logmsg(INFO, Form("Events passing R-cut = %lld (%.2f%%)",count_Rpass,100.*count_Rpass/count_in));

        // Plot saving
        TCanvas *c = new TCanvas("c","",900,700);
        c->SetRightMargin(0.15);
        struct PlotDesc { TH2D* h; const char* tag; TCutG* cut; };
        vector<PlotDesc> plots = {
            {h_th_vs_dp,"th_vs_dp",cut_th_dp},
            {h_phi_vs_dp,"phi_vs_dp",cut_phi_dp},
            {h_phi_vs_y,"phi_vs_y",cut_phi_y},
            {h_th_vs_ph,"th_vs_ph",cut_th_ph},
            {h_y_vs_dp,"y_vs_dp",cut_y_dp},
            {h_th_vs_y,"th_vs_y",cut_th_y}
        };
        for (auto &p : plots){
            p.h->Draw("COLZ");
            if (p.cut) { p.cut->SetLineColor(kRed); p.cut->SetLineWidth(2); p.cut->Draw("same"); }
            c->SaveAs(Form("%s/pid_%s_run%d.png",runOutDir.Data(),p.tag,run));
            c->Clear();
        }

        // Save ROOT outputs
        TFile fout(Form("%s/pid_histos_run%d.root",runOutDir.Data(),run),"RECREATE");
        h_th->Write(); h_ph->Write(); h_dp->Write(); h_y->Write();
        h_th_vs_dp->Write(); h_phi_vs_dp->Write(); h_phi_vs_y->Write(); h_th_vs_ph->Write();
        h_y_vs_dp->Write(); h_th_vs_y->Write(); h_R->Write();
        if (cut_th_dp)  cut_th_dp->Write("cut_th_dp");
        if (cut_phi_dp) cut_phi_dp->Write("cut_phi_dp");
        if (cut_y_dp)   cut_y_dp->Write("cut_y_dp");
        if (cut_th_y)   cut_th_y->Write("cut_th_y");
        if (cut_phi_y)  cut_phi_y->Write("cut_phi_y");
        if (cut_th_ph)  cut_th_ph->Write("cut_th_ph");
        fout.Close();

        delete c;
        delete h_th; delete h_ph; delete h_dp; delete h_y;
        delete h_th_vs_dp; delete h_phi_vs_dp; delete h_phi_vs_y;
        delete h_th_vs_ph; delete h_y_vs_dp; delete h_th_vs_y; delete h_R;

        sw.Stop();
        logmsg(INFO, Form("[pid] Done run %d: %.2fs (CPU=%.2fs)",run,sw.RealTime(),sw.CpuTime()));
    }

    if (fcut){ fcut->Close(); delete fcut; }
    sw_total.Stop();
    logmsg(INFO, Form("All runs done in %.2fs total.", sw_total.RealTime()));
}
