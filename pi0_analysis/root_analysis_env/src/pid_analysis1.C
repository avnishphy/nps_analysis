// ======================================================================
// pid_analysis1.C
// ======================================================================
// PID diagnostics with graphical cuts + R-function with coordinate & area scaling
// Usage:
//   root -l -b -q "src/pid_analysis1.C+()"
//   or
//   root -l -b -q "src/pid_analysis1.C+(\"/path/skim/\",\"/path/out/\",\"config/runlist.txt\",\"config/pid_cuts.root\")"
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

struct PolyInfo {
    double xmin=0, xmax=0, ymin=0, ymax=0;
    double area=0;
    bool valid=false;
};

// Corrected distance from point (x,y) to segment p1-p2
double distanceToSegment_corrected(double x, double y, const Point2D &p1, const Point2D &p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double L2 = dx*dx + dy*dy;
    if (L2 == 0.0) return hypot(x - p1.x, y - p1.y); // degenerate segment

    double di = hypot(x - p1.x, y - p1.y); // distance to vertex i
    double dj = hypot(x - p2.x, y - p2.y); // distance to vertex j
    double L = sqrt(L2);

    // perpendicular distance from point to the line (extended infinitely)
    double d_perp = fabs(dy*(x - p1.x) - dx*(y - p1.y)) / L;

    // conditions for perpendicular to lie within segment
    if ((L2 + di*di - dj*dj > 0) && (L2 + dj*dj - di*di > 0)) {
        return d_perp;
    } else {
        return std::min(di, dj);
    }
}

// Point-in-polygon test (ray-casting)
bool pointInsidePolygon(double x, double y, const vector<Point2D> &poly) {
    bool inside = false;
    int n = (int)poly.size();
    if (n < 3) return false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        bool intersect = ((poly[i].y > y) != (poly[j].y > y)) &&
                         (x < (poly[j].x - poly[i].x) * (y - poly[i].y) / (poly[j].y - poly[i].y + 1e-300) + poly[i].x);
        if (intersect) inside = !inside;
    }
    return inside;
}

// Corrected R-function for a single polygon
double RFunction_raw(double x, double y, const vector<Point2D> &polygon) {
    if (polygon.empty()) return -1e9;
    double min_dist = 1e9;
    int n = (int)polygon.size();
    for (int i = 0; i < n; ++i) {
        const Point2D &p1 = polygon[i];
        const Point2D &p2 = polygon[(i+1)%n];
        double d = distanceToSegment_corrected(x, y, p1, p2); // use corrected distance
        if (d < min_dist) min_dist = d;
    }
    return pointInsidePolygon(x, y, polygon) ? min_dist : -min_dist;
}


// Convert TCutG to polygon points (original coordinates)
vector<Point2D> cutToPolygon(TCutG* cut) {
    vector<Point2D> poly;
    if (!cut) return poly;
    int n = cut->GetN();
    double *x = cut->GetX();
    double *y = cut->GetY();
    for (int i = 0; i < n; ++i) poly.push_back({x[i], y[i]});
    return poly;
}

// Compute polygon info (xmin/xmax/ymin/ymax/area) from polygon points
PolyInfo computePolyInfo(const vector<Point2D> &poly) {
    PolyInfo info;
    int n = (int)poly.size();
    if (n < 3) return info;
    info.valid = true;
    info.xmin = info.xmax = poly[0].x;
    info.ymin = info.ymax = poly[0].y;
    for (int i = 1; i < n; ++i) {
        info.xmin = min(info.xmin, poly[i].x);
        info.xmax = max(info.xmax, poly[i].x);
        info.ymin = min(info.ymin, poly[i].y);
        info.ymax = max(info.ymax, poly[i].y);
    }
    double area = 0.0;
    for (int i = 0; i < n; ++i) {
        int j = (i+1)%n;
        area += poly[i].x * poly[j].y - poly[j].x * poly[i].y;
    }
    info.area = fabs(area) * 0.5;
    return info;
}

// Normalize a point coordinate in plane using xmin,xmax,ymin,ymax
bool normalizePoint(double xin, double yin, double xmin, double xmax, double ymin, double ymax, double &xout, double &yout) {
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    if (dx == 0 || dy == 0) return false;
    xout = (xin - xmin) / dx;
    yout = (yin - ymin) / dy;
    return true;
}

// Normalize polygon coordinates (map polygon into normalized [0,1] box)
vector<Point2D> normalizePolygon(const vector<Point2D> &poly, const PolyInfo &info) {
    vector<Point2D> out;
    if (!info.valid) return out;
    double dx = info.xmax - info.xmin;
    double dy = info.ymax - info.ymin;
    if (dx == 0 || dy == 0) return out;
    out.reserve(poly.size());
    for (auto &p : poly) out.push_back({ (p.x - info.xmin)/dx, (p.y - info.ymin)/dy });
    return out;
}

// R for normalized polygon coordinates
double RFunction_norm(double xnorm, double ynorm, const vector<Point2D> &poly_norm) {
    return RFunction_raw(xnorm, ynorm, poly_norm);
}

// Combine 6 R-functions with coordinate normalization and area scaling
// double R_total_scaled(double ytg, double dtg, double thtg, double phtg,
//                       const vector<Point2D> &poly_phi_dt, const PolyInfo &info_phi_dt,
//                       const vector<Point2D> &poly_phi_th, const PolyInfo &info_phi_th,
//                       const vector<Point2D> &poly_phi_y,  const PolyInfo &info_phi_y,
//                       const vector<Point2D> &poly_th_dt,  const PolyInfo &info_th_dt,
//                       const vector<Point2D> &poly_th_y,   const PolyInfo &info_th_y,
//                       const vector<Point2D> &poly_y_dt,   const PolyInfo &info_y_dt) {

//     // Polygons normalized to unit box
//     vector<pair<vector<Point2D>, PolyInfo>> planes = {
//         { normalizePolygon(poly_phi_dt, info_phi_dt), info_phi_dt },
//         { normalizePolygon(poly_phi_th, info_phi_th), info_phi_th },
//         { normalizePolygon(poly_phi_y,  info_phi_y),  info_phi_y },
//         { normalizePolygon(poly_th_dt, info_th_dt),  info_th_dt },
//         { normalizePolygon(poly_th_y,  info_th_y),   info_th_y },
//         { normalizePolygon(poly_y_dt,  info_y_dt),   info_y_dt }
//     };

//     // Test points in original coordinates
//     vector<pair<double,double>> pts = {
//         { dtg, phtg },  // delta, phi  -> poly_phi_dt
//         { thtg, phtg }, // theta, phi  -> poly_phi_th
//         { ytg,  phtg }, // y, phi      -> poly_phi_y
//         { dtg, thtg },  // delta, theta -> poly_th_dt
//         { ytg,  thtg }, // y, theta    -> poly_th_y
//         { dtg,  ytg }   // delta, y    -> poly_y_dt
//     };

//     vector<double> scaledRs;
//     scaledRs.reserve(6);

//     for (int i=0; i<6; ++i) {
//         const auto &poly_norm = planes[i].first;
//         const PolyInfo &inf = planes[i].second;
//         double x = pts[i].first;
//         double y = pts[i].second;

//         if (!inf.valid || poly_norm.empty()) {
//             scaledRs.push_back(-1e9);
//             continue;
//         }

//         double xnorm, ynorm;
//         bool ok = normalizePoint(x, y, inf.xmin, inf.xmax, inf.ymin, inf.ymax, xnorm, ynorm);
//         if (!ok) { scaledRs.push_back(-1e9); continue; }

//         double rnorm = RFunction_norm(xnorm, ynorm, poly_norm); // uses corrected distance
//         // double scaled = (inf.area > 0) ? rnorm / sqrt(inf.area) : rnorm;
//         double scaled = rnorm;
//         scaledRs.push_back(scaled);
//     }

//     double Rmin = *min_element(scaledRs.begin(), scaledRs.end());
//     return Rmin;
// }

double R_total_scaled(double ytg, double dtg, double thtg, double phtg,
                      const vector<Point2D> &poly_phi_dt, const PolyInfo &info_phi_dt,
                      const vector<Point2D> &poly_phi_th, const PolyInfo &info_phi_th,
                      const vector<Point2D> &poly_phi_y,  const PolyInfo &info_phi_y,
                      const vector<Point2D> &poly_th_dt,  const PolyInfo &info_th_dt,
                      const vector<Point2D> &poly_th_y,   const PolyInfo &info_th_y,
                      const vector<Point2D> &poly_y_dt,   const PolyInfo &info_y_dt,
                      int* minPlane=nullptr) {

    vector<pair<vector<Point2D>, PolyInfo>> planes = {
        { normalizePolygon(poly_phi_dt, info_phi_dt), info_phi_dt },
        { normalizePolygon(poly_phi_th, info_phi_th), info_phi_th },
        { normalizePolygon(poly_phi_y,  info_phi_y),  info_phi_y },
        { normalizePolygon(poly_th_dt, info_th_dt),  info_th_dt },
        { normalizePolygon(poly_th_y,  info_th_y),   info_th_y },
        { normalizePolygon(poly_y_dt,  info_y_dt),   info_y_dt }
    };

    vector<pair<double,double>> pts = {
        { dtg, phtg },  { thtg, phtg }, { ytg,  phtg },
        { dtg, thtg },  { ytg,  thtg }, { dtg,  ytg }
    };

    vector<double> scaledRs(6, -1e9);
    for (int i=0; i<6; ++i) {
        const auto &poly_norm = planes[i].first;
        const PolyInfo &inf = planes[i].second;
        double x = pts[i].first;
        double y = pts[i].second;

        if (!inf.valid || poly_norm.empty()) continue;

        double xnorm, ynorm;
        if (!normalizePoint(x, y, inf.xmin, inf.xmax, inf.ymin, inf.ymax, xnorm, ynorm)) continue;

        double rnorm = RFunction_norm(xnorm, ynorm, poly_norm);
        double scaled = (inf.area > 0) ? rnorm / sqrt(inf.area) : rnorm;
        scaledRs[i] = scaled;
    }

    auto minIt = min_element(scaledRs.begin(), scaledRs.end());
    if (minPlane) *minPlane = distance(scaledRs.begin(), minIt);
    return *minIt;
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
void pid_analysis1(const TString &skimDir_in="output/skimmed/",
                  const TString &outPlotDir_in="output/plots/pid/",
                  const TString &runlistFile="config/runlist_x60_4b.txt",
                  const TString &cutFile="config/pid_cuts_x60_4b.root") {

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
        if (!fcut || fcut->IsZombie()) { logmsg(WARN,"Could not open cuts file"); if(fcut){fcut->Close(); delete fcut;} fcut=nullptr; }
        else logmsg(INFO, "Loaded cuts from " + string(cutFile.Data()));
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
        if (added == 0) {
            pattern = Form("%s/nps_hms_coin_%d_*.root", skimDir.Data(), run);
            added = chain.Add(pattern);
        }
        if (chain.GetNtrees() == 0) {
            logmsg(WARN, Form("No files for run %d in %s", run, skimDir.Data()));
            continue;
        }

        Long64_t nentries = chain.GetEntries();
        if (nentries == 0) { logmsg(WARN, Form("Run %d: zero entries", run)); continue; }
        logmsg(INFO, Form("Run %d: %d trees, %lld entries", run, chain.GetNtrees(), nentries));

        // Enable branches we need
        vector<string> pid_branches = {"H.gtr.th","H.gtr.ph","H.gtr.dp","H.gtr.y"};
        enableBranches(&chain, pid_branches);

        Double_t H_th=0, H_ph=0, H_dp=0, H_y=0;
        if (chain.GetBranch("H.gtr.th")) chain.SetBranchAddress("H.gtr.th", &H_th);
        if (chain.GetBranch("H.gtr.ph")) chain.SetBranchAddress("H.gtr.ph", &H_ph);
        if (chain.GetBranch("H.gtr.dp")) chain.SetBranchAddress("H.gtr.dp", &H_dp);
        if (chain.GetBranch("H.gtr.y"))  chain.SetBranchAddress("H.gtr.y",  &H_y);

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
        
        TH1D *h_R = new TH1D("h_R", Form("R-function (run %d);R-value;counts",run),200,0,0.5);
        TH1D* h_Rmin_plane = new TH1D("h_Rmin_plane", "Plane contributing Rmin;Plane index;counts", 6, 0, 6);


        // Load TCutG polygons
        TCutG *cut_th_dp=nullptr, *cut_phi_dp=nullptr, *cut_y_dp=nullptr;
        TCutG *cut_th_y=nullptr, *cut_phi_y=nullptr, *cut_th_ph=nullptr;

        if (fcut && !fcut->IsZombie()) {
            const char* names[] = {"cut_th_dp","cut_phi_dp","cut_y_dp","cut_th_y","cut_phi_y","cut_th_ph"};
            TCutG** ptrs[] = {&cut_th_dp,&cut_phi_dp,&cut_y_dp,&cut_th_y,&cut_phi_y,&cut_th_ph};
            for (int i=0;i<6;i++){
                TString nrun = Form("%s_run%d", names[i], run);
                TCutG* g = (TCutG*)fcut->Get(nrun);
                if (!g) g = (TCutG*)fcut->Get(names[i]);
                *ptrs[i] = g;
                if (g) logmsg(INFO, Form("Loaded %s", g->GetName()));
            }
        } else logmsg(WARN, "No valid cuts loaded (proceeding without graphical cuts).");

        // Convert TCutG -> polygon points
        auto poly_th_dp  = cutToPolygon(cut_th_dp);
        auto poly_phi_dp = cutToPolygon(cut_phi_dp);
        auto poly_y_dp   = cutToPolygon(cut_y_dp);
        auto poly_th_y   = cutToPolygon(cut_th_y);
        auto poly_phi_y  = cutToPolygon(cut_phi_y);
        auto poly_th_ph  = cutToPolygon(cut_th_ph);

        PolyInfo info_th_dp  = computePolyInfo(poly_th_dp);
        PolyInfo info_phi_dp = computePolyInfo(poly_phi_dp);
        PolyInfo info_y_dp   = computePolyInfo(poly_y_dp);
        PolyInfo info_th_y   = computePolyInfo(poly_th_y);
        PolyInfo info_phi_y  = computePolyInfo(poly_phi_y);
        PolyInfo info_th_ph  = computePolyInfo(poly_th_ph);

        // Counters
        Long64_t count_in=0, count_th_dp=0, count_phi_dp=0, count_phi_y=0, count_y_dp=0, count_th_y=0, count_all=0, count_Rpass=0;
        double Rcut = 0.0;

        // Event loop
        for (Long64_t i=0;i<nentries;i++){
            chain.GetEntry(i);
            ++count_in;

            // fill histos
            h_th->Fill(H_th); h_ph->Fill(H_ph); h_dp->Fill(H_dp); h_y->Fill(H_y);
            h_th_vs_dp->Fill(H_dp,H_th); h_phi_vs_dp->Fill(H_dp,H_ph); h_phi_vs_y->Fill(H_y,H_ph);
            h_th_vs_ph->Fill(H_ph,H_th); h_y_vs_dp->Fill(H_dp,H_y); h_th_vs_y->Fill(H_y,H_th);

            // per-cut counts
            bool pass_thdp  = cut_th_dp  ? cut_th_dp->IsInside(H_dp, H_th)  : true;
            bool pass_phdp  = cut_phi_dp ? cut_phi_dp->IsInside(H_dp, H_ph) : true;
            bool pass_ydp   = cut_y_dp   ? cut_y_dp->IsInside(H_dp, H_y)    : true;
            bool pass_thy   = cut_th_y   ? cut_th_y->IsInside(H_y, H_th)    : true;
            bool pass_phy   = cut_phi_y  ? cut_phi_y->IsInside(H_y, H_ph)   : true;
            bool pass_thph  = cut_th_ph  ? cut_th_ph->IsInside(H_ph, H_th)  : true;

            if (pass_thdp) count_th_dp++;
            if (pass_phdp) count_phi_dp++;
            if (pass_phy)  count_phi_y++;
            if (pass_ydp)  count_y_dp++;
            if (pass_thy)  count_th_y++;
            if (pass_thdp && pass_phdp && pass_phy && pass_ydp && pass_thy && pass_thph) ++count_all;

            // // Compute scaled R-value
            // double Rval = R_total_scaled(
            //     H_y, H_dp, H_th, H_ph,
            //     poly_phi_dp, info_phi_dp,
            //     poly_th_ph,  info_th_ph,
            //     poly_phi_y,  info_phi_y,
            //     poly_th_dp,  info_th_dp,
            //     poly_th_y,   info_th_y,
            //     poly_y_dp,   info_y_dp
            // );


            // h_R->Fill(Rval);
            // if (Rval >= Rcut) ++count_Rpass;
            int minPlane = -1;
            double Rval = R_total_scaled(
                H_y, H_dp, H_th, H_ph,
                poly_phi_dp, info_phi_dp,
                poly_th_ph,  info_th_ph,
                poly_phi_y,  info_phi_y,
                poly_th_dp,  info_th_dp,
                poly_th_y,   info_th_y,
                poly_y_dp,   info_y_dp,
                &minPlane // track which plane
            );

            h_R->Fill(Rval);
            if(minPlane>=0) h_Rmin_plane->Fill(minPlane);
            if (Rval >= Rcut) ++count_Rpass;

        }

        double denom = (count_in>0) ? double(count_in) : 1.0;
        logmsg(INFO, Form("Run %d: total entries = %lld", run, count_in));
        logmsg(INFO, Form("Cuts: th_dp=%.2f%%, phi_dp=%.2f%%, phi_y=%.2f%%, y_dp=%.2f%%, th_y=%.2f%%",
            100.0*count_th_dp/denom, 100.0*count_phi_dp/denom, 100.0*count_phi_y/denom,
            100.0*count_y_dp/denom, 100.0*count_th_y/denom));
        logmsg(INFO, Form("Events passing all cuts = %lld (%.2f%%)", count_all, 100.0*count_all/denom));
        logmsg(INFO, Form("Events passing R-cut (R >= %.6g) = %lld (%.2f%%)", Rcut, count_Rpass, 100.0*count_Rpass/denom));

        // Plot saving
        TCanvas *c = new TCanvas("c","",900,700);
        c->SetRightMargin(0.15);
        struct PlotDesc { TH2D* h; const char* tag; TCutG* cut; };
        vector<PlotDesc> plots = {
            {h_th_vs_dp, "th_vs_dp", cut_th_dp},
            {h_phi_vs_dp, "phi_vs_dp", cut_phi_dp},
            {h_phi_vs_y, "phi_vs_y", cut_phi_y},
            {h_th_vs_ph, "th_vs_ph", cut_th_ph},
            {h_y_vs_dp, "y_vs_dp", cut_y_dp},
            {h_th_vs_y, "th_vs_y", cut_th_y}
        };
        for (auto &p : plots) {
            p.h->Draw("COLZ");
            if (p.cut) { p.cut->SetLineColor(kRed); p.cut->SetLineWidth(2); p.cut->Draw("same"); }
            c->SaveAs(Form("%s/pid_%s_run%d.png", runOutDir.Data(), p.tag, run));
            c->Clear();
        }

        // Save histos
        TFile fout(Form("%s/pid_histos_run%d.root", runOutDir.Data(), run), "RECREATE");
        h_th->Write(); h_ph->Write(); h_dp->Write(); h_y->Write();
        h_th_vs_dp->Write(); h_phi_vs_dp->Write(); h_phi_vs_y->Write();
        h_th_vs_ph->Write(); h_y_vs_dp->Write(); h_th_vs_y->Write(); h_R->Write();
        if(cut_th_dp) cut_th_dp->Write("cut_th_dp");
        if(cut_phi_dp) cut_phi_dp->Write("cut_phi_dp");
        if(cut_y_dp) cut_y_dp->Write("cut_y_dp");
        if(cut_th_y) cut_th_y->Write("cut_th_y");
        if(cut_phi_y) cut_phi_y->Write("cut_phi_y");
        if(cut_th_ph) cut_th_ph->Write("cut_th_ph");
        h_Rmin_plane->Write();

        fout.Close();

        // cleanup
        delete c;
        delete h_th; delete h_ph; delete h_dp; delete h_y;
        delete h_th_vs_dp; delete h_phi_vs_dp; delete h_phi_vs_y;
        delete h_th_vs_ph; delete h_y_vs_dp; delete h_th_vs_y; delete h_R; delete h_Rmin_plane;

        sw.Stop();
        logmsg(INFO, Form("[pid] Done run %d: %.2fs (CPU=%.2fs)", run, sw.RealTime(), sw.CpuTime()));
    }

    if (fcut) { fcut->Close(); delete fcut; fcut=nullptr; }
    sw_total.Stop();
    logmsg(INFO, Form("All runs done in %.2fs total.", sw_total.RealTime()));
}
