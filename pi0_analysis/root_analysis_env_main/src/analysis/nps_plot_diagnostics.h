#ifndef NPS_PLOT_DIAGNOSTICS_H
#define NPS_PLOT_DIAGNOSTICS_H

#include <fstream>
#include <map>
#include <memory>
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TNamed.h"
#include "nps_2d_mass_cut.h"

namespace npsplot {

// Diagnostic qualification only: never replace the stored production flags.
inline std::string ellipse_failure(const nps2d::Params& p, const nps2d::Config& cfg) {
    if (!p.valid || !p.ellipse_valid) return "ellipse fit unavailable";
    if (p.fit_subset_bins < cfg.auto_min_core_bins) return "ellipse fit has too few occupied bins";
    if (p.fit_subset_total_fraction < cfg.auto_min_core_total_fraction)
        return "ellipse fit subset below minimum fraction";
    if (!std::isfinite(p.cov_det) || p.cov_det <= 0.0) return "ellipse covariance invalid";
    return "";
}

class Diagnostics {
    std::vector<std::unique_ptr<TH1>> owned;
public:
    std::map<TH1D*, TH1D*> after;
    std::map<TH2D*, TH2D*> after2d;
    bool stage_pass = false;
    std::string range_source = "per-run pre-cut occupancy";
    std::string mass_status;

    TH1D* clone(TH1D* h, const std::string& suffix) {
        auto* result = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + suffix).c_str()));
        result->SetDirectory(nullptr); result->Reset(); result->Sumw2();
        result->SetFillStyle(0);
        owned.emplace_back(result);
        return result;
    }
    void register_after(TH1D* h) { after[h] = clone(h, "_after"); }
    void register_after(TH2D* h) {
        auto* result = static_cast<TH2D*>(h->Clone((std::string(h->GetName()) + "_after").c_str()));
        result->SetDirectory(nullptr); result->Reset(); result->Sumw2();
        owned.emplace_back(result); after2d[h] = result;
    }
    void fill(TH1D* h, double value, double weight) {
        if (stage_pass && after.count(h)) after.at(h)->Fill(value, weight);
    }
    void fill(TH2D* h, double x, double y, double weight) {
        if (stage_pass && after2d.count(h)) after2d.at(h)->Fill(x, y, weight);
    }

    // The producer chooses display bounds from its pre-cut histogram. Every
    // overlaid selection shares this view; no preparation command is needed.
    std::pair<int, int> display_bins(TH1* h, double cut_lo = NAN, double cut_hi = NAN) const {
        const std::string name = h->GetName();
        if (name.find("h_cut_cluster_x") == 0 || name.find("h_cut_cluster_y") == 0)
            return {1, h->GetNbinsX()}; // Fixed detector coverage, even for sparse runs.
        int first = 1, last = h->GetNbinsX();
        while (first <= last && h->GetBinContent(first) == 0.0) ++first;
        while (last >= first && h->GetBinContent(last) == 0.0) --last;
        if (first > last) { first = 1; last = h->GetNbinsX(); }
        const int pad = std::max(2, (last - first + 1) / 20);
        first -= pad; last += pad;
        for (double cut : {cut_lo, cut_hi}) if (std::isfinite(cut)) {
            first = std::min(first, h->GetXaxis()->FindFixBin(cut) - 2);
            last = std::max(last, h->GetXaxis()->FindFixBin(cut) + 2);
        }
        if (name.find("h_cut_cluster_e") == 0 || name.find("h_cut_mmiss_corr") == 0) {
            first = h->GetXaxis()->FindFixBin(0.4);
            last = std::max(last, first + 1);
        }
        return {std::max(1, std::min(h->GetNbinsX(), first)),
                std::max(1, std::min(h->GetNbinsX(), last))};
    }
    void draw_after(TH1D* before, const char* stage) const {
        const auto found = after.find(before);
        if (found == after.end()) return;
        TH1D* h = found->second;
        h->SetLineColor(kRed + 1); h->SetLineWidth(2); h->SetFillStyle(0);
        h->Draw("HIST SAME");
        auto* legend = new TLegend(0.50, 0.69, 0.94, 0.83);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.03);
        legend->AddEntry(before, Form("Before: %.0f", before->GetEntries()), "l");
        legend->AddEntry(h, Form("After %s: %.0f", stage, h->GetEntries()), "l");
        legend->Draw();
        TLatex text; text.SetNDC(); text.SetTextSize(0.03);
        if (before->GetEntries() > 0.0)
            text.DrawLatex(0.15, 0.79, Form("Retained %.2f%%", 100.0*h->GetEntries()/before->GetEntries()));
    }
    // Keep paired spatial axes identical. NPS cluster maps scale each population
    // independently; other maps retain the shared before color scale.
    // Display copies leave stored axes/counts intact.
    TH2D* draw_2d(TH2D* before, bool selected, bool available = true) const {
        TH2D* source = selected ? after2d.at(before) : before;
        auto* view = static_cast<TH2D*>(source->Clone((std::string(source->GetName()) + "_view").c_str()));
        view->SetDirectory(nullptr); view->SetBit(TObject::kCanDelete);
        const int nx = before->GetNbinsX(), ny = before->GetNbinsY();
        int x0 = nx, x1 = 1, y0 = ny, y1 = 1;
        bool occupied = false;
        double maximum = 0.0;
        const bool fixed = std::string(before->GetName()).find("h_debug_cluster_xy") == 0;
        for (int x = 1; x <= nx; ++x) for (int y = 1; y <= ny; ++y) {
            const double count = before->GetBinContent(x, y);
            maximum = std::max(maximum, fixed ? source->GetBinContent(x, y) : count);
            if (count == 0.0) continue;
            occupied = true;
            x0 = std::min(x0, x); x1 = std::max(x1, x);
            y0 = std::min(y0, y); y1 = std::max(y1, y);
        }
        if (occupied && !fixed) {
            const int px = std::max(2, (x1-x0+1)/20), py = std::max(2, (y1-y0+1)/20);
            view->GetXaxis()->SetRange(std::max(1, x0-px), std::min(nx, x1+px));
            view->GetYaxis()->SetRange(std::max(1, y0-py), std::min(ny, y1+py));
        } else {
            view->GetXaxis()->SetRange(1, nx); view->GetYaxis()->SetRange(1, ny);
        }
        view->SetStats(0); view->SetTitle("");
        view->SetMinimum(0.0); view->SetMaximum(std::max(1.0, maximum));
        view->Draw(available ? "COLZ" : "AXIS");
        return view;
    }
    void write() const {
        for (const auto& h : owned) h->Write();
        TNamed("diagnostic_range_source", range_source.c_str()).Write();
        TNamed("diagnostic_mass_cut_status", mass_status.c_str()).Write();
    }
};

// Reuse the producer's existing mass-cut canvas names, with actual event-filled
// histograms for every selector. Invalid fits still produce an inspectable PDF.
inline void draw_mass_comparison(const std::vector<nps2d::Point>& points,
        const std::vector<int>& decorrelation, const nps2d::Result& result,
        const nps2d::Config& cfg, const std::string& failure) {
    const std::string base = cfg.output_dir + "/" + cfg.tag;
    TFile file((base + ".root").c_str(), result.params.valid ? "UPDATE" : "RECREATE");
    std::ofstream parameters(base + "_params.csv", result.params.valid ? std::ios::app : std::ios::out);
    if (!result.params.valid) parameters << "parameter,value\nellipse_valid,0\nmcd_valid,0\n";
    parameters << "diagnostic_ellipse_valid," << (failure.empty() ? 1 : 0)
               << "\ndiagnostic_fallback_decorrelation," << (failure.empty() ? 0 : 1) << "\n";
    TH2D all("h_mass_compare_all", "Before mass selection;m_{inv} [GeV];M_{miss} [GeV]",
        cfg.n_mpi0_bins, cfg.mpi0_min, cfg.mpi0_max, cfg.n_mmiss_bins, cfg.mmiss_min, cfg.mmiss_max);
    all.SetDirectory(nullptr); all.Sumw2(); all.SetStats(0);
    auto make = [&](const char* name) {
        auto* h = static_cast<TH2D*>(all.Clone(name)); h->SetDirectory(nullptr);
        return std::unique_ptr<TH2D>(h);
    };
    auto decor = make("h_mass_compare_decorrelation");
    auto ellipse = make("h_mass_compare_ellipse");
    auto mcd = make("h_mass_compare_mcd");
    double outside_weight = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        if (!std::isfinite(p.mpi0) || !std::isfinite(p.mmiss) || !std::isfinite(p.weight)) continue;
        if (p.mpi0 < cfg.mpi0_min || p.mpi0 >= cfg.mpi0_max ||
            p.mmiss < cfg.mmiss_min || p.mmiss >= cfg.mmiss_max) {
            outside_weight += p.weight; continue;
        }
        all.Fill(p.mpi0, p.mmiss, p.weight);
        if (i < decorrelation.size() && decorrelation[i]) decor->Fill(p.mpi0,p.mmiss,p.weight);
        if (i < result.pass_ellipse.size() && result.pass_ellipse[i]) ellipse->Fill(p.mpi0,p.mmiss,p.weight);
        if (i < result.pass_mcd.size() && result.pass_mcd[i]) mcd->Fill(p.mpi0,p.mmiss,p.weight);
    }
    parameters << "diagnostic_outside_window_weight," << outside_weight << "\n";
    TCanvas canvas((cfg.tag + "_comparison").c_str(), "Mass selection comparison", 1600, 1200);
    canvas.Divide(2,2);
    const bool good = failure.empty();
    auto* selected = good ? ellipse.get() : decor.get();
    std::vector<std::unique_ptr<TH1D>> projections;
    for (int pad = 1; pad <= 2; ++pad) {
        canvas.cd(pad); gPad->SetRightMargin(0.15); gPad->SetTopMargin(0.14);
        TH2D* background = pad == 1 ? &all : selected;
        background->SetMinimum(0); background->SetMaximum(std::max(1.0, all.GetMaximum()));
        background->SetTitle(pad == 1 ? "Before mass selection" :
            (good ? "After ellipse (diagnostic)" : "After de-correlation (ellipse failed)"));
        background->Draw("COLZ");
        ellipse->SetLineColor(kMagenta+2); mcd->SetLineColor(kGreen+2);
        auto contour = [&](double mx, double my, double xx, double xy, double yy,
                           double det, double d2, int color, int style) {
            if (!(det > 0.0) || !std::isfinite(det)) return;
            nps2d::detail::CovModel model;
            model.mean_x=mx; model.mean_y=my; model.cov_xx=xx; model.cov_xy=xy; model.cov_yy=yy; model.det=det;
            auto* line = nps2d::detail::make_ellipse_line(model,d2,color,style,2);
            if (line) line->Draw("SAME");
        };
        const auto& p = result.params;
        if (p.ellipse_valid) contour(p.mean_mpi0,p.mean_mmiss,p.cov_mpi0_mpi0,p.cov_mpi0_mmiss,
            p.cov_mmiss_mmiss,p.cov_det,p.ellipse_d2_cut,kMagenta+2,good ? 1 : 2);
        if (p.mcd_valid) contour(p.mcd_mean_mpi0,p.mcd_mean_mmiss,p.mcd_cov_mpi0_mpi0,p.mcd_cov_mpi0_mmiss,
            p.mcd_cov_mmiss_mmiss,p.mcd_det,p.mcd_d2_cut,kGreen+2,1);
        double x0 = std::max(cfg.mpi0_min, 0.135 + (0.938-cfg.mmiss_max)/31.95);
        double x1 = std::min(cfg.mpi0_max, 0.135 + (0.938-cfg.mmiss_min)/31.95);
        auto* line = new TLine(x0, 0.938-31.95*(x0-0.135), x1, 0.938-31.95*(x1-0.135));
        line->SetLineStyle(2); line->SetLineWidth(2); line->Draw();
        auto* legend = new TLegend(0.52, 0.65, 0.83, 0.84);
        legend->SetTextSize(0.025); legend->SetFillStyle(0); legend->SetBorderSize(0);
        legend->AddEntry(ellipse.get(), good ? "Ellipse" : "Ellipse FAILED (stored)", "l");
        legend->AddEntry(mcd.get(), "MCD", "l"); legend->Draw();
        TLatex text; text.SetNDC(); text.SetTextSize(0.028);
        text.DrawLatex(0.12, 0.91, "Dashed: reference slope -31.95 (no new cut)");
    }
    for (int axis = 0; axis < 2; ++axis) {
        canvas.cd(axis + 3);
        auto* legend = new TLegend(0.48,0.70,0.89,0.89);
        legend->SetTextSize(0.026); legend->SetBorderSize(0); legend->SetFillStyle(0);
        int i = 0;
        const char* labels[] = {"Before (weighted)", "De-correlation", good ? "Ellipse" : "Ellipse FAILED (stored selector)", "MCD"};
        const int colors[] = {kBlack,kBlue,kMagenta+2,kGreen+2};
        std::vector<TH1D*> drawn;
        double maximum = 0.0;
        for (TH2D* h : {&all,decor.get(),ellipse.get(),mcd.get()}) {
            const std::string name = std::string(h->GetName()) + (axis ? "_y" : "_x");
            TH1D* p = axis ? h->ProjectionY(name.c_str()) : h->ProjectionX(name.c_str());
            p->SetDirectory(nullptr); p->SetLineColor(colors[i]); p->SetLineWidth(2); p->SetStats(0);
            p->SetTitle(axis ? "Missing-mass projection;M_{miss} [GeV];Weighted counts" : "Invariant-mass projection;m_{inv} [GeV];Weighted counts");
            maximum = std::max(maximum, p->GetMaximum()); drawn.push_back(p);
            legend->AddEntry(p, labels[i++], "l"); projections.emplace_back(p);
        }
        drawn.front()->SetMaximum(1.2*std::max(1.0,maximum));
        for (std::size_t j=0; j<drawn.size(); ++j) drawn[j]->Draw(j ? "HIST SAME" : "HIST");
        legend->Draw();
        if (!good) { TLatex text; text.SetNDC(); text.SetTextSize(0.027); text.SetTextColor(kRed+1);
            text.DrawLatex(0.12,0.53,failure.c_str()); }
    }
    canvas.SaveAs((base + ".png").c_str()); canvas.SaveAs((base + ".pdf").c_str());
    file.cd();
    for (TH2D* h : {&all,decor.get(),ellipse.get(),mcd.get()}) h->Write(h->GetName(),TObject::kOverwrite);
    TNamed status("diagnostic_mass_cut_status", good ? "ellipse qualified" : (failure+"; fallback=de-correlation").c_str());
    status.Write(status.GetName(),TObject::kOverwrite);
    canvas.Write("mass_cut_comparison",TObject::kOverwrite);
}
} // namespace npsplot
#endif
