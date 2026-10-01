#pragma once

// ROOT drawing helpers only. These functions never alter yields, acceptance, or fit parameters.
#include "xsec_physics.h"

static void style_current_pad(double left = 0.13, double right = 0.04, double bottom = 0.12, double top = 0.08, bool grid = false) {
    if (!gPad) return;
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetBottomMargin(bottom);
    gPad->SetTopMargin(top);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(grid, grid);
}

static void style_axis(TAxis* ax, double title_size = 0.045, double label_size = 0.04, double title_offset = 1.1) {
    if (!ax) return;
    ax->SetTitleSize(title_size);
    ax->SetLabelSize(label_size);
    ax->SetTitleOffset(title_offset);
}

static void style_hist_axes(TH1* h, double title_size = 0.045, double label_size = 0.04) {
    if (!h) return;
    style_axis(h->GetXaxis(), title_size, label_size, 1.05);
    style_axis(h->GetYaxis(), title_size, label_size, 1.25);
    h->SetTitleSize(title_size);
}

static void style_graph_axes(TGraphErrors* g, double title_size = 0.045, double label_size = 0.04) {
    if (!g) return;
    style_axis(g->GetXaxis(), title_size, label_size, 1.05);
    style_axis(g->GetYaxis(), title_size, label_size, 1.25);
}

static TLegend* make_compact_legend(double x1, double y1, double x2, double y2, double text_size = 0.035) {
    auto leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(text_size);
    return leg;
}

static void draw_pad_message(const std::string& title, const std::string& line1, const std::string& line2 = "") {
    if (!gPad) return;
    gPad->Clear();
    style_current_pad();
    TLatex note;
    note.SetNDC(true);
    note.SetTextFont(42);
    note.SetTextSize(0.052);
    note.DrawLatex(0.12, 0.78, title.c_str());
    note.SetTextSize(0.038);
    note.DrawLatex(0.12, 0.62, line1.c_str());
    if (!line2.empty()) note.DrawLatex(0.12, 0.52, line2.c_str());
}

