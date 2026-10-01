#pragma once

// Ali thesis, printed p.142, Fig.5.3: component-dependent integrated responses
// link generated t' bins to reconstructed (t',phi) bins. The signed LT/TT
// matrices are not probabilities. Our additional count fractions below are
// explicitly conditional on selected MC, and never replace the fit response.
#include "xsec_analysis.h"
#include "xsec_plot_style.h"
#include <TColor.h>
#include <TDirectory.h>
#include <TExec.h>
#include <TBox.h>

inline void ExclPi0XSecAnalysis::make_migration_plots() {
    if (!cfg.diagnostics) return;
    const int ns = static_cast<int>(slices.size());
    const int nb = static_cast<int>(truth_moments.size());
    const int nr = static_cast<int>(migration_response.size());
    if (ns == 0 || nb == 0 || nr == 0) return;
    if (nr != ns * cfg.n_phi || nb < ns)
        die("Migration diagnostics: inconsistent response dimensions");
    for (const auto& row : migration_response)
        if (static_cast<int>(row.size()) != nb)
            die("Migration diagnostics: inconsistent truth-block dimensions");
    if (!fout || fout->IsZombie()) die("Migration diagnostics require an open output ROOT file");

    const auto slice_label = [&](int b) {
        const int it = b / (cfg.n_q2 * cfg.n_xb);
        const int iq = (b / cfg.n_xb) % cfg.n_q2;
        const int ix = b % cfg.n_xb;
        return "t" + std::to_string(it) + "/q" + std::to_string(iq) + "/x" + std::to_string(ix);
    };
    const auto truth_label = [&](int b) {
        if (b < ns) return slice_label(b);
        static const char* guards[] = {"t'<min", "t'>max", "Q^{2}<min", "Q^{2}>max", "x_{B}<min", "x_{B}>max"};
        return "b" + std::to_string(b) + ":" +
               (b - ns < 6 ? std::string(guards[b - ns]) : std::string("guard"));
    };
    const auto display_truth_label = [&](int b) {
        if (b < ns && cfg.n_q2==1 && cfg.n_xb==1 &&
            tprime_edges.size()==static_cast<size_t>(cfg.n_tprime+1)) {
            const int it=b;
            return std::string(Form("b%d:%.2f-%.2f",b,-tprime_edges[it+1],-tprime_edges[it]));
        }
        if (b>=ns) {
            static const char* guard_short[]={"t'<min","t'>max","Q2<min","Q2>max","xB<min","xB>max"};
            return "g"+std::to_string(b-ns)+":"+guard_short[b-ns];
        }
        return "b"+std::to_string(b);
    };
    const auto display_reco_label = [&](int s) {
        if(cfg.n_q2==1 && cfg.n_xb==1 &&
           tprime_edges.size()==static_cast<size_t>(cfg.n_tprime+1))
            return std::string(Form("s%d:%.2f-%.2f",s,-tprime_edges[s+1],-tprime_edges[s]));
        return "s"+std::to_string(s);
    };
    const auto make_hist = [&](const char* name, const char* title, int ny) {
        auto h = std::make_unique<TH2D>(name, title, nb, 0., nb, ny, 0., ny);
        h->SetDirectory(nullptr);
        h->SetStats(false);
        for (int b = 0; b < nb; ++b) h->GetXaxis()->SetBinLabel(b + 1, truth_label(b).c_str());
        return h;
    };

    const char* component_names[] = {"U", "LT", "TT"};
    std::array<std::unique_ptr<TH2D>, 3> response;
    for (int a = 0; a < 3; ++a) {
        const std::string name = std::string("response_") + component_names[a];
        const std::string title = std::string(component_names[a]) +
            " integrated response;Generated truth block;Reconstructed row r;Yield per (#mub/MeV^{2})";
        response[a] = make_hist(name.c_str(), title.c_str(), nr);
    }
    auto counts = make_hist("selected_counts",
        "Selected exclusive MC counts;Generated truth block;Reconstructed slice;Events (unweighted)", ns);
    auto p_reco = make_hist("p_reco_given_truth_selected",
        "P(reco slice | truth block, selected);Generated truth block;Reconstructed slice;Conditional fraction", ns);
    auto p_truth = make_hist("p_truth_given_reco_selected",
        "P(truth block | reco slice, selected);Generated truth block;Reconstructed slice;Conditional fraction", ns);
    std::vector<double> truth_counts(nb, 0.), reco_counts(ns, 0.);
    for (int r = 0; r < nr; ++r) {
        const int s = r / cfg.n_phi;
        for (int b = 0; b < nb; ++b) {
            const auto& cell = migration_response[r][b];
            for (int a = 0; a < 3; ++a)
                response[a]->SetBinContent(b + 1, r + 1, cell.basis[a]);
            const double n = static_cast<double>(cell.events);
            counts->AddBinContent(counts->GetBin(b + 1, s + 1), n);
            truth_counts[b] += n;
            reco_counts[s] += n;
        }
    }
    // Aggregate phi counts first, then normalize using ALL truth blocks and
    // reconstructed slices, including feed-in guards and any fit-excluded rows.
    // Missing/rejected generated events are absent: neither fraction measures
    // acceptance efficiency or an absolute generated-to-detected probability.
    for (int s = 0; s < ns; ++s) {
        for (auto* h : {counts.get(), p_reco.get(), p_truth.get()})
            h->GetYaxis()->SetBinLabel(s + 1, slice_label(s).c_str());
        for (int b = 0; b < nb; ++b) {
            const double n = counts->GetBinContent(b + 1, s + 1);
            p_reco->SetBinContent(b + 1, s + 1, truth_counts[b] > 0. ? n / truth_counts[b] : 0.);
            p_truth->SetBinContent(b + 1, s + 1, reco_counts[s] > 0. ? n / reco_counts[s] : 0.);
        }
    }
    // Preserve original truth-block indices, including empty guard columns, in
    // ROOT. Drawing may suppress empty guards solely to make labels readable.
    {
        TDirectory::TContext restore_directory;
        TDirectory* directory = fout->GetDirectory("migration");
        if (!directory) directory = fout->mkdir("migration");
        if (!directory) die("Cannot create ROOT migration diagnostics directory");
        directory->cd();
        for (auto& h : response) h->Write();
        counts->Write(); p_reco->Write(); p_truth->Write();
        TObjString(
            "Reference: Ali thesis, printed p.142 Fig.5.3, https://www.osti.gov/biblio/1784736\n"
            "signed_tprime=t-tmin; slice=(it*n_q2+iq)*n_xb+ix; reco_row=slice*n_phi+ip; phi varies fastest\n"
            "truth_block=published slice, then six guards: tprime_below, tprime_above, q2_below, q2_above, xb_below, xb_above\n"
            "response_U/LT/TT=raw integrated event basis, yield per coefficient in microbarn/MeV^2; no probability normalization\n"
            "selected_counts=sum_phi cell.events; unweighted selected exclusive MC; includes all fit-excluded rows and guards\n"
            "p_reco_given_truth_selected=counts/column_sum; p_truth_given_reco_selected=counts/row_sum\n"
            "zero denominator=0 stored for unsupported slice/block, not a measured zero probability\n"
            "fractions depend on the sampled selected MC population; not efficiencies; not fitted-cross-section-weighted\n"
            "ROOT preserves all original truth blocks; plots hide only zero-event guard columns"
        ).Write("definitions");
    }
    if (!cfg.write_pdf && !cfg.write_png) return;
    const fs::path plot_dir = fs::path(cfg.out_dir) / "migration";
    fs::create_directories(plot_dir);

    std::vector<int> shown;
    for (int b = 0; b < nb; ++b)
        if (b < ns || truth_counts[b] > 0.) shown.push_back(b);
    const int nx = static_cast<int>(shown.size());
    const auto display_hist = [&](const TH2D& source, const std::string& name, bool phi_rows,
                                  double unit_scale=1.) {
        const int ny = source.GetNbinsY();
        auto h = std::make_unique<TH2D>(name.c_str(), source.GetTitle(), nx, 0., nx, ny, 0., ny);
        h->SetDirectory(nullptr); h->SetStats(false);
        h->GetXaxis()->SetTitle(source.GetXaxis()->GetTitle());
        h->GetYaxis()->SetTitle(source.GetYaxis()->GetTitle());
        h->GetZaxis()->SetTitle(unit_scale==1. ? source.GetZaxis()->GetTitle() :
                                "Weighted yield per (nb/GeV^{2})");
        for (int x = 0; x < nx; ++x) {
            h->GetXaxis()->SetBinLabel(x + 1, display_truth_label(shown[x]).c_str());
            for (int y = 0; y < ny; ++y)
                h->SetBinContent(x + 1, y + 1,unit_scale*source.GetBinContent(shown[x] + 1, y + 1));
        }
        // Slice labels mark groups of n_phi response rows. Their individual row
        // index and angular ordering remain explicit in the canvas annotation.
        if (phi_rows) {
            for (int s = 0; s < ns; ++s)
                h->GetYaxis()->SetBinLabel(s * cfg.n_phi + cfg.n_phi / 2 + 1,
                                           display_reco_label(s).c_str());
            h->GetYaxis()->SetTitle("Reconstructed slice (#phi rows within slice)");
        } else {
            for (int s = 0; s < ns; ++s)
                h->GetYaxis()->SetBinLabel(s + 1,display_reco_label(s).c_str());
        }
        h->GetXaxis()->LabelsOption("v");
        h->GetXaxis()->SetLabelSize(std::min(cfg.n_q2==1 && cfg.n_xb==1 ? .029 : .040,
                                             .95 / std::max(1, nx)));
        h->GetXaxis()->SetTitleOffset(3.0);
        h->GetYaxis()->SetLabelSize(std::min(cfg.n_q2==1 && cfg.n_xb==1 ? .027 : .040,
                                             .85 / std::max(1, ns)));
        h->GetYaxis()->SetTitleOffset(2.6);
        h->GetZaxis()->SetTitleOffset(1.6);
        for (TAxis* axis : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()}) axis->SetTitleSize(.032);
        h->GetZaxis()->SetLabelSize(.030);
        return h;
    };

    // ROOT re-paints pads during every PDF/PNG export. A palette TExec before
    // each COLZ object makes its own palette repeatable across those re-paints.
    // Literal color IDs avoid requiring compiled helper symbols in Cling.
    struct RestorePalette {
        std::vector<int> colors;
        int contours;
        RestorePalette() : contours(gStyle->GetNumberContours()) {
            for (int i = 0; i < gStyle->GetNumberOfColors(); ++i)
                colors.push_back(gStyle->GetColorPalette(i));
        }
        ~RestorePalette() {
            if (!colors.empty()) gStyle->SetPalette(static_cast<int>(colors.size()), colors.data());
            gStyle->SetNumberContours(contours);
        }
    } restore_palette;
    constexpr int ncolors = 101;
    double stops[] = {0., .5, 1.};
    double red[] = {.12, 1., .72}, green[] = {.30, 1., .08}, blue[] = {.72, 1., .12};
    const int first = TColor::CreateGradientColorTable(3, stops, red, green, blue, ncolors);
    std::ostringstream diverging;
    diverging << "{ Int_t colors[" << ncolors << "] = {";
    for (int i = 0; i < ncolors; ++i) diverging << (i ? "," : "") << first + i;
    diverging << "}; gStyle->SetPalette(" << ncolors << ", colors); gStyle->SetNumberContours(" << ncolors << "); }";
    const std::string sequential = "gStyle->SetPalette(57); gStyle->SetNumberContours(101);"; // ROOT kBird.
    const auto draw_matrix = [&](TH2D& h, TExec& palette, bool phi_rows,
                                 bool highlight_diagonal=false) {
        style_current_pad(.20, .21, .23, .13);
        h.Draw("AXIS");
        palette.Draw();
        h.Draw("COLZ SAME");
        TLine separator;
        separator.SetLineColor(kGray + 2); separator.SetLineStyle(3);
        if (phi_rows) {
            for (int s = 1; s < ns; ++s)
                separator.DrawLine(0., s * cfg.n_phi, nx, s * cfg.n_phi);
        }
        if (nx > ns) {
            separator.SetLineColor(kBlack); separator.SetLineStyle(2); separator.SetLineWidth(2);
            separator.DrawLine(ns, 0., ns, h.GetNbinsY());
        }
        if(highlight_diagonal) {
            TBox diagonal;
            diagonal.SetFillStyle(0);
            diagonal.SetLineColor(kBlack);
            diagonal.SetLineWidth(2);
            for(int s=0;s<ns;++s) diagonal.DrawBox(s,s,s+1,s+1);
        }
    };
    const auto canvas_note = [](TCanvas& canvas, const char* title, const char* line1, const char* line2) {
        canvas.cd();
        TLatex note; note.SetNDC(true); note.SetTextFont(42);
        note.SetTextSize(.023); note.DrawLatex(.025, .973, title);
        note.SetTextSize(.017); note.DrawLatex(.025, .041, line1); note.DrawLatex(.025, .018, line2);
    };

    {
        std::array<std::unique_ptr<TH2D>, 3> panels;
        std::array<std::unique_ptr<TExec>, 3> palettes;
        TCanvas canvas("c_migration_response", "Integrated migration responses", 2100, 1100);
        canvas.Divide(3, 1, .006, .065);
        for (int a = 0; a < 3; ++a) {
            // The fit and stored ROOT response stay in native sigcm units.
            // Only the displayed response is per nb/GeV^2 coefficient.
            panels[a] = display_hist(*response[a], std::string("display_response_") + component_names[a], true,1e-9);
            double extent = 0.;
            for (int x = 1; x <= nx; ++x) for (int y = 1; y <= nr; ++y)
                extent = std::max(extent, std::abs(panels[a]->GetBinContent(x, y)));
            if (extent == 0.) extent = 1.;
            panels[a]->SetMinimum(a == 0 ? 0. : -extent);
            panels[a]->SetMaximum(extent);
            palettes[a] = std::make_unique<TExec>((std::string("palette_response_") + component_names[a]).c_str(),
                                                 a == 0 ? sequential.c_str() : diverging.str().c_str());
            canvas.cd(a + 1); draw_matrix(*panels[a], *palettes[a], true);
        }
        canvas_note(canvas, "Integrated U/LT/TT response: yield per nb/GeV^{2} coefficient (not probabilities)",
            "s=(it*n_{Q2}+iq)*n_{xB}+ix; r=s*n_{#phi}+ip. b=published truth block; g=feed-in guard. Full mapping: migration_truth_blocks.csv.",
            "Signed t'=t-t_{min}; #phi varies fastest in rows. Dashed line starts guards; ROOT also retains empty guard columns.");
        canvas.Update();
        write_canvas_pdf_png(&canvas, (plot_dir / "response_components").string());
    }
    {
        auto panel=display_hist(*counts,"display_migration_support",false);
        double highest=0.;
        for(int x=1;x<=nx;++x) for(int y=1;y<=ns;++y) {
            const double value=std::log10(1.+panel->GetBinContent(x,y));
            panel->SetBinContent(x,y,value);
            highest=std::max(highest,value);
        }
        panel->SetTitle(";Generated truth block;Reconstructed slice;log_{10}(1+events)");
        panel->GetZaxis()->SetTitle("log_{10}(1+selected events)");
        panel->SetMinimum(0.);panel->SetMaximum(std::max(1.,highest));
        TExec palette("palette_migration_support",sequential.c_str());
        TCanvas canvas("c_migration_support","Selected MC migration support",1200,950);
        canvas.cd();draw_matrix(*panel,palette,false,true);
        canvas_note(canvas,"Selected exclusive MC: event support by truth and reconstructed bin",
            "Colors show log_{10}(1+selected events); white cells have zero selected events. Diagonal outlines mark corresponding bins.",
            "b/s labels give -t' ranges for single-Q2/xB grids; g=exterior feed-in. Lost events are absent; this is not an efficiency.");
        canvas.Update();
        write_canvas_pdf_png(&canvas,(plot_dir/"selected_migration_support").string());
    }
    {
        std::array<std::unique_ptr<TH2D>,2> panels;
        std::array<std::unique_ptr<TExec>,2> palettes;
        const TH2D* sources[]={p_reco.get(),p_truth.get()};
        TCanvas canvas("c_migration_fractions","Conditional selected-MC migration fractions",1750,1050);
        canvas.Divide(2,1,.008,.065);
        for(int i=0;i<2;++i) {
            panels[i]=display_hist(*sources[i],"display_migration_fraction_"+std::to_string(i),false);
            panels[i]->SetMinimum(0.);panels[i]->SetMaximum(1.);
            palettes[i]=std::make_unique<TExec>(("palette_migration_fraction_"+std::to_string(i)).c_str(),sequential.c_str());
            canvas.cd(i+1);draw_matrix(*panels[i],*palettes[i],false,true);
        }
        canvas_note(canvas,"Selected exclusive MC: conditional migration fractions (not efficiencies)",
            "Left: P(reco | truth, selected), each supported truth column sums to 1. Right: P(truth | reco, selected), each supported reco row sums to 1.",
            "Outlined diagonal means matching bin indices; other cells show migration or exterior feed-in. Zero denominators mean unsupported.");
        canvas.Update();
        write_canvas_pdf_png(&canvas, (plot_dir / "selected_migration_fractions").string());
    }
}
