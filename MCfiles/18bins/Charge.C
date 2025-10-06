#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>

void compare_split_histograms() {
    // Open the ROOT file
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve the histograms
    TH1F* hSplit1 = (TH1F*)f->Get("hrecoTrue_split1");
    TH1F* hSplit2 = (TH1F*)f->Get("hrecoTrue_split2");

    if (!hSplit1 || !hSplit2) {
        std::cerr << "Could not find one or both split histograms!" << std::endl;
        f->Close();
        return;
    }

    // Turn off stats box
    hSplit1->SetStats(0);
    hSplit2->SetStats(0);

    // Style settings
    hSplit1->SetMarkerStyle(20);  // filled circle
    hSplit1->SetMarkerSize(1.0);
    hSplit1->SetMarkerColor(kBlue);
    hSplit1->SetLineColor(kBlue);
    hSplit1->SetLineWidth(2);

    hSplit2->SetMarkerStyle(25);  // open square
    hSplit2->SetMarkerSize(1.0);
    hSplit2->SetMarkerColor(kRed);
    hSplit2->SetLineColor(kRed);
    hSplit2->SetLineWidth(2);
    hSplit2->SetLineStyle(2);     // dashed line

    // Canvas
    TCanvas *c = new TCanvas("c", "hrecoTrue Split Test Comparison", 800, 600);
    hSplit1->SetTitle(";#Delta#phi = #phi_{ee} - #phi_{e+};Entries");
    hSplit1->SetMaximum(std::max(hSplit1->GetMaximum(), hSplit2->GetMaximum()) * 1.2);

    hSplit1->Draw("PE");
    hSplit2->Draw("PE SAME");

    // Legend
    TLegend *legend = new TLegend(0.65, 0.35, 0.88, 0.50);
    legend->AddEntry(hSplit1, "hrecoTrue_split1", "lep");
    legend->AddEntry(hSplit2, "hrecoTrue_split2", "lep");
    legend->Draw();

    // -- Chi2 and P-value computation --
    double chi2 = 0;
    int ndf = 0;
    for (int i = 1; i <= hSplit1->GetNbinsX(); ++i) {
        double x1 = hSplit1->GetBinContent(i);
        double x2 = hSplit2->GetBinContent(i);
        double e1 = hSplit1->GetBinError(i);
        double e2 = hSplit2->GetBinError(i);
        double err2 = e1 * e1 + e2 * e2;

        if (err2 > 0) {
            chi2 += (x1 - x2) * (x1 - x2) / err2;
            ++ndf;
        }
    }

    double pvalue = TMath::Prob(chi2, ndf);

    // -- TLatex annotations --
    TLatex latex;
    latex.SetNDC();

    // CMS + Work in progress
    latex.SetTextFont(62);
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.10, 0.93, "CMS");

    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.16, 0.93, "#it{work in progress}");

    // PbPb 2018
    latex.SetTextSize(0.038);
    latex.SetTextAlign(33);
    latex.DrawLatex(0.89, 0.935, "PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // Chi2/ndf and p-value
    latex.SetTextAlign(13);
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.14, 0.88, Form("#chi^{2}/ndf = %.2f / %d", chi2, ndf));
    latex.DrawLatex(0.14, 0.84, Form("p-value = %.4f", pvalue));

    // Save plot
    c->Update();
    c->SaveAs("Superchic2018hrecoTrue_split_test_comparison.pdf");

    f->Close();
}

