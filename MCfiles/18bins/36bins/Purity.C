#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>
#include <cmath> // for sqrt

void plot_purity() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC36enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histograms for purity
    TH1F* hrecoTrue = (TH1F*)f->Get("hrecoTrue");
    TH1F* hrecoFake = (TH1F*)f->Get("hrecoFake");

    if (!hrecoTrue || !hrecoFake) {
        std::cerr << "Could not find hrecoTrue or hrecoFake!" << std::endl;
        f->Close();
        return;
    }

    // Create histogram for purity
    TH1F* hPurity = (TH1F*)hrecoTrue->Clone("hPurity");
    hPurity->SetTitle(";#Delta#phi;Purity");
    hPurity->SetStats(0);                 // Disable stats box
    hPurity->GetYaxis()->SetRangeUser(0, 1); // Purity between 0 and 1
    hPurity->SetLineColor(kBlue+2);
    hPurity->SetMarkerColor(kBlue+2);
    hPurity->SetMarkerStyle(20);
    hPurity->SetMarkerSize(1.0);

    // ----------- Calculate purity -----------
    for (int i = 1; i <= hPurity->GetNbinsX(); ++i) {
        double nTrue = hrecoTrue->GetBinContent(i);
        double nFake = hrecoFake->GetBinContent(i);
        double denom = nTrue + nFake;

        if (denom > 0) {
            double purity = nTrue / denom;
            double error = sqrt(purity * (1 - purity) / denom); // binomial error
            hPurity->SetBinContent(i, purity);
            hPurity->SetBinError(i, error);
        } else {
            hPurity->SetBinContent(i, 0);
            hPurity->SetBinError(i, 0);
        }
    }

    // ----------- Plot purity -----------
    TCanvas *cPur = new TCanvas("cPur", "Purity", 800, 600);
    hPurity->Draw("E P");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(62);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.10, 0.93, "CMS");
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.16, 0.93, "#it{work in progress}");
    latex.SetTextSize(0.038);
    latex.SetTextAlign(33);
    latex.DrawLatex(0.89, 0.935, "PbPb 2018 #sqrt{#it{s}_{NN}} = 5.02 TeV");

    cPur->SaveAs("PurityAAA.pdf");

    // Close file
    f->Close();
}
