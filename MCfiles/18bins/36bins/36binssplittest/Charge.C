#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>

void plot_two_histograms() {
    // Open first ROOT file (pseudo-data)
    TFile *f1 = TFile::Open("DDS10_pseudodata_25_positron.root");
    if (!f1 || f1->IsZombie()) {
        std::cerr << "Error opening DDS10_pseudodata_25.root!" << std::endl;
        return;
    }

    // Open second ROOT file (response)
    TFile *f2 = TFile::Open("DDS10_response_75_positron.root");
    if (!f2 || f2->IsZombie()) {
        std::cerr << "Error opening DDS10_response_75.root!" << std::endl;
        f1->Close();
        return;
    }

    // Retrieve histograms
    TH1F* hgen = (TH1F*)f1->Get("gen");
    TH1F* hUnfolded = (TH1F*)f2->Get("hUnfoldedPseudodataPurityCorrectedOverAcceptance");

    if (!hgen || !hUnfolded) {
        std::cerr << "Could not find histograms!" << std::endl;
        f1->Close();
        f2->Close();
        return;
    }

    // Style histograms
    hgen->SetStats(0);
    hgen->SetLineColor(kRed);
    hgen->SetLineWidth(2);
    hgen->SetMarkerStyle(21);
    hgen->SetMarkerColor(kRed);
    hgen->SetMarkerSize(1.0);
    hgen->SetTitle("");
    hgen->SetMinimum(0);

    hUnfolded->SetStats(0);
    hUnfolded->SetLineColor(kBlue+1);
    hUnfolded->SetLineWidth(2);
    hUnfolded->SetMarkerStyle(20);
    hUnfolded->SetMarkerColor(kBlue+1);
    hUnfolded->SetMarkerSize(1.0);
    hUnfolded->SetTitle("");
    hUnfolded->SetMinimum(0);

    // Set axis titles
    hgen->GetXaxis()->SetTitleFont(62);
    hgen->GetXaxis()->SetTitleSize(0.04);
    hgen->GetYaxis()->SetTitleFont(62);
    hgen->GetYaxis()->SetTitleSize(0.04);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "gen vs Unfolded", 800, 600);

    // Draw the histogram with the larger maximum first
    double max_gen = hgen->GetMaximum();
    double max_unfold = hUnfolded->GetMaximum();
    
    if (max_gen > max_unfold) {
        hgen->Draw("E P");
        hUnfolded->Draw("E P SAME");
    } else {
        hUnfolded->Draw("E P");
        hgen->Draw("E P SAME");
    }

    // Legend
    TLegend* leg = new TLegend(0.60, 0.70, 0.88, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hgen, "generator Superchic (25%)", "lp");
    leg->AddEntry(hUnfolded, "Unfolded PseudoData (25%)", "lp");
    leg->Draw();

    // CMS-style labels
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(62);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.10, 0.90, "CMS");
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.16, 0.91, "#it{work in progress}");
    latex.SetTextSize(0.038);
    latex.SetTextAlign(33);
    latex.DrawLatex(0.89, 0.935, "PbPb 2018 #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // ---- Compute chi^2 / ndf manually ----
    double chi2 = 0;
    int ndf = 0;

    int nbins = hgen->GetNbinsX();
    for (int i = 1; i <= nbins; i++) {
        double y1 = hgen->GetBinContent(i);
        double y2 = hUnfolded->GetBinContent(i);
        double e1 = hgen->GetBinError(i);
        double e2 = hUnfolded->GetBinError(i);

        double err2 = e1*e1 + e2*e2; // combine errors in quadrature
        if (err2 > 0) {
            chi2 += (y1 - y2)*(y1 - y2)/err2;
            ndf++;
        }
    }

    double chi2ndf = chi2 / ndf;

    // p-value using ROOT's Chi2 probability
    double pvalue = TMath::Prob(chi2, ndf);

    // Draw chi2/ndf and p-value on canvas in "chi2/ndf = X/Y = Z" format
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.SetTextAlign(12); // left-aligned

    latex.DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf = %.2f/%d = %.3f", chi2, ndf, chi2ndf));
    latex.DrawLatex(0.15, 0.75, Form("p-value = %.3f", pvalue));


    // Save output
    c1->SaveAs("gen_vs_Unfolded_chi2ndf.pdf");

    // Close files
    f1->Close();
    f2->Close();
}

