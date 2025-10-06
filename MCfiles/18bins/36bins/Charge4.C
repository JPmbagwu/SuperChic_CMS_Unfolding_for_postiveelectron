#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>

void compare_corrected_vs_forwardfold() {
    // Open the ROOT file
    TFile *f = TFile::Open("FDC36enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) return;

    // Get histograms
    TH1F* hCorrectedData = (TH1F*)f->Get("hDeltaPhiPurityCorrected");
    TH1D* hForwardFolded = (TH1D*)f->Get("measured;1");

    if (!hCorrectedData || !hForwardFolded) {
        std::cerr << "Error: Histogram(s) not found!" << std::endl;
        f->Close();
        return;
    }

    // Style settings
    hCorrectedData->SetStats(0);
    hForwardFolded->SetStats(0);

    hCorrectedData->SetMarkerStyle(20);
    hCorrectedData->SetMarkerSize(1.0);
    hCorrectedData->SetMarkerColor(kRed);
    hCorrectedData->SetLineColor(kRed);
    hCorrectedData->SetLineWidth(2);

    hForwardFolded->SetMarkerStyle(24);
    hForwardFolded->SetMarkerSize(1.0);
    hForwardFolded->SetMarkerColor(kGreen + 2);
    hForwardFolded->SetLineColor(kGreen + 2);
    hForwardFolded->SetLineStyle(2);
    hForwardFolded->SetLineWidth(2);

    // Draw
    TCanvas *c = new TCanvas("c", "2018 Data vs Forward Folded", 800, 600);
    hCorrectedData->SetTitle(";#Delta#phi = #phi_{ee} - #phi_{e^{-}};Entries");
    hCorrectedData->Draw("PE");
    hForwardFolded->Draw("PE SAME");

    // Add legend
    TLegend *legend = new TLegend(0.65, 0.35, 0.88, 0.50);
    legend->AddEntry(hCorrectedData, "Data", "lep");
    legend->AddEntry(hForwardFolded, "Forward Folded Data", "lep");
    legend->Draw();

    // --- Calculate Chi2/ndf and p-value ---
    int ndf = 36; // explicitly set
    double chi2_ndf = hCorrectedData->Chi2Test(hForwardFolded, "UW CHI2"); // gives chi2/ndf
    double chi2_total = chi2_ndf * ndf;
    double pval = TMath::Prob(chi2_total, ndf);

    // Display Chi2/ndf and p-value on plot
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.SetTextAlign(11);
    latex.DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf = %.2f/%d = %.3f", chi2_total, ndf, chi2_ndf));
//    latex.DrawLatex(0.15, 0.80, Form("p-value = %.3f", pval));

    // CMS labels
    latex.SetTextFont(62);
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.10, 0.93, "CMS");

    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.16, 0.93, "#it{work in progress}");

    latex.SetTextSize(0.038);
    latex.SetTextAlign(33);
    latex.DrawLatex(0.89, 0.935, "PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // Save
    c->SaveAs("2018Data_vs_ForwardFold_withChi2.pdf");
    f->Close();
}

