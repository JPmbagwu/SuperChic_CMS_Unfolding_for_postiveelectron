#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>
#include <cmath>

void compare_corrected_vs_measured() {
    // Open the ROOT file
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get histograms
    TH1D* hMeasured = (TH1D*)f->Get("measured;1");
    TH1F* hCorrected = (TH1F*)f->Get("hDeltaPhiPurityCorrected");

    if (!hMeasured || !hCorrected) {
        std::cerr << "Could not find one or both histograms!" << std::endl;
        f->Close();
        return;
    }

    // Style settings
    hMeasured->SetStats(0);
    hCorrected->SetStats(0);

    hMeasured->SetMarkerStyle(20);
    hMeasured->SetMarkerSize(1.0);
    hMeasured->SetMarkerColor(kBlue);
    hMeasured->SetLineColor(kBlue);
    hMeasured->SetLineWidth(2);

    hCorrected->SetMarkerStyle(24);
    hCorrected->SetMarkerSize(1.0);
    hCorrected->SetMarkerColor(kRed);
    hCorrected->SetLineColor(kRed);
    hCorrected->SetLineWidth(2);
    hCorrected->SetLineStyle(2);

    // Canvas for overlay
    TCanvas *c1 = new TCanvas("c1", "Compare Corrected and Measured", 800, 600);
    hMeasured->SetTitle(";#Delta #phi = #phi_{ee} - #phi_{e-};Entries");

    hMeasured->Draw("PE");
    hCorrected->Draw("PE SAME");

    TLegend *legend = new TLegend(0.65, 0.68, 0.88, 0.83);
    legend->AddEntry(hMeasured, "Forward Fold", "lep");
    legend->AddEntry(hCorrected, "2018 Data", "lep");
    legend->Draw();

    // Add TLatex labels
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(62);
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.10, 0.93, "CMS");

    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.16, 0.93, "#it{work in progress}");

    latex.SetTextSize(0.038);
    latex.SetTextAlign(33);
    latex.DrawLatex(0.89, 0.935, "PbPb 2018 #sqrt{#it{s}_{NN}} = 5.02 TeV");

    latex.SetTextAlign(31);
    latex.SetTextSize(0.030);  // This line is optional unless more labels follow

    c1->Update();
    c1->SaveAs("corrected_vs_measured.pdf");

    // Difference plot
    TH1D* hDiff = (TH1D*)hMeasured->Clone("hDiff");
    hDiff->SetTitle(" ;#Delta#phi = #phi_{ee} - #phi_{e+};Difference");
    hDiff->Reset();

    int nBins = hMeasured->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double val1 = hMeasured->GetBinContent(i);
        double err1 = hMeasured->GetBinError(i);
        double val2 = hCorrected->GetBinContent(i);
        double err2 = hCorrected->GetBinError(i);

        double diff = val1 - val2;
        double errDiff = std::sqrt(err1 * err1 + err2 * err2);

        hDiff->SetBinContent(i, diff);
        hDiff->SetBinError(i, errDiff);
    }

    hDiff->SetMarkerStyle(21);
    hDiff->SetMarkerSize(1.2);
    hDiff->SetMarkerColor(kBlack);
    hDiff->SetLineColor(kBlack);
    hDiff->SetStats(0);

    TCanvas *c2 = new TCanvas("c2", "Difference Histogram", 800, 600);
    hDiff->Draw("PE");
    c2->Update();
    c2->SaveAs("correctedData2018_vs_measured2018Data_diff.pdf");

    f->Close();
}

