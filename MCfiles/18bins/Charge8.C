#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>

void plot_UnfoldData_div_hAcceptance() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histograms
    TH1F* hUnfold = (TH1F*)f->Get("UnfoldData2018;1");
    TH1F* hAcc    = (TH1F*)f->Get("hAcceptance;1");

    if (!hUnfold || !hAcc) {
        std::cerr << "Could not find UnfoldData2018 or hAcceptance!" << std::endl;
        f->Close();
        return;
    }

    // Clone to avoid modifying the original
    TH1F* hRatio = (TH1F*)hUnfold->Clone("hRatio");
    hRatio->Divide(hAcc);

    // Style
    hRatio->SetStats(0);
    hRatio->SetLineColor(kBlue+1);
    hRatio->SetLineWidth(2);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(kBlue+1);
    hRatio->SetMarkerSize(1.0);

    hRatio->SetTitle(";#Delta#phi_{reco}; UnfoldData2018 / hAcceptance");
    hRatio->GetXaxis()->SetTitleFont(62);
    hRatio->GetXaxis()->SetTitleSize(0.04);
    hRatio->GetYaxis()->SetTitleFont(62);
    hRatio->GetYaxis()->SetTitleSize(0.04);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "UnfoldData2018 / hAcceptance", 800, 600);
    hRatio->Draw("E P");

    // CMS style labels
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
    latex.DrawLatex(0.89, 0.935, "PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // Save output
    c1->SaveAs("UnfoldData2018_div_hAcceptance.pdf");

    // Close file
    f->Close();
}
