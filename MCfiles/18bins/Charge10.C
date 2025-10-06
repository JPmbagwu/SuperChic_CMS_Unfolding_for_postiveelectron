#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>

void plot_hDeltaPhiPurityCorrected_withBlackErrors() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histogram
    TH1F* h = (TH1F*)f->Get("hDeltaPhiPurityCorrected;1");
    if (!h) {
        std::cerr << "Could not find hDeltaPhiPurityCorrected!" << std::endl;
        f->Close();
        return;
    }

    // Style for histogram line & markers (blue)
    h->SetStats(0);
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(kBlue);
    h->SetMarkerSize(1.0);

    // Axis titles bold & font size
    h->SetTitle(";#Delta#phi_{reco};Entries");
    h->GetXaxis()->SetTitleFont(62);
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleFont(62);
    h->GetYaxis()->SetTitleSize(0.04);

    // Force y-axis to start at 0
    h->SetMinimum(0);

    // Canvas
    TCanvas *c = new TCanvas("c", "hDeltaPhiPurityCorrected", 800, 600);

    // Draw line and markers first (blue)
    h->Draw("HIST P");

    // Draw error bars only in black
    h->SetLineColor(kBlack);
    h->SetMarkerColor(kBlack);
    h->Draw("E SAME");

    // Draw CMS + labels
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
    latex.DrawLatex(0.89, 0.935, "PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    c->Update();
    c->SaveAs("hDeltaPhiPurityCorrected_blackErrors10.pdf");

    f->Close();
}
