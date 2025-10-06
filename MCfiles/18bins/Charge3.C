#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>

void plot_hgenTrue_withBlackErrors() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histogram hgenTrue;1
    TH1F* h = (TH1F*)f->Get("hgenTrue;1");
    if (!h) {
        std::cerr << "Could not find hgenTrue!" << std::endl;
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
    h->SetTitle(";#Delta#phi_{gen};Entries");
    h->GetXaxis()->SetTitleFont(62);
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleFont(62);
    h->GetYaxis()->SetTitleSize(0.04);

    // Create canvas
    TCanvas *c = new TCanvas("c", "hgenTrue", 800, 600);

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
    latex.DrawLatex(0.89, 0.935, "PbPb 2018 #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // Update and save canvas
    c->Update();
    c->SaveAs("hgenTrue_blackErrors.pdf");

    f->Close();
}
