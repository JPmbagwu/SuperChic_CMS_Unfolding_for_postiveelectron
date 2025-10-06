#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>

void plot_hDeltaPhiRaw() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histogram
    TH1F* h = (TH1F*)f->Get("hDeltaPhiRaw;1");
    if (!h) {
        std::cerr << "Could not find hDeltaPhiRaw!" << std::endl;
        f->Close();
        return;
    }

    // Remove stats box
    h->SetStats(0);

    // Optional: axis titles
    h->SetTitle(";#Delta#phi_{reco};Entries");
    h->GetXaxis()->SetTitleFont(62);
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleFont(62);
    h->GetYaxis()->SetTitleSize(0.04);

    // Canvas
    TCanvas *c = new TCanvas("c", "hDeltaPhiRaw", 800, 600);

    // Draw histogram as-is
    h->Draw();

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
    c->SaveAs("hDeltaPhiRaw.pdf");

    f->Close();
}

