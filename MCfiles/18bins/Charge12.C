#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <iostream>

void plot_hResponseMatrix_only() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve hResponse histogram (2D)
    TH2F* hResponse = (TH2F*)f->Get("hResponseMatrix;1");
    if (!hResponse) {
        std::cerr << "Could not find hResponseMatrix!" << std::endl;
        f->Close();
        return;
    }

    // Set axis titles
    hResponse->SetXTitle("#Delta#phi_{gen}");
    hResponse->SetYTitle("#Delta#phi_{reco}");

    // Canvas
    TCanvas *c = new TCanvas("cResponse", "Response Matrix", 800, 600);

    // Style
    gStyle->SetOptStat(1111);  // Show stats box with entries, mean, rms, etc.

    hResponse->Draw("COLZ");

    gPad->Update();

    // Move and resize stats box
    TPaveStats *stat = (TPaveStats*)hResponse->FindObject("stats");
    if (stat) {
        stat->SetX1NDC(0.60);
        stat->SetX2NDC(0.90);
        stat->SetY1NDC(0.28);
        stat->SetY2NDC(0.42);
        stat->SetTextSize(0.03);
        gPad->Modified();
        gPad->Update();
    }

    // Draw CMS + progress + PbPb labels (top area)
    TLatex latex;
    latex.SetNDC();

    latex.SetTextFont(62);      // Bold
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);     // left bottom
    latex.DrawLatex(0.10, 0.93, "CMS");

    latex.SetTextFont(42);      // Normal
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.16, 0.93, "#it{work in progress}");

    latex.SetTextSize(0.038);
    latex.SetTextAlign(33);     // right top
    latex.DrawLatex(0.89, 0.935, "PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // Draw selection cuts text on top-left below CMS text
    latex.SetTextFont(62);      // Bold font
    latex.SetTextAlign(13);     // left bottom
    latex.SetTextSize(0.030);
    latex.DrawLatex(0.18, 0.80, "p_{T}_{ee} < 1.0 GeV");
//    latex.DrawLatex(0.18, 0.80, "p_{T}_{e} > 2.0 GeV");
    latex.DrawLatex(0.18, 0.75, "|#eta_{e}| < 2.4");
    latex.DrawLatex(0.19, 0.70, "7.0 < M_{ee} < 10.0 GeV");
    
    // Save plot
    c->SaveAs("hResponseMatrix_only.pdf");

    f->Close();
}
