#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>

void plot_hgenTrue_hgenMiss() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC36enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histograms
    TH1F* hgenTrue = (TH1F*)f->Get("hgenTrue");
    TH1F* hgenMiss = (TH1F*)f->Get("hgenMiss");

    if (!hgenTrue || !hgenMiss) {
        std::cerr << "Could not find hgenTrue or hgenMiss!" << std::endl;
        f->Close();
        return;
    }

    // ----------- Style common axis settings function -----------
    auto styleAxis = [](TH1* h, const char* xtitle) {
        h->SetStats(0);
        h->SetTitle(Form(";%s;Entries", xtitle));
        h->GetXaxis()->SetTitleFont(62);
        h->GetXaxis()->SetTitleSize(0.04);
        h->GetXaxis()->SetLabelFont(62);
        h->GetXaxis()->SetLabelSize(0.035);

        h->GetYaxis()->SetTitleFont(62);
        h->GetYaxis()->SetTitleSize(0.04);
        h->GetYaxis()->SetLabelFont(62);
        h->GetYaxis()->SetLabelSize(0.035);
        h->GetYaxis()->SetTitle("Entries");
        h->SetMinimum(0); // start y-axis at 0
    };

    // ----------- Plot hgenTrue -----------
    hgenTrue->SetLineColor(kBlue);
    hgenTrue->SetLineWidth(2);
    hgenTrue->SetMarkerStyle(20);
    hgenTrue->SetMarkerColor(kBlue);
    hgenTrue->SetMarkerSize(1.0);
    styleAxis(hgenTrue, "#Delta#phi = #phi_{ee} - #phi_{e-}");

    TCanvas *c1 = new TCanvas("c1", "hgenTrue", 800, 600);
    hgenTrue->Draw("E P");

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextFont(62);
    latex1.SetTextSize(0.04);
    latex1.DrawLatex(0.10, 0.93, "CMS");
    latex1.SetTextFont(42);
    latex1.SetTextSize(0.035);
    latex1.DrawLatex(0.16, 0.93, "#it{work in progress}");
    latex1.SetTextSize(0.038);
    latex1.SetTextAlign(33);
    latex1.DrawLatex(0.89, 0.935, "PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    c1->SaveAs("hgenTrue.pdf");

    // ----------- Plot hgenMiss -----------
    hgenMiss->SetLineColor(kRed);
    hgenMiss->SetLineWidth(2);
    hgenMiss->SetMarkerStyle(21);
    hgenMiss->SetMarkerColor(kRed);
    hgenMiss->SetMarkerSize(1.0);
    styleAxis(hgenMiss, "#Delta#phi = #phi_{ee} - #phi_{e-}");

    TCanvas *c2 = new TCanvas("c2", "hgenMiss", 800, 600);
    hgenMiss->Draw("E P");

    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextFont(62);
    latex2.SetTextSize(0.04);
    latex2.DrawLatex(0.10, 0.93, "CMS");
    latex2.SetTextFont(42);
    latex2.SetTextSize(0.035);
    latex2.DrawLatex(0.16, 0.93, "#it{work in progress}");
    latex2.SetTextSize(0.038);
    latex2.SetTextAlign(33);
    latex2.DrawLatex(0.89, 0.935, "PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");

    c2->SaveAs("hgenMiss.pdf");

    // Close file
    f->Close();
}

