#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>
#include <cmath> // for sqrt

void plot_acceptance() {
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

    // Create histogram for acceptance
    TH1F* hAcceptance = (TH1F*)hgenTrue->Clone("hAcceptance");
    hAcceptance->SetTitle(";#Delta#phi;Acceptance");
    hAcceptance->SetStats(0);                 // Disable stats box
    hAcceptance->GetYaxis()->SetRangeUser(0, 0.4); // y-axis from 0 to 0.4
    hAcceptance->SetLineColor(kGreen+2);
    hAcceptance->SetMarkerColor(kGreen+2);
    hAcceptance->SetMarkerStyle(20);
    hAcceptance->SetMarkerSize(1.0);


    // ----------- Calculate acceptance -----------
    for (int i = 1; i <= hAcceptance->GetNbinsX(); ++i) {
        double nTrue = hgenTrue->GetBinContent(i);
        double nMiss = hgenMiss->GetBinContent(i);
        double total = nTrue + nMiss;

        if (total > 0) {
            double acceptance = nTrue / total;
            double error = sqrt(acceptance * (1 - acceptance) / total); // binomial error
            hAcceptance->SetBinContent(i, acceptance);
            hAcceptance->SetBinError(i, error);
        } else {
            hAcceptance->SetBinContent(i, 0);
            hAcceptance->SetBinError(i, 0);
        }
    }

    // ----------- Plot acceptance -----------
    TCanvas *cAcc = new TCanvas("cAcc", "Acceptance", 800, 600);
    hAcceptance->Draw("E P");

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

    cAcc->SaveAs("AcceptanceAAA1.pdf");

    // Close file
    f->Close();
}
