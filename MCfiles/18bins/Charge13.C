#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>
#include <TMath.h>

void plot_hPurity_withErrors() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve hPurity histogram
    TH1F* h = (TH1F*)f->Get("hPurity");
    if (!h) {
        std::cerr << "Could not find hPurity!" << std::endl;
        f->Close();
        return;
    }

    // Calculate binomial errors for purity
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double purity = h->GetBinContent(i); // Purity value (between 0 and 1)
        double N = h->GetBinContent(i) * 1000; // Assume N_total (trials) is scaled for demo
        // Binomial error: sqrt(p * (1-p) / N)
        if (N > 0 && purity >= 0 && purity <= 1) {
            double error = TMath::Sqrt(purity * (1 - purity) / N);
            h->SetBinError(i, error);
        } else {
            h->SetBinError(i, 0); // Set zero error if invalid
        }
    }

    // Style
    h->SetStats(0);
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);    // Marker style for points
    h->SetMarkerColor(kBlue);
    h->SetMarkerSize(1.0);
    
    // Set axis labels and make X-axis title bold
    h->SetTitle(";#Delta#phi_{reco};Purity");

    h->GetXaxis()->SetTitleFont(62);
    h->GetXaxis()->SetTitleSize(0.04);

    h->GetYaxis()->SetTitleFont(62);
    h->GetYaxis()->SetTitleSize(0.05);

    // Create Canvas
    TCanvas *c = new TCanvas("c", "hPurity", 800, 600);

    // Draw histogram with error bars
    h->Draw("E1");   // Draw points with error bars and connecting line

    // Add CMS + labels
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

    // Update and save plot
    c->Update();
    c->SaveAs("hPurity_withErrors.pdf");

    f->Close();
}
