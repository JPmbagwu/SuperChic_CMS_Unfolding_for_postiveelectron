#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

void compare_unfolded_vs_genTrue() {
    TFile *f = TFile::Open("FDC36enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) return;

    TH1D* hUnfolded = (TH1D*)f->Get("UnfoldSCMC2018;1");
    TH1F* hGenTrue = (TH1F*)f->Get("hgenTrue");

    if (!hUnfolded || !hGenTrue) {
        f->Close();
        return;
    }

    hUnfolded->SetStats(0);
    hGenTrue->SetStats(0);

    // --- Style Unfolded histogram ---
    hUnfolded->SetMarkerStyle(20);
    hUnfolded->SetMarkerSize(1.0);
    hUnfolded->SetMarkerColor(kGreen+2);
    hUnfolded->SetLineColor(kGreen+2);
    hUnfolded->SetLineWidth(2);

    // --- Style GenTrue histogram ---
    hGenTrue->SetMarkerStyle(25);
    hGenTrue->SetMarkerSize(1.0);
    hGenTrue->SetMarkerColor(kMagenta);
    hGenTrue->SetLineColor(kMagenta);
    hGenTrue->SetLineStyle(2);
    hGenTrue->SetLineWidth(2);

    // --- Axis formatting ---
    hUnfolded->SetTitle(";#Delta#phi = #phi_{ee} - #phi_{e-};Entries"); // fixed φe- here
    hUnfolded->GetXaxis()->SetTitleFont(62);
    hUnfolded->GetXaxis()->SetTitleSize(0.04);
    hUnfolded->GetXaxis()->SetLabelFont(62);
    hUnfolded->GetXaxis()->SetLabelSize(0.035);

    hUnfolded->GetYaxis()->SetTitleFont(62);
    hUnfolded->GetYaxis()->SetTitleSize(0.04);
    hUnfolded->GetYaxis()->SetLabelFont(62);
    hUnfolded->GetYaxis()->SetLabelSize(0.035);
    hUnfolded->SetMinimum(0);  // y-axis starts at 0

    // --- Canvas ---
    TCanvas *c = new TCanvas("c", "Unfolded vs Generator-Level True", 800, 600);

    hUnfolded->Draw("PE");
    hGenTrue->Draw("PE SAME");

    // --- Legend ---
    TLegend *legend = new TLegend(0.6, 0.60, 0.88, 0.75);
    legend->AddEntry(hUnfolded, "Reconstructed SuperChic Unfolded", "lep");
    legend->AddEntry(hGenTrue, "SuperChic Generated Level", "lep");
    legend->Draw();

    // --- CMS + labels ---
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

    c->SaveAs("unfolded2018SuperChichrecoTrue_vs_hgenTrue.pdf");
    f->Close();
}

