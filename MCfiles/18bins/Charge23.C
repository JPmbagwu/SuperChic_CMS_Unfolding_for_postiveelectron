#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h> // Include TLatex

void compare_unfolded_vs_genTrue() {
    TFile *f = TFile::Open("FDC10enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) return;

    TH1D* hUnfolded = (TH1D*)f->Get("UnfoldSCMC2018;1");
    TH1F* hGenTrue = (TH1F*)f->Get("hgenTrue");

    if (!hUnfolded || !hGenTrue) {
        f->Close();
        return;
    }

    hUnfolded->SetStats(0);
    hGenTrue->SetStats(0);

    hUnfolded->SetMarkerStyle(20);
    hUnfolded->SetMarkerSize(1.0);
    hUnfolded->SetMarkerColor(kGreen+2);
    hUnfolded->SetLineColor(kGreen+2);
    hUnfolded->SetLineWidth(2);

    hGenTrue->SetMarkerStyle(25);
    hGenTrue->SetMarkerSize(1.0);
    hGenTrue->SetMarkerColor(kMagenta);
    hGenTrue->SetLineColor(kMagenta);
    hGenTrue->SetLineStyle(2);
    hGenTrue->SetLineWidth(2);

    TCanvas *c = new TCanvas("c", "Unfolded vs Generator-Level True", 800, 600);
    hUnfolded->SetTitle(" ;#Delta#phi = #phi_{ee} - #phi_{e+};Entries");
    hUnfolded->Draw("PE");
    hGenTrue->Draw("PE SAME");

    TLegend *legend = new TLegend(0.6, 0.45, 0.88, 0.60);
    legend->AddEntry(hUnfolded, "Reconstructed SuperChic Unfolded", "lep");
    legend->AddEntry(hGenTrue, "SuperChic Generated Level", "lep");
    legend->Draw();

    // Add TLatex annotations
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
    
    latex.SetTextAlign(31);
    latex.SetTextSize(0.030);

    c->SaveAs("unfolded2018SuperchichrecoTrue_vs_hgenTrue.pdf");
    f->Close();
}

