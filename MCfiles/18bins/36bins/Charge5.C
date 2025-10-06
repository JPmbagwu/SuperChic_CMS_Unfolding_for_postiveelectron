#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>  // Include TLatex

void compare_measured_vs_recoTrue() {
    TFile *f = TFile::Open("FDC36enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) return;

    TH1D* hMeasured = (TH1D*)f->Get("measured;2");
    TH1F* hRecoTrue = (TH1F*)f->Get("hrecoTrue");

    if (!hMeasured || !hRecoTrue) {
        f->Close();
        return;
    }

    hMeasured->SetStats(0);
    hRecoTrue->SetStats(0);

    hMeasured->SetMarkerStyle(20);
    hMeasured->SetMarkerSize(1.0);
    hMeasured->SetMarkerColor(kBlue);
    hMeasured->SetLineColor(kBlue);
    hMeasured->SetLineWidth(2);

    hRecoTrue->SetMarkerStyle(24);
    hRecoTrue->SetMarkerSize(1.0);
    hRecoTrue->SetMarkerColor(kRed);
    hRecoTrue->SetLineColor(kRed);
    hRecoTrue->SetLineStyle(2);
    hRecoTrue->SetLineWidth(2);

    TCanvas *c = new TCanvas("c", "Measured vs RecoTrue", 800, 600);
    hMeasured->SetTitle(" ;#Delta#phi = #phi_{ee} - #phi_{e+};Entries");
    hMeasured->Draw("PE");
    hRecoTrue->Draw("PE SAME");

    TLegend *legend = new TLegend(0.65, 0.35, 0.88, 0.50);
    legend->AddEntry(hMeasured, "Forward Fold", "lep");
    legend->AddEntry(hRecoTrue, "2018 Superchic", "lep");
    legend->Draw();

    // Add TLatex CMS-style text
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

    latex.SetTextAlign(31);
    latex.SetTextSize(0.030);

    c->SaveAs("Forward_Fold_vs_2018Superchic.pdf");
    f->Close();
}
