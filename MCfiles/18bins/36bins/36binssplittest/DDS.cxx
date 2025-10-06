#if !defined(__CINT__) && !defined(__CLING__) || defined(__ACLIC__)
#include <TSystem.h>
#include <TROOT.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldInvert.h"
#include <TH1D.h>
#include <TTree.h>
#include <TChain.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <iostream>
#include <fstream>
#include <TClonesArray.h>
#include <cmath>
#include <string>
#include "TRandom.h"
#include "TString.h"
#include "TLeaf.h"
#include "TLegend.h"
#include <TPad.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include <TGaxis.h>
#include <TVector2.h>
#include <TLatex.h>
#endif

using namespace std;

// Debug flag
const bool DEBUG = false;

// Centralized cut definitions
struct Cuts {
    static constexpr float MaxRapidity = 2.4;        // Individual positron |eta| < 2.4
    static constexpr float SingleElePtMin = 0.0;     // Single positron pT > 0.0 GeV/c
    static constexpr float maxDielePt = 1.0;         // Dielectron pT < 1.0 GeV/c
    static constexpr float minMass = 7.0;            // Minimum dielectron mass (GeV/c²)
    static constexpr float maxMass = 10.0;           // Maximum dielectron mass (GeV/c²)
    static constexpr float maxHFpCut = 7.3;          // Max HF+ energy (GeV)
    static constexpr float maxHFmCut = 7.6;          // Max HF- energy (GeV)
    static constexpr float maxVertexZ = 20.0;        // Max vertex z (cm)
    static constexpr float electronMass = 0.000511;  // Positron mass (GeV)
};

// Histogram settings
const float xrangeMax = TMath::Pi();
const float xrangeMin = -TMath::Pi();
const int xbinz = 36;

TChain* CreateChain(const string& inFile, const char* treeName) {
    TChain* chain = new TChain(treeName);
    if (inFile.find(".txt") != string::npos) {
        ifstream file(inFile);
        if (!file.is_open()) {
            cerr << "Error: Could not open file list " << inFile << endl;
            return chain;
        }
        string rootFile;
        while (getline(file, rootFile)) {
            if (!rootFile.empty()) {
                if (DEBUG) cout << "Adding: " << rootFile << endl;
                int added = chain->Add(rootFile.c_str());
                if (added == 0) {
                    cerr << "Error: Could not add file " << rootFile << endl;
                }
            }
        }
    } else {
        if (DEBUG) cout << "Adding single file: " << inFile << endl;
        int added = chain->Add(inFile.c_str());
        if (added == 0) {
            cerr << "Error: Could not add file " << inFile << endl;
        }
    }
    if (DEBUG) cout << "Chain has " << chain->GetEntries() << " entries" << endl;
    return chain;
}

float GetAngle(float vPX1, float vPY1, float vPX2, float vPY2) {
    TVector2 DielectronVect(vPX1, vPY1);
    TVector2 PositronVect(vPX2, vPY2); // MODIFIED FOR e+: Changed ElectronVect to PositronVect

    float DielectronMag = sqrt(DielectronVect.X() * DielectronVect.X() + DielectronVect.Y() * DielectronVect.Y());
    float PositronMag = sqrt(PositronVect.X() * PositronVect.X() + PositronVect.Y() * PositronVect.Y());

    if (DielectronMag == 0 || PositronMag == 0) return 0.0;

    float cosAngle = (DielectronVect.X() * PositronVect.X() + DielectronVect.Y() * PositronVect.Y()) / (DielectronMag * PositronMag);
    float sinAngle = (DielectronVect.X() * PositronVect.Y() - DielectronVect.Y() * PositronVect.X()) / (DielectronMag * PositronMag);

    float tanAngle = atan2(sinAngle, cosAngle);
    return tanAngle;
}

void ManuelPositronResponseMatrix() { // MODIFIED FOR e+: Changed function name to reflect positron analysis
    string inputMC = "/afs/cern.ch/user/j/jmbagwu/RooUnfold/examples/lblfiles/cms/my_file.txt";
    string outfile_pseudo = "DDS10_pseudodata_25_positron.root"; // MODIFIED FOR e+: Added positron to filename
    string outfile_response = "DDS10_response_75_positron.root"; // MODIFIED FOR e+: Added positron to filename
    
    // Check output directories
    if (gSystem->AccessPathName(gSystem->DirName(outfile_pseudo.c_str()))) {
        cerr << "Error: Output directory does not exist for pseudodata: " << gSystem->DirName(outfile_pseudo.c_str()) << endl;
        return;
    }
    if (gSystem->AccessPathName(gSystem->DirName(outfile_response.c_str()))) {
        cerr << "Error: Output directory does not exist for response: " << gSystem->DirName(outfile_response.c_str()) << endl;
        return;
    }
    
    TChain *t1 = CreateChain(inputMC, "ggHiNtuplizer/EventTree");
    if (t1->GetEntries() == 0) {
        cerr << "Error: No entries found in the chain!" << endl;
        return;
    }
    
    if (DEBUG) {
        cout << "Available branches in MC tree:" << endl;
        t1->GetListOfBranches()->Print();
    }
    
    // Create two output files
    TFile f_pseudo(outfile_pseudo.c_str(), "recreate");
    if (!f_pseudo.IsOpen()) {
        cerr << "Error: Could not create pseudodata output file " << outfile_pseudo << endl;
        return;
    }
    
    TFile f_response(outfile_response.c_str(), "recreate");
    if (!f_response.IsOpen()) {
        cerr << "Error: Could not create response output file " << outfile_response << endl;
        f_pseudo.Close();
        return;
    }
    
    // Create tree for response file only
    f_response.cd();
    TTree* SuperEventTree = new TTree("SuperEventTree", "Tree with Event Variables (75% subset)");
    
    // Branch definitions for SuperEventTree
    float maxHFp = 0;
    float maxHFm = 0;
    float maxTower_E = 0;
    float maxTower_eta = 0;
    float maxTower_phi = 0;
    SuperEventTree->Branch("maxHFp", &maxHFp);
    SuperEventTree->Branch("maxHFm", &maxHFm);
    SuperEventTree->Branch("maxTower_E", &maxTower_E);
    SuperEventTree->Branch("maxTower_eta", &maxTower_eta);
    SuperEventTree->Branch("maxTower_phi", &maxTower_phi);
    
    // Set branch addresses
    vector<float> *gen_elePt = nullptr;
    vector<float> *gen_elePhi = nullptr;
    vector<float> *gen_eleEta = nullptr;
    vector<int>   *gen_elePID = nullptr;
    vector<float> *gen_trkvx = nullptr;
    vector<float> *gen_trkvy = nullptr;
    vector<float> *gen_trkvz = nullptr;
    int gen_nEle = 0;
    
    t1->SetBranchAddress("mcPt",   &gen_elePt);
    t1->SetBranchAddress("mcPhi",  &gen_elePhi);
    t1->SetBranchAddress("mcEta",  &gen_eleEta);
    t1->SetBranchAddress("nMC",    &gen_nEle);
    t1->SetBranchAddress("mcPID",  &gen_elePID);
    t1->SetBranchAddress("mcVtx_x", &gen_trkvx);
    t1->SetBranchAddress("mcVtx_y", &gen_trkvy);
    t1->SetBranchAddress("mcVtx_z", &gen_trkvz);
    
    vector<float> *reco_elePt = nullptr;
    vector<float> *reco_elePhi = nullptr;
    vector<float> *reco_eleEta = nullptr;
    vector<int>   *reco_eleCharge = nullptr;
    vector<float> *reco_eleDz = nullptr;
    vector<float> *reco_trkvx = nullptr;
    vector<float> *reco_trkvy = nullptr;
    vector<float> *reco_trkvz = nullptr;
    int reco_nEle = 0;
    int reco_nTrk = 0;
    
    t1->SetBranchAddress("elePt",    &reco_elePt);
    t1->SetBranchAddress("elePhi",   &reco_elePhi);
    t1->SetBranchAddress("eleEta",   &reco_eleEta);
    t1->SetBranchAddress("nEle",     &reco_nEle);
    t1->SetBranchAddress("eleCharge",&reco_eleCharge);
    t1->SetBranchAddress("eleDz",    &reco_eleDz);
    t1->SetBranchAddress("nTrk",     &reco_nTrk);
    t1->SetBranchAddress("trkvx",    &reco_trkvx);
    t1->SetBranchAddress("trkvy",    &reco_trkvy);
    t1->SetBranchAddress("trkvz",    &reco_trkvz);
    
    vector<float> *CaloTower_e = nullptr;
    vector<float> *CaloTower_et = nullptr;
    vector<float> *CaloTower_eta = nullptr;
    vector<float> *CaloTower_hadE = nullptr;
    vector<float> *CaloTower_emE = nullptr;
    vector<float> *CaloTower_phi = nullptr;
    int nTower = 0;
    t1->SetBranchAddress("CaloTower_e",    &CaloTower_e);
    t1->SetBranchAddress("CaloTower_et",   &CaloTower_et);
    t1->SetBranchAddress("CaloTower_eta",  &CaloTower_eta);
    t1->SetBranchAddress("CaloTower_hadE", &CaloTower_hadE);
    t1->SetBranchAddress("CaloTower_emE",  &CaloTower_emE);
    t1->SetBranchAddress("CaloTower_phi",  &CaloTower_phi);
    t1->SetBranchAddress("nTower",         &nTower);
    
    // Initialize histograms
    RooUnfoldResponse response(xbinz, xrangeMin, xrangeMax);
    
    int count_True = 0;
    int count_Fake = 0;
    int count_Miss = 0;
    int Event_count = 0;
    int GEventPass = 0;
    int REventPass = 0;
    
    // Pseudodata histograms for 25% subset
    f_pseudo.cd();
    TH1F *pseudodata = new TH1F("pseudodata", "Reconstructed #Delta#phi (25% subset - PSEUDODATA)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *genTrue = new TH1F("genTrue", "True Generator Level #Delta#phi (25% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *gen = new TH1F("gen", "Generator Level #Delta#phi (25% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    
    // Histograms for 75% subset
    f_response.cd();
    TH1F *hGenMass = new TH1F("hGenMass", "Generator Level Dielectron Mass", 200, 0, 100);
    TH1F *hGenMassCuts = new TH1F("hGenMassCuts", "Generator Level Dielectron Mass After Cuts", 200, 7.0, 100);
    TH1F *hRecoMass = new TH1F("hRecoMass", "Reconstructed Dielectron Mass", 200, 0, 100);
    TH1F *hRecoMassCuts = new TH1F("hRecoMassCuts", "Reconstructed Dielectron Mass After Cuts", 200, 7.0, 100);
    TH1F *hreco = new TH1F("hreco", "Reconstructed #Delta#phi (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hrecoTrue = new TH1F("hrecoTrue", "True Reconstructed #Delta#phi (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hrecoFake = new TH1F("hrecoFake", "Fake Reconstructed #Delta#phi (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hgen = new TH1F("hgen", "Generator Level #Delta#phi (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hgenTrue = new TH1F("hgenTrue", "True Generator Level #Delta#phi (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hgenMiss = new TH1F("hgenMiss", "Missed Generator Level #Delta#phi (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hgenAll = new TH1F("hgenAll", "All Generator Level #Delta#phi (no cuts)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hRDeltaPhi = new TH1F("hRDeltaPhi", "Reconstructed DeltaPhi Explicit (25% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hGDeltaPhi = new TH1F("hGDeltaPhi", "Generator DeltaPhi Explicit (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hEfficiency = new TH1F("hEfficiency", "Efficiency;#Delta#phi;Efficiency (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hPurity = new TH1F("hPurity", "Purity;#Delta#phi;Purity (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hAcceptance = new TH1F("hAcceptance", "Acceptance;#Delta#phi;Acceptance (75% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hUnfoldedPseudodataPurityCorrected = new TH1F("hUnfoldedPseudodataPurityCorrected", "Unfolded Purity-Corrected Pseudodata #Delta#phi (25% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hUnfoldedPseudodataPurityCorrectedOverAcceptance = new TH1F("hUnfoldedPseudodataPurityCorrectedOverAcceptance", "Unfolded Purity-Corrected Pseudodata #Delta#phi / Acceptance (25% subset, 75% acceptance)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hUnfoldedRecoTrue = new TH1F("hUnfoldedRecoTrue", "Unfolded True Reconstructed #Delta#phi (25% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hUnfoldedRecoTrueOverAcceptance = new TH1F("hUnfoldedRecoTrueOverAcceptance", "Unfolded True Reconstructed #Delta#phi / Acceptance (25% subset, 75% acceptance)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    TH1F *hRefolded = new TH1F("hRefolded", "Refolded Purity-Corrected Pseudodata #Delta#phi (25% subset)", xbinz, xrangeMin, xrangeMax); // MODIFIED FOR e+: Updated title
    
    float RDeltaPhi = 0, GDeltaPhi = 0;
    int numberentry = t1->GetEntries();
    
    cout << "Total events: " << numberentry << endl;
    cout << "25% subset (pseudodata): " << 0.25 * numberentry << " events" << endl;
    cout << "75% subset (response): " << 0.75 * numberentry << " events" << endl;
    
    // ======================================================================
    // FIRST 25% OF EVENTS: PSEUDODATA (reco cuts + genTrue + gen)
    // ======================================================================
    int pseudodata_count = 0;
    int pseudodata_true_count = 0;
    int pseudodata_gen_count = 0;
    for (int j = 0; j < 0.25 * numberentry; j++) {
        t1->GetEntry(j);
        
        int passRecoCuts = 0;
        int passGenCuts = 0;
        int TrueEvt = 0;
        
        int iRp = -1, iRn = -1;
        int EleP = -1, EleN = -1;
        
        TLorentzVector RLep1, RLep2, RDiele;
        TLorentzVector GLep1, GLep2, GDiele;
        
        // Generator-level processing
        int idx_electron = -1;
        int idx_positron = -1;
        for (size_t k = 0; k < gen_elePID->size(); ++k) {
            int pid = gen_elePID->at(k);
            if (abs(pid) >= 20) continue;
            
            if (pid == 11 && idx_electron == -1) {
                idx_electron = k;
            }
            else if (pid == -11 && idx_positron == -1) {
                idx_positron = k;
            }
        }
        
        if (idx_electron != -1 && idx_positron != -1) {
            GLep1.SetPtEtaPhiM(gen_elePt->at(idx_electron), gen_eleEta->at(idx_electron),
                               gen_elePhi->at(idx_electron), Cuts::electronMass);
            GLep2.SetPtEtaPhiM(gen_elePt->at(idx_positron), gen_eleEta->at(idx_positron),
                               gen_elePhi->at(idx_positron), Cuts::electronMass);
            GDiele = GLep1 + GLep2;
            
            bool genTwoelectron = (gen_nEle >= 2);
            bool genElectroncharge = (gen_elePID->at(idx_positron) + gen_elePID->at(idx_electron) == 0);
            bool GRapidity = abs(GLep1.Eta()) <= Cuts::MaxRapidity && abs(GLep2.Eta()) <= Cuts::MaxRapidity;
            bool GMassCut = (GDiele.M() > Cuts::minMass && GDiele.M() < Cuts::maxMass);
            bool GpTcut = (GDiele.Pt() < Cuts::maxDielePt);
            bool GVertexCut = (fabs(gen_trkvz->at(idx_positron)) <= Cuts::maxVertexZ &&
                               fabs(gen_trkvz->at(idx_electron)) <= Cuts::maxVertexZ);
            bool GenSingleElectronPtCut = GLep1.Pt() > Cuts::SingleElePtMin && GLep2.Pt() > Cuts::SingleElePtMin;
            
            float GenDielePx = GDiele.Pt() * cos(GDiele.Phi());
            float GenDielePy = GDiele.Pt() * sin(GDiele.Phi());
            float GenPosiPx = GLep2.Pt() * cos(GLep2.Phi()); // MODIFIED FOR e+: Use positron (idx_positron) for angle
            float GenPosiPy = GLep2.Pt() * sin(GLep2.Phi()); // MODIFIED FOR e+: Use positron (idx_positron) for angle
            GDeltaPhi = GetAngle(GenDielePx, GenDielePy, GenPosiPx, GenPosiPy); // MODIFIED FOR e+: Updated for positron
            
            if (GMassCut && GRapidity && GpTcut && genTwoelectron && genElectroncharge && GVertexCut && GenSingleElectronPtCut) {
                passGenCuts = 1;
                GEventPass++;
                f_pseudo.cd();
                gen->Fill(GDeltaPhi);
                pseudodata_gen_count++;
                f_response.cd();
                if (DEBUG && pseudodata_gen_count % 1000 == 0) {
                    cout << "Pseudodata gen event " << pseudodata_gen_count << " (25% subset): GDeltaPhi = " << GDeltaPhi << endl;
                }
            }
        }
        
        // Reconstruction-level processing
        bool Twoelectron = (reco_nEle == 2);
        bool Twotrk = (reco_nTrk == 2);
        
        if (Twoelectron && Twotrk) {
            if (reco_eleCharge->size() < 2) continue;
            bool Electroncharge = (reco_eleCharge->at(0) * reco_eleCharge->at(1) == -1);
            if (!Electroncharge) continue;
            
            if (reco_trkvz->size() == 0) continue;
            
            EleP = ((*reco_eleCharge)[0] == 1) ? 0 : 1; // Positron has charge +1
            EleN = 1 - EleP;
            iRp = EleP;
            iRn = EleN;
            
            RLep1.SetPtEtaPhiM((*reco_elePt)[EleP], (*reco_eleEta)[EleP], (*reco_elePhi)[EleP], Cuts::electronMass); // Positron
            RLep2.SetPtEtaPhiM((*reco_elePt)[EleN], (*reco_eleEta)[EleN], (*reco_elePhi)[EleN], Cuts::electronMass); // Electron
            RDiele = RLep1 + RLep2;
            
            float RecoDielePx = RDiele.Pt() * cos(RDiele.Phi());
            float RecoDielePy = RDiele.Pt() * sin(RDiele.Phi());
            float RecoPosiPx = RLep1.Pt() * cos(RLep1.Phi()); // MODIFIED FOR e+: Use positron (EleP) for angle
            float RecoPosiPy = RLep1.Pt() * sin(RLep1.Phi()); // MODIFIED FOR e+: Use positron (EleP) for angle
            RDeltaPhi = GetAngle(RecoDielePx, RecoDielePy, RecoPosiPx, RecoPosiPy); // MODIFIED FOR e+: Updated for positron
            
            f_response.cd();
            hRecoMass->Fill(RDiele.M());
            hRDeltaPhi->Fill(RDeltaPhi);
            f_pseudo.cd();
            
            maxHFp = 0;
            maxHFm = 0;
            maxTower_E = 0;
            maxTower_eta = 0;
            maxTower_phi = 0;
            for (int t = 0; t < nTower; t++) {
                float tower_eta = CaloTower_eta->at(t);
                float tower_E = CaloTower_e->at(t);
                if (tower_eta > 2.9 && tower_eta < 5.2)
                    maxHFp = TMath::Max(maxHFp, tower_E);
                if (tower_eta < -2.9 && tower_eta > -5.2)
                    maxHFm = TMath::Max(maxHFm, tower_E);
                if (tower_E > maxTower_E) {
                    maxTower_E = tower_E;
                    maxTower_eta = tower_eta;
                    maxTower_phi = CaloTower_phi->at(t);
                }
            }
            
            bool RRapidity = abs(RLep1.Eta()) <= Cuts::MaxRapidity && abs(RLep2.Eta()) <= Cuts::MaxRapidity;
            bool RVertexCut = fabs((*reco_trkvz)[0]) <= Cuts::maxVertexZ;
            bool RMassCut = (RDiele.M() > Cuts::minMass && RDiele.M() < Cuts::maxMass);
            bool RpTCut = (RDiele.Pt() < Cuts::maxDielePt);
            bool HFCuts = (maxHFm <= Cuts::maxHFmCut && maxHFp <= Cuts::maxHFpCut);
            bool RecoSingleElectronPtCut = RLep1.Pt() > Cuts::SingleElePtMin && RLep2.Pt() > Cuts::SingleElePtMin;
            
            if (RMassCut && HFCuts && RVertexCut && RRapidity && RecoSingleElectronPtCut && RpTCut) {
                f_response.cd();
                hRecoMassCuts->Fill(RDiele.M());
                f_pseudo.cd();
                pseudodata->Fill(RDeltaPhi);
                passRecoCuts = 1;
                pseudodata_count++;
                if (DEBUG && pseudodata_count % 1000 == 0) {
                    cout << "Pseudodata event " << pseudodata_count << " (25% subset): RDeltaPhi = " << RDeltaPhi << endl;
                }
            }
        }
        
        // Fill genTrue for true events in 25% subset
        if (passGenCuts == 1 && passRecoCuts == 1) {
            TrueEvt = 1;
            pseudodata_true_count++;
            f_pseudo.cd();
            genTrue->Fill(GDeltaPhi);
            f_response.cd();
            if (DEBUG && pseudodata_true_count % 1000 == 0) {
                cout << "Pseudodata true event " << pseudodata_true_count << " (25% subset): GDeltaPhi = " << GDeltaPhi << endl;
            }
        }
        
        Event_count++;
    }
    
    // Write pseudodata histograms
    f_pseudo.cd();
    pseudodata->Write();
    genTrue->Write();
    gen->Write();
    f_pseudo.Close();
    
    // ======================================================================
    // REMAINING 75% OF EVENTS: RESPONSE MATRIX BUILDING
    // ======================================================================
    int response_true_count = 0;
    int response_fake_count = 0;
    int response_miss_count = 0;
    
    f_response.cd();
    for (int j = 0.25 * numberentry; j < numberentry; j++) {
        t1->GetEntry(j);
        
        int FakeEvt = 0;
        int TrueEvt = 0;
        int MissEvt = 0;
        int passGenCuts = 0;
        int passRecoCuts = 0;
        
        int iRp = -1, iRn = -1;
        int EleP = -1, EleN = -1;
        
        TLorentzVector RLep1, RLep2, RDiele;
        
        int idx_electron = -1;
        int idx_positron = -1;
        
        for (size_t k = 0; k < gen_elePID->size(); ++k) {
            int pid = gen_elePID->at(k);
            if (abs(pid) >= 20) continue;
            
            if (pid == 11 && idx_electron == -1) {
                idx_electron = k;
            }
            else if (pid == -11 && idx_positron == -1) {
                idx_positron = k;
            }
        }
        
        if (idx_electron == -1 || idx_positron == -1) continue;
        
        TLorentzVector GLep1, GLep2, GDiele;
        GLep1.SetPtEtaPhiM(gen_elePt->at(idx_electron), gen_eleEta->at(idx_electron),
                           gen_elePhi->at(idx_electron), Cuts::electronMass);
        GLep2.SetPtEtaPhiM(gen_elePt->at(idx_positron), gen_eleEta->at(idx_positron),
                           gen_elePhi->at(idx_positron), Cuts::electronMass);
        GDiele = GLep1 + GLep2;
        
        bool genTwoelectron = (gen_nEle >= 2);
        bool genElectroncharge = (gen_elePID->at(idx_positron) + gen_elePID->at(idx_electron) == 0);
        bool GRapidity = abs(GLep1.Eta()) <= Cuts::MaxRapidity && abs(GLep2.Eta()) <= Cuts::MaxRapidity;
        bool GMassCut = (GDiele.M() > Cuts::minMass && GDiele.M() < Cuts::maxMass);
        bool GpTcut = (GDiele.Pt() < Cuts::maxDielePt);
        bool GVertexCut = (fabs(gen_trkvz->at(idx_positron)) <= Cuts::maxVertexZ &&
                           fabs(gen_trkvz->at(idx_electron)) <= Cuts::maxVertexZ);
        bool GenSingleElectronPtCut = GLep1.Pt() > Cuts::SingleElePtMin && GLep2.Pt() > Cuts::SingleElePtMin;
        
        float GenDielePx = GDiele.Pt() * cos(GDiele.Phi());
        float GenDielePy = GDiele.Pt() * sin(GDiele.Phi());
        float GenPosiPx = GLep2.Pt() * cos(GLep2.Phi()); // MODIFIED FOR e+: Use positron (idx_positron) for angle
        float GenPosiPy = GLep2.Pt() * sin(GLep2.Phi()); // MODIFIED FOR e+: Use positron (idx_positron) for angle
        GDeltaPhi = GetAngle(GenDielePx, GenDielePy, GenPosiPx, GenPosiPy); // MODIFIED FOR e+: Updated for positron
        
        hgenAll->Fill(GDeltaPhi);
        hGenMass->Fill(GDiele.M());
        
        if (GMassCut && GRapidity && GpTcut && genTwoelectron && genElectroncharge && GVertexCut && GenSingleElectronPtCut) {
            passGenCuts = 1;
            GEventPass++;
            hgen->Fill(GDeltaPhi);
            hGDeltaPhi->Fill(GDeltaPhi);
            hGenMassCuts->Fill(GDiele.M());
        }
        
        bool Twoelectron = (reco_nEle == 2);
        bool Twotrk = (reco_nTrk == 2);
        
        if (Twoelectron && Twotrk) {
            if (reco_eleCharge->size() < 2) continue;
            bool Electroncharge = (reco_eleCharge->at(0) * reco_eleCharge->at(1) == -1);
            if (!Electroncharge) continue;
            
            if (reco_trkvz->size() == 0) continue;
            
            EleP = ((*reco_eleCharge)[0] == 1) ? 0 : 1; // Positron has charge +1
            EleN = 1 - EleP;
            iRp = EleP;
            iRn = EleN;
            
            RLep1.SetPtEtaPhiM((*reco_elePt)[EleP], (*reco_eleEta)[EleP], (*reco_elePhi)[EleP], Cuts::electronMass); // Positron
            RLep2.SetPtEtaPhiM((*reco_elePt)[EleN], (*reco_eleEta)[EleN], (*reco_elePhi)[EleN], Cuts::electronMass); // Electron
            RDiele = RLep1 + RLep2;
            
            float RecoDielePx = RDiele.Pt() * cos(RDiele.Phi());
            float RecoDielePy = RDiele.Pt() * sin(RDiele.Phi());
            float RecoPosiPx = RLep1.Pt() * cos(RLep1.Phi()); // MODIFIED FOR e+: Use positron (EleP) for angle
            float RecoPosiPy = RLep1.Pt() * sin(RLep1.Phi()); // MODIFIED FOR e+: Use positron (EleP) for angle
            RDeltaPhi = GetAngle(RecoDielePx, RecoDielePy, RecoPosiPx, RecoPosiPy); // MODIFIED FOR e+: Updated for positron
            
            hRecoMass->Fill(RDiele.M());
            
            maxHFp = 0;
            maxHFm = 0;
            maxTower_E = 0;
            maxTower_eta = 0;
            maxTower_phi = 0;
            for (int t = 0; t < nTower; t++) {
                float tower_eta = CaloTower_eta->at(t);
                float tower_E = CaloTower_e->at(t);
                if (tower_eta > 2.9 && tower_eta < 5.2)
                    maxHFp = TMath::Max(maxHFp, tower_E);
                if (tower_eta < -2.9 && tower_eta > -5.2)
                    maxHFm = TMath::Max(maxHFm, tower_E);
                if (tower_E > maxTower_E) {
                    maxTower_E = tower_E;
                    maxTower_eta = tower_eta;
                    maxTower_phi = CaloTower_phi->at(t);
                }
            }
            
            bool RRapidity = abs(RLep1.Eta()) <= Cuts::MaxRapidity && abs(RLep2.Eta()) <= Cuts::MaxRapidity;
            bool RVertexCut = fabs((*reco_trkvz)[0]) <= Cuts::maxVertexZ;
            bool RMassCut = (RDiele.M() > Cuts::minMass && RDiele.M() < Cuts::maxMass);
            bool RpTCut = (RDiele.Pt() < Cuts::maxDielePt);
            bool HFCuts = (maxHFm <= Cuts::maxHFmCut && maxHFp <= Cuts::maxHFpCut);
            bool RecoSingleElectronPtCut = RLep1.Pt() > Cuts::SingleElePtMin && RLep2.Pt() > Cuts::SingleElePtMin;
            
            if (RMassCut && HFCuts && RVertexCut && RRapidity && RecoSingleElectronPtCut && RpTCut) {
                hreco->Fill(RDeltaPhi);
                hRecoMassCuts->Fill(RDiele.M());
                REventPass++;
                passRecoCuts = 1;
            }
        }
        
        if (passGenCuts == 1 && passRecoCuts == 0) {
            MissEvt = 1;
            response_miss_count++;
            hgenMiss->Fill(GDeltaPhi);
        }
        if (passGenCuts == 0 && passRecoCuts == 1) {
            FakeEvt = 1;
            response_fake_count++;
            hrecoFake->Fill(RDeltaPhi);
        }
        if (passGenCuts == 1 && passRecoCuts == 1) {
            TrueEvt = 1;
            response_true_count++;
            response.Fill(RDeltaPhi, GDeltaPhi);
            hgenTrue->Fill(GDeltaPhi);
            hrecoTrue->Fill(RDeltaPhi);
            if (DEBUG && response_true_count % 1000 == 0) {
                cout << "Response matrix event " << response_true_count << " (75% subset): RDeltaPhi = " << RDeltaPhi << ", GDeltaPhi = " << GDeltaPhi << endl;
            }
        }
        
        if (MissEvt == 1) {
            count_Miss++;
        }
        if (FakeEvt == 1) {
            count_Fake++;
        }
        if (TrueEvt == 1) {
            count_True++;
        }
        
        Event_count++;
        SuperEventTree->Fill();
    }
    
    // Apply purity correction to pseudodata
    TFile *f_pseudo_read = new TFile(outfile_pseudo.c_str(), "READ");
    if (!f_pseudo_read || f_pseudo_read->IsZombie()) {
        cerr << "Error: Could not open pseudodata file " << outfile_pseudo << endl;
        f_response.Close();
        return;
    }
    
    TH1F *pseudodata_read = (TH1F*)f_pseudo_read->Get("pseudodata");
    if (!pseudodata_read) {
        cerr << "Error: Could not read pseudodata from " << outfile_pseudo << endl;
        f_pseudo_read->Close();
        f_response.Close();
        return;
    }
    
    // Calculate acceptance using 75% subset
    for (int iBinX = 1; iBinX <= hAcceptance->GetNbinsX(); ++iBinX) {
        double genTrue = hgenTrue->GetBinContent(iBinX);
        double genMiss = hgenMiss->GetBinContent(iBinX);
        double totalGenPassing = genTrue + genMiss;
        
        if (totalGenPassing > 0) {
            double acceptance = genTrue / totalGenPassing;
            double error = sqrt(acceptance * (1 - acceptance) / totalGenPassing);
            hAcceptance->SetBinContent(iBinX, acceptance);
            hAcceptance->SetBinError(iBinX, error);
        } else {
            hAcceptance->SetBinContent(iBinX, 0);
            hAcceptance->SetBinError(iBinX, 0);
        }
    }
    
    // Calculate purity using hrecoTrue / (hrecoTrue + hrecoFake)
    f_response.cd();
    for (int iBinX = 1; iBinX <= hPurity->GetNbinsX(); ++iBinX) {
        double trueReco = hrecoTrue->GetBinContent(iBinX);
        double fakeReco = hrecoFake->GetBinContent(iBinX);
        
        double denomPur = trueReco + fakeReco;
        double purity = (denomPur > 0) ? trueReco / denomPur : 0.0;
        double purError = 0.0;
        if (denomPur > 0) {
            purError = sqrt(purity * (1 - purity) / denomPur);
        }
        hPurity->SetBinContent(iBinX, purity);
        hPurity->SetBinError(iBinX, purError);
    }
    
    // Create pseudodataPurityCorrected
    TH1F *pseudodataPurityCorrected = (TH1F*)pseudodata_read->Clone("pseudodataPurityCorrected");
    pseudodataPurityCorrected->SetTitle("Purity-Corrected Pseudodata #Delta#phi (25% subset)"); // MODIFIED FOR e+: Updated title
    for (int iBinX = 1; iBinX <= pseudodataPurityCorrected->GetNbinsX(); ++iBinX) {
        double recoVal = pseudodata_read->GetBinContent(iBinX);
        double purity = hPurity->GetBinContent(iBinX);
        double recoErr = pseudodata_read->GetBinError(iBinX);
        double purityErr = hPurity->GetBinError(iBinX);
        double correctedVal = recoVal * purity;
        double relErr = 0.0;
        if (recoVal > 0 && purity > 0) {
            relErr = sqrt(pow(recoErr/recoVal, 2) + pow(purityErr/purity, 2));
        }
        pseudodataPurityCorrected->SetBinContent(iBinX, correctedVal);
        pseudodataPurityCorrected->SetBinError(iBinX, correctedVal * relErr);
    }
    
    // Unfold pseudodataPurityCorrected
    cout << "Performing Split MC Test..." << endl;
    cout << "Pseudodata (pseudodata) entries: " << pseudodata_read->GetEntries() << endl;
    cout << "Response matrix truth entries: " << hgenTrue->GetEntries() << endl;
    
    RooUnfoldBayes unfold_pseudodata(&response, pseudodataPurityCorrected, 2);
    hUnfoldedPseudodataPurityCorrected = (TH1F*)unfold_pseudodata.Hunfold();
    hUnfoldedPseudodataPurityCorrected->SetName("hUnfoldedPseudodataPurityCorrected");
    hUnfoldedPseudodataPurityCorrected->SetTitle("Unfolded Purity-Corrected Pseudodata #Delta#phi (25% subset)");
    hUnfoldedPseudodataPurityCorrected->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    
    // Refold the unfolded histogram
    hRefolded = (TH1F*)response.ApplyToTruth(hUnfoldedPseudodataPurityCorrected, "hRefolded");
    hRefolded->SetTitle("Refolded Purity-Corrected Pseudodata #Delta#phi (25% subset)");
    hRefolded->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    hRefolded->GetYaxis()->SetTitle("Events");
    
    // Create hUnfoldedPseudodataPurityCorrectedOverAcceptance
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetTitle("Unfolded Purity-Corrected Pseudodata #Delta#phi / Acceptance (25% subset, 75% acceptance)");
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetYaxis()->SetTitle("Unfolded Pseudodata / Acceptance");
    for (int iBinX = 1; iBinX <= hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetNbinsX(); ++iBinX) {
        double unfoldVal = hUnfoldedPseudodataPurityCorrected->GetBinContent(iBinX);
        double acceptance = hAcceptance->GetBinContent(iBinX);
        double unfoldErr = hUnfoldedPseudodataPurityCorrected->GetBinError(iBinX);
        double acceptanceErr = hAcceptance->GetBinError(iBinX);
        if (acceptance > 0 && unfoldVal > 0) {
            double ratio = unfoldVal / acceptance;
            double relErr = sqrt(pow(unfoldErr/unfoldVal, 2) + pow(acceptanceErr/acceptance, 2));
            hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetBinContent(iBinX, ratio);
            hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetBinError(iBinX, ratio * relErr);
        } else {
            hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetBinContent(iBinX, 0);
            hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetBinError(iBinX, 0);
        }
    }
    
    // Unfold hrecoTrue (for true events)
    TH1F *hrecoTrue_read = (TH1F*)f_pseudo_read->Get("pseudodata"); // Using pseudodata as hrecoTrue for consistency
    if (!hrecoTrue_read) {
        cerr << "Error: Could not read pseudodata from " << outfile_pseudo << endl;
        f_pseudo_read->Close();
        f_response.Close();
        return;
    }
    
    RooUnfoldBayes unfold_recoTrue(&response, hrecoTrue_read, 2);
    hUnfoldedRecoTrue = (TH1F*)unfold_recoTrue.Hunfold();
    hUnfoldedRecoTrue->SetName("hUnfoldedRecoTrue");
    hUnfoldedRecoTrue->SetTitle("Unfolded True Reconstructed #Delta#phi (25% subset)");
    hUnfoldedRecoTrue->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    
    // Create hUnfoldedRecoTrueOverAcceptance
    hUnfoldedRecoTrueOverAcceptance->SetTitle("Unfolded True Reconstructed #Delta#phi / Acceptance (25% subset, 75% acceptance)");
    hUnfoldedRecoTrueOverAcceptance->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    hUnfoldedRecoTrueOverAcceptance->GetYaxis()->SetTitle("Unfolded RecoTrue / Acceptance");
    for (int iBinX = 1; iBinX <= hUnfoldedRecoTrueOverAcceptance->GetNbinsX(); ++iBinX) {
        double unfoldVal = hUnfoldedRecoTrue->GetBinContent(iBinX);
        double acceptance = hAcceptance->GetBinContent(iBinX);
        double unfoldErr = hUnfoldedRecoTrue->GetBinError(iBinX);
        double acceptanceErr = hAcceptance->GetBinError(iBinX);
        if (acceptance > 0 && unfoldVal > 0) {
            double ratio = unfoldVal / acceptance;
            double relErr = sqrt(pow(unfoldErr/unfoldVal, 2) + pow(acceptanceErr/acceptance, 2));
            hUnfoldedRecoTrueOverAcceptance->SetBinContent(iBinX, ratio);
            hUnfoldedRecoTrueOverAcceptance->SetBinError(iBinX, ratio * relErr);
        } else {
            hUnfoldedRecoTrueOverAcceptance->SetBinContent(iBinX, 0);
            hUnfoldedRecoTrueOverAcceptance->SetBinError(iBinX, 0);
        }
    }
    
    // Validation: Compare unfolded pseudodata vs generator truth (75% subset)
    TCanvas *cSplitMCValidation = new TCanvas("cSplitMCValidation", "Split MC Test: Unfolded Purity-Corrected Pseudodata vs Generator Truth (75% subset)", 900, 700);
    hUnfoldedPseudodataPurityCorrected->SetLineColor(kBlue);
    hUnfoldedPseudodataPurityCorrected->SetMarkerStyle(20);
    hUnfoldedPseudodataPurityCorrected->SetLineWidth(2);
    hUnfoldedPseudodataPurityCorrected->SetTitle("Split MC Test: Unfolded Purity-Corrected Pseudodata vs Generator Truth (75% subset)");
    hUnfoldedPseudodataPurityCorrected->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    hUnfoldedPseudodataPurityCorrected->GetYaxis()->SetTitle("Events");
    
    hgenTrue->SetLineColor(kRed);
    hgenTrue->SetMarkerStyle(21);
    hgenTrue->SetLineWidth(2);
    
    double max_unfolded = hUnfoldedPseudodataPurityCorrected->GetMaximum();
    double max_genTrue = hgenTrue->GetMaximum();
    double max_y = TMath::Max(max_unfolded, max_genTrue) * 1.2;
    hUnfoldedPseudodataPurityCorrected->SetMaximum(max_y);
    
    hUnfoldedPseudodataPurityCorrected->Draw("E");
    hgenTrue->Draw("E SAME");
    
    TLegend *legValidation = new TLegend(0.65, 0.7, 0.9, 0.85);
    legValidation->AddEntry(hUnfoldedPseudodataPurityCorrected, "Unfolded Purity-Corrected Pseudodata (25% subset)", "lep");
    legValidation->AddEntry(hgenTrue, "Generator Truth (75% subset)", "lep");
    legValidation->Draw();
    
    double chi2_validation = 0;
    int ndf_validation = 0;
    for (int i = 1; i <= hUnfoldedPseudodataPurityCorrected->GetNbinsX(); i++) {
        double diff = hUnfoldedPseudodataPurityCorrected->GetBinContent(i) - hgenTrue->GetBinContent(i);
        double err1 = hUnfoldedPseudodataPurityCorrected->GetBinError(i);
        double err2 = hgenTrue->GetBinError(i);
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2_validation += (diff*diff)/(err*err);
            ndf_validation++;
        }
    }
    
    TPaveText* pt_validation = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt_validation->SetFillColor(0);
    pt_validation->SetBorderSize(1);
    pt_validation->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2_validation, ndf_validation));
    pt_validation->AddText(Form("P-value = %.3f", TMath::Prob(chi2_validation, ndf_validation)));
    pt_validation->Draw();
    
    cSplitMCValidation->Update();
    cSplitMCValidation->SaveAs("DDS10SplitMCValidation_split_25_75_positron.pdf"); // MODIFIED FOR e+: Updated filename
    
    // Refold validation: Compare refolded vs pseudodataPurityCorrected
    TCanvas *cRefoldValidation = new TCanvas("cRefoldValidation", "Refold Validation: Refolded vs Purity-Corrected Pseudodata (25% subset)", 900, 700);
    hRefolded->SetLineColor(kBlue);
    hRefolded->SetMarkerStyle(20);
    hRefolded->SetLineWidth(2);
    hRefolded->SetTitle("Refold Validation: Refolded vs Purity-Corrected Pseudodata (25% subset)");
    hRefolded->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    hRefolded->GetYaxis()->SetTitle("Events");
    
    pseudodataPurityCorrected->SetLineColor(kRed);
    pseudodataPurityCorrected->SetMarkerStyle(21);
    pseudodataPurityCorrected->SetLineWidth(2);
    
    double max_refolded = hRefolded->GetMaximum();
    double max_pseudodata = pseudodataPurityCorrected->GetMaximum();
    max_y = TMath::Max(max_refolded, max_pseudodata) * 1.2;
    hRefolded->SetMaximum(max_y);
    
    hRefolded->Draw("E");
    pseudodataPurityCorrected->Draw("E SAME");
    
    TLegend *legRefold = new TLegend(0.65, 0.7, 0.9, 0.85);
    legRefold->AddEntry(hRefolded, "Refolded Purity-Corrected Pseudodata (25% subset)", "lep");
    legRefold->AddEntry(pseudodataPurityCorrected, "Purity-Corrected Pseudodata (25% subset)", "lep");
    legRefold->Draw();
    
    double chi2_refold = 0;
    int ndf_refold = 0;
    for (int i = 1; i <= hRefolded->GetNbinsX(); i++) {
        double diff = hRefolded->GetBinContent(i) - pseudodataPurityCorrected->GetBinContent(i);
        double err1 = hRefolded->GetBinError(i);
        double err2 = pseudodataPurityCorrected->GetBinError(i);
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2_refold += (diff*diff)/(err*err);
            ndf_refold++;
        }
    }
    
    TPaveText* pt_refold = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt_refold->SetFillColor(0);
    pt_refold->SetBorderSize(1);
    pt_refold->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2_refold, ndf_refold));
    pt_refold->AddText(Form("P-value = %.3f", TMath::Prob(chi2_refold, ndf_refold)));
    pt_refold->Draw();
    
    cRefoldValidation->Update();
    cRefoldValidation->SaveAs("DDS10RefoldValidation_Superchic_positron.pdf"); // MODIFIED FOR e+: Updated filename
    
    // Compare hUnfoldedPseudodataPurityCorrectedOverAcceptance vs gen (25% subset, unscaled)
    TH1F *gen25_read = (TH1F*)f_pseudo_read->Get("gen");
    if (!gen25_read) {
        cerr << "Error: Could not read gen histogram from " << outfile_pseudo << endl;
        f_pseudo_read->Close();
        f_response.Close();
        return;
    }
    
    TCanvas *cUnfoldedPseudodataPurityCorrectedOverAcceptanceVsGen25 = new TCanvas("cUnfoldedPseudodataPurityCorrectedOverAcceptanceVsGen25", "Split MC Test: Unfolded Purity-Corrected/Acceptance vs Generator Truth (25% subset, Unscaled)", 900, 700);
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetLineColor(kBlue);
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetMarkerStyle(20);
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetLineWidth(2);
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetTitle("Split MC Test: Unfolded Purity-Corrected/Acceptance vs Generator Truth (25% subset, Unscaled)");
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetYaxis()->SetTitle("Events");
    
    gen25_read->SetLineColor(kRed);
    gen25_read->SetMarkerStyle(21);
    gen25_read->SetLineWidth(2);
    
    double max_unfolded_purity = hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetMaximum();
    double max_gen25 = gen25_read->GetMaximum();
    max_y = TMath::Max(max_unfolded_purity, max_gen25) * 1.2;
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->SetMaximum(max_y);
    
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->Draw("E");
    gen25_read->Draw("E SAME");
    
    TLegend *legUnfoldedPurityVsGen25 = new TLegend(0.65, 0.7, 0.9, 0.85);
    legUnfoldedPurityVsGen25->AddEntry(hUnfoldedPseudodataPurityCorrectedOverAcceptance, "Unfolded Purity-Corrected Pseudodata / Acceptance (25% subset)", "lep");
    legUnfoldedPurityVsGen25->AddEntry(gen25_read, "Generator Truth (25% subset, All Gen Events)", "lep");
    legUnfoldedPurityVsGen25->Draw();
    
    double chi2_unfolded_purity_vs_gen25 = 0;
    int ndf_unfolded_purity_vs_gen25 = 0;
    for (int i = 1; i <= hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetNbinsX(); i++) {
        double diff = hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetBinContent(i) - gen25_read->GetBinContent(i);
        double err1 = hUnfoldedPseudodataPurityCorrectedOverAcceptance->GetBinError(i);
        double err2 = gen25_read->GetBinError(i);
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2_unfolded_purity_vs_gen25 += (diff*diff)/(err*err);
            ndf_unfolded_purity_vs_gen25++;
        }
    }
    
    TPaveText* pt_unfolded_purity_vs_gen25 = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt_unfolded_purity_vs_gen25->SetFillColor(0);
    pt_unfolded_purity_vs_gen25->SetBorderSize(1);
    pt_unfolded_purity_vs_gen25->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2_unfolded_purity_vs_gen25, ndf_unfolded_purity_vs_gen25));
    pt_unfolded_purity_vs_gen25->AddText(Form("P-value = %.3f", TMath::Prob(chi2_unfolded_purity_vs_gen25, ndf_unfolded_purity_vs_gen25)));
    pt_unfolded_purity_vs_gen25->Draw();
    
    cUnfoldedPseudodataPurityCorrectedOverAcceptanceVsGen25->Update();
    cUnfoldedPseudodataPurityCorrectedOverAcceptanceVsGen25->SaveAs("DDS10UnfoldedPseudodataPurityCorrectedOverAcceptanceVsGen25_Superchic_positron.pdf"); // MODIFIED FOR e+: Updated filename
    
    // Compare hUnfoldedRecoTrueOverAcceptance vs genTrue (25% subset, unscaled)
    TH1F *genTrue25_read = (TH1F*)f_pseudo_read->Get("genTrue");
    if (!genTrue25_read) {
        cerr << "Error: Could not read genTrue histogram from " << outfile_pseudo << endl;
        f_pseudo_read->Close();
        f_response.Close();
        return;
    }
    
    TCanvas *cUnfoldedRecoTrueOverAcceptanceVsGenTrue25 = new TCanvas("cUnfoldedRecoTrueOverAcceptanceVsGenTrue25", "Split MC Test: Unfolded/Acceptance vs Generator Truth (25% subset, True Events, Unscaled)", 900, 700);
    hUnfoldedRecoTrueOverAcceptance->SetLineColor(kBlue);
    hUnfoldedRecoTrueOverAcceptance->SetMarkerStyle(20);
    hUnfoldedRecoTrueOverAcceptance->SetLineWidth(2);
    hUnfoldedRecoTrueOverAcceptance->SetTitle("Split MC Test: Unfolded/Acceptance vs Generator Truth (25% subset, True Events, Unscaled)");
    hUnfoldedRecoTrueOverAcceptance->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
    hUnfoldedRecoTrueOverAcceptance->GetYaxis()->SetTitle("Events");
    
    genTrue25_read->SetLineColor(kRed);
    genTrue25_read->SetMarkerStyle(21);
    genTrue25_read->SetLineWidth(2);
    
    double max_genTrue25 = genTrue25_read->GetMaximum();
    max_y = TMath::Max(hUnfoldedRecoTrueOverAcceptance->GetMaximum(), max_genTrue25) * 1.2;
    hUnfoldedRecoTrueOverAcceptance->SetMaximum(max_y);
    
    hUnfoldedRecoTrueOverAcceptance->Draw("E");
    genTrue25_read->Draw("E SAME");
    
    TLegend *legUnfoldedVsGenTrue25 = new TLegend(0.65, 0.7, 0.9, 0.85);
    legUnfoldedVsGenTrue25->AddEntry(hUnfoldedRecoTrueOverAcceptance, "Unfolded RecoTrue / Acceptance (25% subset)", "lep");
    legUnfoldedVsGenTrue25->AddEntry(genTrue25_read, "Generator Truth (25% subset, True Events)", "lep");
    legUnfoldedVsGenTrue25->Draw();
    
    double chi2_unfolded_vs_genTrue25 = 0;
    int ndf_unfolded_vs_genTrue25 = 0;
    for (int i = 1; i <= hUnfoldedRecoTrueOverAcceptance->GetNbinsX(); i++) {
        double diff = hUnfoldedRecoTrueOverAcceptance->GetBinContent(i) - genTrue25_read->GetBinContent(i);
        double err1 = hUnfoldedRecoTrueOverAcceptance->GetBinError(i);
        double err2 = genTrue25_read->GetBinError(i);
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2_unfolded_vs_genTrue25 += (diff*diff)/(err*err);
            ndf_unfolded_vs_genTrue25++;
        }
    }
    
    TPaveText* pt_unfolded_vs_genTrue25 = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt_unfolded_vs_genTrue25->SetFillColor(0);
    pt_unfolded_vs_genTrue25->SetBorderSize(1);
    pt_unfolded_vs_genTrue25->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2_unfolded_vs_genTrue25, ndf_unfolded_vs_genTrue25));
    pt_unfolded_vs_genTrue25->AddText(Form("P-value = %.3f", TMath::Prob(chi2_unfolded_vs_genTrue25, ndf_unfolded_vs_genTrue25)));
    pt_unfolded_vs_genTrue25->Draw();
    
    cUnfoldedRecoTrueOverAcceptanceVsGenTrue25->Update();
    cUnfoldedRecoTrueOverAcceptanceVsGenTrue25->SaveAs("DDS10UnfoldedRecoTrueOverAcceptanceVsGenTrue25_Superchic_positron.pdf"); // MODIFIED FOR e+: Updated filename
    
    // Write response matrix and histograms
    f_response.cd();
    auto* hResponse = response.HresponseNoOverflow();
    if (hResponse) {
        hResponse->SetTitle("");
        hResponse->GetXaxis()->SetTitle("Superchic 2018 Gen-Level #Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
        hResponse->GetYaxis()->SetTitle("Superchic 2018 Reco-Level #Delta#phi = #phi_{ee} - #phi_{e+}"); // MODIFIED FOR e+: Updated axis label
        hResponse->Write("hResponseMatrix");
        response.Write("responseObject");
    }

    pseudodataPurityCorrected->Write();
    hUnfoldedPseudodataPurityCorrected->Write();
    hUnfoldedPseudodataPurityCorrectedOverAcceptance->Write();
    hUnfoldedRecoTrue->Write();
    hUnfoldedRecoTrueOverAcceptance->Write();
    hRefolded->Write();
    hGenMass->Write();
    hGenMassCuts->Write();
    hRecoMass->Write();
    hRecoMassCuts->Write();
    hreco->Write();
    hrecoTrue->Write();
    hrecoFake->Write();
    hgen->Write();
    hgenTrue->Write();
    hgenMiss->Write();
    hgenAll->Write();
    hEfficiency->Write();
    hPurity->Write();
    hAcceptance->Write();
    hRDeltaPhi->Write();
    hGDeltaPhi->Write();
    cSplitMCValidation->Write();
    cRefoldValidation->Write();
    cUnfoldedPseudodataPurityCorrectedOverAcceptanceVsGen25->Write();
    cUnfoldedRecoTrueOverAcceptanceVsGenTrue25->Write();
    SuperEventTree->Write();
    
    f_pseudo_read->Close();
    f_response.Close();
    
    // Final statistics
    cout << "\n=== SPLIT MC TEST SUMMARY ===" << endl;
    cout << "Pseudodata events (25% subset): " << pseudodata_count << endl;
    cout << "Pseudodata true events (25% subset): " << pseudodata_true_count << endl;
    cout << "Pseudodata gen events (25% subset): " << pseudodata_gen_count << endl;
    cout << "Response matrix true events (75% subset): " << response_true_count << endl;
    cout << "Response matrix fake events (75% subset): " << response_fake_count << endl;
    cout << "Response matrix missed events (75% subset): " << response_miss_count << endl;
    cout << "Refold χ²/ndf (Refolded vs Pseudodata): " << chi2_refold << "/" << ndf_refold << endl;
    cout << "Refold P-value (Refolded vs Pseudodata): " << TMath::Prob(chi2_refold, ndf_refold) << endl;
    cout << "Validation χ²/ndf (Unfolded Purity-Corrected vs Truth 75%): " << chi2_validation << "/" << ndf_validation << endl;
    cout << "Validation P-value (Unfolded Purity-Corrected vs Truth 75%): " << TMath::Prob(chi2_validation, ndf_validation) << endl;
    cout << "Validation χ²/ndf (Unfolded Purity-Corrected/Acceptance vs Gen 25%, Unscaled): " << chi2_unfolded_purity_vs_gen25 << "/" << ndf_unfolded_purity_vs_gen25 << endl;
    cout << "Validation P-value (Unfolded Purity-Corrected/Acceptance vs Gen 25%, Unscaled): " << TMath::Prob(chi2_unfolded_purity_vs_gen25, ndf_unfolded_purity_vs_gen25) << endl;
    cout << "Validation χ²/ndf (Unfolded RecoTrue/Acceptance vs GenTrue 25%, Unscaled): " << chi2_unfolded_vs_genTrue25 << "/" << ndf_unfolded_vs_genTrue25 << endl;
    cout << "Validation P-value (Unfolded RecoTrue/Acceptance vs GenTrue 25%, Unscaled): " << TMath::Prob(chi2_unfolded_vs_genTrue25, ndf_unfolded_vs_genTrue25) << endl;
    
    if (pseudodata_read->GetEntries() == 0) {
        cerr << "ERROR: pseudodata is empty! Check event selection logic." << endl;
    }
    if (hUnfoldedPseudodataPurityCorrected->GetEntries() == 0) {
        cerr << "ERROR: hUnfoldedPseudodataPurityCorrected is empty! Check unfolding." << endl;
    }
    if (hRefolded->GetEntries() == 0) {
        cerr << "ERROR: hRefolded is empty! Check refolding." << endl;
    }
    if (gen25_read->GetEntries() == 0) {
        cerr << "ERROR: gen (25% subset) is empty! Check generator-level selection logic." << endl;
    }
    if (genTrue25_read->GetEntries() == 0) {
        cerr << "ERROR: genTrue (25% subset) is empty! Check generator-level selection logic for true events." << endl;
    }

    cout << "Analysis completed successfully!" << endl;
    cout << "Pseudodata output saved to: " << outfile_pseudo << endl;
    cout << "Response output saved to: " << outfile_response << endl;
}

#ifndef __CINT__
int main() {
    ManuelPositronResponseMatrix(); // MODIFIED FOR e+: Updated function call
    return 0;
}
#endif
