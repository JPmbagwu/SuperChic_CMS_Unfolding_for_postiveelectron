#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
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
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1.h>
#include <TGaxis.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TLatex.h>
#endif

using namespace std;

// Debug flag - set to true for verbose output
const bool DEBUG = false;

// Centralized cut definitions
struct Cuts {
    static constexpr float MaxRapidity = 2.4;        // Individual electron |eta| < 2.4
    static constexpr float SingleElePtMin = 0.0;     // Single electron pT > 0.0 GeV/c
    static constexpr float maxDielePt = 1.0;        // Dielectron pT < 1.0 GeV/c
    static constexpr float minMass = 7.0;            // Minimum dielectron mass (GeV/c²)
    static constexpr float maxMass = 10.0;           // Maximum dielectron mass (GeV/c²)
    static constexpr float maxHFpCut = 7.3;          // Max HF+ energy (GeV)
    static constexpr float maxHFmCut = 7.6;          // Max HF- energy (GeV)
    static constexpr float maxVertexZ = 20.0;        // Max vertex z (cm)
    static constexpr float electronMass = 0.000511;  // Electron mass (GeV)
};

// Histogram settings
const float xrangeMax = TMath::Pi();
const float xrangeMin = -TMath::Pi();
const int xbinz = 18;

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

// Calculate the angle between two 2D vectors given by their (Px, Py) components
float GetAngle(float vPX1, float vPY1, float vPX2, float vPY2) {
    TVector2 DielectronVect(vPX1, vPY1);
    TVector2 ElectronVect(vPX2, vPY2);

    float DielectronMag = sqrt(DielectronVect.X() * DielectronVect.X() + DielectronVect.Y() * DielectronVect.Y());
    float ElectronMag = sqrt(ElectronVect.X() * ElectronVect.X() + ElectronVect.Y() * ElectronVect.Y());

    if (DielectronMag == 0 || ElectronMag == 0) return 0.0; // Avoid division by zero

    float cosAngle = (DielectronVect.X() * ElectronVect.X() + DielectronVect.Y() * ElectronVect.Y()) / (DielectronMag * ElectronMag);
    float sinAngle = (DielectronVect.X() * ElectronVect.Y() - DielectronVect.Y() * ElectronVect.X()) / (DielectronMag * ElectronMag);

    float tanAngle = atan2(sinAngle, cosAngle);
    return tanAngle;
}

void ManuelElectronResponseMatrix() {
    string inputMC = "/afs/cern.ch/user/j/jmbagwu/RooUnfold/examples/lblfiles/cms/my_file.txt";
    string outfile = "FDC10enhanced_superchic2018_dielectrons_full.root";
    
    // Check output directory exists
    if (gSystem->AccessPathName(gSystem->DirName(outfile.c_str()))) {
        cerr << "Error: Output directory does not exist: " << gSystem->DirName(outfile.c_str()) << endl;
        return;
    }
    
    TChain *t1 = CreateChain(inputMC, "ggHiNtuplizer/EventTree");
    if (t1->GetEntries() == 0) {
        cerr << "Error: No entries found in the chain!" << endl;
        return;
    }
    
    // Print available branches to debug
    if (DEBUG) {
        cout << "Available branches in MC tree:" << endl;
        t1->GetListOfBranches()->Print();
    }
    
    // Open output file
    TFile f(outfile.c_str(), "recreate");
    if (!f.IsOpen()) {
        cerr << "Error: Could not create output file " << outfile << endl;
        return;
    }
    f.cd(); // Set current directory to output file
    
    // Create the SuperEventTree
    TTree* SuperEventTree = new TTree("SuperEventTree", "Tree with Event Variables");
    
    // --- Generator level electron branches
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
    
    // --- Reconstructed level electron branches
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
    
    // --- Calorimeter branches
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
    
    // Event variables to store
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
    
    // Load delta phi histogram for electrons (data)
    TFile *dataFile = new TFile("/afs/cern.ch/user/j/jmbagwu/NEWSC/e+/delta_phi_histogram_electron_nozdc01.root", "READ");
    if (!dataFile || dataFile->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the original histogram (read-only)
    TH1F *hDeltaPhiData = (TH1F*)dataFile->Get("hDeltaPhi");
    if (!hDeltaPhiData) {
        std::cerr << "Histogram hDeltaPhi not found in file!" << std::endl;
        dataFile->Close();
        return;
    }

    // Clone histograms in the output file's directory
    f.cd();
    TH1F *hDeltaPhiRaw = (TH1F*)hDeltaPhiData->Clone("hDeltaPhiRaw");
    TH1F *hDeltaPhiPurityCorrected = (TH1F*)hDeltaPhiRaw->Clone("hDeltaPhiPurityCorrected");

    // Close input file
    dataFile->Close();
    delete dataFile;
    
    // Initialize response matrix
    RooUnfoldResponse response(xbinz, xrangeMin, xrangeMax);
    
    // Event counters
    int count_True = 0;
    int count_Fake = 0;
    int count_Miss = 0;
    int Event_count = 0;
    int GEventPass = 0;
    int REventPass = 0;
    
    // Histograms
    TH1F *hGenMass = new TH1F("hGenMass", "Generator Level Dielectron Mass", 200, 0, 100);
    TH1F *hGenMassCuts = new TH1F("hGenMassCuts", "Generator Level Dielectron Mass After Cuts", 200, 7.0, 100);
    TH1F *hRecoMass = new TH1F("hRecoMass", "Reconstructed Dielectron Mass", 200, 0, 100);
    TH1F *hRecoMassCuts = new TH1F("hRecoMassCuts", "Reconstructed Dielectron Mass After Cuts", 200, 7.0, 100);
    
    TH1F *hreco = new TH1F("hreco", "Reconstructed #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hrecoTrue = new TH1F("hrecoTrue", "True Reconstructed #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hrecoTrue_split1 = new TH1F("hrecoTrue_split1", "True Reco Split 1 (even events)", xbinz, xrangeMin, xrangeMax);
    TH1F *hrecoTrue_split2 = new TH1F("hrecoTrue_split2", "True Reco Split 2 (odd events)", xbinz, xrangeMin, xrangeMax);
    TH1F *hrecoFake = new TH1F("hrecoFake", "Fake Reconstructed #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hgen = new TH1F("hgen", "Generator Level #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hgenTrue = new TH1F("hgenTrue", "True Generator Level #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hgenMiss = new TH1F("hgenMiss", "Missed Generator Level #Delta#phi", xbinz, xrangeMin, xrangeMax);
    
    TH1F* hgenAll = new TH1F("hgenAll", "All Generator Level #Delta#phi (no cuts)", xbinz, xrangeMin, xrangeMax);
    
    // Dedicated histograms for RDeltaPhi and GDeltaPhi
    TH1F *hRDeltaPhi = new TH1F("hRDeltaPhi", "Reconstructed DeltaPhi Explicit", xbinz, xrangeMin, xrangeMax);
    TH1F *hGDeltaPhi = new TH1F("hGDeltaPhi", "Generator DeltaPhi Explicit", xbinz, xrangeMin, xrangeMax);
    
    // Efficiency, acceptance, and purity histograms
    TH1F* hEfficiency = new TH1F("hEfficiency", "Efficiency;#Delta#phi;Efficiency", xbinz, xrangeMin, xrangeMax);
    TH1F* hPurity = new TH1F("hPurity", "Purity;#Delta#phi;Purity", xbinz, xrangeMin, xrangeMax);
    TH1F* hAcceptance = new TH1F("hAcceptance", "Acceptance;#Delta#phi;Acceptance", xbinz, xrangeMin, xrangeMax);
    
    float RDeltaPhi = 0, GDeltaPhi = 0;
    int numberentry = t1->GetEntries();
    int split1_count = 0;
    int split2_count = 0;
    
    // Main event loop
    for (int j = 0; j < numberentry; j++) {
        t1->GetEntry(j);
        
        // Reset event flags
        int FakeEvt = 0;
        int TrueEvt = 0;
        int MissEvt = 0;
        int passGenCuts = 0;
        int passRecoCuts = 0;
        
        // Initialize indices to invalid values
        int iRp = -1, iRn = -1;
        int EleP = -1, EleN = -1;
        
        // Declare reconstructed variables at the start of event processing
        TLorentzVector RLep1, RLep2, RDiele;
        
        // --- Generator level pairing
        int idx_electron = -1;  // PID = 11 (e-)
        int idx_positron = -1;  // PID = -11 (e+)
        
        for (size_t k = 0; k < gen_elePID->size(); ++k) {
            int pid = gen_elePID->at(k);
            if (abs(pid) >= 20) continue; // Skip non-electron particles
            
            if (pid == 11 && idx_electron == -1) {
                idx_electron = k;  // Found electron (e-)
            }
            else if (pid == -11 && idx_positron == -1) {
                idx_positron = k;  // Found positron (e+)
            }
        }
        
        if (idx_electron == -1 || idx_positron == -1) continue;
        
        // Consistent assignment: GLep1 = electron, GLep2 = positron
        TLorentzVector GLep1, GLep2, GDiele;
        GLep1.SetPtEtaPhiM(gen_elePt->at(idx_electron), gen_eleEta->at(idx_electron),
                           gen_elePhi->at(idx_electron), Cuts::electronMass);   // Electron (e-)
        GLep2.SetPtEtaPhiM(gen_elePt->at(idx_positron), gen_eleEta->at(idx_positron),
                           gen_elePhi->at(idx_positron), Cuts::electronMass);  // Positron (e+)
        GDiele = GLep1 + GLep2;
        
        // Gen-level cuts
        bool genTwoelectron = (gen_nEle >= 2);
        bool genElectroncharge = (gen_elePID->at(idx_positron) + gen_elePID->at(idx_electron) == 0);
        bool GRapidity = abs(GLep1.Eta()) <= Cuts::MaxRapidity && abs(GLep2.Eta()) <= Cuts::MaxRapidity;
        bool GMassCut = (GDiele.M() > Cuts::minMass && GDiele.M() < Cuts::maxMass);
        bool GpTcut = (GDiele.Pt() < Cuts::maxDielePt);
        bool GVertexCut = (fabs(gen_trkvz->at(idx_positron)) <= Cuts::maxVertexZ &&
                           fabs(gen_trkvz->at(idx_electron)) <= Cuts::maxVertexZ);
        bool GenSingleElectronPtCut = GLep1.Pt() > Cuts::SingleElePtMin && GLep2.Pt() > Cuts::SingleElePtMin;
        
        // Compute DeltaPhi using the positron (GLep2, e+)
        float GenDielePx = GDiele.Pt() * cos(GDiele.Phi());
        float GenDielePy = GDiele.Pt() * sin(GDiele.Phi());
        float GenElePx = GLep2.Pt() * cos(GLep2.Phi());
        float GenElePy = GLep2.Pt() * sin(GLep2.Phi());
        GDeltaPhi = GetAngle(GenDielePx, GenDielePy, GenElePx, GenElePy);
        
        hgenAll->Fill(GDeltaPhi);  // Fill for ALL generated events
        
        hGenMass->Fill(GDiele.M());
        
        if (GMassCut && GRapidity && GpTcut && genTwoelectron && genElectroncharge && GVertexCut && GenSingleElectronPtCut) {
            passGenCuts = 1;
            GEventPass++;
            hgen->Fill(GDeltaPhi);
            hGDeltaPhi->Fill(GDeltaPhi);
            hGenMassCuts->Fill(GDiele.M());
        }
        
        // --- Reco-level selections ---
        bool Twoelectron = (reco_nEle == 2);
        bool Twotrk = (reco_nTrk == 2);
        
        if (Twoelectron && Twotrk) {
            if (reco_eleCharge->size() < 2) continue;
            bool Electroncharge = (reco_eleCharge->at(0) * reco_eleCharge->at(1) == -1);
            if (!Electroncharge) continue;
            
            // Vertex selection
            if (reco_trkvz->size() == 0) continue;
            
            // Assign electrons based on charge
            EleP = ((*reco_eleCharge)[0] == 1) ? 0 : 1;  // Positron (e+)
            EleN = 1 - EleP;                             // Electron (e-)
            
            // Construct TLorentzVectors
            RLep1.SetPtEtaPhiM((*reco_elePt)[EleP], (*reco_eleEta)[EleP], (*reco_elePhi)[EleP], Cuts::electronMass);  // Positron (e+)
            RLep2.SetPtEtaPhiM((*reco_elePt)[EleN], (*reco_eleEta)[EleN], (*reco_elePhi)[EleN], Cuts::electronMass);  // Electron (e-)
            RDiele = RLep1 + RLep2;
            
            // Compute DeltaPhi using the positron (RLep1, e+)
            float RecoDielePx = RDiele.Pt() * cos(RDiele.Phi());
            float RecoDielePy = RDiele.Pt() * sin(RDiele.Phi());
            float RecoElePx = RLep1.Pt() * cos(RLep1.Phi());
            float RecoElePy = RLep1.Pt() * sin(RLep1.Phi());
            RDeltaPhi = GetAngle(RecoDielePx, RecoDielePy, RecoElePx, RecoElePy);
            
            hRecoMass->Fill(RDiele.M());
            
            // Calculate max HF energies and max tower properties
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
            
            // Reco-level cuts
            bool RRapidity = abs(RLep1.Eta()) <= Cuts::MaxRapidity && abs(RLep2.Eta()) <= Cuts::MaxRapidity;
            bool RVertexCut = fabs((*reco_trkvz)[0]) <= Cuts::maxVertexZ;
            bool RMassCut = (RDiele.M() > Cuts::minMass && RDiele.M() < Cuts::maxMass);
            bool RpTCut = (RDiele.Pt() < Cuts::maxDielePt);
            bool HFCuts = (maxHFm <= Cuts::maxHFmCut && maxHFp <= Cuts::maxHFpCut);
            bool RecoSingleElectronPtCut = RLep1.Pt() > Cuts::SingleElePtMin && RLep2.Pt() > Cuts::SingleElePtMin;
            
            if (RMassCut && HFCuts && RVertexCut && RRapidity && RecoSingleElectronPtCut && RpTCut) {
                hreco->Fill(RDeltaPhi);
                hRDeltaPhi->Fill(RDeltaPhi);
                hRecoMassCuts->Fill(RDiele.M());
                REventPass++;
                passRecoCuts = 1;
            }
        }
        
        // --- Event classification
        if (passGenCuts == 1 && passRecoCuts == 0) {
            MissEvt = 1;
        }
        if (passGenCuts == 0 && passRecoCuts == 1) {
            FakeEvt = 1;
        }
        if (passGenCuts == 1 && passRecoCuts == 1) {
            TrueEvt = 1;
        }
        
        if (DEBUG && TrueEvt) {
            cout << "GDeltaPhi is " << GDeltaPhi << endl;
            cout << "RDeltaPhi is " << RDeltaPhi << endl;
            cout << "G lep 1 charge (electron) is " << gen_elePID->at(idx_electron) << endl;
            cout << "G lep 2 charge (positron) is " << gen_elePID->at(idx_positron) << endl;
            if (EleP >= 0 && EleN >= 0 && EleP < reco_eleCharge->size() && EleN < reco_eleCharge->size()) {
                cout << "R lep1 charge (positron) is " << reco_eleCharge->at(EleP) << endl;
                cout << "R lep2 charge (electron) is " << reco_eleCharge->at(EleN) << endl;
                cout << "R lep1 phi (positron) is " << RLep1.Phi() << endl;
                cout << "RDielectron phi is " << RDiele.Phi() << endl;
            }
            cout << "G lep2 phi (positron) is " << GLep2.Phi() << endl;
            cout << "GDielectron phi is " << GDiele.Phi() << endl;
        }
        
        if (DEBUG && TrueEvt) {
            cout << "=== Event " << j << " ===" << endl;
            cout << "Gen: electron idx=" << idx_electron << " PID=" << gen_elePID->at(idx_electron) << endl;
            cout << "Gen: positron idx=" << idx_positron << " PID=" << gen_elePID->at(idx_positron) << endl;
            if (EleP >= 0 && EleN >= 0 && EleP < reco_eleCharge->size() && EleN < reco_eleCharge->size()) {
                cout << "Reco: positron idx=" << EleP << " charge=" << reco_eleCharge->at(EleP) << endl;
                cout << "Reco: electron idx=" << EleN << " charge=" << reco_eleCharge->at(EleN) << endl;
            } else {
                cout << "Reco: No valid positron or electron indices" << endl;
            }
            cout << "GDeltaPhi = " << GDeltaPhi << ", RDeltaPhi = " << RDeltaPhi << endl;
        }
        
        if (MissEvt == 1) {
            count_Miss++;
            hgenMiss->Fill(GDeltaPhi);
        }
        if (FakeEvt == 1) {
            count_Fake++;
            hrecoFake->Fill(RDeltaPhi);
        }
        if (TrueEvt == 1) {
            count_True++;
            response.Fill(RDeltaPhi, GDeltaPhi);
            hrecoTrue->Fill(RDeltaPhi);
            hgenTrue->Fill(GDeltaPhi);
            if (j % 2 == 0) {
                split1_count++;
                hrecoTrue_split1->Fill(RDeltaPhi);
            } else {
                split2_count++;
                hrecoTrue_split2->Fill(RDeltaPhi);
            }
        }
        
        Event_count++;
        SuperEventTree->Fill();
    } // end loop over entries
    
    if (DEBUG) {
        cout << "Split1 entries: " << split1_count << endl;
        cout << "Split2 entries: " << split2_count << endl;
        cout << "True events: " << count_True << endl;
        cout << "Fake events: " << count_Fake << endl;
        cout << "Missed events: " << count_Miss << endl;
    }
    
    // Calculate acceptance: genTrue / (genTrue + genMiss)
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
    
    // Apply purity corrections: purity = true / (true + fake)
    for (int iBinX = 1; iBinX <= hDeltaPhiPurityCorrected->GetNbinsX(); iBinX++) {
        double trueCounts = hrecoTrue->GetBinContent(iBinX);
        double fakeCounts = hrecoFake->GetBinContent(iBinX);
        double denom = trueCounts + fakeCounts;
        
        if (denom > 0) {
            double purity = trueCounts / denom;
            double oldContent = hDeltaPhiRaw->GetBinContent(iBinX);
            double oldError = hDeltaPhiRaw->GetBinError(iBinX);
            
            hDeltaPhiPurityCorrected->SetBinContent(iBinX, oldContent * purity);
            
            double purError = sqrt(purity * (1 - purity) / denom);
            double propagatedError = sqrt(pow(oldContent * purError, 2) + pow(purity * oldError, 2));
            hDeltaPhiPurityCorrected->SetBinError(iBinX, propagatedError);
        } else {
            hDeltaPhiPurityCorrected->SetBinContent(iBinX, 0);
            hDeltaPhiPurityCorrected->SetBinError(iBinX, 0);
        }
    }
    
    // Calculate and fill purity histogram with binomial error propagation
    for (int iBinX = 1; iBinX <= hreco->GetNbinsX(); ++iBinX) {
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
    
    // Diagnostic: Compare hDeltaPhiRaw with hrecoTrue
    TCanvas *cCompareDataMC = new TCanvas("cCompareDataMC", "Data vs MC RecoTrue", 800, 600);
    hDeltaPhiRaw->SetLineColor(kBlue);
    hDeltaPhiRaw->SetMarkerStyle(20);
    hDeltaPhiRaw->SetTitle("Data vs MC Reconstructed #Delta#phi");
    hDeltaPhiRaw->GetXaxis()->SetTitle("#Delta#phi");
    hDeltaPhiRaw->Draw("E");
    hrecoTrue->SetLineColor(kRed);
    hrecoTrue->SetMarkerStyle(21);
    hrecoTrue->Draw("E SAME");
    TLegend *legCompare = new TLegend(0.7, 0.7, 0.9, 0.9);
    legCompare->AddEntry(hDeltaPhiRaw, "Data #Delta#phi", "lep");
    legCompare->AddEntry(hrecoTrue, "MC True #Delta#phi", "lep");
    legCompare->Draw();
    cCompareDataMC->Write();
    
    // --- Split test check on hrecoTrue
    TCanvas* cSplitTest = new TCanvas("cSplitTest", "Split Test Check", 800, 600);
    hrecoTrue_split1->SetLineColor(kBlue);
    hrecoTrue_split2->SetLineColor(kRed);
    hrecoTrue_split1->SetMarkerStyle(20);
    hrecoTrue_split2->SetMarkerStyle(21);
    
    if (hrecoTrue_split1->Integral() > 0) hrecoTrue_split1->Scale(hrecoTrue->Integral() / hrecoTrue_split1->Integral());
    if (hrecoTrue_split2->Integral() > 0) hrecoTrue_split2->Scale(hrecoTrue->Integral() / hrecoTrue_split2->Integral());
    
    hrecoTrue->SetLineColor(kBlack);
    hrecoTrue->SetLineWidth(2);
    hrecoTrue->Draw("HIST");
    hrecoTrue_split1->Draw("HIST SAME");
    hrecoTrue_split2->Draw("HIST SAME");
    
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(hrecoTrue, "2018 Superchic recoTrue", "l");
    leg->AddEntry(hrecoTrue_split1, "Split 1 (recoTrue)", "l");
    leg->AddEntry(hrecoTrue_split2, "Split 2 (recoTrue)", "l");
    leg->Draw();
    
    double chi2 = 0;
    int ndf = 0;
    for (int i = 1; i <= hrecoTrue_split1->GetNbinsX(); i++) {
        double diff = hrecoTrue_split1->GetBinContent(i) - hrecoTrue_split2->GetBinContent(i);
        double err1 = hrecoTrue_split1->GetBinError(i);
        double err2 = hrecoTrue_split2->GetBinError(i);
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2 += (diff*diff)/(err*err);
            ndf++;
        }
    }
    
    TPaveText* pt1 = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt1->SetFillColor(0);
    pt1->SetBorderSize(1);
    pt1->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf));
    pt1->AddText(Form("P-value = %.3f", TMath::Prob(chi2, ndf)));
    pt1->Draw();
    
    cSplitTest->Update();
    
    // --- Unfolding with RooUnfoldBayes for data
    RooUnfoldBayes unfold_bayes_data(&response, hDeltaPhiPurityCorrected, 2);
    TH1* UnfoldData = unfold_bayes_data.Hunfold();
    UnfoldData->SetName("UnfoldData2018");
    UnfoldData->SetTitle("Unfold 2018 Data");
    UnfoldData->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}");
    
    // Forward fold the unfolded data
    TH1* hForwardFoldedData = (TH1*)response.ApplyToTruth(UnfoldData, "hForwardFoldedData");
    hForwardFoldedData->SetTitle("Forward Fold Unfolded 2018 Data");
    hForwardFoldedData->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}");
    
    // Create ratio of forward folded data to original data
    TH1* hRatioData = (TH1*)hForwardFoldedData->Clone("hRatioData");
    hRatioData->Divide(hDeltaPhiPurityCorrected);
    hRatioData->SetTitle("Ratio: Forward Fold Data / 2018 Data");
    hRatioData->GetYaxis()->SetTitle("Ratio");
    
    // --- Monte Carlo closure test
    RooUnfoldBayes unfold_bayes_mc(&response, hrecoTrue, 2);
    TH1* UnfoldMC = unfold_bayes_mc.Hunfold();
    UnfoldMC->SetName("UnfoldSCMC2018");
    UnfoldMC->SetTitle("Unfold 2018 Superchic MC True");
    UnfoldMC->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}");
    
    // Forward fold the unfolded MC
    TH1* hForwardFoldedMC = (TH1*)response.ApplyToTruth(UnfoldMC, "hForwardFoldedMC");
    hForwardFoldedMC->SetTitle("Forward Fold 2018 Superchic MC True");
    hForwardFoldedMC->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}");
    
    // Create ratio of forward folded MC to original MC
    TH1* hRatioMC = (TH1*)hForwardFoldedMC->Clone("hRatioMC");
    hRatioMC->Divide(hrecoTrue);
    hRatioMC->SetTitle("Ratio: Forward Fold 2018 Superchic MC / True Reco 2018 Superchic MC");
    hRatioMC->GetYaxis()->SetTitle("Ratio");
    
    // --- Create UnfoldData / Acceptance plot
    TH1F* hUnfoldDataOverAcceptance = (TH1F*)UnfoldData->Clone("hUnfoldDataOverAcceptance");
    hUnfoldDataOverAcceptance->SetTitle("Unfolded 2018 Data / Acceptance");
    hUnfoldDataOverAcceptance->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e+}");
    hUnfoldDataOverAcceptance->GetYaxis()->SetTitle("Unfolded Data / Acceptance");
    
    // Compute UnfoldData / hAcceptance
    for (int iBinX = 1; iBinX <= hUnfoldDataOverAcceptance->GetNbinsX(); ++iBinX) {
        double unfoldVal = UnfoldData->GetBinContent(iBinX);
        double acceptance = hAcceptance->GetBinContent(iBinX);
        double unfoldErr = UnfoldData->GetBinError(iBinX);
        double acceptanceErr = hAcceptance->GetBinError(iBinX);
        if (acceptance > 0) {
            double ratio = unfoldVal / acceptance;
            double relErr = sqrt(pow(unfoldErr/unfoldVal, 2) + pow(acceptanceErr/acceptance, 2));
            hUnfoldDataOverAcceptance->SetBinContent(iBinX, ratio);
            hUnfoldDataOverAcceptance->SetBinError(iBinX, ratio * relErr);
        } else {
            hUnfoldDataOverAcceptance->SetBinContent(iBinX, 0);
            hUnfoldDataOverAcceptance->SetBinError(iBinX, 0);
        }
    }
    
    // Create canvas for UnfoldData / Acceptance
    TCanvas *cUnfoldDataOverAcceptance = new TCanvas("cUnfoldDataOverAcceptance", "Unfolded Data / Acceptance", 900, 700);
    hUnfoldDataOverAcceptance->SetLineColor(kBlue);
    hUnfoldDataOverAcceptance->SetMarkerStyle(20);
    hUnfoldDataOverAcceptance->SetMarkerColor(kBlue);
    hUnfoldDataOverAcceptance->SetLineWidth(2);
    hUnfoldDataOverAcceptance->Draw("E");
    
    TLegend *legUnfoldAcc = new TLegend(0.65, 0.7, 0.9, 0.85);
    legUnfoldAcc->AddEntry(hUnfoldDataOverAcceptance, "Unfolded Data / Acceptance", "lep");
    legUnfoldAcc->Draw();
    
    // Add CMS and analysis labels
    TLatex latexUnfoldAcc;
    latexUnfoldAcc.SetNDC();
    latexUnfoldAcc.SetTextFont(62);
    latexUnfoldAcc.SetTextSize(0.04);
    latexUnfoldAcc.SetTextAlign(13);
    latexUnfoldAcc.DrawLatex(0.10, 0.93, "CMS");
    latexUnfoldAcc.SetTextFont(42);
    latexUnfoldAcc.SetTextSize(0.035);
    latexUnfoldAcc.DrawLatex(0.16, 0.93, "#it{work in progress}");
    latexUnfoldAcc.SetTextSize(0.038);
    latexUnfoldAcc.SetTextAlign(33);
    latexUnfoldAcc.DrawLatex(0.89, 0.935, "PbPb 2018 #sqrt{#it{s}_{NN}} = 5.02 TeV");
    latexUnfoldAcc.SetTextAlign(31);
    latexUnfoldAcc.SetTextSize(0.030);
    latexUnfoldAcc.DrawLatex(0.84, 0.24, "p_{T}_{ee} < 1.0 GeV");
    latexUnfoldAcc.DrawLatex(0.84, 0.20, "p_{T}_{e} > 0.0 GeV");
    latexUnfoldAcc.DrawLatex(0.84, 0.17, "|#eta_{e}| < 2.4");
    latexUnfoldAcc.DrawLatex(0.84, 0.13, "7.0 GeV < M_{ee} < 10.0 GeV");
    
    cUnfoldDataOverAcceptance->Update();
    cUnfoldDataOverAcceptance->SaveAs("FDC10UnfoldDataOverAcceptance.pdf");
    
    // Bottom-line test plotting
    TCanvas *cDataUnfold = new TCanvas("cDataUnfold", "Bottom Line Test - Data", 900, 700);
    if (!hDeltaPhiPurityCorrected || !hForwardFoldedData) {
        std::cerr << "One or both histograms are null!" << std::endl;
        return;
    }
    
    std::cout << "Corrected Integral: " << hDeltaPhiPurityCorrected->Integral() << std::endl;
    std::cout << "Forward Folded Integral: " << hForwardFoldedData->Integral() << std::endl;
    
    hDeltaPhiPurityCorrected->SetStats(kFALSE);
    hDeltaPhiPurityCorrected->SetLineColor(kBlack);
    hDeltaPhiPurityCorrected->SetMarkerStyle(20);
    hDeltaPhiPurityCorrected->SetLineWidth(2);
    hDeltaPhiPurityCorrected->SetTitle("2018 Data vs Forward Fold with Superchic 2018");
    hDeltaPhiPurityCorrected->GetXaxis()->SetTitle("#Delta#phi");
    
    double maxY = std::max(hDeltaPhiPurityCorrected->GetMaximum(), hForwardFoldedData->GetMaximum());
    if (maxY > 0)
        hDeltaPhiPurityCorrected->SetMaximum(1.2 * maxY);
    
    hDeltaPhiPurityCorrected->Draw("E");
    hForwardFoldedData->SetLineColor(kRed);
    hForwardFoldedData->SetMarkerStyle(21);
    hForwardFoldedData->SetLineWidth(2);
    hForwardFoldedData->Draw("E SAME");
    
    TLegend *legData = new TLegend(0.65, 0.7, 0.9, 0.85);
    legData->AddEntry(hDeltaPhiPurityCorrected, "2018 Data", "lep");
    legData->AddEntry(hForwardFoldedData, "Forward Fold", "lep");
    legData->Draw();
    
    // Calculate chi2 and ndf for hDeltaPhiPurityCorrected vs hForwardFoldedData
    double chi2_data = 0;
    int ndf_data = 0;
    for (int i = 1; i <= hDeltaPhiPurityCorrected->GetNbinsX(); i++) {
        double diff = hDeltaPhiPurityCorrected->GetBinContent(i) - hForwardFoldedData->GetBinContent(i);
        double err1 = hDeltaPhiPurityCorrected->GetBinError(i);
        double err2 = hForwardFoldedData->GetBinError(i);
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2_data += (diff*diff)/(err*err);
            ndf_data++;
        }
    }
    
    TPaveText* pt2 = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt2->SetFillColor(0);
    pt2->SetBorderSize(1);
    pt2->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2_data, ndf_data));
    pt2->AddText(Form("P-value = %.3f", TMath::Prob(chi2_data, ndf_data)));
    pt2->Draw();
    
    cDataUnfold->Update();
    cDataUnfold->SaveAs("FDC10bottom_line_test.pdf");
    
    // Canvas for data ratio
    TCanvas *cDataRatio = new TCanvas("cDataRatio", "2018 Data Ratio ForwardFolded/Corrected", 900, 700);
    hRatioData->SetLineColor(kBlue);
    hRatioData->SetMarkerStyle(20);
    hRatioData->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatioData->Draw("E");
    cDataRatio->Update();
    
    // Canvas for MC closure test
    TCanvas *cMCUnfold = new TCanvas("cMCUnfold", "Bottom Line Test - MC Closure", 900, 700);
    hrecoTrue->SetLineColor(kBlack);
    hrecoTrue->SetMarkerStyle(20);
    hrecoTrue->SetTitle("2018 Superchic MC: recoTrue Split Test Check");
    hrecoTrue->GetXaxis()->SetTitle("#Delta#phi");
    hrecoTrue->Draw("E");
    
    hForwardFoldedMC->SetLineColor(kRed);
    hForwardFoldedMC->SetMarkerStyle(21);
    hForwardFoldedMC->Draw("E SAME");
    
    TLegend *legMC = new TLegend(0.65, 0.7, 0.9, 0.85);
    legMC->AddEntry(hrecoTrue, "Reco True MC", "lep");
    legMC->AddEntry(hForwardFoldedMC, "Forward Folded Unfolded 2018 Superchic MC", "lep");
    legMC->Draw();
    
    cMCUnfold->Update();
    
    // Canvas for MC ratio
    TCanvas *cMCRatio = new TCanvas("cMCRatio", "MC Ratio ForwardFolded/RecoTrue", 900, 700);
    hRatioMC->SetLineColor(kBlue);
    hRatioMC->SetMarkerStyle(20);
    hRatioMC->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatioMC->Draw("E");
    cMCRatio->Update();
    
    // New comparison: Unfolded Data vs Gen MC Truth (fiducial level)
    TCanvas *cUnfoldVsGen = new TCanvas("cUnfoldVsGen", "Unfolded Data vs Gen MC Truth", 900, 700);
    UnfoldData->SetLineColor(kBlue);
    UnfoldData->SetMarkerStyle(20);
    hgen->SetLineColor(kRed);
    hgen->SetMarkerStyle(21);
    UnfoldData->DrawNormalized("E");
    hgen->DrawNormalized("SAME E");

    // Add legend
    TLegend *legUnfoldVsGen = new TLegend(0.65, 0.7, 0.9, 0.85);
    legUnfoldVsGen->AddEntry(UnfoldData, "Unfolded Data", "lep");
    legUnfoldVsGen->AddEntry(hgen, "Gen MC Truth (hgen)", "lep");
    legUnfoldVsGen->Draw();

    // Add χ²/ndf p-value box
    double chi2_unfold = 0;
    int ndf_unfold = 0;
    for (int i = 1; i <= UnfoldData->GetNbinsX(); i++) {
        double content1 = UnfoldData->GetBinContent(i);
        double content2 = hgen->GetBinContent(i);
        double totalIntegral1 = UnfoldData->Integral();
        double totalIntegral2 = hgen->Integral();
        if (totalIntegral1 > 0) content1 /= totalIntegral1;
        if (totalIntegral2 > 0) content2 /= totalIntegral2;
        
        double diff = content1 - content2;
        double err1 = UnfoldData->GetBinError(i) / totalIntegral1;
        double err2 = hgen->GetBinError(i) / totalIntegral2;
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2_unfold += (diff*diff)/(err*err);
            ndf_unfold++;
        }
    }

    TPaveText* pt_unfold = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt_unfold->SetFillColor(0);
    pt_unfold->SetBorderSize(1);
    pt_unfold->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2_unfold, ndf_unfold));
    pt_unfold->AddText(Form("P-value = %.3f", TMath::Prob(chi2_unfold, ndf_unfold)));
    pt_unfold->Draw();

    cUnfoldVsGen->Update();
    cUnfoldVsGen->Write();
    cUnfoldVsGen->SaveAs("FDC10UnfoldVsGenTruth.pdf");
    
    // Canvas for GDeltaPhi vs RDeltaPhi comparison
    TCanvas *cDeltaPhi = new TCanvas("cDeltaPhi", "GDeltaPhi vs RDeltaPhi", 800, 600);
    hgen->SetLineColor(kBlue);
    hgen->SetMarkerStyle(20);
    hgen->SetTitle("Generator vs Reconstructed #Delta#phi");
    hgen->GetXaxis()->SetTitle("#Delta#phi");
    hgen->Draw("E");
    hreco->SetLineColor(kRed);
    hreco->SetMarkerStyle(21);
    hreco->Draw("E SAME");
    TLegend *legDeltaPhi = new TLegend(0.7, 0.7, 0.9, 0.9);
    legDeltaPhi->AddEntry(hgen, "Generator #Delta#phi", "lep");
    legDeltaPhi->AddEntry(hreco, "Reconstructed #Delta#phi", "lep");
    legDeltaPhi->Draw();
    
    // Add chi2 test for hgen vs hreco
    double chi2_deltaPhi = 0;
    int ndf_deltaPhi = 0;
    for (int i = 1; i <= hgen->GetNbinsX(); i++) {
        double diff = hgen->GetBinContent(i) - hreco->GetBinContent(i);
        double err1 = hgen->GetBinError(i);
        double err2 = hreco->GetBinError(i);
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2_deltaPhi += (diff*diff)/(err*err);
            ndf_deltaPhi++;
        }
    }
    TPaveText* pt_deltaPhi = new TPaveText(0.15, 0.65, 0.45, 0.75, "NDC");
    pt_deltaPhi->SetFillColor(0);
    pt_deltaPhi->SetBorderSize(1);
    pt_deltaPhi->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2_deltaPhi, ndf_deltaPhi));
    pt_deltaPhi->AddText(Form("P-value = %.3f", TMath::Prob(chi2_deltaPhi, ndf_deltaPhi)));
    pt_deltaPhi->Draw();
    cDeltaPhi->Update();

    // --- Write all objects to output file
    f.cd();
    
    // Response matrix
    auto* hResponse = response.HresponseNoOverflow();
    if (hResponse) {
        hResponse->SetTitle("");
        hResponse->GetXaxis()->SetTitle("Superchic 2018 Gen-Level #Delta#phi = #phi_{ee} - #phi_{e+}");
        hResponse->GetYaxis()->SetTitle("Superchic 2018 Reco-Level #Delta#phi = #phi_{ee} - #phi_{e+}");
        hResponse->Write("hResponseMatrix");
        response.Write("responseObject");

        TCanvas *c = new TCanvas("cResponse", "Response Matrix", 800, 600);
        gStyle->SetOptStat();
        hResponse->Draw("COLZ");

        gPad->Update();
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
        latex.DrawLatex(0.84, 0.24, "p_{T}_{ee} < 1.0 GeV");
//        latex.DrawLatex(0.84, 0.20, "p_{T}_{e} > 0.0 GeV");
        latex.DrawLatex(0.84, 0.17, "|#eta_{e}| < 2.4");
        latex.DrawLatex(0.84, 0.13, "7.0 GeV < M_{ee} < 10.0 GeV");

        c->Update();
        c->Draw();
        c->Print("FDC10SuperchicResponseMatrix.pdf");
        delete c;
    }

    // Ensure correct axis titles before writing
    hDeltaPhiPurityCorrected->GetXaxis()->SetTitle("#Delta#phi");
    hDeltaPhiPurityCorrected->Write();
    hDeltaPhiRaw->Write();
    hForwardFoldedData->GetXaxis()->SetTitle("#Delta#phi");

    // Data unfolding results
    UnfoldData->Write();
    hForwardFoldedData->Write();
    hRatioData->Write();

    // MC unfolding results
    UnfoldMC->Write();
    hForwardFoldedMC->Write();
    hRatioMC->Write();

    // Split test results
    hrecoTrue->Write();
    hrecoTrue_split1->Write();
    hrecoTrue_split2->Write();
    cSplitTest->Write();

    // Other histograms
    hGenMass->Write();
    hGenMassCuts->Write();
    hRecoMass->Write();
    hRecoMassCuts->Write();
    hreco->Write();
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
    hUnfoldDataOverAcceptance->Write();
    cUnfoldDataOverAcceptance->Write();
    
    // Write comparison canvas
    cDeltaPhi->Write();

    // Event tree
    SuperEventTree->Write();
    
    // Check for empty histograms
    if (hgen->GetEntries() == 0) cerr << "Warning: hgen is empty!" << endl;
    if (hreco->GetEntries() == 0) cerr << "Warning: hreco is empty!" << endl;
    if (hGDeltaPhi->GetEntries() == 0) cerr << "Warning: hGDeltaPhi is empty!" << endl;
    if (hRDeltaPhi->GetEntries() == 0) cerr << "Warning: hRDeltaPhi is empty!" << endl;

    f.Close();

    cout << "Analysis completed successfully!" << endl;
    cout << "Output saved to: " << outfile << endl;
    cout << "True events: " << count_True << endl;
    cout << "Fake events: " << count_Fake << endl;
    cout << "Missed events: " << count_Miss << endl;
    cout << "Generator events passing cuts: " << GEventPass << endl;
    cout << "Reconstructed events passing cuts: " << REventPass << endl;
}

#ifndef __CINT__
int main() {
    ManuelElectronResponseMatrix();
    return 0;
}
#endif
