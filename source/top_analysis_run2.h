//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jan 22 11:35:09 2023 by ROOT version 6.26/10
// from TTree nominal/tree
// found on file: data16_13TeV_merged.root
//////////////////////////////////////////////////////////

#ifndef top_analysis_run2_h
#define top_analysis_run2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TSelector.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TCanvas.h"


// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class top_analysis_run2  : public TSelector{
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
    TH1F * h_t_mass = 0;
    TH1F * h_t_pt = 0;
    TH1F * h_t_eta = 0;
    TH1F * h_t_phi = 0;

    TH1F * h_tbar_mass = 0;
    TH1F * h_tbar_pt = 0;
    TH1F * h_tbar_eta = 0;
    TH1F * h_tbar_phi = 0;


    TH1F * h_ttbar_mass = 0;
    TH1F * h_ttbar_pt = 0;
    TH1F * h_ttbar_eta = 0;
    TH1F * h_ttbar_phi = 0;


     TH1F * h_j1_pt = 0;
     TH1F * h_j2_pt = 0;
     TH1F * h_j1j2_pt = 0;
     TH1F * h_j1j2_m = 0;

     TH1F * h_dilep_m = 0;  

     TH1F * h_lep_pos_pt = 0;
     TH1F * h_lep_neg_pt = 0;

     TH1F *  h_proton_side1_xi = 0;
     TH1F *  h_proton_side2_xi = 0;

   // Declaration of leaf types
   // Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxreco_proton = 33;

   // Declaration of leaf types
   ULong64_t       eventNumber;
   UInt_t          runNumber;
   UInt_t          mcChannelNumber;
   Float_t         mu;
   Float_t         mu_actual;
   Float_t         mu_original_xAOD;
   Float_t         mu_actual_original_xAOD;
   UInt_t          backgroundFlags;
   UInt_t          hasBadMuon;
   Int_t           allhad_2017;
   Int_t           ejets_2017;
   Int_t           mujets_2017;
   Int_t           ee_2017;
   Int_t           mumu_2017;
   Int_t           emu_2017;
   Char_t          HLT_mu50;
   Char_t          HLT_mu26_ivarmedium;
   Char_t          HLT_mu15_L1MU10;
   Char_t          HLT_e26_lhtight_nod0_ivarloose;
   Char_t          HLT_mu14;
   Char_t          HLT_e17_lhloose_nod0;
   Char_t          HLT_2mu4;
   Char_t          HLT_j175;
   Char_t          HLT_j110;
   Char_t          HLT_e140_lhloose_nod0;
   Char_t          HLT_j260;
   Char_t          HLT_j85;
   Char_t          HLT_e60_lhmedium_nod0;
   Char_t          HLT_e15_lhloose_nod0_L1EM12;
   Char_t          HLT_e14_lhtight_nod0;
   Bool_t          reco_fakeEvent;
   vector<float>   *reco_lep_pt;
   vector<float>   *reco_lep_eta;
   vector<float>   *reco_lep_phi;
   vector<float>   *reco_lep_e;
   vector<int>     *reco_lep_pdgid;
   Float_t         reco_dilep_m;
   Float_t         reco_b_pt;
   Float_t         reco_b_m;
   Float_t         reco_b_phi;
   Float_t         reco_b_eta;
   Float_t         reco_bbar_pt;
   Float_t         reco_bbar_m;
   Float_t         reco_bbar_phi;
   Float_t         reco_bbar_eta;
   Float_t         reco_jet1_pt;
   Float_t         reco_jet1_eta;
   Float_t         reco_jet1_phi;
   Float_t         reco_jet1_e;
   Float_t         reco_jet2_pt;
   Float_t         reco_jet2_eta;
   Float_t         reco_jet2_phi;
   Float_t         reco_jet2_e;
   Float_t         reco_met_ex;
   Float_t         reco_met_ey;
   Int_t           reco_nb;
   Int_t           reco_nbbar;
   Int_t           reco_njets;
   Int_t           reco_proton_;
   Double_t        reco_proton_fCoordinates_fPt[kMaxreco_proton];   //[reco_proton_]
   Double_t        reco_proton_fCoordinates_fEta[kMaxreco_proton];   //[reco_proton_]
   Double_t        reco_proton_fCoordinates_fPhi[kMaxreco_proton];   //[reco_proton_]
   Double_t        reco_proton_fCoordinates_fE[kMaxreco_proton];   //[reco_proton_]
   vector<float>   *reco_proton_chi2;
   vector<int>     *reco_proton_side;
   vector<int>     *reco_proton_methodID;
   vector<float>   *reco_proton_xi;
   vector<int>     *reco_proton_ntracks;
   vector<vector<int> > *reco_proton_stations;
   Bool_t          reco_tau_event;
   Int_t           reco_top_reco_method;
   Float_t         reco_t_pt;
   Float_t         reco_t_eta;
   Float_t         reco_t_phi;
   Float_t         reco_t_m;
   Float_t         reco_tbar_pt;
   Float_t         reco_tbar_eta;
   Float_t         reco_tbar_phi;
   Float_t         reco_tbar_m;
   Float_t         reco_ttbar_pt;
   Float_t         reco_ttbar_eta;
   Float_t         reco_ttbar_phi;
   Float_t         reco_ttbar_m;


   // List of branches
  // List of branches
   TBranch        *b_eventNumber;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_mcChannelNumber;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_mu_actual;   //!
   TBranch        *b_mu_original_xAOD;   //!
   TBranch        *b_mu_actual_original_xAOD;   //!
   TBranch        *b_backgroundFlags;   //!
   TBranch        *b_hasBadMuon;   //!
   TBranch        *b_allhad_2017;   //!
   TBranch        *b_ejets_2017;   //!
   TBranch        *b_mujets_2017;   //!
   TBranch        *b_ee_2017;   //!
   TBranch        *b_mumu_2017;   //!
   TBranch        *b_emu_2017;   //!
   TBranch        *b_HLT_mu50;   //!
   TBranch        *b_HLT_mu26_ivarmedium;   //!
   TBranch        *b_HLT_mu15_L1MU10;   //!
   TBranch        *b_HLT_e26_lhtight_nod0_ivarloose;   //!
   TBranch        *b_HLT_mu14;   //!
   TBranch        *b_HLT_e17_lhloose_nod0;   //!
   TBranch        *b_HLT_2mu4;   //!
   TBranch        *b_HLT_j175;   //!
   TBranch        *b_HLT_j110;   //!
   TBranch        *b_HLT_e140_lhloose_nod0;   //!
   TBranch        *b_HLT_j260;   //!
   TBranch        *b_HLT_j85;   //!
   TBranch        *b_HLT_e60_lhmedium_nod0;   //!
   TBranch        *b_HLT_e15_lhloose_nod0_L1EM12;   //!
   TBranch        *b_HLT_e14_lhtight_nod0;   //!
   TBranch        *b_reco_fakeEvent;   //!
   TBranch        *b_reco_t_fCoordinates_fPt;   //!
   TBranch        *b_reco_t_fCoordinates_fEta;   //!
   TBranch        *b_reco_t_fCoordinates_fPhi;   //!
   TBranch        *b_reco_t_fCoordinates_fM;   //!
   TBranch        *b_reco_tbar_fCoordinates_fPt;   //!
   TBranch        *b_reco_tbar_fCoordinates_fEta;   //!
   TBranch        *b_reco_tbar_fCoordinates_fPhi;   //!
   TBranch        *b_reco_tbar_fCoordinates_fM;   //!
   TBranch        *b_reco_t_mass;   //!
   TBranch        *b_reco_ttbar_fCoordinates_fPt;   //!
   TBranch        *b_reco_ttbar_fCoordinates_fEta;   //!
   TBranch        *b_reco_ttbar_fCoordinates_fPhi;   //!
   TBranch        *b_reco_ttbar_fCoordinates_fM;   //!
   TBranch        *b_reco_lep_pt;   //!
   TBranch        *b_reco_lep_eta;   //!
   TBranch        *b_reco_lep_phi;   //!
   TBranch        *b_reco_lep_e;   //!
   TBranch        *b_reco_lep_pdgid;   //!
   TBranch        *b_reco_dilep_m;   //!
   TBranch        *b_reco_b_pt;   //!
   TBranch        *b_reco_b_m;   //!
   TBranch        *b_reco_b_phi;   //!
   TBranch        *b_reco_b_eta;   //!
   TBranch        *b_reco_bbar_pt;   //!
   TBranch        *b_reco_bbar_m;   //!
   TBranch        *b_reco_bbar_phi;   //!
   TBranch        *b_reco_bbar_eta;   //!
   TBranch        *b_reco_jet1_pt;   //!
   TBranch        *b_reco_jet1_eta;   //!
   TBranch        *b_reco_jet1_phi;   //!
   TBranch        *b_reco_jet1_e;   //!
   TBranch        *b_reco_jet2_pt;   //!
   TBranch        *b_reco_jet2_eta;   //!
   TBranch        *b_reco_jet2_phi;   //!
   TBranch        *b_reco_jet2_e;   //!
   TBranch        *b_reco_met_ex;   //!
   TBranch        *b_reco_met_ey;   //!
   TBranch        *b_reco_nb;   //!
   TBranch        *b_reco_nbbar;   //!
   TBranch        *b_reco_njets;   //!
   TBranch        *b_reco_proton_;   //!
   TBranch        *b_reco_proton_fCoordinates_fPt;   //!
   TBranch        *b_reco_proton_fCoordinates_fEta;   //!
   TBranch        *b_reco_proton_fCoordinates_fPhi;   //!
   TBranch        *b_reco_proton_fCoordinates_fE;   //!
   TBranch        *b_reco_proton_chi2;   //!
   TBranch        *b_reco_proton_side;   //!
   TBranch        *b_reco_proton_methodID;   //!
   TBranch        *b_reco_proton_xi;   //!
   TBranch        *b_reco_proton_ntracks;   //!
   TBranch        *b_reco_proton_stations;   //!
   TBranch        *b_reco_tau_event;   //!
   TBranch        *b_reco_top_reco_method;   //!
   TBranch        *b_reco_t_pt;   //!
   TBranch        *b_reco_t_eta;   //!
   TBranch        *b_reco_t_phi;   //!
   TBranch        *b_reco_t_m;   //!
   TBranch        *b_reco_tbar_pt;   //!
   TBranch        *b_reco_tbar_eta;   //!
   TBranch        *b_reco_tbar_phi;   //!
   TBranch        *b_reco_tbar_m;   //!
   TBranch        *b_reco_ttbar_pt;   //!
   TBranch        *b_reco_ttbar_eta;   //!
   TBranch        *b_reco_ttbar_phi;   //!
   TBranch        *b_reco_ttbar_m;   //!

//top_analysis_run2(TTree *tree=0) : fChain(0){};
   virtual ~top_analysis_run2() {};
   //virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    Version() const { return 2; }
   virtual void     Begin(TTree *tree);
   virtual void     SlaveBegin(TTree *tree);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();

   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }

     // Get Output List t osave our histograms in the output file

   virtual TList  *GetOutputList() const { return fOutput; }
   //
   virtual void    define_histograms();
     //
   virtual void    FillOutputList();
     //
   virtual void    WriteHistograms();

   virtual void    PrintHistograms();

   virtual void    SlaveTerminate();

   virtual void    Terminate();

   int nEvents;

   ClassDef(top_analysis_run2,0);
};

#endif
#ifdef top_analysis_run2_cxx




void top_analysis_run2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
    // Set object pointer
   reco_lep_pt = 0;
   reco_lep_eta = 0;
   reco_lep_phi = 0;
   reco_lep_e = 0;
   reco_lep_pdgid = 0;
   reco_proton_chi2 = 0;
   reco_proton_side = 0;
   reco_proton_methodID = 0;
   reco_proton_xi = 0;
   reco_proton_ntracks = 0;
   reco_proton_stations = 0;
   // Set branch addresses and branch pointers
   
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

  
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("mu_actual", &mu_actual, &b_mu_actual);
   fChain->SetBranchAddress("mu_original_xAOD", &mu_original_xAOD, &b_mu_original_xAOD);
   fChain->SetBranchAddress("mu_actual_original_xAOD", &mu_actual_original_xAOD, &b_mu_actual_original_xAOD);
   fChain->SetBranchAddress("backgroundFlags", &backgroundFlags, &b_backgroundFlags);
   fChain->SetBranchAddress("hasBadMuon", &hasBadMuon, &b_hasBadMuon);
   fChain->SetBranchAddress("allhad_2017", &allhad_2017, &b_allhad_2017);
   fChain->SetBranchAddress("ejets_2017", &ejets_2017, &b_ejets_2017);
   fChain->SetBranchAddress("mujets_2017", &mujets_2017, &b_mujets_2017);
   fChain->SetBranchAddress("ee_2017", &ee_2017, &b_ee_2017);
   fChain->SetBranchAddress("mumu_2017", &mumu_2017, &b_mumu_2017);
   fChain->SetBranchAddress("emu_2017", &emu_2017, &b_emu_2017);
   fChain->SetBranchAddress("HLT_mu50", &HLT_mu50, &b_HLT_mu50);
   fChain->SetBranchAddress("HLT_mu26_ivarmedium", &HLT_mu26_ivarmedium, &b_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("HLT_mu15_L1MU10", &HLT_mu15_L1MU10, &b_HLT_mu15_L1MU10);
   fChain->SetBranchAddress("HLT_e26_lhtight_nod0_ivarloose", &HLT_e26_lhtight_nod0_ivarloose, &b_HLT_e26_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("HLT_mu14", &HLT_mu14, &b_HLT_mu14);
   fChain->SetBranchAddress("HLT_e17_lhloose_nod0", &HLT_e17_lhloose_nod0, &b_HLT_e17_lhloose_nod0);
   fChain->SetBranchAddress("HLT_2mu4", &HLT_2mu4, &b_HLT_2mu4);
   fChain->SetBranchAddress("HLT_j175", &HLT_j175, &b_HLT_j175);
   fChain->SetBranchAddress("HLT_j110", &HLT_j110, &b_HLT_j110);
   fChain->SetBranchAddress("HLT_e140_lhloose_nod0", &HLT_e140_lhloose_nod0, &b_HLT_e140_lhloose_nod0);
   fChain->SetBranchAddress("HLT_j260", &HLT_j260, &b_HLT_j260);
   fChain->SetBranchAddress("HLT_j85", &HLT_j85, &b_HLT_j85);
   fChain->SetBranchAddress("HLT_e60_lhmedium_nod0", &HLT_e60_lhmedium_nod0, &b_HLT_e60_lhmedium_nod0);
   fChain->SetBranchAddress("HLT_e15_lhloose_nod0_L1EM12", &HLT_e15_lhloose_nod0_L1EM12, &b_HLT_e15_lhloose_nod0_L1EM12);
   fChain->SetBranchAddress("HLT_e14_lhtight_nod0", &HLT_e14_lhtight_nod0, &b_HLT_e14_lhtight_nod0);
   fChain->SetBranchAddress("reco_fakeEvent", &reco_fakeEvent, &b_reco_fakeEvent);
   

   fChain->SetBranchAddress("reco_lep_pt", &reco_lep_pt, &b_reco_lep_pt);
   fChain->SetBranchAddress("reco_lep_eta", &reco_lep_eta, &b_reco_lep_eta);
   fChain->SetBranchAddress("reco_lep_phi", &reco_lep_phi, &b_reco_lep_phi);
   fChain->SetBranchAddress("reco_lep_e", &reco_lep_e, &b_reco_lep_e);
   fChain->SetBranchAddress("reco_lep_pdgid", &reco_lep_pdgid, &b_reco_lep_pdgid);
   fChain->SetBranchAddress("reco_dilep_m", &reco_dilep_m, &b_reco_dilep_m);
   fChain->SetBranchAddress("reco_b_pt", &reco_b_pt, &b_reco_b_pt);
   fChain->SetBranchAddress("reco_b_m", &reco_b_m, &b_reco_b_m);
   fChain->SetBranchAddress("reco_b_phi", &reco_b_phi, &b_reco_b_phi);
   fChain->SetBranchAddress("reco_b_eta", &reco_b_eta, &b_reco_b_eta);
   fChain->SetBranchAddress("reco_bbar_pt", &reco_bbar_pt, &b_reco_bbar_pt);
   fChain->SetBranchAddress("reco_bbar_m", &reco_bbar_m, &b_reco_bbar_m);
   fChain->SetBranchAddress("reco_bbar_phi", &reco_bbar_phi, &b_reco_bbar_phi);
   fChain->SetBranchAddress("reco_bbar_eta", &reco_bbar_eta, &b_reco_bbar_eta);
   fChain->SetBranchAddress("reco_jet1_pt", &reco_jet1_pt, &b_reco_jet1_pt);
   fChain->SetBranchAddress("reco_jet1_eta", &reco_jet1_eta, &b_reco_jet1_eta);
   fChain->SetBranchAddress("reco_jet1_phi", &reco_jet1_phi, &b_reco_jet1_phi);
   fChain->SetBranchAddress("reco_jet1_e", &reco_jet1_e, &b_reco_jet1_e);
   fChain->SetBranchAddress("reco_jet2_pt", &reco_jet2_pt, &b_reco_jet2_pt);
   fChain->SetBranchAddress("reco_jet2_eta", &reco_jet2_eta, &b_reco_jet2_eta);
   fChain->SetBranchAddress("reco_jet2_phi", &reco_jet2_phi, &b_reco_jet2_phi);
   fChain->SetBranchAddress("reco_jet2_e", &reco_jet2_e, &b_reco_jet2_e);
   fChain->SetBranchAddress("reco_met_ex", &reco_met_ex, &b_reco_met_ex);
   fChain->SetBranchAddress("reco_met_ey", &reco_met_ey, &b_reco_met_ey);
   fChain->SetBranchAddress("reco_nb", &reco_nb, &b_reco_nb);
   fChain->SetBranchAddress("reco_nbbar", &reco_nbbar, &b_reco_nbbar);
   fChain->SetBranchAddress("reco_njets", &reco_njets, &b_reco_njets);
   fChain->SetBranchAddress("reco_proton", &reco_proton_, &b_reco_proton_);
   fChain->SetBranchAddress("reco_proton.fCoordinates.fPt", reco_proton_fCoordinates_fPt, &b_reco_proton_fCoordinates_fPt);
   fChain->SetBranchAddress("reco_proton.fCoordinates.fEta", reco_proton_fCoordinates_fEta, &b_reco_proton_fCoordinates_fEta);
   fChain->SetBranchAddress("reco_proton.fCoordinates.fPhi", reco_proton_fCoordinates_fPhi, &b_reco_proton_fCoordinates_fPhi);
   fChain->SetBranchAddress("reco_proton.fCoordinates.fE", reco_proton_fCoordinates_fE, &b_reco_proton_fCoordinates_fE);
   fChain->SetBranchAddress("reco_proton_chi2", &reco_proton_chi2, &b_reco_proton_chi2);
   fChain->SetBranchAddress("reco_proton_side", &reco_proton_side, &b_reco_proton_side);
   fChain->SetBranchAddress("reco_proton_methodID", &reco_proton_methodID, &b_reco_proton_methodID);
   fChain->SetBranchAddress("reco_proton_xi", &reco_proton_xi, &b_reco_proton_xi);
   fChain->SetBranchAddress("reco_proton_ntracks", &reco_proton_ntracks, &b_reco_proton_ntracks);
   fChain->SetBranchAddress("reco_proton_stations", &reco_proton_stations, &b_reco_proton_stations);
   fChain->SetBranchAddress("reco_tau_event", &reco_tau_event, &b_reco_tau_event);
   fChain->SetBranchAddress("reco_top_reco_method", &reco_top_reco_method, &b_reco_top_reco_method);
   fChain->SetBranchAddress("reco_t_pt", &reco_t_pt, &b_reco_t_pt);
   fChain->SetBranchAddress("reco_t_eta", &reco_t_eta, &b_reco_t_eta);
   fChain->SetBranchAddress("reco_t_phi", &reco_t_phi, &b_reco_t_phi);
   fChain->SetBranchAddress("reco_t_m", &reco_t_m, &b_reco_t_m);
   fChain->SetBranchAddress("reco_tbar_pt", &reco_tbar_pt, &b_reco_tbar_pt);
   fChain->SetBranchAddress("reco_tbar_eta", &reco_tbar_eta, &b_reco_tbar_eta);
   fChain->SetBranchAddress("reco_tbar_phi", &reco_tbar_phi, &b_reco_tbar_phi);
   fChain->SetBranchAddress("reco_tbar_m", &reco_tbar_m, &b_reco_tbar_m);
   fChain->SetBranchAddress("reco_ttbar_pt", &reco_ttbar_pt, &b_reco_ttbar_pt);
   fChain->SetBranchAddress("reco_ttbar_eta", &reco_ttbar_eta, &b_reco_ttbar_eta);
   fChain->SetBranchAddress("reco_ttbar_phi", &reco_ttbar_phi, &b_reco_ttbar_phi);
   fChain->SetBranchAddress("reco_ttbar_m", &reco_ttbar_m, &b_reco_ttbar_m);
   
   Notify();
}

Bool_t top_analysis_run2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef top_analysis_run2_cxx
