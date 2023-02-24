#define top_analysis_run2_cxx
#include "top_analysis_run2.h"
#include "top_analysis_run2_histograms.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TROOT.h"
#include "TKey.h"
#include "TList.h"
#include "TChain.h"
#include "TFile.h"
#include "TProof.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <sstream>
#include <cstdio>
#include <iomanip>


using namespace std;
string name;

void top_analysis_run2::Begin(TTree * )
{
    nEvents=0;
}

void top_analysis_run2::SlaveBegin(TTree * )
{

    TString option = GetOption();
    printf("Starting analysis with process option: %s \n", option.Data());

    name=option;

    define_histograms();

    FillOutputList();
}

vector<float> lepton;  


Bool_t top_analysis_run2::Process(Long64_t entry)
{
    fChain->GetTree()->GetEntry(entry);
    nEvents++; 

  
  int side = -1;  
  double side1 = 0;
  double side2 = 0;
  double reco_proton_methodID_total = 0;
  double reco_proton_ntracks_total = 0;
 


    if (nEvents % 10 == 0) std::cout << "Analysed a total of: " << nEvents << " events out of " << fChain->GetTree()->GetEntries() << " in this sample" << std::endl;
    if(fChain->GetTree()->GetEntries() > 0) {

/*
    reco_lep_pt->resize(0);
    reco_lep_eta->resize(0);
    reco_lep_phi->resize(0);
    reco_lep_e->resize(0);
    reco_lep_pdgid->resize(0);    
  */
   //TLorentzVector ejets_DL1r_2016_p4 = TLorentzVector();
    //Control plots
    h_t_mass->Fill(reco_t_m);
    h_t_pt->Fill(reco_t_pt);
    h_t_eta->Fill(reco_t_eta);
    h_t_phi->Fill(reco_t_phi);

    h_tbar_mass->Fill(reco_tbar_m);
    h_tbar_pt->Fill(reco_tbar_pt);
    h_tbar_eta->Fill(reco_tbar_eta);
    h_tbar_phi->Fill(reco_tbar_phi);

    h_ttbar_mass->Fill(reco_ttbar_m);
    h_ttbar_pt->Fill(reco_ttbar_pt);
    h_ttbar_eta->Fill(reco_ttbar_eta);
    h_ttbar_phi->Fill(reco_ttbar_phi);

    TLorentzVector j1; 
    TLorentzVector j2;
    TLorentzVector j1j2;
    TLorentzVector lep_neg;
    TLorentzVector lep_pos;
    TLorentzVector dilep;

    j1.SetPtEtaPhiE(reco_jet1_pt,reco_jet1_eta,reco_jet1_phi,reco_jet1_e);    
    j2.SetPtEtaPhiE(reco_jet2_pt,reco_jet2_eta,reco_jet2_phi,reco_jet2_e);    
    j1j2 = j1 + j2;

   //Choosing the best combination for leptons  
        for(uint i=0; i< reco_lep_pdgid->size();i++){
       // std::cout<<i<<std::endl;
        if( reco_lep_pdgid->size() == 2 ){
           if( reco_lep_pdgid->at(0) > 0 ){
               lep_neg.SetPtEtaPhiE(reco_lep_pt->at(0),reco_lep_eta->at(0),reco_lep_phi->at(0),reco_lep_e->at(0));
               lep_pos.SetPtEtaPhiE(reco_lep_pt->at(1),reco_lep_eta->at(1),reco_lep_phi->at(1),reco_lep_e->at(1));
            } else{
               lep_neg.SetPtEtaPhiE(reco_lep_pt->at(1),reco_lep_eta->at(1),reco_lep_phi->at(1),reco_lep_e->at(1));
               lep_pos.SetPtEtaPhiE(reco_lep_pt->at(0),reco_lep_eta->at(0),reco_lep_phi->at(0),reco_lep_e->at(0));
               }
            }
         }

   dilep = lep_pos + lep_neg;
    
  
  for(unsigned int i = 0; i < reco_proton_side->size(); i++){  
      side++;
 //if(side == 0 || side ==1 )      std::cout << side << std::endl;
      if(side == 0) side1 = reco_proton_xi->at(i);
      if(side == 1) side2 = reco_proton_xi->at(i);
      reco_proton_methodID_total = reco_proton_methodID->at(i);
      reco_proton_ntracks_total = reco_proton_methodID->at(i);

  }

   h_proton_side1_xi->Fill(side1);
   h_proton_side2_xi->Fill(side2);

   
   
   
   h_lep_pos_pt->Fill(lep_pos.Pt());
   h_lep_neg_pt->Fill(lep_neg.Pt());
   h_dilep_m->Fill(dilep.M());
   h_j1_pt->Fill(j1.Pt());
   h_j2_pt->Fill(j2.Pt());
   h_j1j2_pt->Fill(j1j2.Pt());
   h_j1j2_m->Fill(j1j2.M());

        /*
for(int i = 0; i < ejets_DL1r_2016; i++){ 

   n_ejets_DL1r_2016.push_back(i);

   ejets_DL1r_2016_p4.SetPtEtaPhiE(jet_pt->at(i), jet_eta->at(i), jet_phi->at(i), jet_e->at(i));
   std:cout << "Pt jets: " << ejets_DL1r_2016_p4.Pt() << " njets: " << n_ejets_DL1r_2016.size() << std::endl;
}
*/

}
return kTRUE; 
}

void top_analysis_run2::SlaveTerminate()
{

}

void top_analysis_run2::Terminate()
{


TString path = "/Users/danielernani/Desktop/workspace/MyTopAnalysis/";

TString filename_option = GetOption();
printf("Writting with name option: %s \n", filename_option.Data());
TString output_name = path+"output/"+filename_option+".root";
const char* filename = output_name;
TFile physicsoutput_TTbar(filename,"recreate");
WriteHistograms();
PrintHistograms();
physicsoutput_TTbar.Close();
}
