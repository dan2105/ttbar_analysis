#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <iostream>
#include "TCanvas.h"
#include <TStyle.h>


void top_analysis_run2::define_histograms()
{
    h_t_mass = new TH1F("top_mass","Reco Top mass; M_{t}[GeV]; Events/25 GeV",500,100,300);
    h_t_pt   = new TH1F("top_pt",  "Reco Top pT; pT_{t}[GeV]; Events/25 GeV",20,0,500);
    h_t_eta  = new TH1F("top_eta", "Reco Top eta; #eta_{t}; Events",100,-6,6);
    h_t_phi  = new TH1F("top_phi", "Reco Top phi; #phi_{t}; Events",100,-3.14,3.14);

    h_tbar_mass = new TH1F("tbar_mass", "Reco Top bar mass; M_{#bar{t}}[GeV]; Events/Bin",500,100,300);
    h_tbar_pt   = new TH1F("tbar_pt",   "Reco Top bar pT; pT_{#bar{t}}[[GeV]; Events/Bin",20,0,500);
    h_tbar_eta  = new TH1F("tbar_eta", "Reco Top bar eta; #eta_{#bar{t}}; Events/Bin",100,-6,6);
    h_tbar_phi  = new TH1F("tbar_phi", "Reco Top bar phi; #phi_{#bar{t}}}; Events/Bin",100,-3.14,3.14);

    h_ttbar_mass = new TH1F("ttbar_mass", "Reco tt bar mass; M_{#bar{t}}[GeV]; Events/Bin",40,0,1000);
    h_ttbar_pt   = new TH1F("ttbar_pt",   "Reco tt bar pT; pT_{#bar{t}}[[GeV]; Events/Bin",40,0,1000);
    h_ttbar_eta  = new TH1F("ttbar_eta", "Reco tt bar eta; #eta_{#bar{t}}; Events/Bin",100,-6,6);
    h_ttbar_phi  = new TH1F("ttbar_phi", "Reco tt bar phi; #phi_{#bar{t}}}; Events/Bin",100,-3.14,3.14);

    h_j1_pt = new TH1F("j1_pt","Reco j1 pT; pT(j_{1})[[GeV]; Events/25 GeV",40,0,1000);
    h_j2_pt = new TH1F("j2_pt","Reco j2 pT; pT(j_{2})[[GeV]; Events/25 GeV",40,0,1000);
    h_j1j2_pt = new TH1F("j1j2_pt","Reco j1j2 pT; pT(j_{1}j_{2})[[GeV]; Events/ GeV",80,25,665);
    h_j1j2_m = new TH1F("j1j2_m","Reco j1j2 M; M(j_{1}j_{2})[[GeV]; Events/ GeV",80,25,665);
    
    h_lep_pos_pt = new TH1F("lep_pos_pt","Reco lep_pos_pt M; pT(lep^{+})[GeV]; Events/ GeV",25,0,500);
    h_lep_neg_pt = new TH1F("lep_neg_pt","Reco lep_pos_pt M; pT(lep^{-})[GeV]; Events/ GeV",25,0,500);
    h_dilep_m = new TH1F("lep1lep2_m","Reco lep1lep2 M; M(lep_{1}lep_{2})[[GeV]; Events/ GeV",100,20,120);

    h_proton_side1_xi =new TH1F("proton_side1","Reco proton xi side 1; #xi_{side 1}; Events",100,0,0.2);
    h_proton_side2_xi =new TH1F("proton_side2","Reco proton xi side 2; #xi_{side 2}; Events",100,0,0.2);

}
///////////////////////////////////////////////////////////////////////////////
void top_analysis_run2::FillOutputList()
{

    GetOutputList()->Add(h_t_mass);
    GetOutputList()->Add(h_t_pt);
    GetOutputList()->Add(h_t_eta);
    GetOutputList()->Add(h_t_phi);

    GetOutputList()->Add(h_tbar_mass);
    GetOutputList()->Add(h_tbar_pt);
    GetOutputList()->Add(h_tbar_eta);
    GetOutputList()->Add(h_tbar_phi);

    GetOutputList()->Add(h_ttbar_mass);
    GetOutputList()->Add(h_ttbar_pt);
    GetOutputList()->Add(h_ttbar_eta);
    GetOutputList()->Add(h_ttbar_phi);

    GetOutputList()->Add(h_j1_pt);
    GetOutputList()->Add(h_j2_pt);
    GetOutputList()->Add(h_j1j2_pt);
    GetOutputList()->Add(h_j1j2_m);

    GetOutputList()->Add(h_dilep_m);
    GetOutputList()->Add(h_lep_pos_pt);
    GetOutputList()->Add(h_lep_neg_pt);

    GetOutputList()->Add(h_proton_side1_xi);
    GetOutputList()->Add(h_proton_side2_xi);


}

////////////////////////////////////////////////////////////////////////////////
void top_analysis_run2::WriteHistograms()
{

    h_t_mass->Write();
    h_t_pt->Write();
    h_t_eta->Write();
    h_t_phi->Write();

    h_tbar_mass->Write();
    h_tbar_pt->Write();
    h_tbar_eta->Write();
    h_tbar_phi->Write();

    h_ttbar_mass->Write();
    h_ttbar_pt->Write();
    h_ttbar_eta->Write();
    h_ttbar_phi->Write();

    h_j1_pt->Write();
    h_j2_pt->Write();
    h_j1j2_pt->Write();
    h_j1j2_m->Write();

    h_dilep_m->Write();
    h_lep_pos_pt->Write();
    h_lep_neg_pt->Write();

   h_proton_side1_xi->Write();
   h_proton_side2_xi->Write();
    
}
void top_analysis_run2::PrintHistograms(){

TString path = "/Users/danielernani/Desktop/workspace/MyTopAnalysis";
gStyle->SetOptStat(0);



}
