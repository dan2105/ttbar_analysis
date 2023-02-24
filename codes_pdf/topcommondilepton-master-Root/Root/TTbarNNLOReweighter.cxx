/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/
#include <TMath.h>
#include "TopCommonDilepton/TTbarNNLOReweighter.h"
#include <iostream>
#include "PathResolver/PathResolver.h"
using namespace std;

TTbarNNLOReweighter::TTbarNNLOReweighter(int sampleID, std::string dirName) :
  m_sampleID(sampleID),
  m_weightsFileName(""),
  m_Weights_file(0),
  m_Hist_TopPt_Powheg_Pythia8_Nominal(0),
  m_Hist_TTbarM_Powheg_Pythia8_Nominal(0),
  m_Hist_TTbarY_Powheg_Pythia8_Nominal(0),
  m_Hist_TopY_Powheg_Pythia8_Nominal(0),
  m_Hist_TopPt_PDFvar_Powheg_Pythia8_Nominal(0),
  m_Hist_TopPt_ScaleMin_Powheg_Pythia8_Nominal(0),
  m_Hist_TopPt_ScaleMax_Powheg_Pythia8_Nominal(0),
  m_Hist_TopPt_aMCNLO_Pythia8(0),
  m_Hist_TopPt_Powheg_Herwig7(0),
  m_Hist_TopPt_Powheg_Pythia8_Hdamp(0),
  m_Hist_TopPt_Sherpa(0),
  m_Hist_dPhi_ll(0),
  m_Hist_TopPt_411049_411057_mass173GeV(0),
  m_Hist_TopPt_411046_411054_mass172GeV(0),
  m_Hist_TopPt_Powheg_Herwig713(0),
  m_Hist_TopPt_aMCNLO_Herwig713(0)
{
  
  if(m_weightsDirName==""){
    m_weightsDirName = PathResolverFindCalibDirectory("TopCommonDilepton/");
  }
  else
    m_weightsDirName = dirName;
  
}

TTbarNNLOReweighter::~TTbarNNLOReweighter(){
  if(m_Weights_file){
    m_Weights_file->Close();
    delete m_Weights_file;
  }
}

void TTbarNNLOReweighter::SetInputDirectory(std::string dirName){
  m_weightsDirName = dirName;
}

void TTbarNNLOReweighter::SetInputFile(std::string fileName){
  m_weightsFileName = fileName;
}

void TTbarNNLOReweighter::SetSampleID(int sampleID){
  m_sampleID = sampleID;
}

void TTbarNNLOReweighter::Init(){
  
  if(m_weightsFileName==""){
    if(m_sampleID!=0){
      m_weightsFileName = m_weightsDirName+"/TTbar_NNLO_QCD_NLO_EW_weights.root";
    }
    else{
      std::cout << "TTbarNNLOReweighter::WARNING: No input file nor valid sampleID specified. Not able to get reweighting." << std::endl;
      return;
    }
  }
  
  m_Weights_file = TFile::Open(m_weightsFileName.c_str());
  if(m_Weights_file==0x0){
    std::cout << "TTbarNNLOReweighter::ERROR: Not able to open weight file " << m_weightsFileName << std::endl;
    return;
  }
  
  
  m_Hist_TopPt_Sherpa = (TH1*)m_Weights_file->Get("TopPt_410249_410250_410251_410252");
  m_Hist_TopPt_Powheg_Herwig7 = (TH1*)m_Weights_file->Get("TopPt_410557_410558_410559");
  m_Hist_TopPt_Powheg_Pythia8_Hdamp = (TH1*)m_Weights_file->Get("TopPt_410480_410481_410482");
  m_Hist_TopPt_aMCNLO_Pythia8 = (TH1*)m_Weights_file->Get("TopPt_410464_410465_410466");
  
  if (m_Hist_TopPt_Sherpa) m_Hist_TopPt_Sherpa   ->SetDirectory(0);
  if (m_Hist_TopPt_Powheg_Herwig7) m_Hist_TopPt_Powheg_Herwig7   ->SetDirectory(0);
  if (m_Hist_TopPt_Powheg_Pythia8_Hdamp) m_Hist_TopPt_Powheg_Pythia8_Hdamp  ->SetDirectory(0);
  if (m_Hist_TopPt_aMCNLO_Pythia8) m_Hist_TopPt_aMCNLO_Pythia8   ->SetDirectory(0);
  
  m_Hist_TopPt_Powheg_Pythia8_Nominal = (TH1*)m_Weights_file->Get("TopPt_410470_410471");
  m_Hist_TTbarM_Powheg_Pythia8_Nominal = (TH1*)m_Weights_file->Get("TTbarM_410470_410471");
  m_Hist_TTbarY_Powheg_Pythia8_Nominal = (TH1*)m_Weights_file->Get("TTbarY_410470_410471");
  m_Hist_TopY_Powheg_Pythia8_Nominal = (TH1*)m_Weights_file->Get("TopY_410470_410471");
  m_Hist_TopPt_PDFvar_Powheg_Pythia8_Nominal = (TH1*)m_Weights_file->Get("TopPt_410470_410471_PDFvariation");
  m_Hist_TopPt_ScaleMin_Powheg_Pythia8_Nominal = (TH1*)m_Weights_file->Get("TopPt_410470_410471_scale_min");
  m_Hist_TopPt_ScaleMax_Powheg_Pythia8_Nominal = (TH1*)m_Weights_file->Get("TopPt_410470_410471_scale_max");
  
  
  if (m_Hist_TopPt_Powheg_Pythia8_Nominal)  m_Hist_TopPt_Powheg_Pythia8_Nominal   ->SetDirectory(0);
  if (m_Hist_TTbarM_Powheg_Pythia8_Nominal)  m_Hist_TTbarM_Powheg_Pythia8_Nominal   ->SetDirectory(0);
  if (m_Hist_TTbarY_Powheg_Pythia8_Nominal)  m_Hist_TTbarY_Powheg_Pythia8_Nominal   ->SetDirectory(0);
  if (m_Hist_TopY_Powheg_Pythia8_Nominal)  m_Hist_TopY_Powheg_Pythia8_Nominal   ->SetDirectory(0);
  if (m_Hist_TopPt_PDFvar_Powheg_Pythia8_Nominal)  m_Hist_TopPt_PDFvar_Powheg_Pythia8_Nominal   ->SetDirectory(0);
  if (m_Hist_TopPt_ScaleMin_Powheg_Pythia8_Nominal)  m_Hist_TopPt_ScaleMin_Powheg_Pythia8_Nominal   ->SetDirectory(0);
  if (m_Hist_TopPt_ScaleMax_Powheg_Pythia8_Nominal)  m_Hist_TopPt_ScaleMax_Powheg_Pythia8_Nominal   ->SetDirectory(0);
  
  m_Hist_dPhi_ll = (TH1*)m_Weights_file->Get("dPhi_ll_410472");
  if (m_Hist_dPhi_ll) m_Hist_dPhi_ll ->SetDirectory(0);
 
  m_Hist_TopPt_411049_411057_mass173GeV = (TH1*)m_Weights_file->Get("TopPt_411049_411057_mass173GeV");
  if (m_Hist_TopPt_411049_411057_mass173GeV) m_Hist_TopPt_411049_411057_mass173GeV ->SetDirectory(0);

  m_Hist_TopPt_411046_411054_mass172GeV = (TH1*)m_Weights_file->Get("TopPt_411046_411054_mass172GeV");
  if (m_Hist_TopPt_411046_411054_mass172GeV) m_Hist_TopPt_411046_411054_mass172GeV ->SetDirectory(0);

  m_Hist_TopPt_Powheg_Herwig713 = (TH1*)m_Weights_file->Get("TopPt_411233_411234");
  if (m_Hist_TopPt_Powheg_Herwig713) m_Hist_TopPt_Powheg_Herwig713 ->SetDirectory(0);
 
  m_Hist_TopPt_aMCNLO_Herwig713 = (TH1*)m_Weights_file->Get("TopPt_412116_412117");
  if (m_Hist_TopPt_aMCNLO_Herwig713) m_Hist_TopPt_aMCNLO_Herwig713 ->SetDirectory(0);


  m_Weights_file->Close();
  
  
}

float TTbarNNLOReweighter::GetTopPt_Sherpa(float truth_top_pt){
  if(m_Hist_TopPt_Sherpa==0x0) return 1.;
  return m_Hist_TopPt_Sherpa ->GetBinContent(m_Hist_TopPt_Sherpa ->FindBin( truth_top_pt ));
}

float TTbarNNLOReweighter::GetTopPt_Powheg_Pythia8_Hdamp(float truth_top_pt){
  if(m_Hist_TopPt_Powheg_Pythia8_Hdamp ==0x0) return 1.;
  return m_Hist_TopPt_Powheg_Pythia8_Hdamp  ->GetBinContent(m_Hist_TopPt_Powheg_Pythia8_Hdamp ->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_Powheg_Herwig7(float truth_top_pt){
  if(m_Hist_TopPt_Powheg_Herwig7 ==0x0) return 1.;
  return m_Hist_TopPt_Powheg_Herwig7  ->GetBinContent(m_Hist_TopPt_Powheg_Herwig7 ->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_aMCNLO_Pythia8(float truth_top_pt){
  if(m_Hist_TopPt_aMCNLO_Pythia8 ==0x0) return 1.;
  return m_Hist_TopPt_aMCNLO_Pythia8  ->GetBinContent(m_Hist_TopPt_aMCNLO_Pythia8 ->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_Powheg_Pythia8_Nominal(float truth_top_pt){
  if(m_Hist_TopPt_Powheg_Pythia8_Nominal ==0x0) return 1.;
  return m_Hist_TopPt_Powheg_Pythia8_Nominal  ->GetBinContent(m_Hist_TopPt_Powheg_Pythia8_Nominal->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_PDFvar_Powheg_Pythia8_Nominal(float truth_top_pt){
  if(m_Hist_TopPt_PDFvar_Powheg_Pythia8_Nominal ==0x0) return 1.;
  return m_Hist_TopPt_PDFvar_Powheg_Pythia8_Nominal  ->GetBinContent(m_Hist_TopPt_PDFvar_Powheg_Pythia8_Nominal->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_ScaleMin_Powheg_Pythia8_Nominal(float truth_top_pt){
  if(m_Hist_TopPt_ScaleMin_Powheg_Pythia8_Nominal ==0x0) return 1.;
  return m_Hist_TopPt_ScaleMin_Powheg_Pythia8_Nominal  ->GetBinContent(m_Hist_TopPt_ScaleMin_Powheg_Pythia8_Nominal->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_ScaleMax_Powheg_Pythia8_Nominal(float truth_top_pt){
  if(m_Hist_TopPt_ScaleMax_Powheg_Pythia8_Nominal ==0x0) return 1.;
  return m_Hist_TopPt_ScaleMax_Powheg_Pythia8_Nominal  ->GetBinContent(m_Hist_TopPt_ScaleMax_Powheg_Pythia8_Nominal->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTTbarM_Powheg_Pythia8_Nominal(float truth_TTbar_M){
  if(m_Hist_TTbarM_Powheg_Pythia8_Nominal ==0x0) return 1.;
  return m_Hist_TTbarM_Powheg_Pythia8_Nominal ->GetBinContent(m_Hist_TTbarM_Powheg_Pythia8_Nominal->FindBin( truth_TTbar_M) );
}

float TTbarNNLOReweighter::GetTTbarY_Powheg_Pythia8_Nominal(float truth_TTbar_Y){
  if(m_Hist_TTbarY_Powheg_Pythia8_Nominal ==0x0) return 1.;
  return m_Hist_TTbarY_Powheg_Pythia8_Nominal ->GetBinContent(m_Hist_TTbarY_Powheg_Pythia8_Nominal->FindBin( truth_TTbar_Y ) );
}

float TTbarNNLOReweighter::GetTopY_Powheg_Pythia8_Nominal(float truth_top_Y){
  if( m_Hist_TopY_Powheg_Pythia8_Nominal ==0x0) return 1.;
  return  m_Hist_TopY_Powheg_Pythia8_Nominal ->GetBinContent( m_Hist_TopY_Powheg_Pythia8_Nominal->FindBin( truth_top_Y ) );
}

float TTbarNNLOReweighter::GetDeltaPhi(float truth_dPhi_ll){
  if( m_Hist_dPhi_ll ==0x0) return 1.;
  return m_Hist_dPhi_ll ->GetBinContent( m_Hist_dPhi_ll ->FindBin( truth_dPhi_ll ) );
}

float TTbarNNLOReweighter::GetTopPt_Powheg_Pythia8_TopMass172GeV(float truth_top_pt){
  if(m_Hist_TopPt_411046_411054_mass172GeV ==0x0) return 1.;
  return m_Hist_TopPt_411046_411054_mass172GeV ->GetBinContent( m_Hist_TopPt_411046_411054_mass172GeV ->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_Powheg_Pythia8_TopMass173GeV(float truth_top_pt){
 if(m_Hist_TopPt_411049_411057_mass173GeV ==0x0) return 1.;
 return m_Hist_TopPt_411049_411057_mass173GeV ->GetBinContent( m_Hist_TopPt_411049_411057_mass173GeV ->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_Powheg_Herwig713(float truth_top_pt){
 if(m_Hist_TopPt_Powheg_Herwig713 ==0x0) return 1.;
 return m_Hist_TopPt_Powheg_Herwig713 ->GetBinContent( m_Hist_TopPt_Powheg_Herwig713 ->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_aMCNLO_Herwig713(float truth_top_pt){
 if(m_Hist_TopPt_aMCNLO_Herwig713 ==0x0) return 1.;
 return m_Hist_TopPt_aMCNLO_Herwig713 ->GetBinContent( m_Hist_TopPt_aMCNLO_Herwig713 ->FindBin( truth_top_pt ) );
}

float TTbarNNLOReweighter::GetTopPt_Powheg_Pythia8_Nominal_Average(float truth_top_pt, float truth_antitop_pt){
  if(m_Hist_TopPt_Powheg_Pythia8_Nominal ==0x0) return 1.;
  return TMath::Sqrt( (m_Hist_TopPt_Powheg_Pythia8_Nominal  ->GetBinContent(m_Hist_TopPt_Powheg_Pythia8_Nominal->FindBin( truth_top_pt ) ) ) * (m_Hist_TopPt_Powheg_Pythia8_Nominal  ->GetBinContent(m_Hist_TopPt_Powheg_Pythia8_Nominal->FindBin( truth_antitop_pt ) ) ) );
}


