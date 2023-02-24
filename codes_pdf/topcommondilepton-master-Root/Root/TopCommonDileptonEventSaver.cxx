#include "TopCommonDilepton/TopCommonDileptonEventSaver.h"
#include "TopEvent/Event.h"
#include "TopEvent/EventTools.h"
#include "TopEventSelectionTools/TreeManager.h"
#include "TopParticleLevel/ParticleLevelEvent.h"
#include "TopConfiguration/TopConfig.h"
#include "TopEventSelectionTools/PlotManager.h"

#include "PathResolver/PathResolver.h"

#include "TopConfiguration/ConfigurationSettings.h"
#include "TopCorrections/ScaleFactorRetriever.h"
#include "TopEventSelectionTools/NJetBtagSelector.h"

#include "TopEvent/TopEventMaker.h"
#include "TopEvent/Event.h"

#include "TH1.h"
#include "TH2.h"
#include "TopCommonDilepton/NeutrinoWeighter.h"
#include "TopCommonDilepton/TopDileptonReconstruction.h"

#include "xAODTruth/xAODTruthHelpers.h"

//namespace top {

///-- Always initialise primitive types in the constructor --///
TopCommonDileptonEventSaver::TopCommonDileptonEventSaver() : 

  is_MC(true)
{ 
  top::ConfigurationSettings* configSettings = top::ConfigurationSettings::get();
  //m_NW_reconstructed = false;
  m_doNeutrinoWeighter      = false;
  m_doEllipseMethod         = false;
  m_doSonnenchein           = false;
  m_doNeutrinoWeighter      = (configSettings->value("NeutrinoWeighter") == "True");
  m_doEllipseMethod         = (configSettings->value("EllipseMethod")    == "True");
  m_doSonnenchein           = (configSettings->value("Sonnenchein")      == "True");
  m_store_weights_reco      = false;
  m_store_weights_particle  = false;
  m_store_weights_parton    = false;
  m_store_weights_reco      = (configSettings->value("StoreWeightsReco")     == "True");
  m_store_weights_particle  = (configSettings->value("StoreWeightsParticle") == "True");
  m_store_weights_parton    = (configSettings->value("StoreWeightsParton")   == "True");
  m_skipNominal = (configSettings->value("SkipNominalTree") == "True"); //this is a custom variable for if we ever need to separate ntuples by systematics
  branchFilters().push_back(std::bind(&TopCommonDileptonEventSaver::isBranchStored, this, std::placeholders::_1, std::placeholders::_2));
  m_ttbarNNLO_Reweighter = new TTbarNNLORecursiveRew(0,0,false);
}


void TopCommonDileptonEventSaver::initialize(std::shared_ptr<top::TopConfig> config, TFile* file,const std::vector<std::string>& extraBranches) {

  m_config = config;
  EventSaverFlatNtuple::initialize(config, file, extraBranches);

  if (config->isMC()){
    is_MC = true;
  } else {
    is_MC = false;
  }

  if(is_MC){
    m_ttbarNNLO_Reweighter->SetInputDirectory("./reweightingData/");
    m_ttbarNNLO_Reweighter->Init();
    m_ttbarNNLO_Reweighter->SetSampleID(0);
    m_ttbarNNLO_Reweighter->SetSystVar(0);
    addCustomPartonLevelVariables();
    if ( topConfig()->doTopParticleLevel() ){
      addCustomParticleLevelVariables();
    }
  }
  addCustomRecoLevelVariables(extraBranches);

}


void TopCommonDileptonEventSaver::addCustomPartonLevelVariables(){

  truthTreeManager()->makeOutputVariable(m_parton_weight_nnlo,           "parton_weight_nnlo");
  truthTreeManager()->makeOutputVariable(m_parton_l_pt,                  "parton_l_pt");
  truthTreeManager()->makeOutputVariable(m_parton_l_eta,                 "parton_l_eta");
  truthTreeManager()->makeOutputVariable(m_parton_l_phi,                 "parton_l_phi");
  truthTreeManager()->makeOutputVariable(m_parton_l_m,                   "parton_l_m");
  truthTreeManager()->makeOutputVariable(m_parton_l_pdgid,               "parton_l_pdgid");
  truthTreeManager()->makeOutputVariable(m_parton_lbar_pt,               "parton_lbar_pt");
  truthTreeManager()->makeOutputVariable(m_parton_lbar_eta,              "parton_lbar_eta");
  truthTreeManager()->makeOutputVariable(m_parton_lbar_phi,              "parton_lbar_phi");
  truthTreeManager()->makeOutputVariable(m_parton_lbar_m,                "parton_lbar_m");
  truthTreeManager()->makeOutputVariable(m_parton_lbar_pdgid,            "parton_lbar_pdgid");
  truthTreeManager()->makeOutputVariable(m_parton_t_afterFSR_pt,         "parton_t_pt");
  truthTreeManager()->makeOutputVariable(m_parton_t_afterFSR_eta,        "parton_t_eta");
  truthTreeManager()->makeOutputVariable(m_parton_t_afterFSR_phi,        "parton_t_phi");
  truthTreeManager()->makeOutputVariable(m_parton_t_afterFSR_m,          "parton_t_m");
  truthTreeManager()->makeOutputVariable(m_parton_tbar_afterFSR_pt,      "parton_tbar_pt");
  truthTreeManager()->makeOutputVariable(m_parton_tbar_afterFSR_eta,     "parton_tbar_eta");
  truthTreeManager()->makeOutputVariable(m_parton_tbar_afterFSR_phi,     "parton_tbar_phi");
  truthTreeManager()->makeOutputVariable(m_parton_tbar_afterFSR_m,       "parton_tbar_m");
  truthTreeManager()->makeOutputVariable(m_parton_ttbar_pt,              "parton_ttbar_pt");
  truthTreeManager()->makeOutputVariable(m_parton_ttbar_eta,             "parton_ttbar_eta");
  truthTreeManager()->makeOutputVariable(m_parton_ttbar_phi,             "parton_ttbar_phi");
  truthTreeManager()->makeOutputVariable(m_parton_ttbar_m,               "parton_ttbar_m");
  truthTreeManager()->makeOutputVariable(m_parton_cosphi,                "parton_cosphi");
  truthTreeManager()->makeOutputVariable(m_parton_cos_kaxis_p,           "parton_cos_kaxis_p");
  truthTreeManager()->makeOutputVariable(m_parton_cos_kaxis_m,           "parton_cos_kaxis_m");  
  truthTreeManager()->makeOutputVariable(m_parton_cos_naxis_p,           "parton_cos_naxis_p");
  truthTreeManager()->makeOutputVariable(m_parton_cos_naxis_m,           "parton_cos_naxis_m");
  truthTreeManager()->makeOutputVariable(m_parton_cos_raxis_p,           "parton_cos_raxis_p");
  truthTreeManager()->makeOutputVariable(m_parton_cos_raxis_m,           "parton_cos_raxis_m");
  truthTreeManager()->makeOutputVariable(m_parton_is_tau,                "parton_tau_event");
  if (m_store_weights_parton){
  truthTreeManager()->makeOutputVariable(m_parton_generator_weights,     "parton_generator_weights");
  }
}

void TopCommonDileptonEventSaver::addCustomParticleLevelVariables(){


  particleLevelTreeManager()->makeOutputVariable(m_parton_weight_nnlo,          "particle_weight_nnlo");
  particleLevelTreeManager()->makeOutputVariable(m_particle_lep_pt,             "particle_lep_pt");
  particleLevelTreeManager()->makeOutputVariable(m_particle_lep_eta,            "particle_lep_eta");
  particleLevelTreeManager()->makeOutputVariable(m_particle_lep_phi,            "particle_lep_phi");
  particleLevelTreeManager()->makeOutputVariable(m_particle_lep_m,              "particle_lep_m");
  particleLevelTreeManager()->makeOutputVariable(m_particle_lep_pdgid,          "particle_lep_pdgid");
  particleLevelTreeManager()->makeOutputVariable(m_particle_met_ex,             "particle_met_ex");
  particleLevelTreeManager()->makeOutputVariable(m_particle_met_ey,             "particle_met_ey");
  particleLevelTreeManager()->makeOutputVariable(m_particle_met_met,            "particle_met_met");
  particleLevelTreeManager()->makeOutputVariable(m_particle_met_phi,            "particle_met_phi");
  particleLevelTreeManager()->makeOutputVariable(m_particle_t_pt,               "particle_t_pt");
  particleLevelTreeManager()->makeOutputVariable(m_particle_t_eta,              "particle_t_eta");
  particleLevelTreeManager()->makeOutputVariable(m_particle_t_phi,              "particle_t_phi");
  particleLevelTreeManager()->makeOutputVariable(m_particle_t_m,                "particle_t_m");
  particleLevelTreeManager()->makeOutputVariable(m_particle_tbar_pt,            "particle_tbar_pt");
  particleLevelTreeManager()->makeOutputVariable(m_particle_tbar_eta,           "particle_tbar_eta");
  particleLevelTreeManager()->makeOutputVariable(m_particle_tbar_phi,           "particle_tbar_phi");
  particleLevelTreeManager()->makeOutputVariable(m_particle_tbar_m,             "particle_tbar_m");
  particleLevelTreeManager()->makeOutputVariable(m_particle_ttbar_pt,           "particle_ttbar_pt");
  particleLevelTreeManager()->makeOutputVariable(m_particle_ttbar_eta,          "particle_ttbar_eta");
  particleLevelTreeManager()->makeOutputVariable(m_particle_ttbar_phi,          "particle_ttbar_phi");
  particleLevelTreeManager()->makeOutputVariable(m_particle_ttbar_m,            "particle_ttbar_m");
  particleLevelTreeManager()->makeOutputVariable(m_particle_jet_pt,             "particle_jet_pt");
  particleLevelTreeManager()->makeOutputVariable(m_particle_jet_eta,            "particle_jet_eta");
  particleLevelTreeManager()->makeOutputVariable(m_particle_jet_phi,            "particle_jet_phi");
  particleLevelTreeManager()->makeOutputVariable(m_particle_jet_m,              "particle_jet_m");
  particleLevelTreeManager()->makeOutputVariable(m_particle_jet_btagged,        "particle_jet_btagged");
  particleLevelTreeManager()->makeOutputVariable(m_particle_cosphi,             "particle_cosphi");
  particleLevelTreeManager()->makeOutputVariable(m_particle_cos_kaxis_p,        "particle_cos_kaxis_p");
  particleLevelTreeManager()->makeOutputVariable(m_particle_cos_kaxis_m,        "particle_cos_kaxis_m");  
  particleLevelTreeManager()->makeOutputVariable(m_particle_cos_naxis_p,        "particle_cos_naxis_p");
  particleLevelTreeManager()->makeOutputVariable(m_particle_cos_naxis_m,        "particle_cos_naxis_m");
  particleLevelTreeManager()->makeOutputVariable(m_particle_cos_raxis_p,        "particle_cos_raxis_p");
  particleLevelTreeManager()->makeOutputVariable(m_particle_cos_raxis_m,        "particle_cos_raxis_m");
  particleLevelTreeManager()->makeOutputVariable(m_parton_is_tau,               "particle_truth_tau_event");
  if (m_store_weights_particle){
  particleLevelTreeManager()->makeOutputVariable(m_particle_generator_weights,  "particle_generator_weights");
  }
}


void TopCommonDileptonEventSaver::addCustomRecoLevelVariables(const std::vector<std::string>& extraBranches){


  for (auto systematicTree : treeManagers()) {
    systematicTree->makeOutputVariable(m_passSelection,                 "reco_passSelection");
    systematicTree->makeOutputVariable(m_parton_weight_nnlo,            "reco_weight_nnlo");    
    systematicTree->makeOutputVariable(m_lep_pt,                        "reco_lep_pt");
    systematicTree->makeOutputVariable(m_lep_eta,                       "reco_lep_eta");
    systematicTree->makeOutputVariable(m_lep_phi,                       "reco_lep_phi");
    systematicTree->makeOutputVariable(m_lep_e,                         "reco_lep_e");
    systematicTree->makeOutputVariable(m_lep_pdgid,                     "reco_lep_pdgid");
    if(m_config->isMC() == true){
    systematicTree->makeOutputVariable(m_fakeEvent,                     "reco_fakeEvent");
    systematicTree->makeOutputVariable(m_lep_true_type,                 "reco_lep_true_type");
    systematicTree->makeOutputVariable(m_lep_true_origin,               "reco_lep_true_origin");
    systematicTree->makeOutputVariable(m_lep_true_isPrompt,             "reco_lep_true_isPrompt");
    systematicTree->makeOutputVariable(m_lep_true_isChargeFl,           "reco_lep_true_isChargeFl");
    systematicTree->makeOutputVariable(m_lep_true_pdgid,                "reco_lep_true_pdgid");
    }
    systematicTree->makeOutputVariable(m_dilep_m,                       "reco_dilep_m");
    systematicTree->makeOutputVariable(m_ljet_pt,                       "reco_ljet_pt");
    systematicTree->makeOutputVariable(m_ljet_eta,                      "reco_ljet_eta");
    systematicTree->makeOutputVariable(m_ljet_phi,                      "reco_ljet_phi");
    systematicTree->makeOutputVariable(m_ljet_e,                        "reco_ljet_e");
    systematicTree->makeOutputVariable(m_bjet_pt,                       "reco_bjet_pt");
    systematicTree->makeOutputVariable(m_bjet_eta,                      "reco_bjet_eta");
    systematicTree->makeOutputVariable(m_bjet_phi,                      "reco_bjet_phi");
    systematicTree->makeOutputVariable(m_bjet_e,                        "reco_bjet_e");
    systematicTree->makeOutputVariable(m_met_met,                       "reco_met_met");
    systematicTree->makeOutputVariable(m_met_phi,                       "reco_met_phi");
    if(m_doNeutrinoWeighter || m_doSonnenchein || m_doEllipseMethod){
    systematicTree->makeOutputVariable(m_t_pt,                          "reco_t_pt");
    systematicTree->makeOutputVariable(m_t_eta,                         "reco_t_eta");
    systematicTree->makeOutputVariable(m_t_phi,                         "reco_t_phi");
    systematicTree->makeOutputVariable(m_t_m,                           "reco_t_m");
    systematicTree->makeOutputVariable(m_tbar_pt,                       "reco_tbar_pt");
    systematicTree->makeOutputVariable(m_tbar_eta,                      "reco_tbar_eta");
    systematicTree->makeOutputVariable(m_tbar_phi,                      "reco_tbar_phi");
    systematicTree->makeOutputVariable(m_tbar_m,                        "reco_tbar_m");
    systematicTree->makeOutputVariable(m_ttbar_pt,                      "reco_ttbar_pt");
    systematicTree->makeOutputVariable(m_ttbar_eta,                     "reco_ttbar_eta");
    systematicTree->makeOutputVariable(m_ttbar_phi,                     "reco_ttbar_phi");
    systematicTree->makeOutputVariable(m_ttbar_m,                       "reco_ttbar_m");
    systematicTree->makeOutputVariable(m_top_reco_method,               "reco_top_reco_method");
    systematicTree->makeOutputVariable(m_cosphi,                        "reco_cosphi");
    systematicTree->makeOutputVariable(m_cos_kaxis_p,                   "reco_cos_kaxis_p");
    systematicTree->makeOutputVariable(m_cos_kaxis_m,                   "reco_cos_kaxis_m");
    systematicTree->makeOutputVariable(m_cos_naxis_p,                   "reco_cos_naxis_p");
    systematicTree->makeOutputVariable(m_cos_naxis_m,                   "reco_cos_naxis_m");
    systematicTree->makeOutputVariable(m_cos_raxis_p,                   "reco_cos_raxis_p");
    systematicTree->makeOutputVariable(m_cos_raxis_m,                   "reco_cos_raxis_m");
    }
    systematicTree->makeOutputVariable(m_parton_is_tau,                 "reco_truth_tau_event");
    if (m_store_weights_reco){
    systematicTree->makeOutputVariable(m_generator_weights,             "reco_generator_weights");
    }    
    /// These branches are things that we're using but not storing for now ///s
    //systematicTree->makeOutputVariable(m_lep_cl_eta,                      "reco_lep_cl_eta"); 
    //systematicTree->makeOutputVariable(m_jet_DL1r,                        "reco_jet_DL1r"); 
    //systematicTree->makeOutputVariable(m_jet_jvt,                         "reco_jet_jvt");    
    //systematicTree->makeOutputVariable(m_jet_tagWeightBin_DL1_Continuous, "reco_jet_tagWeightBin_DL1r_Continuous");
    //systematicTree->makeOutputVariable(m_met_ex,                          "reco_met_ex");
    //systematicTree->makeOutputVariable(m_met_ey,                          "reco_met_ey");
    //systematicTree->makeOutputVariable(m_met_sumet,                       "reco_met_sumet");

  } /// End of tree-manager loop

}


void TopCommonDileptonEventSaver::initializePartonLevelEvent(){

  m_parton_weight_nnlo       =   1.;
  m_parton_l_pt              = -99.;         
  m_parton_l_eta             = -99.;
  m_parton_l_phi             = -99.;
  m_parton_l_m               = -99.;
  m_parton_l_pdgid           = -99.;
  m_parton_lbar_pt           = -99.;   
  m_parton_lbar_eta          = -99.;
  m_parton_lbar_phi          = -99.;
  m_parton_lbar_m            = -99.;
  m_parton_lbar_pdgid        = -99.;
  m_parton_t_afterFSR_pt     = -99.;
  m_parton_t_afterFSR_eta    = -99.;
  m_parton_t_afterFSR_phi    = -99.;
  m_parton_t_afterFSR_m      = -99.;
  m_parton_tbar_afterFSR_pt  = -99.;
  m_parton_tbar_afterFSR_eta = -99.;
  m_parton_tbar_afterFSR_phi = -99.;
  m_parton_tbar_afterFSR_m   = -99.;
  m_parton_ttbar_pt          = -99.;
  m_parton_ttbar_eta         = -99.;
  m_parton_ttbar_phi         = -99.;
  m_parton_ttbar_m           = -99.;      
  m_parton_cosphi            = -99.;
  m_parton_cos_kaxis_p       = -99.;
  m_parton_cos_kaxis_m       = -99.;
  m_parton_cos_naxis_p       = -99.;
  m_parton_cos_naxis_m       = -99.;
  m_parton_cos_raxis_p       = -99.;
  m_parton_cos_raxis_m       = -99.;
  m_parton_is_tau            = false;
  m_parton_generator_weights.clear();

}


void TopCommonDileptonEventSaver::initializeParticleLevelBranches(){


  m_particle_lep_pt.clear();
  m_particle_lep_eta.clear();
  m_particle_lep_phi.clear();
  m_particle_lep_m.clear();
  m_particle_lep_pdgid.clear();  
  m_particle_jet_pt.clear();
  m_particle_jet_eta.clear();
  m_particle_jet_phi.clear();
  m_particle_jet_m.clear();
  m_particle_jet_btagged.clear();  
  m_particle_met_ex      = -99.;
  m_particle_met_ey      = -99.;
  m_particle_met_phi     = -99.;
  m_particle_met_met     = -99.;
  m_particle_t_pt        = -99.;
  m_particle_t_eta       = -99.;
  m_particle_t_phi       = -99.;
  m_particle_t_m         = -99.;
  m_particle_tbar_pt     = -99.;
  m_particle_tbar_eta    = -99.;
  m_particle_tbar_phi    = -99.;
  m_particle_tbar_m      = -99.;
  m_particle_ttbar_pt    = -99.;
  m_particle_ttbar_eta   = -99.;
  m_particle_ttbar_phi   = -99.;
  m_particle_ttbar_m     = -99.;
  m_particle_cosphi      = -99.;    
  m_particle_cos_kaxis_p = -99.;
  m_particle_cos_kaxis_m = -99.;
  m_particle_cos_naxis_p = -99.;
  m_particle_cos_naxis_m = -99.;
  m_particle_cos_raxis_p = -99.;
  m_particle_cos_raxis_m = -99.;
  m_particle_generator_weights.clear();
}


void TopCommonDileptonEventSaver::initializeRecoLevelBranches(){

  m_passSelection = false;
  m_lep_pt.clear();
  m_lep_eta.clear();
  m_lep_phi.clear();
  m_lep_e.clear();
  m_lep_pdgid.clear();
  m_fakeEvent = false;
  m_lep_true_type.clear();
  m_lep_true_origin.clear();
  m_lep_true_isPrompt.clear();
  m_lep_true_isChargeFl.clear();
  m_lep_true_pdgid.clear();
  m_dilep_m = -99;
  m_jet_isbtagged_DL1r_85.clear();  
  m_ljet_pt.clear();
  m_ljet_eta.clear();
  m_ljet_phi.clear();
  m_ljet_e.clear();
  m_bjet_pt.clear();
  m_bjet_eta.clear();
  m_bjet_phi.clear();
  m_bjet_e.clear();
  m_met_met         =  -1.;
  m_met_phi         = -99.;
  m_t_pt            = -99.;
  m_t_eta           = -99.;
  m_t_phi           = -99.;
  m_t_m             = -99.;
  m_tbar_pt         = -99.;
  m_tbar_eta        = -99.;
  m_tbar_phi        = -99.;
  m_tbar_m          = -99.;
  m_ttbar_pt        = -99.;
  m_ttbar_eta       = -99.;
  m_ttbar_phi       = -99.;
  m_ttbar_m         = -99.;
  m_top_reco_method =  -1;
  m_cosphi          = -99.;
  m_cos_kaxis_p     = -99.;
  m_cos_kaxis_m     = -99.;
  m_cos_naxis_p     = -99.;
  m_cos_naxis_m     = -99.;
  m_cos_raxis_p     = -99.;
  m_cos_raxis_m     = -99.;
  m_generator_weights.clear();   
}


void TopCommonDileptonEventSaver::saveTruthEvent(){

  // Reset everything
  initializePartonLevelEvent();

  if (!is_MC || !m_config->doTopPartonLevel() || !m_config->doTopPartonHistory()){ // For safety
    return;
  }
  
  TLorentzVector ttbar, top, tbar, lpos, lneg;

  const xAOD::EventInfo* eventInfo(nullptr);
  top::check( evtStore()->retrieve(eventInfo, m_config->sgKeyEventInfo()) , "Failed to retrieve EventInfo" );

  const xAOD::TruthEventContainer * truthEvent(nullptr);
  top::check( evtStore()->retrieve(truthEvent, m_config->sgKeyTruthEvent()) , "Failed to retrieve truth event container" );
  unsigned int truthEventSize = truthEvent->size();
  top::check( truthEventSize==1 , "Failed to retrieve truth PDF info - truth event container size is different from 1 (strange)" );

  ///-- Try to get the tops after FSR --///
  const xAOD::PartonHistoryContainer* topPartonCont = nullptr;
  top::check(evtStore()->event()->retrieve(topPartonCont, "TopPartonHistory"), "FAILURE");
  const xAOD::PartonHistory *topParton = topPartonCont->at(0);

  //Resolve tau event flag
  if(abs(topParton->auxdata<int>("MC_Wdecay1_from_t_pdgId"))    == 16 ||
     abs(topParton->auxdata<int>("MC_Wdecay2_from_t_pdgId"))    == 16 ||
     abs(topParton->auxdata<int>("MC_Wdecay1_from_tbar_pdgId")) == 16 ||
     abs(topParton->auxdata<int>("MC_Wdecay2_from_tbar_pdgId")) == 16 ){
    m_parton_is_tau = 1;
  }

  if (topParton->auxdata<float>("MC_t_afterFSR_pt")>0) {
    m_parton_t_afterFSR_pt     = topParton->auxdata<float>("MC_t_afterFSR_pt")/1000.;
    m_parton_t_afterFSR_eta    = topParton->auxdata<float>("MC_t_afterFSR_eta");
    m_parton_t_afterFSR_phi    = topParton->auxdata<float>("MC_t_afterFSR_phi");
    m_parton_t_afterFSR_m      = topParton->auxdata<float>("MC_t_afterFSR_m")/1000.;

    top.SetPtEtaPhiM(m_parton_t_afterFSR_pt,m_parton_t_afterFSR_eta, m_parton_t_afterFSR_phi, m_parton_t_afterFSR_m);
  }
  if (topParton->auxdata<float>("MC_tbar_afterFSR_pt")>0) {
    m_parton_tbar_afterFSR_pt  = topParton->auxdata<float>("MC_tbar_afterFSR_pt")/1000.;
    m_parton_tbar_afterFSR_eta = topParton->auxdata<float>("MC_tbar_afterFSR_eta");
    m_parton_tbar_afterFSR_phi = topParton->auxdata<float>("MC_tbar_afterFSR_phi");
    m_parton_tbar_afterFSR_m   = topParton->auxdata<float>("MC_tbar_afterFSR_m")/1000.;

    tbar.SetPtEtaPhiM(m_parton_tbar_afterFSR_pt, m_parton_tbar_afterFSR_eta, m_parton_tbar_afterFSR_phi, m_parton_tbar_afterFSR_m);
  }

  if (top.Pt()>0 && tbar.Pt()>0){
  	ttbar = top + tbar;	
  	m_parton_ttbar_pt  = ttbar.Pt();
  	m_parton_ttbar_eta = ttbar.Eta();
  	m_parton_ttbar_phi = ttbar.Phi();
  	m_parton_ttbar_m   = ttbar.M();
  }
 

  if (topParton->auxdata<int>("MC_Wdecay1_from_t_pdgId") == -11 || 
      topParton->auxdata<int>("MC_Wdecay1_from_t_pdgId") == -13 || 
      topParton->auxdata<int>("MC_Wdecay1_from_t_pdgId") == -15 ){
    m_parton_lbar_pt    = topParton->auxdata<float>("MC_Wdecay1_from_t_pt")/1000.;
    m_parton_lbar_eta   = topParton->auxdata<float>("MC_Wdecay1_from_t_eta");
    m_parton_lbar_phi   = topParton->auxdata<float>("MC_Wdecay1_from_t_phi");
    m_parton_lbar_m     = topParton->auxdata<float>("MC_Wdecay1_from_t_m")/1000.;
    m_parton_lbar_pdgid = topParton->auxdata<int>("MC_Wdecay1_from_t_pdgId");
  }
  if (topParton->auxdata<int>("MC_Wdecay2_from_t_pdgId") == -11 || 
      topParton->auxdata<int>("MC_Wdecay2_from_t_pdgId") == -13 || 
      topParton->auxdata<int>("MC_Wdecay2_from_t_pdgId") == -15 ){
    m_parton_lbar_pt    = topParton->auxdata<float>("MC_Wdecay2_from_t_pt")/1000.;
    m_parton_lbar_eta   = topParton->auxdata<float>("MC_Wdecay2_from_t_eta");
    m_parton_lbar_phi   = topParton->auxdata<float>("MC_Wdecay2_from_t_phi");
    m_parton_lbar_m     = topParton->auxdata<float>("MC_Wdecay2_from_t_m")/1000.;
    m_parton_lbar_pdgid = topParton->auxdata<int>("MC_Wdecay2_from_t_pdgId");
  }
  if (topParton->auxdata<int>("MC_Wdecay1_from_tbar_pdgId") == 11 || 
      topParton->auxdata<int>("MC_Wdecay1_from_tbar_pdgId") == 13 || 
      topParton->auxdata<int>("MC_Wdecay1_from_tbar_pdgId") == 15 ){
    m_parton_l_pt    = topParton->auxdata<float>("MC_Wdecay1_from_tbar_pt")/1000.;
    m_parton_l_eta   = topParton->auxdata<float>("MC_Wdecay1_from_tbar_eta");
    m_parton_l_phi   = topParton->auxdata<float>("MC_Wdecay1_from_tbar_phi");
    m_parton_l_m     = topParton->auxdata<float>("MC_Wdecay1_from_tbar_m")/1000.;
    m_parton_l_pdgid = topParton->auxdata<int>("MC_Wdecay1_from_tbar_pdgId");
  }
  if (topParton->auxdata<int>("MC_Wdecay2_from_tbar_pdgId") == 11 || 
      topParton->auxdata<int>("MC_Wdecay2_from_tbar_pdgId") == 13 || 
      topParton->auxdata<int>("MC_Wdecay2_from_tbar_pdgId") == 15 ){
    m_parton_l_pt    = topParton->auxdata<float>("MC_Wdecay2_from_tbar_pt")/1000.;
    m_parton_l_eta   = topParton->auxdata<float>("MC_Wdecay2_from_tbar_eta");
    m_parton_l_phi   = topParton->auxdata<float>("MC_Wdecay2_from_tbar_phi");
    m_parton_l_m     = topParton->auxdata<float>("MC_Wdecay2_from_tbar_m")/1000.;
    m_parton_l_pdgid = topParton->auxdata<int>("MC_Wdecay2_from_tbar_pdgId");
  }

  lneg.SetPtEtaPhiM(m_parton_l_pt,    m_parton_l_eta,    m_parton_l_phi,    m_parton_l_m);
  lpos.SetPtEtaPhiM(m_parton_lbar_pt, m_parton_lbar_eta, m_parton_lbar_phi, m_parton_lbar_m);

  m_parton_cosphi      = cos_phi(lpos, lneg, top, tbar, ttbar);
  m_parton_cos_kaxis_p = cos_theta_helicity(   top, top,  ttbar, lpos, +1);
  m_parton_cos_kaxis_m = cos_theta_helicity(   top, tbar, ttbar, lneg, -1);
  m_parton_cos_naxis_p = cos_theta_transverse( top, top,  ttbar, lpos, +1);
  m_parton_cos_naxis_m = cos_theta_transverse( top, tbar, ttbar, lneg, -1);
  m_parton_cos_raxis_p = cos_theta_raxis(      top, top,  ttbar, lpos, +1);
  m_parton_cos_raxis_m = cos_theta_raxis(      top, tbar, ttbar, lneg, -1);


  ///-- Get MC generator weights --///
  //const xAOD::EventInfo* eventInfo(nullptr);
  //top::check(evtStore()->retrieve(eventInfo, m_config->sgKeyEventInfo()), "Failed to retrieve EventInfo for loading of MCGenWeights!");
  m_parton_weight_nnlo = m_ttbarNNLO_Reweighter->GetWeight(top.Pt(), tbar.Pt(), ttbar.M(), ttbar.Pt());

  std::vector<float> temp_generator_weights = eventInfo->mcEventWeights();

  ///-- WARNING: THIS IS HARDCODED FOR THE FS NOMINAL PP8 SAMPLE ONLY! --///
  for(uint i = 0; i < temp_generator_weights.size(); ++i){
    if(i == 1   || // "muR = 1.0, muF = 2.0"
       i == 2   || // "muR = 1.0, muF = 0.5"
       i == 3   || // "muR = 2.0, muF = 1.0"
       i == 4   || // "muR = 0.5, muF = 1.0"
       i == 115 || // "PDF set = 90901"
       i == 116 || // "PDF set = 90902"
       i == 117 || // "PDF set = 90903"
       i == 118 || // "PDF set = 90904"
       i == 119 || // "PDF set = 90905"
       i == 120 || // "PDF set = 90906"
       i == 121 || // "PDF set = 90907"
       i == 122 || // "PDF set = 90908"
       i == 123 || // "PDF set = 90909"
       i == 124 || // "PDF set = 90910"
       i == 125 || // "PDF set = 90911"
       i == 126 || // "PDF set = 90912"
       i == 127 || // "PDF set = 90913"
       i == 128 || // "PDF set = 90914"
       i == 129 || // "PDF set = 90915"
       i == 130 || // "PDF set = 90916"
       i == 131 || // "PDF set = 90917"
       i == 132 || // "PDF set = 90918"
       i == 133 || // "PDF set = 90919"
       i == 134 || // "PDF set = 90920"
       i == 135 || // "PDF set = 90921"
       i == 136 || // "PDF set = 90922"
       i == 137 || // "PDF set = 90923"
       i == 138 || // "PDF set = 90924"
       i == 139 || // "PDF set = 90925"
       i == 140 || // "PDF set = 90926"
       i == 141 || // "PDF set = 90927"
       i == 142 || // "PDF set = 90928"
       i == 143 || // "PDF set = 90929"
       i == 144 || // "PDF set = 90930"
       i == 193 || // "Var3cUp"
       i == 194 || // "Var3cDown"
       i == 198 || // "isr:muRfac=1.0_fsr:muRfac=2.0"
       i == 199    // "isr:muRfac=1.0_fsr:muRfac=0.5"
       )
     m_parton_generator_weights.push_back(temp_generator_weights[i]);
  }

  EventSaverFlatNtuple::saveTruthEvent();    
  
}


void TopCommonDileptonEventSaver::saveParticleLevelEvent(const top::ParticleLevelEvent& plEvent){ // This is particle level
  
  // Reset everything
  initializeParticleLevelBranches();
  

  if (!is_MC || !m_config->doTopParticleLevel()){ // For safety?
    return;
  }

  TLorentzVector lep_pos, lep_neg, nu, nubar, b, bbar, top, tbar, ttbar, jet1, jet2;

  // Electrons
  for (const auto & elPtr : * plEvent.m_electrons) {
    m_particle_lep_pt.push_back(elPtr->pt()/1000.);
    m_particle_lep_eta.push_back(elPtr->eta());
    m_particle_lep_phi.push_back(elPtr->phi());
    m_particle_lep_m.push_back(elPtr->m()/1000.);
    m_particle_lep_pdgid.push_back(11.*elPtr->charge()*-1.);
  }
  
  // Muons
  for (const auto & muPtr : * plEvent.m_muons) {
    m_particle_lep_pt.push_back(muPtr->pt()/1000.);
    m_particle_lep_eta.push_back(muPtr->eta());
    m_particle_lep_phi.push_back(muPtr->phi());
    m_particle_lep_m.push_back(muPtr->m()/1000.);
    m_particle_lep_pdgid.push_back(13.*muPtr->charge()*-1.);
  }

  m_particle_met_ex  = plEvent.m_met->mpx()/1000.;
  m_particle_met_ey  = plEvent.m_met->mpy()/1000.;
  m_particle_met_met = plEvent.m_met->met()/1000.;
  m_particle_met_phi = plEvent.m_met->phi();

  //Grab info for pseudo-top reconstruction
  bool leptonsSet = false;
  if (m_particle_lep_pt.size() == 2){
    if (m_particle_lep_pdgid[0] < 0 && m_particle_lep_pdgid[1] > 0){
      lep_pos.SetPtEtaPhiM(m_particle_lep_pt[0], m_particle_lep_eta[0], m_particle_lep_phi[0], m_particle_lep_m[0]);
      lep_neg.SetPtEtaPhiM(m_particle_lep_pt[1], m_particle_lep_eta[1], m_particle_lep_phi[1], m_particle_lep_m[1]);
      leptonsSet = true;
   } else if (m_particle_lep_pdgid[0] > 0 && m_particle_lep_pdgid[1] < 0){
      lep_pos.SetPtEtaPhiM(m_particle_lep_pt[1], m_particle_lep_eta[1], m_particle_lep_phi[1], m_particle_lep_m[1]);
      lep_neg.SetPtEtaPhiM(m_particle_lep_pt[0], m_particle_lep_eta[0], m_particle_lep_phi[0], m_particle_lep_m[0]);
      leptonsSet = true;
    } else {
      std::cout << "WARNING: Particle same-sign lepton event!" << std::endl;
    }
  }

  // Get the particle level neutrinos
  const xAOD::TruthParticleContainer* neutrinos{nullptr};
  top::check(evtStore()->retrieve(neutrinos, "TruthNeutrinos"), "Failed to retrieve Truth Neutrinos");
  
  ///-- Here we rely on the ordering of the Neutrinos, i.e. that the first two are the hard ones and the subsequent ones come from further tau decay or from jets --///
  bool nu_set    = false;
  bool nubar_set = false;
  if (neutrinos != nullptr) {
    for (const auto& neutrino : *neutrinos) {
      if(!nu_set || !nubar_set){
        if( neutrino->pdgId() > 0 && !nu_set){
          nu.SetPtEtaPhiE(neutrino->pt()/1000., neutrino->eta(), neutrino->phi(), neutrino->e()/1000.);
          nu_set               = true;

        }
        if( neutrino->pdgId() < 0 && !nubar_set){
          nubar.SetPtEtaPhiE(neutrino->pt()/1000., neutrino->eta(), neutrino->phi(), neutrino->e()/1000.);
      	  nubar_set              = true;
        }
      }
    }
  }


  // Grab a little jet info just for the pseudo-top reconstruction
  std::vector<float> tmp_particle_bjet_pt;
  std::vector<float> tmp_particle_bjet_eta;
  std::vector<float> tmp_particle_bjet_phi;
  std::vector<float> tmp_particle_bjet_m;
  std::vector<float> tmp_particle_bjet_ghosts;
  int tmp_particle_bjet_n = 0;
  std::vector<float> tmp_particle_jet_pt;
  std::vector<float> tmp_particle_jet_eta;
  std::vector<float> tmp_particle_jet_phi;
  std::vector<float> tmp_particle_jet_m;
  std::vector<float> tmp_particle_jet_ghosts;
  int tmp_particle_jet_n = 0;
  for (const auto & jetPtr : * plEvent.m_jets) {
    tmp_particle_jet_pt.push_back(        jetPtr->pt()/1000.);
    tmp_particle_jet_eta.push_back(       jetPtr->eta());
    tmp_particle_jet_phi.push_back(       jetPtr->phi());
    tmp_particle_jet_m.push_back(         jetPtr->m()/1000.);
    tmp_particle_jet_ghosts.push_back(    jetPtr->auxdata<int>( "GhostBHadronsFinalCount" ));
    ++tmp_particle_jet_n;
    if( jetPtr->auxdata<int>( "GhostBHadronsFinalCount" ) > 0){
      tmp_particle_bjet_pt.push_back(     jetPtr->pt()/1000.);
      tmp_particle_bjet_eta.push_back(    jetPtr->eta());
      tmp_particle_bjet_phi.push_back(    jetPtr->phi());
      tmp_particle_bjet_m.push_back(      jetPtr->m()/1000.);
      tmp_particle_bjet_ghosts.push_back( jetPtr->auxdata<int>( "GhostBHadronsFinalCount" ));
      ++tmp_particle_bjet_n;
    }
  }


  // Do the pseudo-top reconstruction
  if(leptonsSet == true && tmp_particle_jet_n >= 2){

    if(tmp_particle_bjet_n > 1){ //if 2 bjets, use them both                                             
      jet1.SetPtEtaPhiM(tmp_particle_bjet_pt[0], tmp_particle_bjet_eta[0], tmp_particle_bjet_phi[0], tmp_particle_bjet_m[0]);
      jet2.SetPtEtaPhiM(tmp_particle_bjet_pt[1], tmp_particle_bjet_eta[1], tmp_particle_bjet_phi[1], tmp_particle_bjet_m[1]);
    }
    else if(tmp_particle_bjet_n > 0){ // use 1 bjet              
      jet1.SetPtEtaPhiM(tmp_particle_bjet_pt[0], tmp_particle_bjet_eta[0], tmp_particle_bjet_phi[0], tmp_particle_bjet_m[0]);
      if (tmp_particle_jet_pt[0] != tmp_particle_bjet_pt[0]){ // make sure other isnt also bjet            
	      jet2.SetPtEtaPhiM(tmp_particle_jet_pt[0], tmp_particle_jet_eta[0], tmp_particle_jet_phi[0], tmp_particle_jet_m[0]);
      }
      else { // use different jet       
	      jet2.SetPtEtaPhiM(tmp_particle_jet_pt[1], tmp_particle_jet_eta[1], tmp_particle_jet_phi[1], tmp_particle_jet_m[1]);
      }
    }
    else{ // don't use any bjets because none exist
      jet1.SetPtEtaPhiM(tmp_particle_jet_pt[0], tmp_particle_jet_eta[0], tmp_particle_jet_phi[0], tmp_particle_jet_m[0]);
      jet2.SetPtEtaPhiM(tmp_particle_jet_pt[1], tmp_particle_jet_eta[1], tmp_particle_jet_phi[1], tmp_particle_jet_m[1]);
    }

    TLorentzVector Wp, Wm, top_a, tbar_a, top_b, tbar_b;

    Wp     = lep_pos + nu;
    Wm     = lep_neg + nubar;

    top_a  = Wp + jet1;
    top_b  = Wp + jet2;

    tbar_a = Wm + jet2;
    tbar_b = Wm + jet1;

    // loose mass decision for combinatorics
    double diff_a = fabs(top_a.M() - 172.5) + fabs(tbar_a.M() - 172.5);
    double diff_b = fabs(top_b.M() - 172.5) + fabs(tbar_b.M() - 172.5);
    if(diff_a < diff_b){ // A is the correct combination                                               
      top  = top_a;
      tbar = tbar_a;
    } else if (diff_a > diff_b){ // B is the right combination                                         
      top  = top_b;
      tbar = tbar_b;
    } else {
      top  = top_a;
      tbar = tbar_a;
    }
    
    m_particle_t_pt       = top.Pt();
    m_particle_t_eta      = top.Eta();
    m_particle_t_phi      = top.Phi();
    m_particle_t_m        = top.M();

    m_particle_tbar_pt    = tbar.Pt();
    m_particle_tbar_eta   = tbar.Eta();
    m_particle_tbar_phi   = tbar.Phi();
    m_particle_tbar_m     = tbar.M();

    m_particle_ttbar_pt   = ttbar.Pt();
    m_particle_ttbar_eta  = ttbar.Eta();
    m_particle_ttbar_phi  = ttbar.Phi();
    m_particle_ttbar_m    = ttbar.M();   
  }

  EventSaverFlatNtuple::saveParticleLevelEvent(plEvent);

}


void TopCommonDileptonEventSaver::saveEvent(const top::Event& event){ //This is reco level


  if (m_config->saveOnlySelectedEvents() && !event.m_saveEvent){
    return;
  }


  // Reset eveything
  initializeRecoLevelBranches();


  const xAOD::EventInfo* eventInfo(nullptr);
  top::check( evtStore()->retrieve(eventInfo, m_config->sgKeyEventInfo()) , "Failed to retrieve EventInfo" );

  if (event.m_hashValue == m_config->nominalHashValue()){
    //ONLY SAVE WEIGHTS IF TREE IS NOMINAL
    if(m_config->isMC() == true){
    std::vector<float> temp_generator_weights = eventInfo->mcEventWeights();

    ///-- WARNING: THIS IS HARDCODED FOR THE FS NOMINAL PP8 SAMPLE ONLY! --///
    for(uint i = 0; i < temp_generator_weights.size(); ++i){
      if(i == 1   || // "muR = 1.0, muF = 2.0"
        i == 2   || // "muR = 1.0, muF = 0.5"
        i == 3   || // "muR = 2.0, muF = 1.0"
        i == 4   || // "muR = 0.5, muF = 1.0"
        i == 115 || // "PDF set = 90901"
        i == 116 || // "PDF set = 90902"
        i == 117 || // "PDF set = 90903"
        i == 118 || // "PDF set = 90904"
        i == 119 || // "PDF set = 90905"
        i == 120 || // "PDF set = 90906"
        i == 121 || // "PDF set = 90907"
        i == 122 || // "PDF set = 90908"
        i == 123 || // "PDF set = 90909"
        i == 124 || // "PDF set = 90910"
        i == 125 || // "PDF set = 90911"
        i == 126 || // "PDF set = 90912"
        i == 127 || // "PDF set = 90913"
        i == 128 || // "PDF set = 90914"
        i == 129 || // "PDF set = 90915"
        i == 130 || // "PDF set = 90916"
        i == 131 || // "PDF set = 90917"
        i == 132 || // "PDF set = 90918"
        i == 133 || // "PDF set = 90919"
        i == 134 || // "PDF set = 90920"
        i == 135 || // "PDF set = 90921"
        i == 136 || // "PDF set = 90922"
        i == 137 || // "PDF set = 90923"
        i == 138 || // "PDF set = 90924"
        i == 139 || // "PDF set = 90925"
        i == 140 || // "PDF set = 90926"
        i == 141 || // "PDF set = 90927"
        i == 142 || // "PDF set = 90928"
        i == 143 || // "PDF set = 90929"
        i == 144 || // "PDF set = 90930"
        i == 193 || // "Var3cUp"
        i == 194 || // "Var3cDown"
        i == 198 || // "isr:muRfac=1.0_fsr:muRfac=2.0"
        i == 199    // "isr:muRfac=1.0_fsr:muRfac=0.5"
        )
      m_generator_weights.push_back(temp_generator_weights[i]);
    }
    }
  }



  // Electrons
  for (const auto* const electron : event.m_electrons) {

    unsigned int j(0);

    m_lep_pt.push_back(electron->pt()/1000.);
    m_lep_eta.push_back(electron->eta());
    //m_lep_cl_eta.push_back(electron->caloCluster()->etaBE(2));
    m_lep_phi.push_back(electron->phi());
    m_lep_e.push_back(electron->e()/1000.);
    if(electron->charge() > 0){
      m_lep_pdgid.push_back(-11);
    } else {
      m_lep_pdgid.push_back(11);
    }

    if(m_config->isMC() == true){
      static SG::AuxElement::Accessor<int> typeel("truthType");
      static SG::AuxElement::Accessor<int> origel("truthOrigin");
      static SG::AuxElement::Accessor<int> firstEgMotherTruthType("firstEgMotherTruthType");
      static SG::AuxElement::Accessor<int> firstEgMotherTruthOrigin("firstEgMotherTruthOrigin");
      static SG::AuxElement::Accessor<int> firstEgMotherPdgId("firstEgMotherPdgId");
      auto elTruth = xAOD::TruthHelpers::getTruthParticle(*electron);
      bool eltruth_good = false;
      if(elTruth){
	      if(elTruth->isAvailable<int>("pdgId")){
	        eltruth_good = true;
	        m_lep_true_pdgid.push_back(elTruth->pdgId());
	      }
      }
      if(!eltruth_good){
	      m_lep_true_pdgid.push_back(-999);
      }

      if (typeel.isAvailable(*electron)) { 
	      m_lep_true_type.push_back(typeel(*electron)); 
	      if(electron->auxdata<int>("truthType") != 2){ 
          m_fakeEvent = true;
        }
      } 
      else { 
        m_lep_true_type.push_back(-99.0); 
      }

      if (origel.isAvailable(*electron)){ 
        m_lep_true_origin.push_back(origel(*electron)); 
      } else {
        m_lep_true_origin.push_back(-99.0); 
      }

      if (typeel.isAvailable(*electron) && 
        origel.isAvailable(*electron) && 
        firstEgMotherTruthType.isAvailable(*electron) && 
        firstEgMotherTruthOrigin.isAvailable(*electron) && 
        firstEgMotherPdgId.isAvailable(*electron)){
	        m_lep_true_isPrompt.push_back( (isPromptElectron( typeel(*electron), origel(*electron), firstEgMotherTruthType(*electron), firstEgMotherTruthOrigin(*electron), firstEgMotherPdgId(*electron), electron->charge())).first );
          m_lep_true_isChargeFl.push_back( (isPromptElectron( typeel(*electron), origel(*electron), firstEgMotherTruthType(*electron), firstEgMotherTruthOrigin(*electron), firstEgMotherPdgId(*electron), electron->charge())).second );
      }
      else { 
        m_lep_true_isPrompt.push_back(-99.0); 
      }
    }

    ++j;
  }


  // Muons
  for (const auto* const muon : event.m_muons) {

    unsigned int k(0);

    m_lep_pt.push_back(muon->pt()/1000.);
    m_lep_eta.push_back(muon->eta());
    //m_lep_cl_eta.push_back(-99.0);
    m_lep_phi.push_back(muon->phi());
    m_lep_e.push_back(muon->e()/1000.);
    //m_lep_charge.push_back(muon->charge());
    if(muon->charge() > 0){
      m_lep_pdgid.push_back(-13);
    } else {
      m_lep_pdgid.push_back(13);
    }

    if(m_config->isMC() == true){
      static SG::AuxElement::Accessor<int> typemu("truthType");
      static SG::AuxElement::Accessor<int> origmu("truthOrigin");
      const xAOD::TrackParticle* part = muon->primaryTrackParticle();
      bool mutruth_good = false;
      if(part){
      	const xAOD::TruthParticle* muTruth = xAOD::TruthHelpers::getTruthParticle(*part);
	      if(muTruth){
	        if(muTruth->isAvailable<int>("pdgId")){
	          mutruth_good = true;
	          m_lep_true_pdgid.push_back(muTruth->pdgId());
	        }
	      }
      }
      if(!mutruth_good){
      	m_lep_true_pdgid.push_back(-999);
      }

      if (typemu.isAvailable(*muon)) {
	      m_lep_true_type.push_back(typemu(*muon));
	      if(muon->primaryTrackParticle()->auxdata<int>("truthType") != 6){
          m_fakeEvent = true;
        }
      } else {
        m_lep_true_type.push_back(-99.0);
      }

      if (origmu.isAvailable(*muon)) {
        m_lep_true_origin.push_back(origmu(*muon));
      } else {
        m_lep_true_origin.push_back(-99.0);
      }

      if( typemu.isAvailable(*muon) && origmu.isAvailable(*muon)){ 
        m_lep_true_isPrompt.push_back(isPromptMuon(typemu(*muon),origmu(*muon))); 
      } else { 
        m_lep_true_isPrompt.push_back(-99.0);
      }
    }

    ++k;
  }

  std::vector<TLorentzVector> bjets;    
  std::vector<TLorentzVector> ljets;      

  for (const auto* const jetPtr : event.m_jets){

    TLorentzVector jet = TLorentzVector();
    jet.SetPtEtaPhiE(jetPtr->pt()/1000., jetPtr->eta(), jetPtr->phi(), jetPtr->e()/1000.);

    if(jetPtr->isAvailable<char>("isbtagged_DL1r_FixedCutBEff_85")){
      if(jetPtr->auxdataConst<char>("isbtagged_DL1r_FixedCutBEff_85")){
        m_bjet_pt.push_back( jet.Pt());
        m_bjet_eta.push_back(jet.Eta());
        m_bjet_phi.push_back(jet.Phi());
        m_bjet_e.push_back(  jet.E());
        bjets.push_back(jet);      
      } else {
        m_ljet_pt.push_back( jet.Pt());
        m_ljet_eta.push_back(jet.Eta());
        m_ljet_phi.push_back(jet.Phi());
        m_ljet_e.push_back(  jet.E());
        ljets.push_back(jet);              
      }
    }
  }

  // Let's do the MET
  float m_met_ex    = event.m_met->mpx()/1000.;
  float m_met_ey    = event.m_met->mpy()/1000.;
//float m_met_sumet = event.m_met->sumet()/1000.;
  m_met_met         = event.m_met->met()/1000.;
  m_met_phi         = event.m_met->phi();

  ///-- For top reconstruction 
  TLorentzVector lepton_pos;
  TLorentzVector lepton_neg;
  TLorentzVector b;
  TLorentzVector bbar;


  if(m_lep_pdgid.size() == 2){
    if(m_lep_pdgid[0] > 0){
      lepton_neg.SetPtEtaPhiE(m_lep_pt[0], m_lep_eta[0], m_lep_phi[0], m_lep_e[0]);
      lepton_pos.SetPtEtaPhiE(m_lep_pt[1], m_lep_eta[1], m_lep_phi[1], m_lep_e[1]);    
    } else {
      lepton_neg.SetPtEtaPhiE(m_lep_pt[1], m_lep_eta[1], m_lep_phi[1], m_lep_e[1]);
      lepton_pos.SetPtEtaPhiE(m_lep_pt[0], m_lep_eta[0], m_lep_phi[0], m_lep_e[0]);      
    }
  }  

  TLorentzVector dilep = lepton_pos + lepton_neg;
  m_dilep_m = dilep.M();

  if (event.m_electrons.size() == 2 || event.m_muons.size() == 2){
    if (fabs(m_dilep_m - 91.0) < 10 || m_dilep_m < 15){
      m_passSelection = false;
    } else {
      m_passSelection = true;
    }
  } else if (event.m_electrons.size() == 1 && event.m_muons.size() == 1){
    m_passSelection = true;
  } else {
    m_passSelection = false;
    std::cout << "ERROR: You really shouldn't be getting to this point!" << std::endl;
  }


  bool use_1_bjet_region = true;

  bool b_set    = false;
  bool bbar_set = false;

  for(uint i = 0; i < bjets.size(); ++i){
    if(!b_set){
      b_set = true;
      b = bjets[i];
      continue;
    }
    else if(!bbar_set) {
      bbar_set = true;
      bbar = bjets[i];
      continue;
    }   
  }

  if(use_1_bjet_region && (!bbar_set || !b_set) && ljets.size() > 0){
    if(!b_set){
      b_set = true;
      b = ljets[0];
    }
    else if(!bbar_set) {
      bbar_set = true;
      bbar = ljets[0];
    }   
  }


  if(m_lep_pdgid.size() == 2 && b_set && bbar_set){

    TopDileptonReconstruction topDilepRecoEM = TopDileptonReconstruction();
    TopDileptonReconstruction topDilepRecoNW = TopDileptonReconstruction();    

    if(m_doEllipseMethod){
      topDilepRecoEM.RunEM();
    }
    if(m_doNeutrinoWeighter){
      topDilepRecoNW.RunNW();
    }


    TLorentzVector top, tbar, ttbar;

    uint nsmears = 20;
    std::vector<float> mtop_smears;
    std::vector<float> mtbar_smears;
    std::vector<float> mWpos_smears;
    std::vector<float> mWneg_smears;

    TRandom3 random = TRandom3(top::EventSaverFlatNtuple::eventNumber());    


    ///-- Setup the smears --///
    for( uint i = 0; i < nsmears; ++i ){
      double mtop  = random.Gaus(172.5,  1.480);
      double mtbar = random.Gaus(172.5,  1.480);
      double mWpos = random.Gaus(80.379, 2.085);
      double mWneg = random.Gaus(80.379, 2.085);

      mtop_smears.push_back(mtop);
      mtbar_smears.push_back(mtbar);
      mWpos_smears.push_back(mWpos);
      mWneg_smears.push_back(mWneg);
    }


    ///-- Run the EM first --///
    for( uint i = 0; i < nsmears; ++i ){ 
      m_top_reco_method = 0;

      float mtop  = mtop_smears[i];
      float mtbar = mtbar_smears[i];
      float mWpos = mWpos_smears[i];
      float mWneg = mWneg_smears[i];      

      topDilepRecoEM.Reconstruct(lepton_pos, lepton_neg, b, bbar, m_met_ex, m_met_ey, mtop, mtbar, mWpos, mWneg);
      topDilepRecoEM.Reconstruct(lepton_pos, lepton_neg, bbar, b, m_met_ex, m_met_ey, mtop, mtbar, mWpos, mWneg);    
      top  = topDilepRecoEM.GetEMTop();
      tbar = topDilepRecoEM.GetEMTbar();
    }


    ///-- Now run NW only if the EM completely failed --///
    if(top.Pt() < 0.1 and tbar.Pt() < 0.1){
      m_top_reco_method = 1;

      for( uint i = 0; i < nsmears; ++i ){ 
        float mtop  = mtop_smears[i];
        float mtbar = mtbar_smears[i];
        float mWpos = mWpos_smears[i];
        float mWneg = mWneg_smears[i];  

        topDilepRecoNW.Reconstruct(lepton_pos, lepton_neg, b, bbar, m_met_ex, m_met_ey, mtop, mtbar, mWpos, mWneg);
        topDilepRecoNW.Reconstruct(lepton_pos, lepton_neg, bbar, b, m_met_ex, m_met_ey, mtop, mtbar, mWpos, mWneg);    
        top  = topDilepRecoNW.GetNWTop();
        tbar = topDilepRecoNW.GetNWTbar();
      }
    }


    ///-- Now run rudimentary reconstruction if both failed --/// 
    if(top.Pt() < 0.1 and tbar.Pt() < 0.1){
      m_top_reco_method = 2;

      /// Pair leptons and bs by closest DR , only runs if EM and NW failed
      if (lepton_pos.DeltaR(b) < lepton_pos.DeltaR(bbar)){
        top  = lepton_pos + b;
        tbar = lepton_neg + bbar;
      } else {
        top  = lepton_pos + bbar;
        tbar = lepton_neg + b;
      }
    }


    /// Fill Branches ///
    if(top.Pt() > 0 and tbar.Pt() > 0){
      ttbar = top + tbar;
      m_t_pt        = top.Pt();
      m_t_eta       = top.Eta();
      m_t_phi       = top.Phi();
      m_t_m         = top.M();
      m_tbar_pt     = tbar.Pt();
      m_tbar_eta    = tbar.Eta();
      m_tbar_phi    = tbar.Phi();
      m_tbar_m      = tbar.M();
      m_ttbar_pt    = ttbar.Pt();
      m_ttbar_eta   = ttbar.Eta();
      m_ttbar_phi   = ttbar.Phi();
      m_ttbar_m     = ttbar.M();
      m_cosphi      = cos_phi(lepton_pos, lepton_neg, top, tbar, ttbar);
      m_cos_kaxis_p = cos_theta_helicity(   top, top,  ttbar, lepton_pos, +1);
      m_cos_kaxis_m = cos_theta_helicity(   top, tbar, ttbar, lepton_neg, -1);
      m_cos_naxis_p = cos_theta_transverse( top, top,  ttbar, lepton_pos, +1);
      m_cos_naxis_m = cos_theta_transverse( top, tbar, ttbar, lepton_neg, -1);
      m_cos_raxis_p = cos_theta_raxis(      top, top,  ttbar, lepton_pos, +1);
      m_cos_raxis_m = cos_theta_raxis(      top, tbar, ttbar, lepton_neg, -1);
    }

  }

  ///-- We want to set all the original variables too and the call fill on the TTree --///
  EventSaverFlatNtuple::saveEvent(event);
}


int TopCommonDileptonEventSaver::isBranchStored(top::TreeManager const *treeManager, const std::string &variableName)
{
  /* Remove branches from output n-tuple using pattern matching */

  // Could we replace all of the if statements in here with a check against a variable 'blacklist' which can be read-in dynamically?
  // In future this should be replaced by the new FilterBranch feature of AnalysisTop!

  if(false){ // turn off, this is just for some brute-force debugging
    std::cout << "treeManager->name() = " << treeManager->name() << "  ,  variableName = " << variableName << std::endl;
  }

  
  // individual lepton (debugging) SFs
  if (variableName.find("weight_indiv_SF") != std::string::npos){
    return 0;
  }

  if (variableName.find("weight_globalLeptonTriggerSF_") != std::string::npos){
    return 0;
  }
  if (variableName.find("weight_oldTriggerSF") != std::string::npos){
    return 0;
  }

  //kill weight_mc since it's also the first entry in mc_generator_weights
  //if ((variableName.find("weight_mc") != std::string::npos) && (treeManager->name() != "truth")){
  //return 0;
  //}

  if (variableName.find("met_met") != std::string::npos){
    if(!(variableName.find("d_") != std::string::npos)){
      return 0;
    }
  }
  
  if (variableName.find("met_phi") != std::string::npos){
    if(!(variableName.find("d_") != std::string::npos)){
      return 0;
    }
  }

  // we are overwriting these variables with our own, to swap from vector<char> to vector<int>
  if ( variableName.find("jet_isbtagged_DL1r_85") != std::string::npos
      || variableName.find("jet_isTrueHS") != std::string::npos
      || variableName.find("jet_passfjvt") != std::string::npos)
    {
      if(!(variableName.find("d_jet") != std::string::npos)){
	return 0;
      }
    }

  // we are overwrriting these variables as a simple name change for consistency with the variables above which undergo a type change
  if(variableName.find("jet_") != std::string::npos){
    if(!(variableName.find("d_jet") != std::string::npos)){
      if(!(variableName.find("jet_passfjvt") != std::string::npos)){
	if(!(variableName.find("failJvt_jet") != std::string::npos)){
	  return 0;  
	}
      }
    }
  }

  
  // we are also overwrriting the trigger variables just o change their type from vector<char> to vector<int>
  if(variableName.find("trigMatch") != std::string::npos){
    if(!(variableName.find("d_") != std::string::npos)){
      return 0;
    }
  }

  // electron/muon variables which have been combined into 'lepton' variables
  if(variableName.find("el_pt") != std::string::npos
     || variableName.find("el_eta") != std::string::npos
     || variableName.find("el_phi") != std::string::npos
     || variableName.find("el_e") != std::string::npos
     || variableName.find("el_charge") != std::string::npos
     || variableName.find("el_topoetcone20") != std::string::npos
     || variableName.find("el_ptvarcone20") != std::string::npos
     || variableName.find("el_d0sig") != std::string::npos
     || variableName.find("el_delta_z0_sintheta") != std::string::npos
     || variableName.find("el_cl_eta") != std::string::npos
     || variableName.find("el_CF") != std::string::npos
     || variableName.find("el_ECIDS") != std::string::npos
     || variableName.find("el_ECIDSResult") != std::string::npos
     || variableName.find("el_true_type") != std::string::npos
     || variableName.find("el_true_origin") != std::string::npos
     || variableName.find("el_true_isPrompt") != std::string::npos
     || variableName.find("el_true_firstEgMotherTruthType") != std::string::npos
     || variableName.find("el_true_firstEgMotherTruthOrigin") != std::string::npos
     || variableName.find("el_true_firstEgMotherPdgId") != std::string::npos
     || variableName.find("el_true_isChargeFl") != std::string::npos
     ){return 0;}
  if(variableName.find("mu_pt") != std::string::npos
     || variableName.find("mu_eta") != std::string::npos
     || variableName.find("mu_phi") != std::string::npos
     || variableName.find("mu_e") != std::string::npos
     || variableName.find("mu_charge") != std::string::npos
     || variableName.find("mu_topoetcone20") != std::string::npos
     || variableName.find("mu_ptvarcone20") != std::string::npos
     || variableName.find("mu_d0sig") != std::string::npos
     || variableName.find("mu_delta_z0_sintheta") != std::string::npos
     || variableName.find("mu_true_type") != std::string::npos
     || variableName.find("mu_true_origin") != std::string::npos
     || variableName.find("mu_true_isPrompt") != std::string::npos
     ){return 0;}

  // uncessesary / renamed truth record info
  if (variableName.find("MC_") != std::string::npos){
    if (!(variableName.find("mc_generator_weight") != std::string::npos)){ // don't remove mc_weights by accident
      return 0;
    }
  }

  // jet b-tagging variables (no longer needed)
  if (variableName.find("jet_ip3dsv1") != std::string::npos
      || variableName.find("jet_mv2c00") != std::string::npos
      || variableName.find("jet_mv2c20") != std::string::npos
      //|| variableName.find("jet_DL1") != std::string::npos
      || variableName.find("jet_DL1_") != std::string::npos
      || variableName.find("jet_DL1r") != std::string::npos
      || variableName.find("jet_MV2c10mu") != std::string::npos
      || variableName.find("jet_MV2c100") != std::string::npos
      || variableName.find("jet_MV2cl100") != std::string::npos
      || variableName.find("jet_MV2c10rnn") != std::string::npos
      || variableName.find("jet_MV2r") != std::string::npos
      || variableName.find("jet_MV2rmu")!= std::string::npos){
    return 0;
  }

  // unused top-tagging variables
  if (variableName.find("ljet_isTopTagged") != std::string::npos
      || variableName.find("ljet_isWTagged") != std::string::npos
      || variableName.find("ljet_isZTagged") != std::string::npos
      || variableName == "ljet_sd12"){
    return 0;
  }

  // unused KLFitter variables
  if (variableName.find("klfitter") != std::string::npos){
    return 0;
  }

  // This gets used in the rare case that we want to save the systematic variations but not the nominal (i.e when the grid jobs take up too much RAM)
  if(m_skipNominal && treeManager->name() == "nominal"){
    return 0;
  }

  return -1;
}


float TopCommonDileptonEventSaver::cos_theta_helicity(TLorentzVector top, TLorentzVector parent_t, TLorentzVector ttbar, TLorentzVector lep, float sign){

  TVector3 boost_to_ttbar = ttbar.BoostVector();
  boost_to_ttbar *= -1.;
    
  parent_t.Boost(boost_to_ttbar);
  top.Boost(boost_to_ttbar);
  lep.Boost(boost_to_ttbar);
    
  TVector3 boost_to_parent_t = parent_t.BoostVector();
  boost_to_parent_t *= -1.;

  lep.Boost(boost_to_parent_t);
  
  TVector3 k_vector = top.Vect().Unit();
  k_vector *= sign;
  float theta = lep.Vect().Unit().Dot(k_vector);  

  //-- If we have a nan move to this to the overflow bins --//
  if (isnan(theta)){
    return -55.;
  } else {
    return theta;
  }
}


float TopCommonDileptonEventSaver::cos_theta_transverse(TLorentzVector top, TLorentzVector parent_t, TLorentzVector ttbar, TLorentzVector lep, float sign){

  TVector3 boost_to_ttbar = ttbar.BoostVector();
  boost_to_ttbar *= -1.;
    
  parent_t.Boost(boost_to_ttbar);
  top.Boost(boost_to_ttbar);
  lep.Boost(boost_to_ttbar);
    
  TVector3 boost_to_parent_t = parent_t.BoostVector();
  boost_to_parent_t *= -1.;

  lep.Boost(boost_to_parent_t);
  
  TVector3 k_vector = top.Vect().Unit();
  k_vector *= sign;
  TVector3 p_vector(0,0,1);

  float y = p_vector.Dot(k_vector);
  float r = pow((1. - y*y),0.5);

  TVector3 n_vector = (1./r)*(p_vector.Cross(k_vector)); ///-- Should this be Unit Vector? --///

  if (sign == 1){
    ///-- We're in the a axis --///
    if(y > 0) n_vector *= 1.;
    if(y < 0) n_vector *= -1.;
  } else if (sign == -1){
    ///-- We're in the b axis --///
    if(y > 0) n_vector *= -1.;
    if(y < 0) n_vector *= 1.;
  }

  float theta = lep.Vect().Unit().Dot(n_vector);

  //-- If we have a nan move to this to the overflow bins --//
  if (isnan(theta)){
    return -55.;
  } else {
    return theta;
  }
}


float TopCommonDileptonEventSaver::cos_theta_raxis(TLorentzVector top, TLorentzVector parent_t, TLorentzVector ttbar, TLorentzVector lep, float sign){

  TVector3 boost_to_ttbar = ttbar.BoostVector();
  boost_to_ttbar *= -1.;
    
  parent_t.Boost(boost_to_ttbar);
  top.Boost(boost_to_ttbar);
  lep.Boost(boost_to_ttbar);
    
  TVector3 boost_to_parent_t = parent_t.BoostVector();
  boost_to_parent_t *= -1.;

  lep.Boost(boost_to_parent_t);
  
  TVector3 k_vector = top.Vect().Unit();
  k_vector *= sign;
  TVector3 p_vector(0,0,1);

  float y = p_vector.Dot(k_vector);
  float r = pow((1. - y*y),0.5);

  TVector3 r_vector = (1./r)*(p_vector - y*k_vector);
  if (sign == 1){
    ///-- We're in the a axis --///
    if(y > 0) r_vector *= 1.;
    if(y < 0) r_vector *= -1.;
  } else if (sign == -1){
    ///-- We're in the b axis --///
    if(y > 0) r_vector *= -1.;
    if(y < 0) r_vector *= 1.;
  }

  float theta = lep.Vect().Unit().Dot(r_vector);

  //-- If we have a nan move to this to the overflow bins --//
  if (isnan(theta)){
    return -55.;
  } else {
    return theta;
  }
}


float TopCommonDileptonEventSaver::cos_phi(TLorentzVector lep_pos, TLorentzVector lep_neg, TLorentzVector top, TLorentzVector tbar, TLorentzVector ttbar){

  // Boost everything to ttbar
  TVector3 boost_to_ttbar = ttbar.BoostVector();
  boost_to_ttbar *= -1;

  // Boost leptons to their parent top frames
  TVector3 boost_to_top = top.BoostVector();
  boost_to_top *= -1;
  lep_pos.Boost(boost_to_top); 

  TVector3 boost_to_tbar = tbar.BoostVector();
  boost_to_tbar *= -1;
  lep_neg.Boost(boost_to_tbar); 

  //Take cosine of deltaphi
  float D = lep_pos.Vect().Unit().Dot(lep_neg.Vect().Unit());

  return D;
}
