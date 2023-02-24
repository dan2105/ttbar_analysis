#include "TopCommonDilepton/TopDileptonReconstruction.h"

//amespace top {

TopDileptonReconstruction::~TopDileptonReconstruction(){}

TopDileptonReconstruction::TopDileptonReconstruction(){
  Reset();
}

void TopDileptonReconstruction::Reset(){

  m_NW_tops.clear();
  m_EM_tops.clear();
  m_SN_tops.clear();
  m_NW_tbars.clear();
  m_EM_tbars.clear();
  m_SN_tbars.clear();
  m_NW_nus.clear();
  m_EM_nus.clear();
  m_SN_nus.clear();
  m_NW_nubars.clear();
  m_EM_nubars.clear();
  m_SN_nubars.clear();
  m_NW_Wposs.clear();
  m_EM_Wposs.clear();
  m_SN_Wposs.clear();
  m_NW_Wnegs.clear();
  m_EM_Wnegs.clear();
  m_SN_Wnegs.clear();

  m_NW_weights.clear();
  NWhighestWeight = 0;
  EMlowestMass    = 999999999999.;
  SNlowestMass    = 999999999999.;

  m_NW_highestWeightTop   = TLorentzVector();
  m_NW_highestWeightTbar  = TLorentzVector();
  m_EM_averageTop         = TLorentzVector();  
  m_EM_averageTbar        = TLorentzVector();
  m_EM_lowestMassTop      = TLorentzVector();  
  m_EM_lowestMassTbar     = TLorentzVector();     
  m_SN_averageTop         = TLorentzVector();  
  m_SN_averageTbar        = TLorentzVector();  
  m_SN_lowestMassTop      = TLorentzVector();  
  m_SN_lowestMassTbar     = TLorentzVector(); 

  m_NW_highestWeightNu    = TLorentzVector();
  m_NW_highestWeightNubar = TLorentzVector();
  m_EM_averageNu          = TLorentzVector();  
  m_EM_averageNubar       = TLorentzVector();
  m_EM_lowestMassNu       = TLorentzVector();  
  m_EM_lowestMassNubar    = TLorentzVector();     
  m_SN_averageNu          = TLorentzVector();  
  m_SN_averageNubar       = TLorentzVector();  
  m_SN_lowestMassNu       = TLorentzVector();  
  m_SN_lowestMassNubar    = TLorentzVector(); 

  m_NW_highestWeightWpos  = TLorentzVector();
  m_NW_highestWeightWneg  = TLorentzVector();
  m_EM_averageWpos        = TLorentzVector();  
  m_EM_averageWneg        = TLorentzVector();
  m_EM_lowestMassWpos     = TLorentzVector();  
  m_EM_lowestMassWneg     = TLorentzVector();    
  m_SN_averageWpos        = TLorentzVector();  
  m_SN_averageWneg        = TLorentzVector();
  m_SN_lowestMassWpos     = TLorentzVector();  
  m_SN_lowestMassWneg     = TLorentzVector();  
}

void TopDileptonReconstruction::Reconstruct(TLorentzVector lepton_pos, 
                                            TLorentzVector lepton_neg, 
                                            TLorentzVector b, 
                                            TLorentzVector bbar, 
                                            double met_ex, 
                                            double met_ey, 
                                            double mtop,
                                            double mtbar,
                                            double mWpos,
                                            double mWneg){

  if(doNW){
    ReconstructNW(lepton_pos, lepton_neg, b, bbar, met_ex, met_ey, mtop, mtbar, mWpos, mWneg);
  }

  if(doEM){
    ReconstructEM(lepton_pos, lepton_neg, b, bbar, met_ex, met_ey, mtop, mtbar, mWpos, mWneg);
    if(m_EM_tops.size() > 0){
      m_EM_averageTop   = Average(m_EM_tops);
      m_EM_averageTbar  = Average(m_EM_tbars);
      m_EM_averageNu    = Average(m_EM_nus);
      m_EM_averageNubar = Average(m_EM_nubars);
      m_EM_averageWpos  = Average(m_EM_Wposs);
      m_EM_averageWneg  = Average(m_EM_Wnegs);      
    }

  }

  if(doSN){
    ReconstructSN(lepton_pos, lepton_neg, b, bbar, met_ex, met_ey, mtop, mtbar, mWpos, mWneg);
    if(m_SN_tops.size() > 0){
      m_SN_averageTop   = Average(m_SN_tops);
      m_SN_averageTbar  = Average(m_SN_tbars);
      m_SN_averageNu    = Average(m_SN_nus);
      m_SN_averageNubar = Average(m_SN_nubars);
      m_SN_averageWpos  = Average(m_SN_Wposs);
      m_SN_averageWneg  = Average(m_SN_Wnegs);
    }
  }

}


TLorentzVector TopDileptonReconstruction::Average(std::vector<TLorentzVector> vecs){

  double average_m = 0;
  TVector3 three_vecs = TVector3();
  TLorentzVector blank;
  if(vecs.size() == 0){
    return blank;
  }

  for(uint i = 0; i < vecs.size(); ++i){
    TVector3 tmp_vec = TVector3();
    tmp_vec.SetXYZ(vecs[i].X(), vecs[i].Y(), vecs[i].Z());
    average_m += vecs[i].M();
    three_vecs = three_vecs + tmp_vec;
  }

  three_vecs = three_vecs*(1./vecs.size());
  average_m = average_m/vecs.size();

  TLorentzVector new_vec = TLorentzVector();
  new_vec.SetPtEtaPhiM(three_vecs.Pt(), three_vecs.Eta(), three_vecs.Phi(), average_m);
  return new_vec;
}


void TopDileptonReconstruction::ReconstructNW(TLorentzVector lepton_pos, 
                                              TLorentzVector lepton_neg, 
                                              TLorentzVector b, 
                                              TLorentzVector bbar, 
                                              double met_ex, 
                                              double met_ey, 
                                              double mtop,
                                              double mtbar,
                                              double mWpos,
                                              double mWneg){

  uint nsmears = 10;
  bool sol_nu_success = false;
  bool sol_nubar_success = false;
  std::vector<TLorentzVector> solution_nu;
  std::vector<TLorentzVector> solution_nubar;
  double nu_eta_mean    = lepton_pos.Eta()*0.72 - 0.018*lepton_pos.Eta()*lepton_pos.Eta()*lepton_pos.Eta();
  double nubar_eta_mean = lepton_neg.Eta()*0.72 - 0.018*lepton_neg.Eta()*lepton_neg.Eta()*lepton_neg.Eta();

  double width = 1.16;
  TRandom3 random = TRandom3(0);

  for(uint i = 0; i < nsmears; ++i){

    if (!sol_nu_success){
      double nu_eta = random.Gaus(nu_eta_mean, width);
      solution_nu = NW_solveForNeutrinoEta(&lepton_pos, &b, nu_eta, mtop, mWpos);
      if (solution_nu.size() > 0) {
        sol_nu_success = true ;
      }
    }

    if (!sol_nubar_success){
      double nubar_eta = random.Gaus(nubar_eta_mean, width);
      solution_nubar = NW_solveForNeutrinoEta(&lepton_neg, &bbar, nubar_eta, mtbar, mWneg);
      if (solution_nubar.size() > 0) {
        sol_nubar_success = true ;
      }
    }

    if( sol_nu_success && sol_nubar_success) break;
  }

  if(solution_nu.size()    == 0) return;
  if(solution_nubar.size() == 0) return;

  TLorentzVector top, tbar, ttbar, Wpos, Wneg, nu, nubar;

  for(uint index_nu = 0; index_nu < solution_nu.size(); ++index_nu){
    for(uint index_nubar = 0; index_nubar < solution_nubar.size(); ++index_nubar){

      nu    = solution_nu.at(index_nu);
      nubar = solution_nubar.at(index_nubar);
      Wpos  = lepton_pos + nu;
      Wneg  = lepton_neg + nubar;  
      top   = Wpos + b;
      tbar  = Wneg + bbar;

      m_NW_tops.push_back(top);
      m_NW_tbars.push_back(tbar);
      m_NW_nus.push_back(nu);
      m_NW_nubars.push_back(nubar);
      m_NW_Wposs.push_back(Wpos);
      m_NW_Wnegs.push_back(Wneg);

      double weight = NW_get_weight(nu, nubar, met_ex, met_ey);
      m_NW_weights.push_back(weight);
      if(weight > NWhighestWeight && weight > 0.00001 && top.E() > 0 && tbar.E() > 0){
      	NWhighestWeight = weight;
      	m_NW_highestWeightTop   = top;
      	m_NW_highestWeightTbar  = tbar;
      	m_NW_highestWeightNu    = nu;
      	m_NW_highestWeightNubar = nubar;
      	m_NW_highestWeightWpos  = Wpos;
      	m_NW_highestWeightWneg  = Wneg;
      }
    }
  }
  return;
}


void TopDileptonReconstruction::ReconstructEM(TLorentzVector lepton_pos, 
                                              TLorentzVector lepton_neg, 
                                              TLorentzVector b, 
                                              TLorentzVector bbar, 
                                              double met_ex, 
                                              double met_ey, 
                                              double mtop,
                                              double mtbar,
                                              double mWpos,
                                              double mWneg){

  double mNu    = 0;
  double mNubar = 0;

  std::vector<TMatrixD> nu_results = EM_getNeutrinoEllipse(b, lepton_pos, mtop,  mWpos, mNu);
  if(nu_results.size()==0) return;
  TMatrixD Nu_perp     = nu_results[0]; // Nperp
  TMatrixD H_nu        = nu_results[1]; // H
  TMatrixD Hperp_nu    = nu_results[2]; // Hperp
  TMatrixD HperpInv_nu = nu_results[3]; // HperpInv

  std::vector<TMatrixD> nubar_results = EM_getNeutrinoEllipse(bbar, lepton_neg, mtbar, mWneg, mNubar);  
  if(nubar_results.size()==0) return;
  TMatrixD Nubar_perp     = nubar_results[0];
  TMatrixD H_nubar        = nubar_results[1];
  TMatrixD Hperp_nubar    = nubar_results[2];
  TMatrixD HperpInv_nubar = nubar_results[3];

  if (Nu_perp == 0 || Nubar_perp == 0) return;

  // calculate second ellipse for nu
  TMatrixD Gamma(3,3);
  Gamma.Zero();
  Gamma[0][0] = -1;
  Gamma[1][1] = -1;
  Gamma[2][2] =  1;
  Gamma[0][2] = met_ex;
  Gamma[1][2] = met_ey;

  // second ellipse matrix for nu
  TMatrixD n_perp = TMatrixD(Gamma, TMatrixD::kTransposeMult, TMatrixD(Nubar_perp, TMatrixD::kMult, Gamma));

  std::vector<TVectorD> nu_perp = EM_intersect_ell_ell(Nu_perp, n_perp);
  std::vector<TVectorD> nubar_perp;

  if (nu_perp.size() == 0){

    // calculate n_perp ellipse eccentricity
    TMatrixD n_perp22(2,2);
    n_perp22[0][0] = n_perp[0][0];
    n_perp22[0][1] = n_perp[0][1];
    n_perp22[1][0] = n_perp[1][0];
    n_perp22[1][1] = n_perp[1][1];

    TMatrixDEigen n_p22_eig(n_perp22);

    TVectorD eval_re = n_p22_eig.GetEigenValuesRe();

    //Double_t ecc = sqrt(1. - eval_re[1]/eval_re[0]);
    // end - calulate n_perp ellipse eccentricity

    TMatrixD Xp = TMatrixD(Hperp_nu,TMatrixD::kTransposeMult,TMatrixD(n_perp, TMatrixD::kMult, Hperp_nu));

    TMatrixD D(3,3);
    D.Zero();
    D[0][1] = -1;
    D[1][0] =  1;

    TMatrixD XpD = TMatrixD(Xp, TMatrixD::kMult, D);
    TMatrixD Mp  = TMatrixD(TMatrixD::kTransposed, XpD) + XpD;

    TMatrixD U(3,3);
    U.Zero();
    U[0][0] =  1;
    U[1][1] =  1;
    U[2][2] = -1;

    std::vector<TVectorD> t_sols = EM_intersect_ell_ell(Mp, U);
    std::vector<std::pair<Double_t,TVectorD>> kt_sols;

    if (t_sols.size() == 0) return;

    for (uint i = 0; i < t_sols.size(); i++){
      TVectorD tX(3);
      for (int j=0; j<3; j++){
        TVectorD Xp_col = TMatrixDColumn(Xp,j);
        tX[j] = t_sols.at(i)*Xp_col;
      }
      Double_t k = fabs(tX*t_sols.at(i));
      kt_sols.push_back(std::make_pair(k, t_sols.at(i)));
    }
    std::sort(kt_sols.begin(), kt_sols.end(), [this] (std::pair<Double_t,TVectorD> kv1, std::pair<Double_t,TVectorD> kv2) {return EM_cmp(kv1, kv2);});
    //sort(kt_sols.begin(), kt_sols.end(), EM_cmp);
    nu_perp.push_back(Hperp_nu*kt_sols.at(0).second);

    TMatrixD D2 = TMatrixD(D,TMatrixD::kMult,D);
    TVectorD Np_nup = Nu_perp*nu_perp.at(0);
    TVectorD D2_Np_nup = D2*Np_nup;

    TVectorD line_perp(3);
    line_perp[0] = nu_perp.at(0)[1]*D2_Np_nup[2] - nu_perp.at(0)[2]*D2_Np_nup[1];
    line_perp[1] = nu_perp.at(0)[2]*D2_Np_nup[0] - nu_perp.at(0)[0]*D2_Np_nup[2];
    line_perp[2] = nu_perp.at(0)[0]*D2_Np_nup[1] - nu_perp.at(0)[1]*D2_Np_nup[0];

    std::vector<Double_t> kv;
    std::vector<TVectorD> nu_perp_p_tmp = EM_intersect_ell_line(n_perp, line_perp, kv);

    Double_t dist = 1e9;
    Int_t idist = -999;

    if (nu_perp_p_tmp.size() == 0) return;

    for (uint i=0; i<nu_perp_p_tmp.size(); i++){
      Double_t d = sqrt(pow(nu_perp_p_tmp.at(i)[0] - nu_perp.at(0)[0],2) + pow(nu_perp_p_tmp.at(i)[1] - nu_perp.at(0)[1],2));
      if (d < dist) {dist = d; idist = i;}
    }
    if (idist == -999) return;

    //if (ecc > 0.99) return;

    TVectorD nubar_perp_tmp = Gamma*nu_perp_p_tmp.at(idist);
    nubar_perp.push_back(nubar_perp_tmp);

    } else{
    for (uint i = 0; i < nu_perp.size(); i++){
      nubar_perp.push_back(Gamma*nu_perp.at(i));
    }
  }

  // H*HperpInv
  TMatrixD S_nu    = TMatrixD(H_nu,    TMatrixD::kMult, HperpInv_nu);
  TMatrixD S_nubar = TMatrixD(H_nubar, TMatrixD::kMult, HperpInv_nubar);

  assert (nu_perp.size() == nubar_perp.size());

  for (uint i=0; i < nu_perp.size(); i++){

    TVectorD p_nu_tmp    = S_nu    * nu_perp.at(i);
    TVectorD p_nubar_tmp = S_nubar * nubar_perp.at(i);

    TLorentzVector nu;
    TLorentzVector nubar;
    TLorentzVector top;
    TLorentzVector tbar;
    TLorentzVector Wpos;
    TLorentzVector Wneg;    

    nu.SetPxPyPzE(   p_nu_tmp[0],    p_nu_tmp[1],    p_nu_tmp[2],    sqrt(p_nu_tmp[0]*p_nu_tmp[0]       + p_nu_tmp[1]*p_nu_tmp[1]       + p_nu_tmp[2]*p_nu_tmp[2]));
    nubar.SetPxPyPzE(p_nubar_tmp[0], p_nubar_tmp[1], p_nubar_tmp[2], sqrt(p_nubar_tmp[0]*p_nubar_tmp[0] + p_nubar_tmp[1]*p_nubar_tmp[1] + p_nubar_tmp[2]*p_nubar_tmp[2]));    

    top  = nu    + lepton_pos + b;
    tbar = nubar + lepton_neg + bbar;

    Wpos = lepton_pos + nu;
    Wneg = lepton_neg + nubar;    

    m_EM_nus.push_back(nu);
    m_EM_nubars.push_back(nubar);
    m_EM_tops.push_back(top);
    m_EM_tbars.push_back(tbar);
    m_EM_Wposs.push_back(Wpos);
    m_EM_Wnegs.push_back(Wneg);  

    TLorentzVector ttbar;
    ttbar = top + tbar;

    double mttbar = ttbar.M();

    if(mttbar < EMlowestMass){
      EMlowestMass = mttbar;
      m_EM_lowestMassTop   = top;
      m_EM_lowestMassTbar  = tbar;
      m_EM_lowestMassWneg  = Wneg;
      m_EM_lowestMassWpos  = Wpos;
      m_EM_lowestMassNu    = nu;
      m_EM_lowestMassNubar = nubar;
    }
  }

  return;
}


void TopDileptonReconstruction::ReconstructSN(TLorentzVector lepton_pos, 
                                              TLorentzVector lepton_neg,
                                              TLorentzVector b,
                                              TLorentzVector bbar,
                                              double met_ex,
                                              double met_ey,
                                              double mtop,
                                              double mtbar,
                                              double mWpos,
                                              double mWneg){

  TVector3 bjet_v    = b.Vect();
  TVector3 lep_v     = lepton_neg.Vect();
  TVector3 bjetbar_v = bbar.Vect();
  TVector3 lepbar_v  = lepton_pos.Vect();

  //  hope I got these the right way around...
  double mW  = mWpos;
  double mW_ = mWneg;  

  //Accessing TLorentzVectors is a bit slow. Replace with doubles for frequently accessed components
  double b_E     = b.E();
  double b_Px    = b.Px();
  double b_Py    = b.Py(); 
  double b_Pz    = b.Pz();
  double mb      = b.M();   

  double bbar_E  = bbar.E();
  double bbar_Px = bbar.Px();
  double bbar_Py = bbar.Py(); 
  double bbar_Pz = bbar.Pz();  
  double mbbar   = bbar.M();   

  double lpos_E  = lepton_pos.E();
  double lpos_Px = lepton_pos.Px();
  double lpos_Py = lepton_pos.Py();
  double lpos_Pz = lepton_pos.Pz();
  double lpos_M  = lepton_pos.M();
  double lpos_M2 = lpos_M * lpos_M;
  double lpos_E2 = lpos_E * lpos_E;

  double lneg_E  = lepton_neg.E();
  double lneg_Px = lepton_neg.Px();
  double lneg_Py = lepton_neg.Py();
  double lneg_Pz = lepton_neg.Pz();
  double lneg_M  = lepton_neg.M();
  double lneg_M2 = lneg_M * lneg_M;
  double lneg_E2 = lneg_E * lneg_E;

  double mnu     = 0.;
  double mnubar  = 0.;

  double mW2     = mW     * mW;
  double mW2_    = mW_    * mW_;

  double mb2     = mb     * mb;
  double mbbar2  = mbbar  * mbbar;

  double mtop2   = mtop   * mtop;
  double mtbar2  = mtbar  * mtbar; 

  double mnu2    = mnu    * mnu;
  double mnubar2 = mnubar * mnubar;


  Double_t a1 = (b_E    + lpos_E) * (mW2  - lpos_M2 - mnu2)    - lpos_E * (mtop2  - mb2    - lpos_M2 - mnu2)    + 2. * b_E    * lpos_E2 - 2. * lpos_E * (bjet_v.Dot(lepbar_v));
  Double_t b1 = (bbar_E + lneg_E) * (mW2_ - lneg_M2 - mnubar2) - lneg_E * (mtbar2 - mbbar2 - lneg_M2 - mnubar2) + 2. * bbar_E * lneg_E2 - 2. * lneg_E * (bjetbar_v.Dot(lep_v));
                                                                    
  Double_t a2 = 2. * (b_E    * lpos_Px - lpos_E * b_Px);
  Double_t a3 = 2. * (b_E    * lpos_Py - lpos_E * b_Py);
  Double_t a4 = 2. * (b_E    * lpos_Pz - lpos_E * b_Pz);
  Double_t b2 = 2. * (bbar_E * lneg_Px - lneg_E * bbar_Px);
  Double_t b3 = 2. * (bbar_E * lneg_Py - lneg_E * bbar_Py);
  Double_t b4 = 2. * (bbar_E * lneg_Pz - lneg_E * bbar_Pz);


  Double_t c22  = (mW2  - lpos_M2 - mnu2)    * (mW2  - lpos_M2 - mnu2)    - 4. * (lpos_E2 - lpos_Pz * lpos_Pz) * (a1 / a4) * (a1 / a4) - 4. * (mW2  - lpos_M2 - mnu2)    * lpos_Pz * a1 / a4;
  Double_t dp22 = (mW2_ - lneg_M2 - mnubar2) * (mW2_ - lneg_M2 - mnubar2) - 4. * (lneg_E2 - lneg_Pz * lneg_Pz) * (b1 / b4) * (b1 / b4) - 4. * (mW2_ - lneg_M2 - mnubar2) * lneg_Pz * b1 / b4;

  Double_t c21  = 4. * (mW2  - lpos_M2 - mnu2)    * (lpos_Px - lpos_Pz * a2 / a4) - 8. * (lpos_E2 - lpos_Pz * lpos_Pz) * a1 * a2 / (a4 * a4) - 8. * lpos_Px * lpos_Pz * a1 / a4;
  Double_t dp21 = 4. * (mW2_ - lneg_M2 - mnubar2) * (lneg_Px - lneg_Pz * b2 / b4) - 8. * (lneg_E2 - lneg_Pz * lneg_Pz) * b1 * b2 / (b4 * b4) - 8. * lneg_Px * lneg_Pz * b1 / b4;
  Double_t c11  = 4. * (mW2  - lpos_M2 - mnu2)    * (lpos_Py - lpos_Pz * a3 / a4) - 8. * (lpos_E2 - lpos_Pz * lpos_Pz) * a1 * a3 / (a4 * a4) - 8. * lpos_Py * lpos_Pz * a1 / a4;
  Double_t dp11 = 4. * (mW2_ - lneg_M2 - mnubar2) * (lneg_Py - lneg_Pz * b3 / b4) - 8. * (lneg_E2 - lneg_Pz * lneg_Pz) * b1 * b3 / (b4 * b4) - 8. * lneg_Py * lneg_Pz * b1 / b4;

  Double_t c20  = -4. * (lpos_E2 - lpos_Px * lpos_Px) - 4. * (lpos_E2 - lpos_Pz * lpos_Pz) * (a2 / a4) * (a2 / a4) - 8. * lpos_Px * lpos_Pz * a2 / a4;
  Double_t dp20 = -4. * (lneg_E2 - lneg_Px * lneg_Px) - 4. * (lneg_E2 - lneg_Pz * lneg_Pz) * (b2 / b4) * (b2 / b4) - 8. * lneg_Px * lneg_Pz * b2 / b4;
  Double_t c00  = -4. * (lpos_E2 - lpos_Py * lpos_Py) - 4. * (lpos_E2 - lpos_Pz * lpos_Pz) * (a3 / a4) * (a3 / a4) - 8. * lpos_Py * lpos_Pz * a3 / a4;
  Double_t dp00 = -4. * (lneg_E2 - lneg_Py * lneg_Py) - 4. * (lneg_E2 - lneg_Pz * lneg_Pz) * (b3 / b4) * (b3 / b4) - 8. * lneg_Py * lneg_Pz * b3 / b4;

  Double_t c10  = -8. * (lpos_E2 - lpos_Pz * lpos_Pz) * a2 * a3 / (a4 * a4) + 8. * lpos_Px * lpos_Py - 8. * lpos_Px * lpos_Pz * a3 / a4 - 8. * lpos_Py * lpos_Pz * a2 / a4;
  Double_t dp10 = -8. * (lneg_E2 - lneg_Pz * lneg_Pz) * b2 * b3 / (b4 * b4) + 8. * lneg_Px * lneg_Py - 8. * lneg_Px * lneg_Pz * b3 / b4 - 8. * lneg_Py * lneg_Pz * b2 / b4;

  c22 = c22 * a4 * a4;
  c21 = c21 * a4 * a4;
  c20 = c20 * a4 * a4;
  c11 = c11 * a4 * a4;
  c10 = c10 * a4 * a4;
  c00 = c00 * a4 * a4;

  Double_t d22 = dp22 + met_ex * met_ex * dp20 + met_ey * met_ey * dp00 + met_ex * met_ey * dp10 + met_ex * dp21 + met_ey * dp11;
  Double_t d21 = -1. * dp21 - 2.* met_ex * dp20 - met_ey * dp10;
  Double_t d20 = dp20;
  Double_t d11 = -1. * dp11 - 2. * met_ey * dp00 - met_ex * dp10;
  Double_t d10 = dp10;
  Double_t d00 = dp00;

  d22 = d22 * b4 * b4;
  d21 = d21 * b4 * b4;
  d20 = d20 * b4 * b4;
  d11 = d11 * b4 * b4;
  d10 = d10 * b4 * b4;
  d00 = d00 * b4 * b4;


  Double_t h4 = c00 * c00 * d22 * d22 + c11 * d22 * (c11 * d00 - c00 * d11) + c00 * c22 * (d11 * d11 - 2. * d00 * d22) + c22 * d00 * (c22 * d00 - c11 * d11);
  Double_t h3 = c00 * d21 * (2. * c00 * d22 - c11 * d11) + c00 * d11 * (2. * c22 * d10 + c21 * d11) + c22 * d00 * (2. * c21 * d00 -      c11 * d10) -      c00 * d22 * (c11 * d10 + c10 * d11) - 2. * c00 * d00 * (     c22 * d21 + c21 * d22) - d00 * d11 * (c11 * c21 + c10 * c22) + c11 * d00 *( c11 * d21 + 2 * c10 * d22);
  Double_t h2 = c00 * c00 * (2. * d22 * d20 + d21 * d21) - c00 * d21 * (     c11 * d10 + c10 * d11) + c11 * d20 * (     c11 * d00 -      c00 * d11) +      c00 * d10 * (c22 * d10 - c10 * d22) +      c00 * d11 * (2. * c21 * d10 + c20 * d11) + ( 2. * c22 * c20 + c21 * c21) * d00 * d00 - 2. * c00 * d00 * (c22 * d20 + c21 * d21 + c20 * d22) + c10 * d00 * (2. * c11 * d21 + c10 * d22) - d00 * d10 * (c11 * c21 + c10 * c22) -d00 * d11 * (c11 * c20 + c10 * c21);
  Double_t h1 = c00 * d21 * (2. * c00 * d20 - c10 * d10) - c00 * d20 * (     c11 * d10 + c10 * d11) + c00 * d10 * (     c21 * d10 + 2. * c20 * d11) - 2. * c00 * d00 * (c21 * d20 + c20 * d21) +      c10 * d00 * (2. * c11 * d20 + c10 * d21) + c20 * d00 * (2 * c21 * d00 - c10 * d11) - d00 * d10 * (c11 * c20 + c10 * c21);
  Double_t h0 = c00 * c00 * d20 * d20 + c10 * d20 * (c10 * d00 - c00 * d10) + c20 * d10 * (c00 * d10 - c10 * d00) + c20 * d00 * (c20 * d00 - 2. * c00 * d20);

  ROOT::Math::Polynomial poly = ROOT::Math::Polynomial(h0,h1,h2,h3,h4);
  std::vector<Double_t> roots = poly.FindRealRoots();

  if (roots.size() == 0) return;

  // nu and nubar momenta
  std::vector<TVector3> p_nu;
  std::vector<TVector3> p_nubar;

  Double_t p_nu_x, p_nu_y, p_nu_z;
  Double_t p_nubar_x, p_nubar_y, p_nubar_z;

  for (std::vector<Double_t>::iterator it = roots.begin() ; it!=roots.end() ; ++it){

    p_nu_x = *it;

    // calculate p_nu_y
    Double_t c0 = c00;
    Double_t c1 = c11 + c10 * p_nu_x;
    Double_t c2 = c22 + c21 * p_nu_x + c20 * p_nu_x * p_nu_x;

    Double_t d0 = d00;
    Double_t d1 = d11 + d10 * p_nu_x;
    Double_t d2 = d22 + d21 * p_nu_x + d20 * p_nu_x * p_nu_x;

    p_nu_y = (c0 * d2 - c2 * d0)/(c1 * d0 - c0 * d1);

    // calculate p_nu_z
    p_nu_z = (-1. * a1 - a2 * p_nu_x - a3 * p_nu_y) / a4;

    // calculate nubar momentum components
    p_nubar_x = met_ex - p_nu_x;
    p_nubar_y = met_ey - p_nu_y;
    p_nubar_z = (-1. * b1 - b2 * p_nubar_x - b3 * p_nubar_y) / b4;

    TVector3 p_nu_xyz(   p_nu_x,    p_nu_y,    p_nu_z);
    TVector3 p_nubar_xyz(p_nubar_x, p_nubar_y, p_nubar_z);

    TLorentzVector nu, nubar;
    nu.SetPtEtaPhiM(p_nu_xyz.Pt(),       p_nu_xyz.Eta(),    p_nu_xyz.Phi(),    0.);
    nubar.SetPtEtaPhiM(p_nubar_xyz.Pt(), p_nubar_xyz.Eta(), p_nubar_xyz.Phi(), 0.);  

    TLorentzVector top, tbar, Wpos, Wneg;
    Wpos = lepton_pos + nu;
    Wneg = lepton_neg + nubar;
    top  = Wpos + b;
    tbar = Wneg + bbar;

    m_SN_tops.push_back(top);
    m_SN_tbars.push_back(tbar);
    m_SN_nus.push_back(nu);
    m_SN_nubars.push_back(nubar); 
    m_SN_Wposs.push_back(Wpos);
    m_SN_Wnegs.push_back(Wneg);

    TLorentzVector ttbar;
    ttbar = top + tbar;

    double mttbar = ttbar.M();

    if(mttbar < SNlowestMass){
      SNlowestMass = mttbar;
      m_SN_lowestMassTop   = top;
      m_SN_lowestMassTbar  = tbar;
      m_SN_lowestMassWneg  = Wneg;
      m_SN_lowestMassWpos  = Wpos;
      m_SN_lowestMassNu    = nu;
      m_SN_lowestMassNubar = nubar;
    }
  }
  return;
}


std::vector<TLorentzVector> TopDileptonReconstruction::NW_solveForNeutrinoEta(TLorentzVector* lepton, 
                                                                              TLorentzVector* bJet, 
                                                                              double nu_eta, 
                                                                              double mtop, 
                                                                              double mW) {

  double nu_cosh = cosh(nu_eta);
  double nu_sinh = sinh(nu_eta);

  double Wmass2  = mW*mW;
  double bmass   = 4.5;
  double Elprime = lepton->E() * nu_cosh - lepton->Pz() * nu_sinh;
  double Ebprime = bJet->E()   * nu_cosh - bJet->Pz()   * nu_sinh;

  double A = (lepton->Py() * Ebprime - bJet->Py() * Elprime) / (bJet->Px() * Elprime - lepton->Px() * Ebprime);
  double B = (Elprime * (mtop * mtop - Wmass2 - bmass * bmass - 2. * lepton->Dot(*bJet)) - Ebprime * Wmass2) / (2. * (lepton->Px() * Ebprime - bJet->Px() * Elprime));

  double par1 = (lepton->Px() * A + lepton->Py()) / Elprime;
  double C    = A * A + 1. - par1 * par1;
  double par2 = (Wmass2 / 2. + lepton->Px() * B) / Elprime;
  double D    = 2. * (A * B - par2 * par1);
  double F    = B * B - par2 * par2;
  double det  = D * D - 4. * C * F;

  std::vector<TLorentzVector> sol;

  ///-- 0 solutions case --///
  if (det < 0.0){
    return std::move(sol);
  }

  ///-- Only one real solution case --///
  if (det == 0.) {
    double py1 = -D / (2. * C);
    double px1 = A * py1 + B;
    double pT2_1 = px1 * px1 + py1 * py1;
    double pz1 = sqrt(pT2_1) * nu_sinh;

    TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));

    if (!TMath::IsNaN(a1.E()) )
      sol.push_back(a1);
    return std::move(sol);
  }

  ///-- 2 solutions case --///
  if(det > 0){
    double tmp   = sqrt(det) / (2. * C);
    double py1   = -D / (2. * C) + tmp;
    double py2   = -D / (2. * C) - tmp;
    double px1   = A * py1 + B;
    double px2   = A * py2 + B;
    double pT2_1 = px1 * px1 + py1 * py1;
    double pT2_2 = px2 * px2 + py2 * py2;
    double pz1   = sqrt(pT2_1) * nu_sinh;
    double pz2   = sqrt(pT2_2) * nu_sinh;
    TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));
    TLorentzVector a2(px2, py2, pz2, sqrt(pT2_2 + pz2 * pz2));

    if (!TMath::IsNaN(a1.E()) && !TMath::IsNaN(a2.E())){
      sol.push_back(a1);
      sol.push_back(a2);
    }
    return std::move(sol);
  }
  
  ///-- Should never reach this point --///
  return std::move(sol);
}


double TopDileptonReconstruction::NW_get_weight(TLorentzVector nu1, TLorentzVector nu2, double met_ex, double met_ey){

    /*double dx = met_ex - (nu1 + nu2).Px();
    double dy = met_ey - (nu2 + nu2).Py();
    double m_sigma_met_ex  = 2*met_ex;
    double m_sigma_met_ey  = 2*met_ey;

    double numerator_x = -dx*dx;
    double numerator_y = -dy*dy;

    double denominator_x = m_sigma_met_ex*m_sigma_met_ex;
    double denominator_y = m_sigma_met_ey*m_sigma_met_ey;

    double exp_x = exp(numerator_x/denominator_x);
    double exp_y = exp(numerator_y/denominator_y);

    return exp_x*exp_y;*/

    double met = sqrt(met_ex*met_ex + met_ey*met_ey);
    double dt  = met - (nu1 + nu2).Pt();
    double sigma = 2*met;
    double numerator = -dt*dt;
    double denominator = sigma*sigma;
    return exp(numerator/denominator);
}


std::vector<TMatrixD> TopDileptonReconstruction::EM_getNeutrinoEllipse(TLorentzVector& bjet, TLorentzVector& lepton, double& mtop, double& mW, double& mNu){
  
  double mt2  = mtop * mtop;
  double mW2  = mW   * mW;
  double mnu2 = mNu  * mNu;

  std::vector<TMatrixD> empty_vectors;

  Double_t bjet_beta     = bjet.Beta();
  Double_t bjet_beta2    = bjet_beta * bjet_beta;

  Double_t lepton_beta   = lepton.Beta();
  Double_t lepton_beta2  = lepton_beta * lepton_beta;

  Double_t cos = ROOT::Math::VectorUtil::CosTheta(lepton,bjet);
  Double_t s   = sqrt(1. - cos * cos);
  
  TMatrixD Ab(4,4);
  TMatrixD Al(4,4);
  TMatrixD Htilde(3,3);
  TMatrixD H(3,3);
  TMatrixD Hperp(3,3);
  TMatrixD HperpInv(3,3);
  TMatrixD Nperp(3,3);
  TMatrixD neutrino_ellipse(3,3); 

  Ab.Zero();
  Al.Zero();
  Htilde.Zero();
  H.Zero();
  Hperp.Zero();
  HperpInv.Zero();
  Nperp.Zero();
  neutrino_ellipse.Zero();

  //- W surface 
  Double_t x0p = -(0.5 / bjet.E())   * (mt2 - mW2  - bjet.M2());
  Double_t x0  = -(0.5 / lepton.E()) * (mW2 - mnu2 - lepton.M2());
  Double_t Sx  = (1. / lepton_beta2) * (x0 * lepton_beta - lepton.P() * (1.- lepton_beta2));
  Double_t epsilon2 = (1. - lepton_beta2) * ( mW2 - mnu2);

  //- lepton ellipsoid
  Al[0][0] = 1. - lepton_beta2;
  Al[1][0] = 0;
  Al[2][0] = 0;
  Al[3][0] = Sx * lepton_beta2;
    
  Al[0][1] = 0;
  Al[1][1] = 1;
  Al[2][1] = 0;
  Al[3][1] = 0;
    
  Al[0][2] = 0;
  Al[1][2] = 0;
  Al[2][2] = 1;
  Al[3][2] = 0;
    
  Al[0][3] = Sx * lepton_beta2;
  Al[1][3] = 0;
  Al[2][3] = 0;
  Al[3][3] = mW2 - x0 * x0 - epsilon2;

  //- bjet ellipsoid
  Ab[0][0] = 1 - cos * cos * bjet_beta2;
  Ab[1][0] = -cos * s * bjet_beta2;
  Ab[2][0] = 0;
  Ab[3][0] = cos * x0p * bjet_beta;
    
  Ab[0][1] = -cos * s * bjet_beta2;
  Ab[1][1] = 1 - s * s * bjet_beta2;
  Ab[2][1] = 0;
  Ab[3][1] = s * x0p * bjet_beta;
    
  Ab[0][2] = 0;
  Ab[1][2] = 0;
  Ab[2][2] = 1;
  Ab[3][2] = 0;
    
  Ab[0][3] = cos * x0p * bjet_beta;
  Ab[1][3] = s   * x0p * bjet_beta;
  Ab[2][3] = 0;
  Ab[3][3] = mW2 - x0p * x0p;

  //Neutrino solution
  Double_t Sy = (1. / s) * (x0p / bjet_beta - cos * Sx);
  Double_t omega = (1. / s) * (lepton_beta / bjet_beta - cos); //only the positive slope
  Double_t Omega = sqrt(std::max(0., omega * omega + 1. - lepton_beta2));
  double Omega2 = Omega * Omega;
  Double_t x1 = Sx - (1. / Omega2) * (Sx + omega * Sy);
  Double_t y1 = Sy - (1. / Omega2) * omega * (Sx + omega * Sy);
  Double_t Z2 = x1 * x1 * Omega2 - (Sy - omega * Sx) * (Sy - omega * Sx) - (mW2 - x0 * x0 - epsilon2);
  double Z = sqrt(std::max(0., Z2));

  Htilde[0][0] = Z / Omega;
  Htilde[0][1] = 0;
  Htilde[0][2] = x1 - lepton.P();
  
  Htilde[1][0] = omega * Z / Omega;
  Htilde[1][1] = 0;
  Htilde[1][2] = y1;
    
  Htilde[2][0] = 0;
  Htilde[2][1] = Z;
  Htilde[2][2] = 0;

  //lab system transform, rotate Htilde to H
  TMatrixD R(3,3);
  R.Zero();

  TMatrixD Rz = EM_rotationMatrix(2, - lepton.Phi());
  TMatrixD Ry = EM_rotationMatrix(1, 0.5 * M_PI - lepton.Theta());

  double bJetP[3] = {bjet.Px(), bjet.Py(), bjet.Pz()};
  TMatrixD bJet_xyz(3, 1, bJetP);
  TMatrixD rM(Ry, TMatrixD::kMult, TMatrixD(Rz, TMatrixD::kMult,bJet_xyz));
  double* rA = rM.GetMatrixArray();
  double phi = -TMath::ATan2(rA[2], rA[1]);
  TMatrixD Rx = EM_rotationMatrix(0, phi);
  R = TMatrixD(Rz, TMatrixD::kTransposeMult, TMatrixD(Ry, TMatrixD::kTransposeMult, Rx.T()));
  H = TMatrixD(R, TMatrixD::kMult, Htilde);

  //calculate Hperp
  double Hvalues[9]={H[0][0], H[0][1], H[0][2], H[1][0], H[1][1], H[1][2], 0, 0, 1};
  TArrayD Harray(9, Hvalues);
  Hperp.SetMatrixArray(Harray.GetArray());

  //calculate Nperp
  double Hperp_Det = Hperp[0][0] * Hperp[1][1] * Hperp[2][2] - 
                     Hperp[0][0] * Hperp[1][2] * Hperp[2][1] -
                     Hperp[0][1] * Hperp[1][0] * Hperp[2][2] +
                     Hperp[0][1] * Hperp[1][2] * Hperp[2][0] +
                     Hperp[0][2] * Hperp[1][0] * Hperp[2][1] -
                     Hperp[0][2] * Hperp[1][1] * Hperp[2][0];

  if (Hperp_Det == 0) return empty_vectors;

  HperpInv = Hperp;
  HperpInv.Invert();

  TMatrixD U(3,3);
  U.Zero();
  U[0][0] =  1;
  U[1][1] =  1;
  U[2][2] = -1;
  Nperp = TMatrixD(HperpInv, TMatrixD::kTransposeMult, TMatrixD(U, TMatrixD::kMult, HperpInv));

  std::vector<TMatrixD> return_vectors;
  return_vectors.push_back(Nperp);
  return_vectors.push_back(H);
  return_vectors.push_back(Hperp);
  return_vectors.push_back(HperpInv);

  return return_vectors;
}


std::vector<TVectorD> TopDileptonReconstruction::EM_intersect_ell_ell(TMatrixD A, TMatrixD B){

    TMatrixD TMP = TMatrixD(A);

    Double_t DetA(0.), DetB(0.);

    for (int j=0; j<3; j++){
        DetA += A[0][j] * EM_cofactor(A,0,j);
        DetB += B[0][j] * EM_cofactor(B,0,j);
    }

    if (fabs(DetB) > fabs(DetA)) {A=B; B=TMP;}

    TMatrixD AinvB = TMatrixD(A,TMatrixD::kInvMult,B);

    TMatrixDEigen EIG(AinvB);

    TVectorD eigval_re = EIG.GetEigenValuesRe();
    TVectorD eigval_im = EIG.GetEigenValuesIm();

    int i_eig = 0;

    for (int i=0; i<3; i++){
        if (eigval_im[i]!=0) continue;
        i_eig = i;
    }

    std::vector<TVectorD> lines = EM_factor_degenerate(B - eigval_re[i_eig]*A);

    std::vector<TVectorD> points;
    std::vector<Double_t> kvect;

    for (uint i=0; i<lines.size(); i++){
        TVectorD line = lines.at(i);
	std::vector<Double_t> kv_tmp;
	std::vector<TVectorD> p_tmp = EM_intersect_ell_line(A,line,kv_tmp);
        points.insert(points.end(),p_tmp.begin(),p_tmp.end());
        kvect.insert(kvect.end(),kv_tmp.begin(),kv_tmp.end());
    }

    return points;
}

bool TopDileptonReconstruction::EM_cmp(std::pair<Double_t,TVectorD> kv1, std::pair<Double_t,TVectorD> kv2){
    return kv1.first < kv2.first;
}


std::vector<TVectorD> TopDileptonReconstruction::EM_intersect_ell_line(TMatrixD E, TVectorD L, std::vector<Double_t> &kv){

    //cross product L X E
    TMatrixD LXE(3,3);
    for (int j = 0; j < 3; j++){
        LXE[0][j] = L[1]*E[2][j] - L[2]*E[1][j];
        LXE[1][j] = L[2]*E[0][j] - L[0]*E[2][j];
        LXE[2][j] = L[0]*E[1][j] - L[1]*E[0][j];
    }

    TMatrixDEigen EIG(LXE);

    TMatrixD eig_v = EIG.GetEigenVectors();

    std::vector<TVectorD> sols;
    std::vector<std::pair<Double_t,TVectorD>> ksols;

    for (int j=0; j<3; j++){
        TVectorD v = TMatrixDColumn(eig_v,j);;

        v*=(1./sqrt(v.Norm2Sqr()));
        Double_t lv = L*v;
        TVectorD ve(3);
        for (int jj=0; jj<3; jj++){
            TVectorD E_col = TMatrixDColumn(E,jj);
            ve[jj] = v*E_col;
        }
        Double_t vev = ve*v;
        Double_t k = pow(lv,2) + pow(vev,2);

        v*=(1./v[2]);

        ksols.push_back(std::make_pair(k,v));
    }

    if (ksols.size() < 2) return sols;

    sort(ksols.begin(), ksols.end(), [this] (std::pair<Double_t,TVectorD> kv1, std::pair<Double_t,TVectorD> kv2) {return EM_cmp(kv1, kv2);});
    for (int i=0; i<2; i++){
        if (ksols.at(i).first > 1e-12) continue;
        sols.push_back(ksols.at(i).second);
        kv.push_back(ksols.at(i).first);
    }
    return sols;
}


Double_t TopDileptonReconstruction::EM_cofactor(TMatrixD A, int row, int col){

    int di,dj;
    TMatrixD M(2,2);

    di = 0;
    for (int i=0; i<2; i++){
        if (i == row) di = 1;
        dj = 0;
        for (int j=0; j<2; j++){
            if (j == col) dj = 1;
            M[i][j] = A[i+di][j+dj];
        }
    }

    return pow(-1,row+col) * (M[0][0]*M[1][1] - M[1][0]*M[0][1]);
}


TMatrixD TopDileptonReconstruction::EM_rotationMatrix(int axis, double angle){
 
  TMatrixD r(3,3);
  r.Zero();
  if (axis!=0 && axis!=1 && axis!=2) return r;
  
  for( int i=0; i<3; i++ ){
    r[i][i]=cos(angle);
  }

  for( int i=-1; i<=1; i++ ){
    double row  = (axis - i)%3; if(row<0) row+=3;
    double col  = (axis + i)%3; if(col<0) col+=3;
    r[row][col] = i * sin(angle) + (1 - i*i);
  }

  //r.Print();

  return r;
}


std::vector<TVectorD> TopDileptonReconstruction::EM_factor_degenerate(TMatrixD G){

    TVectorD Lp(3);
    TVectorD Lm(3);

    if (G[0][0] == 0. && G[1][1] == 0.){
        Lp[0] = G[0][1]; Lp[1] = 0; Lp[2] = G[1][2];
        Lm[0] = 0; Lm[1] = G[0][1]; Lm[2] = G[0][2] - G[1][2];
        return std::vector<TVectorD> {Lp,Lm};
    }

    bool swapXY = (fabs(G[0][0]) > fabs(G[1][1]));

    if (swapXY) {
        TMatrixD TMP = G;
        G[0][1] = TMP[1][0];
        G[1][0] = TMP[0][1];
        G[0][0] = TMP[1][1];
        G[1][1] = TMP[0][0];
        G[0][2] = TMP[1][2];
        G[2][0] = TMP[2][1];
        G[1][2] = TMP[0][2];
        G[2][1] = TMP[2][0];
    }

    //G*=(1./G[1][1]); ???
    Double_t g22 = EM_cofactor(G,2,2);

    if (g22 == 0. && G[1][1] !=0.){
        Double_t g00 = EM_cofactor(G,0,0);
        if (g00 > 0.) {return (std::vector<TVectorD>) 0;}

        Lp[0] = G[0][1]; Lp[1] = G[1][1]; Lp[2] = G[1][2] + sqrt(-g00);
        Lm[0] = G[0][1]; Lm[1] = G[1][1]; Lm[2] = G[1][2] - sqrt(-g00);
    } else {
      if(g22 > 0.) {return (std::vector<TVectorD>) 0;}

        Double_t x0 = EM_cofactor(G,0,2)/g22;
        Double_t y0 = EM_cofactor(G,1,2)/g22;

        Lp[0] = G[0][1] + sqrt(-g22); Lp[1] = G[1][1]; Lp[2] = -G[1][1]*y0 - (G[0][1] + sqrt(-g22))*x0;
        Lm[0] = G[0][1] - sqrt(-g22); Lm[1] = G[1][1]; Lm[2] = -G[1][1]*y0 - (G[0][1] - sqrt(-g22))*x0;

        TVectorD LpTMP = Lp;
        TVectorD LmTMP = Lm;

        Lp[0] = LpTMP[(int)swapXY]; Lp[1] = LpTMP[(int)!swapXY]; Lp[2] = LpTMP[2];
        Lm[0] = LmTMP[(int)swapXY]; Lm[1] = LmTMP[(int)!swapXY]; Lm[2] = LmTMP[2];

    }

    return std::vector<TVectorD> {Lp,Lm};
}

//} // End of namespace

