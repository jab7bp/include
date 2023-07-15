#ifndef GMNA_H
#define GMNA_H

//////////////////////////////
////Static Detector Parameters
//Tracks
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
//Trigger TDC
const Int_t maxTDCTrigChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Double_t tdiffwidecut = 50; // Set at 50ns from nominal 510ns (passed by user) from GMn
//HCal - Note that actual measurement of vertical length is 381.6cm, indicating that the MC figures are correct
const Int_t maxHCalChan = 288; // Total HCal channels
const Int_t maxHCalRows = 24; // Total HCal rows
const Int_t maxHCalCols = 12; // Total HCal cols
const Int_t maxClusters = 10; // Total HCal clusters with information saved
const Double_t HCALHeight = 0.365; // Height of HCal above beamline in m
const Double_t HCalblk_l = 0.15; // Width and height of HCAL blocks in m
const Double_t HCalblk_l_h_MC = 0.15494; // Width and height of all HCAL blocks in m from MC database
const Double_t HCalblk_l_v_MC = 0.15875; // Width and height of all HCAL blocks in m from MC database
//const Double_t posHCalXi = -2.20; // Distance from beam center to top of HCal in m
//const Double_t posHCalXf = 1.47; // Distance from beam center to bottom of HCal in m
//const Double_t posHCalYi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
//const Double_t posHCalYf = 0.853; // Distance from beam center to beam side of HCal in m
const Double_t posHCalXi = -2.165; // Distance from beam center to top of HCal in m
const Double_t posHCalXf = 1.435; // Distance from beam center to bottom of HCal in m
const Double_t posHCalYi = -0.9; // Distance from beam center to opposite-beam side of HCal in m
const Double_t posHCalYf = 0.9; // Distance from beam center to beam side of HCal in m
const Double_t posHCalXi_MC = -2.355005; // Distance from beam center to top of HCal in m from MC database
const Double_t posHCalXf_MC = 1.454995; // Distance from beam center to bottom of HCal in m from MC database
const Double_t posHCalYi_MC = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
const Double_t posHCalYf_MC = 0.92964; // Distance from beam center to beam side of HCal in m from MC database
const Double_t HCalSampFrac = 0.077;  //Re-evaluated with MC GEn settings using second to outermost shower column for kin2

//BBCal
const Int_t maxBBCalShChan = 189; // Total BBCal Shower Channels
const Int_t maxBBCalShRows = 27;
const Int_t maxBBCalShCols = 7;
const Int_t maxBBCalPSChan = 52; // Total BBCal Preshower Channels
const Int_t maxBBCalPSRows = 26;
const Int_t maxBBCalPSCols = 2;
//Beamline
const Int_t chargeConvert = 3318; // See D.Flay Doc DB sbs.jlab.org/DocDB/0001/000164/002/dflay_bcm-ana-update_02-21-22.pdf p.8
const Int_t clockActual = 103700; // Needed to convert the 104kHz clock to the actual counting rate
//SBS Magnet
const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
const Double_t maxSBSfield = 1.26; //Tesla
const Double_t SBSfield = 1.0; //fraction of max field. TODO: should be variable per run
const Double_t SBSdist = 2.25; //m
const Double_t dipGap = 1.22; //m
const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;
//GEMs
const Double_t GEMpitch = 10*TMath::DegToRad();

///////////////
///Physics/Math
// const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;
const UInt_t us = 1000000; //For conversion to seconds used by reporting time delays

////////////////////////////
////Static Target/Scattering Chamber Parameters
const Double_t l_tgt = 0.15; // Length of the target (m)
const Double_t rho_tgt = 0.0723; // Density of target (g/cc)
const Double_t rho_Al = 2.7; // Density of aluminum windows (g/cc)
const Double_t celldiameter = 1.6*2.54; //cm, right now this is a guess
const Double_t Ztgt = 1.0;
const Double_t Atgt = 1.0;
const Double_t Mmol_tgt = 1.008; //g/mol
const Double_t dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const Double_t dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const Double_t uwallthick_LH2 = 0.0145; //cm
const Double_t dwallthick_LH2 = 0.015; //cm
const Double_t cellthick_LH2 = 0.02; //cm, this is a guess;
const Double_t Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch 
const Double_t tdiffmax = 20; //ns. Maximum allowed deviation from HCal trig = BBCal trig = t_trig (from config file)

//////////////////////////////
////BG fit functions
//BG Gaussian Fit
const Double_t sbs4mag0prise = -1.2;
const Double_t sbs4mag0nfall = 0.48;
const Double_t sbs4mag30prise = -1.2;
const Double_t sbs4mag30nfall = 0.48;
const Double_t sbs4mag50prise = -1.77;
const Double_t sbs4mag50nfall = 0.51;

Double_t g4BGFit0(double *x, double *par){
  Double_t amp = par[0];
  Double_t loc = par[1];
  Double_t sigma = par[2];
  Double_t xpos = x[0];

  if(xpos>sbs4mag0prise && xpos<sbs4mag0nfall) { 
    TF1::RejectPoint();
    return 0;
  }

  return amp * TMath::Exp(-0.5 * pow(((xpos - loc) / sigma), 2.));
}

Double_t g4BGFit30(double *x, double *par){
  Double_t amp = par[0];
  Double_t loc = par[1];
  Double_t sigma = par[2];
  Double_t xpos = x[0];

  if(xpos>sbs4mag30prise && xpos<sbs4mag30nfall) { 
    TF1::RejectPoint();
    return 0;
  }

  return amp * TMath::Exp(-0.5 * pow(((xpos - loc) / sigma), 2.));
}

Double_t g4BGFit50(double *x, double *par){
  Double_t amp = par[0];
  Double_t loc = par[1];
  Double_t sigma = par[2];
  Double_t xpos = x[0];

  if(xpos>sbs4mag50prise && xpos<sbs4mag50nfall){
    TF1::RejectPoint();
    return 0;
  }

  return amp * TMath::Exp(-0.5 * pow(((xpos - loc) / sigma), 2.));
}

#endif // GMNA_H