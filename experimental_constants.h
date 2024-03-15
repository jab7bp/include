#ifndef EXPERIMENTAL_CONSTANTS_H
#define EXPERIMENTAL_CONSTANTS_H


double E_beam; //Electron beam energy (electron energy) in GeV
double BB_dist;
double BB_theta;
double HCal_dist;
double HCal_theta;

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
const double MN = (Mp + Mn)/2.0;
const double Me = 0.00051; //Mass of electron [GeV]

const double alpha = 1.0/137.0;
const double mu_p = 2.79;
const double mu_n = 1.91;
const double delta = 0.710;

//SBS Magnet
const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
const Double_t maxSBSfield = 1.26; //Tesla
const Double_t SBSdist = 2.25; //m
const Double_t dipGap = 1.22; //m
const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;

//Static Detector Parameters
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
// const double hcal_height = -0.2897; // Height of HCal above beamline
double hcal_height;

const Double_t sampfrac = 0.077; 	//Estimate of the sampling fraction from MC
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable max number of tracks per event
const Int_t kNtdc = 1000; // Reasonable max number of tdc signals per event
const Int_t max_clus = 10;

//Values have ben updated -- Known as of Oct 1, 2023......
const Double_t Xi = -2.655; //Distance from beam center to top of HCal in meters, from database
const Double_t Xf = 1.155; //Distance from beam center to bottom of HCal in meters, from database
const Double_t Yi = -0.92964; //Distance from beam center to opposite-beam side of HCal in meters, from MC database
const Double_t Yf = 0.92964; //Distance from beam center to beam side of HCal in meters, from MC database

const Double_t HCalblk_l_h = 0.15494; //Horizontal length of all HCal blocks in meters, from MC database
const Double_t HCalblk_l_v = 0.15875; //Vertical length of all HCal blocks in meters, from MC database
double fiducial_xi, fiducial_xf, fiducial_yi, fiducial_yf;
double hcal_x_fmin, hcal_x_fmax, hcal_y_fmin, hcal_y_fmax;

//PRevious values:
// const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
// const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
// const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
// const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m

//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double cell_diameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0;
const double Atgt = 1.0;
const double Mmol_tgt = 1.008; //g/mol

//For energy-loss correction to beam energy:
const double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const double uwallthick_LH2 = 0.0145; //cm
const double dwallthick_LH2 = 0.015; //cm
const double cellthick_LH2 = 0.02; //cm, this is a guess;
const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch


#endif