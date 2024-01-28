#ifndef SYST_CLASSES_H
#define SYST_CLASSES_H

// class syst {
class cl_syst {
	public:

	int kine;

	double sbsfieldscale_Frac;
	TString systematic_varVal_str = "";
	TString systematic_var_str = "";

	double val_min, val_max;

	double systematic_varVal = 0.0;
	double w2_mult_min, w2_mult_max, invMass_min, invMass_max;
	double SHPS_sigma_mult_min, SHPS_sigma_mult_max;
	double Ep_min;
	double PS_min, HCal_clus_e_min;
	int gem_track_nhits_min;
	double fcut_mult;
	double fcut_mult_min, fcut_mult_max;
};


#endif