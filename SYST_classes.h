#ifndef SYST_CLASSES_H
#define SYST_CLASSES_H

// class syst {
class cl_syst {
	public:

	int kine;
	int pass;

	double sbsfieldscale_Frac;
	TString systematic_varVal_str = "";
	TString systematic_var_str = "";

	double val_min, val_max;

	double systematic_varVal = 0.0;
	double w2_mult_min, w2_mult_max, invMass_min, invMass_max;
	double SHPS_sigma_mult_min, SHPS_sigma_mult_max;
	double SH_PS_min, SH_PS_max;
	double Ep_min, Ep_max;
	double PS_min, PS_max, HCal_clus_e_min, HCal_clus_e_max;
	int gem_track_nhits_min;
	double fcut_mult;
	double fcut_mult_min, fcut_mult_max;
	double Ep_sigma_mult_min, Ep_sigma_mult_max;
	double HCal_XY_multiplier = 1.0;
	double ADC_time_mult = 1.0;
	double ADC_diff_time_mult = 1.0;
	double Ep_sig_mult = 4.0;
	double HCal_nblocks = 1;
};


#endif