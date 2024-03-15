#ifndef FF_CLASSES_H
#define FF_CLASSES_H

// class kine_ff_extract {
class cl_SBSkine {
	public:


//Experimental Kinematics
	int kine;
	int sbsfieldscale;
	int pass;
	TString run_target;
	TString magmod;

	double E_beam;

	double BB_theta_rad;
	double BB_theta_deg;
	double HCal_theta_rad;
	double HCal_theta_deg;

	double Q2_calculated;
	double Q2_data_mean;
	double Q2_data_std_dev;

	double e_scat_momentum;

//Booleans
	bool extracted = false;

	bool use_dy_anticut_BGsub; //Background source: dy_anticut used for BGsub
	bool use_inelastic_BGsub; //Background source: fit to simc inelastics BG (always !use_dy_anticut_BGsub)
	bool use_collected; //Use some histograms collected from others' analysis?

//Calibrations
	double dx_p_mean;
	double dx_p_sigma;
	double dx_n_mean;
	double dx_n_sigma;
	double dy_mean;
	double dy_sigma;

	double dy_mult; //Multiplier used on the dy_sigma to define the dy_anticut size/region
	double dy_mult_minus, dy_mult_plus; // Used for assymmetrical dy cut

//Data histograms
	TH1D *h_data_dx; //Final dx plot to be sent to BGsub, etc. PRE-BGSUB, etc. (includes wcut, fcut, etc.)
	TH1D *h_data_dy; //Similar to h_data_dx but for 'dy'
	TH1D *h_data_Q2; //Momentum-transfer, Q^2, from reconstructed kinematic variables

	TH2D *h2_data_dxdy; //histogram with wcut and fcut that is used for ProjectionY() to get the h_data_dx plot


	//dyAnticut 
	TH1D *h_data_dx_from_dyAnticut; //dx projection from the 2D dxdy with a dyAnticut
	TH1D *h_data_dyAnticut_BG_scaled; //h_data_dx_from_dyAnticut SCALED to match original total_dx BG
	TH1D *h_data_dx_with_dyAnticut_BGsub; //dx plot AFTER dyAnticut BGsub

	TH2D *h2_data_dxdy_with_dyAnticut; //2D dxdy with the dyAnticut applied

	//inelastics
	TH1D *h_simc_dx_inelastics_BG;  //histogram made from fit to simc inelastics BG
	TH1D *h_simc_dx_inelastics_BG_scaled; //histo from simc inelastics scaled to total_dx BG
	TH1D *h_data_dx_with_inelastics_BGsub; //dx plot AFTER inelastics BGsub

	////
	TH1D *h_data_dx_with_select_BGsub; //SELECTED (dyAnticut or inelastics) dx plot WITH BGsub
	Int_t N_data_dx_with_select_BGsub; //Number of entries in BG sub (selected) data dx
	TH1D *h_data_dx_with_select_BGsub_DoubleBinWidth; 
	TH1D *h_data_dx_with_select_BGsub_DoubleBinWidth_CloneForPlotting; 
	TH1D *h_simc_dx_scaled_to_BGsub;
	TH1D *h_simc_dx_scaled_to_BGsub_CloneForPlotting;
	TH1D *h_simc_dx_scaled_to_BGsub_DoubleBinWidth;
	TH1D *h_simc_dx_scaled_to_BGsub_DoubleBinWidth_CloneForPlotting;

//simc histograms
	TH1D *h_simc_dx; //Scaled and chi-square minimized sum of h_simc_dx_p and h_simc_dx_n
	TH1D *h_simc_dx_p; //Scaled and chi-square minimized simc proton data;
	TH1D *h_simc_dx_n; //Scaled anc chi-square minimized simc neutron data;
	TH1D *h_simc_dx_p_DoubleBinWidth, *h_simc_dx_n_DoubleBinWidth;
	TH1D *h_simc_dx_n_BeforeScaling; //raw h_simc_dx_n, saved from original and un-modified
	TH1D *h_simc_dx_p_BeforeScaling; //raw h_simc_dx_p, saved from original and un-modified

	TH1D *h_mc_p_sigma; //Weighted cross-section for proton -- taken from the simc output tree
	TH1D *h_mc_n_sigma; //Weighted cross-section for neutron -- taken from the simc output tree
	TH1D *h_mc_np_sigma_ratio; //Ratio of the weighted cross-sections from simc

//Residual
	TH1D *h_data_simc_dx_residual;

//Fit Functions
	TF1 *tf_dx_total; //Total fit function for the data: p_gaus + n_gaus + BGpol
	
	//dyAnticut
	TF1 *tf_dyAnticut_BG; //Fit to the dx projection of the 2D dxdy with dyAnticut
	TF1 *tf_dyAnticut_BG_scaled; //Fit to the dx proj dyAnticut SCALED to mach BGpol from tf_dx_total

	//inelastics
	TF1 *tf_inelastics_BG; //Fit to the simc inelastics BG
	TF1 *tf_inelastics_BG_scaled; //Fit to the simc inelastics BG SCALED to total_dx BG



//------------------------------------------------------------------------
//------------------------------------------------------------------------
//DATA and SIMC Comparisons - yields, scale factors, etc. 
	double yield_n;  //Integral of SCALED (and chi-square min.) h_simc_dx_n
	double yield_n_error;
	double yield_n_error_Nnorm;
	double yield_p;  //Integral of SCALED (and chi-square min.) h_simc_dx_p
	double yield_p_error;  //Integral of SCALED (and chi-square min.) h_simc_dx_p
	double yield_p_error_Nnorm;
	double yield_np_ratio; //Ratio of SCALED yields
	double yield_np_ratio_error;
	double yield_total;

	double scale_n;  //Scale factor used to chi-square min h_simc_dx_n and data
	double scale_n_error;
	double scale_n_error_Nnorm;
	double scale_p;  //Scale factor used to chi-square min h_simc_dx_p and data
	double scale_p_error;
	double scale_p_error_Nnorm;
	double scale_np_ratio;  //Ratio of scale factors
	double scale_np_ratio_error;

	double data_simc_yield_ratio; //ratio of N in simc_dx p + n yields and Data_dx_with_select_BGsub

	double minuit_chiSquared; //Minimum chi-squared found in TMinuit minimization
	double minuit_chiSquared_NDF; //Redcued chi-squared: chiSquared/NDF (NDF = 500 - 2);


//Form Factor and Cross-section variables
	double eps_n = 0.0, eps_p = 0.0, tau_n = 0.0, tau_p = 0.0, G_D = 0.0;
	double mott_CS;
	double sigma_R_p_param = 0.0;
	double sigma_R_p_param_err = 0.0;

	double sigma_R_n_param = 0.0;
	double sigma_R_n_param_err = 0.0;

//Form Factors - Theoretical
	//general variables to hold any parameterization
	double GEp_param = 0.0;
	double GEp_param_err = 0.0;

	//-----------

	double GMp_param = 0.0;
	double GMp_param_err = 0.0;

	//-----------

	double GEn_param = 0.0;
	double GEn_param_err = 0.0;

	double GEn_param_norm = 0.0;
	double GEn_param_norm_err = 0.0;

	//-----------

	double GMn_param = 0.0;
	double GMn_param_err = 0.0;

////////////////////
	//bosted
	double GEp_bosted = 0.0;
	double GEp_bosted_err = 0.0;

	//-----------

	double GMp_bosted = 0.0;
	double GMp_bosted_err = 0.0;

	//-----------

	double GEn_bosted = 0.0;
	double GEn_bosted_err = 0.0;

	double GEn_bosted_norm = 0.0;
	double GEn_bosted_norm_err = 0.0;

	//-----------

	double GMn_bosted = 0.0;
	double GMn_bosted_err = 0.0;

	//-----------

	//Experimental Form Factors

	double GMn_exper = 0.0;
	double GMn_exper_err = 0.0;

	double GMn_exper_norm = 0.0;
	double GMn_exper_norm_err = 0.0;


};

class cl_nTPE {
	public:

	//Experimental Kinematics
	int kine;
	int sbsfieldscale;
	int pass;
	TString run_target;
	TString magmod;

	double E_beam;

	double BB_theta_rad;
	double BB_theta_deg;
	double HCal_theta_rad;
	double HCal_theta_deg;

	double Q2_calculated;
	double Q2_data_mean;

	double e_scat_momentum;

	double SBS8_yield_n, SBS8_yield_n_error, SBS8_yield_p, SBS8_yield_p_error;
	double SBS8_yield_ratio, SBS8_yield_ratio_error;
	double SBS8_scale_n, SBS8_scale_n_error, SBS8_scale_p, SBS8_scale_p_error;
	double SBS8_scale_ratio, SBS8_scale_ratio_error;
	double SBS8_GMn, SBS8_GMn_error, SBS8_GEn_, SBS8_GEn_error;
	double SBS8_sigma_R, SBS8_sigma_R_error;
	double SBS8_tau_n, SBS8_tau_p, SBS8_eps_n, SBS8_eps_p;
	double SBS8_RS_n, SBS8_RS_n_error, SBS8_RS_p, SBS8_RS_p_error;
	double SBS8_sigma_T, SBS8_sigma_T_error, SBS8_sigma_L, SBS8_sigma_L_error;

	double SBS9_yield_n, SBS9_yield_n_error, SBS9_yield_p, SBS9_yield_p_error;
	double SBS9_yield_ratio, SBS9_yield_ratio_error;
	double SBS9_scale_n, SBS9_scale_n_error, SBS9_scale_p, SBS9_scale_p_error;
	double SBS9_scale_ratio, SBS9_scale_ratio_error;
	double SBS9_GMn, SBS9_GMn_error, SBS9_GEn_, SBS9_GEn_error;
	double SBS9_sigma_R, SBS9_sigma_R_error;
	double SBS9_tau_n, SBS9_tau_p, SBS9_eps_n, SBS9_eps_p;
	double SBS9_RS_n, SBS9_RS_n_error, SBS9_RS_p, SBS9_RS_p_error;
	double SBS9_sigma_T, SBS9_sigma_T_error, SBS9_sigma_L, SBS9_sigma_L_error;


};

#endif