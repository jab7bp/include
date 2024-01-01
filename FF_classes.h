#ifndef FF_CLASSES_H
#define FF_CLASSES_H

class kine_ff_extract {
	public:

//Experimental Kinematics
	int kine;
	int sbsfieldscale;
	TString run_target;

	double Ebeam;

	double BB_theta_rad;
	double BB_theta_deg;
	double HCal_theta_rad;
	double HCal_theta_deg;

//Calibrations
	double dx_p_mean;
	double dx_p_sigma;
	double dx_n_mean;
	double dx_n_sigma;
	double dy_mean;
	double dy_sigma;


//Data histograms
	TH1D *h_data_dx; //Final dx plot to be sent to BGsub, etc. PRE-BGSUB, etc. (includes wcut, fcut, etc.)
	TH1D *h_data_dy; //Similar to h_data_dx but for 'dy'
	TH1D *h_data_Q2; //Momentum-transfer, Q^2, from reconstructed kinematic variables

	TH2D *h_data_dxdy; //histogram with wcut and fcut that is used for ProjectionY() to get the h_data_dx plot


//simc histograms
	TH1D *h_simc_dx; //Scaled and chi-square minimized sum of h_simc_dx_p and h_simc_dx_n
	TH1D *h_simc_dx_p; //Scaled and chi-square minimized simc proton data;
	TH1D *h_simc_dx_n; //Scaled anc chi-square minimized simc neutron data;


//Fit Functions


//DATA and SIMC Comparisons - yields, scale factors, etc. 
	double yield_n;  //Integral of SCALED (and chi-square min.) h_simc_dx_n
	double yield_p;  //Integral of SCALED (and chi-square min.) h_simc_dx_p

	double scale_n;  //Scale factor used to chi-square min h_simc_dx_n and data
	double scale_p;  //Scale factor used to chi-square min h_simc_dx_p and data
	double scale_np_ratio;  //Ratio of scale factors



//Form Factors - Theoretical
	double GEp;
	double GEp_err;

	double GMp;
	double GMp_err;

	double GEn;
	double GEn_err;

	double GMn_theory;
	double GMN_theory_err;

	//Experimental Form Factors

	double GMn_exp;
	double GMn_exp_error;


};

class cl_SBSkine {
	public:


//Experimental values:
	int kine, sbsfieldscale;
	TString run_target;

//Data experimental values:
	double Q2_data_mean;
	double dx_p_mean, dx_p_sigma, dx_n_mean, dx_n_sigma, dy_mean, dy_sigma;

//Analysis/extraction values:
	double dy_mult;
	double scale_p, scale_n, scale_np_ratio;
	double minuit_chiSquared, minuit_chiSquared_NDF;

	//BOOLEANS
	bool use_inelastic_BGsub, use_dy_anticut_BGsub;

	//Fits:
	TF1 *tf_dx_total;
	TF1 *tf_dyAnticut_BG, *tf_inelastics_BG;
	TF1 *tf_dyAnticut_BG_scaled, *tf_inelastics_BG_scaled;

//Analysis histograms:
	//DATA
	TH1D *h_data_dx, *h_data_dy;
	TH1D *h_data_dx_from_dyAnticut;
	TH1D *h_data_dyAnticut_BG_scaled, *h_data_inelastics_BG_scaled;
	TH1D *h_data_dx_with_dyAnticut_BGsub, *h_data_dx_with_inelastics_BGsub;

	TH2D *h2_data_dxdy;
	TH2D *h2_data_dxdy_with_dyAnticut;

	//simc
	TH1D *h_simc_dx, *h_simc_dx_p, *h_simc_dx_n;
	TH1D *h_simc_dx_p_BeforeScaling, *h_simc_dx_n_BeforeScaling;
	TH1D *h_mc_p_sigma, *h_mc_n_sigma, *h_mc_np_sigma_ratio;
		//BG
	TH1D *h_simc_dx_inelastics_BG, *h_simc_dx_inelastics_BG_scaled;

//Form Factor input variables:

//Form Factor calculations:

//Form Factor outputs:
	double GMn_final;
	double GMn_norm_final;

};

#endif