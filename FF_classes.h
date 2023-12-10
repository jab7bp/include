#ifndef FF_CLASSES_H
#define FF_CLASSES_H

class kine_ff_extract {
	public:

	TString run_target;
	int kine, sbsfieldscale;

	TH1D *h_data;
	TH1D *h_data_BeforeScaling;
	TH1D *h_data_dx_dyAntiCut;
	TH1D *h_data_dx_BGsub;
	TH1D *h_data_dx_BGsub_selected;
	TH1D *h_data_dx_dyAntiCut_BG_scaled;

	TH1D *h_simu_dx_inelastic_BGhisto;
	TH1D *h_simu_dx_inelastic_BGhisto_scaled;
	TH1D *h_simu_dx_inelastic_BGsub;

	TH2D *h_data_dxdy;
	TH2D *h_data_dxdy_dyAntiCut;
	TH2D *h_data_dxdy_dyCut;

	TH1D *h_simc;
	TH1D *h_simc_dx_p;
	TH1D *h_simc_dx_n;
	TH1D *h_simc_dx_BeforeScaling;
	TH1D *h_simc_dx_p_BeforeScaling;
	TH1D *h_simc_dx_n_BeforeScaling;

	TH1D *h_simc_dx_p_BGsub, *h_simc_dx_n_BGsub, *h_simc_dx_BGsub;

	double dx_p_mean, dx_p_sigma, dx_n_mean, dx_n_sigma, dy_mean, dy_sigma;
	TF1 *tf_data_dx_total, *tf_data_dx_pinit, *tf_data_dx_ninit, *tf_data_dx_BGinit, *tf_data_dx_dyAntiCut_BG;
	TF1 *tf_data_dx_dyAntiCut_BG_scaled, *tf_inelasticBG, *tf_inelasticBG_scaled;
	TF1 *tf_data_dx_BGsub_p, *tf_data_dx_BGsub_n, *tf_data_dx_BGsub_total;

	double data_dx_BGsub_total_integral;
	double data_dx_BGsub_p_fit_integral, data_dx_BGsub_n_fit_integral, data_dx_BGsub_total_fit_integral;

	double data_dx_total_par[11];
	int dx_dyAnticut_BG_polN;

	double dy_mult;
	double antiCut_BGscale_factor;

	double simc_dx_total_integral;
	double p_yield;
	double n_yield;
	double R;
	double Rcorr;

	double BGsub_min_chiSquared, BGsub_min_chiSquared_ndf, BGsub_min_chiSquared_err, BGsub_min_chiSquared_ndf_err;

	//To calculate the approximate ratio of n/p in the data histogram we can form the relation between the total integral and n/p ratios
	//       ( total data dx integral )	            ( total sim dx integral)
	//	-------------------------------------  = -------------------------------
	// ( R_data = data_n_yield/data_p_yield )    ( R = simc_n_yield/simc_p_yield)

	// --> (I_data/R_data) = (MC_total_integral/R);

	Double_t R_sim; //R_sim = (MC total integral)/(R)
	Double_t R_data_calc; //R_data = (I_data)*(R_sim)
	Double_t data_p_yield_calc;
	Double_t data_n_yield_calc;

	double p_scale;
	double n_scale;
	double scale_R;
	double scale_Rcorr;
	double np_scale_ratio;
	double prescale_ratio;

	//---- BG sub ratio values
	double BGsub_simc_n_yield;
	double BGsub_simc_n_yield_error;
	double BGsub_simc_p_yield;
	double BGsub_simc_p_yield_error;
	double BGsub_np_ratio;
	double BGsub_np_ratio_error;
	double BGsub_np_ratio_corr;
	double BGsub_np_ratio_corr_error;
	double BGsub_eff_ratio_error;

	double BGsub_np_scale_factor_ratio;
	double BGsub_np_scale_ratio;

	double BGsub_global_scale_p;
	double BGsub_global_scale_p_error;
	double BGsub_global_scale_n;
	double BGsub_global_scale_n_error;
	double BGsub_np_ratio_scaled_error;
	double BGsub_np_ratio_scaled_corr_error;
	double BGsub_simc_p_scaled_yield_error;
	double BGsub_simc_n_scaled_yield_error;

	double p_eff;
	double p_eff_err;
	double n_eff;
	double n_eff_err;

	double simc_p_yield_error;
	double simc_n_yield_error;
	double eff_ratio_error;
	double np_ratio_error;
	double np_ratio_corr_error;
	double np_scale_ratio_error;
	double global_scale_p_error;
	double global_scale_n_error;

	double GEn_param;
	double GEp_param;
	double GMp_param;
	double GMn_param;
	double GMn;
	double GMn_norm;

	double tau_p;
	double tau_n;
	double eps_p;
	double eps_n;

	double theta_deg;
	double theta_rad;
	double Q2_mean;
	double Ebeam_mean;
	double Enucl_mean;
	double Ee_scat_mean;

	double nTPE;

	double sigma_r_p;
	double G_D;

////---------- SAME VARIABLES BUT USING Q2, EBEAM, ETC FROM DATA HISTOGRAMS
	double p_yield_data;
	double n_yield_data;
	double R_data;
	double Rcorr_data;

	double p_scale_data;
	double n_scale_data;
	double scale_R_data;
	double scale_Rcorr_data;

	double p_eff_data;
	double p_eff_err_data;
	double n_eff_data;
	double n_eff_err_data;

	double GEn_param_data;
	double GEp_param_data;
	double GMp_param_data;
	double GMn_param_data;
	double GMn_data;

	double tau_p_data;
	double tau_n_data;
	double eps_p_data;
	double eps_n_data;

	double theta_deg_data;
	double theta_rad_data;
	double Q2_mean_data;
	double Ebeam_mean_data;
	double Enucl_mean_data;
	double Ee_scat_mean_data;

	double nTPE_data;

	double sigma_r_p_data;

	double G_D_data;

};

#endif