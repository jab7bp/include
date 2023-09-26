#ifndef FF_CLASSES_H
#define FF_CLASSES_H

class kine_ff_extract {
	public:

	TH1D *h_data;
	TH1D *h_simc;

	double p_yield;
	double n_yield;
	double R;
	double Rcorr;

	double p_scale;
	double n_scale;
	double scale_R;
	double scale_Rcorr;

	double p_eff;
	double p_eff_err;
	double n_eff;
	double n_eff_err;

	double GEn_param;
	double GEp_param;
	double GMp_param;
	double GMn_param;
	double GMn;

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