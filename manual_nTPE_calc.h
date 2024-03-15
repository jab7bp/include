#ifndef MANUAL_NTPE_CALC_H
#define MANUAL_NTPE_CALC_H

#include "/work/halla/sbs/jboyd/include/FF_classes.h"
#include "/work/halla/sbs/jboyd/include/experimental_constants.h"
#include "/work/halla/sbs/jboyd/include/beam_variables.h"
#include "/work/halla/sbs/jboyd/include/calc_errors.h"
#include "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/extract_GMn/include/calc_FF_bosted.h"
#include "/work/halla/sbs/jboyd/analysis/gmn/world_data/parameterizations/kelly_param_only.h"
#include "/work/halla/sbs/jboyd/analysis/gmn/world_data/jboyd_data_points.h"

bool use_theoretical_Q2 = false;
TString parameterization = "kelly";

cl_nTPE nTPE;

// General values/constants:
int sbsfieldscale = 70;
TString run_target = "LD2";

double SBS8_eps_err = 0.001, SBS8_sigma_R_err, SBS9_eps_err = 0.001, SBS9_sigma_R_err;

double Sn_TPE, R_sigma;
double Sn_approx;
double Sn_OPE, TPE;
double mu_GEn_GMn_proposal = 0.55;
double GEn_GMn_proposal;
double Sn_OPE_proposal = 0.126;
double GEn_average, GMn_average;
double Q_sq_average;
double tau_n_average, eps_n_average;
double GEn_OPE, GEn_OPE_error, GMn_OPE, GMn_OPE_error;
double GEn_GMn_ratio_OPE, GEn_GMn_ratio_OPE_error;

double SBS8_data_Q2, SBS8_data_Q2_err, SBS9_data_Q2, SBS9_data_Q2_err;

double R_corr_eps1, R_corr_eps2, A, B, Sn, Sp;
double RS_p_average;

//----------------------------------------
// +-+-+-+-+-+-+-+- SBS8 +-+-+-+-+-+-+-+-

// Experimental setup, SBS8:
int SBS8_kine = 8;
double SBS8_E_beam;
double SBS8_theta_e, SBS8_theta_e_deg;
double SBS8_e_scat_momentum, SBS8_Q_sq;
double SBS8_mott_CS_numerator = 0, SBS8_mott_CS_denominator = 0;
double SBS8_mott_CS, SBS8_G_D;

// Analaysis/Extraction values, SBS8:
double SBS8_yield_n, SBS8_yield_n_error, SBS8_yield_p, SBS8_yield_p_error;
double SBS8_scale_factor_n, SBS8_scale_factor_n_error, SBS8_scale_factor_p, SBS8_scale_factor_p_error;
double SBS8_GMn_calc, SBS8_GMn_calc_stat_error;
double SBS8_GMn_selected, SBS8_GMn_selected_syst_error, SBS8_GMn_selected_stat_error;
double SBS8_sigma_R_n_exp, SBS8_sigma_R_n_exp_error;
double SBS8_tau_n, SBS8_tau_p, SBS8_eps_n, SBS8_eps_p;
double SBS8_RS_n, SBS8_RS_n_error, SBS8_RS_p, SBS8_RS_p_error;
double SBS8_sigma_T_n, SBS8_sigma_T_n_error, SBS8_sigma_L_n, SBS8_sigma_L_n_error;
double SBS8_sigma_T_p, SBS8_sigma_T_p_error, SBS8_sigma_L_p, SBS8_sigma_L_p_error;

// Parameterizations, SBS8:
double SBS8_GEp_param, SBS8_GEp_param_error, SBS8_GMp_param, SBS8_GMp_param_error;
double SBS8_GEn_param, SBS8_GEn_param_error, SBS8_GMn_param, SBS8_GMn_param_error;
double SBS8_sigma_R_n_param, SBS8_sigma_R_n_param_error;

//----------------------------------------
// +-+-+-+-+-+-+-+- SBS9 +-+-+-+-+-+-+-+-

// Experimental setup, SBS9:
int SBS9_kine = 9;
double SBS9_E_beam;
double SBS9_theta_e, SBS9_theta_e_deg;
double SBS9_e_scat_momentum, SBS9_Q_sq;
double SBS9_mott_CS_numerator = 0, SBS9_mott_CS_denominator = 0;
double SBS9_mott_CS, SBS9_G_D;

// Analaysis/Extraction values, SBS9:
double SBS9_yield_n, SBS9_yield_n_error, SBS9_yield_p, SBS9_yield_p_error;
double SBS9_scale_factor_n, SBS9_scale_factor_n_error, SBS9_scale_factor_p, SBS9_scale_factor_p_error;
double SBS9_GMn, SBS9_GMn_syst_error, SBS9_GMn_stat_error, SBS9_GEn, SBS9_GEn_error;
double SBS9_GMn_calc, SBS9_GMn_calc_stat_error;
double SBS9_GMn_selected, SBS9_GMn_selected_syst_error, SBS9_GMn_selected_stat_error;
double SBS9_sigma_R_n_exp, SBS9_sigma_R_n_exp_error;
double SBS9_tau_n, SBS9_tau_p, SBS9_eps_n, SBS9_eps_p;
double SBS9_RS_n, SBS9_RS_n_error, SBS9_RS_p, SBS9_RS_p_error;
double SBS9_sigma_T_n, SBS9_sigma_T_n_error, SBS9_sigma_L_n, SBS9_sigma_L_n_error;
double SBS9_sigma_T_p, SBS9_sigma_T_p_error, SBS9_sigma_L_p, SBS9_sigma_L_p_error;

// Parameterizations, SBS9:
double SBS9_GEp_param, SBS9_GEp_param_error, SBS9_GMp_param, SBS9_GMp_param_error;
double SBS9_GEn_param, SBS9_GEn_param_error, SBS9_GMn_param, SBS9_GMn_param_error;
double SBS9_sigma_R_n_param, SBS9_sigma_R_n_param_error;

double SBS8_sigma_R_n_stat_err, SBS8_sigma_R_n_syst_err;
double SBS9_sigma_R_n_stat_err, SBS9_sigma_R_n_syst_err;

//----------------------------------------

void calc_nTPE(){

	// FILL THESE VALUES FROM THE ANALYSIS/EXTRACTION
	// MANUALLY FILL THIS IN FROM AVERAGES OF Y-PROJECTIONS AND STD. DEVs.
	
//SBS8
	SBS8_sigma_R_err = jboyd_GMn_SBS8_syst_error/100.0;
	// SBS8_sigma_R_err = sqrt( pow(SBS8_GMn_syst_error, 2) + pow(SBS8_GMn_stat_error, 2));

	SBS8_data_Q2 = 4.4071;
	SBS8_data_Q2_err = 0.3100;

//SBS9
	SBS9_sigma_R_err = jboyd_GMn_SBS9_syst_error/100.0;
	// SBS9_sigma_R_err = sqrt( pow(jboyd_GMn_SBS9_syst_error, 2) + pow(jboyd_GMn_SBS9_stat_error, 2));

	SBS9_data_Q2 = 4.3927; 
	SBS9_data_Q2_err = 0.1857;

/////////////
	//Choose beteween using theoretical Q2 or data-based

// We start with some general kinematic calculations for each kine, SBS8 and SBS9:

//SBS8
	SBS8_E_beam = lookup_beam_energy_from_kine(SBS8_kine);
	SBS8_theta_e = lookup_BB_angle_by_kine(SBS8_kine, "rad");
	SBS8_theta_e_deg = lookup_BB_angle_by_kine(SBS8_kine, "deg");

	SBS8_e_scat_momentum = (Mp*SBS8_E_beam)/(Mp + SBS8_E_beam*(1 - cos(SBS8_theta_e)));

	if( use_theoretical_Q2 ){
		SBS8_Q_sq = (2*MN*SBS8_E_beam*SBS8_E_beam*(1 - cos(SBS8_theta_e)))/(MN + SBS8_E_beam*(1 - cos(SBS8_theta_e)));
	}	
	else{
		SBS8_Q_sq = SBS8_data_Q2;
	}

	SBS8_mott_CS_numerator = (alpha*alpha*pow(cos( SBS8_theta_e/2), 2)*SBS8_e_scat_momentum);
	SBS8_mott_CS_denominator = 4*SBS8_E_beam*SBS8_E_beam*SBS8_E_beam*pow( sin(SBS8_theta_e/2), 4);
	SBS8_mott_CS = SBS8_mott_CS_numerator/SBS8_mott_CS_denominator;

	SBS8_G_D = pow( 1 + (SBS8_Q_sq/delta), -2);

	SBS8_tau_p = (SBS8_Q_sq)/(4*Mp*Mp);
	SBS8_tau_n = (SBS8_Q_sq)/(4*Mn*Mn);
	SBS8_eps_p = pow( 1.0 + 2.0*(1.0 + SBS8_tau_p)*pow( tan(SBS8_theta_e/2.0), 2 ), -1);
	SBS8_eps_n = pow( 1.0 + 2.0*(1.0 + SBS8_tau_n)*pow( tan(SBS8_theta_e/2.0), 2 ), -1);

//SBS9
	SBS9_E_beam = lookup_beam_energy_from_kine(SBS9_kine);
	SBS9_theta_e = lookup_BB_angle_by_kine(SBS9_kine, "rad");
	SBS9_theta_e_deg = lookup_BB_angle_by_kine(SBS9_kine, "deg");

	SBS9_e_scat_momentum = (Mp*SBS9_E_beam)/(Mp + SBS9_E_beam*(1 - cos(SBS9_theta_e)));

	if( use_theoretical_Q2 ){
		SBS9_Q_sq = (2*MN*SBS9_E_beam*SBS9_E_beam*(1 - cos(SBS9_theta_e)))/(MN + SBS9_E_beam*(1 - cos(SBS9_theta_e)));
	}
	else{
		SBS9_Q_sq = SBS9_data_Q2;
	}

	SBS9_mott_CS_numerator = (alpha*alpha*pow(cos( SBS9_theta_e/2), 2)*SBS9_e_scat_momentum);
	SBS9_mott_CS_denominator = 4*SBS9_E_beam*SBS9_E_beam*SBS9_E_beam*pow( sin(SBS9_theta_e/2), 4);
	SBS9_mott_CS = SBS9_mott_CS_numerator/SBS9_mott_CS_denominator;

	SBS9_G_D = pow( 1 + (SBS9_Q_sq/delta), -2);

	SBS9_tau_p = (SBS9_Q_sq)/(4*Mp*Mp);
	SBS9_tau_n = (SBS9_Q_sq)/(4*Mn*Mn);
	SBS9_eps_p = pow( 1.0 + 2.0*(1.0 + SBS9_tau_p)*pow( tan(SBS9_theta_e/2.0), 2 ), -1);
	SBS9_eps_n = pow( 1.0 + 2.0*(1.0 + SBS9_tau_n)*pow( tan(SBS9_theta_e/2.0), 2 ), -1);

///--------------
	tau_n_average = (SBS8_tau_n + SBS9_tau_n)/2.0;
	eps_n_average = (SBS8_eps_n + SBS9_eps_n)/2.0;

/////////////////////////////////////////////////////////////////////
//  PARAMETERIZATION CALCULATIONS

	if( parameterization == "bosted" ){
		//SBS8
		SBS8_GEp_param = calc_FF_bosted(1, SBS8_Q_sq, false);
		SBS8_GEp_param_error = calc_FF_bosted_error(1, SBS8_Q_sq, false);

		SBS8_GMp_param = calc_FF_bosted(2, SBS8_Q_sq, false);
		SBS8_GMp_param_error = calc_FF_bosted_error(2, SBS8_Q_sq, false);

		SBS8_GEn_param = calc_FF_bosted(3, SBS8_Q_sq, false);
		SBS8_GEn_param_error = calc_FF_bosted_error(3, SBS8_Q_sq, false);

		SBS8_GMn_param = calc_FF_bosted(4, SBS8_Q_sq, false);
		SBS8_GMn_param_error = calc_FF_bosted_error(4, SBS8_Q_sq, false);

		SBS8_sigma_R_n_param = SBS8_eps_n*(SBS8_GEn_param*SBS8_GEn_param) + SBS8_tau_n*(SBS8_GMn_param*SBS8_GMn_param);
		SBS8_sigma_R_n_param_error = calc_sigma_R_N_Err(SBS8_eps_n, SBS8_GEn_param, SBS8_GEn_param_error, SBS8_tau_n, SBS8_GMn_param, SBS8_GMn_param_error);
		

		//SBS9
		SBS9_GEp_param = calc_FF_bosted(1, SBS9_Q_sq, false);
		SBS9_GEp_param_error = calc_FF_bosted_error(1, SBS9_Q_sq, false);

		SBS9_GMp_param = calc_FF_bosted(2, SBS9_Q_sq, false);
		SBS9_GMp_param_error = calc_FF_bosted_error(2, SBS9_Q_sq, false);

		SBS9_GEn_param = calc_FF_bosted(3, SBS9_Q_sq, false);
		SBS9_GEn_param_error = calc_FF_bosted_error(3, SBS9_Q_sq, false);

		SBS9_GMn_param = calc_FF_bosted(4, SBS9_Q_sq, false);
		SBS9_GMn_param_error = calc_FF_bosted_error(4, SBS9_Q_sq, false);		

		SBS9_sigma_R_n_param = SBS9_eps_n*(SBS9_GEn_param*SBS9_GEn_param) + SBS9_tau_n*(SBS9_GMn_param*SBS9_GMn_param);
		SBS9_sigma_R_n_param_error = calc_sigma_R_N_Err(SBS9_eps_n, SBS9_GEn_param, SBS9_GEn_param_error, SBS9_tau_n, SBS9_GMn_param, SBS9_GMn_param_error);
		
	}
	else{ //default is Kelly
		//SBS8
		SBS8_GEp_param = calc_kelly_parameterization(1, SBS8_Q_sq, false);
		SBS8_GEp_param_error = calc_kelly_errors(1, SBS8_Q_sq, false);

		SBS8_GMp_param = calc_kelly_parameterization(2, SBS8_Q_sq, false);
		SBS8_GMp_param_error = calc_kelly_errors(2, SBS8_Q_sq, false);

		SBS8_GEn_param = calc_kelly_parameterization(3, SBS8_Q_sq, false);
		SBS8_GEn_param_error = calc_kelly_errors(3, SBS8_Q_sq, false);

		SBS8_GMn_param = calc_kelly_parameterization(4, SBS8_Q_sq, false);
		SBS8_GMn_param_error = calc_kelly_errors(4, SBS8_Q_sq, false);

		SBS8_sigma_R_n_param = SBS8_eps_n*(SBS8_GEn_param*SBS8_GEn_param) + SBS8_tau_n*(SBS8_GMn_param*SBS8_GMn_param);
		SBS8_sigma_R_n_param_error = calc_sigma_R_N_Err(SBS8_eps_n, SBS8_GEn_param, SBS8_GEn_param_error, SBS8_tau_n, SBS8_GMn_param, SBS8_GMn_param_error);
		

		//SBS9
		SBS9_GEp_param = calc_kelly_parameterization(1, SBS9_Q_sq, false);
		SBS9_GEp_param_error = calc_kelly_errors(1, SBS9_Q_sq, false);

		SBS9_GMp_param = calc_kelly_parameterization(2, SBS9_Q_sq, false);
		SBS9_GMp_param_error = calc_kelly_errors(2, SBS9_Q_sq, false);

		SBS9_GEn_param = calc_kelly_parameterization(3, SBS9_Q_sq, false);
		SBS9_GEn_param_error = calc_kelly_errors(3, SBS9_Q_sq, false);

		SBS9_GMn_param = calc_kelly_parameterization(4, SBS9_Q_sq, false);
		SBS9_GMn_param_error = calc_kelly_errors(4, SBS9_Q_sq, false);		

		SBS9_sigma_R_n_param = SBS9_eps_n*(SBS9_GEn_param*SBS9_GEn_param) + SBS9_tau_n*(SBS9_GMn_param*SBS9_GMn_param);
		SBS9_sigma_R_n_param_error = calc_sigma_R_N_Err(SBS9_eps_n, SBS9_GEn_param, SBS9_GEn_param_error, SBS9_tau_n, SBS9_GMn_param, SBS9_GMn_param_error);
		
	}


/////////////////////////////////////////////////////////////////////
	///  Calculation of GMn using the scale factor ratio 
// Reduce cross-section in simc is calculated from the Bosted and Riordan Paramterizations for GMn and GEn

	//SBS8
	SBS8_GMn_calc =  sqrt( (1/SBS8_tau_n)*( ( (SBS8_scale_factor_ratio)*(SBS8_sigma_R_n_param) ) - (SBS8_eps_n*SBS8_GEn_param*SBS8_GEn_param) ) );
	//Normalize:
	SBS8_GMn_calc = SBS8_GMn_calc/(mu_n*SBS8_G_D);
	SBS8_GMn_calc_stat_error = calc_GMn_err( SBS8_tau_n, SBS8_eps_n, SBS8_scale_factor_n, SBS8_scale_factor_n_error, SBS8_scale_factor_p, SBS8_scale_factor_p_error, SBS8_sigma_R_n_param, SBS8_sigma_R_n_param_error, SBS8_GEn_param, SBS8_GEn_param_error);

	//SBS9
	SBS9_GMn_calc =  sqrt( (1/SBS9_tau_n)*( ( (SBS9_scale_factor_ratio)*(SBS9_sigma_R_n_param) ) - (SBS9_eps_n*SBS9_GEn_param*SBS9_GEn_param) ) );
	//Normalize:
	SBS9_GMn_calc = SBS9_GMn_calc/(mu_n*SBS9_G_D);
	SBS9_GMn_calc_stat_error = calc_GMn_err( SBS9_tau_n, SBS9_eps_n, SBS9_scale_factor_n, SBS9_scale_factor_n_error, SBS9_scale_factor_p, SBS9_scale_factor_p_error, SBS9_sigma_R_n_param, SBS9_sigma_R_n_param_error, SBS9_GEn_param, SBS9_GEn_param_error);


/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////
	/// nTPE Calculations:

	// Longitudinal and Transverse Cross-sections:
	// For GMn we use the experimental value 
	// For GEn we use the simc method: Bosted
	bool use_GMn_exp = false;

	if( use_GMn_exp	){
		SBS8_GMn_selected = jboyd_GMn_SBS8;
		SBS8_GMn_selected_stat_error = jboyd_GMn_SBS8_stat_error;
		SBS8_GMn_selected_syst_error = jboyd_GMn_SBS8_syst_error;

		SBS9_GMn_selected = jboyd_GMn_SBS9;
		SBS9_GMn_selected_stat_error = jboyd_GMn_SBS9_stat_error;
		SBS9_GMn_selected_syst_error = jboyd_GMn_SBS9_syst_error;
	}
	else{
		SBS8_GMn_selected = SBS8_GMn_calc;
		SBS8_GMn_selected_stat_error = SBS8_GMn_calc_stat_error;
		SBS8_GMn_selected_syst_error = 0.0;

		SBS9_GMn_selected = SBS9_GMn_calc;
		SBS9_GMn_selected_stat_error = SBS9_GMn_calc_stat_error;
		SBS9_GMn_selected_syst_error = 0.0;		
	}


	//SBS8
	//Reduced C.S, neutron
	SBS8_sigma_R_n_exp = SBS8_eps_n*(SBS8_GEn_param*SBS8_GEn_param) + SBS8_tau_n*(pow(SBS8_GMn_selected*mu_n*SBS8_G_D, 2));

	SBS8_sigma_T_n = SBS8_tau_n*(pow(SBS8_GMn_selected*mu_n*SBS8_G_D, 2));
	SBS8_sigma_L_n = SBS8_GEn_param*SBS8_GEn_param;

	SBS8_sigma_T_p = SBS8_tau_p*(pow(SBS8_GMp_param*mu_p, 2));
	SBS8_sigma_L_p = SBS8_GEp_param*SBS8_GEp_param;

	//SBS9
	//Reduced C.S, neutron
	SBS9_sigma_R_n_exp = SBS9_eps_n*(SBS9_GEn_param*SBS9_GEn_param) + SBS9_tau_n*(pow(SBS9_GMn_selected*mu_n*SBS8_G_D, 2));

	SBS9_sigma_T_n = SBS9_tau_n*(pow(SBS9_GMn_selected*mu_n*SBS9_G_D, 2));
	SBS9_sigma_L_n = SBS9_GEn_param*SBS9_GEn_param;

	SBS9_sigma_T_p = SBS9_tau_p*(pow(SBS9_GMp_param*mu_p, 2));
	SBS9_sigma_L_p = SBS9_GEp_param*SBS9_GEp_param;

	// R_sigma is the Scale Factor Ratio corrected SIMC Reduced C.S.
	R_sigma = ( SBS8_scale_factor_ratio*SBS8_sigma_R_n_param )/(SBS9_scale_factor_ratio*SBS9_sigma_R_n_param);

	// Rosebluth Slopes:
	// S^{n}:
	SBS8_RS_n = (SBS8_sigma_L_n)/(SBS8_sigma_T_n);
	SBS9_RS_n = (SBS9_sigma_L_n)/(SBS9_sigma_T_n);
	// S^{p}:
	SBS8_RS_p = (SBS8_sigma_L_p)/(SBS8_sigma_T_p);
	SBS9_RS_p = (SBS9_sigma_L_p)/(SBS9_sigma_T_p);
	RS_p_average = (SBS8_RS_p + SBS9_RS_p)/2.0;

	Sn_TPE = ( R_sigma*(SBS9_sigma_T_n/SBS8_sigma_T_n) - 1 )/( SBS8_eps_n - ((R_sigma*(SBS9_sigma_T_n/SBS8_sigma_T_n))*SBS9_eps_n) );

	Q_sq_average = (SBS8_Q_sq + SBS9_Q_sq)/2.0;

	GEn_OPE = calc_kelly_parameterization( 3, Q_sq_average, true );
	GMn_OPE = calc_kelly_parameterization( 4, Q_sq_average, true );
	GEn_GMn_ratio_OPE = GEn_OPE/GMn_OPE;
	Sn_OPE = (1.0/tau_n_average)*pow(((1/mu_n)*GEn_GMn_ratio_OPE), 2);

	GEn_average = (SBS8_GEn_param + SBS9_GEn_param)/2.0;
	GMn_average = (SBS8_GMn_selected + SBS9_GMn_selected)/2.0;

	Sn_approx = ( (R_sigma*(SBS9_sigma_T_n/SBS8_sigma_T_n)) - 1 )/( (R_sigma*(SBS9_sigma_T_n/SBS8_sigma_T_n))*(SBS9_eps_n - SBS8_eps_n));

	R_corr_eps1 = (SBS8_eps_n*SBS8_sigma_L_n + SBS8_sigma_T_n)/(SBS8_eps_p*SBS8_sigma_L_p + SBS8_sigma_T_p);
	R_corr_eps2 = (SBS9_eps_n*SBS9_sigma_L_n + SBS9_sigma_T_n)/(SBS9_eps_p*SBS9_sigma_L_p + SBS9_sigma_T_p);
	A = R_corr_eps1/R_corr_eps2;
	B = (1 + SBS9_eps_p*SBS9_RS_p)/(1 + SBS8_eps_p*SBS8_RS_p);
	Sn = (B - A)/((A*SBS9_eps_n) - (B*SBS8_eps_n));

	GEn_GMn_proposal = mu_GEn_GMn_proposal/mu_n;

	TPE = Sn - (pow(GEn_GMn_proposal, 2)/tau_n_average);

}

double eps[2], eps_err[2], RS[2], RS_err[2];

double min_sigma_R, max_sigma_R;

TGraph *tg_nTPE_sigma_R;
TGraphErrors *tge_nTPE_sigma_R;

//double SBS8_eps = 0.798389481, double SBS8_sigma_R = 0.00385,  double SBS9_eps = 0.513874804,  double SBS9_sigma_R = 0.00375

void plot_RS_v_epsilson(){
	
	double SBS8_eps = SBS8_eps_n;
	double SBS8_sigma_R = SBS8_sigma_R_n_exp;

	double SBS9_eps = SBS9_eps_n;
	double SBS9_sigma_R = SBS9_sigma_R_n_exp;

	int num_points = 2;
	eps[0] = SBS8_eps; 
	eps[1] = SBS9_eps;

	eps_err[0] = SBS8_eps_err;
	eps_err[1] = SBS9_eps_err;

	RS[0] = SBS8_sigma_R;
	RS[1] = SBS9_sigma_R;

	RS_err[0] = SBS8_sigma_R_err;
	RS_err[1] = SBS9_sigma_R_err;

	double min_sigma_R, max_sigma_R;

	if( SBS8_sigma_R > SBS9_sigma_R ){
		min_sigma_R = SBS9_sigma_R;
		max_sigma_R = SBS8_sigma_R;
	}
	else{
		min_sigma_R = SBS8_sigma_R;
		max_sigma_R = SBS9_sigma_R;
	}
	double y_lower_lim = 0.8*min_sigma_R;
	double y_upper_lim = 1.03*max_sigma_R;

	tg_nTPE_sigma_R = new TGraph(num_points, eps, RS);
	tg_nTPE_sigma_R->SetTitle("nTPE - Reduced C.S., #sigma_{R} vs Virtual Photon Polarization, #epsilon");
	tg_nTPE_sigma_R->GetXaxis()->SetLimits(0.0, 1.0);
	tg_nTPE_sigma_R->GetXaxis()->SetTitle("#epsilon");

	tg_nTPE_sigma_R->GetYaxis()->SetTitle("#tauG_{M}^{2} + #epsilonG_{E}^{2}");
	tg_nTPE_sigma_R->GetYaxis()->SetRangeUser(y_lower_lim, y_upper_lim);
	tg_nTPE_sigma_R->GetYaxis()->SetMaxDigits(2);
	tg_nTPE_sigma_R->SetMarkerStyle(7);

	tge_nTPE_sigma_R = new TGraphErrors(num_points, eps, RS, eps_err, RS_err);
	tge_nTPE_sigma_R->SetTitle("nTPE - Reduced C.S., #sigma_{R} vs Virtual Photon Polarization, #epsilon");
	tge_nTPE_sigma_R->GetXaxis()->SetLimits(0.0, 1.0);
	tge_nTPE_sigma_R->GetXaxis()->SetTitle("#epsilon");

	tge_nTPE_sigma_R->GetYaxis()->SetTitle("#tauG_{M}^{2} + #epsilonG_{E}^{2}");
	tge_nTPE_sigma_R->GetYaxis()->SetRangeUser(y_lower_lim, y_upper_lim);
	tge_nTPE_sigma_R->GetYaxis()->SetMaxDigits(2);
	tge_nTPE_sigma_R->SetMarkerStyle(7);

	TCanvas *c_nTPE_sigma_R = new TCanvas("c_nTPE_sigma_R", "c_nTPE_sigma_R", 600, 500);
	c_nTPE_sigma_R->SetGrid();
	tg_nTPE_sigma_R->Draw("APN");
	tge_nTPE_sigma_R->Draw("same");

	TF1 *tf_nTPE_sigma_R = new TF1("tf_nTPE_sigma_R", "pol1", 0.0, 1.0);
	tf_nTPE_sigma_R->SetLineStyle(6);
	tf_nTPE_sigma_R->SetLineColor(1);
	tg_nTPE_sigma_R->Fit("tf_nTPE_sigma_R", "R+");

	double y_intercept = tf_nTPE_sigma_R->GetParameter(0);
	double slope = tf_nTPE_sigma_R->GetParameter(1);

	if( y_intercept > y_upper_lim ){
		tg_nTPE_sigma_R->GetYaxis()->SetRangeUser(y_lower_lim, y_intercept);
	}
	if( y_intercept < y_lower_lim ){
		tg_nTPE_sigma_R->GetYaxis()->SetRangeUser(y_intercept, y_upper_lim);
	}

	double GMn = sqrt(y_intercept/tau_n_average);
	double GEn = sqrt(slope);

	TLine *tl_eps_1 = new TLine(SBS8_eps, y_lower_lim, SBS8_eps, y_upper_lim);
	tl_eps_1->SetLineColor(17);
	tl_eps_1->SetLineStyle(3);
	tl_eps_1->Draw("same");

	TLine *tl_eps_2 = new TLine(SBS9_eps, y_lower_lim, SBS9_eps, y_upper_lim);
	tl_eps_2->SetLineColor(17);
	tl_eps_2->SetLineStyle(3);
	tl_eps_2->Draw("same");


	TText *txt_sbs8 = new TText(SBS8_eps - 0.0005, 1.01*y_lower_lim, "SBS8");
	// txt_sbs8->SetTextAngle(90);
	txt_sbs8->SetTextSize(0.035f);
	txt_sbs8->Draw("same");

	TText *txt_sbs9 = new TText(SBS9_eps - 0.0005, 1.01*y_lower_lim, "SBS9");
	// txt_sbs9->SetTextAngle(90);
	txt_sbs9->SetTextSize(0.035f);
	txt_sbs9->Draw("same");

	TArrow *ta_tauGM2 = new TArrow(0.025, y_lower_lim, 0.025, y_intercept, 0.02, "<|>");
	ta_tauGM2->SetFillColor(1);
	ta_tauGM2->Draw();

	TLatex txt_tauGMN2;
	txt_tauGMN2.SetTextSize(0.05f);
	txt_tauGMN2.SetTextAlign(13);
	txt_tauGMN2.DrawLatex(0.03, (y_intercept + y_lower_lim)/2.0, Form("#tauG_{M}^{2} #rightarrow G_{M} = %0.3f", GMn));

	TArrow *ta_GE2 = new TArrow(1.03*SBS8_eps, min_sigma_R, 1.03*SBS8_eps, max_sigma_R, 0.02, "<|>");
	ta_GE2->SetFillColor(1);
	ta_GE2->Draw();

	//check for negative slope
	double min_sigma_R_line_y_val;
	if( slope < 0 ){
		min_sigma_R_line_y_val = max_sigma_R;
	}
	else{
		min_sigma_R_line_y_val = min_sigma_R;
	}

	TLine *tl_GEn_horiz = new TLine(SBS9_eps, min_sigma_R_line_y_val, SBS8_eps, min_sigma_R_line_y_val);
	tl_GEn_horiz->SetLineStyle(3);
	tl_GEn_horiz->SetLineColor(17);
	tl_GEn_horiz->Draw("same");

	TLatex txt_GE2;
	txt_GE2.SetTextSize(0.05f);
	txt_GE2.SetTextAlign(13);
	txt_GE2.DrawLatex(1.035*SBS8_eps, 1.02*min_sigma_R, "G_{E}^{2}" );

	cout << "----------------------------------------" << endl;
	cout << Form("Intercept ( #tau (G_{M})^{2} ) = %E", tf_nTPE_sigma_R->GetParameter(0)) << endl;
	cout << Form("Slope ( (G_{E})^{2} ) = %E", tf_nTPE_sigma_R->GetParameter(1) ) << endl;

	double OPE_at_x1 = y_intercept + pow( (calc_kelly_parameterization(3, Q_sq_average, false)), 2);

	if( slope > 0 ){
		TLine *tl_OPE = new TLine(0.0, y_intercept, 1.0, OPE_at_x1);
		tl_OPE->SetLineColor(17);
		tl_OPE->SetLineStyle(2);
		tl_OPE->Draw("same");

		TText *txt_pol_trans = new TText(0.1, 1.004*y_intercept, "Polarization Transfer Slope");
		double OPE_rise = OPE_at_x1 - y_intercept;
		double OPE_line_angle = atan(1100*OPE_rise)*TMath::RadToDeg();
		txt_pol_trans->SetTextAngle(OPE_line_angle);
		txt_pol_trans->SetTextSize(0.025f);
		txt_pol_trans->SetTextColor(16);
		txt_pol_trans->Draw("same");		

		TLine *tl_flat = new TLine(0.0, y_intercept, 1.0, y_intercept);
		tl_flat->SetLineColor(1);
		tl_flat->SetLineStyle(7);
		tl_flat->Draw("same");
	}

	double GEn2_from_diff = fabs(SBS8_sigma_R_n_exp - SBS9_sigma_R_n_exp);
	double GEn_from_diff = sqrt(GEn2_from_diff);

	double GMn2_from_int = tf_nTPE_sigma_R->GetParameter(0)/tau_n_average;
	double GMn_from_int = sqrt(GMn2_from_int);

	cout << "---------------------" << endl;
	cout << "SBS8 tau_n = " << SBS8_tau_n << endl;
	cout << "SBS8 eps_n = " << SBS8_eps_n << endl;
	cout << "SBS8 G_D = " << SBS8_G_D << endl;
	cout << "---" << endl;
	cout << "SBS9 tau_n = " << SBS9_tau_n << endl;
	cout << "SBS9 eps_n = " << SBS9_eps_n << endl;
	cout << "SBS9 G_D = " << SBS9_G_D << endl;
	cout << "------------" << endl;
	cout << "GEn = " << GEn_from_diff << endl;
	cout << "GMn = " << GMn_from_int << endl;
	cout << "GMn norm. = " << GMn_from_int/(mu_n*( (SBS8_G_D + SBS9_G_D)/2.0)) << endl;



}


#endif
