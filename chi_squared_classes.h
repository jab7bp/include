#ifndef CHI_SQUARED_CLASSES_H
#define CHI_SQUARED_CLASSES_H


// class cl_chi2_histo {
	
// public:

// 	TH1D *data_dx_histo, *sim_dx_histo, *sim_dx_p_histo, *sim_dx_n_histo, *sim_chi2_dx_cutdy_histo;
//     TH1D *data_dx_BGsub_histo;

//     TH2D *h2_chi2_dxdy_cutdy,*h2_chi2_dxdy_cutdy_SigmaMult;
// 	double dx_global_scale_p, dx_global_scale_n, ChiSquared_low_x, ChiSquared_high_x;
//     double dx_global_scale_p_error, dx_global_scale_n_error;
//     double dx_global_dy_cut_min, dx_global_dy_cut_max, ChiSquared_dy_cut_low_x, ChiSquared_dy_cut_high_x;
//     double dx_trough_min, dx_trough_max, sim_dx_trough_integral;
//     double dy_global_sigma_mult;
//     double dy_cut_pn, dy_cut_pn_sigma;
// 	vector<double> *chiSquaredByScaleValues, *chiSquaredDyCutValues, *chiSquaredDyCutValuesSigmaMult;

//     double dx_global_BGsub_scale_p, dx_global_BGsub_scale_n;
//     double global_BGsub_scale_p, global_BGsub_scale_n;
// };


class cl_chi2_minimization{
    
	public:
		//BGsub stuff
	    TH1D *h_data_dx = NULL, *h_data_dx_with_BGsub = NULL;
	    TH1D *h_simc_dx = NULL, *h_simc_dx_p = NULL, *h_simc_dx_n = NULL;
	    TH1D *h_simc_dx_p_BeforeScaling = NULL, *h_simc_dx_n_BeforeScaling = NULL;

	    double chiSquared_low_x = -99999.9, chiSquared_high_x = 99999.0;
	    double scale_p = 0.0, scale_n = 0.0;
	    double scale_p_error = 0.0, scale_n_error = 0.0;
	    double global_scale_p, global_scale_n;

	    //NON BGsub minimization stuff:
	    TH1D *h_simc_chi2_cutdy = NULL;
	    TH2D *h2_chi2_dxdy_cutdy = NULL;
	    TH2D *h2_chi2_dxdy_cutdy_SigmaMult = NULL;

	    double dy_cut_pn = 0.0;
	    double dy_cut_pn_sigma = 0.0;
	    double dy_cut_sigmaMult = 1.0;

	    //axis limits:
	    double dy_cut_min = -0.5;
	    double dy_cut_max = 0.5;

	    double dx_trough_min = -0.3;
	    double dx_trough_max = 0.3;
	    double simc_dx_trough_integral = 0.0;

	    double dx_Xaxis_min;
	    double dx_Xaxis_max;

	    vector<double> *chiSquared_dyCut_values = {};
	    vector<double> *chiSquared_dyCut_values_sigmaMult = {};

	    vector<double> chiSquared_values = {};
	    vector<double> scale_n_values = {};
	    vector<double> scale_p_values = {};
	    // vector<vector<double>> *scale_n_p_values{};

	void Reset(){

		cout << "*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" << endl;
		cout << endl;
		cout << "       Resetting chi-squared minimization class" << endl;
		cout << endl;
		cout << "*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" << endl << endl;

		//BGsub stuff
	    TH1D *h_data_dx = NULL, *h_data_dx_with_BGsub = NULL;
	    TH1D *h_simc_dx = NULL, *h_simc_dx_p = NULL, *h_simc_dx_n = NULL;

	    double chiSquared_low_x = -99999.9, chiSquared_high_x = 99999.0;
	    double scale_p = 0.0, scale_n = 0.0;
	    double scale_p_error = 0.0, scale_n_error = 0.0;
	    double scale_p_temp = 0.0, scale_n_temp = 0.0;


	    //NON BGsub minimization stuff:
	    TH1D *h_simc_chi2_cutdy = NULL;
	    TH2D *h2_chi2_dxdy_cutdy = NULL;
	    TH2D *h2_chi2_dxdy_cutdy_SigmaMult = NULL;

	    double dy_cut_pn = 0.0;
	    double dy_cut_pn_sigma = 0.0;
	    double dy_cut_sigmaMult = 1.0;
	    double dy_cut_min = -0.5;
	    double dy_cut_max = 0.5;

	    double dx_trough_min = -0.3;
	    double dx_trough_max = 0.3;
	    double simc_dx_trough_integral = 0.0;

	    vector<double> *chiSquared_dyCut_values = {};
	    vector<double> *chiSquared_dyCut_values_sigmaMult = {};

	    vector<double> chiSquared_values = {};
	    vector<double> scale_n_values = {};
	    vector<double> scale_p_values = {};
	    // vector<vector<double>> *scale_n_p_values{};
	}

};

#endif