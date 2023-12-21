#ifndef CHI_SQUARED_FUNCS_H
#define CHI_SQUARED_FUNCS_H

class cl_chi2_histo {
	
public:

	TH1D *data_dx_histo, *sim_dx_histo, *sim_dx_p_histo, *sim_dx_n_histo, *sim_chi2_dx_cutdy_histo;
    TH1D *data_dx_BGsub_histo;

    TH2D *h2_chi2_dxdy_cutdy,*h2_chi2_dxdy_cutdy_SigmaMult;
	double dx_global_scale_p, dx_global_scale_n, ChiSquared_low_x, ChiSquared_high_x;
    double dx_global_scale_p_error, dx_global_scale_n_error;
    double dx_global_dy_cut_min, dx_global_dy_cut_max, ChiSquared_dy_cut_low_x, ChiSquared_dy_cut_high_x;
    double dx_trough_min, dx_trough_max, sim_dx_trough_integral;
    double dy_global_sigma_mult;
    double dy_cut_pn, dy_cut_pn_sigma;
	vector<double> *chiSquaredByScaleValues, *chiSquaredDyCutValues, *chiSquaredDyCutValuesSigmaMult;

    double dx_global_BGsub_scale_p, dx_global_BGsub_scale_n;
    double global_BGsub_scale_p, global_BGsub_scale_n;
};

cl_chi2_histo cl_dx_chi2_data_simc, cl_dx_chi2_data_simc_BGsub;

double CalculateChiSquaredInRangeByScale(const cl_chi2_histo& cl_dx_chi2_data_simc) {
    
	double scale_p = cl_dx_chi2_data_simc.dx_global_scale_p;
	double scale_n = cl_dx_chi2_data_simc.dx_global_scale_n;
	TH1D *hin_dx_cutdy = cl_dx_chi2_data_simc.data_dx_histo;
	TH1D *h_simc_dx_p = cl_dx_chi2_data_simc.sim_dx_p_histo;
	TH1D *h_simc_dx_n = cl_dx_chi2_data_simc.sim_dx_n_histo;
	double ChiSquared_low_x = cl_dx_chi2_data_simc.ChiSquared_low_x;
	double ChiSquared_high_x = cl_dx_chi2_data_simc.ChiSquared_high_x; 
	vector<double> *chiSquaredByScaleValues = cl_dx_chi2_data_simc.chiSquaredByScaleValues;
    double chiSquared = 0.0;

    int low_bin = hin_dx_cutdy->GetXaxis()->FindBin(ChiSquared_low_x);
    int high_bin = hin_dx_cutdy->GetXaxis()->FindBin(ChiSquared_high_x);

    // for (int bin = 1; bin <= hin_dx_cutdy->GetNbinsX(); ++bin) {
    for (int bin = low_bin; bin <= high_bin; ++bin) {
        double content_data = hin_dx_cutdy->GetBinContent(bin);
        double content_expected = scale_p * h_simc_dx_p->GetBinContent(bin) + scale_n * h_simc_dx_n->GetBinContent(bin);

        double error_data = hin_dx_cutdy->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    // Store the chi-squared value
    chiSquaredByScaleValues->push_back(chiSquared);

    return chiSquared;
}

void MinimizationFunctionInRangeByScale(int &npar, double *gin, double &result, double *param, int flag){
	// global_scale_p = param[0];
	// global_scale_n = param[1];
	cl_dx_chi2_data_simc.dx_global_scale_p = param[0];
	cl_dx_chi2_data_simc.dx_global_scale_n = param[1];

	//BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
	result = CalculateChiSquaredInRangeByScale(cl_dx_chi2_data_simc);
}

//-------------------------------------------
//-------------------------------------------

double CalculateChiSquaredInRangeByScale_BGsub(const cl_chi2_histo& cl_dx_chi2_data_simc_BGsub) {
    
    double scale_p = cl_dx_chi2_data_simc_BGsub.dx_global_BGsub_scale_p;
    double scale_n = cl_dx_chi2_data_simc_BGsub.dx_global_BGsub_scale_n;
    TH1D *hin_dx_BGsub_cutdy = cl_dx_chi2_data_simc_BGsub.data_dx_BGsub_histo;
    TH1D *h_simc_dx_p = cl_dx_chi2_data_simc_BGsub.sim_dx_p_histo;
    TH1D *h_simc_dx_n = cl_dx_chi2_data_simc_BGsub.sim_dx_n_histo;

    double ChiSquared_low_x = cl_dx_chi2_data_simc_BGsub.ChiSquared_low_x;
    double ChiSquared_high_x = cl_dx_chi2_data_simc_BGsub.ChiSquared_high_x; 

    vector<double> *chiSquaredByScaleValues = cl_dx_chi2_data_simc_BGsub.chiSquaredByScaleValues;
    double chiSquared = 0.0;

    int low_bin = hin_dx_BGsub_cutdy->GetXaxis()->FindBin(ChiSquared_low_x);
    int high_bin = hin_dx_BGsub_cutdy->GetXaxis()->FindBin(ChiSquared_high_x);

    // for (int bin = 1; bin <= hin_dx_BGsub_cutdy->GetNbinsX(); ++bin) {
    for (int bin = low_bin; bin <= high_bin; ++bin) {
        double content_data = hin_dx_BGsub_cutdy->GetBinContent(bin);
        double content_expected = scale_p * h_simc_dx_p->GetBinContent(bin) + scale_n * h_simc_dx_n->GetBinContent(bin);

        double error_data = hin_dx_BGsub_cutdy->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    // Store the chi-squared value
    chiSquaredByScaleValues->push_back(chiSquared);

    return chiSquared;
}

void MinimizationFunctionInRangeByScale_BGsub(int &npar, double *gin, double &result, double *param, int flag){
    // global_scale_p = param[0];
    // global_scale_n = param[1];
    cl_dx_chi2_data_simc_BGsub.dx_global_BGsub_scale_p = param[0];
    cl_dx_chi2_data_simc_BGsub.dx_global_BGsub_scale_n = param[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByScale_BGsub(cl_dx_chi2_data_simc_BGsub);
}

double CalculateChiSquaredInRangeByScaleWithNormalization(const cl_chi2_histo& cl_dx_chi2_data_simc) {
    
    double scale_p = cl_dx_chi2_data_simc.dx_global_scale_p;
    double scale_n = cl_dx_chi2_data_simc.dx_global_scale_n;
    TH1D *hin_dx_cutdy = cl_dx_chi2_data_simc.data_dx_histo;
    TH1D *h_simc_dx_p = cl_dx_chi2_data_simc.sim_dx_p_histo;
    TH1D *h_simc_dx_n = cl_dx_chi2_data_simc.sim_dx_n_histo;
    double ChiSquared_low_x = cl_dx_chi2_data_simc.ChiSquared_low_x;
    double ChiSquared_high_x = cl_dx_chi2_data_simc.ChiSquared_high_x; 
    vector<double> *chiSquaredByScaleValues = cl_dx_chi2_data_simc.chiSquaredByScaleValues;
    double chiSquared = 0.0;

    int low_bin = hin_dx_cutdy->GetXaxis()->FindBin(ChiSquared_low_x);
    int high_bin = hin_dx_cutdy->GetXaxis()->FindBin(ChiSquared_high_x);

    // for (int bin = 1; bin <= hin_dx_cutdy->GetNbinsX(); ++bin) {

    TH1D *h_simc_dx_added = new TH1D("h_simc_dx_added", "h_simc_dx_added", h_simc_dx_p->GetNbinsX(), h_simc_dx_p->GetXaxis()->GetXmin(), h_simc_dx_p->GetXaxis()->GetXmax());
    
    for( int i = 1; i <= h_simc_dx_added->GetNbinsX(); i++){
        h_simc_dx_added->SetBinContent(i, scale_p*h_simc_dx_p->GetBinContent(i) + scale_n*h_simc_dx_n->GetBinContent(i) );
    }

    //Normalize histograms
    hin_dx_cutdy->Scale(1.0/hin_dx_cutdy->Integral());
    h_simc_dx_added->Scale(1.0/h_simc_dx_added->Integral());

    for (int bin = low_bin; bin <= high_bin; ++bin) {
        double content_data = hin_dx_cutdy->GetBinContent(bin);
        // double content_expected = scale_p * h_simc_dx_p->GetBinContent(bin) + scale_n * h_simc_dx_n->GetBinContent(bin);
        // double content_expected = h_simc_dx_added->GetBinContent(bin);
        double content_expected = h_simc_dx_added->GetBinContent(bin);

        double error_data = hin_dx_cutdy->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    // Store the chi-squared value
    chiSquaredByScaleValues->push_back(chiSquared);

    return chiSquared;
}

void MinimizationFunctionInRangeByScaleWithNormalization(int &npar, double *gin, double &result, double *param, int flag){
    // global_scale_p = param[0];
    // global_scale_n = param[1];
    cl_dx_chi2_data_simc.dx_global_scale_p = param[0];
    cl_dx_chi2_data_simc.dx_global_scale_n = param[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByScaleWithNormalization(cl_dx_chi2_data_simc);
}

double CalculateChiSquaredInRangeByDyCut(const cl_chi2_histo& cl_dx_chi2_data_simc) {
    TH2D *h2_dxdy_trough = cl_dx_chi2_data_simc.h2_chi2_dxdy_cutdy;
    int h2_dxdy_xbins = h2_dxdy_trough->GetNbinsX();
    int h2_dxdy_ybins = h2_dxdy_trough->GetNbinsY();

    double dy_cut_min = cl_dx_chi2_data_simc.dx_global_dy_cut_min;
    int dy_cut_min_bin = h2_dxdy_trough->GetXaxis()->FindBin(dy_cut_min);

    double dy_cut_max = cl_dx_chi2_data_simc.dx_global_dy_cut_max;
    int dy_cut_max_bin = h2_dxdy_trough->GetXaxis()->FindBin(dy_cut_max); 

    for( int xbin = 0; xbin < h2_dxdy_xbins; xbin++ ){
        for( int ybin = 0; ybin < h2_dxdy_ybins; ybin++ ){
            if( (xbin < dy_cut_min_bin) || (xbin > dy_cut_max_bin) ){
                h2_dxdy_trough->SetBinContent(xbin, ybin, 0.0);
            }
        }
    }

    double trough_min = cl_dx_chi2_data_simc.dx_trough_min;
    double trough_max = cl_dx_chi2_data_simc.dx_trough_max;
    int data_trough_min_bin = h2_dxdy_trough->GetYaxis()->FindBin(trough_min);
    int data_trough_max_bin = h2_dxdy_trough->GetYaxis()->FindBin(trough_max);
    TH1D *h_dxdy_ProjY = (TH1D*)h2_dxdy_trough->ProjectionY();
    double h2_dxdy_ProjY_trough_integral = h_dxdy_ProjY->Integral(data_trough_min_bin, data_trough_max_bin);
    h_dxdy_ProjY->Delete();

    double h_sim_trough_integral = cl_dx_chi2_data_simc.sim_dx_trough_integral;

    double integral_difference = abs(h2_dxdy_ProjY_trough_integral - h_sim_trough_integral);
    // cout << "----------------------------------------------" << endl;
    // cout << "dy_min: " << dy_cut_min << ", dy_max: " << dy_cut_max << endl;
    // cout << "Data trough integral: " << h2_dxdy_ProjY_trough_integral << endl;
    // cout << "sim trough integral: " << h2_sim_trough_integral << endl;
    // cout << "Integral difference: " << integral_difference << endl;
    // cout << "----------------------------------------------" << endl;
    return integral_difference;

}

void MinimizationFunctionInRangeByDyCut(int &npar, double *gin, double &result, double *parDyCut, int flag){
    // global_scale_p = param[0];
    // global_scale_n = param[1];
    cl_dx_chi2_data_simc.dx_global_dy_cut_min = parDyCut[0];
    cout << "par[0]: " << parDyCut[0] << ", par[1]: " << parDyCut[1] << endl;
    cl_dx_chi2_data_simc.dx_global_dy_cut_max = parDyCut[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByDyCut(cl_dx_chi2_data_simc);
}

double CalculateChiSquaredInRangeDyCutSigmaMult(const cl_chi2_histo& cl_dx_chi2_data_simc) {
    TH2D *h2_dxdy_SigmaMult_trough = (TH2D*)cl_dx_chi2_data_simc.h2_chi2_dxdy_cutdy_SigmaMult->Clone("h2_dxdy_SigmaMult_trough");
    int h2_dxdy_SigmaMult_xbins = h2_dxdy_SigmaMult_trough->GetNbinsX();
    int h2_dxdy_SigmaMult_ybins = h2_dxdy_SigmaMult_trough->GetNbinsY();

    double dy_cut_minimizer_sigma = cl_dx_chi2_data_simc.dy_cut_pn_sigma;
    double dy_cut_minimizer_pn = cl_dx_chi2_data_simc.dy_cut_pn;
    double dy_cut_minimizer_mult = cl_dx_chi2_data_simc.dy_global_sigma_mult;

    double dy_cut_SigmaMult_min = dy_cut_minimizer_pn - dy_cut_minimizer_mult*dy_cut_minimizer_sigma;
    int dy_cut_SigmaMult_min_bin = h2_dxdy_SigmaMult_trough->GetXaxis()->FindBin(dy_cut_SigmaMult_min);

    double dy_cut_SigmaMult_max = dy_cut_minimizer_pn + dy_cut_minimizer_mult*dy_cut_minimizer_sigma;
    int dy_cut_SigmaMult_max_bin = h2_dxdy_SigmaMult_trough->GetXaxis()->FindBin(dy_cut_SigmaMult_max); 

    for( int xbin = 0; xbin < h2_dxdy_SigmaMult_xbins; xbin++ ){
        for( int ybin = 0; ybin < h2_dxdy_SigmaMult_ybins; ybin++ ){
            if( (xbin <= dy_cut_SigmaMult_min_bin) || (xbin >= dy_cut_SigmaMult_max_bin) ){
                h2_dxdy_SigmaMult_trough->SetBinContent(xbin, ybin, 0.0);
            }
        }
    }

    double trough_min = cl_dx_chi2_data_simc.dx_trough_min;
    double trough_max = cl_dx_chi2_data_simc.dx_trough_max;
    int data_trough_min_bin = h2_dxdy_SigmaMult_trough->GetYaxis()->FindBin(trough_min);
    int data_trough_max_bin = h2_dxdy_SigmaMult_trough->GetYaxis()->FindBin(trough_max);
    TH1D *h_dxdy_ProjY = h2_dxdy_SigmaMult_trough->ProjectionY("h_dxdy_ProjY");
    double h2_dxdy_SigmaMult_ProjY_trough_integral = h_dxdy_ProjY->Integral(201, 231);

    double h_sim_trough_integral = cl_dx_chi2_data_simc.sim_dx_trough_integral;

    double SigmaMult_integral_difference = abs(h2_dxdy_SigmaMult_ProjY_trough_integral - h_sim_trough_integral);
    // cout << "----------------------------------------------" << endl;
    // cout << "dy_min: " << dy_cut_min << ", dy_max: " << dy_cut_max << endl;
    // cout << "Data trough integral: " << h2_dxdy_SigmaMult_ProjY_trough_integral << endl;
    // cout << "sim trough integral: " << h_sim_trough_integral << endl;
    // cout << "SigmaMult_integral_difference: " << SigmaMult_integral_difference << endl;
    // cout << "----------------------------------------------" << endl;

    h2_dxdy_SigmaMult_trough->Delete();
    h_dxdy_ProjY->Delete();
    return SigmaMult_integral_difference;

}

void MinimizationFunctionInRangeDyCutSigmaMult(int &npar, double *gin, double &result, double *parSigmaMult, int flag){
    // global_scale_p = param[0];
    // global_scale_n = param[1];
    cl_dx_chi2_data_simc.dy_global_sigma_mult = parSigmaMult[0];
    // cout << "Sigma Mult: " << parSigmaMult[0] << endl;

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeDyCutSigmaMult(cl_dx_chi2_data_simc);
}


double CalculateChiSquared(const cl_chi2_histo& cl_dx_chi2_data_simc) {
    
	double scale_p = cl_dx_chi2_data_simc.dx_global_scale_p;
	double scale_n = cl_dx_chi2_data_simc.dx_global_scale_n;
	TH1D *hin_dx_cutdy = cl_dx_chi2_data_simc.data_dx_histo;
	TH1D *h_simc_dx_p = cl_dx_chi2_data_simc.sim_dx_p_histo;
	TH1D *h_simc_dx_n = cl_dx_chi2_data_simc.sim_dx_n_histo;
	double ChiSquared_low_x = cl_dx_chi2_data_simc.ChiSquared_low_x;
	double ChiSquared_high_x = cl_dx_chi2_data_simc.ChiSquared_high_x; 
	vector<double> *chiSquaredByScaleValues = cl_dx_chi2_data_simc.chiSquaredByScaleValues;
    double chiSquared = 0.0;

    for (int bin = 1; bin <= hin_dx_cutdy->GetNbinsX(); ++bin) {

        double content_data = hin_dx_cutdy->GetBinContent(bin);
        double content_expected = scale_p * h_simc_dx_p->GetBinContent(bin) + scale_n * h_simc_dx_n->GetBinContent(bin);

        double error_data = hin_dx_cutdy->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    // Store the chi-squared value
    chiSquaredByScaleValues->push_back(chiSquared);

    return chiSquared;
}

void MinimizationFunction(int &npar, double *gin, double &result, double *param, int flag){
	// global_scale_p = param[0];
	// global_scale_n = param[1];

	cl_dx_chi2_data_simc.dx_global_scale_p = param[0];
	cl_dx_chi2_data_simc.dx_global_scale_n = param[1];

	// result = CalculateChiSquared(global_scale_p, global_scale_n);
	//BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
	result = CalculateChiSquared(cl_dx_chi2_data_simc);
}

double SimpleCalculateChiSquared(TH1D* observed, TH1D* expected) {
    if (!observed || !expected) {
        std::cerr << "One or both histograms are nullptr." << std::endl;
        return -1.0; // Return a negative value to indicate an error
    }

    if (observed->GetNbinsX() != expected->GetNbinsX()) {
        std::cerr << "Histogram bin counts do not match!" << std::endl;
        return -1.0;
    }

    double chi2 = 0.0;
    for (int i = 1; i <= observed->GetNbinsX(); ++i) {
        double obsContent = observed->GetBinContent(i);
        double expContent = expected->GetBinContent(i);

        if (expContent > 0.0) { // To avoid division by zero
            chi2 += std::pow(obsContent - expContent, 2) / expContent;
        }
    }

    return chi2;
}

double SimpleCalculateChiSquared_InRange(TH1D* observed, TH1D* expected, double min_x, double max_x) {
    
    double chiSquared = 0.0;

    int low_bin = observed->GetXaxis()->FindBin(min_x);
    int high_bin = expected->GetXaxis()->FindBin(max_x);

    // for (int bin = 1; bin <= hin_dx_BGsub_cutdy->GetNbinsX(); ++bin) {
    for (int bin = low_bin; bin <= high_bin; ++bin) {
        double content_data = observed->GetBinContent(bin);
        double content_expected = expected->GetBinContent(bin);

        double error_data = observed->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }


    return chiSquared;
}
#endif