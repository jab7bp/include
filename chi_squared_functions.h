#ifndef CHI_SQUARED_FUNCS_H
#define CHI_SQUARED_FUNCS_H

class cl_chi2_histo {
	
public:

	TH1D *data_dx_histo, *sim_dx_histo, *sim_dx_p_histo, *sim_dx_n_histo;
	double dx_global_scale_p, dx_global_scale_n, ChiSquared_low_x, ChiSquared_high_x;
	vector<double> *chiSquaredValues;

};

cl_chi2_histo cl_dx_chi2_data_simc;

double CalculateChiSquaredInRange(const cl_chi2_histo& cl_dx_chi2_data_simc) {
    
	double scale_p = cl_dx_chi2_data_simc.dx_global_scale_p;
	double scale_n = cl_dx_chi2_data_simc.dx_global_scale_n;
	TH1D *hin_dx_cutdy = cl_dx_chi2_data_simc.data_dx_histo;
	TH1D *h_simc_dx_p = cl_dx_chi2_data_simc.sim_dx_p_histo;
	TH1D *h_simc_dx_n = cl_dx_chi2_data_simc.sim_dx_n_histo;
	double ChiSquared_low_x = cl_dx_chi2_data_simc.ChiSquared_low_x;
	double ChiSquared_high_x = cl_dx_chi2_data_simc.ChiSquared_high_x; 
	vector<double> *chiSquaredValues = cl_dx_chi2_data_simc.chiSquaredValues;
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
    chiSquaredValues->push_back(chiSquared);

    return chiSquared;
}

void MinimizationFunctionInRange(int &npar, double *gin, double &result, double *param, int flag){
	// global_scale_p = param[0];
	// global_scale_n = param[1];
	cl_dx_chi2_data_simc.dx_global_scale_p = param[0];
	cl_dx_chi2_data_simc.dx_global_scale_n = param[1];

	//BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
	result = CalculateChiSquaredInRange(cl_dx_chi2_data_simc);
}

double CalculateChiSquared(const cl_chi2_histo& cl_dx_chi2_data_simc) {
    
	double scale_p = cl_dx_chi2_data_simc.dx_global_scale_p;
	double scale_n = cl_dx_chi2_data_simc.dx_global_scale_n;
	TH1D *hin_dx_cutdy = cl_dx_chi2_data_simc.data_dx_histo;
	TH1D *h_simc_dx_p = cl_dx_chi2_data_simc.sim_dx_p_histo;
	TH1D *h_simc_dx_n = cl_dx_chi2_data_simc.sim_dx_n_histo;
	double ChiSquared_low_x = cl_dx_chi2_data_simc.ChiSquared_low_x;
	double ChiSquared_high_x = cl_dx_chi2_data_simc.ChiSquared_high_x; 
	vector<double> *chiSquaredValues = cl_dx_chi2_data_simc.chiSquaredValues;
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
    chiSquaredValues->push_back(chiSquared);

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
#endif