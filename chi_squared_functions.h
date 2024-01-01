#ifndef CHI_SQUARED_FUNCS_H
#define CHI_SQUARED_FUNCS_H

cl_chi2_minimization cl_dx_chi2_data_simc;
cl_chi2_minimization cl_dx_chi2_inelastics, cl_dx_chi2_dyAnticut, cl_dx_chi2_BGsub_final;
cl_chi2_minimization cl_dx_chi2_data_simc_BGsub_global;


double get_chiSquared_low_high_x( int kine, TString low_high, double mean, double sigma ){
    //using sigma multiplier of 2:
    double sigma_mult = 1.0;
    double chiSquared_low_high_x;

    if( low_high == "low" ){
        chiSquared_low_high_x = mean - ( sigma_mult * sigma );
    }
    if( low_high == "high" ){
        chiSquared_low_high_x = mean + ( sigma_mult * sigma );
    }

    return chiSquared_low_high_x;
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

double CalculateChiSquaredInRangeByScale(const cl_chi2_minimization& cl_dx_chi2_data_simc) {
    
	double minimizer_scale_p = cl_dx_chi2_data_simc.scale_p;
	double minimizer_scale_n = cl_dx_chi2_data_simc.scale_n;
	TH1D *hin_dx_cutdy = cl_dx_chi2_data_simc.h_data_dx;
	TH1D *h_simc_dx_p = cl_dx_chi2_data_simc.h_simc_dx_p;
	TH1D *h_simc_dx_n = cl_dx_chi2_data_simc.h_simc_dx_n;
	double ChiSquared_low_x = cl_dx_chi2_data_simc.chiSquared_low_x;
	double ChiSquared_high_x = cl_dx_chi2_data_simc.chiSquared_high_x; 
	vector<double> chiSquared_values = cl_dx_chi2_data_simc.chiSquared_values;
    double chiSquared = 0.0;

    int low_bin = hin_dx_cutdy->GetXaxis()->FindBin(ChiSquared_low_x);
    int high_bin = hin_dx_cutdy->GetXaxis()->FindBin(ChiSquared_high_x);

    // for (int bin = 1; bin <= hin_dx_cutdy->GetNbinsX(); ++bin) {
    for (int bin = low_bin; bin <= high_bin; ++bin) {
        double content_data = hin_dx_cutdy->GetBinContent(bin);
        double content_expected = minimizer_scale_p*h_simc_dx_p->GetBinContent(bin) + minimizer_scale_n*h_simc_dx_n->GetBinContent(bin);

        double error_data = hin_dx_cutdy->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    // Store the chi-squared value
    chiSquared_values.push_back(chiSquared);

    cout << "data_simc - scale_n: " << minimizer_scale_n << ", scale_p: " << minimizer_scale_p << endl;

    return chiSquared;
}

void MinimizationFunctionInRangeByScale(int &npar, double *gin, double &result, double *param, int flag){
	// global_scale_p = param[0];
	// global_scale_n = param[1];
	cl_dx_chi2_data_simc.scale_p = param[0];
	cl_dx_chi2_data_simc.scale_n = param[1];

	//BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
	result = CalculateChiSquaredInRangeByScale(cl_dx_chi2_data_simc);
}

//-------------------------------------------
//-------------------------------------------

double CalculateChiSquaredInRangeByScale_BGsub(const cl_chi2_minimization& cl_dx_chi2_data_simc_BGsub_global) {
    
    double minimizer_scale_p = cl_dx_chi2_data_simc_BGsub_global.global_scale_p;
    double minimizer_scale_n = cl_dx_chi2_data_simc_BGsub_global.global_scale_n;
    TH1D *minimizer_h_data_dx_with_BGsub = cl_dx_chi2_data_simc_BGsub_global.h_data_dx_with_BGsub;
    TH1D *minimizer_h_simc_dx_p = cl_dx_chi2_data_simc_BGsub_global.h_simc_dx_p;
    TH1D *minimizer_h_simc_dx_n = cl_dx_chi2_data_simc_BGsub_global.h_simc_dx_n;

    double ChiSquared_low_x = cl_dx_chi2_data_simc_BGsub_global.chiSquared_low_x;
    double ChiSquared_high_x = cl_dx_chi2_data_simc_BGsub_global.chiSquared_high_x; 

    double chiSquared = 0.0;

    int low_bin = minimizer_h_data_dx_with_BGsub->GetXaxis()->FindBin(ChiSquared_low_x);
    int high_bin = minimizer_h_data_dx_with_BGsub->GetXaxis()->FindBin(ChiSquared_high_x);

    // for (int bin = 1; bin <= hin_dx_BGsub_cutdy->GetNbinsX(); ++bin) {
    for (int bin = low_bin; bin <= high_bin; ++bin) {
        double content_data = minimizer_h_data_dx_with_BGsub->GetBinContent(bin);
        double content_expected = ( minimizer_scale_p * minimizer_h_simc_dx_p->GetBinContent(bin) ) + ( minimizer_scale_n * minimizer_h_simc_dx_n->GetBinContent(bin) );

        double error_data = minimizer_h_data_dx_with_BGsub->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    // Store the chi-squared value
    // (cl_dx_chi2_data_simc_BGsub_global.chiSquared_values).push_back(chiSquared);

    return chiSquared;
}

void MinimizationFunctionInRangeByScale_BGsub(int &npar, double *gin, double &result, double *param, int flag){
    // global_dyAnticut_scale_p = param[0];
    // global_dyAnticut_scale_n = param[1];

    cl_dx_chi2_data_simc_BGsub_global.global_scale_p = param[0];
    cl_dx_chi2_data_simc_BGsub_global.global_scale_n = param[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByScale_BGsub(cl_dx_chi2_data_simc_BGsub_global);;
}

double CalculateChiSquaredIn_CustomRange_ByScale_BGsub(const cl_chi2_minimization& cl_dx_chi2_data_simc_BGsub_global, double low_x, double high_x) {
    
    double minimizer_scale_p = cl_dx_chi2_data_simc_BGsub_global.global_scale_p;
    double minimizer_scale_n = cl_dx_chi2_data_simc_BGsub_global.global_scale_n;
    TH1D *minimizer_h_data_dx_with_BGsub = (TH1D*)(cl_dx_chi2_data_simc_BGsub_global.h_data_dx_with_BGsub)->Clone("minimizer_h_data_dx_with_BGsub");;
    TH1D *minimizer_h_simc_dx_p = (TH1D*)(cl_dx_chi2_data_simc_BGsub_global.h_simc_dx_p)->Clone("minimizer_h_simc_dx_p");
    TH1D *minimizer_h_simc_dx_n = (TH1D*)(cl_dx_chi2_data_simc_BGsub_global.h_simc_dx_n)->Clone("minimizer_h_simc_dx_n");

    double nBins = minimizer_h_simc_dx_p->GetNbinsX();
    double h_xmin = minimizer_h_simc_dx_p->GetXaxis()->GetXmin();
    double h_xmax = minimizer_h_simc_dx_p->GetXaxis()->GetXmax();

    TH1D *minimizer_h_simc_dx = new TH1D("minimizer_h_simc_dx", "minimizer_h_simc_dx", nBins, h_xmin, h_xmax);
    minimizer_h_simc_dx->Add( minimizer_h_simc_dx_p, minimizer_h_simc_dx_n, minimizer_scale_p, minimizer_scale_n);

    // //Normalize both integrals:
    // minimizer_h_data_dx_with_BGsub->Scale(1.0/(minimizer_h_data_dx_with_BGsub->Integral()));
    // minimizer_h_simc_dx->Scale(1.0/(minimizer_h_simc_dx->Integral()));

    double chiSquared = 0.0;

    int low_bin = minimizer_h_data_dx_with_BGsub->GetXaxis()->FindBin(low_x);
    int high_bin = minimizer_h_data_dx_with_BGsub->GetXaxis()->FindBin(high_x);

    // for (int bin = 1; bin <= hin_dx_BGsub_cutdy->GetNbinsX(); ++bin) {
    for (int bin = low_bin; bin <= high_bin; ++bin) {
        double content_data = minimizer_h_data_dx_with_BGsub->GetBinContent(bin);
        double content_expected = minimizer_h_simc_dx->GetBinContent(bin);
        // double content_expected = ( minimizer_scale_p * minimizer_h_simc_dx_p->GetBinContent(bin) ) + ( minimizer_scale_n * minimizer_h_simc_dx_n->GetBinContent(bin) );

        double error_data = minimizer_h_data_dx_with_BGsub->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    // Store the chi-squared value
    // (cl_dx_chi2_data_simc_BGsub_global.chiSquared_values).push_back(chiSquared);

    delete minimizer_h_data_dx_with_BGsub;
    delete minimizer_h_simc_dx_p;
    delete minimizer_h_simc_dx_n;
    delete minimizer_h_simc_dx;

    return chiSquared;
}

void MinimizationFunctionIn_CustomRange_ByScale_BGsub(int &npar, double *gin, double &result, double *param, int flag){
    // global_dyAnticut_scale_p = param[0];
    // global_dyAnticut_scale_n = param[1];

    cl_dx_chi2_data_simc_BGsub_global.global_scale_p = param[0];
    cl_dx_chi2_data_simc_BGsub_global.global_scale_n = param[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    // result = CalculateChiSquaredIn_CustomRange_ByScale_BGsub(cl_dx_chi2_data_simc_BGsub_global, -1.0, -0.7) 
    //         // + CalculateChiSquaredIn_CustomRange_ByScale_BGsub(cl_dx_chi2_data_simc_BGsub_global, -0.6, -0.2)
    //         + CalculateChiSquaredIn_CustomRange_ByScale_BGsub(cl_dx_chi2_data_simc_BGsub_global, -0.15, 0.2);

    double Xaxis_min = cl_dx_chi2_data_simc_BGsub_global.dx_Xaxis_min;
    double Xaxis_max = cl_dx_chi2_data_simc_BGsub_global.dx_Xaxis_max;

    result = CalculateChiSquaredIn_CustomRange_ByScale_BGsub(cl_dx_chi2_data_simc_BGsub_global, Xaxis_min, Xaxis_max);
}

double CalculateChiSquaredInRangeByScale_BGsub_SimpleCalc( TString BG_type = "" ) {

    cl_chi2_minimization cl_dx_to_minimize;
    
    bool inelastics = false;
    bool dyAnticut = false;

    if( BG_type == "inelastics" ){
        inelastics = true;
        dyAnticut = false;
        cl_dx_to_minimize = cl_dx_chi2_inelastics;
    }

    if( BG_type == "dyAnticut" ){
        inelastics = false;
        dyAnticut = true;
        cl_dx_to_minimize = cl_dx_chi2_dyAnticut;
    }

    double minimizer_scale_p = cl_dx_to_minimize.scale_p;
    double minimizer_scale_n = cl_dx_to_minimize.scale_n;
    TH1D *h_data_dx_with_BGsub= cl_dx_to_minimize.h_data_dx_with_BGsub;
    TH1D *h_simc_dx_p = (TH1D*)cl_dx_to_minimize.h_simc_dx_p_BeforeScaling->Clone("h_simc_dx_p");
    TH1D *h_simc_dx_n = (TH1D*)cl_dx_to_minimize.h_simc_dx_n_BeforeScaling->Clone("h_simc_dx_n");

    double chiSquared_low_x = cl_dx_to_minimize.chiSquared_low_x;
    double chiSquared_high_x = cl_dx_to_minimize.chiSquared_high_x; 

    // vector<vector<double>> *scale_n_p_values = cl_to_minimize.scale_n_p_values;

    double chiSquared = 0.0;

    h_simc_dx_p->Scale( minimizer_scale_p );
    h_simc_dx_n->Scale( minimizer_scale_n );

    TH1D *h_simc_dx = (TH1D*)h_simc_dx_p->Clone("h_simc_dx");
    h_simc_dx->Add(h_simc_dx_n, 1);

    chiSquared = SimpleCalculateChiSquared_InRange( h_data_dx_with_BGsub, h_simc_dx, chiSquared_low_x, chiSquared_high_x );

    // Store the chi-squared value
    cl_dx_to_minimize.chiSquared_values.push_back(chiSquared);
    cl_dx_to_minimize.scale_n_values.push_back(minimizer_scale_n);
    cl_dx_to_minimize.scale_p_values.push_back(minimizer_scale_p);

    cl_dx_to_minimize.scale_n = minimizer_scale_n;
    cl_dx_to_minimize.scale_p = minimizer_scale_p;

    // Store the neutron and proton scales
    // scale_n_p_values->push_back( {scale_n, scale_p} );

    cout << "Chi-squared output: " << endl;
    cout << "chi_squared: " << chiSquared << endl;
    cout << "n: " << minimizer_scale_n << ", p: " << minimizer_scale_p << endl;
    cout << "----------------------------------" << endl;

    if( BG_type == "inelastics" ){
        cl_dx_chi2_inelastics = cl_dx_to_minimize;
    }

    if( BG_type == "dyAnticut" ){
        cl_dx_chi2_dyAnticut = cl_dx_to_minimize;
    }

    return chiSquared;
}

void MinimizationFunctionInRangeByScale_dyAnticut_BGsub_SimpleCalc(int &npar, double *gin, double &result, double *param, int flag){
    // global_dyAnticut_scale_p = param[0];
    // global_dyAnticut_scale_n = param[1];

    cl_dx_chi2_dyAnticut.scale_p = param[0];
    cl_dx_chi2_dyAnticut.scale_n = param[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByScale_BGsub_SimpleCalc("dyAnticut");;
}

void MinimizationFunctionInRangeByScale_inelastics_BGsub_SimpleCalc(int &npar, double *gin, double &result, double *param, int flag){
    // global_inelastics_scale_p = param[0];
    // global_inelastics_scale_n = param[1];

    cl_dx_chi2_inelastics.scale_p = param[0];
    cl_dx_chi2_inelastics.scale_n = param[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByScale_BGsub_SimpleCalc("inelastics");;
}

// double CalculateChiSquaredInRangeByScale_dyAnticut_BGsub( cl_chi2_minimization& cl_dx_chi2_dyAnticut ) {
    
//     double minimizer_scale_p = cl_dx_chi2_dyAnticut.scale_p;
//     double minimizer_scale_n = cl_dx_chi2_dyAnticut.scale_n;
//     TH1D *h_data_dx_with_BGsub= cl_dx_chi2_dyAnticut.h_data_dx_with_BGsub;
//     TH1D *h_simc_dx_p = cl_dx_chi2_dyAnticut.h_simc_dx_p_BeforeScaling;
//     TH1D *h_simc_dx_n = cl_dx_chi2_dyAnticut.h_simc_dx_n_BeforeScaling;

//     double chiSquared_low_x = cl_dx_chi2_dyAnticut.chiSquared_low_x;
//     double chiSquared_high_x = cl_dx_chi2_dyAnticut.chiSquared_high_x; 

//     // vector<vector<double>> *scale_n_p_values = cl_to_minimize.scale_n_p_values;

//     double chiSquared = 0.0;

//     int low_bin = h_data_dx_with_BGsub->GetXaxis()->FindBin(chiSquared_low_x);
//     int high_bin = h_data_dx_with_BGsub->GetXaxis()->FindBin(chiSquared_high_x);

//     // for (int bin = 1; bin <= hin_dx_BGsub_cutdy->GetNbinsX(); ++bin) {
//     for (int bin = low_bin; bin <= high_bin; ++bin) {
//         double content_data = h_data_dx_with_BGsub->GetBinContent(bin);
//         double content_expected = ( minimizer_scale_p*h_simc_dx_p->GetBinContent(bin) ) + ( minimizer_scale_n*h_simc_dx_n->GetBinContent(bin) );

//         double error_data = h_data_dx_with_BGsub->GetBinError(bin);

//         if (error_data > 0.0) {
//             double residual = content_data - content_expected;
//             chiSquared += (residual * residual) / (error_data * error_data);
//         }
//     }

//     TH1D *h_simc_dx = (TH1D*)h_simc_dx_p->Clone("h_simc_dx");
//     h_simc_dx->Add(h_simc_dx_n, 1);

//     // Store the chi-squared value
//     (cl_dx_chi2_dyAnticut.chiSquared_values).push_back(chiSquared);
//     (cl_dx_chi2_dyAnticut.scale_n_values).push_back(minimizer_scale_n);
//     (cl_dx_chi2_dyAnticut.scale_p_values).push_back(minimizer_scale_p);

//     // Store the neutron and proton scales
//     // scale_n_p_values->push_back( {scale_n, scale_p} );

//     cout << "dyAnticut Chi-squared output: " << endl;
//     cout << "chi_squared: " << chiSquared << endl;
//     cout << "n: " << minimizer_scale_n << ", p: " << minimizer_scale_p << endl;
//     cout << "----------------------------------" << endl;

//     return chiSquared;
// }

// void MinimizationFunctionInRangeByScale_BGsub_dyAnticut(int &npar, double *gin, double &result, double *param, int flag){
//     // global_dyAnticut_scale_p = param[0];
//     // global_dyAnticut_scale_n = param[1];

//     cl_dx_chi2_dyAnticut.scale_p = param[0];
//     cl_dx_chi2_dyAnticut.scale_n = param[1];

//     //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
//     result = CalculateChiSquaredInRangeByScale_dyAnticut_BGsub(cl_dx_chi2_dyAnticut);
// }

// double CalculateChiSquaredInRangeByScale_inelastics_BGsub( cl_chi2_minimization& cl_dx_chi2_inelastics ) {
    
//     double minimizer_scale_p = cl_dx_chi2_inelastics.scale_p;
//     double minimizer_scale_n = cl_dx_chi2_inelastics.scale_n;
//     TH1D *h_data_dx_with_BGsub= cl_dx_chi2_inelastics.h_data_dx_with_BGsub;
//     TH1D *h_simc_dx_p = cl_dx_chi2_inelastics.h_simc_dx_p_BeforeScaling;
//     TH1D *h_simc_dx_n = cl_dx_chi2_inelastics.h_simc_dx_n_BeforeScaling;

//     double chiSquared_low_x = cl_dx_chi2_inelastics.chiSquared_low_x;
//     double chiSquared_high_x = cl_dx_chi2_inelastics.chiSquared_high_x; 

//     // vector<double> chiSquared_values = {};
//     // vector<double> scale_n_values = {};
//     // vector<double> scale_p_values = {};
//     // vector<vector<double>> *scale_n_p_values = cl_to_minimize.scale_n_p_values;

//     double chiSquared = 0.0;

//     int low_bin = h_data_dx_with_BGsub->GetXaxis()->FindBin(chiSquared_low_x);
//     int high_bin = h_data_dx_with_BGsub->GetXaxis()->FindBin(chiSquared_high_x);

//     // for (int bin = 1; bin <= hin_dx_BGsub_cutdy->GetNbinsX(); ++bin) {
//     for (int bin = low_bin; bin <= high_bin; ++bin) {
//         double content_data = h_data_dx_with_BGsub->GetBinContent(bin);
//         double content_expected = ( minimizer_scale_p*h_simc_dx_p->GetBinContent(bin) ) + ( minimizer_scale_n*h_simc_dx_n->GetBinContent(bin) );

//         double error_data = h_data_dx_with_BGsub->GetBinError(bin);

//         if (error_data > 0.0) {
//             double residual = content_data - content_expected;
//             chiSquared += (residual * residual) / (error_data * error_data);
//         }
//     }

//     TH1D *h_simc_dx = (TH1D*)h_simc_dx_p->Clone("h_simc_dx");
//     h_simc_dx->Add(h_simc_dx_n, 1);

//     // Store the chi-squared value
//     (cl_dx_chi2_inelastics.chiSquared_values).push_back(chiSquared);
//     (cl_dx_chi2_inelastics.scale_n_values).push_back(minimizer_scale_n);
//     (cl_dx_chi2_inelastics.scale_p_values).push_back(minimizer_scale_p);

//     // Store the neutron and proton scales
//     // scale_n_p_values->push_back( {scale_n, scale_p} );

//     cout << "inelastics Chi-squared output: " << endl;
//     cout << "chi_squared: " << chiSquared << endl;
//     cout << "n: " << minimizer_scale_n << ", p: " << minimizer_scale_p << endl;
//     cout << "----------------------------------" << endl;

//     return chiSquared;
// }

// void MinimizationFunctionInRangeByScale_BGsub_inelastics(int &npar, double *gin, double &result, double *param, int flag){
//     // global_inelastics_scale_p = param[0];
//     // global_inelastics_scale_n = param[1];

//     cl_dx_chi2_inelastics.scale_p = param[0];
//     cl_dx_chi2_inelastics.scale_n = param[1];

//     //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
//     result = CalculateChiSquaredInRangeByScale_inelastics_BGsub(cl_dx_chi2_inelastics);
// }

double CalculateChiSquaredInRangeByScale_ALL_BGsub( cl_chi2_minimization& cl_dx_chi2_BGsub_final ) {
    
    double minimizer_scale_p = cl_dx_chi2_BGsub_final.scale_p;
    double minimizer_scale_n = cl_dx_chi2_BGsub_final.scale_n;
    TH1D *h_data_dx_with_BGsub = (TH1D*)cl_dx_chi2_BGsub_final.h_data_dx_with_BGsub->Clone("h_data_dx_with_BGsub");
    TH1D *h_simc_dx_p = (TH1D*)cl_dx_chi2_BGsub_final.h_simc_dx_p_BeforeScaling->Clone("h_simc_dx_p");
    TH1D *h_simc_dx_n = (TH1D*)cl_dx_chi2_BGsub_final.h_simc_dx_n_BeforeScaling->Clone("h_simc_dx_n");

    double chiSquared_low_x = cl_dx_chi2_BGsub_final.chiSquared_low_x;
    double chiSquared_high_x = cl_dx_chi2_BGsub_final.chiSquared_high_x; 

    // vector<double> chiSquared_values = {};
    // vector<double> scale_n_values = {};
    // vector<double> scale_p_values = {};
    // vector<vector<double>> *scale_n_p_values = cl_to_minimize.scale_n_p_values;

    double chiSquared = 0.0;

    int low_bin = h_data_dx_with_BGsub->GetXaxis()->FindBin(chiSquared_low_x);
    int high_bin = h_data_dx_with_BGsub->GetXaxis()->FindBin(chiSquared_high_x);

    // for (int bin = 1; bin <= hin_dx_BGsub_cutdy->GetNbinsX(); ++bin) {
    for (int bin = low_bin; bin <= high_bin; ++bin) {
        double content_data = h_data_dx_with_BGsub->GetBinContent(bin);
        double content_expected = ( minimizer_scale_p*h_simc_dx_p->GetBinContent(bin) ) + ( minimizer_scale_n*h_simc_dx_n->GetBinContent(bin) );

        double error_data = h_data_dx_with_BGsub->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    TH1D *h_simc_dx = (TH1D*)h_simc_dx_p->Clone("h_simc_dx");
    h_simc_dx->Add(h_simc_dx_n, 1);

    // Store the chi-squared value
    (cl_dx_chi2_BGsub_final.chiSquared_values).push_back(chiSquared);
    (cl_dx_chi2_BGsub_final.scale_n_values).push_back(minimizer_scale_n);
    (cl_dx_chi2_BGsub_final.scale_p_values).push_back(minimizer_scale_p);

    // Store the neutron and proton scales
    // scale_n_p_values->push_back( {scale_n, scale_p} );

    cout << "inelastics Chi-squared output: " << endl;
    cout << "chi_squared: " << chiSquared << endl;
    cout << "n: " << minimizer_scale_n << ", p: " << minimizer_scale_p << endl;
    cout << "----------------------------------" << endl;

    return chiSquared;
}

void MinimizationFunctionInRangeByScale_BGsub_ALL(int &npar, double *gin, double &result, double *param, int flag){
    // global_inelastics_scale_p = param[0];
    // global_inelastics_scale_n = param[1];

    cl_dx_chi2_BGsub_final.scale_p = param[0];
    cl_dx_chi2_BGsub_final.scale_n = param[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByScale_ALL_BGsub(cl_dx_chi2_BGsub_final);
}

double CalculateChiSquaredInRangeByScaleWithNormalization(const cl_chi2_minimization& cl_dx_chi2_data_simc) {
    
    double minimizer_scale_p = cl_dx_chi2_data_simc.scale_p;
    double minimizer_scale_n = cl_dx_chi2_data_simc.scale_n;
    TH1D *hin_dx_cutdy = cl_dx_chi2_data_simc.h_data_dx;
    TH1D *h_simc_dx_p = cl_dx_chi2_data_simc.h_simc_dx_p;
    TH1D *h_simc_dx_n = cl_dx_chi2_data_simc.h_simc_dx_n;
    double ChiSquared_low_x = cl_dx_chi2_data_simc.chiSquared_low_x;
    double ChiSquared_high_x = cl_dx_chi2_data_simc.chiSquared_high_x; 
    vector<double> chiSquared_values = cl_dx_chi2_data_simc.chiSquared_values;
    double chiSquared = 0.0;

    int low_bin = hin_dx_cutdy->GetXaxis()->FindBin(ChiSquared_low_x);
    int high_bin = hin_dx_cutdy->GetXaxis()->FindBin(ChiSquared_high_x);

    // for (int bin = 1; bin <= hin_dx_cutdy->GetNbinsX(); ++bin) {

    TH1D *h_simc_dx_added = new TH1D("h_simc_dx_added", "h_simc_dx_added", h_simc_dx_p->GetNbinsX(), h_simc_dx_p->GetXaxis()->GetXmin(), h_simc_dx_p->GetXaxis()->GetXmax());
    
    for( int i = 1; i <= h_simc_dx_added->GetNbinsX(); i++){
        h_simc_dx_added->SetBinContent(i, minimizer_scale_p*h_simc_dx_p->GetBinContent(i) + minimizer_scale_n*h_simc_dx_n->GetBinContent(i) );
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
    chiSquared_values.push_back(chiSquared);

    return chiSquared;
}

void MinimizationFunctionInRangeByScaleWithNormalization(int &npar, double *gin, double &result, double *param, int flag){
    // global_scale_p = param[0];
    // global_scale_n = param[1];
    cl_dx_chi2_data_simc.scale_p = param[0];
    cl_dx_chi2_data_simc.scale_n = param[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByScaleWithNormalization(cl_dx_chi2_data_simc);
}

double CalculateChiSquaredInRangeByDyCut(const cl_chi2_minimization& cl_dx_chi2_data_simc) {
    TH2D *h2_dxdy_trough = cl_dx_chi2_data_simc.h2_chi2_dxdy_cutdy;
    int h2_dxdy_xbins = h2_dxdy_trough->GetNbinsX();
    int h2_dxdy_ybins = h2_dxdy_trough->GetNbinsY();

    double dy_cut_min = cl_dx_chi2_data_simc.dy_cut_min;
    int dy_cut_min_bin = h2_dxdy_trough->GetXaxis()->FindBin(dy_cut_min);

    double dy_cut_max = cl_dx_chi2_data_simc.dy_cut_max;
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

    double h_sim_trough_integral = cl_dx_chi2_data_simc.simc_dx_trough_integral;

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
    cl_dx_chi2_data_simc.dy_cut_min = parDyCut[0];
    cout << "par[0]: " << parDyCut[0] << ", par[1]: " << parDyCut[1] << endl;
    cl_dx_chi2_data_simc.dy_cut_max = parDyCut[1];

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeByDyCut(cl_dx_chi2_data_simc);
}

double CalculateChiSquaredInRangeDyCutSigmaMult(const cl_chi2_minimization& cl_dx_chi2_data_simc) {
    TH2D *h2_dxdy_SigmaMult_trough = (TH2D*)cl_dx_chi2_data_simc.h2_chi2_dxdy_cutdy_SigmaMult->Clone("h2_dxdy_SigmaMult_trough");
    int h2_dxdy_SigmaMult_xbins = h2_dxdy_SigmaMult_trough->GetNbinsX();
    int h2_dxdy_SigmaMult_ybins = h2_dxdy_SigmaMult_trough->GetNbinsY();

    double dy_cut_minimizer_sigma = cl_dx_chi2_data_simc.dy_cut_pn_sigma;
    double dy_cut_minimizer_pn = cl_dx_chi2_data_simc.dy_cut_pn;
    double dy_cut_minimizer_mult = cl_dx_chi2_data_simc.dy_cut_sigmaMult;

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

    double h_sim_trough_integral = cl_dx_chi2_data_simc.simc_dx_trough_integral;

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
    cl_dx_chi2_data_simc.dy_cut_sigmaMult = parSigmaMult[0];
    // cout << "Sigma Mult: " << parSigmaMult[0] << endl;

    //BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
    result = CalculateChiSquaredInRangeDyCutSigmaMult(cl_dx_chi2_data_simc);
}


double CalculateChiSquared(const cl_chi2_minimization& cl_dx_chi2_data_simc) {
    
	double minimizer_scale_p = cl_dx_chi2_data_simc.scale_p;
	double minimizer_scale_n = cl_dx_chi2_data_simc.scale_n;
	TH1D *hin_dx_cutdy = cl_dx_chi2_data_simc.h_data_dx;
	TH1D *h_simc_dx_p = cl_dx_chi2_data_simc.h_simc_dx_p;
	TH1D *h_simc_dx_n = cl_dx_chi2_data_simc.h_simc_dx_n;
	double ChiSquared_low_x = cl_dx_chi2_data_simc.chiSquared_low_x;
	double ChiSquared_high_x = cl_dx_chi2_data_simc.chiSquared_high_x; 
	vector<double> chiSquared_values = cl_dx_chi2_data_simc.chiSquared_values;
    double chiSquared = 0.0;

    for (int bin = 1; bin <= hin_dx_cutdy->GetNbinsX(); ++bin) {

        double content_data = hin_dx_cutdy->GetBinContent(bin);
        double content_expected = minimizer_scale_p*h_simc_dx_p->GetBinContent(bin) + minimizer_scale_n*h_simc_dx_n->GetBinContent(bin);

        double error_data = hin_dx_cutdy->GetBinError(bin);

        if (error_data > 0.0) {
            double residual = content_data - content_expected;
            chiSquared += (residual * residual) / (error_data * error_data);
        }
    }

    // Store the chi-squared value
    chiSquared_values.push_back(chiSquared);

    return chiSquared;
}

void MinimizationFunction(int &npar, double *gin, double &result, double *param, int flag){
	// global_scale_p = param[0];
	// global_scale_n = param[1];

	cl_dx_chi2_data_simc.scale_p = param[0];
	cl_dx_chi2_data_simc.scale_n = param[1];

	// result = CalculateChiSquared(global_scale_p, global_scale_n);
	//BEFORE THIS ACTIVATION YOU MUST HAVE CREATED A CLASS CALLED: cl_dx_chi2_data_simc
	result = CalculateChiSquared(cl_dx_chi2_data_simc);
}



void collect_chiSquared_values( const cl_chi2_minimization& cl_dx_chi2_data_simc, kine_ff_extract& SBSkine ){

    global_BGsub_scale_p = cl_dx_chi2_data_simc.scale_p;
    SBSkine.BGsub_scale_p = cl_dx_chi2_data_simc.scale_p;

    global_BGsub_scale_p_error = cl_dx_chi2_data_simc.scale_p_error;
    SBSkine.BGsub_scale_p_error = cl_dx_chi2_data_simc.scale_p_error;

    global_BGsub_scale_n = cl_dx_chi2_data_simc.scale_n;
    SBSkine.BGsub_scale_n = cl_dx_chi2_data_simc.scale_n;

    global_BGsub_scale_n_error = cl_dx_chi2_data_simc.scale_n_error;
    SBSkine.BGsub_scale_n_error = cl_dx_chi2_data_simc.scale_n_error;

    SBSkine.BGsub_scale_np_ratio = global_BGsub_scale_n/global_BGsub_scale_p;

}



#endif