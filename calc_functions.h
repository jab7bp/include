#ifndef calc_functions
#define calc_functions

// double calc_mean(double *arr_for_mean){
// 	int n = (sizeof(*(&arr_for_mean))/sizeof(arr_for_mean[0]));

// 	double sum = 0.0, mean;

// 	for(int i = 0; i < n; i++){
// 		sum += arr_for_mean[i];
// 	}
// 	mean = sum/n;
// 	return mean;
// }

// double calc_StdDev(double *arr_for_stdDev){
// 	int n = sizeof(&arr_for_stdDev)/sizeof(&arr_for_stdDev[0]);
// 	cout << "n: " << n << endl;

// 	return 0;
// }

template<typename T>
double VectorMean(std::vector<T> const& v){
	if(v.empty()){
		return 0;
	}
	return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

// Function to calculate the standard deviation of a vector
double VectorStandardDeviation(const std::vector<double>& data) {
    int n = data.size();
    if (n <= 1) {
        std::cerr << "Error: Vector size must be greater than 1 for standard deviation calculation." << std::endl;
        return 0.0; // or handle the error as appropriate for your application
    }

    // Calculate mean
    double mean = TMath::Mean(n, &data[0]);

    // Calculate sum of squared differences
    double sumSquaredDiff = 0.0;
    for (int i = 0; i < n; ++i) {
        double diff = data[i] - mean;
        sumSquaredDiff += diff * diff;
    }

    // Calculate standard deviation
    double standardDeviation = TMath::Sqrt(sumSquaredDiff / (n - 1));

    return standardDeviation;
}

double VectorRootSumOfSquares( const std::vector<double>& data ){
	int n = data.size();
	double RSS = 0.0;

	for( int i = 0; i < n; i++ ){
		RSS += pow( data[i], 2);
	}

	return sqrt(RSS)/n;
}

double VectorStandardErrorOfMean( const std::vector<double>& data ){

	int n = data.size();

	double sum = 0.0;

	for( int i = 0; i < n; i++ ){
		sum += pow( data[i], 2);
	}

	double SEM = sum/(sqrt(n));
	return SEM;
}

double calc_gaus_error(double par_0, double par_0_err, double par_2, double par_2_err, double bin_width ){
	double gaus_error = 0.0;

	gaus_error = ( sqrt( TMath::Pi()*( pow(par_2, 2)*pow(par_0_err, 2) + pow(par_0, 2)*pow(par_2_err, 2) ) ) )/bin_width;

	return gaus_error;
}

void list_files(TString directory, vector<TString> &filenames, TString ext){
	const char *dirname = directory.Data();

	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();

	if(files){
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while( (file = (TSystemFile*)next() ) ){
			fname = file->GetName();
			if( !file->IsDirectory() && fname.BeginsWith( ext.Data() )) {
				filenames.push_back( Form("%s/%s", directory.Data(), fname.Data()) );
				// cout << fname.Data() << endl;
			}
		}
	}
}

double calc_luminosity( double I_beam, TString run_targ ){

	double luminosity = 0.0;

	double e_charge = 1.60E-19;
	double N_avogadro = 6.02E23;

	double targ_len = 1.50E01;
	double dens_LH2 = 7.23E-02;
	double mol_m_LH2 = 2.02E00;
	double atom_dens_LH2 = (dens_LH2/mol_m_LH2)*(N_avogadro)*2.0;

	double dens_LD2 = 1.67E-01;
	double mol_m_LD2 = 4.03E00;
	double atom_dens_LD2 = (dens_LD2/mol_m_LD2)*(N_avogadro)*2.0;

//pdbforce: (I_bamn/qe)*len_tar*ld2_dens*(N_A/LD2_mass)
	double lumi_LH2 = (I_beam/e_charge)*targ_len*dens_LH2*(N_avogadro/mol_m_LH2)*(1E-6);
	double lumi_LD2 = (I_beam/e_charge)*targ_len*dens_LD2*(N_avogadro/mol_m_LD2)*(1E-6);

//adr: (a_dens_LH2*tar_len/qe)*(E-6)*I_beam
	// double lumi_LH2 = (atom_dens_LH2*targ_len/e_charge)*(1E-6)*I_beam;
	// double lumi_LD2 = (atom_dens_LD2*targ_len/e_charge)*(1E-6)*I_beam;

	if( run_targ == "LH2" ){
		return lumi_LH2;
	}

	//run_targ == "LD2"
	else{
		return lumi_LD2;
	}

}

double CalculateChiSquaredTwoHistograms( TH1D *h_observed, TH1D *h_expected ){
	if (h_observed->GetNbinsX() != h_expected->GetNbinsX()) {
        std::cerr << "Histograms have different numbers of bins!" << std::endl;
        return 0;
    }

    double chiSquared = 0.0;
    double Ndf = h_observed->GetNbinsX() - 2;

    for (int i = 1; i <= h_observed->GetNbinsX(); ++i) {
        double Oi = h_observed->GetBinContent(i);
        double Ei = h_expected->GetBinContent(i);

        if (Ei > 0) { // Avoid division by zero
            double term = (Oi - Ei) * (Oi - Ei) / Ei;
            chiSquared += term;
        }
    }

    return chiSquared/Ndf;
}

double CalculateChiSquaredInRangeTwoHistograms( TH1D *h_observed, TH1D *h_expected, double low_x, double high_x){
	if (h_observed->GetNbinsX() != h_expected->GetNbinsX()) {
        std::cerr << "Histograms have different numbers of bins!" << std::endl;
        return 0;
    }

    int low_bin = h_observed->GetXaxis()->FindBin(low_x);
    int high_bin = h_observed->GetXaxis()->FindBin(high_x);

    double chiSquared = 0.0;
    double Ndf = h_observed->GetNbinsX() - 2;

    for (int i = low_bin; i <= high_bin; ++i) {
        double Oi = h_observed->GetBinContent(i);
        double Ei = h_expected->GetBinContent(i);

        if (Ei > 0) { // Avoid division by zero
            double term = (Oi - Ei) * (Oi - Ei) / Ei;
            chiSquared += term;
        }
    }

    return chiSquared/Ndf;
}

double CalculateChiSquaredWithErrors(TH1* observed, TH1* expected) {
    if (observed->GetNbinsX() != expected->GetNbinsX()) {
        std::cerr << "Histograms have different numbers of bins!" << std::endl;
        return 0 ;
    }

    double chiSquared = 0.0;

    for (int i = 1; i <= observed->GetNbinsX(); ++i) {
        double Oi = observed->GetBinContent(i);
        double Ei = expected->GetBinContent(i);
        double sigma_i = observed->GetBinError(i);
        double tau_i = expected->GetBinError(i);

        double term = (Oi - Ei) * (Oi - Ei) / (sigma_i * sigma_i + tau_i * tau_i);
        chiSquared += term;
    }

    return chiSquared;
}

double calculateHCalDetEff(int kine, TString nucleon, bool return_error = false ){

	double hcal_det_eff = 0.0;
	double hcal_det_eff_err = 0.0;

	double p_nucl;
	if( kine == 8 ){ p_nucl = 3.22; }
	if( kine == 9 ){ p_nucl = 3.21; }

	double a0 = 85.689454;
	double a0_err = 0.41106;

	double a1 = 7.5699056;
	double a1_err = 0.421948;

	double a2 = -1.8512754;
	double a2_err = 0.145954;

	double a3 = 0.19593452;
	double a3_err = 0.0206016;

	double a4 = -0.007713635;
	double a4_err = 0.00101772;

//------------------
	double b0 = 79.3477;
	double b0_err = 0.733834;

	double b1 = 13.8929;
	double b1_err = 0.796093;

	double b2 = -3.60264;
	double b2_err = 0.284625;

	double b3 = 0.386701;
	double b3_err = 0.0408757;

	double b4 = -0.0152329;
	double b4_err = 0.00203378;

	if( nucleon == "n" ){
		hcal_det_eff = a0 + (a1*p_nucl) + (a2*pow(p_nucl, 2)) + (a3*pow(p_nucl, 3) ) + (a4*pow(p_nucl, 4));
		hcal_det_eff_err = sqrt( pow(a0_err, 2) + pow( p_nucl*a1_err, 2) + pow( pow(p_nucl, 2)*a2_err, 2) + pow( pow(p_nucl, 3)*a3_err, 3) + pow( pow(p_nucl, 4)*a4_err, 4));
	}

	if( nucleon == "p" ){
		hcal_det_eff = b0 + (b1*p_nucl) + (b2*pow(p_nucl, 2)) + (b3*pow(p_nucl, 3) ) + (b4*pow(p_nucl, 4));
		hcal_det_eff_err = sqrt( pow(b0_err, 2) + pow( p_nucl*b1_err, 2) + pow( pow(p_nucl, 2)*b2_err, 2) + pow( pow(p_nucl, 3)*b3_err, 3) + pow( pow(p_nucl, 4)*b4_err, 4));
	}
	
	if( return_error ){
		return hcal_det_eff_err/100.0;
	}

	//!return_error
	else{
		return hcal_det_eff/100.0;		
	}

}

int FindMinimumDifferenceBinInRange(TH1D* hist1, TH1D* hist2, double xMin, double xMax) {
    int numBins = hist1->GetNbinsX();
    int binXMin = hist1->GetXaxis()->FindBin(xMin);
    int binXMax = hist1->GetXaxis()->FindBin(xMax);

    double minDifference = std::numeric_limits<double>::max();
    int minBin = -1;

    for (int bin = binXMin; bin <= binXMax; ++bin) {
        double difference = std::abs(hist1->GetBinContent(bin) - hist2->GetBinContent(bin));
        if (difference < minDifference) {
            minDifference = difference;
            minBin = bin;
        }
    }

    return minBin; // Return the bin index that corresponds to the minimum difference
}


#endif