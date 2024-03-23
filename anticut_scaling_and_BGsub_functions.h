#ifndef ANTICUT_SCALING_AND_BGSUB_FUNCTIONS_H
#define ANTICUT_SCALING_AND_BGSUB_FUNCTIONS_H

	// if( kine == -1 || sbsfieldscale == -1 || run_target == "" ){
	// 	cout << "------- Error -------" << endl;
	// 	cout << "Function make_dyAnticut_histo requires valid inputs for kine, sbsfieldscale, and run_target" << endl;
	// 	exit(1);
	// }

	// const int total_parN = BG_polN + 6;
	// double parTotalDx[total_parN];

	// vector<double> parPolN_vec = {};

	// if( BG_polN == 0 ){
	// 	cout << "No value specified for BG_polN for dy anticut BG fit. Defaulting to: pol2" << endl;
	// 	BG_polN = 4;
	// 	SBSkine.dx_dyAnticut_BG_polN = BG_polN;
	// }

TF1 *fit_and_get_inelastic_BG_fit_from_file( cl_SBSkine& SBSkine, bool anticut_reject = false, bool inelastic_BG_wcut = false, double fit_x_min = -2.0, double fit_x_max = 1.0){

	int kine = SBSkine.kine;
	int sbsfieldscale = SBSkine.sbsfieldscale;
	TString magmod = SBSkine.magmod;
	TString file_date = "";

	Double_t dy_mean = SBSkine.dy_mean;
	Double_t dy_sigma = SBSkine.dy_sigma;

	double dy_mult = SBSkine.dy_mult;
	double dy_mult_minus = SBSkine.dy_mult_minus;
	double dy_mult_plus = SBSkine.dy_mult_plus;

	Double_t dy_min, dy_max;
	if( SBSkine.pass == 1 ){
		dy_min = dy_mean - (dy_mult)*dy_sigma;
		dy_max = dy_mean + (dy_mult)*dy_sigma;		
	}
	else{
		dy_min = dy_mean - (dy_mult_minus)*dy_sigma;
		dy_max = dy_mean + (dy_mult_plus)*dy_sigma;		
	}


	TString I_beam_str = "";

	if( kine == 4 ){
		I_beam_str = "0175";
		file_date = "21_12_2023";
	}
	if( kine == 8 ){
		I_beam_str = "0500";
		file_date = "21_12_2023";
	}

	if( kine == 9 ){
		I_beam_str = "1200";
		file_date = "06_01_2024";
	}

	TString rootfile_dir = "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/MC/rootfiles";
	TString filename = Form("MC_BACKGROUND_SBS%i_LD2_mag%imod%s_%suA_dxdy_%s.root", kine, sbsfieldscale, magmod.Data(), I_beam_str.Data(), file_date.Data());

	TFile *BG_infile = new TFile(Form("%s/%s", rootfile_dir.Data(), filename.Data()), "READ" );

	TString h_BG_dx_name = "";
	TH1D *h_BG_dx;
	TH2D *h_BG_dxdy;

	TF1 *tf_dx_inelasticsBG_fit = new TF1("tf_dx_inelasticsBG_fit", "pol4", fit_x_min, fit_x_max);

	if( anticut_reject == false && inelastic_BG_wcut == false ){
		h_BG_dx_name = "h_BG_dx";
		h_BG_dx = (TH1D*)(BG_infile->Get("h_BG_dx"))->Clone("h_BG_dx");
		h_BG_dx->GetXaxis()->SetRangeUser(fit_x_min, fit_x_max);

		for( int par = 0; par < 4; par++ ){
			tf_dx_inelasticsBG_fit->SetParName( par, Form("inelastic_BG_dx_p%i", par));
		}

		// h_BG_dx->Draw();
		h_BG_dx->Fit("tf_dx_inelasticsBG_fit", "RMSE0+");
	}

	if( anticut_reject == true ){
		h_BG_dxdy = (TH2D*)(BG_infile->Get("h_BG_dxdy"))->Clone("h_BG_dxdy");
		h_BG_dxdy->GetXaxis()->SetRangeUser(dy_min, dy_max);
		h_BG_dx = (TH1D*)h_BG_dxdy->ProjectionY();
		h_BG_dx->GetXaxis()->SetRangeUser(fit_x_min, fit_x_max);

		h_BG_dx->SetName("h_BG_dx");
		h_BG_dx->SetTitle("h_BG_dx_with_dyAnticut");

		for( int par = 0; par < 4; par++ ){
			tf_dx_inelasticsBG_fit->SetParName( par, Form("inelastic_BG_with_dyAnticut_dx_p%i", par));
		}

		// h_BG_dx->Draw();
		h_BG_dx->Fit("tf_dx_inelasticsBG_fit", "RMSE0+");
	}

	cout << "Finished extracting inelastic fit from BG file: " << endl;
	cout << filename.Data() << endl;

	return tf_dx_inelasticsBG_fit;

}

vector<double> get_pN_polN_vec_from_TF1_fit_to_file( cl_SBSkine& SBSkine, bool anticut_reject = false, bool inelastic_BG_wcut = false, double fit_x_min = -2.0, double fit_x_max = 1.0){

	TF1 *tf_bg_fit = fit_and_get_inelastic_BG_fit_from_file( SBSkine, anticut_reject, inelastic_BG_wcut, fit_x_min, fit_x_max);
	int polN = tf_bg_fit->GetNpar();

	vector<double> pN_polN_vec = {};

	for( int par = 0; par < polN; par++ ){
		pN_polN_vec.push_back( tf_bg_fit->GetParameter(par) );
	}

	return pN_polN_vec;
}

vector<double> get_pN_polN_vec( int kine = 8, bool anticut_reject = false, bool inelastic_BG_wcut = false ){


	vector<double> pN_polN_vec = {};

	if( kine == 4 ){
	//STANDARD
		if( !anticut_reject ){
			pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};		

			if( inelastic_BG_wcut ){
				// Dec 22, 2023
				pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};				
			}
		}		
	//Anti-cut  REJECT SBS8
		if( anticut_reject ){
			pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};				

			if( inelastic_BG_wcut ){
				// Dec 22, 2023
				pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};				
			}
		}
	}

	if( kine == 8 ){
	//STANDARD
		if( !anticut_reject ){
			pN_polN_vec = {3.60894e-06, -6.05927e-07, -7.98384e-07, -1.28239e-07, -4.61731e-08};

			if( inelastic_BG_wcut ){
				pN_polN_vec = {1.71886e-08, -1.22842e-08, -8.24278e-09, 4.47680e-09, 1.95109e-09};	
				// Dec 22, 2023		
			}
		}
	//Anti-cut  REJECT SBS8
		if( anticut_reject ){
			pN_polN_vec = {1.94233e-07, -1.14126e-07, -8.25939e-08, 4.49891e-08, 2.18278e-08};				

			if( inelastic_BG_wcut ){
				pN_polN_vec = {1.71886e-08, -1.22842e-08, -8.24278e-09, 4.47680e-09, 1.95109e-09};	
				// Dec 22, 2023			
			}
		}
	}
	if( kine == 9 ){
	//Standard SBS9
		if( !anticut_reject ){
			pN_polN_vec = {4.52486e-06, -4.66733e-07, -9.19219e-07, -4.20985e-07, -1.74717e-07};
			// Dec 22, 2023		
			
			if( inelastic_BG_wcut ){
				pN_polN_vec = {9.52164e-08, -6.16508e-08, -4.70501e-08, 2.08224e-08, 9.76397e-09};	
				// Dec 22, 2023		
			}

		}
	
	//Anti-cut  REJECT	SBS9
		if( anticut_reject ){
			pN_polN_vec = {7.42661e-07, -2.95458e-07, -2.85281e-07, 5.32731e-08, 3.51582e-08};	
			// Dec 22, 2023

			if( inelastic_BG_wcut ){
				pN_polN_vec = {9.52164e-08, 9.52164e-08, -4.70501e-08, 2.08224e-08, 9.76397e-09};	
				// Dec 22, 2023			
			}					
		}
	}

	return pN_polN_vec;
}

TH1D *subtract_histograms( TH1D *h_primary, TH1D *h_to_sub, TString cust_histo_name = "", TString cust_histo_ID = "" ){

	int n_bins_primary = h_primary->GetNbinsX();
	int n_bins_secondary = h_to_sub->GetNbinsX();

	if( n_bins_primary != n_bins_secondary ){
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
		cout << "To take the difference between two histograms they must have the same nBinsX" << endl;
		cout << "N bins - h_primary: " << n_bins_primary << ", h_to_sub: " << n_bins_secondary << endl;
		cout << "h_primary name: " << h_primary->GetName() << endl;
		cout << "h_to_sub name: " << h_to_sub->GetName() << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	double histo_x_min = h_primary->GetXaxis()->GetXmin();
	double histo_x_max = h_primary->GetXaxis()->GetXmax();

	TString histo_name = "h_difference";
	if( cust_histo_name != "" ){
		histo_name = cust_histo_name.Data();
	}

	TString histo_ID = h_primary->GetName();
	if( cust_histo_ID != "" ){
		histo_ID = cust_histo_ID.Data();
	}

	TH1D *h_difference = new TH1D(histo_name.Data(), histo_ID.Data(), n_bins_primary, histo_x_min, histo_x_max );

	for( int bin = 1; bin <= n_bins_primary; bin++ ){
		h_difference->SetBinContent(bin, h_primary->GetBinContent(bin) - h_to_sub->GetBinContent(bin) );
	}	

	return h_difference;

}

TH1D *subtract_histograms_positive_values_only( TH1D *h_primary, TH1D *h_to_sub, TString cust_histo_name = "", TString cust_histo_ID = "" ){

	int n_bins_primary = h_primary->GetNbinsX();
	int n_bins_secondary = h_to_sub->GetNbinsX();

	if( n_bins_primary != n_bins_secondary ){
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
		cout << "To take the difference between two histograms they must have the same nBinsX" << endl;
		cout << "N bins - h_primary: " << n_bins_primary << ", h_to_sub: " << n_bins_secondary << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	double histo_x_min = h_primary->GetXaxis()->GetXmin();
	double histo_x_max = h_primary->GetXaxis()->GetXmax();

	TString histo_name = "h_difference";
	if( cust_histo_name != "" ){
		histo_name = cust_histo_name.Data();
	}

	TString histo_ID = h_primary->GetName();
	if( cust_histo_ID != "" ){
		histo_ID = cust_histo_ID.Data();
	}

	TH1D *h_pos_difference = new TH1D(histo_name.Data(), histo_ID.Data(), n_bins_primary, histo_x_min, histo_x_max );

	double difference = 0.0;

	for( int bin = 1; bin <= n_bins_primary; bin++ ){
		difference = 0.0;
		difference = h_primary->GetBinContent(bin) - h_to_sub->GetBinContent(bin);
		if( difference >= 0.0 ){
			h_pos_difference->SetBinContent(bin, h_primary->GetBinContent(bin) - h_to_sub->GetBinContent(bin) );
		}
	}	

	return h_pos_difference;

}

//Given a dx plot this function will fit the total dx: p_gaus + n_gaus + BG_pol4
TF1 *get_dx_total_fit( TH1D *h_dx, cl_SBSkine& SBSkine, TString identifier = "", int BG_polN = 4 ){

	if( identifier != "" ){ identifier.Prepend("_"); }

	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	const int total_parN = BG_polN + 7;
	double totalDx_par[total_parN];

	// VALUES FOR THE TOTAL DX FIT FOR H_DX
	Double_t dx_p_mean = SBSkine.dx_p_mean;
	Double_t dx_p_sigma = SBSkine.dx_p_sigma;
	Double_t dx_n_mean = SBSkine.dx_n_mean;
	Double_t dx_n_sigma = SBSkine.dx_n_sigma;
	Double_t dy_mean = SBSkine.dy_mean;
	Double_t dy_sigma = SBSkine.dy_sigma;

	double dy_mult = SBSkine.dy_mult;
	double dy_mult_minus = SBSkine.dy_mult_minus;
	double dy_mult_plus = SBSkine.dy_mult_plus;

	Double_t dy_min, dy_max;
	if( SBSkine.pass == 1 ){
		dy_min = dy_mean - (dy_mult)*dy_sigma;
		dy_max = dy_mean + (dy_mult)*dy_sigma;		
	}
	else{
		dy_min = dy_mean - (dy_mult_minus)*dy_sigma;
		dy_max = dy_mean + (dy_mult_plus)*dy_sigma;		
	}

	TAxis *xaxis = h_dx->GetXaxis();
	int xbins = xaxis->GetNbins();
	Double_t xaxis_min = xaxis->GetXmin();
	Double_t xaxis_max = xaxis->GetXmax();

//Let's draw the histogram and do the total fit with p_gaus + n_gaus + BG_polN:
	if( false ){
		TCanvas *c_get_total_dx = new TCanvas( Form("c_get_total_dx%s", identifier.Data()), Form("c_get_total_dx%s", identifier.Data()), 600, 500);
		h_dx->Draw();		
	}


	TF1 *tf_get_total_dx_fit = new TF1("tf_get_total_dx_fit", fit_full_dx_peaks_with_pol4BG, plot_xmin, plot_xmax, BG_polN + 7 );
	tf_get_total_dx_fit->SetNpx( h_dx->GetNbinsX() );

//Name the parameters
	tf_get_total_dx_fit->SetParName(0, Form("get_total_dx%s_p_norm", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(1, Form("get_total_dx%s_p_mean", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(2, Form("get_total_dx%s_p_sigma", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(3, Form("get_total_dx%s_n_norm", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(4, Form("get_total_dx%s_n_mean", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(5, Form("get_total_dx%s_n_sigma", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(6, Form("get_total_dx%s_BG_pol_p0", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(7, Form("get_total_dx%s_BG_pol_p1", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(8, Form("get_total_dx%s_BG_pol_p2", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(9, Form("get_total_dx%s_BG_pol_p3", identifier.Data() ));
	tf_get_total_dx_fit->SetParName(10, Form("get_total_dx%s_BG_pol_p4", identifier.Data() ));	

	cout << "Finished naming get_total total fit parameters. " << endl;
	cout << "Initializing tf_get_total_dx_fit parLimits....";

	tf_get_total_dx_fit->SetParLimits(0, 0.8*h_dx->GetMaximum(), 1.1*h_dx->GetMaximum());
	tf_get_total_dx_fit->SetParLimits(1, 1.1*dx_p_mean, 0.9*dx_p_mean);
	tf_get_total_dx_fit->SetParLimits(2, 0.7*dx_p_sigma, 1.1*dx_p_sigma);
	tf_get_total_dx_fit->SetParLimits(3, 0.25*h_dx->GetMaximum(), 0.6*h_dx->GetMaximum());
	tf_get_total_dx_fit->SetParLimits(4, dx_n_mean - 0.08, dx_n_mean + 0.08);
	tf_get_total_dx_fit->SetParLimits(5, 0.8*dx_n_sigma, 1.25*dx_n_sigma);	
	// tf_get_total_dx_fit->SetParLimits(6, 200, 500);

	cout << " Finished. " << endl << endl;
	cout << "Fitting h_dx" << identifier.Data() << endl;

	h_dx->Fit("tf_get_total_dx_fit", "R0+");

	return tf_get_total_dx_fit;
}

TF1 *get_tf_from_dx_total_fit( TF1 *tf_dx, TString plot_type = "" ){
	int start_par, end_par;
	TF1 *tf_return;

	TF1 *tf_p_gaus = new TF1("tf_p_gaus", "gaus", -2.0, 1.0);
	tf_p_gaus->SetLineColor(kRed);

	TF1 *tf_n_gaus = new TF1("tf_n_gaus", "gaus", -2.0, 1.0);
	tf_n_gaus->SetLineColor(kBlue);

	TF1 *tf_BGpol4 = new TF1("tf_BGpol4", "pol4", -2.0, 1.0);
	tf_BGpol4->SetLineColor(kViolet);

	if( plot_type == "proton" ){
		start_par = 0;
		end_par = 3;
		tf_return = tf_p_gaus;
	}
	if( plot_type == "neutron" ){
		start_par = 3;
		end_par = 6;
		tf_return = tf_n_gaus;
	}
	if( plot_type == "BG" ){
		start_par = 6;
		end_par = 11;
		tf_return = tf_BGpol4;
	}

	int return_par_counter = 0;

	for(int par = start_par; par < end_par; par++){
		tf_return->FixParameter( return_par_counter, tf_dx->GetParameter(par) );
		return_par_counter++;
	}

	return tf_return;
}

TH1D *make_dyAnticut_dx_histo( TH2D *h_dxdy, cl_SBSkine& SBSkine){
	int kine = SBSkine.kine;
//Variables for setting the dy anticut region:
	Double_t dx_p_mean = SBSkine.dx_p_mean;
	Double_t dx_p_sigma = SBSkine.dx_p_sigma;
	Double_t dx_n_mean = SBSkine.dx_n_mean;
	Double_t dx_n_sigma = SBSkine.dx_n_sigma;
	Double_t dy_mean = SBSkine.dy_mean;
	Double_t dy_sigma = SBSkine.dy_sigma;

	double dy_mult = SBSkine.dy_mult;
	double dy_mult_minus = SBSkine.dy_mult_minus;
	double dy_mult_plus = SBSkine.dy_mult_plus;

	Double_t dy_min, dy_max;
	if( SBSkine.pass == 1 ){
		dy_min = dy_mean - (dy_mult)*dy_sigma;
		dy_max = dy_mean + (dy_mult)*dy_sigma;		
	}
	else{
		dy_min = dy_mean - (dy_mult_minus)*dy_sigma;
		dy_max = dy_mean + (dy_mult_plus)*dy_sigma;		
	}

	//Histogram limits:
	TAxis *xaxis = h_dxdy->GetXaxis();
	int xbins = xaxis->GetNbins();
	Double_t xaxis_min = xaxis->GetXmin();
	Double_t xaxis_max = xaxis->GetXmax();

	TAxis *yaxis = h_dxdy->GetYaxis();
	int ybins = yaxis->GetNbins();
	Double_t yaxis_min = yaxis->GetXmin();
	Double_t yaxis_max = yaxis->GetXmax();

	//dyAnticut dxdy histogram
	TH2D *h_dxdy_dyAntiCut = new TH2D("h_dxdy_dyAntiCut", Form("SBS%i, dxdy plot with p and n peak regions set to zero", kine ), xbins, xaxis_min, xaxis_max, ybins, yaxis_min, yaxis_max);

//Let's set the region within the p and n peak regions to 0.
//We can use the values provided for dy_min and dy_max (these should correspond to some sigma parameter scaling of the dy_plot)
	for( int binx = 1; binx <= xbins; binx++ ){
		Double_t binCenter = xaxis->GetBinCenter(binx);

		//We replicate the original values outide the dy region:
		if( binCenter < dy_min || binCenter > dy_max ){
			for( int biny = 1; biny <= ybins; biny++ ){
				h_dxdy_dyAntiCut->SetBinContent( binx, biny, h_dxdy->GetBinContent( binx, biny ));
			}
		}

		//We set the values inside the dy region to zero:
		if( binCenter > dy_min && binCenter < dy_max ){
			for( int biny = 1; biny <= ybins; biny++ ){
				h_dxdy_dyAntiCut->SetBinContent( binx, biny, 0.0);
			}
		}
	}

	//We now have the dxdy histogram with the dyAnticut.
	//For dx we need the Y projection:

	TH1D *h_dx_dyAnticut = (TH1D*)h_dxdy_dyAntiCut->ProjectionY();

	TCanvas *c_dx_dyAnticut_BG = new TCanvas(Form("c_dyAnticut_BG_SBS%i", kine), Form("c_dyAnticut_BG_SBS%i", kine), 600, 500);
	h_dx_dyAnticut->Draw();

	return h_dx_dyAnticut;
}

TH2D *make_dyAnticut_dxdy_histo( TH2D *h_dxdy, cl_SBSkine& SBSkine){	
	int kine = SBSkine.kine;
//Variables for setting the dy anticut region:
	Double_t dx_p_mean = SBSkine.dx_p_mean;
	Double_t dx_p_sigma = SBSkine.dx_p_sigma;
	Double_t dx_n_mean = SBSkine.dx_n_mean;
	Double_t dx_n_sigma = SBSkine.dx_n_sigma;
	Double_t dy_mean = SBSkine.dy_mean;
	Double_t dy_sigma = SBSkine.dy_sigma;

	double dy_mult = SBSkine.dy_mult;
	double dy_mult_minus = SBSkine.dy_mult_minus;
	double dy_mult_plus = SBSkine.dy_mult_plus;

	Double_t dy_min, dy_max;
	if( SBSkine.pass == 1 ){
		dy_min = dy_mean - (dy_mult)*dy_sigma;
		dy_max = dy_mean + (dy_mult)*dy_sigma;		
	}
	else{
		dy_min = dy_mean - (dy_mult_minus)*dy_sigma;
		dy_max = dy_mean + (dy_mult_plus)*dy_sigma;		
	}

	//Histogram limits:
	TAxis *xaxis = h_dxdy->GetXaxis();
	int xbins = xaxis->GetNbins();
	Double_t xaxis_min = xaxis->GetXmin();
	Double_t xaxis_max = xaxis->GetXmax();

	TAxis *yaxis = h_dxdy->GetYaxis();
	int ybins = yaxis->GetNbins();
	Double_t yaxis_min = yaxis->GetXmin();
	Double_t yaxis_max = yaxis->GetXmax();

	//dyAnticut dxdy histogram
	TH2D *h_dxdy_dyAntiCut = new TH2D("h_dxdy_dyAntiCut", Form("SBS%i, dxdy plot with p and n peak regions set to zero", kine ), xbins, xaxis_min, xaxis_max, ybins, yaxis_min, yaxis_max);

//Let's set the region within the p and n peak regions to 0.
//We can use the values provided for dy_min and dy_max (these should correspond to some sigma parameter scaling of the dy_plot)
	for( int binx = 1; binx <= xbins; binx++ ){
		Double_t binCenter = xaxis->GetBinCenter(binx);

		//We replicate the original values outide the dy region:
		if( binCenter < dy_min || binCenter > dy_max ){
			for( int biny = 1; biny <= ybins; biny++ ){
				h_dxdy_dyAntiCut->SetBinContent( binx, biny, h_dxdy->GetBinContent( binx, biny ));
			}
		}

		//We set the values inside the dy region to zero:
		if( binCenter > dy_min && binCenter < dy_max ){
			for( int biny = 1; biny <= ybins; biny++ ){
				h_dxdy_dyAntiCut->SetBinContent( binx, biny, 0.0);
			}
		}
	}

	//We now have the dxdy histogram with the dyAnticut.

	TCanvas *c_dxdy_dyAnticut_BG = new TCanvas(Form("c_dyAnticut_BG_SBS%i", kine), Form("c_dyAnticut_BG_SBS%i", kine), 600, 500);
	h_dxdy_dyAntiCut->Draw("colz");

	return h_dxdy_dyAntiCut;
}

TF1 *fit_dyAnticut_BG( TH1D *h_dx_dyAntiCut = NULL ){
	if( h_dx_dyAntiCut == NULL ){
		cout << "ERROR IN fit_dyAnticut_BG...." << endl;
		cout << "h_dx_dyAnticut must be sent to function and CAN NOT BE NULL" << endl;
		exit(1);
	}

	TF1 *tf_dx_dyAntiCut;

	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	TAxis *xaxis = h_dx_dyAntiCut->GetXaxis();
	int xbins = xaxis->GetNbins();
	Double_t xaxis_min = xaxis->GetXmin();
	Double_t xaxis_max = xaxis->GetXmax();	

	TCanvas *c_dx_dyAntiCut_fit = new TCanvas("c_dx_dyAntiCut_fit", "c_dx_dyAntiCut_fit", 600, 500);
	h_dx_dyAntiCut->Draw();

	tf_dx_dyAntiCut = new TF1("tf_dx_dyAntiCut", background_pol4, xaxis_min, xaxis_max, 5);
	tf_dx_dyAntiCut->SetNpx( h_dx_dyAntiCut->GetNbinsX() );

	h_dx_dyAntiCut->Fit("tf_dx_dyAntiCut", "RMSE+");
	h_dx_dyAntiCut->GetXaxis()->SetRangeUser(plot_xmin, plot_xmax);

	double antiCutBG_max_x = FindMaxTF1ValueInXRange(-2.5, 0, tf_dx_dyAntiCut).x;
	double antiCutBG_max_y = FindMaxTF1ValueInXRange(-2.5, 0, tf_dx_dyAntiCut).y;

	TPaveText *tpt_dyAntiCut_fit = new TPaveText(0.60, 0.63, 0.89, 0.70, "NDC");
	tpt_dyAntiCut_fit->SetFillColor(0);
	tpt_dyAntiCut_fit->SetBorderSize(1);
	tpt_dyAntiCut_fit->AddText(Form("x-value of Anti-Cut BG max: %0.2f", antiCutBG_max_x));
	c_dx_dyAntiCut_fit->cd();
	tpt_dyAntiCut_fit->Draw("same");

	return tf_dx_dyAntiCut;

}


TF1 *scale_dyAnticut_by_TF1_eval( cl_SBSkine& SBSkine, double total_dx_x_eval = -99999.0, double x_eval_value = -2.0, TF1 *tf_dyAnticut_BG = NULL, TString identifier = "", double custom_scaling = 1.0 ){
	if( identifier != "" ){ identifier.Prepend("_"); }

	int kine = SBSkine.kine;
	double dyAnticut_prescale_factor = 0.48; //larger value subtracts more BG from data
	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	TF1 *tf_scale_dyAnticut_by_TF1_eval = new TF1(Form("tf_scale_dyAnticut_by_TF1_eval%s", identifier.Data()), "pol4", plot_xmin, plot_xmax );
	tf_scale_dyAnticut_by_TF1_eval->SetNpx( tf_dyAnticut_BG->GetNpx() );

//We need to create a background histogram using the inelastic fits
//Get the fit parameters from the function get_pN_polN_vec:
	vector< double> pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};
	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec[BG_par] = tf_dyAnticut_BG->GetParameter(BG_par);
	}

	vector< double> pN_polN_vec_scaled = {0.0, 0.0, 0.0, 0.0, 0.0};

	//Get the evaluations of the fit functions at the specifiec evaluation value:
	double dyAnticut_eval = tf_dyAnticut_BG->Eval(x_eval_value);

	double dyAnticut_BG_scale_factor = total_dx_x_eval/dyAnticut_eval;

	cout << "------- scale_dyAnticut_by_TF1_eval-------" << endl;
	cout << "total_dx_x_eval: " << total_dx_x_eval << ", dyAnticut_eval: " << dyAnticut_eval << endl;
	cout << "dyAnticut_BG_scale_factor: " << dyAnticut_BG_scale_factor << endl;
	cout << "------------------------------------------" << endl;

//Now we actually scale it up:
	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec_scaled[BG_par] = custom_scaling*dyAnticut_prescale_factor*pN_polN_vec[BG_par]*dyAnticut_BG_scale_factor;

		tf_scale_dyAnticut_by_TF1_eval->FixParameter(BG_par, pN_polN_vec_scaled[BG_par]);
	}
	return tf_scale_dyAnticut_by_TF1_eval;
}

TF1 *scale_dyAnticut_by_TF1_eval_MaxMin( cl_SBSkine& SBSkine, double total_dx_x_eval = -99999.0, double low_x_eval_value = -2.0, double high_x_eval_value = 1.0, TF1 *tf_dyAnticut_BG = NULL, TString identifier = "", bool systematics = false, bool bg_sub_systematics_analysis = false, double bg_sub_syst_val = 1.0, double custom_scaling = 1.0 ){
	if( identifier != "" ){ identifier.Prepend("_"); }

	int kine = SBSkine.kine;

	double dyAnticut_prescale_factor;

	if( systematics && !bg_sub_systematics_analysis ){
		if( kine == 4 ){
			dyAnticut_prescale_factor = 0.40; //0.35867812; //Larger number subtracts more background from data.
		}
		if( kine == 8 ){
			dyAnticut_prescale_factor = 0.32; //0.32276110; //0.43; //0.32276110; //Larger number subtracts more background from data.
		}
		if( kine == 9 ){
			dyAnticut_prescale_factor = 0.40; //0.260803; //0.48; //0.218666; //0.11132874; //Larger number subtracts more background from data.
		}		
	}

	else if( systematics && bg_sub_systematics_analysis ){
		dyAnticut_prescale_factor = bg_sub_syst_val;
	}
	else{
		if( kine == 4 ){
			if( SBSkine.sbsfieldscale == 30 ){
				dyAnticut_prescale_factor = 0.6; //0.35867812; //Larger number subtracts more background from data.
			}
			if( SBSkine.sbsfieldscale == 50 ){
				dyAnticut_prescale_factor = 0.0; //0.35867812; //Larger number subtracts more background from data.
			}
			
		}
		if( kine == 8 ){
			dyAnticut_prescale_factor = 0.42; //0.322276110; //0.43; //0.32276110; //Larger number subtracts more background from data.
		}
		if( kine == 9 ){
			dyAnticut_prescale_factor = 0.60; //0.260803; //0.48; //0.218666; //0.11132874; //Larger number subtracts more background from data.
		}		
	}

	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	TF1 *tf_scale_dyAnticut_by_TF1_eval = new TF1(Form("tf_scale_dyAnticut_by_TF1_eval%s", identifier.Data()), "pol4", plot_xmin, plot_xmax );
	tf_scale_dyAnticut_by_TF1_eval->SetNpx( tf_dyAnticut_BG->GetNpx() );

//We need to create a background histogram using the inelastic fits
//Get the fit parameters from the function get_pN_polN_vec:
	vector< double> pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};
	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec[BG_par] = tf_dyAnticut_BG->GetParameter(BG_par);
	}

	vector< double> pN_polN_vec_scaled = {0.0, 0.0, 0.0, 0.0, 0.0};

	//Get the evaluations of the fit functions at the specifiec evaluation value:
	double low_x_dyAnticut_eval = tf_dyAnticut_BG->Eval(low_x_eval_value);
	double high_x_dyAnticut_eval = tf_dyAnticut_BG->Eval(high_x_eval_value);

	double dyAnticut_eval;

	//pick the lower evaluation value:
	if( low_x_dyAnticut_eval > high_x_dyAnticut_eval ){
		dyAnticut_eval = high_x_dyAnticut_eval;
	}
	if( low_x_dyAnticut_eval < high_x_dyAnticut_eval ){
		dyAnticut_eval = low_x_dyAnticut_eval;
	}

	double dyAnticut_BG_scale_factor = total_dx_x_eval/dyAnticut_eval;

	cout << "------- scale_dyAnticut_by_TF1_eval-------" << endl;
	cout << "total_dx_x_eval: " << total_dx_x_eval << ", dyAnticut_eval: " << dyAnticut_eval << endl;
	cout << "dyAnticut_BG_scale_factor: " << dyAnticut_BG_scale_factor << endl;
	cout << "------------------------------------------" << endl;

//Now we actually scale it up:
	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec_scaled[BG_par] = custom_scaling*dyAnticut_prescale_factor*pN_polN_vec[BG_par]*dyAnticut_BG_scale_factor;

		tf_scale_dyAnticut_by_TF1_eval->FixParameter(BG_par, pN_polN_vec_scaled[BG_par]);
	}
	return tf_scale_dyAnticut_by_TF1_eval;
}

//This function is used to scale the dyAnticut BG function to fit the original fit from the total dx BG
//Therefore, we need to provide this function with the original fit from the total dx.
TF1 *scale_dyAnticut_by_ScalarCoeff( cl_SBSkine& SBSkine, double total_dx_fit_Scalar_Coefficient = -99999.0, TF1 *tf_dyAnticut_BG = NULL ){

	int kine = SBSkine.kine;
	double dyAnticut_prescale_factor;
	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	double parTotalDyAnticutDx[11];

//We need to create a background histogram using the inelastic fits
//Get the fit parameters from the function get_pN_polN_vec:
	vector< double> pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};

	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec[BG_par] = tf_dyAnticut_BG->GetParameter(BG_par);
	}

	vector< double> pN_polN_vec_scaled = {0.0, 0.0, 0.0, 0.0, 0.0};

//Get the ratio of the Scalar Coefficients:

	double dyAnticut_BG_scale_factor = total_dx_fit_Scalar_Coefficient/pN_polN_vec[0];

	cout << "--------------scale_dyAnticut_by_ScalarCoeff--------------" << endl;
	cout << "Data Scalar Coefficient: " << total_dx_fit_Scalar_Coefficient << endl;
	cout << "dyAnticut Scalar: " << pN_polN_vec[0] << endl;
	cout << "dyAnticut BG Scale Factor: " << dyAnticut_BG_scale_factor << endl;
	cout << "----------------------------" << endl;

//Create a new TF1 function with scaled up parameters:

	TF1 *tf_dyAnticut_BG_scaled = new TF1("tf_dyAnticut_BG_scaled", "pol4", plot_xmin, plot_xmax);
	tf_dyAnticut_BG_scaled->SetNpx(tf_dyAnticut_BG->GetNpx());

	for( size_t BGsub_par = 0; BGsub_par < pN_polN_vec.size(); BGsub_par++ ){

		pN_polN_vec_scaled[BGsub_par] = dyAnticut_BG_scale_factor*pN_polN_vec[BGsub_par];

		cout << "Par " << BGsub_par << "- Pre-scale: " << pN_polN_vec[BGsub_par] << ", Post-scale: " << pN_polN_vec_scaled[BGsub_par] << endl;
		
		tf_dyAnticut_BG_scaled->FixParameter(BGsub_par, pN_polN_vec_scaled[BGsub_par]);

	}

	return tf_dyAnticut_BG_scaled;
}


TH1D *make_inelastic_dx_histo( int kine, bool anticut_reject = false, bool inelastic_BG_wcut = false ){
	//We need to create a background histogram using the inelastic fits
	//Get the fit parameters from the function get_pN_polN_vec:
	// get_pN_polN_vec(kine = 8,anticut_reject, inelastic_BG_wcut)
	vector< double> pN_polN_vec;

	if( kine == 9 ){
		pN_polN_vec = get_pN_polN_vec( kine, true, false );		
	}
	else{
		pN_polN_vec = get_pN_polN_vec( kine, false, false );		
	}

	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	TF1 *tf_make_dx_inelastics_histo = new TF1("tf_make_dx_inelastics_histo", "pol4", plot_xmin, plot_xmax );
	tf_make_dx_inelastics_histo->SetNpx(500);

	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		tf_make_dx_inelastics_histo->FixParameter( BG_par, pN_polN_vec[BG_par] );
	}

	TH1D *h_make_dx_inelastics_histo = (TH1D*)tf_make_dx_inelastics_histo->GetHistogram();

	return h_make_dx_inelastics_histo;
}

TH1D *make_inelastic_dx_histo_from_BGfile( cl_SBSkine& SBSkine, bool anticut_reject = false, bool inelastic_BG_wcut = false ){
	
	int kine = SBSkine.kine;

	//We need to create a background histogram using the inelastic fits
	//Get the fit parameters from the function get_pN_polN_vec:
	// get_pN_polN_vec(kine = 8,anticut_reject, inelastic_BG_wcut)
	vector< double> pN_polN_vec;
	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	if( kine == 9 ){
		pN_polN_vec = get_pN_polN_vec_from_TF1_fit_to_file( SBSkine, false, false, plot_xmin, plot_xmax );		
	}
	else{
		pN_polN_vec = get_pN_polN_vec_from_TF1_fit_to_file( SBSkine, false, false, plot_xmin, plot_xmax );		
	}

	TF1 *tf_make_dx_inelastics_histo = new TF1("tf_make_dx_inelastics_histo", "pol4", plot_xmin, plot_xmax );
	tf_make_dx_inelastics_histo->SetNpx(500);

	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		tf_make_dx_inelastics_histo->FixParameter( BG_par, pN_polN_vec[BG_par] );
	}

	TH1D *h_make_dx_inelastics_histo = (TH1D*)tf_make_dx_inelastics_histo->GetHistogram();

	return h_make_dx_inelastics_histo;
}

TF1 *fit_inelastics_BG( TH1D *h_dx_inelastics = NULL ){
	if( h_dx_inelastics == NULL ){
		cout << "ERROR IN fit_inelastics...." << endl;
		cout << "h_dx_inelastics must be sent to function and CAN NOT BE NULL" << endl;
		exit(1);
	}

	TF1 *tf_dx_inelastics_BG;

	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	TAxis *xaxis = h_dx_inelastics->GetXaxis();
	int xbins = xaxis->GetNbins();
	Double_t xaxis_min = xaxis->GetXmin();
	Double_t xaxis_max = xaxis->GetXmax();	

	if( false ){
		TCanvas *c_dx_inelastics_fit = new TCanvas("c_dx_inelastics_fit", "c_dx_inelastics_fit", 600, 500);
		h_dx_inelastics->Draw();		
	}

	TF1 *tf_dx_inelastics = new TF1("tf_dx_inelastics", background_pol4, xaxis_min, xaxis_max, 5);
	tf_dx_inelastics->SetNpx( h_dx_inelastics->GetNbinsX() );

	h_dx_inelastics->Fit("tf_dx_inelastics", "R0+");
	h_dx_inelastics->GetXaxis()->SetRangeUser(plot_xmin, plot_xmax);

	double antiCutBG_max_x = FindMaxTF1ValueInXRange(-2.5, 0, tf_dx_inelastics).x;
	double antiCutBG_max_y = FindMaxTF1ValueInXRange(-2.5, 0, tf_dx_inelastics).y;

	TPaveText *tpt_inelastics_fit = new TPaveText(0.60, 0.63, 0.89, 0.70, "NDC");
	tpt_inelastics_fit->SetFillColor(0);
	tpt_inelastics_fit->SetBorderSize(1);
	tpt_inelastics_fit->AddText(Form("x-value of Anti-Cut BG max: %0.2f", antiCutBG_max_x));
	tpt_inelastics_fit->Draw("same");

	return tf_dx_inelastics;

}


TF1 *scale_inelastics_by_TF1_eval( cl_SBSkine& SBSkine, double total_dx_x_eval = -99999.0, double x_eval_value = -2.0, TF1 *tf_inelastics_BG = NULL, TString identifier = "", double custom_scaling = 1.0 ){
	if( identifier != "" ){ identifier.Prepend("_"); }

	int kine = SBSkine.kine;
	double inelastics_prescale_factor = 0.40; //higher number should remove more BG
	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	TF1 *tf_scale_inelastics_by_TF1_eval = new TF1(Form("tf_scale_inelastics_by_TF1_eval%s", identifier.Data()), "pol4", plot_xmin, plot_xmax );
	tf_scale_inelastics_by_TF1_eval->SetNpx( tf_inelastics_BG->GetNpx() );

//We need to create a background histogram using the inelastic fits
//Get the fit parameters from the function get_pN_polN_vec:
	vector< double> pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};
	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec[BG_par] = tf_inelastics_BG->GetParameter(BG_par);
	}

	vector< double> pN_polN_vec_scaled = {0.0, 0.0, 0.0, 0.0, 0.0};

	//Get the evaluations of the fit functions at the specifiec evaluation value:
	double inelastics_eval = tf_inelastics_BG->Eval(x_eval_value);

	double inelastics_BG_scale_factor = total_dx_x_eval/inelastics_eval;

	cout << "------- scale_inelastics_by_TF1_eval-------" << endl;
	cout << "total_dx_x_eval: " << total_dx_x_eval << ", inelastics_eval: " << inelastics_eval << endl;
	cout << "inelastics_BG_scale_factor: " << inelastics_BG_scale_factor << endl;
	cout << "------------------------------------------" << endl;

//Now we actually scale it up:
	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec_scaled[BG_par] = custom_scaling*inelastics_prescale_factor*pN_polN_vec[BG_par]*inelastics_BG_scale_factor;

		tf_scale_inelastics_by_TF1_eval->FixParameter(BG_par, pN_polN_vec_scaled[BG_par]);
	}
	return tf_scale_inelastics_by_TF1_eval;
}

TF1 *scale_inelastics_by_TF1_eval_MaxMin( cl_SBSkine& SBSkine, double total_dx_x_eval = -99999.0, double low_x_eval_value = -2.0, double high_x_eval_value = 1.0, TF1 *tf_inelastics_BG = NULL, TString identifier = "", bool bg_sub_systematics_analysis = false, double bg_sub_syst_val = 1.0, double custom_scaling = 1.0 ){
	
	if( identifier != "" ){ identifier.Prepend("_"); }

	int kine = SBSkine.kine;
	double inelastics_prescale_factor;

	if( bg_sub_systematics_analysis ){
		inelastics_prescale_factor = bg_sub_syst_val;
	}
	if( kine == 4 ){
		inelastics_prescale_factor = 0.50437188; //Larger number subtracts more background from data.
	}
	if( kine == 8 ){
		inelastics_prescale_factor = 0.50437188; //Larger number subtracts more background from data.
	}
	if( kine == 9 ){
		inelastics_prescale_factor = 0.504; //0.315359; //0.308942; //0.347967; //0.31546166; //Larger number subtracts more background from data.
	}

	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	cout << endl;
	cout << "Scaling inelastics by TF1 Eval max/min: " << endl;
	cout << "inelastics prescale_factor: " << inelastics_prescale_factor << endl;
	cout << endl;

	TF1 *tf_scale_inelastics_by_TF1_eval = new TF1(Form("tf_scale_inelastics_by_TF1_eval%s", identifier.Data()), "pol4", plot_xmin, plot_xmax );
	tf_scale_inelastics_by_TF1_eval->SetNpx( tf_inelastics_BG->GetNpx() );

//We need to create a background histogram using the inelastic fits
//Get the fit parameters from the function get_pN_polN_vec:
	vector< double> pN_polN_vec = {0.0, 0.0, 0.0, 0.0, 0.0};
	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec[BG_par] = tf_inelastics_BG->GetParameter(BG_par);
	}

	vector< double> pN_polN_vec_scaled = {0.0, 0.0, 0.0, 0.0, 0.0};

	//Get the evaluations of the fit functions at the specifiec evaluation value:
	double low_x_inelastics_eval = tf_inelastics_BG->Eval(low_x_eval_value);
	double high_x_inelastics_eval = tf_inelastics_BG->Eval(high_x_eval_value);

	double inelastics_eval;

	//pick the lower evaluation value:
	if( low_x_inelastics_eval > high_x_inelastics_eval ){
		inelastics_eval = high_x_inelastics_eval;
	}
	if( low_x_inelastics_eval < high_x_inelastics_eval ){
		inelastics_eval = low_x_inelastics_eval;
	}

	double inelastics_BG_scale_factor = total_dx_x_eval/inelastics_eval;

	cout << "------- scale_inelastics_by_TF1_eval-------" << endl;
	cout << "total_dx_x_eval: " << total_dx_x_eval << ", inelastics_eval: " << inelastics_eval << endl;
	cout << "inelastics_BG_scale_factor: " << inelastics_BG_scale_factor << endl;
	cout << "------------------------------------------" << endl;

//Now we actually scale it up:
	for( int BG_par = 0; BG_par < 5; BG_par++ ){
		pN_polN_vec_scaled[BG_par] = custom_scaling*inelastics_prescale_factor*pN_polN_vec[BG_par]*inelastics_BG_scale_factor;

		tf_scale_inelastics_by_TF1_eval->FixParameter(BG_par, pN_polN_vec_scaled[BG_par]);
	}
	return tf_scale_inelastics_by_TF1_eval;
}


//This function is used to scale the inelastic BG function to fit the original fit from the total dx BG
//Therefore, we need to provide this function with the original fit from the total dx.
TF1 *scale_inelastics_by_ScalarCoeff( cl_SBSkine& SBSkine, double total_dx_fit_Scalar_Coefficient = -99999.0, TF1 *tf_inelastics_BG = NULL ){

	int kine = SBSkine.kine;
	double plot_xmin = -2.0;
	double plot_xmax = 1.0;

	double parTotalMCinelasticsDx[11];


//We need to create a background histogram using the inelastic fits
//Get the fit parameters from the function get_pN_polN_vec:
	vector< double> pN_polN_vec = get_pN_polN_vec( kine, false, false );
	vector< double> pN_polN_vec_scaled = {0.0, 0.0, 0.0, 0.0, 0.0};

//Get the ratio of the Scalar Coefficients:

	double MCinelastics_BG_scale_factor = total_dx_fit_Scalar_Coefficient/pN_polN_vec[0];

	cout << "----------------------------" << endl;
	cout << "Data Scalar Coefficient: " << total_dx_fit_Scalar_Coefficient << endl;
	cout << "MC inelastics Scalar: " << pN_polN_vec[0] << endl;
	cout << "MC inelastics BG Scale Factor: " << MCinelastics_BG_scale_factor << endl;
	cout << "----------------------------" << endl;

//Create a new TF1 function with scaled up parameters:

	TF1 *tf_inelastics_BG_scaled = new TF1("tf_inelastics_BG_scaled", "pol4", plot_xmin, plot_xmax);
	tf_inelastics_BG_scaled->SetNpx( tf_inelastics_BG->GetNpx() );

	for( size_t BGsub_par = 0; BGsub_par < pN_polN_vec.size(); BGsub_par++ ){

		pN_polN_vec_scaled[BGsub_par] = MCinelastics_BG_scale_factor*pN_polN_vec[BGsub_par];

		cout << "Par " << BGsub_par << "- Pre-scale: " << pN_polN_vec[BGsub_par] << ", Post-scale: " << pN_polN_vec_scaled[BGsub_par] << endl;
		
		tf_inelastics_BG_scaled->FixParameter(BGsub_par, pN_polN_vec_scaled[BGsub_par]);

	}

	return tf_inelastics_BG_scaled;
}

TH1D *make_simc_dx_histo( TH1D *h_simc_dx_n, TH1D *h_simc_dx_p, TString identifier = "", double scale_n = 1.0, double scale_p = 1.0){

	if( identifier != "" ){
		identifier.Prepend("_");
	}

	h_simc_dx_n->Scale(scale_n);
	h_simc_dx_p->Scale(scale_p);

	TH1D *h_simc_dx_histo = (TH1D*)h_simc_dx_p->Clone(Form("h_simc_dx_histo%s", identifier.Data()) );

	h_simc_dx_histo->Add( h_simc_dx_n, 1 );

	return h_simc_dx_histo;
}

double get_data_dx_n_maxVal( TH1D *h_data, double n_mean, double n_sigma ){
	double original_range_min = h_data->GetXaxis()->GetXmin();
	double original_range_max = h_data->GetXaxis()->GetXmax();

	double new_range_min = n_mean - n_sigma;
	double new_range_max = n_mean + n_sigma;

	double n_maxVal = 0.0;

	h_data->GetXaxis()->SetRangeUser( new_range_min, new_range_max );

	n_maxVal = h_data->GetMaximum();

	h_data->GetXaxis()->SetRangeUser( original_range_min, original_range_max );

	return n_maxVal;
}

#endif