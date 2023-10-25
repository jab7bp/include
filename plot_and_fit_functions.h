#ifndef PLOT_AND_FIT_FUNCS_H
#define PLOT_AND_FIT_FUNCS_H

TLegend *tl_dx_mag0, *tl_dx_mag, *tl_dy_mag, *tl_dx_mag_dyBG;
TF1 *tf_data_dx_p, *tf_data_dx_n;
TF1 *tf_dx_bg, *fit_dx_dyBG, *fit_polBG_with_reject, *fit_dx, *fit_dx_origBG;
TF1 *tf_BG_dx_pol4, *tf_BG_dx_pol5, *tf_BG_dx_pol6, *tf_BG_dx_pol8;
TF1 *tf_BG_dx_wcut_pol4, *tf_BG_dx_wcut_pol5, *tf_BG_dx_wcut_pol6, *tf_BG_dx_wcut_pol8;
Double_t par_BG_reject[5], par_dy_bg_reject[5], parPoly[9], parPol8[9], parPol2[3], parPol3[4], parPol4[5], parPol5[6], parPol6[7], parPol7[8], parElastics[6];
Double_t dx_total_fit_par[10];

Double_t background_pol2(Double_t *x, Double_t *parPol2){
	return parPol2[0] + parPol2[1]*x[0] + parPol2[2]*x[0]*x[0];
}
Double_t background_pol3(Double_t *x, Double_t *parPol3){
	return parPol3[0] + parPol3[1]*x[0] + parPol3[2]*x[0]*x[0] + parPol3[3]*x[0]*x[0]*x[0];
}
Double_t background_pol4(Double_t *x, Double_t *parPol4){
	return parPol4[0] + parPol4[1]*x[0] + parPol4[2]*x[0]*x[0] + parPol4[3]*x[0]*x[0]*x[0] + parPol4[4]*x[0]*x[0]*x[0]*x[0];
}
Double_t background_pol5(Double_t *x, Double_t *parPol5){
	return parPol5[0] + parPol5[1]*x[0] + parPol5[2]*x[0]*x[0] + parPol5[3]*x[0]*x[0]*x[0] + parPol5[4]*x[0]*x[0]*x[0]*x[0] + parPol5[5]*x[0]*x[0]*x[0]*x[0]*x[0];
}
Double_t background_pol6(Double_t *x, Double_t *parPol6){
	return parPol6[0] + parPol6[1]*x[0] + parPol6[2]*x[0]*x[0] + parPol6[3]*x[0]*x[0]*x[0] + parPol6[4]*x[0]*x[0]*x[0]*x[0] + parPol6[5]*x[0]*x[0]*x[0]*x[0]*x[0] + parPol6[6]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
}
Double_t background_pol8(Double_t *x, Double_t *parPol8){
	return parPol8[0] + parPol8[1]*x[0] + parPol8[2]*pow(x[0], 2) + parPol8[3]*pow(x[0], 3) + parPol8[4]*pow(x[0], 4) + parPol8[5]*pow(x[0], 5) + parPol8[6]*pow(x[0], 6) + parPol8[7]*pow(x[0], 7) + parPol8[8]*pow(x[0], 8);
}

Double_t pol8_bin_to_y(Int_t bin, Double_t *parPol8){
	return parPol8[0] + parPol8[1]*bin + parPol8[2]*pow(bin, 2) + parPol8[3]*pow(bin, 3) + parPol8[4]*pow(bin, 4) + parPol8[5]*pow(bin, 5) + parPol8[6]*pow(bin, 6) + parPol8[7]*pow(bin, 7) + parPol8[8]*pow(bin, 8);
}

Double_t yieldPeaks(Double_t *x, Double_t *par){
	return par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2)) + par[3]*exp((-0.5)*pow(((x[0] -  par[4])/par[5]),2));
}

Double_t fullFitElastics(Double_t *x, Double_t *parElastics){
	return yieldPeaks(x, parElastics);
}

Double_t polBG_with_reject( Double_t *x, Double_t *par ){
	// cout << "Fitting polynomial BG while rejecting: -" << dy_pn_sigma << " thru " << dy_pn_sigma << endl;
	if( x[0] > -1.5 && x[0] < -0.5 ){
		TF1::RejectPoint();
		return 0;
	}
	if( x[0] > -0.25 && x[0] < 0.5 ){
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t fit_gaus(Double_t *x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

Double_t fit_full_dx_peaks_with_pol2BG( Double_t *x, Double_t *dx_total_fit_par ){
	Double_t p_peak = fit_gaus( x, dx_total_fit_par );
	Double_t n_peak = fit_gaus( x, &dx_total_fit_par[3] );
	Double_t dx_BG = background_pol2(x, &dx_total_fit_par[6] );

	return p_peak + n_peak + dx_BG;
}

Double_t fit_full_dx_peaks_with_pol3BG( Double_t *x, Double_t *dx_total_fit_par ){
	Double_t p_peak = fit_gaus( x, dx_total_fit_par );
	Double_t n_peak = fit_gaus( x, &dx_total_fit_par[3] );
	Double_t dx_BG = background_pol3(x, &dx_total_fit_par[6] );

	return p_peak + n_peak + dx_BG;
}

Double_t fit_full_dx_peaks_with_pol4BG( Double_t *x, Double_t *dx_total_fit_par ){
	Double_t p_peak = fit_gaus( x, dx_total_fit_par );
	Double_t n_peak = fit_gaus( x, &dx_total_fit_par[3] );
	Double_t dx_BG = background_pol4(x, &dx_total_fit_par[6] );

	return p_peak + n_peak + dx_BG;
}

Double_t fit_dy(Double_t *x, Double_t *par){
	double total_fit = 0.0, dy_gaus = 0.0, BG_dy_gaus = 0.0;

	dy_gaus = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	BG_dy_gaus = par[3]*exp((-0.5)*pow(((x[0] -  par[4])/par[5]),2));

	total_fit = dy_gaus + BG_dy_gaus;

	return total_fit;
}
Double_t dyBackground(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t dyYieldPeaks(Double_t *x, Double_t *par){
	return par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
}


Double_t fulldyFunction(Double_t *x, Double_t *par){
	return dyYieldPeaks(x,par) + dyBackground(x, &par[3]);
}

// Double_t fit_dy_polyBG(Double_t *x, Double_t *par){
// 	double total_fit = 0.0, dy_gaus = 0.0, polyBG = 0.0;

// 	dy_gaus = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
// 	polyBG = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2)) + par[3]*exp((-0.5)*pow(((x[0] -  par[4])/par[5]),2));

// 	total_fit = dy_gaus + polyBG;

// 	return total_fit;
// }

Double_t fit_dxdy_mag0(Double_t * x, Double_t *par){
	double total_fit = 0.0, p_gaus = 0.0, n_gaus = 0.0, BG_gaus = 0.0, pn_gaus = 0.0;

		pn_gaus = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
		BG_gaus = par[3]*exp((-0.5)*pow(((x[0] -  par[4])/par[5]),2));

		total_fit = pn_gaus + BG_gaus;

	return total_fit;	
}

Double_t fit_dxdy(Double_t * x, Double_t *par){
	double total_fit = 0.0, p_gaus = 0.0, n_gaus = 0.0, BG_gaus = 0.0, pn_gaus = 0.0;
	
	p_gaus = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	n_gaus = par[3]*exp((-0.5)*pow(((x[0] -  par[4])/par[5]),2));
	BG_gaus = par[6]*exp((-0.5)*pow(((x[0] -  par[7])/par[8]),2));

	total_fit = p_gaus + n_gaus + BG_gaus;


	return total_fit;
}

Double_t plot_dxdy_polyBG_fit(Double_t *par, bool plot_sim, bool dyBG, TH1D *h_dx_final_histo){
	Double_t fit_x[2400], gaus1_y[2400], gaus2_y[2400], polyBG_y[2400];
	Int_t fit_n = 550;
	int fit_cnt = 0;

	for(double i = -4; i < 1.5; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		gaus2_y[fit_cnt] = par[3]*exp((-0.5)*pow(((i -  par[4])/par[5]),2));
		// polyBG_y[fit_cnt] = par[10]*pow(i, 4) + par[9]*pow(i, 3) + par[8]*pow(i, 2) + par[7]*i + par[6];
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	gr_gaus1->SetLineColor(1);
	gr_gaus1->SetLineStyle(7);
	TGraph *gr_gaus2 = new TGraph(fit_n, fit_x, gaus2_y);
	gr_gaus2->SetLineColor(6);
	gr_gaus2->SetLineStyle(6);

	gr_gaus1->Draw("same");
	gr_gaus2->Draw("same");

	tl_dx_mag = new TLegend(0.50, 0.65, 0.75, 0.88,"",  "NDC");
	tl_dx_mag->AddEntry(h_dx_final_histo, "Data");
	tl_dx_mag->AddEntry(fit_dx, "Total data fit");
	tl_dx_mag->AddEntry(gr_gaus1, "p Fit");
	tl_dx_mag->AddEntry(gr_gaus2, "n Fit");

	tl_dx_mag->Draw("same");
	
	return 1;

}

Double_t plot_dy_polyBG_fit(Double_t *par){
	Double_t fit_x[2400], gaus1_y[2400], polyBG_y[2400];
	Int_t fit_n = 250;
	int fit_cnt = 0;

	for(double i = -1.25; i < 1.25; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		// polyBG_y[fit_cnt] = par[3]*pow(i, 4) + par[4]*pow(i, 3) + par[5]*pow(i, 2) + par[6]*i + par[7];
		polyBG_y[fit_cnt] = par[3] + par[4]*i + par[5]*i*i + par[6]*i*i*i + par[7]*i*i*i*i;
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	gr_gaus1->SetLineColor(3);
	TGraph *gr_polyBG_y = new TGraph(fit_n, fit_x, polyBG_y);

	gr_polyBG_y->SetLineColor(7);
	gr_gaus1->Draw("same");
	gr_polyBG_y->Draw("same");

	return 1;

}

Double_t plot_dxdy_fit_mag0(Double_t *par){
	Double_t fit_x[1800], gaus1_y[1800], gaus2_y[1800], gaus3_y[1800];
	Int_t fit_n = 500;
	int fit_cnt = 0;

	for(double i = -3; i < 3; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		gaus2_y[fit_cnt] = par[3]*exp((-0.5)*pow(((i -  par[4])/par[5]),2));
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	gr_gaus1->SetLineColor(3);
	TGraph *gr_gaus2 = new TGraph(fit_n, fit_x, gaus2_y);
	gr_gaus2->SetLineColor(7);

	gr_gaus1->Draw("same");
	gr_gaus2->Draw("same");

	tl_dx_mag0 = new TLegend(0.15, 0.75, 0.35, 0.85, "", "NDC");
	tl_dx_mag0->AddEntry(gr_gaus1, "pn Fit");
	tl_dx_mag0->AddEntry(gr_gaus2, "BG fit");
	tl_dx_mag0->Draw("same");
	
	return 1;

}

Double_t plot_dxdy_fit(Double_t *par){
	Double_t fit_x[1800], gaus1_y[1800], gaus2_y[1800], gaus3_y[1800];
	Int_t fit_n = 500;
	int fit_cnt = 0;

	for(double i = -3; i < 3; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		gaus2_y[fit_cnt] = par[3]*exp((-0.5)*pow(((i -  par[4])/par[5]),2));
		gaus3_y[fit_cnt] = par[6]*exp((-0.5)*pow(((i -  par[7])/par[8]),2));
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	gr_gaus1->SetLineColor(3);
	TGraph *gr_gaus2 = new TGraph(fit_n, fit_x, gaus2_y);
	gr_gaus2->SetLineColor(6);
	TGraph *gr_gaus3 = new TGraph(fit_n, fit_x, gaus3_y);
	gr_gaus3->SetLineColor(7);

	gr_gaus1->Draw("same");
	gr_gaus2->Draw("same");
	gr_gaus3->Draw("same");

	tl_dx_mag_dyBG = new TLegend(0.15, 0.75, 0.35, 0.85,"",  "NDC");
	tl_dx_mag_dyBG->AddEntry(gr_gaus1, "p Fit");
	tl_dx_mag_dyBG->AddEntry(gr_gaus2, "n Fit");
	tl_dx_mag_dyBG->AddEntry(gr_gaus3, "dyBG fit");
	tl_dx_mag_dyBG->Draw("same");
	
	return 1;

}
Double_t plot_dy_fit(Double_t *par){
	Double_t fit_x[1800], gaus1_y[1800], gaus2_y[1800];
	Int_t fit_n = 500;
	int fit_cnt = 0;

	for(double i = -3; i < 3; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		gaus2_y[fit_cnt] = par[3]*exp((-0.5)*pow(((i -  par[4])/par[5]),2));
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	gr_gaus1->SetLineColor(3);
	TGraph *gr_gaus2 = new TGraph(fit_n, fit_x, gaus2_y);
	gr_gaus2->SetLineColor(7);

	gr_gaus1->Draw("same");
	gr_gaus2->Draw("same");
	
	return 1;
}

Double_t plot_pol4(Double_t *par){
	Double_t fit_x[2400], polyBG_y[2400];
	Int_t fit_n = 450;
	int fit_cnt = 0;

	for(double i = -2.5; i < 2.0; i+=.01){
		fit_x[fit_cnt] = i;
		polyBG_y[fit_cnt] = par[0] + par[1]*i + par[2]*i*i + par[3]*i*i*i + par[4]*i*i*i*i;
		fit_cnt++;
	}

	TGraph *gr_polyBG_y = new TGraph(fit_n, fit_x, polyBG_y);

	gr_polyBG_y->SetLineColor(1);
	gr_polyBG_y->Draw("same");

	return 1;

}


Int_t cnt_ellipse_steps(double h, double k, double a, double b){
	Int_t cnt = 0;
	double y = 0.0;
	for(double x = -abs(h-a); x <= abs(h+a); x+=0.001){
		y = k + (b)*sqrt( 1 - ((pow( (x - h), 2))/(a*a)) );
		cnt++;
	}
	// 0.001
	for(double x = abs(h+a); x >= -abs(h-a); x-=0.001){
		y = k - (b)*sqrt( 1 - ( (pow( (x - h), 2))/(a*a)) );
		cnt++;
	}
	return cnt;
}

Int_t set_tcg_ellipse(double h, double k, double a, double b, TCutG *tcg){
	Int_t cnt = 0;
	double y = 0.0;
	cout << "-------------------start----------------------" << endl;
	for(double x = -abs(h-a); x <= abs(h+a); x+=0.001){
		y = k + (b)*sqrt( 1 - ((pow( (x - h), 2))/(a*a)) );
		tcg->SetPoint(cnt, x, y);
		cnt++;
		// cout << "x: " << x << ", y: " << y << endl;
	}
	
	for(double x = abs(h+a); x >= -abs(h-a); x-=0.001){
		y = k - (b)*sqrt( 1 - ( (pow( (x - h), 2))/(a*a)) );
		tcg->SetPoint(cnt, x, y);
		// cout << "x: " << x << ", y: " << y << endl;
		cnt++;
	}
	cout << "-------------------end----------------------" << endl;
	return cnt;
}

void dyAnticut_and_dxBG_fit( kine_ff_extract& SBSkine ){

	double BG_prescale_factor = 0.48;

	cout << "Using dxdy dyAntiCut and finding BG fit " << endl;

	TH1D *h_dx_dyAntiCut = SBSkine.h_data_dx_dyAntiCut;
	int BG_polN = SBSkine.dx_dyAnticut_BG_polN;

	cout << "BG_polN: " << BG_polN;

	TAxis *xaxis = h_dx_dyAntiCut->GetXaxis();
	int xbins = xaxis->GetNbins();
	Double_t xaxis_min = xaxis->GetXmin();
	Double_t xaxis_max = xaxis->GetXmax();

	TF1 *tf_dx_dyAntiCut, *tf_data_dx_dyAntiCut_BG_scaled;


	if( BG_polN == 2){
		tf_dx_dyAntiCut = new TF1("tf_dx_dyAntiCut", background_pol2, xaxis_min, xaxis_max, 3);
		tf_data_dx_dyAntiCut_BG_scaled = new TF1("tf_data_dx_dyAntiCut_BG_scaled", background_pol2, xaxis_min, xaxis_max, 3);
	}
	if( BG_polN == 3){
		tf_dx_dyAntiCut = new TF1("tf_dx_dyAntiCut", background_pol3, xaxis_min, xaxis_max, 4);
		tf_data_dx_dyAntiCut_BG_scaled = new TF1("tf_data_dx_dyAntiCut_BG_scaled", background_pol3, xaxis_min, xaxis_max, 4);
	}
	if( BG_polN == 4){
		tf_dx_dyAntiCut = new TF1("tf_dx_dyAntiCut", background_pol4, xaxis_min, xaxis_max, 5);
		tf_data_dx_dyAntiCut_BG_scaled = new TF1("tf_data_dx_dyAntiCut_BG_scaled", background_pol4, xaxis_min, xaxis_max, 5);
	}
	tf_dx_dyAntiCut->SetNpx(500);
	tf_data_dx_dyAntiCut_BG_scaled->SetNpx(500);

	h_dx_dyAntiCut->Fit("tf_dx_dyAntiCut", "R+");

	//Now we need to scale up the BG fit to the scale of the original BG from the total fit:
	//We get our scale factor by looking at the pure scalar p0 coefficient:
	// double antiCut_BGscale_factor = SBSkine.tf_data_dx_total->GetParameter(6)/tf_dx_dyAntiCut->GetParameter(0);

	double antiCut_BGscale_factor = BG_prescale_factor*SBSkine.h_data->GetBinContent(SBSkine.h_data->FindBin(-1.75))/tf_dx_dyAntiCut->Eval(-1.75);

	cout << "-----------------------------------------------" << endl;
	cout << "Scaling dyAntiCut by: " << antiCut_BGscale_factor << endl;
	cout << "-----------------------------------------------" << endl;

	tf_data_dx_dyAntiCut_BG_scaled->FixParameter( 0, antiCut_BGscale_factor*tf_dx_dyAntiCut->GetParameter(0) );
	tf_data_dx_dyAntiCut_BG_scaled->FixParameter( 1, antiCut_BGscale_factor*tf_dx_dyAntiCut->GetParameter(1) );
	tf_data_dx_dyAntiCut_BG_scaled->FixParameter( 2, antiCut_BGscale_factor*tf_dx_dyAntiCut->GetParameter(2) );
	
	if( BG_polN >= 3 ){
		tf_data_dx_dyAntiCut_BG_scaled->FixParameter( 3, antiCut_BGscale_factor*tf_dx_dyAntiCut->GetParameter(3) );
	}
	if( BG_polN >= 4 ){
		tf_data_dx_dyAntiCut_BG_scaled->FixParameter( 4, antiCut_BGscale_factor*tf_dx_dyAntiCut->GetParameter(4) );
	}

	SBSkine.tf_data_dx_dyAntiCut_BG = tf_dx_dyAntiCut;
	SBSkine.tf_data_dx_dyAntiCut_BG_scaled = tf_data_dx_dyAntiCut_BG_scaled;

}

void fit_dyAnticut_and_BGsub( kine_ff_extract& SBSkine){

	int kine = SBSkine.kine;
	int sbsfieldscale = SBSkine.sbsfieldscale;
	TString run_target = SBSkine.run_target;

	TCanvas *c_anticut = new TCanvas("c_anticut", "c_anticut", 600, 500);

	int BG_polN = SBSkine.dx_dyAnticut_BG_polN;

	if( BG_polN == 0 ){
		cout << "No value specified for BG_polN for dy anticut BG fit. Defaulting to: pol2" << endl;
		BG_polN = 2;
		SBSkine.dx_dyAnticut_BG_polN = BG_polN;
	}

	Double_t dx_p_mean = SBSkine.dx_p_mean;
	Double_t dx_p_sigma = SBSkine.dx_p_sigma;
	Double_t dx_n_mean = SBSkine.dx_n_mean;
	Double_t dx_n_sigma = SBSkine.dx_n_sigma;
	Double_t dy_mean = SBSkine.dy_mean;
	Double_t dy_sigma = SBSkine.dy_sigma;

	double dx_plot_min = -2;
	double dx_plot_max = 1.0;

	double dy_mult = 3.0;
	Double_t dy_min = dy_mean - (dy_mult)*dy_sigma;
	Double_t dy_max = dy_mean + (dy_mult)*dy_sigma;

	TH1D *h_dx_dyAntiCut, *h_dx_dyBGsub;
	TH1D *h_dx_final_form = (TH1D*)SBSkine.h_data->Clone("h_dx_final_form");
	TH2D *h_dxdy = SBSkine.h_data_dxdy;
	TH2D *h_dxdy_dyAntiCut;
	TH2D *h_dxdy_dyCut;

	TAxis *xaxis = h_dxdy->GetXaxis();
	int xbins = xaxis->GetNbins();
	Double_t xaxis_min = xaxis->GetXmin();
	Double_t xaxis_max = xaxis->GetXmax();

	TAxis *yaxis = h_dxdy->GetYaxis();
	int ybins = yaxis->GetNbins();
	Double_t yaxis_min = yaxis->GetXmin();
	Double_t yaxis_max = yaxis->GetXmax();

//Let's first do a total fit of the final data_dx plot (best selections and final, etc.)
	h_dx_final_form->Draw();

	TF1 *dx_total_fit;

	if( BG_polN == 2 ){
		dx_total_fit = new TF1("dx_total_fit", fit_full_dx_peaks_with_pol2BG, dx_plot_min, dx_plot_max, BG_polN + 7);		
	}
	if( BG_polN == 3 ){
		dx_total_fit = new TF1("dx_total_fit", fit_full_dx_peaks_with_pol3BG, dx_plot_min, dx_plot_max, BG_polN + 7);		
	}
	if( BG_polN == 4 ){
		dx_total_fit = new TF1("dx_total_fit", fit_full_dx_peaks_with_pol3BG, dx_plot_min, dx_plot_max, BG_polN + 7);		
	}
	//Set reasonable parlimits
	dx_total_fit->SetParName(0, "dx_p_norm");
	dx_total_fit->SetParName(1, "dx_p_mean");
	dx_total_fit->SetParName(2, "dx_p_sigma");
	dx_total_fit->SetParName(3, "dx_n_norm");
	dx_total_fit->SetParName(4, "dx_n_mean");
	dx_total_fit->SetParName(5, "dx_n_sigma");
	dx_total_fit->SetParName(6, "dx_BG_pol_p0");
	dx_total_fit->SetParName(7, "dx_BG_pol_p1");
	dx_total_fit->SetParName(8, "dx_BG_pol_p2");
	if( BG_polN >= 3 ){
		dx_total_fit->SetParName(9, "dx_BG_pol_p3");	
	}
	if( BG_polN >= 4 ){
		dx_total_fit->SetParName(10, "dx_BG_pol_p4");
	}

	dx_total_fit->SetParLimits(0, 0.9*h_dx_final_form->GetMaximum(), 1.1*h_dx_final_form->GetMaximum());
	dx_total_fit->SetParLimits(1, 1.1*dx_p_mean, 0.9*dx_p_mean);
	dx_total_fit->SetParLimits(2, 0.7*dx_p_sigma, 1.1*dx_p_sigma);
	dx_total_fit->SetParLimits(3, 0.25*h_dx_final_form->GetMaximum(), 0.6*h_dx_final_form->GetMaximum());
	dx_total_fit->SetParLimits(4, dx_n_mean - 0.01, dx_n_mean + 0.01);
	dx_total_fit->SetParLimits(5, 0.8*dx_n_sigma, 1.25*dx_n_sigma);	
	dx_total_fit->SetParLimits(6, 200, 500);
	// dx_total_fit->SetParLimits(7, -140, -120);
	// dx_total_fit->SetParLimits(8, -100, -90);
	// dx_total_fit->SetParLimits(9, 0, 10);

	h_dx_final_form->Fit("dx_total_fit", "R+");

	SBSkine.tf_data_dx_total = dx_total_fit;

	TF1 *tf_dx_pinit = new TF1("tf_dx_pinit", fit_gaus, dx_plot_min, dx_plot_max, 3);
	tf_dx_pinit->FixParameter(0, dx_total_fit->GetParameter(0) );
	tf_dx_pinit->FixParameter(1, dx_total_fit->GetParameter(1) );
	tf_dx_pinit->FixParameter(2, dx_total_fit->GetParameter(2) );
	tf_dx_pinit->SetLineColor(kViolet);
	tf_dx_pinit->Draw("same");
	SBSkine.tf_data_dx_pinit = tf_dx_pinit;

	TF1 *tf_dx_ninit = new TF1("tf_dx_ninit", fit_gaus, dx_plot_min, dx_plot_max, 3);
	tf_dx_ninit->FixParameter(0, dx_total_fit->GetParameter(3) );
	tf_dx_ninit->FixParameter(1, dx_total_fit->GetParameter(4) );
	tf_dx_ninit->FixParameter(2, dx_total_fit->GetParameter(5) );
	tf_dx_ninit->SetLineColor(kMagenta);
	tf_dx_ninit->Draw("same");
	SBSkine.tf_data_dx_ninit = tf_dx_ninit;

	double BG_scale_factor = 0.5;

	TF1 *tf_dx_BGinit;

	if( BG_polN == 2 ){
		tf_dx_BGinit = new TF1("tf_dx_BGinit", background_pol2, dx_plot_min, dx_plot_max, BG_polN + 1);
	}
	if( BG_polN == 3 ){
		tf_dx_BGinit = new TF1("tf_dx_BGinit", background_pol3, dx_plot_min, dx_plot_max, BG_polN + 1);
	}
	if( BG_polN == 4 ){
		tf_dx_BGinit = new TF1("tf_dx_BGinit", background_pol4, dx_plot_min, dx_plot_max, BG_polN + 1);
	}
	tf_dx_BGinit->SetNpx(h_dx_final_form->GetNbinsX());
	tf_dx_BGinit->FixParameter(0, BG_scale_factor*dx_total_fit->GetParameter(6) );
	tf_dx_BGinit->FixParameter(1, BG_scale_factor*dx_total_fit->GetParameter(7) );
	tf_dx_BGinit->FixParameter(2, BG_scale_factor*dx_total_fit->GetParameter(8) );
	if( BG_polN >= 3 ){
		tf_dx_BGinit->FixParameter(3, dx_total_fit->GetParameter(9) );
	}
	if( BG_polN >= 4 ){
		tf_dx_BGinit->FixParameter(4, dx_total_fit->GetParameter(10) );
	}
	tf_dx_BGinit->SetLineStyle(6);
	tf_dx_BGinit->Draw("same");
	SBSkine.tf_data_dx_BGinit = tf_dx_BGinit;

	// TH1D *h_data_dx_BG_from_fit = (TH1D*)tf_dx_BGinit->GetHistogram();

//Let's set the region within the p and n peak regions to 0.
//We can use the values provided for dy_min and dy_max (these should correspond to some sigma parameter scaling of the dy_plot)
	h_dxdy_dyAntiCut = new TH2D("h_dxdy_dyAntiCut", "dxdy plot with p and n peak regions set to zero", xbins, xaxis_min, xaxis_max, ybins, yaxis_min, yaxis_max);

	for( int binx = 1; binx <= xbins; binx++ ){
		Double_t binCenter = xaxis->GetBinCenter(binx);

		if( binCenter < dy_min || binCenter > dy_max ){
			for( int biny = 1; biny <= ybins; biny++ ){
				h_dxdy_dyAntiCut->SetBinContent( binx, biny, h_dxdy->GetBinContent( binx, biny ));
			}
		}

		if( binCenter > dy_min && binCenter < dy_max ){
			for( int biny = 1; biny <= ybins; biny++ ){
				h_dxdy_dyAntiCut->SetBinContent( binx, biny, 0.0);
			}
		}
	}

	SBSkine.h_data_dxdy_dyAntiCut = h_dxdy_dyAntiCut;

	h_dx_dyAntiCut = (TH1D*)h_dxdy_dyAntiCut->ProjectionY();
	SBSkine.h_data_dx_dyAntiCut = h_dx_dyAntiCut;


//Now that this histogram is filled we can send it to the anti-cut BG fit function....
	dyAnticut_and_dxBG_fit(SBSkine);
	// h_dxdy_dyAntiCut->Draw();

	//BACKGROUND SUBTRACTION

	TH1D *h_data_dx_dyAntiCut_BG_scaled = (TH1D*)SBSkine.tf_data_dx_dyAntiCut_BG_scaled->GetHistogram();
	SBSkine.h_data_dx_dyAntiCut_BG_scaled = h_data_dx_dyAntiCut_BG_scaled;

	h_dx_dyBGsub = new TH1D("h_dx_dyBGsub", "Data dx with dyAntiCut BG subtracted", ybins, yaxis_min, yaxis_max);
	cout << "BG fit nbins: " << h_data_dx_dyAntiCut_BG_scaled->GetNbinsX() << endl;
	cout << "dx final form nbins: " << h_dx_final_form->GetNbinsX() << endl;
	for( int bin = 1; bin <= h_data_dx_dyAntiCut_BG_scaled->GetNbinsX(); bin++ ){
		h_dx_dyBGsub->SetBinContent(bin, h_dx_final_form->GetBinContent(bin) - h_data_dx_dyAntiCut_BG_scaled->GetBinContent(bin) );
	}
	h_dx_dyBGsub->GetXaxis()->SetRangeUser(yaxis_min, yaxis_max);
	h_dx_dyBGsub->SetMinimum(0);
	SBSkine.h_data_dx_BGsub = h_dx_dyBGsub;

}





#endif