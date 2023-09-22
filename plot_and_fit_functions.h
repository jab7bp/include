#ifndef PLOT_AND_FIT_FUNCS_H
#define PLOT_AND_FIT_FUNCS_H

TLegend *tl_dx_mag0, *tl_dx_mag, *tl_dy_mag, *tl_dx_mag_dyBG;
TF1 *tf_data_dx_p, *tf_data_dx_n;
TF1 *tf_dx_bg, *fit_dx_dyBG, *fit_polBG_with_reject, *fit_dx, *fit_dx_origBG;
TF1 *tf_BG_dx_pol4, *tf_BG_dx_pol5, *tf_BG_dx_pol6, *tf_BG_dx_pol8;
TF1 *tf_BG_dx_wcut_pol4, *tf_BG_dx_wcut_pol5, *tf_BG_dx_wcut_pol6, *tf_BG_dx_wcut_pol8;
Double_t par_BG_reject[5], par_dy_bg_reject[5], parPoly[9], parPol8[9], parPol4[5], parPol5[6], parPol6[7], parPol7[8], parElastics[6];

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

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
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







#endif