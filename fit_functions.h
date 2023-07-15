#ifndef FIT_FUNCTIONS_H
#define FIT_FUNCTIONS_H

// TGraph *gr_gaus1, *gr_polyBG;

// double dy_reject_min = -0.20;
// double dy_reject_max = 0.175;
double dy_reject_min = 0.4;
double dy_reject_max = 1.21;
double expo_and_polN_N;
// int polN = 6;
// int dypolN = 8;

Double_t fit_det_eff( Double_t *x, Double_t *par){
	double det_eff_log = 0.0;

	det_eff_log = par[0] - log(par[1]/x[0]);

	return det_eff_log;
}

Double_t fit_expo_and_polN( Double_t *x, Double_t *par){
	double expo = 0.0;
	double polN = 0.0;

	if( x[0] > dy_reject_max ){
		expo = par[0] - (par[0] - par[1])*exp(-par[2]*x[0]);	
		polN = 0.0;	
	}

	else{
		for( int N = 0; N <= expo_and_polN_N; N++ ){
			polN += par[N+3]*pow(x[0], N);
		}	
		expo = 0.0;	
	}

	return expo + polN;

}

Double_t fit_expo_and_polN_with_reject( Double_t *x, Double_t *par){
	double expo = 0.0;
	double polN = 0.0;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		expo = par[0] - (par[0] - par[1])*exp(-par[2]*x[0]);

		for( int N = 0; N <= expo_and_polN_N; N++ ){
			polN += par[N+3]*pow(x[0], N);
		}		
		return expo + polN;
	}
}

Double_t fit_expo(Double_t *x, Double_t *par){
	// cout << "Fiting with exponential" << endl;
	double expo = 0.0;
	expo = par[0] - (par[0] - par[1])*exp(-par[2]*x[0]);
	return expo;
}

Double_t fit_expo_with_reject(Double_t *x, Double_t *par){
	double expo = 0.0;

	// double reject_min = 0.50217680;
	// // double reject_max = 1.09;
	// double reject_max = 1.15;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		expo = par[0] - (par[0] - par[1])*exp(-par[2]*x[0]);
		return expo;
	}
}

Double_t fitPol0(Double_t *x, Double_t *par){
	return par[0];
}

Double_t fitPol1(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0];
}

Double_t fitPol2(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t fitPol3(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

Double_t fitPol4(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t fitPol5(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5);
}

Double_t fitPol6(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6);
}

Double_t fitPol7(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7);
}

Double_t fitPol8(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7) + par[8]*pow(x[0], 8);
}

Double_t fitPol10(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7) + par[8]*pow(x[0], 8) + par[9]*pow(x[0], 9) + par[10]*pow(x[0], 10);
}

Double_t fitPol3_with_reject(Double_t *x, Double_t *par){
	
	// double reject_min = 0.50217680;
	// // double reject_max = 1.09;
	// double reject_max = 1.05;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
	}
}

Double_t fitPol4_with_reject(Double_t *x, Double_t *par){
	
	// double reject_min = 0.50217680;
	// double reject_max = 1.075;
	// double reject_max = 1.2;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
	}
}

Double_t dyfitPol4_with_reject(Double_t *x, Double_t *par){


	if( x[0] > (dy_reject_min) && x[0] < (dy_reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];		
	}

}

Double_t fitPol5_with_reject(Double_t *x, Double_t *par){
	
	// double reject_min = 0.50217680;
	// // double reject_max = 1.09;
	// double reject_max = 1.05;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0] + par[5]*x[0]*x[0]*x[0]*x[0]*x[0];
	}
}

Double_t fitPol6_with_reject(Double_t *x, Double_t *par){

	// double reject_min = 0.50217680;
	// double reject_max = 1.09;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6);
	}
}

Double_t dyfitPol6_with_reject(Double_t *x, Double_t *par){

	// double reject_min = -0.20;
	// double reject_max = 0.15;

	if( x[0] > (dy_reject_min) && x[0] < (dy_reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6);
	}
}

Double_t fitPol8_with_reject(Double_t *x, Double_t *par){

	// double reject_min = 0.50217680;
	// double reject_max = 1.09;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7) + par[8]*pow(x[0], 8);
	}
}

Double_t dyfitPol8_with_reject(Double_t *x, Double_t *par){

	// double reject_min = -0.20;
	// double reject_max = 0.15;

	if( x[0] > (dy_reject_min) && x[0] < (dy_reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7) + par[8]*pow(x[0], 8);
	}
}

Double_t dyfitPol10_with_reject(Double_t *x, Double_t *par){

	// double reject_min = -0.20;
	// double reject_max = 0.15;

	if( x[0] > (dy_reject_min) && x[0] < (dy_reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7) + par[8]*pow(x[0], 8) + par[9]*pow(x[0], 9) + par[10]*pow(x[0], 10);
	}
}

Double_t W2_BG_fit_with_reject( Double_t *x, Double_t *par ){
	if( x[0] > 0.55 && x[0] < 1.2 ){
		TF1::RejectPoint();
		return 0;
	}
	else{
		return exp(par[0] +par[1]*x[0]);
	}
}

Double_t W2_BG_fit_NO_reject( Double_t *x, Double_t *par ){
	return exp(par[0] +par[1]*x[0]);
}

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

Double_t fullFitFunction_reject( Double_t *x, Double_t *par ){
	// return fit_gaus(x, par) + W2_BG_fit_NO_reject(x, &par[3]);
	return fit_gaus(x, par) + fitPol6_with_reject(x, &par[3]);
	// return fit_gaus(x, par) + fitPol8_with_reject(x, &par[3]);
}

Double_t fullFitFunction( Double_t *x, Double_t *par ){
	// return fit_gaus(x, par) + W2_BG_fit_NO_reject(x, &par[3]);

	if( polN == 99 ){
		return fit_gaus(x, par) + fit_expo(x, &par[3]);		
	}
	if( polN == 3 ){
		return fit_gaus(x, par) + fitPol3(x, &par[3]);		
	}
	if( polN == 4 ){
		return fit_gaus(x, par) + fitPol4(x, &par[3]);		
	}
	// return fit_gaus(x, par) + fitPol5(x, &par[3]);
	if( polN == 5 ){
		return fit_gaus(x, par) + fitPol5(x, &par[3]);		
	}
	if( polN == 6 ){
		return fit_gaus(x, par) + fitPol6(x, &par[3]);		
	}
	if( polN == 8 ){
		return fit_gaus(x, par) + fitPol8(x, &par[3]);		
	}
	// return fit_gaus(x, par) + fitPol10(x, &par[3]);
}

Double_t fullBGCutFitFunction( Double_t *x, Double_t *par ){

	return fit_gaus(x, par) + fitPol8(x, &par[3]);		

}

Double_t fullBGSubFitFunction( Double_t *x, Double_t *par ){

	return fit_gaus(x, par) + fitPol8(x, &par[3]);		

}

Double_t fullInterpolFitFunction( Double_t *x, Double_t *par ){
	// return fit_gaus(x, par) + W2_BG_fit_NO_reject(x, &par[3]);

	if( interpolN == 99 ){
		return fit_gaus(x, par) + fit_expo(x, &par[3]);		
	}
	if( interpolN == 3 ){
		return fit_gaus(x, par) + fitPol3(x, &par[3]);		
	}
	if( interpolN == 4 ){
		return fit_gaus(x, par) + fitPol4(x, &par[3]);		
	}
	// return fit_gaus(x, par) + fitPol5(x, &par[3]);
	if( interpolN == 5 ){
		return fit_gaus(x, par) + fitPol5(x, &par[3]);		
	}
	if( interpolN == 6 ){
		return fit_gaus(x, par) + fitPol6(x, &par[3]);		
	}
	if( interpolN == 8 ){
		return fit_gaus(x, par) + fitPol8(x, &par[3]);		
	}
	// return fit_gaus(x, par) + fitPol10(x, &par[3]);
}

Double_t fullDXFitFunction( Double_t *x, Double_t *par ){

	if( dxpolN == 2){
		return fit_gaus(x, par) + fitPol2(x, &par[3]);		
	}

	if( dxpolN == 3){
		return fit_gaus(x, par) + fitPol2(x, &par[3]);		
	}
	// return fitPol3(x, par) + fitPol0(x, &par[4]) + fitPol2(x, &par[5]);
	// if( dypolN == 4 ){
	// 	return fit_gaus(x, par) + fitPol4(x, &par[3]);		
	// }
	// if( dypolN == 6 ){
	// 	return fit_gaus(x, par) + fitPol6(x, &par[3]);		
	// }
	// if( dypolN == 8 ){
	// 	return fit_gaus(x, par) + fitPol8(x, &par[3]);		
	// }

}

Double_t fullDYFitFunction( Double_t *x, Double_t *par ){

	if( dypolN == 2 ){
		return fit_gaus(x, par) + fitPol2(x, &par[3]);
	}

	if( dypolN == 3 ){
		return fit_gaus(x, par) + fitPol3(x, &par[3]);
	}

	// return fitPol3(x, par) + fitPol0(x, &par[4]) + fitPol2(x, &par[5]);
	// if( dypolN == 4 ){
	// 	return fit_gaus(x, par) + fitPol4(x, &par[3]);		
	// }
	// if( dypolN == 6 ){
	// 	return fit_gaus(x, par) + fitPol6(x, &par[3]);		
	// }
	// if( dypolN == 8 ){
	// 	return fit_gaus(x, par) + fitPol8(x, &par[3]);		
	// }

}

Double_t plot_gaus(Double_t *par){
	Double_t fit_x[2400], gaus1_y[2400];
	Int_t fit_n = 800;
	int fit_cnt = 0;

	for(double i = -2.5; i < 2.5; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	gr_gaus1->SetLineColor(1);
	gr_gaus1->SetLineStyle(7);

	gr_gaus1->Draw("same");
	
	return 1;
}

Double_t plot_sig_and_BG(Double_t *par){
	Double_t fit_x[180], gaus1_y[180], polyBG_y[180];
	Int_t fit_n = 180;
	int fit_cnt = 0;

	for(double i = 0.4; i < 1.4; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		polyBG_y[fit_cnt] = par[11]*pow(i, 8) + par[10]*pow(i, 7) + par[9]*pow(i, 6) + par[8]*pow(i, 5) + par[7]*pow(i, 4) + par[6]*pow(i, 3) + par[5]*pow(i, 2) + par[4]*i + par[3];
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	// cout << "Max element of sig: " << gr_gaus1->GetMaximum();
	gr_gaus1->SetLineColor(1);
	gr_gaus1->SetLineStyle(7);
	TGraph *gr_polyBG = new TGraph(fit_n, fit_x, polyBG_y);
	gr_polyBG->SetLineColor(7);

	gr_gaus1->Draw("same");
	gr_polyBG->Draw("same");
	
	return 1;

}

Double_t plot_polN(Double_t *par, Int_t polOrd){
	const int fit_n = 180;
	Double_t fit_x[fit_n], polN_y[fit_n];
	int fit_cnt = 0;

	for(double i = 0.4; i < 1.4; i += 0.01 ){
		fit_x[fit_cnt] = i;
		for( int ord = 0; ord <= polOrd; ord ++){
			polN_y[fit_cnt] += par[ord]*pow(i, ord);
		}
		fit_cnt++;
	}

	TGraph *gr_polN = new TGraph(fit_n, fit_x, polN_y);
	gr_polN->SetLineColor(3);
	gr_polN->Draw("same");
	return 1;

}


#endif