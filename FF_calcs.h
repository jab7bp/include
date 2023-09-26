#ifndef FF_CALCS_H
#define FF_CALCS_H

double calculate_Q2_by_kine( int kine ){

	double theta_deg = lookup_BB_angle_by_kine(kine, "deg");
	double theta_rad = lookup_BB_angle_by_kine(kine, "rad");

	double Ebeam = lookup_beam_energy_from_kine(kine);
	double Escat_e = (Mp*Ebeam)/(Mp + Ebeam*( 1 - cos(theta_rad) ));

	double Q2 = 4.0*Escat_e*Ebeam*pow( sin(theta_rad/2.0), 2 );

	return Q2;
}

double calculate_mott( int kine ){
	double theta_deg = lookup_BB_angle_by_kine(kine, "deg");
	double theta_rad = lookup_BB_angle_by_kine(kine, "rad");
	double Ebeam = lookup_beam_energy_from_kine(kine);
	double Escat_e = (Mp*Ebeam)/(Mp + Ebeam*( 1 - cos(theta_rad) ));

	double mott = ( pow(alpha, 2)*pow( cos(theta_rad/2.0), 2 )*(Escat_e) )/( 4.0*pow(Ebeam, 3)*pow( sin(theta_rad/2.0), 4 ) );	

	return mott;
}

double calc_GMp_param( double Q2 ){
	double tau_p = Q2/(4.0*Mp*Mp);

	double aMp1, bMp1, bMp2, bMp3;
	aMp1 = 1.09;
	bMp1 = 12.31;
	bMp2 = 25.57;
	bMp3 = 30.61;

	double calc_GMp_param = (mu_p)*( 1 + aMp1*tau_p )/( 1 + bMp1*tau_p + bMp2*pow(tau_p, 2) + bMp3*pow(tau_p, 3) );

	return calc_GMp_param;
}

double calc_GEp_param( double Q2 ){
	double tau_p = Q2/(4.0*Mp*Mp);

	double aEp1, bEp1, bEp2, bEp3;
	aEp1 = -0.19;
	bEp1 = 11.12;
	bEp2 = 15.16;
	bEp3 = 21.25;

	double calc_GEp_param = ( 1 + aEp1*tau_p )/( 1 + bEp1*tau_p + bEp2*pow(tau_p, 2) + bEp3*pow(tau_p, 3) );

	return calc_GEp_param;
}

double calc_GEn_param( double Q2 ){
	double aEn1, bEn1, bEn2;
	aEn1 = -0.1;
	bEn1 = 2.83;
	bEn2 = 0.43;

	double calc_GEn_param = ( aEn1 )/( pow( (1 + bEn1*Q2), 2 ) ) - ( aEn1 )/( pow( (1 + bEn2*Q2), 2 ) );

	return calc_GEn_param;
}

double calc_GMn_param( double Q2 ){
	double tau_n = Q2/(4.0*Mn*Mn);

	double aMn1, bMn1, bMn2, bMn3;
	aMn1 = 8.28;
	bMn1 = 21.3;
	bMn2 = 77.0;
	bMn3 = 238.0;

	double calc_GMn_param = (mu_n)*( 1 + aMn1*tau_n )/( 1 + bMn1*tau_n + bMn2*pow(tau_n, 2) + bMn3*pow(tau_n, 3) );

	return calc_GMn_param;
}

double calc_reduced_CS( double eps, double GE, double tau, double GM ){
	double sigma_r = eps*pow(GE, 2) + tau*pow(GM, 2);
	return sigma_r;
}

#endif