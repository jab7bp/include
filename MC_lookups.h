#ifndef MC_LOOKUPS_H
#define MC_LOOKUPS_H

double lookup_MC_dxdy(int kine, int sbsfieldscale, TString run_targ, TString param){
	double lookup_val = -1;

	//dxdy vector:
// 	kine, sbsfieldscale, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma

	vector<vector<double>> MC_LD2_dxdy = {
//    kine, mag,    dx_p,      dx_p_sig,       dx_n,       dx_n_sig,       dy,         dy_sig
		{4, 0, },
		{4, 30, },
		{4, 50, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{8, 70, -8.79633e-01, 1.62141e-01, -1.12461e-02, 1.46882e-01, -2.27440e-02, 2.76191e-03},
		// {8, 70, -8.73151e-01, 1.74421e-01, -2.41303e-02, 1.52362e-01, 6.31940e-02, 2.26323e-01},
		{9, 70, -8.55327e-01, 1.49493e-01, -5.06166e-02, 2.00000e-01, 6.87480e-02, 1.79564e-01}
		// {8, 70, -0.849403, 0.14, -0.00315522, 0.14, -0.0463545, 0.23}
		// {8, 70, -0.852848, 0.170124, 2.96250e-04, 0.160995, -2.56391e-02, 1.95379e-01}

	};

	if( run_targ == "LD2" ){
		for( size_t i = 0; i < MC_LD2_dxdy.size(); i++){
			if( MC_LD2_dxdy[i][0] == kine && MC_LD2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = MC_LD2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = MC_LD2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = MC_LD2_dxdy[i][4];}
				if( param == "dx_n_sigma" ){lookup_val = MC_LD2_dxdy[i][5];}
				if( param == "dy" ){lookup_val = MC_LD2_dxdy[i][6];}
				if( param == "dy_sigma" ){lookup_val = MC_LD2_dxdy[i][7];}

			}
		}
	}

	//kine, mag, dx_p, dx_p_sigma, dy, dy_sigma
	vector<vector<double>> MC_LH2_dxdy = {
		{4, 0, },
		{4, 30, },
		{8, 70, -1.137277, 0.0855973, 0.0054094, 0.059083}, 
		{9, 70, -1.137277, 0.0855973, 0.0054094, 0.059083}, 
	};

	if( run_targ == "LH2" ){
		for( size_t i = 0; i < MC_LH2_dxdy.size(); i++){
			if( MC_LD2_dxdy[i][0] == kine && MC_LH2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = MC_LH2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = MC_LH2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = -1;}
				if( param == "dx_n_sigma" ){lookup_val = -1;}
				if( param == "dy" ){lookup_val = MC_LH2_dxdy[i][4];}
				if( param == "dy_sigma" ){lookup_val = MC_LH2_dxdy[i][5];}

			}
		}	
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;

}

double lookup_pass2_MC_dxdy(int kine, int sbsfieldscale, TString run_targ, TString param){
	double lookup_val = -1;

	//dxdy vector:
// 	kine, sbsfieldscale, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma

	vector<vector<double>> MC_LD2_dxdy = {
//    kine, mag,    dx_p,      dx_p_sig,       dx_n,       dx_n_sig,       dy,         dy_sig
		{4, 0, },
		{4, 30, -6.85956e-01, 1.83169e-01, -2.69853e-02, 1.89613e-01, -1.73561e-02, 2.16250e-01},
		{4, 50, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{8, 70, -8.55567e-01, 1.63631e-01, -8.74682e-03, 1.43677e-01, -1.23321e-02, 2.05283e-01},
		// {8, 70, -8.73151e-01, 1.74421e-01, -2.41303e-02, 1.52362e-01, 6.31940e-02, 2.26323e-01},
		{9, 70, -8.55327e-01, 1.49493e-01, -5.06166e-02, 2.00000e-01, 6.87480e-02, 1.79564e-01}
		// {8, 70, -0.849403, 0.14, -0.00315522, 0.14, -0.0463545, 0.23}
		// {8, 70, -0.852848, 0.170124, 2.96250e-04, 0.160995, -2.56391e-02, 1.95379e-01}

	};

	if( run_targ == "LD2" ){
		for( size_t i = 0; i < MC_LD2_dxdy.size(); i++){
			if( MC_LD2_dxdy[i][0] == kine && MC_LD2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = MC_LD2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = MC_LD2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = MC_LD2_dxdy[i][4];}
				if( param == "dx_n_sigma" ){lookup_val = MC_LD2_dxdy[i][5];}
				if( param == "dy" ){lookup_val = MC_LD2_dxdy[i][6];}
				if( param == "dy_sigma" ){lookup_val = MC_LD2_dxdy[i][7];}

			}
		}
	}

	//kine, mag, dx_p, dx_p_sigma, dy, dy_sigma
	vector<vector<double>> MC_LH2_dxdy = {
		{4, 0, },
		{4, 30, },
		{8, 70, -1.137277, 0.0855973, 0.0054094, 0.059083}, 
		{9, 70, -1.137277, 0.0855973, 0.0054094, 0.059083}, 
	};

	if( run_targ == "LH2" ){
		for( size_t i = 0; i < MC_LH2_dxdy.size(); i++){
			if( MC_LD2_dxdy[i][0] == kine && MC_LH2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = MC_LH2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = MC_LH2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = -1;}
				if( param == "dx_n_sigma" ){lookup_val = -1;}
				if( param == "dy" ){lookup_val = MC_LH2_dxdy[i][4];}
				if( param == "dy_sigma" ){lookup_val = MC_LH2_dxdy[i][5];}

			}
		}	
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;

}

double lookup_MC_cut( int kine, int sbsfield, TString run_targ, TString param ){
	double lookup_val = -1;

	int target_int;

	if( run_targ == "LH2" ){
		target_int = 0;
	}
	if( run_targ == "LD2" ){
		target_int = 1;
	}

	//Cuts vector:
// 	[0]target_int, [1]kine, [2]sbs_field [3]PS_clus_e_cut, [4]SH_PS_mean, [5]SH_PS_sigma, [6]HCal_clus_e_cut, [7]Ep, [8]Ep_sigma, [9]W2, [10]W2_sigma, [11]W, [12]W_sigma, 
	vector<vector<double>> MC_cuts = {
	//tar, kine, mag, PS,  SH_PS,  SH_PS_sig,   HCal,    Ep,        Ep_sigma,      W2,        W2_sig,          W,        W_sig
		{1, 8, 70, 0.200, 3.33556, 3.31747e-01, 0.02, 9.55274e-01, 7.66726e-02, 9.08361e-01, 3.08407e-01, 9.73196e-01, 1.58584e-01},
		{1, 9, 70, 0.150, 3.24000, 0.417381,    0.02, 0.960000,    0.0750000,   0.979146,    3.22005e-01,    1.005,       0.195985},
		{0, 9, 70, 0.150, 1.64755, 0.219328, 0.01, 0.718334, 0.0908741}
	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-run_02

	for( size_t i = 0; i < MC_cuts.size(); i++){
		if( MC_cuts[i][0] == target_int && MC_cuts[i][1] == kine && MC_cuts[i][2] == sbsfield ){
			if( param == "PS_clus_e_cut" ){lookup_val = MC_cuts[i][3];}
			if( param == "SH_PS_mean" ){lookup_val = MC_cuts[i][4];}
			if( param == "SH_PS_sigma" ){lookup_val = MC_cuts[i][5];}
			if( param == "HCal_clus_e_cut" ){lookup_val = MC_cuts[i][6];}
			if( param == "Ep" ){lookup_val = MC_cuts[i][7];}
			if( param == "Ep_sigma" ){lookup_val = MC_cuts[i][8];}
			if( param == "W2" ){lookup_val = MC_cuts[i][9];}
			if( param == "W2_sigma" ){lookup_val = MC_cuts[i][10];}
			if( param == "W" ){lookup_val = MC_cuts[i][11];}
			if( param == "W_sigma" ){lookup_val = MC_cuts[i][12];}

		}
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;
}


double lookup_simc_dxdy(int kine, int sbsfieldscale, TString run_targ, TString param){
	double lookup_val = -1;

	//dxdy vector:
// 	kine, sbsfieldscale, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma

	vector<vector<double>> simc_LD2_dxdy = {
//    kine, mag,    dx_p,      dx_p_sig,       dx_n,       dx_n_sig,       dy,         dy_sig
		{4, 0, },
		{4, 30, },
		{4, 50, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{8, 70, -9.32040e-01, 1.58298e-01, -5.89213e-02, 1.52057e-01, -0.0440112, 0.23},
		{9, 70, -8.80616e-01, 1.67226e-01, -5.09922e-03, 1.71203e-01, -2.80750e-02, 1.81892e-01}
		// {8, 70, -0.913432, 0.14, -0.0622097, 0.14, -0.0322545, 0.219831}
		// {8, 70, -0.852848, 0.170124, 2.96250e-04, 0.160995, -2.56391e-02, 1.95379e-01}

	};

	if( run_targ == "LD2" ){
		for( size_t i = 0; i < simc_LD2_dxdy.size(); i++){
			if( simc_LD2_dxdy[i][0] == kine && simc_LD2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = simc_LD2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = simc_LD2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = simc_LD2_dxdy[i][4];}
				if( param == "dx_n_sigma" ){lookup_val = simc_LD2_dxdy[i][5];}
				if( param == "dy" ){lookup_val = simc_LD2_dxdy[i][6];}
				if( param == "dy_sigma" ){lookup_val = simc_LD2_dxdy[i][7];}

			}
		}
	}

	//kine, mag, dx_p, dx_p_sigma, dy, dy_sigma
	vector<vector<double>> simc_LH2_dxdy = {
		{4, 0, },
		{4, 30, },
		{8, 70, -1.137277, 0.0855973, 0.0054094, 0.059083}, 
	};

	if( run_targ == "LH2" ){
		for( size_t i = 0; i < simc_LH2_dxdy.size(); i++){
			if( simc_LD2_dxdy[i][0] == kine && simc_LH2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = simc_LH2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = simc_LH2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = -1;}
				if( param == "dx_n_sigma" ){lookup_val = -1;}
				if( param == "dy" ){lookup_val = simc_LH2_dxdy[i][4];}
				if( param == "dy_sigma" ){lookup_val = simc_LH2_dxdy[i][5];}

			}
		}	
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;

}

double lookup_pass2_simc_dxdy(int kine, int sbsfieldscale, TString run_targ, TString param){
	double lookup_val = -1;

	//dxdy vector:
// 	kine, sbsfieldscale, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma

	vector<vector<double>> simc_LD2_dxdy = {
//    kine, mag,    dx_p,      dx_p_sig,       dx_n,       dx_n_sig,       dy,         dy_sig
		{4, 0, },
		{4, 30, -6.85956e-01, 1.83169e-01, -2.69853e-02, 1.89613e-01, -1.73561e-02, 2.16250e-01},
		{4, 50, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{8, 70, -8.55567e-01, 1.63631e-01, -8.74682e-03, 1.43677e-01, -1.23321e-02, 2.05283e-01},
		{9, 70, -8.67616e-01, 1.52463e-01, -4.13039e-03, 1.42609e-01, -2.95006e-02, 1.65635e-01}
		// {8, 70, -0.913432, 0.14, -0.0622097, 0.14, -0.0322545, 0.219831}
		// {8, 70, -0.852848, 0.170124, 2.96250e-04, 0.160995, -2.56391e-02, 1.95379e-01}

	};

	if( run_targ == "LD2" ){
		for( size_t i = 0; i < simc_LD2_dxdy.size(); i++){
			if( simc_LD2_dxdy[i][0] == kine && simc_LD2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = simc_LD2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = simc_LD2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = simc_LD2_dxdy[i][4];}
				if( param == "dx_n_sigma" ){lookup_val = simc_LD2_dxdy[i][5];}
				if( param == "dy" ){lookup_val = simc_LD2_dxdy[i][6];}
				if( param == "dy_sigma" ){lookup_val = simc_LD2_dxdy[i][7];}

			}
		}
	}

	//kine, mag, dx_p, dx_p_sigma, dy, dy_sigma
	vector<vector<double>> simc_LH2_dxdy = {
		{4, 0, },
		{4, 30, },
		{8, 70, -1.137277, 0.0855973, 0.0054094, 0.059083}, 
	};

	if( run_targ == "LH2" ){
		for( size_t i = 0; i < simc_LH2_dxdy.size(); i++){
			if( simc_LD2_dxdy[i][0] == kine && simc_LH2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = simc_LH2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = simc_LH2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = -1;}
				if( param == "dx_n_sigma" ){lookup_val = -1;}
				if( param == "dy" ){lookup_val = simc_LH2_dxdy[i][4];}
				if( param == "dy_sigma" ){lookup_val = simc_LH2_dxdy[i][5];}

			}
		}	
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;

}

double lookup_simc_dxdy_with_portField(int kine, int sbsfieldscale, TString run_targ, TString portField, TString param){
	double lookup_val = -1;
	double portField_val = portField.Atof();
	// cout << "portfield val: " << portField_val << endl;

	//dxdy vector:
// 	kine, sbsfieldscale, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma

	vector<vector<double>> simc_LD2_dxdy = {
//    kine, mag, portField, dx_p,      dx_p_sig,       dx_n,       dx_n_sig,       dy,         dy_sig
		{4, 0, },
		{4, 30, 387.0, -6.88539e-01, 1.85003e-01, -2.93900e-02, 1.94680e-01, -1.83491e-02, 2.17775e-01},
		{4, 30, 382.0, -6.77375e-01, 1.73755e-01, -3.24282e-02, 2.00758e-01, -2.34307e-02, 2.14616e-01},
		{4, 30, 380.0, -6.73803e-01, 1.80103e-01, -2.59548e-02, 1.93379e-01, -2.01166e-02, 2.11148e-01},
		{4, 30, 378.0, -6.69676e-01, 1.86114e-01, -3.08786e-02, 1.94588e-01, -2.16924e-02, 2.05162e-01},
		{4, 30, 376.0, -6.68969e-01, 1.80815e-01, -3.31734e-02, 1.94083e-01, -2.14107e-02, 2.08351e-01},
		{4, 50, 630.0, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{4, 50, 635.0, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{4, 50, 640.0, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{4, 50, 645.0, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{4, 50, 650.0, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{4, 50, 655.0, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
		{8, 70, 635.0, -8.51918e-01, 1.57751e-01, -7.21366e-03, 1.44332e-01, -1.37180e-02, 1.99796e-01},
		{8, 70, 637.0, -8.57597e-01, 1.58685e-01, -7.50604e-03, 1.40883e-01, -1.43177e-02, 2.02419e-01},
		{8, 70, 653.0, -8.76133e-01, 1.59905e-01, -7.56828e-03, 1.44230e-01, -3.08164e-02, 2.26097e-01},
		{8, 70, 660.0, -8.81911e-01, 1.55295e-01, -8.48764e-03, 1.46549e-01, -2.65032e-02, 1.94552e-01},
		{9, 70, 653.0, -8.63129e-01, 1.41754e-01, -3.08180e-03, 1.39902e-01, -2.94542e-02, 1.63196e-01},
		{9, 70, 657.0, -8.67735e-01, 1.43540e-01, -3.23057e-03, 1.45355e-01, -3.07056e-02, 1.67445e-01},
		{9, 70, 666.0, -8.80075e-01, 1.48067e-01, -4.13151e-03, 1.43590e-01, -2.97798e-02, 1.65274e-01},
		{9, 70, 670.0, -8.90092e-01, 1.44202e-01, -9.96894e-04, 1.36398e-01, -2.52938e-02, 1.60053e-01}

		// {8, 70, -0.913432, 0.14, -0.0622097, 0.14, -0.0322545, 0.219831}
		// {8, 70, -0.852848, 0.170124, 2.96250e-04, 0.160995, -2.56391e-02, 1.95379e-01}

	};

	if( run_targ == "LD2" ){
		for( size_t i = 0; i < simc_LD2_dxdy.size(); i++){
			if( simc_LD2_dxdy[i][0] == kine && simc_LD2_dxdy[i][1] == sbsfieldscale && simc_LD2_dxdy[i][2] == portField_val ){
				if( param == "dx_p" ){lookup_val = simc_LD2_dxdy[i][3];}
				if( param == "dx_p_sigma" ){lookup_val = simc_LD2_dxdy[i][4];}
				if( param == "dx_n" ){lookup_val = simc_LD2_dxdy[i][5];}
				if( param == "dx_n_sigma" ){lookup_val = simc_LD2_dxdy[i][6];}
				if( param == "dy" ){lookup_val = simc_LD2_dxdy[i][7];}
				if( param == "dy_sigma" ){lookup_val = simc_LD2_dxdy[i][8];}

			}
		}
	}

	//kine, mag, dx_p, dx_p_sigma, dy, dy_sigma
	vector<vector<double>> simc_LH2_dxdy = {
		{4, 0, },
		{4, 30, },
		{8, 70, -1.137277, 0.0855973, 0.0054094, 0.059083}, 
	};

	if( run_targ == "LH2" ){
		for( size_t i = 0; i < simc_LH2_dxdy.size(); i++){
			if( simc_LD2_dxdy[i][0] == kine && simc_LH2_dxdy[i][1] == sbsfieldscale ){
				if( param == "dx_p" ){lookup_val = simc_LH2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = simc_LH2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = -1;}
				if( param == "dx_n_sigma" ){lookup_val = -1;}
				if( param == "dy" ){lookup_val = simc_LH2_dxdy[i][4];}
				if( param == "dy_sigma" ){lookup_val = simc_LH2_dxdy[i][5];}

			}
		}	
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;

}

double lookup_simc_cut( int kine, int sbsfield, TString run_targ, TString param ){
	double lookup_val = -1;

	int target_int;

	if( run_targ == "LH2" ){
		target_int = 0;
	}
	if( run_targ == "LD2" ){
		target_int = 1;
	}

	//Cuts vector:
// 	[0]target_int, [1]kine, [2]sbs_field [3]PS_clus_e_cut, [4]SH_PS_mean, [5]SH_PS_sigma, [6]HCal_clus_e_cut, [7]Ep, [8]Ep_sigma, [9]W2, [10]W2_sigma, [11]W, [12]W_sigma, 
	vector<vector<double>> simc_cuts = {
	//targ, kine,  mag,  PS,   SHPS,   SHPS_sig, HCal_e, Ep,       Ep_sigma,     W2,      W2_sig,     W,      W_sig
		{1,  4,    30, 0.150, 2.24380, 0.287984, 0.0200, 1.044270, 0.1116450, 0.927754, 0.194505, 0.968430, 0.100385},
		{1,  8,    70, 0.150, 3.25473, 0.419763, 0.0200, 0.956448, 0.0780319, 0.909982, 0.306327, 0.961529, 0.185036},
		{1,  9,    70, 0.150, 1.56300, 0.291986, 0.0200, 1.037640, 0.1315590, 0.956325, 0.280737, 1.003870, 0.164814}
	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-run_02

	for( size_t i = 0; i < simc_cuts.size(); i++){
		if( simc_cuts[i][0] == target_int && simc_cuts[i][1] == kine && simc_cuts[i][2] == sbsfield ){
			if( param == "PS_clus_e_cut" ){lookup_val = simc_cuts[i][3];}
			if( param == "SH_PS_mean" ){lookup_val = simc_cuts[i][4];}
			if( param == "SH_PS_sigma" ){lookup_val = simc_cuts[i][5];}
			if( param == "HCal_clus_e_cut" ){lookup_val = simc_cuts[i][6];}
			if( param == "Ep" ){lookup_val = simc_cuts[i][7];}
			if( param == "Ep_sigma" ){lookup_val = simc_cuts[i][8];}
			if( param == "W2" ){lookup_val = simc_cuts[i][9];}
			if( param == "W2_sigma" ){lookup_val = simc_cuts[i][10];}
			if( param == "W" ){lookup_val = simc_cuts[i][11];}
			if( param == "W_sigma" ){lookup_val = simc_cuts[i][12];}

		}
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;
}

double lookup_pass2_simc_cut( int kine, int sbsfield, TString run_targ, TString param ){
	double lookup_val = -1;

	int target_int;

	if( run_targ == "LH2" ){
		target_int = 0;
	}
	if( run_targ == "LD2" ){
		target_int = 1;
	}

	//Cuts vector:
// 	[0]target_int, [1]kine, [2]sbs_field [3]PS_clus_e_cut, [4]SH_PS_mean, [5]SH_PS_sigma, [6]HCal_clus_e_cut, [7]Ep, [8]Ep_sigma, [9]W2, [10]W2_sigma, [11]W, [12]W_sigma, 
	vector<vector<double>> simc_cuts = {
	//targ, kine,  mag,  PS,   SHPS,   SHPS_sig,    HCal_e, Ep,       Ep_sigma,        W2,      W2_sig,          W,      W_sig
		{1,  4,    30, 0.150, 3.43217, 4.24507e-01, 0.0050, 1.04396, 1.12542e-01, 9.23760e-01, 1.87922e-01, 9.68462e-01, 1.02889e-01},
		{1,  4,    50, 0.150, 3.43217, 4.24507e-01, 0.0050, 1.04396, 1.12542e-01, 9.23760e-01, 1.87922e-01, 9.68462e-01, 1.02889e-01},
		{1,  8,    70, 0.150, 3.43217, 4.24507e-01, 0.0050, 1.00237, 8.69793e-02, 9.59612e-01, 2.80758e-01, 9.97310e-01, 1.51970e-01},
		{1,  9,    70, 0.150, 1.78846, 2.52708e-01, 0.0050, 1.07721, 1.32326e-01, 9.33520e-01, 2.45559e-01, 9.77031e-01, 1.31265e-01}
	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-run_02

	for( size_t i = 0; i < simc_cuts.size(); i++){
		if( simc_cuts[i][0] == target_int && simc_cuts[i][1] == kine && simc_cuts[i][2] == sbsfield ){
			if( param == "PS_clus_e_cut" ){lookup_val = simc_cuts[i][3];}
			if( param == "SH_PS_mean" ){lookup_val = simc_cuts[i][4];}
			if( param == "SH_PS_sigma" ){lookup_val = simc_cuts[i][5];}
			if( param == "HCal_clus_e_cut" ){lookup_val = simc_cuts[i][6];}
			if( param == "Ep" ){lookup_val = simc_cuts[i][7];}
			if( param == "Ep_sigma" ){lookup_val = simc_cuts[i][8];}
			if( param == "W2" ){lookup_val = simc_cuts[i][9];}
			if( param == "W2_sigma" ){lookup_val = simc_cuts[i][10];}
			if( param == "W" ){lookup_val = simc_cuts[i][11];}
			if( param == "W_sigma" ){lookup_val = simc_cuts[i][12];}

		}
	}


	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;
}
double lookup_parsed_dx_parlimits( int kine, TString run_targ, int sbsfield, TString param){
	double lookup_val = -1;

	int target_int;

	if( run_targ == "LH2" ){
		target_int = 0;
	}
	if( run_targ == "LD2" ){
		target_int = 1;
	}

	//Cuts vector:
// 	[0]target_int, [1]kine, [2]sbs_field [3]PS_clus_e_cut, [4]SH_PS_mean, [5]SH_PS_sigma, [6]HCal_clus_e_cut, [7]Ep, [8]Ep_sigma, [9]W2, [10]W2_sigma, [11]W, [12]W_sigma, 
	vector<vector<double>> dx_parlimits = {
		{1, 8, 70, }
		
	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-run_02

	for( size_t i = 0; i < dx_parlimits.size(); i++){
		if( dx_parlimits[i][0] == target_int && dx_parlimits[i][1] == kine && dx_parlimits[i][2] == sbsfield ){
			if( param == "par1min" ){lookup_val = dx_parlimits[i][3];}
			if( param == "par1max" ){lookup_val = dx_parlimits[i][4];}
			if( param == "par2min" ){lookup_val = dx_parlimits[i][5];}
			if( param == "par2max" ){lookup_val = dx_parlimits[i][6];}
			if( param == "par4min" ){lookup_val = dx_parlimits[i][7];}
			if( param == "par4max" ){lookup_val = dx_parlimits[i][8];}
			if( param == "par5min" ){lookup_val = dx_parlimits[i][9];}
			if( param == "par5max" ){lookup_val = dx_parlimits[i][10];}
			if( param == "par6min" ){lookup_val = dx_parlimits[i][11];}
			if( param == "par6max" ){lookup_val = dx_parlimits[i][12];}
			if( param == "par7min" ){lookup_val = dx_parlimits[i][13];}
			if( param == "par7max" ){lookup_val = dx_parlimits[i][14];}

		}
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}

	return lookup_val;	
}

double lookup_parsed_dy_parlimits( int kine, TString run_targ, int sbsfield, TString param){
	double lookup_val = -1;

	int target_int;

	if( run_targ == "LH2" ){
		target_int = 0;
	}
	if( run_targ == "LD2" ){
		target_int = 1;
	}

	//Cuts vector:
// 	[0]target_int, [1]kine, [2]sbs_field [3]PS_clus_e_cut, [4]SH_PS_mean, [5]SH_PS_sigma, [6]HCal_clus_e_cut, [7]Ep, [8]Ep_sigma, [9]W2, [10]W2_sigma, [11]W, [12]W_sigma, 
	vector<vector<double>> dy_parlimits = {
		{1, 8, 70, }
		
	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-run_02

	for( size_t i = 0; i < dy_parlimits.size(); i++){
		if( dy_parlimits[i][0] == target_int && dy_parlimits[i][1] == kine && dy_parlimits[i][2] == sbsfield ){
			if( param == "par1min" ){lookup_val = dy_parlimits[i][3];}
			if( param == "par1max" ){lookup_val = dy_parlimits[i][4];}
			if( param == "par2min" ){lookup_val = dy_parlimits[i][5];}
			if( param == "par2max" ){lookup_val = dy_parlimits[i][6];}
			if( param == "par4min" ){lookup_val = dy_parlimits[i][7];}
			if( param == "par4max" ){lookup_val = dy_parlimits[i][8];}
			if( param == "par5min" ){lookup_val = dy_parlimits[i][9];}
			if( param == "par5max" ){lookup_val = dy_parlimits[i][10];}
			if( param == "par6min" ){lookup_val = dy_parlimits[i][11];}
			if( param == "par6max" ){lookup_val = dy_parlimits[i][12];}
			if( param == "par7min" ){lookup_val = dy_parlimits[i][13];}
			if( param == "par7max" ){lookup_val = dy_parlimits[i][14];}

		}
	}

	if( lookup_val == -1 ){ 
		std::cout << "Lookup for " << param.Data() << " may not exist. Returned value of -1. " << std::endl;
		cout << "Sleeping to allow for a ctrl-c catch..." << endl;
		sleep(20);
		cout << "Exiting...." << endl;
		exit(1);	
	}
	
	return lookup_val;	
}

double lookup_simc_cut_byPass(int pass, TString run_target, int kine, int sbsfieldscale, TString selection ){
	double return_val = 0.0;
	cout << "lookup_simc_cut_byPass: " << pass << ", " << run_target.Data() << ", " << kine << ", " << sbsfieldscale << ", " << selection.Data() << endl;
	if( pass == 1 ){
		return_val = lookup_simc_cut( kine, sbsfieldscale, run_target, selection );
	}
	if( pass == 2 ){
		return_val = lookup_pass2_simc_cut( kine, sbsfieldscale, run_target, selection );
	}
	return return_val;
}

double lookup_simc_dxdy_by_kine_and_mag_byPass(int pass, TString run_target, int kine, int sbsfieldscale, TString selection ){
	double return_val = 0.0;
	if( pass == 1 ){
		return_val = lookup_simc_dxdy( kine, sbsfieldscale, run_target, selection );
	}
	if( pass == 2 ){
		return_val = lookup_pass2_simc_dxdy( kine, sbsfieldscale, run_target, selection );
	}
	return return_val;
}

double lookup_simc_dxdy_by_kine_and_mag_with_portField_byPass(int pass, TString run_target, int kine, int sbsfieldscale, TString portField, TString selection ){
	double return_val = 0.0;
	// cout << "lookup_simc_dxdy_by_kine_and_mag_with_portField_byPass: " << pass << ", " << run_target.Data() << ", " << kine << ", " << sbsfieldscale << ", portField: " << portField.Data() << ", selection: " << selection.Data() << endl;

	if( pass == 1 ){
		return_val = lookup_simc_dxdy_with_portField( kine, sbsfieldscale, run_target, portField, selection );
	}
	if( pass == 2 ){
		return_val = lookup_simc_dxdy_with_portField( kine, sbsfieldscale, run_target, portField, selection );
	}
	return return_val;
}

#endif
