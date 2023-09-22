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
		{8, 70, -8.73151e-01, 1.74421e-01, -2.41303e-02, 1.52362e-01, 6.31940e-02, 2.26323e-01},
		{9, 70, -8.55327e-01, 1.49493e-01, -5.06166e-02, 2.00000e-01, 6.87480e-02, 1.79564e-01}
		// {8, 70, -0.849403, 0.14, -0.00315522, 0.14, -0.0463545, 0.23}
		// {8, 70, -0.852848, 0.170124, 2.96250e-04, 0.160995, -2.56391e-02, 1.95379e-01}

	};

	if( run_targ == "LD2" ){
		for( size_t i = 0; i < MC_LD2_dxdy.size(); i++){
			if( MC_LD2_dxdy[i][0] == kine && MC_LD2_dxdy[i][1] == param ){
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
			if( MC_LD2_dxdy[i][0] == kine && MC_LH2_dxdy[i][1] == param ){
				if( param == "dx_p" ){lookup_val = MC_LH2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = MC_LH2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = -1;}
				if( param == "dx_n_sigma" ){lookup_val = -1;}
				if( param == "dy" ){lookup_val = MC_LH2_dxdy[i][4];}
				if( param == "dy_sigma" ){lookup_val = MC_LH2_dxdy[i][5];}

			}
		}	
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
			if( simc_LD2_dxdy[i][0] == kine && simc_LD2_dxdy[i][1] == param ){
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
			if( simc_LD2_dxdy[i][0] == kine && simc_LH2_dxdy[i][1] == param ){
				if( param == "dx_p" ){lookup_val = simc_LH2_dxdy[i][2];}
				if( param == "dx_p_sigma" ){lookup_val = simc_LH2_dxdy[i][3];}
				if( param == "dx_n" ){lookup_val = -1;}
				if( param == "dx_n_sigma" ){lookup_val = -1;}
				if( param == "dy" ){lookup_val = simc_LH2_dxdy[i][4];}
				if( param == "dy_sigma" ){lookup_val = simc_LH2_dxdy[i][5];}

			}
		}	
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
		{1,  8,    70, 0.150, 3.25473, 0.419763, 0.0200, 0.956448, 0.0780319, 0.909982, 0.306327, 0.961529, 0.185036},
		{1,  9,    70, 0.150, 1.56300, 0.291986, 0.0200, 1.037640, 0.1315590, 9.563250, 0.280737, 1.003870, 0.164814}
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
	return lookup_val;	
}

#endif
