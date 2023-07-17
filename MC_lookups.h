#ifndef MC_LOOKUPS_H
#define MC_LOOKUPS_H

double lookup_MC_dxdy(int kine, int sbsfieldscale, TString param){
	double lookup_val = -1;

	//dxdy vector:
// 	kine, sbsfieldscale, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma

	vector<vector<double>> MC_dxdy = {
		{4, 0, },
		{4, 30, },
		{4, 50, -1.1, 0.15, 0.0, 0.11, 0.0, 0.20}, 
	};

	for( size_t i = 0; i < MC_dxdy.size(); i++){
		if( MC_dxdy[i][0] == kine && MC_dxdy[i][1] == param ){
			if( param == "dx_p" ){lookup_val = MC_dxdy[i][2];}
			if( param == "dx_p_sigma" ){lookup_val = MC_dxdy[i][3];}
			if( param == "dx_n" ){lookup_val = MC_dxdy[i][4];}
			if( param == "dx_n_sigma" ){lookup_val = MC_dxdy[i][5];}
			if( param == "dy" ){lookup_val = MC_dxdy[i][6];}
			if( param == "dy_sigma" ){lookup_val = MC_dxdy[i][7];}

		}
	}
	return lookup_val;

}

double lookup_MC_cut( int kine, int sbsfieldscale, TString param ){
	double lookup_val = -1;

	//dxdy vector:
// 	W, W_sigma, W2, W2_sigma

	vector<vector<double>> MC_cut = {
		{4, 0, },
		{4, 30, },
		{4, 50, 0.9772, 0.15, 0.9605, 0.15}, 
	};

	for( size_t i = 0; i < MC_cut.size(); i++){
		if( MC_cut[i][0] == kine && MC_cut[i][1] == param ){
			if( param == "dx_p" ){lookup_val = MC_cut[i][2];}
			if( param == "dx_p_sigma" ){lookup_val = MC_cut[i][3];}
			if( param == "dx_n" ){lookup_val = MC_cut[i][4];}
			if( param == "dx_n_sigma" ){lookup_val = MC_cut[i][5];}
			if( param == "dy" ){lookup_val = MC_cut[i][6];}
			if( param == "dy_sigma" ){lookup_val = MC_cut[i][7];}

		}
	}
	return lookup_val;

}


#endif
