#ifndef BEAM_VARIABLES_H
#define BEAM_VARIABLES_H
#include <algorithm>

TString target;
double sbsdist = 2.25;

double lookup_beam_energy( int runnum ){
	double beam_energy = 0;

	if( runnum >= 11436 && runnum <= 11616 ){
		beam_energy = 3.728;
	}
	else if( runnum >= 11989 && runnum <= 12073 ){
		beam_energy = 7.906;
	}
	else if( runnum >= 12313 && runnum <= 13063 ){
		beam_energy = 9.91;
	}
	else if( runnum >= 13239 && runnum <= 13407 ){
		beam_energy = 5.965;
	}
	else if( runnum >= 13444 && runnum <= 13620 ){
		beam_energy = 5.965;
	}
	else if( runnum >= 13656 && runnum <= 13799 ){
		beam_energy = 4.013;
	}
	else{ beam_energy = -1.111; }

	return beam_energy;
}

int lookup_kinematic( int runnum ){
	int kinematic = 0;

	if( runnum >= 11436 && runnum <= 11616 ){
		kinematic = 4;
	}
	else if( runnum >= 11989 && runnum <= 12073 ){
		kinematic = 7;
	}
	else if( runnum >= 12313 && runnum <= 13063 ){
		kinematic = 11;
	}
	else if( runnum >= 13239 && runnum <= 13407 ){
		kinematic = 14;
	}
	else if( runnum >= 13444 && runnum <= 13620 ){
		kinematic = 8;
	}
	else if( runnum >= 13656 && runnum <= 13799 ){
		kinematic = 9;
	}
	else{ kinematic = -1; }

	return kinematic;
}

//Angle is in degrees
double lookup_BB_angle( int runnum ){
	double BB_angle_lookup;
	if( lookup_kinematic(runnum) == 1 ){
		BB_angle_lookup = 51.0;
	}
	if( lookup_kinematic(runnum) == 4 ){
		BB_angle_lookup = 36.0;
	}
	if( lookup_kinematic(runnum) == 7 ){
		BB_angle_lookup = 40.0;
	}
	if( lookup_kinematic(runnum) == 8 ){
		BB_angle_lookup = 26.5;
	}
	if( lookup_kinematic(runnum) == 9 ){
		BB_angle_lookup = 49.0;
	}
	if( lookup_kinematic(runnum) == 11 ){
		BB_angle_lookup = 42.0;
	}
	if( lookup_kinematic(runnum) == 14 ){
		BB_angle_lookup = 46.5;
	}
	return BB_angle_lookup;
}

double lookup_SBS_angle( int runnum ){
	double SBS_angle_lookup;
	if( lookup_kinematic(runnum) == 1 ){
		SBS_angle_lookup = 33.5;
	}
	if( lookup_kinematic(runnum) == 4 ){
		SBS_angle_lookup = 31.9;
	}
	if( lookup_kinematic(runnum) == 7 ){
		SBS_angle_lookup = 16.1;
	}
	if( lookup_kinematic(runnum) == 8 ){
		SBS_angle_lookup = 29.9;
	}
	if( lookup_kinematic(runnum) == 9 ){
		SBS_angle_lookup = 22.5;
	}
	if( lookup_kinematic(runnum) == 11 ){
		SBS_angle_lookup = 13.3;
	}
	if( lookup_kinematic(runnum) == 14 ){
		SBS_angle_lookup = 17.3;
	}
	return SBS_angle_lookup;
}

double lookup_HCal_angle( int runnum ){
	double HCal_angle_lookup;
	if( lookup_kinematic(runnum) == 1 ){
		HCal_angle_lookup = 33.5;
	}
	if( lookup_kinematic(runnum) == 4 ){
		HCal_angle_lookup = 31.9;
	}
	if( lookup_kinematic(runnum) == 7 ){
		HCal_angle_lookup = 16.1;
	}
	if( lookup_kinematic(runnum) == 8 ){
		HCal_angle_lookup = 29.4;
	}
	if( lookup_kinematic(runnum) == 9 ){
		HCal_angle_lookup = 22.0;
	}
	if( lookup_kinematic(runnum) == 11 ){
		HCal_angle_lookup = 13.3;
	}
	if( lookup_kinematic(runnum) == 14 ){
		HCal_angle_lookup = 17.3;
	}
	return HCal_angle_lookup;
}

double lookup_BB_dist( int runnum ){
	double BB_dist_lookup;
	if( lookup_kinematic(runnum) == 1 ){
		BB_dist_lookup = 1.85;
	}
	if( lookup_kinematic(runnum) == 4 ){
		BB_dist_lookup = 1.80;
	}
	if( lookup_kinematic(runnum) == 7 ){
		BB_dist_lookup = 1.85;
	}
	if( lookup_kinematic(runnum) == 8 ){
		BB_dist_lookup = 2.00;
	}
	if( lookup_kinematic(runnum) == 9 ){
		BB_dist_lookup = 1.55;
	}
	if( lookup_kinematic(runnum) == 11 ){
		BB_dist_lookup = 1.55;
	}
	if( lookup_kinematic(runnum) == 14 ){
		BB_dist_lookup = 1.85;
	}
	return BB_dist_lookup;
}

double lookup_HCal_dist( int runnum ){
	double HCal_dist_lookup;
	if( lookup_kinematic(runnum) == 1 ){
		HCal_dist_lookup = 13.5;
	}
	if( lookup_kinematic(runnum) == 4 ){
		HCal_dist_lookup = 11.0;
	}
	if( lookup_kinematic(runnum) == 7 ){
		HCal_dist_lookup = 14.0;
	}
	if( lookup_kinematic(runnum) == 8 ){
		HCal_dist_lookup = 11.0;
	}
	if( lookup_kinematic(runnum) == 9 ){
		HCal_dist_lookup = 11.0;
	}
	if( lookup_kinematic(runnum) == 11 ){
		HCal_dist_lookup = 14.5;
	}
	if( lookup_kinematic(runnum) == 14 ){
		HCal_dist_lookup = 14.0;
	}
	return HCal_dist_lookup;
}

TString lookup_target( int runnum ){

	vector<TString> targets = {"LH2", "LD2"};
	TString target_lookup = "XXX";
	int target_found = 0;

	vector<int> runs_LH2 = {11436, 11500, 11547, 11548, 11573, 11587, 11588, 11589, 11590, 11592, 11616, 11989, 11990, 11991, 11992, 11993, 11994, 12000, 12008, 12022, 12035, 12049, 12052, 12053, 12058, 12064, 12072, 12313, 12320, 12335, 12336, 12340, 12345, 12355, 12358, 12363, 12367, 12368, 12369, 12370, 12380, 12382, 12400, 12401, 12414, 12415, 12427, 12428, 12429, 12471, 12472, 12491, 12496, 12497, 12524, 12525, 12526, 12548, 12561, 12563, 12575, 12576, 12587, 12588, 12625, 12626, 12657, 12658, 12675, 12676, 12677, 12961, 12692, 12704, 12705, 12716, 12730, 12731, 12747, 12758, 12765, 12776, 12777, 12809, 12810, 12812, 12822, 12823, 12895, 12910, 12923, 12930, 12931, 12959, 12960, 12973, 12974, 13027, 13028, 13041, 13042, 13043, 13056, 13057, 13239, 13240, 13241, 13242, 13243, 13244, 13312, 13313, 13320, 13321, 13345, 13346, 13348, 13349, 13351, 13352, 13375, 13376, 13377, 13378, 13379, 13396, 13397, 13405, 13444, 13450, 13451, 13452, 13459, 13460, 13461, 13463, 13464, 13465, 13466, 13482, 13483, 13484, 13485, 13486, 13487, 13488, 13489, 13490, 13537, 13538, 13539, 13540, 13541, 13542, 13543, 13573, 13574, 13575, 13576, 13577, 13578, 13580, 13656, 13657, 13663, 13676, 13683, 13696, 13697, 13719, 13720, 13730, 13747, 13768, 13769 ,13795, 13796};

	vector<int> runs_LD2 = {11449, 11451, 11452, 11456, 11493, 11494, 11495, 11496, 11551, 11554, 11562, 11563, 11564, 11565, 11568, 11570, 11571, 11579, 11580, 11581, 11582, 11583, 11586, 11593, 11594, 11595, 11996, 11997, 11998, 11999, 12001, 12002, 12004, 12006, 12013, 12014, 12017, 12019, 12021, 12029, 12030, 12038, 12039, 12040, 12041, 12042, 12043, 12044, 12045, 12046, 12047, 12048, 12050, 12051, 12055, 12056, 12057, 12059, 12060, 12062, 12063, 12065, 12066, 12068, 12069, 12070, 12071, 12073, 12314, 12315, 12316, 12318, 12319, 12321, 12322, 12323, 12333, 12337, 12338, 12339, 12341, 12342, 12343, 12344, 12346, 12347, 12348, 12349, 12356, 12357, 12359, 12360, 12362, 12364, 12366, 12371, 12372, 12373, 12374, 12375, 12376, 12377, 12378, 12379, 12383, 12384, 12385, 12387, 12388, 12394, 12395, 12396, 12403, 12404, 12405, 12406, 12407, 12408,12409, 12410, 12411, 12412, 12416, 12423, 12424, 12425, 12473, 12474, 12478, 12479, 12480, 12481, 12482, 12485, 12486, 12498, 12499, 12500, 12501, 12502, 12512, 12513, 12514, 12515, 12516, 12517, 12521, 12523, 12527, 12532, 12534, 12535, 12536, 12537, 12538, 12539, 12541, 12542, 12543, 12544, 12547, 12549, 12550, 12551, 12552, 12553, 12554, 12556, 12558, 12559, 12565, 12566, 12567, 12568, 12569, 12570, 12571, 12572, 12574, 12577, 12580, 12581, 12582, 12583, 12584, 12585, 12586, 12589, 12590, 12592, 12593, 12597, 12598, 12599, 12614, 12618, 12619, 12620, 12623, 12627, 12628, 12638, 12639, 12641, 12642, 12643, 12644, 12645, 12647, 12653, 12654, 12655, 12656, 12659, 12660, 12661, 12662, 12664, 12666, 12671, 12673, 12674, 12678, 12679, 12680, 12681, 12682, 12683, 12684, 12686, 12687, 12688, 12693, 12694, 12695, 12697, 12699, 12700, 12701, 12702, 12706, 12707, 12708, 12709, 12710, 12711, 12713, 12714, 12715, 12721, 12722, 12726, 12728, 12729, 12732, 12733, 12734, 12735, 12737, 12738, 12739, 12740, 12742, 12743, 12767, 12768, 12769, 12770, 12771, 12772, 12773, 12778, 12779, 12781, 12782, 12796, 12797, 12798, 12799, 12801, 12802, 12813, 12814, 12815, 12816, 12817, 12819, 12820, 12824, 12825, 12826, 12827, 12828, 12829, 12830, 12894, 12896, 12897, 12902, 12903, 12904, 12906, 12907, 12908, 12909, 12911, 12912, 12913, 12914, 12915, 12916, 12917, 12918, 12919, 12920, 12921, 12932, 12933, 12934, 12935, 12936, 12939, 12940, 12946, 12957, 12961, 12962, 12963, 12965, 12968, 12969, 12970, 12971, 12972, 12975, 12976, 12977, 12979, 12980, 12982, 12983, 12984, 12986, 12987, 12988, 12989, 13000, 13029, 13030, 13034, 13035, 13036, 13038, 13039, 13040, 13044, 13045, 13048, 13049, 13050, 13051, 13052, 13053, 13054, 13055, 13058, 13059, 13060, 13061, 13062, 13063, 13305, 13306, 13307, 13308, 13309, 13314, 13315, 13316, 13317, 13318, 13319, 13322, 13323, 13324, 13325, 13342, 13344, 13353, 13357, 13358, 13359, 13360, 13361, 13362, 13363, 13364, 13368, 13369, 13370, 13371, 13372, 13373, 13381, 13382, 13384, 13385, 13387, 13388, 13389, 13390, 13391, 13392, 13393, 13394, 13395, 13398, 13399, 13400, 13402, 13403, 13406, 13407, 13453, 13454, 13455, 13468, 13470, 13472, 13473, 13474, 13475, 13476, 13477, 13478, 13479, 13491, 13492, 13493, 13494, 13495, 13496, 13497, 13502, 13503, 13504, 13505, 13544, 13545, 13546, 13547, 13548, 13549, 13550, 13551, 13552, 13554, 13556, 13557, 13558, 13559, 13560, 13561, 13562, 13563, 13564, 13565, 13566, 13567,13568, 13569, 13570, 13571, 13581, 13582, 13583, 13584, 13585, 13586, 13587, 13588, 13589, 13590, 13591, 13592, 13593, 13596, 13597, 13608, 13609, 13610, 13612, 13613, 13614, 13615, 13616, 13617, 13618, 13619, 13620, 13660, 13661, 13662, 13664, 13665, 13666, 13677, 13678, 13679, 13680, 13681, 13682, 13684, 13685, 13686, 13687, 13688, 13689, 13694, 13695, 13698, 13699, 13700, 13710, 13711, 13712, 13714, 13715, 13716, 13717, 13721, 13723, 13724, 13727, 13728, 13729, 13731, 13732, 13734, 13736, 13737, 13746, 13748, 13749, 13753, 13754, 13755,13756, 13757, 13758, 13760, 13761, 13764, 13765, 13766, 13767, 13770, 13771, 13773, 13775, 13776, 13777, 13778, 13779, 13793, 13797, 13798, 13799};


	for(size_t i = 0; i < runs_LH2.size(); i++){
		if( runs_LH2[i] == runnum ){
			target_lookup = targets[0];
			target_found = 1;
			break;

		}
	}

	for(size_t i = 0; i < runs_LD2.size(); i++){
		if( runs_LD2[i] == runnum ){
			target_lookup = targets[1];
			target_found = 1;
			break;
		}
	}

	return target_lookup;
}

double lookup_cut(int runnum, TString param){
	double lookup_val = -1;
	//Cuts vector:
// 	[0] runnum, [1] Ep, [2]Ep_sigma, [3]PS_clus_e_cut, [4]SH_PS_clus_e_cut, [5]SH_PS_sigma, [6]HCal_clus_e_cut, [7]W2, [8]W2_sigma, [9]W, [10]W_sigma, 
	vector<vector<double>> cuts = {
		{11449, 0.974283, 0.0661193, 0.150, 1.93, 0.183409, 0.03, 9.81425e-01, 2.27171e-01, 1.005, 1.41647e-01},
		{11493, 9.64468e-01, 8.41696e-02, 0.150, 1.49721, 2.43353e-01, .0075, 9.51695e-01, 1.88384e-01 },
		{11494, 0.964638, 0.258728, 0.150, 1.7, 0.02, 0.959381, .0723608, 0.996746, 0.129869 },
		{13479, 9.61859e-01, 7.52286e-02, 0.25, 2.80079, 4.23989e-01, 3.28841e-02, 9.17414e-01, 4.85251e-01},
		{13544, 9.57190e-01, 7.38953e-02, 0.25, 2.74275, 3.96009e-01, 1.84462e-02, 1.02312, 2.68668e-01},
		{13566, 1.02859, 0.0765564, 0.15, 2.95883, 0.417844, 0.015, 0.997588, 0.29556, 1, 0.16742 },
		{13585, 0.956980, 0.068656, 0.25, 2.76844, 3.96853e-01, .025, 9.86933e-01, 4.72408e-01, 0.955, 0.22863400 }

	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-02
	for( size_t i = 0; i < cuts.size(); i++){
		if( cuts[i][0] == runnum ){
			if( param == "Ep" ){lookup_val = cuts[i][1];}
			if( param == "Ep_sigma" ){lookup_val = cuts[i][2];}
			if( param == "PS_clus_e_cut" ){lookup_val = cuts[i][3];}
			if( param == "SH_PS_clus_e_cut" ){lookup_val = cuts[i][4];}
			if( param == "SH_PS_sigma" ){lookup_val = cuts[i][5];}
			if( param == "HCal_clus_e_cut" ){lookup_val = cuts[i][6];}
			if( param == "W2" ){lookup_val = cuts[i][7];}
			if( param == "W2_sigma" ){lookup_val = cuts[i][8];}
			if( param == "W" ){lookup_val = cuts[i][9];}
			if( param == "W_sigma" ){lookup_val = cuts[i][10];}

		}
	}
	return lookup_val;
}

double lookup_calibrate_fit(TString target, int kine, int sbsfield, TString param){
	double lookup_val = -1;
	int target_int;
	if( target == "LH2"){
		target_int = 0;
	}
	if( target == "LD2" ){
		target_int = 1;
	}
	//Cuts vector:
// [0]target_int, [1]kine, [2]sbsfield, [3]PS_min, [4]PS_max, [5]PS_sigma_min [6]PS_sigma_max, [7] SH_PS_min, [8]SH_PS_max, [9]SH_PS_sigma_min, [10]SH_PS_sigma_max
// [11] Ep_min, [12]Ep_max, [13]Ep_sigma_min, [14] Ep_sigma_max, [15]W2_min, [16]W2_max, 
//[17]  W2_sigma_min, [18] W2_sigma_max, [19]W_min, [20]W_max, [21]W_sigma_min, [22]W_sigma_max
	vector<vector<double>> cuts = {
		{1, 8, 0, 0.8, 1.2, 0.01, 0.4, 1.07, 1.13, 0.01, 0.8, 1.9, 3.7, 0.01, 0.5, 0.77, 0.84, 0.01, 0.27, 0.87, 0.94, 0.01, 0.4},
		{1, 8, 50, 0.8, 1.2, 0.01, 0.4, 1.07, 1.13, 0.01, 0.8, 1.9, 3.7, 0.01, 0.5, 0.87, 0.925, 0.01, 0.26, 0.95, 1.05, 0.01, 0.35},
		{1, 8, 70, 0.9, 1.25, 0.01, 0.4, 0.95, 1.2, 0.01, 0.6, 1.9, 3.7, 0.01, 0.5, 0.85, .949, 0.01, 0.3, 0.95, 1.05, 0.01, 0.35},
		{1, 9, 70, 0.975, 0.99, 0.01, 0.4, 0.55, 0.625, 0.01, 0.6, 1.56, 1.66, 0.01, 0.5 }

	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-02
	for( size_t i = 0; i < cuts.size(); i++){
		if( cuts[i][0] == target_int && cuts[i][1] == kine && cuts[i][2] == sbsfield ){
			if( param == "Ep_min" ){lookup_val = cuts[i][3];}
			if( param == "Ep_max" ){lookup_val = cuts[i][4];}
			if( param == "Ep_min_sigma" ){lookup_val = cuts[i][5];}
			if( param == "Ep_max_sigma" ){lookup_val = cuts[i][6];}
			if( param == "PS_min" ){lookup_val = cuts[i][7];}
			if( param == "PS_max" ){lookup_val = cuts[i][8];}
			if( param == "PS_min_sigma" ){lookup_val = cuts[i][9];}
			if( param == "PS_max_sigma" ){lookup_val = cuts[i][10];}
			if( param == "SH_PS_min" ){lookup_val = cuts[i][11];}
			if( param == "SH_PS_max" ){lookup_val = cuts[i][12];}
			if( param == "SH_PS_min_sigma" ){lookup_val = cuts[i][13];}
			if( param == "SH_PS_max_sigma" ){lookup_val = cuts[i][14];}
			if( param == "W2_min" ){lookup_val = cuts[i][15];}
			if( param == "W2_max" ){lookup_val = cuts[i][16];}
			if( param == "W2_min_sigma" ){lookup_val = cuts[i][17];}
			if( param == "W2_max_sigma" ){lookup_val = cuts[i][18];}
			if( param == "W_min" ){lookup_val = cuts[i][19];}
			if( param == "W_max" ){lookup_val = cuts[i][20];}
			if( param == "W_min_sigma" ){lookup_val = cuts[i][21];}
			if( param == "W_max_sigma" ){lookup_val = cuts[i][22];}

		}
	}

	vector<vector<double>> default_cuts = {
		{-1, -1, -1, 0.8, 1.2, 0.01, 0.4, 0.7, 1.2, 0.01, 0.8, 1.9, 3.7, 0.01, 0.5, 0.85, 0.98, 0.01, 0.4, 0.80, 0.985, 0.01, 0.4}

	};
	if( lookup_val == -1){
		for( size_t i = 0; i < default_cuts.size(); i++){
			if( param == "Ep_min" ){lookup_val = default_cuts[i][3];}
			if( param == "Ep_max" ){lookup_val = default_cuts[i][4];}
			if( param == "Ep_min_sigma" ){lookup_val = default_cuts[i][5];}
			if( param == "Ep_max_sigma" ){lookup_val = default_cuts[i][6];}
			if( param == "PS_min" ){lookup_val = default_cuts[i][7];}
			if( param == "PS_max" ){lookup_val = default_cuts[i][8];}
			if( param == "PS_min_sigma" ){lookup_val = default_cuts[i][9];}
			if( param == "PS_max_sigma" ){lookup_val = default_cuts[i][10];}
			if( param == "SH_PS_min" ){lookup_val = default_cuts[i][11];}
			if( param == "SH_PS_max" ){lookup_val = default_cuts[i][12];}
			if( param == "SH_PS_min_sigma" ){lookup_val = default_cuts[i][13];}
			if( param == "SH_PS_max_sigma" ){lookup_val = default_cuts[i][14];}
			if( param == "W2_min" ){lookup_val = default_cuts[i][15];}
			if( param == "W2_min_sigma" ){lookup_val = default_cuts[i][16];}
			if( param == "W2_max" ){lookup_val = default_cuts[i][17];}
			if( param == "W2_max_sigma" ){lookup_val = default_cuts[i][18];}
			if( param == "W_min" ){lookup_val = default_cuts[i][19];}
			if( param == "W_min_sigma" ){lookup_val = default_cuts[i][20];}
			if( param == "W_max" ){lookup_val = default_cuts[i][21];}
			if( param == "W_max_sigma" ){lookup_val = default_cuts[i][22];}
		}	
	}
	
	return lookup_val;
}

double lookup_parsed_fit(TString target, int kine, TString param){
	double lookup_val = -1;
	int target_int;
	if( target == "LH2"){
		target_int = 0;
	}
	if( target == "LD2" ){
		target_int = 1;
	}
	//Cuts vector:
// [0]target_int, [1]kine, [2]sbsfield, [3]PS_min, [4]PS_min_sigma, [5]PS_max [6]PS_max_sigma, [7] SH_PS_min, [8]SH_PS_min_sigma, [9]SH_PS_max, [10]SH_PS_max_sigma
// [11] Ep_min, [12] Ep_min_sigma, [13]Ep_max, [14] Ep_max_sigma, [15]W2_min, [16] W2_min_sigma, 
//[17] W2_max, [18] W2_max_sigma, [19]W_min, [20]W_min_sigma, [21]W_max, [22]W_max_sigma
	vector<vector<double>> parsed_fit_cuts = {
		{1, 9, -1, 0.95, 1.02, 0.01, 0.12, 0.41, 0.51, 0.01, 0.8, 1.45, 1.55, 0.01, 0.5, 0.85, 0.97, 0.01, 0.4, 0.80, 0.995, 0.01, 0.4}

	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-02
	for( size_t i = 0; i < parsed_fit_cuts.size(); i++){
		if( parsed_fit_cuts[i][0] == target_int && parsed_fit_cuts[i][1] == kine ){
			if( param == "Ep_min" ){lookup_val = parsed_fit_cuts[i][3];}
			if( param == "Ep_max" ){lookup_val = parsed_fit_cuts[i][4];}
			if( param == "Ep_min_sigma" ){lookup_val = parsed_fit_cuts[i][5];}
			if( param == "Ep_max_sigma" ){lookup_val = parsed_fit_cuts[i][6];}
			if( param == "PS_min" ){lookup_val = parsed_fit_cuts[i][7];}
			if( param == "PS_max" ){lookup_val = parsed_fit_cuts[i][8];}
			if( param == "PS_min_sigma" ){lookup_val = parsed_fit_cuts[i][9];}
			if( param == "PS_max_sigma" ){lookup_val = parsed_fit_cuts[i][10];}
			if( param == "SH_PS_min" ){lookup_val = parsed_fit_cuts[i][11];}
			if( param == "SH_PS_max" ){lookup_val = parsed_fit_cuts[i][12];}
			if( param == "SH_PS_min_sigma" ){lookup_val = parsed_fit_cuts[i][13];}
			if( param == "SH_PS_max_sigma" ){lookup_val = parsed_fit_cuts[i][14];}
			if( param == "W2_min" ){lookup_val = parsed_fit_cuts[i][15];}
			if( param == "W2_min_sigma" ){lookup_val = parsed_fit_cuts[i][16];}
			if( param == "W2_max" ){lookup_val = parsed_fit_cuts[i][17];}
			if( param == "W2_max_sigma" ){lookup_val = parsed_fit_cuts[i][18];}
			if( param == "W_min" ){lookup_val = parsed_fit_cuts[i][19];}
			if( param == "W_min_sigma" ){lookup_val = parsed_fit_cuts[i][20];}
			if( param == "W_max" ){lookup_val = parsed_fit_cuts[i][21];}
			if( param == "W_max_sigma" ){lookup_val = parsed_fit_cuts[i][22];}

		}
	}

	vector<vector<double>> default_parsed_fit_cuts = {
		{-1, -1, -1, 0.8, 1.2, 0.01, 0.4, 0.8, 1.2, 0.01, 0.8, 1.35, 1.6, 0.01, 0.5, 0.85, 0.98, 0.01, 0.4, 0.80, 0.985, 0.01, 0.4},

	};
	if( lookup_val == -1){
		for( size_t i = 0; i < default_parsed_fit_cuts.size(); i++){
			if( param == "Ep_min" ){lookup_val = default_parsed_fit_cuts[i][3];}
			if( param == "Ep_max" ){lookup_val = default_parsed_fit_cuts[i][4];}
			if( param == "Ep_min_sigma" ){lookup_val = default_parsed_fit_cuts[i][5];}
			if( param == "Ep_max_sigma" ){lookup_val = default_parsed_fit_cuts[i][6];}
			if( param == "PS_min" ){lookup_val = default_parsed_fit_cuts[i][7];}
			if( param == "PS_max" ){lookup_val = default_parsed_fit_cuts[i][8];}
			if( param == "PS_min_sigma" ){lookup_val = default_parsed_fit_cuts[i][9];}
			if( param == "PS_max_sigma" ){lookup_val = default_parsed_fit_cuts[i][10];}
			if( param == "SH_PS_min" ){lookup_val = default_parsed_fit_cuts[i][11];}
			if( param == "SH_PS_max" ){lookup_val = default_parsed_fit_cuts[i][12];}
			if( param == "SH_PS_min_sigma" ){lookup_val = default_parsed_fit_cuts[i][13];}
			if( param == "SH_PS_max_sigma" ){lookup_val = default_parsed_fit_cuts[i][14];}
			if( param == "W2_min" ){lookup_val = default_parsed_fit_cuts[i][15];}
			if( param == "W2_min_sigma" ){lookup_val = default_parsed_fit_cuts[i][16];}
			if( param == "W2_max" ){lookup_val = default_parsed_fit_cuts[i][17];}
			if( param == "W2_max_sigma" ){lookup_val = default_parsed_fit_cuts[i][18];}
			if( param == "W_min" ){lookup_val = default_parsed_fit_cuts[i][19];}
			if( param == "W_min_sigma" ){lookup_val = default_parsed_fit_cuts[i][20];}
			if( param == "W_max" ){lookup_val = default_parsed_fit_cuts[i][21];}
			if( param == "W_max_sigma" ){lookup_val = default_parsed_fit_cuts[i][22];}
		}	
	}
	
	return lookup_val;
}

double lookup_parsed_cut(TString target, int kinematic, int mag, TString param){
	double lookup_val = -1;

	int target_int = -1;

	//Targets: LH2 == 0, LD2 == 1
	if( target == "LH2"){
		target_int = 0;
	}
	if( target == "LD2"){
		target_int = 1;
	}

	//Cuts vector:
// [0] run_target, [1]Kine, [2] sbsfield, [3]Ep_mean, [4]Ep_sigma [5]PS_min, [6]SH_PS_mean, [7]SH_PS_sigma, [8]HCal_min, [9]W2_mean, [10]W2_sigma, [11]W_mean, [12]W_sigma, 
//OLD:  runnum, PS_clus_e_cut, SH_PS_clus_e_cut, SH_PS_sigma, HCal_clus_e_cut, Ep, Ep_sigma, W2, W2_sigma, W, W_sigma, 

	vector<vector<double>> cuts = {
		{1, 4, 0, 0.956527, 0.064244, 0.15, 2.01496, 0.158437, 0.026, 0.911756, 0.219805, 0.966038, 0.110807 },
		{1, 4, 30, 0.958483, 0.0666202, 0.15, 2.00752, 0.157567, 0.026, 0.959671, 0.198965, 0.985, 0.103403 },
		{1, 4, 50, 0.965375, 0.0673284, 0.15, 2.01962, 0.159043, 0.026, 0.98, 0.182398, 0.985, 0.0913345 },
		{1, 7, 1.30, 0.6, 0, 0.2, 1.10, 0.04, 0.75, 1.3},
		{1, 8, 0, 0.966829, 0.0604167, 0.18, 3.4578, 0.267153, 0.028, 0.84, 0.27, 0.94, 0.150904 },
		{1, 8, 50, 0.963503, 0.0588214, 0.18, 3.4242, 0.262852, 0.028, 0.925, 0.26, 0.987409, 0.131511 },
		{1, 8, 70, 1.00182, 0.0695257, 0.18, 3.55784, 0.320457, 0.026, 0.949, 0.265948, 0.994452, 0.123653 },
		// {1, 9, 70,},
		{1, 11, },
		{1, 14, },
		{0, 4, },
		{0, 7, },
		{0, 11, },
		{0, 14, },
		{0, 8, },
		{0, 9, }

	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-02
	for( size_t i = 0; i < cuts.size(); i++){
		if( cuts[i][0] == target_int && cuts[i][1] == kinematic && cuts[i][2] == mag){
			if( param == "Ep_mean" ){lookup_val = cuts[i][3];}
			if( param == "Ep_sigma" ){lookup_val = cuts[i][4];}
			if( param == "Ep_min" ){lookup_val = cuts[i][3] - cuts[i][4];}
			if( param == "Ep_max" ){lookup_val = cuts[i][3] + cuts[i][4];}
			if( param == "PS_min" ){lookup_val = cuts[i][5];}
			if( param == "SH_PS_mean" ){lookup_val = cuts[i][6];}
			if( param == "SH_PS_sigma" ){lookup_val = cuts[i][7];}
			if( param == "SH_PS_min" ){lookup_val = cuts[i][6] - cuts[i][7];}
			if( param == "SH_PS_max" ){lookup_val = cuts[i][6] + cuts[i][7];}
			if( param == "HCal_min" ){lookup_val = cuts[i][8];}
			if( param == "W2_min" ){lookup_val = cuts[i][9] - cuts[i][10];}
			if( param == "W2_max" ){lookup_val = cuts[i][9] + cuts[i][11];}
			if( param == "W2_mean" ){lookup_val = cuts[i][9];}
			if( param == "W2_sigma" ){lookup_val = cuts[i][10];}
			if( param == "W_mean" ){lookup_val = cuts[i][11];}
			if( param == "W_sigma" ){lookup_val = cuts[i][12];}

		}
	}
	return lookup_val;
}

double lookup_pre_parsed_cut(TString target, int kinematic, TString param){
	double lookup_val = -1;

	int target_int = -1;

	//Targets: LH2 == 0, LD2 == 1
	if( target == "LH2"){
		target_int = 0;
	}
	if( target == "LD2"){
		target_int = 1;
	}

	//Cuts vector:
// run_target, Kine, Ep_min, Ep_max, PS_min, SH_PS_min, HCal_min, W2_min, W2_max, W_min, W_max
	vector<vector<double>> cuts = {
		{1, 4, 0.75, 1.2, 0.140, 1.25, 0.025, 0.6, 1.3, 9.51695e-01, 1.88384e-01, 0.955, 0.22863400},
		{1, 7, 0.75, 1.30, 0.2, 1.10, 0.04, 0.6, 1.3},
		{1, 8, 0.78, 1.22, 0.175, 2.20, 0.025, 0.6, 1.3},
		// {1, 9, 0.21, 1.45, 0.241244, 0.006, 1.02, 0.118335 },
		{1, 9, 0.725, 1.4, 0.21, 0.85, 0.03, 0.6, 1.25, 0.75, 1.2 },
		// {13660, 0.81, 1.7, 0.187523, 0.006, 1.03634, 0.123049, 0.98, 0.314041, 0.985, 0.179382 },
		{1, 11, },
		{1, 14, },
		{0, 4, },
		{0, 7, },
		{0, 11, },
		{0, 14, },
		{0, 8, },
		{0, 9, }

	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-02
	for( size_t i = 0; i < cuts.size(); i++){
		if( cuts[i][0] == target_int && cuts[i][1] == kinematic ){
			if( param == "Ep_min" ){lookup_val = cuts[i][2];}
			if( param == "Ep_max" ){lookup_val = cuts[i][3];}
			if( param == "PS_min" ){lookup_val = cuts[i][4];}
			if( param == "SH_PS_min" ){lookup_val = cuts[i][5];}
			if( param == "HCal_min" ){lookup_val = cuts[i][6];}
			if( param == "W2_min" ){lookup_val = cuts[i][7];}
			if( param == "W2_max" ){lookup_val = cuts[i][8];}
			if( param == "W2_mean" ){lookup_val = cuts[i][9];}
			if( param == "W2_sigma" ){lookup_val = cuts[i][10];}
			if( param == "W_mean" ){lookup_val = cuts[i][11];}
			if( param == "W_sigma" ){lookup_val = cuts[i][12];}

		}
	}
	return lookup_val;
}

double lookup_dxdy(int runnum, TString param){
	double lookup_val = -1;
	//Cuts vector:
// 	runnum, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma, dx_pn_max
	vector<vector<double>> dxdy_vec = {
		{11493,},
		{11494,},
		{13479,},
		{13544,},
		{13566, -0.787787, 0.176623, 0.0352742, 0.169769, -0.00455544, 0.258216, 1.16945},
		{13585, 0.9, 0.09, 0.62, 0.15, 0.9, 0.09, 0.62, 0.15, 0.65}

	};

	for( size_t i = 0; i < dxdy_vec.size(); i++){
		if( dxdy_vec[i][0] == runnum ){
			if( param == "dx_p" ){lookup_val = dxdy_vec[i][1];}
			if( param == "dx_p_sigma" ){lookup_val = dxdy_vec[i][2];}
			if( param == "dx_n" ){lookup_val = dxdy_vec[i][3];}
			if( param == "dx_n_sigma" ){lookup_val = dxdy_vec[i][4];}
			if( param == "dy" ){lookup_val = dxdy_vec[i][5];}
			if( param == "dy_sigma" ){lookup_val = dxdy_vec[i][6];}
			if( param == "dx_pn_max" ){lookup_val = dxdy_vec[i][7];}
		}
	}
	return lookup_val;
}


double lookup_parsed_dxdy(TString target, int kine, int sbsfield, TString param){
	double lookup_val = -1;
	int target_int;
	if( target == "LH2" ){
		target_int = 0;
	}
	if( target == "LD2" ){
		target_int = 1;
	}
	//Cuts vector:
// 	runnum, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma, dx_pn_max
	vector<vector<double>> dxdy_vec = {
	
		{1, 4, -0.787787, 0.176623, 0.0352742, 0.169769, -0.00455544, 0.258216, 1.16945},
		{1, 8, -0.792278, 0.204108, 0.0789937, 0.170958, -0.0126515, 0.260425, 1.24634}

	};

	for( size_t i = 0; i < dxdy_vec.size(); i++){
		if( dxdy_vec[i][0] == target_int && dxdy_vec[i][1] == kine && dxdy_vec[i][2] == sbsfield ){
			if( param == "dx_p" ){lookup_val = dxdy_vec[i][2];}
			if( param == "dx_p_sigma" ){lookup_val = dxdy_vec[i][3];}
			if( param == "dx_n" ){lookup_val = dxdy_vec[i][4];}
			if( param == "dx_n_sigma" ){lookup_val = dxdy_vec[i][5];}
			if( param == "dy" ){lookup_val = dxdy_vec[i][6];}
			if( param == "dy_sigma" ){lookup_val = dxdy_vec[i][7];}
			if( param == "dx_pn_max" ){lookup_val = dxdy_vec[i][8];}
		}
	}
	return lookup_val;
}

double lookup_errors(int runnum, TString param){
	double lookup_val = -1;
	//Errors vector:
// 	runnum, Ep, Ep_sigma, W2, W2_sigma, W, W_sigma, SH_PS, PS_SH_sigma

	vector<vector<double>> errors = {
		{11494, 6.69871E-05, 5.73715E-5, 4.05074E-3, 9.22625E-3, 1.91140E-3, 4.81047E-3 },
		{11595, 1.33104E-4, 1.22936E0-4, 4.56945E-3, 9.42378E-3, 2.68183E-3, 4.83629E-3},
		{13655, 4.74332e-5, 4.14841e-05, -1, -1, -1, -1},
		{13585, 3.68758E-5, 3.54258E-5, -1, -1, -1, -1, 2.12500e-0, 2.12500e-0 }
	};

	for( size_t i = 0; i < errors.size(); i++){
		if( errors[i][0] == runnum ){
			if( param == "Ep" ){lookup_val = errors[i][1];}
			if( param == "Ep_sigma" ){lookup_val = errors[i][2];}
			if( param == "W2" ){lookup_val = errors[i][3];}
			if( param == "W2_sigma" ){lookup_val = errors[i][4];}
			if( param == "W" ){lookup_val = errors[i][5];}
			if( param == "W_sigma" ){lookup_val = errors[i][6];}

		}
	}
	return lookup_val;

}

double lookup_XTALK_cut(int runnum, TString param){
	double lookup_val = -1;
	//Cuts vector:
// 	runnum, PS_clus_e_cut, SH_PS_clus_e_cut, SH_PS_sigma, HCal_clus_e_cut, Ep, Ep_sigma, W2, W2_sigma, W, W_sigma, 
	vector<vector<double>> xtalk_cuts = {
		// { 13566, 0.72, 3.21765, 0.448265, 0.07, 1.02202, 0.0749792 }// PRIM
		// { 13566, 0.25, 3.21723, 0.450056, 0.08, 1.02126, 0.0751095, 0.0751095, 0.99109, 0.291603, 1, 0.168658 } //XTALK_ON
		{13566, 1.02859, 0.0765564, 0.15, 2.95883, 0.417844, 0.015, 1.00188, 0.291855, 1.03, 0.178035 } //STD
	};
	//	9.38724e-01 4.57268e-01
//13585: HCal_clus_e_cut = 2.39514e-02 --> sigma: 1.08801e-02
	for( size_t i = 0; i < xtalk_cuts.size(); i++){
		if( xtalk_cuts[i][0] == runnum ){
			if( param == "Ep" ){lookup_val = xtalk_cuts[i][1];}
			if( param == "Ep_sigma" ){lookup_val = xtalk_cuts[i][2];}
			if( param == "PS_clus_e_cut" ){lookup_val = xtalk_cuts[i][3];}
			if( param == "SH_PS_clus_e_cut" ){lookup_val = xtalk_cuts[i][4];}
			if( param == "SH_PS_sigma" ){lookup_val = xtalk_cuts[i][5];}
			if( param == "HCal_clus_e_cut" ){lookup_val = xtalk_cuts[i][6];}
			if( param == "W2" ){lookup_val = xtalk_cuts[i][7];}
			if( param == "W2_sigma" ){lookup_val = xtalk_cuts[i][8];}
			if( param == "W" ){lookup_val = xtalk_cuts[i][9];}
			if( param == "W_sigma" ){lookup_val = xtalk_cuts[i][10];}

		}
	}
	return lookup_val;
}

int lookup_parsed_runs_cnt(TString target, int kine, int sbsfield){
	int return_var = -1, target_int = -1;

	//Targets: LH2 == 0, LD2 == 1
	if( target == "LH2"){
		target_int = 0;
	}
	if( target == "LD2"){
		target_int = 1;
	}

	vector<vector<int>> parsed_run_cnts = {
		{1, 4, 0, 6},
		{1, 4, 30, 17},
		{1, 4, 50, 3},
		{1, 7, 85, 42},
		{1, 8, 0, 10},
		{1, 8, 50, 6},
		{1, 8, 70, 49},
		{1, 8, 100, 12},
		{1, 9, 70, 68},
		{1, 11, 0, 3},		
		{1, 11, 100, 300},
		{1, 14, 70, 52},
		{0, 4, 30, 5},
		{0, 4, 0, 3},
		{0, 4, 50, 3},
		{0, 7, 85, 16},
		{0, 11, 100, 75},
		{0, 11, 0, 2},
		{0, 14, 70, 20},
		{0, 14, 0, 4},
		{0, 8, 70, 19},
		{0, 8, 0, 7},
		{0, 8, 100, 7},
		{0, 8, 50, 1},
		{0, 9, 70, 15}
	};

	for( size_t i = 0; i < parsed_run_cnts.size(); i++){
		if( parsed_run_cnts[i][0] == target_int && parsed_run_cnts[i][1] == kine && parsed_run_cnts[i][2] == sbsfield ){
			return_var = parsed_run_cnts[i][3];
		}
	}
	return return_var;

}

int lookup_parsed_runnums(TString target, int kine, int sbsfield, int run_index){
	int run_cnt = lookup_parsed_runs_cnt(target, kine, sbsfield);
	int lookup_run = -1;

	int target_int = -1;

	//Targets: LH2 == 0, LD2 == 1
	if( target == "LH2"){
		target_int = 0;
	}
	if( target == "LD2"){
		target_int = 1;
	}

	vector<vector<int>> parsed_runnums = {
		{1, 4, 0, 6, 11579, 11580, 11581, 11582, 11583, 11586},
		{1, 4, 30, 17, 11449, 11451, 11452, 11456, 11493, 11494, 11495, 11496, 11551, 11554, 11562, 11563, 11564, 11565, 11568, 11570, 11571},
		{1, 4, 50, 3, 11593, 11594, 11595},
		{1, 7, 85, 42, 11996, 11997, 11998, 11999, 12001, 12002, 12004, 12006, 12013, 12014, 12017, 12019, 12021, 12029, 12030, 12038, 12039, 
			12040, 12041, 12042, 12043, 12044, 12045, 12046, 12047, 12048, 12050, 12051, 12055, 12056, 12057, 12059, 12060, 12062, 12063, 12065, 12066, 
			12068, 12069, 12070, 12071, 12073},
		{1, 8, 0, 10, 13468, 13470, 13472, 13473, 13474, 13475, 13476, 13477, 13478, 13479},
		{1, 8, 50, 6, 13581, 13582, 13583, 13584, 13585, 13586},
		{1, 8, 70, 49, 13453, 13454, 13455, 13491, 13492, 13493, 13494, 13495, 13496, 13497, 13502, 13503, 13504, 13505, 13558, 13559, 13560, 13561, 13562,
			13563, 13564, 13565, 13566, 13567, 13568, 13569, 13570, 13571, 13587, 13588, 13589, 13590, 13591, 13592, 13593, 13596, 13597, 13608, 13609, 13610, 
			13612, 13613, 13614, 13615, 13616, 13617, 13618, 13619, 13620},
		{1, 8, 100, 12, 13544, 13545, 13546, 13547, 13548, 13549, 13550, 13551, 13552, 13554, 13556, 13557},
		{1, 9, 70, 68, 13660, 13661, 13662, 13664, 13665, 13666, 13677, 13678, 13679, 13680, 13681, 13682, 13684, 13685, 13686, 13687, 13688, 13689,
			13694, 13695, 13698, 13699, 13700, 13710, 13711, 13712, 13714, 13715, 13716, 13717, 13721, 13723, 13724, 13727, 13728, 13729, 13731, 13732,
			13734, 13736, 13737, 13746, 13748, 13749, 13753, 13754, 13755, 13756, 13757, 13758, 13760, 13761, 13764, 13765, 13766, 13767, 13770, 13771,
			13773, 13775, 13776, 13777, 13778, 13779, 13793, 13797, 13798, 13799},
		{1, 11, 0, 3, 12726, 12728, 12729},		
		{1, 11, 100, 300, 12314, 12315, 12316, 12318, 12319, 12321, 12322, 12323, 12333, 12337, 12338, 12339, 12341, 12342, 12343, 12344, 12346, 12347, 
			12348, 12349, 12356, 12357, 12359, 12360, 12362, 12364, 12366, 12371, 12372, 12373, 12374, 12375, 12376, 12377, 12378, 12379, 12383, 12384, 
			12385, 12387, 12388, 12394, 12395, 12396, 12403, 12404, 12405, 12406, 12407, 12408, 12409, 12410, 12411, 12412, 12416, 12423, 12424, 12425, 
			12473, 12474, 12478, 12479, 12480, 12481, 12482, 12485, 12486, 12498, 12499, 12500, 12501, 12502, 12512, 12513, 12514, 12515, 12516, 12517, 
			12521, 12523, 12527, 12532, 12534, 12535, 12536, 12537, 12538, 12539, 12541, 12542, 12543, 12544, 12547, 12549, 12550, 12551, 12552, 12553, 
			12554, 12556, 12558, 12559, 12565, 12566, 12567, 12568, 12569, 12570, 12571, 12572, 12574, 12577, 12580, 12581, 12582, 12583, 12584, 12585, 
			12586, 12589, 12590, 12592, 12593, 12597, 12598, 12599, 12614, 12618, 12619, 12620, 12623, 12627, 12628, 12638, 12639, 12641, 12642, 12643, 
			12644, 12645, 12647, 12653, 12654, 12655, 12656, 12659, 12660, 12661, 12662, 12664, 12666, 12671, 12673, 12674, 12678, 12679, 12680, 12681, 
			12682, 12683, 12684, 12686, 12687, 12688, 12693, 12694, 12695, 12697, 12699, 12700, 12701, 12702, 12706, 12707, 12708, 12709, 12710, 12711, 
			12713, 12714, 12715, 12721, 12722, 12732, 12733, 12734, 12735, 12737, 12738, 12739, 12740, 12742, 12743, 12767, 12768, 12769, 12770, 12771, 
			12772, 12773, 12778, 12779, 12781, 12782, 12796, 12797, 12798, 12799, 12801, 12802, 12813, 12814, 12815, 12816, 12817, 12819, 12820, 12824,
			12825, 12826, 12827, 12828, 12829, 12830, 12894, 12896, 12897, 12902, 12903, 12904, 12906, 12907, 12908, 12909, 12911, 12912, 12913, 12914, 
			12915, 12916, 12917, 12918, 12919, 12920, 12921, 12932, 12933, 12934, 12935, 12936, 12939, 12940, 12946, 12957, 12961, 12962, 12963, 12965, 
			12968, 12969, 12970, 12971, 12972, 12975, 12976, 12977, 12979, 12980, 12982, 12983, 12984, 12986, 12987, 12988, 12989, 13000, 13029, 13030, 
			13034, 13035, 13036, 13038, 13039, 13040, 13044, 13045, 13048, 13049, 13050, 13051, 13052, 13053, 13054, 13055, 13058, 13059, 13060, 13061, 13062, 13063},
		{1, 14, 70, 52, 13305, 13306, 13307, 13308, 13309, 13314, 13315, 13316, 13317, 13318, 13319, 13322, 13323, 13324, 13325, 13342, 13344, 13353, 13357,
			13358, 13359, 13360, 13361, 13362, 13363, 13364, 13368, 13369, 13370, 13371, 13372, 13373, 13381, 13382, 13384, 13385, 13387, 13388, 13389, 
			13390, 13391, 13392, 13393, 13394, 13395, 13398, 13399, 13400, 13402, 13403, 13406, 13407},
		{0, 4, 30, 5, 11436, 11500, 11547, 11548, 11616},
		{0, 4, 0, 3, 11573, 11587, 11588},
		{0, 4, 50, 3, 11589, 11590, 11592},
		{0, 7, 85, 16, 11989, 11990, 11991, 11992, 11993, 11994, 12000, 12008, 12022, 12035, 12049, 12052, 12053, 12058, 12064, 12072},
		{0, 11, 100, 75, 12313, 12320, 12335, 12336, 12340, 12345, 12355, 12358, 12363, 12367, 12368, 12369, 12370, 12380, 12382, 12400, 12401, 12414, 
			12415, 12427, 12428, 12429, 12471, 12472, 12491, 12496, 12497, 12524, 12525, 12526, 12548, 12561, 12563, 12575, 12576, 12587, 12588, 12625, 
			12626, 12657, 12658, 12675, 12676, 12677, 12961, 12692, 12704, 12705, 12716, 12747, 12758, 12765, 12776, 12777, 12809, 12810, 12812, 12822, 
			12823, 12895, 12910, 12923, 12930, 12931, 12959, 12960, 12973, 12974, 13027, 13028, 13041, 13042, 13043, 13056, 13057},
		{0, 11, 0, 2, 12730, 12731},
		{0, 14, 70, 20, 13239, 13240, 13241, 13242, 13243, 13244, 13312, 13313, 13320, 13321, 13345, 13346, 13348, 13349, 13351, 13352, 13379, 13396, 13397, 13405},
		{0, 14, 0, 4, 13375, 13376, 13377, 13378},
		{0, 8, 70, 19, 13444, 13450, 13451, 13452, 13482, 13483, 13484, 13485, 13486, 13487, 13488, 13489, 13490, 13573, 13574, 13575, 13576, 13577, 13578},
		{0, 8, 0, 7, 13459, 13460, 13461, 13463, 13464, 13465, 13466},
		{0, 8, 100, 7, 13537, 13538, 13539, 13540, 13541, 13542, 13543},
		{0, 8, 50, 1, 13580},
		{0, 9, 70, 15, 13656, 13657, 13663, 13676, 13683, 13696, 13697, 13719, 13720, 13730, 13747, 13768, 13769, 13795, 13796}
	};

	for( size_t i = 0; i < parsed_runnums.size(); i++){
		if( parsed_runnums[i][0] == target_int && parsed_runnums[i][1] == kine && parsed_runnums[i][2] == sbsfield ){
			// for(size_t run = parsed_runnums[i][4]; run < parsed_runnums[i].size(); run++){
			// 	runnum_vector.push_back(parsed_runnums[i][run]);
			// 	cout << parsed_runnums[i][run] << endl;
			// }
			// break;
			lookup_run = parsed_runnums[i][run_index + 4];
		}
	}
	
	return lookup_run;
}



double lookup_run_info( int runnum, TString lookup_var){
	double return_var = -99999;
	//runnum, beam_current, sbs_current, bb_current, bbcal_thresh, sbs_field, target (0 = LH2, 1 = LD2)
	
	vector<vector<double>> all_run_info ={	{11436,6.5,630,750,-500,30,0}
		,{11500,3.5,630,750,30,0}
		,{11547,3.5,630,750,-425,30,0}
		,{11548,3.5,630,750,-425,30,0}
		,{11573,3.5,0,750,-425,0,0}
		,{11587,3.5,0,750,-425,0,0}
		,{11588,3.5,0,750,-425,0,0}
		,{11589,3.5,1050,750,-425,50,0}
		,{11590,3.5,1050,750,-425,50,0}
		,{11592,3.5,1050,750,-425,50,0}
		,{11616,3.5,630,750,-425,30,0}
		,{11989,1,1785,750,-380,85,0}
		,{11990,2,1785,750,-380,85,0}
		,{11991,4,1785,750,-380,85,0}
		,{11992,4,1785,750,-380,85,0}
		,{11993,8,1785,750,-380,85,0}
		,{11994,8,1785,750,-380,85,0}
		,{12000,8,1785,750,-380,85,0}
		,{12008,8,1785,750,-380,85,0}
		,{12022,10,1785,750,-450,85,0}
		,{12035,10,1785,750,-450,85,0}
		,{12049,10,1785,750,-450,85,0}
		,{12052,10,1785,750,-450,85,0}
		,{12053,10,1785,750,-450,85,0}
		,{12058,10,1785,750,-450,85,0}
		,{12064,10,1785,750,-450,85,0}
		,{12072,10,1785,750,-450,85,0}
		,{12313,15,2100,750,400,100,0}
		,{12320,15,2100,750,400,100,0}
		,{12335,15,2100,750,450,100,0}
		,{12336,15,2100,750,450,100,0}
		,{12340,8,2100,750,450,100,0}
		,{12345,15,2100,750,450,100,0}
		,{12355,15,2100,750,450,100,0}
		,{12358,15,2100,750,450,100,0}
		,{12363,15,2100,750,450,100,0}
		,{12367,15,2100,750,450,100,0}
		,{12368,15,2100,750,450,100,0}
		,{12369,15,2100,750,450,100,0}
		,{12370,15,2100,750,450,100,0}
		,{12380,15,2100,750,450,100,0}
		,{12382,15,2100,750,450,100,0}
		,{12400,15,2100,750,450,100,0}
		,{12401,15,2100,750,450,100,0}
		,{12414,8,2100,750,450,100,0}
		,{12415,8,2100,750,450,100,0}
		,{12427,15,2100,750,450,100,0}
		,{12428,15,2100,750,450,100,0}
		,{12429,15,2100,750,450,100,0}
		,{12471,15,2100,750,450,100,0}
		,{12472,15,2100,750,450,100,0}
		,{12491,15,2100,750,450,100,0}
		,{12496,15,2100,750,450,100,0}
		,{12497,15,2100,750,450,100,0}
		,{12524,15,2100,750,451,100,0}
		,{12525,15,2100,750,451,100,0}
		,{12526,15,2100,750,451,100,0}
		,{12548,8,2100,750,450,100,0}
		,{12561,12,2100,750,470,100,0}
		,{12563,12,2100,750,470,100,0}
		,{12575,12,2100,750,470,100,0}
		,{12576,12,2100,750,470,100,0}
		,{12587,12,2100,750,470,100,0}
		,{12588,12,2100,750,470,100,0}
		,{12625,12,2100,750,470,100,0}
		,{12626,12,2100,750,470,100,0}
		,{12657,12,2100,750,470,100,0}
		,{12658,12,2100,750,470,100,0}
		,{12675,15,2100,750,470,100,0}
		,{12676,12,2100,750,470,100,0}
		,{12677,15,2100,750,470,100,0}
		,{12961,15,2100,750,470,100,0}
		,{12692,15,2100,750,470,100,0}
		,{12704,15,2100,750,480,100,0}
		,{12705,15,2100,750,480,100,0}
		,{12716,15,2100,750,480,100,0}
		,{12730,15,0,750,480,0,0}
		,{12731,15,0,750,480,0,0}
		,{12747,15,2100,750,480,100,0}
		,{12758,15,2100,750,480,100,0}
		,{12765,15,2100,750,480,100,0}
		,{12776,15,2100,750,480,100,0}
		,{12777,15,2100,750,480,100,0}
		,{12809,15,2100,750,470,100,0}
		,{12810,15,2100,750,470,100,0}
		,{12812,15,2100,750,480,100,0}
		,{12822,15,2100,750,480,100,0}
		,{12823,15,2100,750,480,100,0}
		,{12895,15,2100,750,480,100,0}
		,{12910,15,2100,750,480,100,0}
		,{12923,15,2100,750,480,100,0}
		,{12930,15,2100,750,480,100,0}
		,{12931,15,2100,750,480,100,0}
		,{12959,11,2100,750,480,100,0}
		,{12960,11,2100,750,480,100,0}
		,{12973,15,2100,750,480,100,0}
		,{12974,15,2100,750,480,100,0}
		,{13027,15,2100,750,480,100,0}
		,{13028,15,2100,750,480,100,0}
		,{13041,15,2100,750,480,100,0}
		,{13042,15,2100,750,480,100,0}
		,{13043,15,2100,750,480,100,0}
		,{13056,15,2100,750,480,100,0}
		,{13057,15,2100,750,480,100,0}
		,{13239,10,1470,750,-500,70,0}
		,{13240,1.7,1470,750,-500,70,0}
		,{13241,9.6,1470,750,-500,70,0}
		,{13242,9.6,1470,750,-500,70,0}
		,{13243,9.6,1470,750,-500,70,0}
		,{13244,9.6,1470,750,-500,70,0}
		,{13312,12,1466,750,500,70,0}
		,{13313,12,1466,750,500,70,0}
		,{13320,12,1466,750,500,70,0}
		,{13321,12,1466,750,500,70,0}
		,{13345,12,1466,750,500,70,0}
		,{13346,12,1466,750,500,70,0}
		,{13348,12,1466,750,500,70,0}
		,{13349,12,1466,750,500,70,0}
		,{13351,15,1466,750,500,70,0}
		,{13352,15,1466,750,500,70,0}
		,{13375,15,0,750,525,0,0}
		,{13376,15,0,750,525,0,0}
		,{13377,15,0,750,525,0,0}
		,{13378,15,0,750,525,0,0}
		,{13379,15,1466,750,525,70,0}
		,{13396,15,1466,750,525,70,0}
		,{13397,15,1466,750,525,70,0}
		,{13405,15,1466,750,525,70,0}
		,{13444,1,1470,750,-385,70,0}
		,{13450,4,1470,750,-385,70,0}
		,{13451,5,1470,750,-385,70,0}
		,{13452,5,1470,750,-405,70,0}
		,{13459,8,0,750,-459,0,0}
		,{13460,8,0,750,-459,0,0}
		,{13461,8,0,750,-459,0,0}
		,{13463,8,0,750,-459,0,0}
		,{13464,8,0,750,-459,0,0}
		,{13465,8,0,750,-459,0,0}
		,{13466,8,0,750,-459,0,0}
		,{13482,8,1470,750,-477,70,0}
		,{13483,8,1470,750,-477,70,0}
		,{13484,8,1470,750,-477,70,0}
		,{13485,8,1470,750,-477,70,0}
		,{13486,8,1470,750,-477,70,0}
		,{13487,8,1470,750,-477,70,0}
		,{13488,8,1470,750,-477,70,0}
		,{13489,8,1470,750,-477,70,0}
		,{13490,8,1470,750,-477,70,0}
		,{13537,8,2100,750,-452,100,0}
		,{13538,8,2100,750,-452,100,0}
		,{13539,8,2100,750,-452,100,0}
		,{13540,8,2100,750,-452,100,0}
		,{13541,8,2100,750,-452,100,0}
		,{13542,8,2100,750,-452,100,0}
		,{13543,8,2100,750,-452,100,0}
		,{13573,8,1470,750,477,70,0}
		,{13574,8,1470,750,477,70,0}
		,{13575,8,1470,750,477,70,0}
		,{13576,8,1470,750,477,70,0}
		,{13577,8,1470,750,477,70,0}
		,{13578,8,1470,750,477,70,0}
		,{13580,8,1050,750,447,50,0}
		,{13656,3,1470,750,-428,70,0}
		,{13657,3,1470,750,-428,70,0}
		,{13663,10,1470,750,-512,70,0}
		,{13676,10,1470,750,-512,70,0}
		,{13683,15,1470,750,-554,70,0}
		,{13696,15,1470,750,-554,70,0}
		,{13697,15,1470,750,-554,70,0}
		,{13719,15,1470,750,-554,70,0}
		,{13720,15,1470,750,-554,70,0}
		,{13730,15,1470,750,-554,70,0}
		,{13747,15,1470,750,-554,70,0}
		,{13768,15,1470,750,-560,70,0}
		,{13769,15,1470,750,-560,70,0}
		,{13795,15,1470,750,-553,70,0}
		,{13796,15,1470,750,-553,70,0} 
		,{11449,3,630,750,-500,30,1}
		,{11451,3,630,750,-475,30,1}
		,{11452,3,630,750,-500,30,1}
		,{11456,1.75,630,750,-425,30,1}
		,{11493,1,630,750,-425,30,1}
		,{11494,1,630,750,-425,30,1}
		,{11495,1,630,750,-425,30,1}
		,{11496,1,630,750,-425,30,1}
		,{11551,1.75,630,750,-425,30,1}
		,{11554,1.75,630,750,-425,30,1}
		,{11562,1.75,630,750,-425,30,1}
		,{11563,1.75,630,750,-425,30,1}
		,{11564,1.75,630,750,-425,30,1}
		,{11565,1.75,630,750,-425,30,1}
		,{11568,1.75,630,750,-425,30,1}
		,{11570,1.75,630,750,-425,30,1}
		,{11571,1.75,630,750,-425,30,1}
		,{11579,1.75,0,750,-425,0,1}
		,{11580,1.75,0,750,-425,0,1}
		,{11581,1.75,0,750,-425,0,1}
		,{11582,1.75,0,750,-425,0,1}
		,{11583,1.75,0,750,-425,0,1}
		,{11586,1.75,0,750,-425,0,1}
		,{11593,1.75,1050,750,-425,50,1}
		,{11594,1.75,1050,750,-425,50,1}
		,{11595,1.75,1050,750,-425,50,1}
		,{11996,4,1785,750,-380,85,1}
		,{11997,4,1785,750,-380,85,1}
		,{11998,4,1785,750,-380,85,1}
		,{11999,4,1785,750,-380,85,1}
		,{12001,4,1785,750,-380,85,1}
		,{12002,4,1785,750,-380,85,1}
		,{12004,4,1785,750,-380,85,1}
		,{12006,4,1785,750,-380,85,1}
		,{12013,6,1785,750,-450,85,1}
		,{12014,4,1785,750,-450,85,1}
		,{12017,4,1785,750,-450,85,1}
		,{12019,5,1785,750,-450,85,1}
		,{12021,5,1785,750,-450,85,1}
		,{12029,5,1785,750,-450,85,1}
		,{12030,5,1785,750,-450,85,1}
		,{12038,5,1785,750,-450,85,1}
		,{12039,5,1785,750,-450,85,1}
		,{12040,5,1785,750,-450,85,1}
		,{12041,5,1785,750,-450,85,1}
		,{12042,5,1785,750,-450,85,1}
		,{12043,5,1785,750,-450,85,1}
		,{12044,5,1785,750,-450,85,1}
		,{12045,5,1785,750,-450,85,1}
		,{12046,5,1785,750,-450,85,1}
		,{12047,5,1785,750,-450,85,1}
		,{12048,5,1785,750,-450,85,1}
		,{12050,5,1785,750,-450,85,1}
		,{12051,5,1785,750,-450,85,1}
		,{12055,5,1785,750,-450,85,1}
		,{12056,5,1785,750,-450,85,1}
		,{12057,5,1785,750,-450,85,1}
		,{12059,5,1785,750,-450,85,1}
		,{12060,5,1785,750,-450,85,1}
		,{12062,5,1785,750,-450,85,1}
		,{12063,5,1785,750,-450,85,1}
		,{12065,5,1785,750,-450,85,1}
		,{12066,5,1785,750,-450,85,1}
		,{12068,5,1785,750,-450,85,1}
		,{12069,5,1785,750,-450,85,1}
		,{12070,5,1785,750,-450,85,1}
		,{12071,5,1785,750,-450,85,1}
		,{12073,5,1785,750,-450,85,1}
		,{12314,8,2100,750,400,100,1}
		,{12315,8,2100,750,450,100,1}
		,{12316,8,2100,750,450,100,1}
		,{12318,8,2100,750,450,100,1}
		,{12319,8,2100,750,450,100,1}
		,{12321,8,2100,750,450,100,1}
		,{12322,8,2100,750,450,100,1}
		,{12323,8,2100,750,450,100,1}
		,{12333,8,2100,750,450,100,1}
		,{12337,8,2100,750,450,100,1}
		,{12338,8,2100,750,450,100,1}
		,{12339,8,2100,750,450,100,1}
		,{12341,4,2100,750,450,100,1}
		,{12342,4,2100,750,450,100,1}
		,{12343,4,2100,750,450,100,1}
		,{12344,4,2100,750,450,100,1}
		,{12346,8,2100,750,450,100,1}
		,{12347,8,2100,750,450,100,1}
		,{12348,8,2100,750,450,100,1}
		,{12349,8,2100,750,450,100,1}
		,{12356,8,2100,750,450,100,1}
		,{12357,8,2100,750,450,100,1}
		,{12359,8,2100,750,450,100,1}
		,{12360,8,2100,750,450,100,1}
		,{12362,8,2100,750,450,100,1}
		,{12364,8,2100,750,450,100,1}
		,{12366,8,2100,750,450,100,1}
		,{12371,8,2100,750,450,100,1}
		,{12372,8,2100,750,450,100,1}
		,{12373,8,2100,750,450,100,1}
		,{12374,8,2100,750,450,100,1}
		,{12375,8,2100,750,450,100,1}
		,{12376,8,2100,750,450,100,1}
		,{12377,8,2100,750,450,100,1}
		,{12378,8,2100,750,450,100,1}
		,{12379,8,2100,750,450,100,1}
		,{12383,8,2100,750,450,100,1}
		,{12384,8,2100,750,450,100,1}
		,{12385,8,2100,750,450,100,1}
		,{12387,8,2100,750,450,100,1}
		,{12388,8,2100,750,450,100,1}
		,{12394,8,2100,750,450,100,1}
		,{12395,8,2100,750,450,100,1}
		,{12396,8,2100,750,450,100,1}
		,{12403,8,2100,750,450,100,1}
		,{12404,8,2100,750,450,100,1}
		,{12405,8,2100,750,450,100,1}
		,{12406,8,2100,750,450,100,1}
		,{12407,8,2100,750,450,100,1}
		,{12408,8,2100,750,450,100,1}
		,{12409,8,2100,750,450,100,1}
		,{12410,8,2100,750,450,100,1}
		,{12411,8,2100,750,450,100,1}
		,{12412,8,2100,750,450,100,1}
		,{12416,8,2100,750,450,100,1}
		,{12423,2,2100,750,450,100,1}
		,{12424,4,2100,750,450,100,1}
		,{12425,4,2100,750,450,100,1}
		,{12473,8,2100,750,450,100,1}
		,{12474,8,2100,750,450,100,1}
		,{12478,8,2100,750,450,100,1}
		,{12479,8,2100,750,450,100,1}
		,{12480,8,2100,750,450,100,1}
		,{12481,8,2100,750,450,100,1}
		,{12482,8,2100,750,450,100,1}
		,{12485,8,2100,750,450,100,1}
		,{12486,8,2100,750,450,100,1}
		,{12498,8,2100,750,450,100,1}
		,{12499,8,2100,750,450,100,1}
		,{12500,8,2100,750,450,100,1}
		,{12501,8,2100,750,450,100,1}
		,{12502,8,2100,750,450,100,1}
		,{12512,8,2100,750,450,100,1}
		,{12513,8,2100,750,450,100,1}
		,{12514,8,2100,750,450,100,1}
		,{12515,8,2100,750,450,100,1}
		,{12516,8,2100,750,451,100,1}
		,{12517,8,2100,750,451,100,1}
		,{12521,8,2100,750,451,100,1}
		,{12523,8,2100,750,451,100,1}
		,{12527,8,2100,750,451,100,1}
		,{12532,8,2100,750,451,100,1}
		,{12534,8,2100,750,451,100,1}
		,{12535,8,2100,750,451,100,1}
		,{12536,8,2100,750,451,100,1}
		,{12537,8,2100,750,451,100,1}
		,{12538,8,2100,750,451,100,1}
		,{12539,8,2100,750,451,100,1}
		,{12541,8,2100,750,451,100,1}
		,{12542,8,2100,750,451,100,1}
		,{12543,8,2100,750,451,100,1}
		,{12544,8,2100,750,451,100,1}
		,{12547,8,2100,750,450,100,1}
		,{12549,8,2100,750,450,100,1}
		,{12550,8,2100,750,450,100,1}
		,{12551,8,2100,750,450,100,1}
		,{12552,8,2100,750,450,100,1}
		,{12553,8,2100,750,450,100,1}
		,{12554,8,2100,750,450,100,1}
		,{12556,8,2100,750,450,100,1}
		,{12558,8,2100,750,450,100,1}
		,{12559,10,2100,750,470,100,1}
		,{12565,10,2100,750,470,100,1}
		,{12566,10,2100,750,470,100,1}
		,{12567,10,2100,750,470,100,1}
		,{12568,10,2100,750,470,100,1}
		,{12569,10,2100,750,470,100,1}
		,{12570,10,2100,750,470,100,1}
		,{12571,10,2100,750,470,100,1}
		,{12572,10,2100,750,470,100,1}
		,{12574,10,2100,750,470,100,1}
		,{12577,10,2100,750,470,100,1}
		,{12580,10,2100,750,470,100,1}
		,{12581,10,2100,750,470,100,1}
		,{12582,10,2100,750,470,100,1}
		,{12583,10,2100,750,470,100,1}
		,{12584,10,2100,750,470,100,1}
		,{12585,10,2100,750,470,100,1}
		,{12586,10,2100,750,470,100,1}
		,{12589,10,2100,750,470,100,1}
		,{12590,10,2100,750,470,100,1}
		,{12592,10,2100,750,470,100,1}
		,{12593,10,2100,750,470,100,1}
		,{12597,10,2100,750,470,100,1}
		,{12598,10,2100,750,470,100,1}
		,{12599,10,2100,750,470,100,1}
		,{12614,10,2100,750,470,100,1}
		,{12618,10,2100,750,470,100,1}
		,{12619,10,2100,750,470,100,1}
		,{12620,10,2100,750,470,100,1}
		,{12623,10,2100,750,470,100,1}
		,{12627,10,2100,750,470,100,1}
		,{12628,10,2100,750,470,100,1}
		,{12638,10,2100,750,470,100,1}
		,{12639,10,2100,750,470,100,1}
		,{12641,10,2100,750,470,100,1}
		,{12642,10,2100,750,460,100,1}
		,{12643,10,2100,750,460,100,1}
		,{12644,10,2100,750,460,100,1}
		,{12645,10,2100,750,460,100,1}
		,{12647,10,2100,750,470,100,1}
		,{12653,10,2100,750,470,100,1}
		,{12654,1,2100,750,300,100,1}
		,{12655,10,2100,750,470,100,1}
		,{12656,10,2100,750,470,100,1}
		,{12659,10,2100,750,470,100,1}
		,{12660,10,2100,750,470,100,1}
		,{12661,10,2100,750,470,100,1}
		,{12662,10,2100,750,470,100,1}
		,{12664,10,2100,750,470,100,1}
		,{12666,10,2100,750,470,100,1}
		,{12671,10,2100,750,470,100,1}
		,{12673,10,2100,750,470,100,1}
		,{12674,11,2100,750,470,100,1}
		,{12678,11,2100,750,470,100,1}
		,{12679,11,2100,750,470,100,1}
		,{12680,11,2100,750,470,100,1}
		,{12681,11,2100,750,470,100,1}
		,{12682,11,2100,750,470,100,1}
		,{12683,11,2100,750,470,100,1}
		,{12684,11,2100,750,470,100,1}
		,{12686,11,2100,750,470,100,1}
		,{12687,11,2100,750,470,100,1}
		,{12688,11,2100,750,470,100,1}
		,{12693,11,2100,750,470,100,1}
		,{12694,11,2100,750,470,100,1}
		,{12695,11,2100,750,470,100,1}
		,{12697,11,2100,750,470,100,1}
		,{12699,11,2100,750,470,100,1}
		,{12700,11,2100,750,480,100,1}
		,{12701,11,2100,750,480,100,1}
		,{12702,11,2100,750,480,100,1}
		,{12706,11,2100,750,480,100,1}
		,{12707,11,2100,750,480,100,1}
		,{12708,11,2100,750,480,100,1}
		,{12709,11,2100,750,480,100,1}
		,{12710,11,2100,750,480,100,1}
		,{12711,11,2100,750,480,100,1}
		,{12713,11,2100,750,480,100,1}
		,{12714,11,2100,750,480,100,1}
		,{12715,11,2100,750,480,100,1}
		,{12721,11,2100,750,480,100,1}
		,{12722,11,2100,750,480,100,1}
		,{12726,11,0,750,480,0,1}
		,{12728,11,0,750,480,0,1}
		,{12729,11,0,750,480,0,1}
		,{12732,11,2100,750,480,100,1}
		,{12733,11,2100,750,480,100,1}
		,{12734,11,2100,750,480,100,1}
		,{12735,11,2100,750,480,100,1}
		,{12737,11,2100,750,480,100,1}
		,{12738,11,2100,750,480,100,1}
		,{12739,11,2100,750,480,100,1}
		,{12740,11,2100,750,480,100,1}
		,{12742,11,2100,750,480,100,1}
		,{12743,11,2100,750,480,100,1}
		,{12767,11,2100,750,480,100,1}
		,{12768,11,2100,750,480,100,1}
		,{12769,11,2100,750,480,100,1}
		,{12770,11,2100,750,480,100,1}
		,{12771,11,2100,750,480,100,1}
		,{12772,11,2100,750,480,100,1}
		,{12773,11,2100,750,480,100,1}
		,{12778,11,2100,750,480,100,1}
		,{12779,11,2100,750,480,100,1}
		,{12781,11,2100,750,480,100,1}
		,{12782,11,2100,750,480,100,1}
		,{12796,11,2100,750,480,100,1}
		,{12797,11,2100,750,480,100,1}
		,{12798,11,2100,750,480,100,1}
		,{12799,11,2100,750,480,100,1}
		,{12801,11,2100,750,480,100,1}
		,{12802,11,2100,750,480,100,1}
		,{12813,11,2100,750,480,100,1}
		,{12814,11,2100,750,480,100,1}
		,{12815,11,2100,750,480,100,1}
		,{12816,11,2100,750,480,100,1}
		,{12817,11,2100,750,480,100,1}
		,{12819,11,2100,750,480,100,1}
		,{12820,11,2100,750,480,100,1}
		,{12824,11,2100,750,480,100,1}
		,{12825,11,2100,750,480,100,1}
		,{12826,11,2100,750,480,100,1}
		,{12827,11,2100,750,480,100,1}
		,{12828,11,2100,750,480,100,1}
		,{12829,11,2100,750,480,100,1}
		,{12830,11,2100,750,480,100,1}
		,{12894,11,2100,750,480,100,1}
		,{12896,11,2100,750,480,100,1}
		,{12897,11,2100,750,480,100,1}
		,{12902,11,2100,750,480,100,1}
		,{12903,11,2100,750,480,100,1}
		,{12904,11,2100,750,480,100,1}
		,{12906,11,2100,750,480,100,1}
		,{12907,11,2100,750,480,100,1}
		,{12908,11,2100,750,480,100,1}
		,{12909,11,2100,750,480,100,1}
		,{12911,11,2100,750,480,100,1}
		,{12912,11,2100,750,480,100,1}
		,{12913,11,2100,750,480,100,1}
		,{12914,11,2100,750,480,100,1}
		,{12915,11,2100,750,480,100,1}
		,{12916,11,2100,750,480,100,1}
		,{12917,11,2100,750,480,100,1}
		,{12918,11,2100,750,480,100,1}
		,{12919,11,2100,750,480,100,1}
		,{12920,11,2100,750,480,100,1}
		,{12921,11,2100,750,480,100,1}
		,{12932,11,2100,750,480,100,1}
		,{12933,11,2100,750,480,100,1}
		,{12934,11,2100,750,480,100,1}
		,{12935,11,2100,750,480,100,1}
		,{12936,11,2100,750,480,100,1}
		,{12939,11,2100,750,480,100,1}
		,{12940,11,2100,750,480,100,1}
		,{12946,11,2100,750,480,100,1}
		,{12957,11,2100,750,480,100,1}
		,{12961,11,2100,750,480,100,1}
		,{12962,11,2100,750,480,100,1}
		,{12963,11,2100,750,480,100,1}
		,{12965,11,2100,750,480,100,1}
		,{12968,11,2100,750,480,100,1}
		,{12969,11,2100,750,480,100,1}
		,{12970,11,2100,750,480,100,1}
		,{12971,11,2100,750,480,100,1}
		,{12972,11,2100,750,480,100,1}
		,{12975,11,2100,750,480,100,1}
		,{12976,11,2100,750,480,100,1}
		,{12977,11,2100,750,480,100,1}
		,{12979,11,2100,750,480,100,1}
		,{12980,11,2100,750,480,100,1}
		,{12982,11,2100,750,480,100,1}
		,{12983,11,2100,750,480,100,1}
		,{12984,11,2100,750,480,100,1}
		,{12986,11,2100,750,480,100,1}
		,{12987,11,2100,750,480,100,1}
		,{12988,11,2100,750,480,100,1}
		,{12989,11,2100,750,480,100,1}
		,{13000,11,2100,750,480,100,1}
		,{13029,11,2100,750,480,100,1}
		,{13030,11,2100,750,480,100,1}
		,{13034,11,2100,750,480,100,1}
		,{13035,11,2100,750,480,100,1}
		,{13036,11,2100,750,480,100,1}
		,{13038,11,2100,750,480,100,1}
		,{13039,11,2100,750,480,100,1}
		,{13040,11,2100,750,480,100,1}
		,{13044,11,2100,750,480,100,1}
		,{13045,11,2100,750,480,100,1}
		,{13048,11,2100,750,480,100,1}
		,{13049,11,2100,750,480,100,1}
		,{13050,11,2100,750,480,100,1}
		,{13051,11,2100,750,480,100,1}
		,{13052,11,2100,750,480,100,1}
		,{13053,11,2100,750,480,100,1}
		,{13054,11,2100,750,480,100,1}
		,{13055,11,2100,750,480,100,1}
		,{13058,11,2100,750,480,100,1}
		,{13059,11,2100,750,480,100,1}
		,{13060,11,2100,750,480,100,1}
		,{13061,11,2100,750,480,100,1}
		,{13062,11,2100,750,480,100,1}
		,{13063,11,2100,750,480,100,1}
		,{13305,8,1470,750,500,70,1}
		,{13306,8,1470,750,500,70,1}
		,{13307,8,1470,750,500,70,1}
		,{13308,8,1470,750,500,70,1}
		,{13309,8,1470,750,500,70,1}
		,{13314,8,1466,750,500,70,1}
		,{13315,8,1466,750,500,70,1}
		,{13316,8,1466,750,500,70,1}
		,{13317,8,1466,750,500,70,1}
		,{13318,8,1466,750,500,70,1}
		,{13319,8,1466,750,500,70,1}
		,{13322,8,1466,750,500,70,1}
		,{13323,8,1466,750,500,70,1}
		,{13324,8,1466,750,500,70,1}
		,{13325,8,1466,750,500,70,1}
		,{13342,8,1466,750,500,70,1}
		,{13344,8,1466,750,500,70,1}
		,{13353,10,1466,750,500,70,1}
		,{13357,10,1466,750,525,70,1}
		,{13358,10,1466,750,525,70,1}
		,{13359,10,1466,750,525,70,1}
		,{13360,10,1466,750,525,70,1}
		,{13361,10,1466,750,525,70,1}
		,{13362,10,1466,750,525,70,1}
		,{13363,10,1466,750,525,70,1}
		,{13364,10,1466,750,525,70,1}
		,{13368,10,1466,750,525,70,1}
		,{13369,10,1466,750,525,70,1}
		,{13370,10,1466,750,525,70,1}
		,{13371,10,1466,750,525,70,1}
		,{13372,10,1466,750,525,70,1}
		,{13373,10,1466,750,525,70,1}
		,{13381,10,1466,750,525,70,1}
		,{13382,10,1466,750,525,70,1}
		,{13384,10,1466,750,525,70,1}
		,{13385,10,1466,750,525,70,1}
		,{13387,10,1466,750,525,70,1}
		,{13388,10,1466,750,525,70,1}
		,{13389,10,1466,750,525,70,1}
		,{13390,10,1466,750,525,70,1}
		,{13391,10,1466,750,525,70,1}
		,{13392,10,1466,750,525,70,1}
		,{13393,10,1466,750,525,70,1}
		,{13394,10,1466,750,525,70,1}
		,{13395,10,1466,750,525,70,1}
		,{13398,10,1466,750,525,70,1}
		,{13399,10,1466,750,525,70,1}
		,{13400,10,1466,750,525,70,1}
		,{13402,10,1466,750,525,70,1}
		,{13403,10,1466,750,525,70,1}
		,{13406,10,1466,750,525,70,1}
		,{13407,10,1466,750,525,70,1}
		,{13453,3,1470,750,-405,70,1}
		,{13454,3,1470,750,-405,70,1}
		,{13455,3,1470,750,-425,70,1}
		,{13468,5,0,750,-477,0,1}
		,{13470,5,0,750,-477,0,1}
		,{13472,5,0,750,-477,0,1}
		,{13473,5,0,750,-477,0,1}
		,{13474,5,0,750,-477,0,1}
		,{13475,5,0,750,-477,0,1}
		,{13476,5,0,750,-477,0,1}
		,{13477,5,0,750,-477,0,1}
		,{13478,5,0,750,-477,0,1}
		,{13479,5,0,750,-477,0,1}
		,{13491,5,1466,750,-477,70,1}
		,{13492,5,1466,750,-477,70,1}
		,{13493,5,1466,750,-477,70,1}
		,{13494,5,1466,750,-477,70,1}
		,{13495,5,1466,750,-477,70,1}
		,{13496,5,1466,750,-477,70,1}
		,{13497,5,1466,750,-477,70,1}
		,{13502,5,1466,750,-477,70,1}
		,{13503,5,1466,750,-477,70,1}
		,{13504,5,1466,750,-477,70,1}
		,{13505,5,1466,750,-477,70,1}
		,{13544,5,2100,750,-462,100,1}
		,{13545,5,2100,750,-462,100,1}
		,{13546,5,2100,750,-462,100,1}
		,{13547,5,2100,750,-462,100,1}
		,{13548,5,2100,750,-462,100,1}
		,{13549,5,2100,750,-462,100,1}
		,{13550,5,2100,750,-462,100,1}
		,{13551,5,2100,750,-462,100,1}
		,{13552,5,2100,750,-462,100,1}
		,{13554,4.5,2100,750,-462,100,1}
		,{13556,5,2100,750,-462,100,1}
		,{13557,5,2100,750,-462,100,1}
		,{13558,5,1470,750,477,70,1}
		,{13559,5,1470,750,477,70,1}
		,{13560,5,1470,750,477,70,1}
		,{13561,5,1470,750,477,70,1}
		,{13562,5,1470,750,477,70,1}
		,{13563,5,1470,750,477,70,1}
		,{13564,5,1470,750,477,70,1}
		,{13565,5,1470,750,477,70,1}
		,{13566,5,1470,750,477,70,1}
		,{13567,5,1470,750,477,70,1}
		,{13568,5,1470,750,477,70,1}
		,{13569,5,1470,750,477,70,1}
		,{13570,5,1470,750,477,70,1}
		,{13571,5,1470,750,477,70,1}
		,{13581,5,1050,750,458,50,1}
		,{13582,5,1050,750,458,50,1}
		,{13583,5,1050,750,458,50,1}
		,{13584,5,1050,750,458,50,1}
		,{13585,5,1050,750,458,50,1}
		,{13586,5,1050,750,458,50,1}
		,{13587,5,1466,750,473,70,1}
		,{13588,5,1466,750,473,70,1}
		,{13589,5,1466,750,473,70,1}
		,{13590,5,1466,750,473,70,1}
		,{13591,5,1466,750,478,70,1}
		,{13592,5,1466,750,478,70,1}
		,{13593,5,1466,750,478,70,1}
		,{13596,5,1470,750,478,70,1}
		,{13597,5,1470,750,478,70,1}
		,{13608,5.5,1470,750,478,70,1}
		,{13609,5.5,1470,750,478,70,1}
		,{13610,5.5,1470,750,478,70,1}
		,{13612,5.5,1470,750,478,70,1}
		,{13613,5.5,1470,750,478,70,1}
		,{13614,5.5,1470,750,478,70,1}
		,{13615,5.5,1470,750,478,70,1}
		,{13616,5.5,1470,750,478,70,1}
		,{13617,5.5,1470,750,478,70,1}
		,{13618,5.5,1470,750,478,70,1}
		,{13619,5.5,1470,750,478,70,1}
		,{13620,5.5,1470,750,478,70,1}
		,{13660,1.5,1470,750,-428,70,1}
		,{13661,4.4,1470,750,-512,70,1}
		,{13662,5,1470,750,-512,70,1}
		,{13664,5,1470,750,-512,70,1}
		,{13665,5,1470,750,-512,70,1}
		,{13666,5,1470,750,-512,70,1}
		,{13677,10,1470,750,-585,70,1}
		,{13678,10,1470,750,-585,70,1}
		,{13679,10,1470,750,-585,70,1}
		,{13680,10,1470,750,-585,70,1}
		,{13681,10,1470,750,-585,70,1}
		,{13682,10,1470,750,-585,70,1}
		,{13684,10,1470,750,-585,70,1}
		,{13685,10,1470,750,-585,70,1}
		,{13686,10,1470,750,-585,70,1}
		,{13687,10,1470,750,-585,70,1}
		,{13688,10,1470,750,-585,70,1}
		,{13689,10,1470,750,-585,70,1}
		,{13694,10,1470,750,-585,70,1}
		,{13695,10,1470,750,-585,70,1}
		,{13698,10,1470,750,-585,70,1}
		,{13699,10,1470,750,-585,70,1}
		,{13700,10,1470,750,-585,70,1}
		,{13710,10,1470,750,-585,70,1}
		,{13711,12,1470,750,-609,70,1}
		,{13712,12,1470,750,-609,70,1}
		,{13714,12,1470,750,-609,70,1}
		,{13715,12,1470,750,-609,70,1}
		,{13716,12,1470,750,-609,70,1}
		,{13717,12,1470,750,-609,70,1}
		,{13721,12,1470,750,-609,70,1}
		,{13723,12,1470,750,-609,70,1}
		,{13724,12,1470,750,-609,70,1}
		,{13727,12,1470,750,-609,70,1}
		,{13728,12,1470,750,-609,70,1}
		,{13729,12,1470,750,-609,70,1}
		,{13731,12,1470,750,-609,70,1}
		,{13732,12,1470,750,-609,70,1}
		,{13734,12,1470,750,-609,70,1}
		,{13736,12,1470,750,-609,70,1}
		,{13737,12,1470,750,-609,70,1}
		,{13746,12,1470,750,-609,70,1}
		,{13748,12,1470,750,-609,70,1}
		,{13749,12,1470,750,-609,70,1}
		,{13753,12,1470,750,-609,70,1}
		,{13754,12,1470,750,-609,70,1}
		,{13755,12,1470,750,-609,70,1}
		,{13756,12,1470,750,-609,70,1}
		,{13757,12,1470,750,-609,70,1}
		,{13758,12,1470,750,-609,70,1}
		,{13760,12,1470,750,-609,70,1}
		,{13761,12,1470,750,-609,70,1}
		,{13764,12,1470,750,-609,70,1}
		,{13765,12,1470,750,-609,70,1}
		,{13766,12,1470,750,-609,70,1}
		,{13767,12,1470,750,-609,70,1}
		,{13770,12,1470,750,-608,70,1}
		,{13771,12,1470,750,-607,70,1}
		,{13773,12,1470,750,-607,70,1}
		,{13775,12,1470,750,-607,70,1}
		,{13776,12,1470,750,-607,70,1}
		,{13777,12,1470,750,-607,70,1}
		,{13778,12,1470,750,-607,70,1}
		,{13779,12,1470,750,-607,70,1}
		,{13793,12,1470,750,-607,70,1}
		,{13797,12,1470,750,-607,70,1}
		,{13798,12,1470,750,-607,70,1}
		,{13799,12,1470,750,-607,70,1}
	};
	
	for( size_t i = 0; i < all_run_info.size(); i++){
		if( all_run_info[i][0] == runnum ){
			if( !(strncmp(lookup_var, "beam_current", 12)) ){
				return_var = all_run_info[i][1];
			}
			if( !(strncmp(lookup_var, "sbs_current", 11)) ){
				return_var = all_run_info[i][2];
			}
			if( !(strncmp(lookup_var, "bb_current", 10)) ){
				return_var = all_run_info[i][3];
			}
			if( !(strncmp(lookup_var, "bbcal_thresh", 12)) ){
				return_var = all_run_info[i][4];
			}
			if( !(strncmp(lookup_var, "sbs_field", 9)) ){
				return_var = 0.01*all_run_info[i][5];
			}
			if( !(strncmp(lookup_var, "target", 6)) ){
				return_var = all_run_info[i][6];
			}
		}
	}
	
	return return_var;
}
#endif
