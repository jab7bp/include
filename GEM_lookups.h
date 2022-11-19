#ifndef GEM_LOOKUPS_H
#define GEM_LOOKUPS_H

#include <algorithm>

int config;
int gem_type;
int nModules_total;
int nModules_layer;
int global_module_number;
int numAPVs_U;
int numAPVs_V;

int lookup_config(int runnum){
	if( runnum >= 11180 && runnum <= 12073 ){
		config = 0;
	}
	else if( runnum >= 12078 && runnum <= 13086 ){
		config = 1;
	}
	else if( runnum >= 13095 && runnum <= 13799 ){
		config = 2;
	}
	else{ config = -99999; }

	return config;
}

int lookup_GEM_type(int runnum, int layer){
	int config_num = lookup_config(runnum);

	if( config_num == 0 ){
		if( layer == 0 ){
			gem_type = 0;
		}
		if( layer == 1 ){
			gem_type = 2;
		}
		if( layer == 2 ){
			gem_type = 0;
		}
		if( layer == 3 ){
			gem_type = 2;
		}
		if( layer == 4 ){
			gem_type = 1;
		}
	}

	if( config_num == 1 ){
		if( layer == 0 ){
			gem_type = 0;
		}
		if( layer == 1 ){
			gem_type = 2;
		}
		if( layer == 2 ){
			gem_type = 0;
		}
		if( layer == 3 ){
			gem_type = 0;
		}
		if( layer == 4 ){
			gem_type = 1;
		}
	}

	if( config_num == 2 ){
		if( layer == 0 ){
			gem_type = 0;
		}
		if( layer == 1 ){
			gem_type = 0;
		}
		if( layer == 2 ){
			gem_type = 0;
		}
		if( layer == 3 ){
			gem_type = 0;
		}
		if( layer == 4 ){
			gem_type = 1;
		}
	}

	return gem_type;
}

int lookup_GEM_type_from_global_mod_num(int runnum, int mod){
	int config_num = lookup_config(runnum);

	if( config_num == 0 ){
		if( (mod == 0) || (mod == 4) ){
			gem_type = 0;
		}
		if( (mod == 8) || (mod == 9) || (mod == 10) || ( mod == 11 ) ){
			gem_type = 1;
		}
		if( (mod == 1) || (mod == 2) || (mod == 3) || (mod == 5) || (mod == 6) || (mod == 7) ){
			gem_type = 2;
		}
		if( mod > 11 ){
			gem_type = -99900;
		}
	} 
	if( config_num == 1 ){
		if( (mod == 0) || (mod == 4) ){
			gem_type = 0;
		}
		if( (mod == 6) || (mod == 7) || (mod == 8) || ( mod == 9 ) ){
			gem_type = 1;
		}
		if( (mod == 1) || (mod == 2) || (mod == 3) ){
			gem_type = 2;
		}
		if( mod > 9 ){
			gem_type = -99901;
		}
	} 
	if( config_num == 2 ){
		if( (mod == 0) || (mod == 1) || (mod == 2) || (mod == 4) ){
			gem_type = 0;
		}
		if( (mod == 4) || (mod == 5) || (mod == 6) || ( mod == 7 ) ){
			gem_type = 1;
		}
		if( mod > 7 ){
			gem_type = -99900;
		}
	} 
	return gem_type;
}

int lookup_nModules_total(int runnum){
	int config_num = lookup_config(runnum);

	if( config_num == 0 ){
		nModules_total = 12;
	}
	if( config_num == 1 ){
		nModules_total = 10;
	}
	if( config_num == 2 ){
		nModules_total = 8;
	}
	return nModules_total;
}

int lookup_nModules_layer(int runnum, int layer){
	int config_num = lookup_config(runnum);

	if( config_num == 0 ){
		if( layer == 0 ){
			nModules_layer = 1;
		}
		if( layer == 1 ){
			nModules_layer = 3;
		}
		if( layer == 2 ){
			nModules_layer = 1;
		}
		if( layer == 3 ){
			nModules_layer = 3;
		}
		if( layer == 4 ){
			nModules_layer = 4;
		}
	}
	if( config_num == 1 ){
		if( layer == 0 ){
			nModules_layer = 1;
		}
		if( layer == 1 ){
			nModules_layer = 3;
		}
		if( layer == 2 ){
			nModules_layer = 1;
		}
		if( layer == 3 ){
			nModules_layer = 1;
		}
		if( layer == 4 ){
			nModules_layer = 4;
		}
	}
	if( config_num == 2 ){
		if( layer == 0 ){
			nModules_layer = 1;
		}
		if( layer == 1 ){
			nModules_layer = 1;
		}
		if( layer == 2 ){
			nModules_layer = 1;
		}
		if( layer == 3 ){
			nModules_layer = 1;
		}
		if( layer == 4 ){
			nModules_layer = 4;
		}
	}
	return nModules_layer;
}

int lookup_global_mod_num(int runnum, int layer, int module_on_layer){
	int config_num = lookup_config(runnum);

	if( config_num == 0){
		if( layer == 0 ){
			global_module_number = 0;
			if( module_on_layer > 0 ){
				global_module_number = -99900;
			}
		}
		if( layer == 1 ){
			global_module_number = module_on_layer + 1;
			if( module_on_layer > 2 ){
				global_module_number = -99901;
			}
		}
		if( layer == 2 ){
			global_module_number = 4;
			if( module_on_layer > 0 ){
				global_module_number = -99902;
			}
		}
		if( layer == 3 ){
			global_module_number = module_on_layer + 5;
			if( module_on_layer > 2 ){
				global_module_number = -99903;
			}
		}
		if( layer == 4 ){
			global_module_number = module_on_layer + 8;
			if( module_on_layer > 3 ){
				global_module_number = -99904;
			}
		}
	}

	if( config_num == 1){
		if( layer == 0 ){
			global_module_number = 0;
			if( module_on_layer > 0 ){
				global_module_number = -99910;
			}
		}
		if( layer == 1 ){
			global_module_number = module_on_layer + 1;
			if( module_on_layer > 2 ){
				global_module_number = -99911;
			}
		}
		if( layer == 2 ){
			global_module_number = 4;
			if( module_on_layer > 0 ){
				global_module_number = -99912;
			}
		}
		if( layer == 3 ){
			global_module_number = 5;
			if( module_on_layer > 0 ){
				global_module_number = -99913;
			}
		}
		if( layer == 4 ){
			global_module_number = module_on_layer + 6;
			if( module_on_layer > 3 ){
				global_module_number = -99914;
			}
		}
	}

	if( config_num == 2){
		if( layer == 0 ){
			global_module_number = 0;
			if( module_on_layer > 0 ){
				global_module_number = -99920;
			}
		}
		if( layer == 1 ){
			global_module_number = 1;
			if( module_on_layer > 0 ){
				global_module_number = -99921;
			}
		}
		if( layer == 2 ){
			global_module_number = 2;
			if( module_on_layer > 0 ){
				global_module_number = -99922;
			}
		}
		if( layer == 3 ){
			global_module_number = 3;
			if( module_on_layer > 0 ){
				global_module_number = -99923;
			}
		}
		if( layer == 4 ){
			global_module_number = module_on_layer + 4;
			if( module_on_layer > 3 ){
				global_module_number = -99924;
			}
		}
	}

	return global_module_number;
}

int lookup_nAPVs(int gemType, int UorV){
	int nAPVs = -99996;
	if( gemType == 0 ){
		numAPVs_U = 30;
		numAPVs_V = 30;
	}
	if( gemType == 1 ){
		numAPVs_U = 10;
		numAPVs_V = 12;
	}
	if( gemType == 2 ){
		numAPVs_U = 10;
		numAPVs_V = 8;
	}

	if( UorV == 0 ){
		nAPVs = numAPVs_U;
	}
	if( UorV == 1){
		nAPVs = numAPVs_V;
	}

	return nAPVs;
}

double APV_strip_nums(int APV, TString minmax){
	int minmax_output = -91000;
	int apv_min_max_strips[30][2] =	{{0, 127},
									{128, 255},
									{256, 383},
									{384, 511},
									{512, 639},
									{640, 767},
									{768, 895},
									{896, 1023},
									{1024, 1151},
									{1152, 1279},
									{1280, 1407},
									{1408, 1535},
									{1536, 1663},
									{1664, 1791},
									{1792, 1919},
									{1920, 2047},
									{2048, 2175},
									{2176, 2303},
									{2304, 2431},
									{2432, 2559},
									{2560, 2687},
									{2688, 2815},
									{2816, 2943},
									{2944, 3071},
									{3072, 3199},
									{3200, 3327},
									{3328, 3455},
									{3456, 3583},
									{3584, 3711},
									{3712, 3839}};
	
	int apv_min = apv_min_max_strips[APV][0];
	int apv_max = apv_min_max_strips[APV][1];

	if(minmax == "min"){ minmax_output = apv_min; }
	else if (minmax == "max"){ minmax_output = apv_max; }
	
	return minmax_output;
}

int strip_to_APV(int strip){
	int APV;
	if( strip >= 0 && strip <= 127 ){ APV = 0; }
	if( strip >= 128 && strip <= 255 ){ APV = 1; }
	if( strip >= 256 && strip <= 383 ){ APV = 2; }
	if( strip >= 384 && strip <= 511 ){ APV = 3; }
	if( strip >= 512 && strip <= 639 ){ APV = 4; }
	if( strip >= 640 && strip <= 767 ){ APV = 5; }
	if( strip >= 768 && strip <= 895 ){ APV = 6; }
	if( strip >= 896 && strip <= 1023 ){ APV = 7; }
	if( strip >= 1024 && strip <= 1151 ){ APV = 8; }
	if( strip >= 1152 && strip <= 1279 ){ APV = 9; }
	if( strip >= 1280 && strip <= 1407 ){ APV = 10; }
	if( strip >= 1408 && strip <= 1535 ){ APV = 11; }
	if( strip >= 1536 && strip <= 1663 ){ APV = 12; }
	if( strip >= 1664 && strip <= 1791 ){ APV = 13; }
	if( strip >= 1792 && strip <= 1919 ){ APV = 14; }
	if( strip >= 1920 && strip <= 2047 ){ APV = 15; }
	if( strip >= 2048 && strip <= 2175 ){ APV = 16; }
	if( strip >= 2176 && strip <= 2303 ){ APV = 17; }
	if( strip >= 2304 && strip <= 2431 ){ APV = 18; }
	if( strip >= 2432 && strip <= 2559 ){ APV = 19; }
	if( strip >= 2560 && strip <= 2687 ){ APV = 20; }
	if( strip >= 2688 && strip <= 2815 ){ APV = 21; }
	if( strip >= 2816 && strip <= 2943 ){ APV = 22; }
	if( strip >= 2944 && strip <= 3071 ){ APV = 23; }
	if( strip >= 3072 && strip <= 3199 ){ APV = 24; }
	if( strip >= 3200 && strip <= 3327 ){ APV = 25; }
	if( strip >= 3328 && strip <= 3455 ){ APV = 26; }
	if( strip >= 3456 && strip <= 3583 ){ APV = 27; }
	if( strip >= 3584 && strip <= 3711 ){ APV = 28; }
	if( strip >= 3712 && strip <= 3839 ){ APV = 29; }
	return APV;
}

int UVa_UV_APV_strip_to_channel(int strip) {
	int UV_APV_channel;
	//Mapping for UVa UV GEM APV strip to channel
	const int _mapped_strip_uva_uv[128] = {
	     31,  15, 127, 111,  27,  11, 123, 107,  23,   7,
	    119, 103,  19,   3, 115,  99,  30,  14, 126, 110,
	     26,  10, 122, 106,  22,   6, 118, 102,  18,   2,
	    114,  98,  29,  13, 125, 109,  25,   9, 121, 105,
	     21,   5, 117, 101,  17,   1, 113,  97,  28,  12,
	    124, 108,  24,   8, 120, 104,  20,   4, 116, 100,
	     16,   0, 112,  96,  32,  48,  64,  80,  36,  52,
	     68,  84,  40,  56,  72,  88,  44,  60,  76,  92,
	     33,  49,  65,  81,  37,  53,  69,  85,  41,  57,
	     73,  89,  45,  61,  77,  93,  34,  50,  66,  82,
	     38,  54,  70,  86,  42,  58,  74,  90,  46,  62,
	     78,  94,  35,  51,  67,  83,  39,  55,  71,  87,
	     43,  59,  75,  91,  47,  63,  79,  95
	};

	if(strip >= 128){
		std::cout << std::endl << " WARNING!!! APVs only have strips from 0 - 127. Cannot have strips >= 128." << std::endl;
		UV_APV_channel = -99990;
	}
	else{
		auto chan = std::find(_mapped_strip_uva_uv, _mapped_strip_uva_uv + 128, strip);
		UV_APV_channel = std::distance(_mapped_strip_uva_uv, chan);
	}
	return UV_APV_channel;	
}

int UVa_UV_APV_channel_to_strip(int channel){
	int UVa_UV_APV_strip;
	//Mapping for UVa UV GEM APV strip to channel
	const int _mapped_strip_uva_uv[128] = {
	     31,  15, 127, 111,  27,  11, 123, 107,  23,   7,
	    119, 103,  19,   3, 115,  99,  30,  14, 126, 110,
	     26,  10, 122, 106,  22,   6, 118, 102,  18,   2,
	    114,  98,  29,  13, 125, 109,  25,   9, 121, 105,
	     21,   5, 117, 101,  17,   1, 113,  97,  28,  12,
	    124, 108,  24,   8, 120, 104,  20,   4, 116, 100,
	     16,   0, 112,  96,  32,  48,  64,  80,  36,  52,
	     68,  84,  40,  56,  72,  88,  44,  60,  76,  92,
	     33,  49,  65,  81,  37,  53,  69,  85,  41,  57,
	     73,  89,  45,  61,  77,  93,  34,  50,  66,  82,
	     38,  54,  70,  86,  42,  58,  74,  90,  46,  62,
	     78,  94,  35,  51,  67,  83,  39,  55,  71,  87,
	     43,  59,  75,  91,  47,  63,  79,  95
	};

	if(channel >= 128){
		std::cout << std::endl << " WARNING!!! APVs only have channels from 0 - 127. Cannot have channels >= 128." << std::endl;
		UVa_UV_APV_strip = -99991;
	}

	else{
		UVa_UV_APV_strip = _mapped_strip_uva_uv[channel];
	}
	return UVa_UV_APV_strip;
}

int UVa_XY_APV_strip_to_channel(int strip){
	int UVa_XY_APV_channel;
	//Mapping for UVa XY APV strip to channel
	const int _mapped_strip_uva_xy[128] = {
	     1,  33, 65,  97,  9,  41, 73, 105, 17,  49,
	    81, 113, 25,  57, 89, 121,  3,  35, 67,  99,
	    11,  43, 75, 107, 19,  51, 83, 115, 27,  59,
	    91, 123,  5,  37, 69, 101, 13,  45, 77, 109,
	    21,  53, 85, 117, 29,  61, 93, 125,  7,  39,
	    71, 103, 15,  47, 79, 111, 23,  55, 87, 119,
	    31,  63, 95, 127,  0,  32, 64,  96,  8,  40,
	    72, 104, 16,  48, 80, 112, 24,  56, 88, 120,
	     2,  34, 66,  98, 10,  42, 74, 106, 18,  50,
	    82, 114, 26,  58, 90, 122,  4,  36, 68, 100,
	    12,  44, 76, 108, 20,  52, 84, 116, 28,  60,
	    92, 124,  6,  38, 70, 102, 14,  46, 78, 110,
	    22,  54, 86, 118, 30,  62, 94, 126
	};

	if(strip >= 128){
		std::cout << std::endl << " WARNING!!! APVs only have strips from 0 - 127. Cannot have strips >= 128." << std::endl;
		UVa_XY_APV_channel = -99993;
	}
	else{
		auto chan = std::find(_mapped_strip_uva_xy, _mapped_strip_uva_xy + 128, strip);
		UVa_XY_APV_channel = std::distance(_mapped_strip_uva_xy, chan);
	}
	return UVa_XY_APV_channel;
}

int UVa_XY_APV_channel_to_strip(int channel){
	int UVa_XY_APV_strip;
	//Mapping for UVa XY GEM APV strip to channel
	const int _mapped_strip_uva_xy[128] = {
	     1,  33, 65,  97,  9,  41, 73, 105, 17,  49,
	    81, 113, 25,  57, 89, 121,  3,  35, 67,  99,
	    11,  43, 75, 107, 19,  51, 83, 115, 27,  59,
	    91, 123,  5,  37, 69, 101, 13,  45, 77, 109,
	    21,  53, 85, 117, 29,  61, 93, 125,  7,  39,
	    71, 103, 15,  47, 79, 111, 23,  55, 87, 119,
	    31,  63, 95, 127,  0,  32, 64,  96,  8,  40,
	    72, 104, 16,  48, 80, 112, 24,  56, 88, 120,
	     2,  34, 66,  98, 10,  42, 74, 106, 18,  50,
	    82, 114, 26,  58, 90, 122,  4,  36, 68, 100,
	    12,  44, 76, 108, 20,  52, 84, 116, 28,  60,
	    92, 124,  6,  38, 70, 102, 14,  46, 78, 110,
	    22,  54, 86, 118, 30,  62, 94, 126
	};

	if(channel >= 128){
		std::cout << std::endl << " WARNING!!! APVs only have channels from 0 - 127. Cannot have channels >= 128." << std::endl;
		UVa_XY_APV_strip = -99991;
	}

	else{
		UVa_XY_APV_strip = _mapped_strip_uva_xy[channel];
	}
	return UVa_XY_APV_strip;
}

int INFN_XY_APV_strip_to_channel(int strip){
	int INFN_XY_APV_channel;
	//Mapping for UVa XY APV strip to channel
	const int _mapped_strip_infn_xy[128] = {
	    0,  32,  64,  96,   8,  40,  72, 104,  16,  48,
	    80, 112,  24,  56,  88, 120,   1,  33,  65,  97,
	     9,  41,  73, 105,  17,  49,  81, 113,  25,  57,
	    89, 121,   2,  34,  66,  98,  10,  42,  74, 106,
	    18,  50,  82, 114,  26,  58,  90, 122,   3,  35,
	    67,  99,  11,  43,  75, 107,  19,  51,  83, 115,
	    27,  59,  91, 123,   4,  36,  68, 100,  12,  44,
	    76, 108,  20,  52,  84, 116,  28,  60,  92, 124,
	     5,  37,  69, 101,  13,  45,  77, 109,  21,  53,
	    85, 117,  29,  61,  93, 125,   6,  38,  70, 102,
	    14,  46,  78, 110,  22,  54,  86, 118,  30,  62,
	    94, 126,   7,  39,  71, 103,  15,  47,  79, 111,
	    23,  55,  87, 119,  31,  63,  95, 127
	};


	if(strip >= 128){
		std::cout << std::endl << " WARNING!!! APVs only have strips from 0 - 127. Cannot have strips >= 128." << std::endl;
		INFN_XY_APV_channel = -99993;
	}
	else{
		auto chan = std::find(_mapped_strip_infn_xy, _mapped_strip_infn_xy + 128, strip);
		INFN_XY_APV_channel = std::distance(_mapped_strip_infn_xy, chan);
	}
	return INFN_XY_APV_channel;
}

int INFN_XY_APV_channel_to_strip(int channel){
	int INFN_XY_APV_strip;
	//Mapping for UVa XY GEM APV strip to channel
	const int _mapped_strip_infn_xy[128] = {
	    0,  32,  64,  96,   8,  40,  72, 104,  16,  48,
	    80, 112,  24,  56,  88, 120,   1,  33,  65,  97,
	     9,  41,  73, 105,  17,  49,  81, 113,  25,  57,
	    89, 121,   2,  34,  66,  98,  10,  42,  74, 106,
	    18,  50,  82, 114,  26,  58,  90, 122,   3,  35,
	    67,  99,  11,  43,  75, 107,  19,  51,  83, 115,
	    27,  59,  91, 123,   4,  36,  68, 100,  12,  44,
	    76, 108,  20,  52,  84, 116,  28,  60,  92, 124,
	     5,  37,  69, 101,  13,  45,  77, 109,  21,  53,
	    85, 117,  29,  61,  93, 125,   6,  38,  70, 102,
	    14,  46,  78, 110,  22,  54,  86, 118,  30,  62,
	    94, 126,   7,  39,  71, 103,  15,  47,  79, 111,
	    23,  55,  87, 119,  31,  63,  95, 127
	};

	if(channel >= 128){
		std::cout << std::endl << " WARNING!!! APVs only have channels from 0 - 127. Cannot have channels >= 128." << std::endl;
		INFN_XY_APV_strip = -99991;
	}

	else{
		INFN_XY_APV_strip = _mapped_strip_infn_xy[channel];
	}
	return INFN_XY_APV_strip;
}

int GEM_channel_to_strip(int gemType, int channel){
	int lookup_strip = -99997;
	if( gemType == 0 ){
		lookup_strip = UVa_UV_APV_channel_to_strip(channel);
	}
	if( gemType == 1 ){
		lookup_strip = UVa_XY_APV_channel_to_strip(channel);
	}
	if( gemType == 2 ){
		lookup_strip = INFN_XY_APV_channel_to_strip(channel);
	}
	return lookup_strip;
}

int GEM_strip_to_channel(int gemType, int strip){
	int lookup_channel = -99998;
	if( gemType == 0 ){
		lookup_channel = UVa_UV_APV_strip_to_channel(strip);
	}
	if( gemType == 1 ){
		lookup_channel = UVa_XY_APV_strip_to_channel(strip);
	}
	if( gemType == 2 ){
		lookup_channel = INFN_XY_APV_strip_to_channel(strip);
	}
	return lookup_channel;
}

#endif

