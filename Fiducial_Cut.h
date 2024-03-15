#ifndef FIDUCIAL_CUT_H
#define FIDUCIAL_CUT_H

#include "/work/halla/sbs/jboyd/include/experimental_constants.h"


bool fiducial_cut_expected_pos( double kine_dx_pn_max, double kine_dx_p_sigma, double kine_dx_n_sigma, double kine_dy_sigma, double fid_cut_multiplier, double expected_x, double expected_y, double kine_hcal_x_fmin, double kine_hcal_x_fmax, double kine_hcal_y_fmin, double kine_hcal_y_fmax ){
	
	bool passed_test;
	double sigma_multiplier = 0.5;

	// Let us start with a simple acceptance area cut using the fiducial area:

	// // X-dim
	// if( (expected_x > kine_hcal_x_fmax) || (expected_x < kine_hcal_x_fmin) ){
	// 	return false;
	// }


	// Min Edge (Top of HCal) (Direction that proton bends towards)
	// Take the expected position and subtract dn_pn_max from it.
	// If this falls off the edge (is less than x_fmin) then we cut it. 
	// double adjusted_min_x_position = expected_x - kine_dx_pn_max + fid_cut_multiplier*(kine_dx_n_sigma - kine_dx_p_sigma);
	// bool min_x_edge_exceeeded = adjusted_min_x_position < kine_hcal_x_fmin;

	fiducial_xi = kine_hcal_x_fmin + fid_cut_multiplier*kine_dx_pn_max;
	fiducial_xf = kine_hcal_x_fmax + fid_cut_multiplier*kine_dx_n_sigma;
	fiducial_yi = kine_hcal_y_fmin + sigma_multiplier*kine_dy_sigma;
	fiducial_yf = kine_hcal_y_fmax - sigma_multiplier*kine_dy_sigma;

	// We can set some max/min values for these so that the fiducial region doesn't go off of HCal
	// This can happen if given too large of multiplier values:
	if( fiducial_xi > kine_hcal_x_fmin ){
		fiducial_xi = kine_hcal_x_fmin;
	}
	if( fiducial_xf < kine_hcal_x_fmax ){
		fiducial_xf = kine_hcal_x_fmax;
	}
	if( fiducial_yi < kine_hcal_y_fmin ){
		fiducial_yi = kine_hcal_y_fmin;
	}
	if( fiducial_yf > kine_hcal_y_fmax ){
		fiducial_yf = kine_hcal_y_fmax;
	}

	// // Y-dim
	bool min_y_edge_exceeded = expected_y - fiducial_yi < 0;
	bool max_y_edge_exceeded = expected_y - fiducial_yf > 0;
	if( min_y_edge_exceeded || max_y_edge_exceeded ){
		return false;
	}

	// if( ((expected_y + sigma_multiplier*kine_dy_sigma) > kine_hcal_y_fmax) || ((expected_y - sigma_multiplier*kine_dy_sigma) < kine_hcal_y_fmin) ){
	// 	return false;
	// }

	// X-dim
	bool min_x_edge_exceeeded = expected_x - fiducial_xi < 0;
	// bool min_x_edge_exceeeded = (expected_x - kine_dx_pn_max - fid_cut_multiplier*kine_dx_p_sigma ) < kine_hcal_x_fmin;


	// Max Edge (Bottom of HCal) (Opposite direction of proton bending)
	// Here we only give a safety factor for the width of the neutron peak:
	// double adjusted_max_x_position = expected_x + fid_cut_multiplier*kine_dx_n_sigma;
	// bool max_x_edge_exceeeded = adjusted_max_x_position > kine_hcal_x_fmax;

	bool max_x_edge_exceeeded = expected_x - fiducial_xf > 0;
	// bool max_x_edge_exceeeded = (expected_x - fid_cut_multiplier*kine_dx_n_sigma ) > kine_hcal_x_fmax;

	if( min_x_edge_exceeeded || max_x_edge_exceeeded ){
		passed_test = false;
	}
	else{
		passed_test = true;
	}

	return passed_test;

}

bool fiducial_cut_HCal_x_pos( double kine_dx_pn_max, double kine_dx_p_sigma, double kine_dx_n_sigma, double fid_cut_multiplier, double hcal_hit_x, double hcal_hit_y, double kine_hcal_x_fmin, double kine_hcal_x_fmax, double kine_hcal_y_fmin, double kine_hcal_y_fmax ){

	bool passed_test;

	// Let us start with a simple acceptance area cut using the fiducial area:

	// X-dim
	if( (hcal_hit_x > kine_hcal_x_fmax) || (hcal_hit_x < kine_hcal_x_fmin) ){
		return false;
	}
	// Y-dim
	if( (hcal_hit_y > kine_hcal_y_fmax) || (hcal_hit_y < kine_hcal_y_fmin) ){
		return false;
	}

	// Min Edge (Top of HCal) (Direction that proton bends towards)
	// In worst case, it is a neutron and the proton missed. So, we require a proton to be able to land "above" it
	// Take the HCal hit position and subtract dn_pn_max from it.
	// If this falls off the edge (is less than x_fmin) then we cut it. 
	double adjusted_min_x_position = hcal_hit_x - fid_cut_multiplier*kine_dx_pn_max;
	bool min_x_edge_exceeeded = adjusted_min_x_position < kine_hcal_x_fmin;

	// bool min_x_edge_exceeeded = (hcal_hit_x - fid_cut_multiplier*kine_dx_pn_max) < kine_hcal_x_fmin;

	// Max Edge (Bottom of HCal) (Opposite direction of proton bending)
	// In worst case this was a bended up proton. So, the neutron missed.
	// Therefore, we need to account for the neutron having landed "below" it. 
	// So, we add the dx_pn_max distance and see if it falls off. If it does, then we cut it:

	double adjusted_max_x_position = hcal_hit_x + fid_cut_multiplier*kine_dx_pn_max;
	// double adjusted_max_x_position = expected_x + fid_cut_multiplier*kine_dx_n_sigma;

	bool max_x_edge_exceeeded = adjusted_max_x_position > kine_hcal_x_fmax;

	// bool max_x_edge_exceeeded = (hcal_hit_x + fid_cut_multiplier*kine_dx_pn_max) > kine_hcal_x_fmax;

	if( min_x_edge_exceeeded || max_x_edge_exceeeded ){
		passed_test = false;
	}
	else{
		passed_test = true;
	}

	return passed_test;
}

#endif