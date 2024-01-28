#ifndef SETUP_SYSTEMATIC_VARIABLE_H
#define SETUP_SYSTEMATIC_VARIABLE_H

// Function to select/sort/print/etc systematic variables
// This function will re-define global variables so everything sent to the function should be within the scope

void SetupSystematicVariable(cl_syst& syst_vars){
	int kine = syst_vars.kine;

	//Should be pre-defined
	double val_min = syst_vars.val_min;
	double val_max = syst_vars.val_max;
	TString systematic_var_str = syst_vars.systematic_var_str;
	TString systematic_varVal_str = syst_vars.systematic_varVal_str;

	//To be set depending on which systematic variable is used
	//These should get default values though as only the systematic variable gets set, not all others. 

	double w2_mult_min, w2_mult_max, invMass_min, invMass_max;
	double SHPS_sigma_mult_min, SHPS_sigma_mult_max, Ep_min, PS_min, fcut_mult, HCal_clus_e_min;
	int gem_track_nhits_min;

	//Set Defaults by kinematic:
	if( kine == 4){
		w2_mult_min = 1.5;
		w2_mult_max = 1.5;
		SHPS_sigma_mult_min = 3.5;
		SHPS_sigma_mult_max = 3.5;
		Ep_min = 0.5;
		gem_track_nhits_min = 3;
		fcut_mult = 2;
		PS_min = 0.15;
		HCal_clus_e_min = 0.005;
	}
	if( kine == 8 ){
		w2_mult_min = 1.5;
		w2_mult_max = 1.5;
		SHPS_sigma_mult_min = 2.0;
		SHPS_sigma_mult_max = 2.0;
		Ep_min = 0.5;
		gem_track_nhits_min = 3;
		fcut_mult = 2;
		PS_min = 0.15;
		HCal_clus_e_min = 0.005;
	}
	if( kine == 9 ){
		w2_mult_min = 1.5;
		w2_mult_max = 1.5;
		SHPS_sigma_mult_min = 3.0;
		SHPS_sigma_mult_max = 3.0;
		Ep_min = 0.5;
		gem_track_nhits_min = 3;
		fcut_mult = 2;
		PS_min = 0.15;
		HCal_clus_e_min = 0.005;
	}
	//Options for systematic_varVal:
	// bb_gem_track_nhits, bb_ps_e, sbs_hcal_e, fcut_mult

	cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
	cout << "Systematic varVal str: " << systematic_varVal_str.Data() << endl;
	cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;


// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//					INVARIANT MASS
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
	//inv_mass_fixed_mean_vary_sigma_mult
	if( systematic_var_str == "inv_mass_fixed_mean_vary_sigma_mult" ){

		w2_mult_min = val_min;
		w2_mult_max = val_max;

		systematic_varVal_str = Form("%0.3f", w2_mult_min);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.w2_mult_min = w2_mult_min;
		syst_vars.w2_mult_max = w2_mult_max;		
		cout << "w2_mult_min: " << w2_mult_min << " ---  w2_mult_max " << w2_mult_max << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;	
	}

	//inv_mass_fixed_min_vary_max
	if( systematic_var_str == "inv_mass_fixed_min_vary_max" ){
		invMass_min = val_min;
		invMass_max = val_max;

		systematic_varVal_str = Form("%0.3f_to_%0.3f", invMass_min, invMass_max);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.invMass_min = invMass_min;
		syst_vars.invMass_max = invMass_max;

		cout << "invMass_min: " << invMass_min << " --- invMass_max: " << invMass_max << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
	}
	//inv_mass_fixed_moving_range
	if( systematic_var_str == "inv_mass_fixed_moving_range" ){
		invMass_min = val_min;
		invMass_max = val_max;

		systematic_varVal_str = Form("%0.3f_to_%0.3f", invMass_min, invMass_max);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.invMass_min = invMass_min;
		syst_vars.invMass_max = invMass_max;

		cout << "invMass_min: " << invMass_min << " --- invMass_max: " << invMass_max << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
	}


// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//					Shower and Pre-Shower
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
// VAR and DIR: SH_PS_sigma_mult
	if( systematic_var_str == "SH_PS_sigma_mult" ){
		SHPS_sigma_mult_min = val_min;
		SHPS_sigma_mult_max = val_max;

		systematic_varVal_str = Form("%0.3f", SHPS_sigma_mult_min);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.SHPS_sigma_mult_min = SHPS_sigma_mult_min;
		syst_vars.SHPS_sigma_mult_max = SHPS_sigma_mult_max;
		
		cout << "SHPS_sigma_mult_min: " << SHPS_sigma_mult_min << " ---  SHPS_sigma_mult_max " << SHPS_sigma_mult_max << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;	
	}
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//					E over P
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
// VAR and DIR: Ep_min
	if( systematic_var_str == "Ep_min" ){
		Ep_min = val_min;

		systematic_varVal_str = Form("%0.3f", Ep_min);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.Ep_min = Ep_min;
		
		cout << "Ep_min: " << Ep_min  << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;	
	}

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//					Pre-Shower
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// VAR and DIR: PS_min
	if( systematic_var_str == "PS_min" ){
		PS_min = val_min;

		systematic_varVal_str = Form("%0.3f", PS_min);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.PS_min = PS_min;
		
		cout << "PS_min: " << PS_min  << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;	
	}
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//					Fiducial Cut
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
// VAR and DIR: fcut_mult
	if( systematic_var_str == "fcut_mult" ){
		fcut_mult = val_min;

		systematic_varVal_str = Form("%0.3f", fcut_mult);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.fcut_mult = fcut_mult;
		
		cout << "fcut_mult: " << fcut_mult  << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;	
	}
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//					HCal Cluster Energy
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
// VAR and DIR: HCal_clus_e_min
	if( systematic_var_str == "HCal_clus_e_min" ){
		HCal_clus_e_min = val_min;

		systematic_varVal_str = Form("%0.3f", HCal_clus_e_min);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.HCal_clus_e_min = HCal_clus_e_min;
		
		cout << "HCal_clus_e_min: " << HCal_clus_e_min  << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;	
	}
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//					GEM tracks 
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// VAR and DIR: gem_track_nhits_min
	if( systematic_var_str == "gem_track_nhits_min" ){
		gem_track_nhits_min = val_min;

		systematic_varVal_str = Form("%i", gem_track_nhits_min);
		systematic_varVal_str.ReplaceAll(".", "");

		syst_vars.gem_track_nhits_min = gem_track_nhits_min;
		
		cout << "gem_track_nhits_min: " << gem_track_nhits_min  << endl;
		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;	
	}
//*********************************************************

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//			CRUCIAL!!!!
// MUST SET THE VARIABLES IN THE SYST_VAR CLASS WITH ALL NEW VALUES:

	syst_vars.systematic_varVal_str = systematic_varVal_str;

}
#endif