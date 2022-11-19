# include

Includes multiple files that can be included as header files. 

These headers will provide functions that can be used for retrieving experiment data/info such as GEM information and beam condtions, etc.

The two main "include" files are:

*GEM_lookup.h* --> Functions for retrieving GEM information
  lookup_config: returns Experiment Run Config number (0, 1, or 2)
  lookup_GEM_type: returns GEM type if given a runnum and GEM layer; 0 = UV, 1 = UVa XY, 2 = INFN XY
  lookup_GEM_type_from_global_mod_num: provide a runnum and a global module number (0 thru 11 depending on the config) and it turns a GEM Type (see above)
  lookup_nModules_total: returns the total number of GEM modules on the setup
  lookup_nModules_layer: returns the number of GEM modules on a given layer
  lookup_global_mod_num: Input a runnum, layer, and module on layer and you get back the global module number
  lookup_nAPVs: Provide a GEM Type and specificy U or V and this returns the number of APVs
  APV_strip_nums: Provide an APV adnd either "min" or "max". The return is either the first ("min") or last ("max") strip on that APV
  strip_to_APV: Input a strip and this returns what APV that corresponds to
  
  ---------------------------------------------
  GEM_channel_to_strip: This is a channel-to-strip mapping for GEM APVs. Input a gemType (0: UV, 1: UVa XY, 2: INFN XY) and a channel and get back the strip
  GEM_strip_channel: This is a strip-to-channel mapping for GEM APVs. Input a gemType (0: UV, 1: UVa XY, 2: INFN XY) and a strip and get back the channel
  
  THE FOLLOWING FUNCTIONS ARE INDIVIDIAULS SCRIPTS FOR STRIP<-->CHANNEL MAPPING
  THESE ARE ALL REFERENCED BY EITHER 'GEM_channel_to_strip' or 'GEM_strip_to_channel'
  UVa_UV_APV_strip_to_channel: strip-to-channel mapping for UVa UV APVs
  UVa_UV_APV_channel_to_straip: channel-to-strip mapping for UVa UV APVs
  UVa_XY_APV_strip_to_channel: strip-to-channel mapping for UVa XY APVs
  UVa_XY_APV_channel_to_strip: channel-to-strip mapping for UVa XY APVs
  INFN_XY_APV_strip_to_channel: strip-to-channel mapping for INFN XY APVs
  INFN_XY_APV_channel_to_strip: channel-to-strip mapping for INFN XY APVs
  ---------------------------------------------
  
  *beam_variables.h* --> Functions for retrieving beam info/variables
  
  This contains multiple functions that were created and used before a more general function was created. This one is:
  
  lookup_run_info
  
  Inputs: runnum and lookup_var
  lookup_var can be: 
    "beam_current"
    "sbs_current"
    "bb_current"
    "bbcal_thresh"
    "sbs_field"
    "target"
    
    The return value from this function will be the "lookup_var" requested. For target 0 = LH2 and 1 = LD2. 
    
 lookup_target can also be used and returns the target as a string of either "LH2" or "LD2"
