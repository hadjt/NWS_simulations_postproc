

### Dictionaries to name and rename variables, attributes etc.
long_name_dict = {}

long_name_dict["SST"] = "Sea Surface Temperature"
long_name_dict["NBT"] = "Near Bed Temperature"
long_name_dict["DFT"] = "Difference between Sea Surface and Near Bed Temperature"
long_name_dict["SSS"] = "Sea Surface Salinity"
long_name_dict["NBS"] = "Near Bed Salinity"
long_name_dict["DFS"] = "Difference between Sea Surface and Near Bed Salinity"
long_name_dict["SSH"] = "Sea Surface Height above Geoid"
long_name_dict["MLD"] = "Mixed Layer Depth using the Kara approach"
long_name_dict["PEA"] = "Potential Energy Anomaly"
long_name_dict["DMU"] = "Eastward Ocean Barotropic current"
long_name_dict["DMV"] = "Northward Ocean Barotropic current"
long_name_dict["DMUV"] = "Barotropic current speed on T grid"
long_name_dict["RegAveSST"] = "Regional Mean Sea Surface Temperature"
long_name_dict["RegAveNBT"] = "Regional Mean Near Bed Temperature"
long_name_dict[
    "RegAveDFT"
] = "Regional Mean Difference between Sea Surface and Near Bed Temperature"
long_name_dict["RegAveSSS"] = "Regional Mean Sea Surface Salinity"
long_name_dict["RegAveNBS"] = "Regional Mean Near Bed Salinity"
long_name_dict[
    "RegAveDFS"
] = "Regional Mean Difference between Sea Surface and Near Bed Salinity"
long_name_dict["RegAveSSH"] = "Regional Mean Sea Surface Height above Geoid"
long_name_dict["RegAvePEA"] = "Regional Mean Potential Energy Anomaly"
long_name_dict["reg_id"] = "Region ID number"
long_name_dict["cnt"] = "Regional count - number of grid boxes in region"
long_name_dict["mask"] = "Region Mask"

standard_name_dict = {}
standard_name_dict["SST"] = "sea_surface_temperature"
standard_name_dict["SSH"] = "sea_surface_height_above_geoid"
standard_name_dict["MLD"] = "ocean_mixed_layer_thickness"


unit_dict = {}
unit_dict["RegAveSST"] = "degC"
unit_dict["RegAveNBT"] = "degC"
unit_dict["RegAveDFT"] = "degC"
unit_dict["RegAveSSS"] = "1e-3"
unit_dict["RegAveNBS"] = "1e-3"
unit_dict["RegAveDFS"] = "1e-3"
unit_dict["RegAveSSH"] = "m"
unit_dict["RegAvePEA"] = "J/m3"
unit_dict["SST"] = "degC"
unit_dict["NBT"] = "degC"
unit_dict["DFT"] = "degC"
unit_dict["SSS"] = "1e-3"
unit_dict["NBS"] = "1e-3"
unit_dict["DFS"] = "1e-3"
unit_dict["SSH"] = "m"
unit_dict["MLD"] = "m"
unit_dict["PEA"] = "J/m3"
unit_dict["DMU"] = "m/s"
unit_dict["DMV"] = "m/s"
unit_dict["DMUV"] = "m/s"
unit_dict["reg_id"] = "1"
unit_dict["cnt"] = "1"
unit_dict["mask"] = "1"

orig_var_name_dict = {}
orig_var_name_dict["RegAveSST"] = "reg_sst_ave"
orig_var_name_dict["RegAveNBT"] = "reg_nbt_ave"
orig_var_name_dict["RegAveDFT"] = "reg_dft_ave"
orig_var_name_dict["RegAveSSS"] = "reg_sss_ave"
orig_var_name_dict["RegAveNBS"] = "reg_nbs_ave"
orig_var_name_dict["RegAveDFS"] = "reg_dfs_ave"
orig_var_name_dict["RegAveSSH"] = "reg_ssh_ave"
orig_var_name_dict["RegAvePEA"] = "reg_pea_ave"
orig_var_name_dict["SST"] = "votemper_top"
orig_var_name_dict["NBT"] = "votemper_bot"
orig_var_name_dict["DFT"] = "votemper_dif"
orig_var_name_dict["SSS"] = "vosaline_top"
orig_var_name_dict["NBS"] = "vosaline_bot"
orig_var_name_dict["DFS"] = "vosaline_dif"
orig_var_name_dict["SSH"] = "sossheig"
orig_var_name_dict["MLD"] = "mldkara"
orig_var_name_dict["PEA"] = "pea"
orig_var_name_dict["DMU"] = "vobtcrtx"
orig_var_name_dict["DMV"] = "vobtcrty"

ens_stat_long_name_format_dict = {}
# ens_stat_long_name_format_dict[var + ens_stat]% (var.upper(),perlab,long_name_dict[var])
ens_stat_long_name_format_dict["_ensmean"] = "Ensemble mean %s for %s: %s"
ens_stat_long_name_format_dict["_ensvar"] = "Ensemble variability %s for %s: %s"
ens_stat_long_name_format_dict["_intvar"] = "Interannual variability %s for %s: %s"
ens_stat_long_name_format_dict["_ensstd"] = "Ensemble standard deviation %s for %s: %s"

ens_stat_diff_long_name_format_dict = {}
# ens_stat_long_name_format_dict[var + ens_stat]% (var.upper(),perlab,long_name_dict[var])
ens_stat_diff_long_name_format_dict["_ensmean"] = "Change in ensemble mean %s for %s: %s"
ens_stat_diff_long_name_format_dict["_ensvar"] = "Change in ensemble variability %s for %s: %s"
ens_stat_diff_long_name_format_dict["_intvar"] = "Change in interannual variability %s for %s: %s"
ens_stat_diff_long_name_format_dict["_ensstd"] = "Change in ensemble standard deviation %s for %s: %s"
ens_stat_diff_long_name_format_dict["_projensmean"] = "Projected mean change in %s for %s: %s"
ens_stat_diff_long_name_format_dict["_projensstd"] = "Uncertainty (standard deviation) on projected change in %s for %s: %s"



ens_stat_cell_methods_dict = {}
ens_stat_cell_methods_dict[
    "_ensmean"
] = "time: mean within years time: mean over years realization: mean"
ens_stat_cell_methods_dict[
    "_ensvar"
] = "time: mean within years time: mean over years realization: variance"
ens_stat_cell_methods_dict[
    "_intvar"
] = "time: mean within years time: variance over years realization: mean"
ens_stat_cell_methods_dict[
    "_ensstd"
] = "time: mean within years time: mean over years realization: standard_deviation"
ens_stat_cell_methods_dict[
    "_projensmean"
] = "time: mean within years time: mean over years realization: mean"
ens_stat_cell_methods_dict[
    "_projensstd"
] = "time: mean within years time: mean over years realization: standard_deviation"




ens_stat_units_format_dict = {}
# ens_stat_units_format_dict[var + ens_stat]%unit_dict[var]
ens_stat_units_format_dict["_ensmean"] = "%s"
ens_stat_units_format_dict["_ensvar"] = "(%s)^2"
ens_stat_units_format_dict["_intvar"] = "(%s)^2"
ens_stat_units_format_dict["_ensstd"] = "%s"
ens_stat_units_format_dict["_projensmean"] = "%s"
ens_stat_units_format_dict["_projensstd"] = "%s"

var_dict = {}
var_dict["T"] = [
    "SST", "NBT", "DFT", "SSS", "NBS",
    "DFS", "SSH", "MLD", "PEA","DMUV"
    ]
var_dict["U"] = ["DMU"]
var_dict["V"] = ["DMV"]
var_dict["R"] = [
    "RegAveSST","RegAveNBT","RegAveDFT","RegAveSSS",
    "RegAveNBS","RegAveDFS","RegAveSSH","RegAvePEA",
    ]
