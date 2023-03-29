from netCDF4 import Dataset
import numpy as np

import NWS_simulations_postproc_config as NWS_config

HCCP_path_in_pattern = NWS_config.HCCP_path_in_pattern

# Standardised Lat/Lon values for output variables

eg_file_lst = [
    HCCP_path_in_pattern % ("r001i1p00000") + "/19900101_Monthly2D_grid_%s.nc" % ss
    for ss in ["T", "U", "V"]
]
lon_str_lst, lat_str_lst = ["nav_lon", "nav_lon_grid_U", "nav_lon_grid_V"], [
    "nav_lat",
    "nav_lat_grid_U",
    "nav_lat_grid_V",
]

lonlat_dict = {}
#lonlat_fixed_dict = {}
lonlat_out_nc_dict = {}

for eg_file, lon_str, lat_str, grid_val in zip(
    eg_file_lst, lon_str_lst, lat_str_lst, ["T", "U", "V"]
):

    rootgrp = Dataset(eg_file, "r", format="NETCDF4")
    lonlat_dict[lon_str] = rootgrp.variables[lon_str][:, :]
    lonlat_dict[lat_str] = rootgrp.variables[lat_str][:, :]
    rootgrp.close()

    nav_lon_st = lonlat_dict[lon_str][250, :].copy()
    nav_lat_st = lonlat_dict[lat_str][:, 50].copy()

    nav_lon_ma_st = np.ma.masked_equal(lonlat_dict[lon_str].copy(), 0).mean(axis=0)
    nav_lat_ma_st = np.ma.masked_equal(lonlat_dict[lat_str].copy(), 0).mean(axis=1)

    if np.corrcoef(nav_lon_st, nav_lon_ma_st)[0, 1] < 0.99999:
        pdb.set_trace()

    lonlat_out_nc_dict["lon_" + grid_val] = nav_lon_st
    lonlat_out_nc_dict["lat_" + grid_val] = nav_lat_st
