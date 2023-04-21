from netCDF4 import Dataset
import os

# local directory (where code is stored)
localdir=os.path.dirname(os.path.realpath(__file__))

# Input data as output from NEMO - used for running run_CEDA_regional_means run_CEDA_monthly within the Met Office.
PDCtrl_path_in = localdir + "/"
NWSPPE_path_in_pattern = localdir + "/"

# Region Mask
# Download from github. Edit location of NWSPPE_mask_file if not stored with this code.
NWSPPE_mask_file = localdir + "/NWSClim_Wakelin_regmask.nc"
rootgrp = Dataset(NWSPPE_mask_file, "r", format="NETCDF4")
reg_mask_mat = rootgrp.variables["mask"][:,:]
rootgrp.close()

# Output data - Change to a local directory if reprocessing climatologies and ensemble statistics.
NWSPPE_output_dir = localdir + "/"
