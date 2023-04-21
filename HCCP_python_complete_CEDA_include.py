
##################################################################################
########### Remove MO file structure info                               ##########
##################################################################################


PDCtrl_path_in = (
    "/scratch/hadjt/UKCP/tmp_cray_data/UKCP18_gc30_present_Ctrl_amm7_CO6_river_clim/"
)
HCCP_path_in_pattern = "/scratch/hadjt/HCCP_UKCP_PPE/Results/amm7_ensemble/ap977_ar095_au084/HCCP_CO9_ap977_ar095_au084_%s_02/"


# Region Mask
mask_file = "/data/cr1/hadjt/data/reffiles/nemo_CO6_CMEMS_NWS_RAN_v5_comb_mask.nc"
rootgrp = Dataset(mask_file, "r", format="NETCDF4")
mask = rootgrp.variables["mask"][:, :]
rootgrp.close()

file_output_freq_ann = True

if file_output_freq_ann:
    HCCP_output_dir = "/project/shelf_sea_projection/HCCP_UKCP_PPE/proc_files/CEDA/"
else:
    HCCP_output_dir = (
        "/project/shelf_sea_projection/HCCP_UKCP_PPE/proc_files/CEDA_python_mon/"
    )

if not os.path.exists(HCCP_output_dir):
    os.makedirs(HCCP_output_dir)


##################################################################################
##################################################################################

