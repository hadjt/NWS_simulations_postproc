# -*- coding: utf-8 -*-
from netCDF4 import Dataset,num2date
import pdb
import numpy as np
import matplotlib.pyplot as plt
import sys
from shutil import copyfile,move
from datetime import datetime,timedelta
import os
import glob

##################################################################################
##################################################################################
########### Remove MO file structure info                               ##########
##################################################################################
##################################################################################


PDCtrl_path_in = '/scratch/hadjt/UKCP/tmp_cray_data/UKCP18_gc30_present_Ctrl_amm7_CO6_river_clim/'
HCCP_path_in_pattern = '/scratch/hadjt/HCCP_UKCP_PPE/Results/amm7_ensemble/ap977_ar095_au084/HCCP_CO9_ap977_ar095_au084_%s_02/'


# Region Mask
mask_file = '/data/cr1/hadjt/data/reffiles/nemo_CO6_CMEMS_NWS_RAN_v5_comb_mask.nc'
rootgrp = Dataset(mask_file,'r', format='NETCDF4')
mask = rootgrp.variables['mask'][:,:]
rootgrp.close()

file_output_freq_ann = True

if file_output_freq_ann:
    HCCP_output_dir='/project/shelf_sea_projection/HCCP_UKCP_PPE/proc_files/CEDA/'
else:
    HCCP_output_dir='/project/shelf_sea_projection/HCCP_UKCP_PPE/proc_files/CEDA_python_mon/'

if not os.path.exists(HCCP_output_dir): os.makedirs(HCCP_output_dir)


##################################################################################
##################################################################################


HCCP_ens_mat_12 = ['r001i1p00000', 'r001i1p00605', 'r001i1p00834', 'r001i1p01113', 'r001i1p01554','r001i1p01649', 'r001i1p01843', 'r001i1p01935', 'r001i1p02123', 'r001i1p02242','r001i1p02491', 'r001i1p02868']  #

monstr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
monmxx = ['m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12']

reg_dim = 'region'
lat_dim = 'lat'
lon_dim = 'lon'
time_dim = 'time'
bnd_dim = 'bnds'

### Dictionaries to name and rename variables, attributes etc.
long_name_dict = {}

long_name_dict['SST'] = 'Sea Surface Temperature'
long_name_dict['NBT'] = 'Near Bed Temperature'
long_name_dict['DFT'] = 'Difference between Sea Surface and Near Bed Temperature'
long_name_dict['SSS'] = 'Sea Surface Salinity'
long_name_dict['NBS'] = 'Near Bed Salinity'
long_name_dict['DFS'] = 'Difference between Sea Surface and Near Bed Salinity'
long_name_dict['SSH'] = 'Sea Surface Height above Geoid'
long_name_dict['MLD'] = 'Mixed Layer Depth using the Kara approach'
long_name_dict['PEA'] = 'Potential Energy Anomaly'
long_name_dict['DMU'] = 'Eastward Ocean Barotropic current'
long_name_dict['DMV'] = 'Northward Ocean Barotropic current'
long_name_dict['DMUV'] = 'Barotropic current speed on T grid'
long_name_dict['RegAveSST'] = 'Regional Mean Sea Surface Temperature'
long_name_dict['RegAveNBT'] = 'Regional Mean Near Bed Temperature'
long_name_dict['RegAveDFT'] = 'Regional Mean Difference between Sea Surface and Near Bed Temperature'
long_name_dict['RegAveSSS'] = 'Regional Mean Sea Surface Salinity'
long_name_dict['RegAveNBS'] = 'Regional Mean Near Bed Salinity'
long_name_dict['RegAveDFS'] = 'Regional Mean Difference between Sea Surface and Near Bed Salinity'
long_name_dict['RegAveSSH'] = 'Regional Mean Sea Surface Height above Geoid'
long_name_dict['RegAvePEA'] = 'Regional Mean Potential Energy Anomaly'
long_name_dict['reg_id'] = 'Region ID number'
long_name_dict['cnt'] = 'Regional count - number of grid boxes in region'
long_name_dict['mask'] = 'Region Mask'

standard_name_dict = {}
standard_name_dict['SST'] = 'sea_surface_temperature'
standard_name_dict['SSH'] = 'sea_surface_height_above_geoid'
standard_name_dict['MLD'] = 'ocean_mixed_layer_thickness'

unit_dict = {}
unit_dict['RegAveSST'] = 'degC'
unit_dict['RegAveNBT'] = 'degC'
unit_dict['RegAveDFT'] = 'degC'
unit_dict['RegAveSSS'] = '1e-3'
unit_dict['RegAveNBS'] = '1e-3'
unit_dict['RegAveDFS'] = '1e-3'
unit_dict['RegAveSSH'] = 'm'
unit_dict['RegAvePEA'] = 'J/m3'
unit_dict['SST'] = 'degC'
unit_dict['NBT'] = 'degC'
unit_dict['DFT'] = 'degC'
unit_dict['SSS'] = '1e-3'
unit_dict['NBS'] = '1e-3'
unit_dict['DFS'] = '1e-3'
unit_dict['SSH'] = 'm'
unit_dict['MLD'] = 'm'
unit_dict['PEA'] = 'J/m3'
unit_dict['DMU'] = 'm/s'
unit_dict['DMV'] = 'm/s'
unit_dict['DMUV'] = 'm/s'
unit_dict['reg_id'] = '1'
unit_dict['cnt'] = '1'
unit_dict['mask'] = '1'

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
ens_stat_long_name_format_dict['_ensmean'] = 'Ensemble mean %s for %s: %s'
ens_stat_long_name_format_dict['_ensvar'] = 'Ensemble variability %s for %s: %s'
ens_stat_long_name_format_dict['_intvar'] = 'Interannual variability %s for %s: %s'
ens_stat_long_name_format_dict['_ensstd'] = 'Ensemble standard deviation %s for %s: %s'

ens_stat_cell_methods_dict = {}
ens_stat_cell_methods_dict['_ensmean'] = 'time: mean within years time: mean over years realization: mean'
ens_stat_cell_methods_dict['_ensvar'] = 'time: mean within years time: mean over years realization: variance'
ens_stat_cell_methods_dict['_intvar'] = 'time: mean within years time: variance over years realization: mean'
ens_stat_cell_methods_dict['_ensstd'] = 'time: mean within years time: mean over years realization: standard_deviation'

ens_stat_units_format_dict = {}
#ens_stat_units_format_dict[var + ens_stat]%unit_dict[var]
ens_stat_units_format_dict['_ensmean'] = '%s'
ens_stat_units_format_dict['_ensvar'] = "(%s)^2"
ens_stat_units_format_dict['_intvar'] = "(%s)^2"
ens_stat_units_format_dict['_ensstd'] = '%s'

add_DMUV = True

var_dict = {}
var_dict['T'] = ['SST', 'NBT', 'DFT', 'SSS', 'NBS', 'DFS', 'SSH', 'MLD', 'PEA']
if add_DMUV: var_dict['T'] = ['SST', 'NBT', 'DFT', 'SSS', 'NBS', 'DFS', 'SSH', 'MLD', 'PEA','DMUV']
var_dict['U'] = ['DMU']
var_dict['V'] = ['DMV']
var_dict['R'] = ['RegAveSST','RegAveNBT','RegAveDFT','RegAveSSS','RegAveNBS','RegAveDFS','RegAveSSH','RegAvePEA']

# Read in and calculate a standard lon and lat array.
grid_list = ['T','U','V']

# Standardised Lat/Lon values for output variables

eg_file_lst = [HCCP_path_in_pattern%('r001i1p00000') + '/19900101_Monthly2D_grid_%s.nc'%ss for ss in grid_list]
lon_str_lst,lat_str_lst = ['nav_lon','nav_lon_grid_U','nav_lon_grid_V'],['nav_lat','nav_lat_grid_U','nav_lat_grid_V']

lonlat_dict = {}
lonlat_fixed_dict = {}
lonlat_out_nc_dict = {}

for eg_file,lon_str,lat_str, grid_val in zip(eg_file_lst,lon_str_lst,lat_str_lst,['T','U','V']):

    rootgrp = Dataset(eg_file,'r', format='NETCDF4')
    lonlat_dict[lon_str] = rootgrp.variables[lon_str][:,:]
    lonlat_dict[lat_str] = rootgrp.variables[lat_str][:,:]
    rootgrp.close()

    nav_lon_st = lonlat_dict[lon_str][250,:].copy()
    nav_lat_st = lonlat_dict[lat_str][:,50].copy()

    nav_lon_ma_st = np.ma.masked_equal(lonlat_dict[lon_str].copy(),0).mean(axis = 0)
    nav_lat_ma_st =  np.ma.masked_equal(lonlat_dict[lat_str].copy(),0).mean(axis = 1)

    if np.corrcoef(nav_lon_st , nav_lon_ma_st)[0,1]<0.99999:
        pdb.set_trace()

    lonlat_fixed_dict[lon_str] = lonlat_dict[lon_str].copy()
    lonlat_fixed_dict[lat_str] = lonlat_dict[lat_str].copy()
    lonlat_fixed_dict[lon_str][:,:] = nav_lon_st
    lonlat_fixed_dict[lat_str].T[:,:] = nav_lat_st.T

    lonlat_out_nc_dict['lon_' + grid_val] = nav_lon_st
    lonlat_out_nc_dict['lat_' + grid_val] = nav_lat_st

lonlat_keys_lst = lonlat_fixed_dict.keys()

out_Conventions = "CF-1.8" # "CF-1.6"

out_institution = "Met Office Hadley Centre, Exeter, UK."

out_title = "Marine climate projections for the North West European Shelf Seas.\nA GCM Perturbed Parameter Ensemble and present day control simulation downscaled with a shelf seas model.\n"

out_history = 'Model output post processed by NWS_simulations_postproc.py'

out_references = 'Tinker et al. 2023, A set of climate projections for the NW European Shelf Seas, in prep.'

out_comment_monthly_mean = 'These are monthly mean values calculated from every timestep, as output by NEMO shelf version %s.'

out_comment_regmean = "Regional mean time series are output by NEMO shelf version %s, on the Wakelin et al. (2012) region mask (see region_refererence). NEMO reads in the region mask, and every time step it averages the variable within each mask region. The regional mean is then averaged over the month, and output by NEMO. This file includes the regional mean time-series for a selection of variables (as RegAveSST, RegAveSSS etc. for the regional mean SST and SSS respectively), for each month, and each of the 14 regions. There are variables for the region id (reg_id) and the number of grid boxes within each region (cnt). The region mask included (mask, with the associated longitude and latitude variables). This methodology is described by Tinker et al. (2019; see region_methology_reference)."



out_region_names = "['Shelf','Southern North Sea','Central North Sea','Northern North Sea','English Channel','Skagerrak/Kattegat','Norwegian Trench','Shetland Shelf','Irish Shelf','Irish Sea','Celtic Sea','Armorican Shelf','NE Atlantic (S)','NE Atlantic (N)']"

out_region_notes = 'The region Shelf is the combination of Southern North Sea; Central North Sea; Northern North Sea; English Channel; Shetland Shelf; Irish Shelf; Irish Sea; and Celtic Sea. It does not include the Norwegian Trench; Skagerrak/Kattegat; Armorican Shelf; NE Atlantic (S); NE Atlantic (N)'

out_region_refererence = 'Adapted from Wakelin, S. L., Holt, J., Blackford, J., Allen, I., Butenschön, M., and Artioli, Y.: Modeling the carbon fluxes of the northwest European continental shelf: Validation and budgets, 117, C05020, https://doi.org/10.1029/2011JC007402, 2012.'.encode('utf8')  #

out_region_methology_reference = 'Tinker, J., Renshaw, R., Barciela, R., and Wood, R.: Regional mean time series for the Northwest European Shelf seas. In: Copernicus Marine Service Ocean State Report, Issue 3, J. Oper. Oceanogr., 12, s26-s30, https://doi.org/10.1080/1755876X.2019.1633075, 2019.'

out_source_PDCtrl = "Underlying GCM: HadGEM3 GC3.\nDownscaling shelf seas model: NEMO shelf version 3.6.\nPresent Day Control Simulation representing the year 2000.\nModel simulations developed and run by Dr. Jonathan Tinker"

out_source_NWSPPE = "Underlying GCM: HadGEM3 GC3.\nDownscaling shelf seas model: NEMO shelf version 4.0.4.\nRepresentative Concentration Pathway: RCP8.5.\nModel simulations developed and run by Dr. Jonathan Tinker"



out_comment_grid = 'Note that NEMO uses the Arakawa "C" grid, where the T, U and V grids are offset. Most variables are on the T grid, apart from the Eastward and Northward components of the ocean barotropic current, which are on the U and V grids respectively. We also provide the barotropic current speed on T grid, where we transform the U and V velocity components onto the T grid, before calculating their magnitude. We therefore separate the variables on the T, U and V grids into separate files. This file contains variables on the %s grid.'


out_comment_EnsStat = 'These ensemble statistics are calculated from the climatological mean and standard deviation for each NWSPPE ensemble member, which in turn are based on the monthly mean values calculated from every timestep, as output by NEMO shelf version 4.0.4. For a given climatological period (e.g. 2000-2019, 2079-2098), and month, season or year, the climatological mean and standard deviation of all 12 NWSPPE ensmeble members are loaded. The mean of these climatological means gives the Ensemble Mean (ensmean). The variance of the climatological means gives the Ensemble Variance (ensvar), and its square-root gives the Ensemble Standard Deviation (ensstd). The mean of climatological standard deviation squared (i.e. the climatological variance) gives the Interannual Variance (intvar). Each variable (SST, SSS, NBT, etc.) has 4 statistic associated with it, and the NetCDF variable names join the variable name (e.g. SST) with the statistic name (e.g. ensmean), separated by an underscore. For example the ensmeble mean SST is SST_ensmean and the NBT interannual variability is NBT_intvar. This is all captured in the "cell_methods", where "realization" denotes the ensemble members. The Ensemble Variance and Interannual Variance can be combined into the Total Variance following Tinker et al. 2016; see variance_separation_references).'


out_variance_separation_references = 'Tinker, J., Lowe, J., Pardaens, A., Holt, J., and Barciela, R.: Uncertainty in climate projections for the 21st century northwest European shelf seas, Prog. Oceanogr., 148, 56-73, https://doi.org/10.1016/j.pocean.2016.09.003, 2016.'











# compress variables? smaller, but file size varies, so spot incomplete files
compress_variables = False

def create_ann_seas_lst(yrmat):
    """
    Create a list of the months and years that make up annual and seasonal means, for a given list of years

    Jonathan Tinker 29/03/2023
    """
    ann_lst = []
    djf_lst = []
    jja_lst = []
    mam_lst = []
    son_lst = []
    for yr in yrmat:

        tmp_ann_lst = []
        tmp_djf_lst = []
        tmp_jja_lst = []
        tmp_mam_lst = []
        tmp_son_lst = []
        for mi in range(12):tmp_ann_lst.append('%04i%02i'%(yr,mi+1))
        for mi in [2,3,4]:tmp_mam_lst.append('%04i%02i'%(yr,mi+1))
        for mi in [5,6,7]:tmp_jja_lst.append('%04i%02i'%(yr,mi+1))
        for mi in [8,9,10]:tmp_son_lst.append('%04i%02i'%(yr,mi+1))
        tmp_djf_lst.append('%04i12'%(yr-1))
        for mi in [0,1]:tmp_djf_lst.append('%04i%02i'%(yr,mi+1))

        ann_lst.append(tmp_ann_lst)
        djf_lst.append(tmp_djf_lst)
        mam_lst.append(tmp_mam_lst)
        jja_lst.append(tmp_jja_lst)
        son_lst.append(tmp_son_lst)

    return ann_lst,djf_lst,mam_lst,jja_lst,son_lst

def create_mon_lst(yrmat):
    """    
    Create a list of the months and years that for each month, for a given list of years

    Jonathan Tinker 29/03/2023
    """
    mon_lst = []
    for mi in range(12):
        tmp_mon_lst = []
        for yr in yrmat:tmp_mon_lst.append(['%04i%02i'%(yr,mi+1)])
        mon_lst.append(tmp_mon_lst)

    return mon_lst


def run_CEDA_regional_means():
    """
    Read in the regional mean files output from NEMO COx and create a single regional mean file for CEDA, with standardised and corrected meta data.

    Jonathan Tinker 29/03/2023
    """
    grid_val = 'R'
    mnmat = np.arange(1,12+1)

    # Create the directory if missing.
    if not os.path.exists(HCCP_output_dir):
        os.makedirs(HCCP_output_dir)

    #cycle through the ensemble members
    for ens in ['PDCtrl']+ HCCP_ens_mat_12:
        print(ens, datetime.now())

        # input path and years for NWSPPE and PDCtrl differ.
        if ens =='PDCtrl':

            path_in = PDCtrl_path_in

            yrmat = np.arange(2050,2250+1)
            fname_ens = ens
            file_out= 'NWSClim_PDCtrl_2050-2250_regmean.nc'
            path_out = HCCP_output_dir + 'NWSClim/%s/regmean/'%ens
        else:
            path_in = HCCP_path_in_pattern%(ens)

            yrmat = np.arange(1990,2098+1)
            fname_ens = 'NWSPPE_%s'%ens
            file_out= 'NWSClim_%s_1990-2098_regmean.nc'%fname_ens
            path_out = HCCP_output_dir + 'NWSClim/NWSPPE/%s/regmean/'%ens

        if not os.path.exists(path_out):
            os.makedirs(path_out)


        # Open output file, create variables and attributes.
        rootgrp_out = Dataset(path_out +file_out, 'w', format='NETCDF4')

        region_dim_nc = rootgrp_out.createDimension(reg_dim, 14)
        lat_dim_nc = rootgrp_out.createDimension(lat_dim, 375)
        lon_dim_nc = rootgrp_out.createDimension(lon_dim, 297)
        time_counter_dim_nc = rootgrp_out.createDimension(time_dim, None)
        axis_nbounds_dim_nc = rootgrp_out.createDimension(bnd_dim, 2)

        # create time variables, and add attributes
        time_var = rootgrp_out.createVariable('time', 'f8',(time_dim))
        timebnds_var = rootgrp_out.createVariable('time_bounds', 'f8',(time_dim,bnd_dim))
        time_var.setncattr('long_name',"Time axis")
        time_var.setncattr('calendar',"360_day")
        time_var.setncattr('units',"seconds since 1950-01-01 00:00:00")
        time_var.setncattr('time_origin',"1950-01-01 00:00:00")
        time_var.setncattr('bounds',"time_bounds")

        tmpvardict={}
        for var in var_dict[grid_val]:
            tmpvardict[var] = rootgrp_out.createVariable(var, 'f4', (time_dim, reg_dim),fill_value = 1.e+20,zlib  = compress_variables)

        for var in ['reg_id','cnt']:
            tmpvardict[var] = rootgrp_out.createVariable(var, 'f4', ( reg_dim),fill_value = 1.e+20,zlib  = compress_variables)

        for var in ['mask']:
            tmpvardict[var] = rootgrp_out.createVariable(var, 'f4', ( lat_dim,lon_dim),fill_value = 1.e+20,zlib  = compress_variables)

        for var in var_dict[grid_val] + ['reg_id','cnt']:
            tmpvardict[var].setncattr('long_name',long_name_dict[var])
            if var in standard_name_dict.keys():                           tmpvardict[var].setncattr('standard_name',standard_name_dict[var])
            tmpvardict[var].setncattr('units',unit_dict[var])

        for var in var_dict[grid_val] :
            tmpvardict[var].setncattr('online_operation',"average")
            tmpvardict[var].setncattr('interval_write',"1 month")
            tmpvardict[var].setncattr('interval_operation',"1 time-step")
            tmpvardict[var].setncattr('cell_methods',"Area: mean time: mean")
        # Add mask
        for var in ['mask']:
            tmpvardict[var].setncattr('long_name',long_name_dict[var])
            tmpvardict[var].setncattr('units',unit_dict[var])

        #Add lon and lats
        lon_var = rootgrp_out.createVariable('lon', 'f4',(lon_dim))
        lat_var = rootgrp_out.createVariable('lat', 'f4',(lat_dim))

        lon_var.setncattr('long_name',"Longitude")
        lon_var.setncattr('standard_name',"longitude")
        lon_var.setncattr('units',"degrees_east")
        lon_var.setncattr('nav_model',"grid_T")
        lat_var.setncattr('long_name',"Latitude")
        lat_var.setncattr('standard_name',"latitude")
        lat_var.setncattr('units',"degrees_north")
        lat_var.setncattr('nav_model',"grid_T")

        lon_var[:] = lonlat_out_nc_dict['lon_T']
        lat_var[:] = lonlat_out_nc_dict['lat_T']

        tmpvardict['mask'][:,:] = mask[2]

        # Add Global Attributes

        if ens == 'PDCtrl':
            out_source = out_source_PDCtrl
            tmp_nemo_version = '3.6'
        else:
            out_source = out_source_NWSPPE
            tmp_nemo_version = '4.0.4'


        rootgrp_out.setncattr('Conventions',out_Conventions)
        rootgrp_out.setncattr('institution',out_institution)
        rootgrp_out.setncattr('title',out_title)
        rootgrp_out.setncattr('source',out_source)

        rootgrp_out.setncattr('region_names', out_region_names)
        rootgrp_out.setncattr('region_notes', out_region_notes)
        rootgrp_out.setncattr('region_refererence', out_region_refererence)
        rootgrp_out.setncattr('region_methology_reference', out_region_methology_reference)

        rootgrp_out.setncattr('history',out_history)
        rootgrp_out.setncattr('references',out_references)
        rootgrp_out.setncattr('comments',out_comment_regmean%tmp_nemo_version)


        # loop throuhg year and month with a counter.
        counter = 0
        for yr in yrmat:
            if (yr %10 )==0 :print(yr, datetime.now())
            for mn in mnmat:

                #Open files
                file_in = glob.glob('%s%04i%02i01_Region_mean_timeseries_Monthly*.nc'%(path_in,yr,mn))
                rootgrp_in = Dataset(file_in[0], 'r', format='NETCDF4')

                # time: Check origin, 1950 or 1900
                orig_time_counter_origin = rootgrp_in.variables['time_counter'].time_origin
                if orig_time_counter_origin == "1900-01-01 00:00:00":
                    calendar_time_offset = (50.*360.*86400.)
                elif orig_time_counter_origin == "1950-01-01 00:00:00":
                    calendar_time_offset = 0.
                else:
                    print('time_counter origin needs checking')

                #PDCtrl overflows when second counter gets over 2**32, so work around.
                if ens =='PDCtrl':

                    tmp_calc_time_counter = (yr - 1950.)*360.*86400. + (mn-1.)*30.*86400 + 15.*86400
                    tmp_calc_time_counter_bounds = np.zeros((2))
                    tmp_calc_time_counter_bounds[0] = (yr - 1950.)*360.*86400. + (mn-1.)*30.*86400.
                    tmp_calc_time_counter_bounds[1] = (yr - 1950.)*360.*86400. + (mn-1.)*30.*86400. + (30.*86400.)

                    # test work around.
                    # confirm that the time calculation matches that from NEMO, when accounting from the 2**32 overflow
                    if not ((tmp_calc_time_counter_bounds % (2**32))== ((rootgrp_in.variables['time_counter_bounds'][:,:] - calendar_time_offset)% (2**32))).all():
                        print('Calculated time bounds (instantaneous values) not offset by 2**32 - check')
                        pdb.set_trace()

                    time_var[counter:counter+1] = tmp_calc_time_counter
                    timebnds_var[counter:counter+1,:] = tmp_calc_time_counter_bounds
                else:
                    time_var[counter:counter+1] = rootgrp_in.variables['time_counter'][:] - calendar_time_offset
                    timebnds_var[counter:counter+1,:] = rootgrp_in.variables['time_counter_bounds'][:,:] - calendar_time_offset

                # check which regions are the Shelf, and Wakelin regions.
                tmp_cnt =  rootgrp_in.variables['cnt'][:].ravel()
                tmp_reg_id =  rootgrp_in.variables['reg_id'][:].ravel()
                tmp_mask_id =  rootgrp_in.variables['mask_id'][:].ravel()
                #pdb.set_trace()

                reg_ind =  ((tmp_mask_id == 2)  & (tmp_reg_id == 1))| ((tmp_mask_id == 3) &  (tmp_reg_id > 0))

                tmpvardict['cnt'][:] = tmp_cnt[reg_ind]
                tmpvardict['reg_id'][:] = tmp_reg_id[reg_ind]

                # fill the variables with data.
                for var in var_dict[grid_val]:
                    tmp_orig_var_name_dict = orig_var_name_dict[var]

                    tmpvardict[var][counter:counter+1,:] = rootgrp_in.variables[tmp_orig_var_name_dict][:].ravel()[reg_ind]

                counter+=1
                rootgrp_in.close()

        rootgrp_out.close()


def run_CEDA_monthly(file_output_freq_ann = True):
    """
    read in the monthly mean files output from NEMO COx and recreate them for CEDA, with standardised and corrected meta data.

    Jonathan Tinker 29/03/2023
    """
    mnmat = np.arange(1,12+1)

    # Create the directory if missing.
    if not os.path.exists(HCCP_output_dir):
        os.makedirs(HCCP_output_dir)


    #cycle through the ensemble members
    for ens in ['PDCtrl'] + HCCP_ens_mat_12:
        print(ens, datetime.now())

        # input path and years for NWSPPE and PDCtrl differ.
        if ens =='PDCtrl':
            path_in = PDCtrl_path_in

            yrmat = np.arange(2050,2250+1)
            fname_ens = ens
            path_out = HCCP_output_dir + 'NWSClim/%s/annual/'%ens
        else:
            path_in = HCCP_path_in_pattern%(ens)
            fname_ens = 'NWSPPE_%s'%ens

            yrmat = np.arange(1990,2098+1)
            path_out = HCCP_output_dir + 'NWSClim/NWSPPE/%s/annual/'%ens

        #output path, create it missing.
        if not os.path.exists(path_out):
            os.makedirs(path_out)

        #Cycle through years, months and grid (T, U, V)
        for yr in yrmat:
            if (yr %10 )==0 :print(yr, datetime.now())
            for mi,mn in enumerate(mnmat):

                # input file names
                tmp_file_in_dict = {}
                for grid_val in ['T','U','V']:
                    tmp_file_in_dict[grid_val] = glob.glob('%s%04i%02i01_Monthly2D_grid_%s*.nc'%(path_in,yr,mn,grid_val))[0]

                for grid_val in ['T','U','V']:

                    #  output and input file names
                    if file_output_freq_ann == True:
                        fname_date = '%04i'%(yr)
                        mi_ind = mi
                    else:
                        fname_date = '%04i%02i'%(yr,mn)
                        mi_ind = 0

                    file_out = 'NWSClim_%s_%s_grid%s.nc'%(fname_ens,fname_date,grid_val)
                    file_in = tmp_file_in_dict[grid_val]

                    # Open input and output file names, and then open them
                    rootgrp_in = Dataset(file_in, 'r', format='NETCDF4')


                    # Write file with basic dimensions, variables, and attributes, then close.
                    # then reopen with append, and add data - this allows annual files of 12 months to be easily written

                    # if first month of the year, or if monthly output write file.

                    if (file_output_freq_ann == False)| (file_output_freq_ann & (mi==0)):

                        rootgrp_out = Dataset(path_out +file_out, 'w', format='NETCDF4')

                        # Add dimensions, variables and attributes to output files.
                        lat_dim_nc = rootgrp_out.createDimension(lat_dim, 375)
                        lon_dim_nc = rootgrp_out.createDimension(lon_dim, 297)
                        time_counter_dim_nc = rootgrp_out.createDimension(time_dim, None)
                        axis_nbounds_dim_nc = rootgrp_out.createDimension(bnd_dim, 2)

                        time_var = rootgrp_out.createVariable('time', 'f8',(time_dim))
                        timebnds_var = rootgrp_out.createVariable('time_bounds', 'f8',(time_dim,bnd_dim))

                        #add output variable, and use a dictionary as the handle.
                        tmpvardict={}
                        for var in var_dict[grid_val]: tmpvardict[var] = rootgrp_out.createVariable(var, 'f4', (time_dim, lat_dim, lon_dim),fill_value = 1.e+20,zlib  = compress_variables)

                        # add time variables
                        time_var.setncattr('long_name',"Time axis")
                        time_var.setncattr('calendar',"360_day")
                        time_var.setncattr('units',"seconds since 1950-01-01 00:00:00")
                        time_var.setncattr('time_origin',"1950-01-01 00:00:00")
                        time_var.setncattr('bounds',"time_bounds")

                        # add lon/lat variables
                        lon_var = rootgrp_out.createVariable('lon', 'f4',(lon_dim))
                        lat_var = rootgrp_out.createVariable('lat', 'f4',(lat_dim))

                        # add attributes to lon/lat
                        lon_var.setncattr('long_name',"Longitude")
                        lon_var.setncattr('standard_name',"longitude")
                        lon_var.setncattr('units',"degrees_east")
                        lon_var.setncattr('nav_model',"grid_%s"%grid_val)
                        lat_var.setncattr('long_name',"Latitude")
                        lat_var.setncattr('standard_name',"latitude")
                        lat_var.setncattr('units',"degrees_north")
                        lat_var.setncattr('nav_model',"grid_%s"%grid_val)

                        #Add attributes to variables
                        for var in var_dict[grid_val]:
                            tmpvardict[var].setncattr('long_name',long_name_dict[var])
                            if var in standard_name_dict.keys():                           tmpvardict[var].setncattr('standard_name',standard_name_dict[var])
                            tmpvardict[var].setncattr('units',unit_dict[var])
                            #tmpvardict[var].setncattr('coordinates',"ensemble")
                            tmpvardict[var].setncattr('online_operation',"average")
                            tmpvardict[var].setncattr('interval_write',"1 month")
                            tmpvardict[var].setncattr('interval_operation',"1 time-step")
                            tmpvardict[var].setncattr('cell_methods',"time: mean")


                        if ens == 'PDCtrl':
                            out_source = out_source_PDCtrl
                            tmp_nemo_version = '3.6'
                        else:
                            out_source = out_source_NWSPPE
                            tmp_nemo_version = '4.0.4'

                        rootgrp_out.setncattr('Conventions',out_Conventions)
                        rootgrp_out.setncattr('institution',out_institution)
                        rootgrp_out.setncattr('title',out_title)
                        rootgrp_out.setncattr('source',out_source)

                        rootgrp_out.setncattr('history',out_history)
                        rootgrp_out.setncattr('references',out_references)
                        rootgrp_out.setncattr('comments',out_comment_monthly_mean%(tmp_nemo_version))
                        rootgrp_out.setncattr('comments_grid',out_comment_grid%(grid_val))

                        # Add lon/lat
                        lon_var[:] = lonlat_out_nc_dict['lon_' + grid_val]
                        lat_var[:] = lonlat_out_nc_dict['lat_' + grid_val]

                        rootgrp_out.close()


                    if add_DMUV and grid_val == 'T':
                        rootgrp_U_in = Dataset(tmp_file_in_dict['U'], 'r', format='NETCDF4')
                        rootgrp_V_in = Dataset(tmp_file_in_dict['V'], 'r', format='NETCDF4')

                        tmp_DMU = rootgrp_U_in.variables[orig_var_name_dict["DMU"]][:]
                        tmp_DMV = rootgrp_V_in.variables[orig_var_name_dict["DMV"]][:]
                        DMU_T = np.ma.zeros((1, 375, 297))*np.ma.masked
                        DMV_T = np.ma.zeros((1, 375, 297))*np.ma.masked
                        DMUV_T = np.ma.zeros((1, 375, 297))*np.ma.masked

                        # Interpolate DMU and DMV from U and V grid, to T grid.
                        # Confirmed with Jeff Polton
                        DMU_T[0, 1:, 1:] = (tmp_DMU[0, 1:,    :-1] + tmp_DMU[0, 1:, 1:] )/2.
                        DMV_T[0, 1:, 1:] = (tmp_DMV[0,  :-1, 1:  ] + tmp_DMV[0, 1:, 1:] )/2.

                        # Calculate the barotropic current speed
                        DMUV_T = np.sqrt(DMU_T**2 + DMV_T**2)

                        rootgrp_U_in.close()
                        rootgrp_V_in.close()

                    #pdb.set_trace()
                    rootgrp_out = Dataset(path_out +file_out, 'a', format='NETCDF4')

                    # Add dimensions, variables and attributes to output files.
                    time_var = rootgrp_out.variables['time']
                    timebnds_var = rootgrp_out.variables['time_bounds']

                    #add output variable, and use a dictionary as the handle.
                    tmpvardict={}
                    for var in var_dict[grid_val]: tmpvardict[var] = rootgrp_out.variables[var]
                    # add attributes to time.

                    # PDCtrl and NWSPPE have different time origins (1900 and 1950)
                    # Also PDCtrl has a time overflow problem, when seconds are > 2**32.
                    orig_time_counter_origin = rootgrp_in.variables['time_counter'].time_origin
                    if orig_time_counter_origin == "1900-01-01 00:00:00":
                        calendar_time_offset = (50.*360.*86400.)
                    elif orig_time_counter_origin == "1950-01-01 00:00:00":
                        calendar_time_offset = 0.
                    else:
                        print('time_counter origin needs checking')

                    # Implement a time overflow work around for when PDCtrl seconds are > 2**32.
                    if ens =='PDCtrl':
                        tmp_calc_time_counter = (yr - 1950.)*360.*86400. + (mn-1.)*30.*86400 + 15.*86400
                        tmp_calc_time_counter_bounds = np.zeros((2))
                        tmp_calc_time_counter_bounds[0] = (yr - 1950.)*360.*86400. + (mn-1.)*30.*86400.
                        tmp_calc_time_counter_bounds[1] = (yr - 1950.)*360.*86400. + (mn-1.)*30.*86400. + (30.*86400.)


                        if not ((tmp_calc_time_counter_bounds % (2**32))== ((rootgrp_in.variables['time_counter_bounds'][:,:] - calendar_time_offset)% (2**32))).all():
                            print('Calculated time bounds (instantaneous values) not offset by 2**32 - check')
                            pdb.set_trace()()


                        time_var[mi_ind] = tmp_calc_time_counter
                        timebnds_var[mi_ind,:] = tmp_calc_time_counter_bounds

                    else:
                        time_var[mi_ind] = rootgrp_in.variables['time_counter'][:] - calendar_time_offset
                        timebnds_var[mi_ind,:] = rootgrp_in.variables['time_counter_bounds'][:,:] - calendar_time_offset



                    # Loop through the variables, finding the NEMO original nc var name, and outputing with the new standardised names
                    for var in var_dict[grid_val]:

                        if var == 'DMUV': continue
                        tmp_orig_var_name_dict = orig_var_name_dict[var]
                        #PDCtrl and NWSPPE have different nc var names for MLD, so work around.
                        if (ens == 'PDCtrl'):
                            if (var == 'MLD'):
                                tmp_orig_var_name_dict = 'Mixed_Layer_Depth_Kara_et_al_definition'

                        tmpvardict[var][mi_ind,:,:] = rootgrp_in.variables[tmp_orig_var_name_dict][:]

                    #pdb.set_trace()

                    if add_DMUV and grid_val == 'T':

                        # add SST mask, so setting 15 additional sea points to land
                        DMUV_T.mask = tmpvardict['SST'][:].mask|DMUV_T.mask
                        tmpvardict['DMUV'][mi_ind,:,:] = DMUV_T

                    rootgrp_in.close()
                    rootgrp_out.close()






def run_CEDA_ens_climatologies(skip_existing = False, Test = False, Testvar = None, file_output_freq_ann = True):
    """
    Read in the CEDA reprocessed monthly mean files, and create climatological means and standard deviations, 
    for a near present day period (2000-2019) and a end of century period (2079-2098). This is repeated for all ensemble member for the NWSPPE, 
    and is produced for monthly, seasonal and annual means. The climatologiccal means and standard deviations are in different files, and the variable names are unchanged 
    (although the processing and statisitic is recorded in the variable attributes, cell_methods)
            
    Jonathan Tinker 29/03/2023
    """



    import warnings
    warnings.filterwarnings("ignore")

    for yrmat in [np.arange(2000,2019+1),np.arange(2079,2098+1)]:

        yrstr = '%04i-%04i'%(yrmat[0],yrmat[-1])

        ann_lst,djf_lst,mam_lst,jja_lst,son_lst = create_ann_seas_lst(yrmat)
        mon_lst = create_mon_lst(yrmat)

        date_name_a_s_mat = ['ann','djf','mam','jja','son']
        date_name_m_mat =[ii.lower() for ii in monmxx]
        date_name_mat = date_name_m_mat + date_name_a_s_mat

        for ens in HCCP_ens_mat_12:


            input_dir = HCCP_output_dir + 'NWSClim/NWSPPE/%s/annual/'%ens
            path_out = HCCP_output_dir + 'NWSClim/NWSPPE/%s/clim/'%ens
            if not os.path.exists(path_out): os.makedirs(path_out)

            #for date_lst,date_name in zip([ann_lst,djf_lst,mam_lst,jja_lst,son_lst],date_name_mat):
            for date_lst,date_name in zip(mon_lst + [ann_lst,djf_lst,mam_lst,jja_lst,son_lst],date_name_mat):

                #How many months (files) in each period.
                if date_name.capitalize() in monstr:
                    tmp_nfiles = 1
                    monthnamestr = date_name
                elif date_name in monmxx:
                    tmp_nfiles = 1
                    monthnamestr = np.array(monstr)[date_name == np.array(monmxx)][0]
                elif date_name.upper() in ['DJF', 'MAM', 'JJA', 'SON']:
                    tmp_nfiles = 3
                elif date_name.upper() in ['ANN']:
                    tmp_nfiles = 12


                for grid_val in ['T','U','V']:


                    output_file_pref = 'NWSClim_NWSPPE_%s_clim_%s_%04i-%04i_grid%s'%(ens,date_name,yrmat[0],yrmat[-1],grid_val)

                    tmpfname_clim = path_out+ output_file_pref + '_mean.nc'
                    tmpfname_clim_stdev = path_out+output_file_pref + '_stddev.nc'


                    print(datetime.now(),tmpfname_clim )
                    file_pref = 'NWSClim_NWSPPE_%s_'%(ens)
                    file_suff = '_grid%s.nc'%(grid_val)

                    if skip_existing:
                        if os.path.isfile(tmpfname_clim) & os.path.isfile(tmpfname_clim_stdev):
                            return

                    time_bnd = {}
                    if Test: test_lst = []

                    for yi,tmp_date_lst in enumerate(date_lst): # year within climatological period

                        # cycle through the months within period for the climatological period
                        int_cnt = 0
                        for td,tmp_date in enumerate(tmp_date_lst): #

                            if file_output_freq_ann:

                                tmp_fname = input_dir+file_pref+tmp_date[:4]+file_suff
                                mi = int(tmp_date[4:])-1

                            else:
                                tmp_fname = input_dir+file_pref+tmp_date+file_suff
                                mi = 0
                            ###############################################
                            if not os.path.isfile(tmp_fname):
                                continue
                            ###############################################

                            rootgrp = Dataset(tmp_fname, 'r', format='NETCDF4')
                            var = rootgrp.variables
                            if td ==  0:
                                tmp_key_mat = list(var.keys())
                                tmp_key_mat.remove('lon')
                                tmp_key_mat.remove('lat')

                                int_var_mat = {}
                                for tmp_key in tmp_key_mat:
                                    int_var_mat[tmp_key] = var[tmp_key][mi:mi+1].astype('double')
                                int_cnt=+1
                            else:
                                for tmp_key in tmp_key_mat:
                                    int_var_mat[tmp_key][:] += var[tmp_key][mi:mi+1].astype('double')
                                int_cnt += 1

                            # identify time and time bound var names, and take a copy of the time bounds.
                                # later iterations, take the min and max of the time.
                            if (td == 0) & (yi == 0):
                                tbnd_key_lst = [ss for ss in tmp_key_mat if (ss[-7:] == '_bounds') & (ss[:4] == 'time')]
                                time_key_lst = [ss for ss in tmp_key_mat if (ss[-7:] != '_bounds') & (ss[:4] == 'time')]
                                for tmp_tbnd_key in tbnd_key_lst:  time_bnd[tmp_tbnd_key] = var[tmp_tbnd_key][mi:mi+1,:]
                            else:

                                for tmp_tbnd_key in tbnd_key_lst:
                                    time_bnd[tmp_tbnd_key][0,0] = np.minimum(time_bnd[tmp_tbnd_key][0,0],var[tmp_tbnd_key][mi,0])
                                    time_bnd[tmp_tbnd_key][0,1] = np.maximum(time_bnd[tmp_tbnd_key][0,1],var[tmp_tbnd_key][mi,1])

                            rootgrp.close()

                        #if didn't count any files, stop
                        if int_cnt == 0:
                            print('no files this season')
                            pdb.set_trace()

                        #if didn't count any files, stop
                        if int_cnt != tmp_nfiles:
                            print('not enough files for ',date_name ,'Expected ',tmp_nfiles, 'got', int_cnt)
                            pdb.set_trace()

                        # internal variable average matrix
                            #monthly, seasonal or annual mean.

                        int_var_ave_mat = {}
                        for tmp_key in tmp_key_mat:
                            int_var_ave_mat[tmp_key] = int_var_mat[tmp_key][:]/int_cnt
                        if Test: test_lst.append(int_var_ave_mat[Testvar][0,120,120])

                        # if the first year, take a copy, and a copy of the square
                        if yi ==  0:
                            var_mat = {}
                            var_mat_stdev = {}
                            for tmp_key in tmp_key_mat:
                                var_mat[tmp_key] = int_var_ave_mat[tmp_key][:]
                                var_mat_stdev[tmp_key] = int_var_ave_mat[tmp_key][:]**2
                            tot_cnt = 1

                        # otherwise increment the sum, and the sum of the square, and the counter
                        else:
                            for tmp_key in tmp_key_mat:
                                var_mat[tmp_key][:] += int_var_ave_mat[tmp_key][:]
                                var_mat_stdev[tmp_key][:] += int_var_ave_mat[tmp_key][:]**2
                            tot_cnt += 1




                    if Test:
                        for tmp_key in tmp_key_mat:
                            tmp_mean=var_mat[Testvar][:]/tot_cnt
                            tmp_std=np.sqrt(var_mat_stdev[Testvar][:]/tot_cnt - tmp_mean**2)
                            print (tmp_mean[0,120,120], np.array(test_lst).mean())
                            print (tmp_std[0,120,120], np.array(test_lst).std())

                            pdb.set_trace()


                    for tmp_fill_fname in  [tmpfname_clim,tmpfname_clim_stdev]:

                        rootgrp_out = Dataset(tmp_fill_fname, 'w', format='NETCDF4_CLASSIC')
                        lat_dim_nc = rootgrp_out.createDimension(lat_dim, 375)
                        lon_dim_nc = rootgrp_out.createDimension(lon_dim, 297)
                        time_counter_dim_nc = rootgrp_out.createDimension(time_dim, None)
                        axis_nbounds_dim_nc = rootgrp_out.createDimension(bnd_dim, 2)

                        time_var = rootgrp_out.createVariable('time', 'f8',(time_dim))
                        timebnds_var = rootgrp_out.createVariable('time_bounds', 'f8',(time_dim,bnd_dim))

                        #add output variable, and use a dictionary as the handle.
                        tmpvardict={}
                        for var in var_dict[grid_val]: tmpvardict[var] = rootgrp_out.createVariable(var, 'f4', (time_dim, lat_dim, lon_dim),fill_value = 1.e+20,zlib  = compress_variables)

                        # add time variables
                        time_var.setncattr('long_name',"Time axis")
                        time_var.setncattr('calendar',"360_day")
                        time_var.setncattr('units',"seconds since 1950-01-01 00:00:00")
                        time_var.setncattr('time_origin',"1950-01-01 00:00:00")
                        time_var.setncattr('bounds',"time_bounds")

                        # add lon/lat variables
                        lon_var = rootgrp_out.createVariable('lon', 'f4',(lon_dim))
                        lat_var = rootgrp_out.createVariable('lat', 'f4',(lat_dim))

                        # add attributes to lon/lat
                        lon_var.setncattr('long_name',"Longitude")
                        lon_var.setncattr('standard_name',"longitude")
                        lon_var.setncattr('units',"degrees_east")
                        lon_var.setncattr('nav_model',"grid_%s"%grid_val)
                        lat_var.setncattr('long_name',"Latitude")
                        lat_var.setncattr('standard_name',"latitude")
                        lat_var.setncattr('units',"degrees_north")
                        lat_var.setncattr('nav_model',"grid_%s"%grid_val)


                        #Add attributes to variables
                        for var in var_dict[grid_val]:
                            tmpvardict[var].setncattr('long_name',long_name_dict[var])
                            if var in standard_name_dict.keys():                           tmpvardict[var].setncattr('standard_name',standard_name_dict[var])
                            tmpvardict[var].setncattr('units',unit_dict[var])
                            tmpvardict[var].setncattr('online_operation',"average")
                            tmpvardict[var].setncattr('interval_write',"1 month")
                            tmpvardict[var].setncattr('interval_operation',"1 time-step")
                            tmpvardict[var].setncattr('cell_methods',"time: mean")




                        out_source = out_source_NWSPPE

                        rootgrp_out.setncattr('Conventions',out_Conventions)
                        rootgrp_out.setncattr('institution',out_institution)
                        rootgrp_out.setncattr('title',out_title)
                        rootgrp_out.setncattr('source',out_source)


                        rootgrp_out.setncattr('history',out_history)
                        rootgrp_out.setncattr('references',out_references)
                        # Add lon/lat
                        lon_var[:] = lonlat_out_nc_dict['lon_' + grid_val]
                        lat_var[:] = lonlat_out_nc_dict['lat_' + grid_val]


                        rootgrp_out.close()

                    #pdb.set_trace()
                    rootgrp_out = Dataset(tmpfname_clim, 'a', format='NETCDF4_CLASSIC')
                    rootgrp_stdev_out = Dataset(tmpfname_clim_stdev, 'a', format='NETCDF4_CLASSIC')


                    # cycle through the variables, calcuate the mean and standard devation, and write them out.
                    for tmp_key in tmp_key_mat:
                        tmp_mean=var_mat[tmp_key][:]/tot_cnt
                        tmp_std=np.sqrt(var_mat_stdev[tmp_key][:]/tot_cnt - tmp_mean**2)
                        #tmp_variance=(var_mat_stdev[tmp_key][:]/tot_cnt - tmp_mean**2)
                        out_var = rootgrp_out.variables[tmp_key]
                        out_var[:] = tmp_mean
                        out_stdev_clim = rootgrp_stdev_out.variables[tmp_key]
                        out_stdev_clim[:] = tmp_std

                        #update the cell_methods
                        if tmp_key not in ['time', 'time_bounds', 'lon', 'lat']:
                            out_var.setncattr('cell_methods','time: mean within years time: mean over years')
                            out_stdev_clim.setncattr('cell_methods','time: mean within years time: standard_deviation over years')

                    # update the time bounds, and the time to reflect the climatologies
                    for tmp_tbnd_key in tbnd_key_lst:
                        rootgrp_out.variables[tmp_tbnd_key][:] = time_bnd[tmp_tbnd_key][:]
                        rootgrp_stdev_out.variables[tmp_tbnd_key][:] = time_bnd[tmp_tbnd_key][:]

                    for tmp_time_key in time_key_lst:
                        rootgrp_stdev_out.variables[tmp_time_key][:] = rootgrp_out.variables[tmp_time_key][:]

                    #update the comment global attribute
                    if (date_name.capitalize() in monstr) | (date_name in monmxx):
                        out_comment_mean = "These are climatological means for each ensemble member. The %s monthly means (as output by NEMO shelf version 4.0.4.) within the current climatological period (%s) are selected, and their mean is taken to give these climatological means."%(monthnamestr,yrstr)

                        out_comment_stdev = "These are climatological standard deviations for each ensemble member. The %s monthly means (as output by NEMO shelf version 4.0.4.) within the current climatological period (%s) are selected, and their standard deviation is taken to give these climatological standard deviation."%(monthnamestr,yrstr)

                    elif date_name.upper() in ['DJF', 'MAM', 'JJA', 'SON']:

                        out_comment_mean = "These are climatological means for each ensemble member. The monthly means (as output by NEMO shelf version 4.0.4.) are averaged into seasonal means for %s. The seasonal means within the current climatological period (%s) are selected, and their mean is taken to give these climatological means." %(date_name.upper(),yrstr)

                        out_comment_stdev = "These are climatological standard deviations for each ensemble member. The monthly means (as output by NEMO shelf version 4.0.4.) are averaged into seasonal means for %s. The seasonal means within the current climatological period (%s) are selected, and their standard deviation is taken to give these climatological standard deviation." %(date_name.upper(),yrstr)

                    elif date_name.upper() in ['ANN']:

                        out_comment_mean = "These are climatological means for each ensemble member. The monthly means (as output by NEMO shelf version 4.0.4.) are averaged into annual means. The annual means within the current climatological period (%s) are selected, and their mean is taken to give these climatological means."%(yrstr)

                        out_comment_stdev = "These are climatological standard deviations for each ensemble member. The monthly means (as output by NEMO shelf version 4.0.4.) are averaged into annual means. The annual means within the current climatological period (%s) are selected, and their standard deviation is taken to give these climatological standard deviation."%(yrstr)
                    else:
                        print("check date_name",date_name)
                        #pdb.set_trace()

                    rootgrp_out.setncattr('comments',out_comment_mean)
                    rootgrp_stdev_out.setncattr('comments',out_comment_stdev)

                    rootgrp_out.setncattr('comments_grid',out_comment_grid%(grid_val))
                    rootgrp_stdev_out.setncattr('comments_grid',out_comment_grid%(grid_val))


                    rootgrp_out.close()
                    rootgrp_stdev_out.close()


def run_CEDA_ens_stats( Test = False, Testvar = None,yrmat_1 = np.arange(2000,2019+1),yrmat_2 = np.arange(2079,2098+1)):
    """
    Read in the CEDA climatological means and standard deviations, for ensemble member for the NWSPPE, 
    for the near present day period (2000-2019) and end of century period (2079-2098), for the monthly, seasonal and annual means, 
    and calculate ensemble statistics. Ensemble statistics include ensmeble mean,  ensmeble standard deviations,  ensmeble variance and interannual variance, 
    for the present day, end of century, and the difference between them. All statistics are in the same file, with the variable name modified with a suffix to indentify the statistic (and the processing and statistic is recorded in the variable attributes, cell_methods)

    Jonathan Tinker 29/03/2023
    """

    print ('start: ',datetime.now())


    path_out = '%s/NWSClim/EnsStats/'%HCCP_output_dir

    #if the output dir is missing, create it
    if not os.path.exists(path_out): os.makedirs(path_out)

    import warnings
    warnings.filterwarnings("ignore")

    date_name_mat = ['ann','djf','mam','jja','son']+ [ii.lower() for ii in monmxx]

    ens_stat_lst = ['_ensmean','_ensvar','_intvar','_ensstd']


    # Labeling output date strings for filesnames
    yrmat_1_str = '%04i-%04i'%(yrmat_1[0],yrmat_1[-1])
    yrmat_2_str = '%04i-%04i'%(yrmat_2[0],yrmat_2[-1])
    diff_yrmat_str = '%sminus%s'%(yrmat_2_str, yrmat_1_str)

    for grid_val in ['T','U','V',]:

        var_mat = var_dict[grid_val]
        #ncvar_mat = var_dict[grid_val]


        # create processing and projections dictionaries
        # proc dict to increment the sum and the sum of the squares
        proc_dat = {}
        proj_dat = {}
        proj_dat['yrmat_1'] = yrmat_1
        proj_dat['yrmat_2'] = yrmat_2

        #other dictionary for lon, lat, time etc.
        other_dat = {}

        for seas in date_name_mat:
            proc_dat[seas] = {}
            proj_dat[seas] = {}
            other_dat[seas] = {}
            #for var, ncvar in zip(var_mat,ncvar_mat):
            for var in var_mat:
                proc_dat[seas][var] = {}
                proj_dat[seas][var] = {}


            if Test:
                Test_mean = {}
                Test_std = {}

        # loop though date types (months, seasons, annuals)
        #   open files for present day, and future climatologies (means and std devs)
        #   sum and sum squares of the variables.
        for di,seas in enumerate(date_name_mat):
            proc_cnt = 0

            if Test:
                Test_mean[seas] = np.ma.zeros(12)
                Test_std[seas] = np.ma.zeros(12)


            print (datetime.now(), ' Loading data clim_%s_%s_grid%s'%( seas,yrmat_1_str,grid_val))
            for ei,ens in enumerate(HCCP_ens_mat_12):
                input_dir = HCCP_output_dir + 'NWSClim/NWSPPE/%s/clim/'%ens

                # file names
                tmpfname_clim_1 =  ('%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_mean.nc'%(input_dir, ens,seas,yrmat_1_str,grid_val))
                tmpfname_clim_2 =  ('%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_mean.nc'%(input_dir, ens,seas,yrmat_2_str,grid_val))
                tmpfname_clim_stdev_1 =  ('%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_stddev.nc'%(input_dir, ens,seas,yrmat_1_str,grid_val))
                tmpfname_clim_stdev_2 =  ('%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_stddev.nc'%(input_dir, ens,seas,yrmat_2_str,grid_val))


                #open files
                rootgrp_clim_1 = Dataset(tmpfname_clim_1, 'r', format='NETCDF4')
                rootgrp_clim_2 = Dataset(tmpfname_clim_2, 'r', format='NETCDF4')
                rootgrp_clim_stdev_1 = Dataset(tmpfname_clim_stdev_1, 'r', format='NETCDF4')
                rootgrp_clim_stdev_2 = Dataset(tmpfname_clim_stdev_2, 'r', format='NETCDF4')

                #find variables
                var_clim_1 = rootgrp_clim_1.variables
                var_clim_stdev_1 = rootgrp_clim_stdev_1.variables
                var_clim_2 = rootgrp_clim_2.variables
                var_clim_stdev_2 = rootgrp_clim_stdev_2.variables

                #Time variables (and lon/lat)
                #   note mean, min and max of the variables,
                v_name_mat = rootgrp_clim_1.variables.keys()
                for v_name in v_name_mat:
                    if v_name in var_mat: continue
                    #if v_name in ncvar_mat: continue
                    if ei == 0:
                        other_dat[seas][v_name+'_pres_sum'] = var_clim_1[v_name][:].astype('double')
                        other_dat[seas][v_name+'_fut_sum'] = var_clim_2[v_name][:].astype('double')
                        other_dat[seas][v_name+'_pres_min'] = var_clim_1[v_name][:].astype('double')
                        other_dat[seas][v_name+'_fut_min'] = var_clim_2[v_name][:].astype('double')
                        other_dat[seas][v_name+'_pres_max'] = var_clim_1[v_name][:].astype('double')
                        other_dat[seas][v_name+'_fut_max'] = var_clim_2[v_name][:].astype('double')
                    else:
                        other_dat[seas][v_name+'_pres_sum'] += var_clim_1[v_name][:].astype('double')
                        other_dat[seas][v_name+'_fut_sum'] += var_clim_2[v_name][:].astype('double')
                        other_dat[seas][v_name+'_pres_min'] = np.minimum(other_dat[seas][v_name+'_pres_min'],var_clim_1[v_name][:].astype('double'))
                        other_dat[seas][v_name+'_fut_min'] = np.minimum(other_dat[seas][v_name+'_fut_min'],var_clim_2[v_name][:].astype('double'))
                        other_dat[seas][v_name+'_pres_max'] = np.minimum(other_dat[seas][v_name+'_pres_max'],var_clim_1[v_name][:].astype('double'))
                        other_dat[seas][v_name+'_fut_max'] = np.minimum(other_dat[seas][v_name+'_fut_max'],var_clim_2[v_name][:].astype('double'))

                # Variables
                # sum the values and their squares.
                #   First iteration copy them, and following iterations add
                #   Note that we are working in variance, rather than standard deviation, so need to convert the climatology std files into variances... we do this a this stage, but squaring var_clim_stdev_1  etc.
                #for var, ncvar in zip(var_mat,ncvar_mat):
                #for var, ncvar in zip(var_mat,ncvar_mat):
                for var in var_mat:
                    if ei == 0:
                        proc_dat[seas][var]['pres_mean_sum'] = var_clim_1[var][0,:,:].astype('double')
                        proc_dat[seas][var]['fut_mean_sum'] = var_clim_2[var][0,:,:].astype('double')
                        proc_dat[seas][var]['diff_mean_sum'] = var_clim_2[var][0,:,:].astype('double') - var_clim_1[var][0,:,:].astype('double')
                        proc_dat[seas][var]['pres_mean_ss'] = (var_clim_1[var][0,:,:].astype('double'))**2
                        proc_dat[seas][var]['fut_mean_ss'] = (var_clim_2[var][0,:,:].astype('double'))**2
                        proc_dat[seas][var]['diff_mean_ss'] = (var_clim_2[var][0,:,:].astype('double') - var_clim_1[var][0,:,:].astype('double'))**2
                        # we need the interannual variance of each ensmeble member, so at this step need to convert from Std Dev by squaring
                        proc_dat[seas][var]['pres_var_sum'] = (var_clim_stdev_1[var][0,:,:].astype('double')**2)
                        proc_dat[seas][var]['fut_var_sum'] = (var_clim_stdev_2[var][0,:,:].astype('double')**2)
                        proc_dat[seas][var]['diff_var_sum'] = (var_clim_stdev_2[var][0,:,:].astype('double')**2) - (var_clim_stdev_1[var][0,:,:].astype('double')**2)
                    else:
                        #pdb.set_trace()
                        proc_dat[seas][var]['pres_mean_sum'] += var_clim_1[var][0,:,:].astype('double')
                        proc_dat[seas][var][ 'fut_mean_sum'] += var_clim_2[var][0,:,:].astype('double')
                        proc_dat[seas][var]['diff_mean_sum'] += var_clim_2[var][0,:,:].astype('double') - var_clim_1[var][0,:,:].astype('double')
                        proc_dat[seas][var]['pres_mean_ss'] += (var_clim_1[var][0,:,:].astype('double'))**2
                        proc_dat[seas][var][ 'fut_mean_ss'] += (var_clim_2[var][0,:,:].astype('double'))**2
                        proc_dat[seas][var]['diff_mean_ss'] += (var_clim_2[var][0,:,:].astype('double') - var_clim_1[var][0,:,:].astype('double'))**2
                        # we need the interannual variance of each ensmeble member, so at this step need to convert from Std Dev by squaring
                        proc_dat[seas][var]['pres_var_sum'] += (var_clim_stdev_1[var][0,:,:].astype('double')**2)
                        proc_dat[seas][var]['fut_var_sum'] += (var_clim_stdev_2[var][0,:,:].astype('double')**2)
                        proc_dat[seas][var]['diff_var_sum'] += (var_clim_stdev_2[var][0,:,:].astype('double')**2) - (var_clim_stdev_1[var][0,:,:].astype('double')**2)
                # increment ensemble counter
                proc_cnt += 1

                if Test:
                    Test_mean[seas][ei] = var_clim_1[Testvar][0,120,120]
                    Test_std[seas][ei] = var_clim_stdev_1[Testvar][0,120,120]

                #close input files
                rootgrp_clim_1.close()
                rootgrp_clim_stdev_1.close()
                rootgrp_clim_2.close()
                rootgrp_clim_stdev_2.close()


        # process ensemble statistics, from sum, and sum of square
        for seas in date_name_mat:
            # calc the mean for other variables (time, lon, lat (even time_bounds))
            for v_name in v_name_mat:
                #if v_name in ncvar_mat: continue
                if v_name in var_mat: continue
                other_dat[seas][v_name + '_pres_mean'] = other_dat[seas][v_name + '_pres_sum']/proc_cnt
                other_dat[seas][v_name +  '_fut_mean'] = other_dat[seas][v_name +  '_fut_sum']/proc_cnt
            #Variables
            #   Calc the ens_mean, int_var, ens_var (ens_std inferred later)
            #for var, ncvar in zip(var_mat,ncvar_mat):
            for var in var_mat:
                for pername, perlab in zip(['pres', 'fut', 'diff'],[yrmat_1_str,yrmat_2_str,diff_yrmat_str]):

                    #Ens_mean -
                    proj_dat[seas][var]['%s_ensmean'%pername] = proc_dat[seas][var]['%s_mean_sum'%pername]/proc_cnt

                    # int_var
                    #   note we have already converted the clim std into variance when we read the climatology files in.
                    proj_dat[seas][var]['%s_intvar'%pername] = proc_dat[seas][var]['%s_var_sum'%pername]/proc_cnt

                    #ens_var
                    proj_dat[seas][var]['%s_ensvar'%pername] = proc_dat[seas][var]['%s_mean_ss'%pername]/proc_cnt - proj_dat[seas][var]['%s_ensmean'%pername]**2

                    #ens_stddev
                    proj_dat[seas][var]['%s_ensstd'%pername] = np.sqrt(proj_dat[seas][var]['%s_ensvar'%pername])

        if Test:
            pdb.set_trace()
            Test_mean['dec'].mean(),proj_dat['dec'][Testvar]['pres_ensmean'][120,120]

            Test_mean['dec'].var(),proj_dat['dec'][Testvar]['pres_ensvar'][120,120]
            ((Test_std['dec'])**2).mean(),proj_dat['dec'][Testvar]['pres_intvar'][120,120]

        # Write output files.
        #loop through seasons, and then pres, fut and diff.

        for seas in date_name_mat:
            for pername, perlab in zip(['pres', 'fut', 'diff'],[yrmat_1_str,yrmat_2_str,diff_yrmat_str]):
                print(datetime.now(),' Writing files: ',seas,perlab,grid_val)

                # output file name.
                tmpfname_out =  '%sNWSClim_NWSPPE_EnsStats_clim_%s_%s_grid%s_stats.nc' %(path_out,seas,perlab,grid_val)

                # open file, and add dimensions, and variables.
                rootgrp_out = Dataset(tmpfname_out, 'w', format='NETCDF4')


                #Add dimensions
                lat_dim_nc = rootgrp_out.createDimension(lat_dim, 375)
                lon_dim_nc = rootgrp_out.createDimension(lon_dim, 297)
                time_counter_dim_nc = rootgrp_out.createDimension(time_dim, None)
                axis_nbounds_dim_nc = rootgrp_out.createDimension(bnd_dim, 2)


                #Add time variables and attributes
                time_var = rootgrp_out.createVariable('time', 'f8',(time_dim))
                timebnds_var = rootgrp_out.createVariable('time_bounds', 'f8',(time_dim,bnd_dim))

                time_var.setncattr('long_name',"Time axis")
                time_var.setncattr('calendar',"360_day")
                time_var.setncattr('units',"seconds since 1950-01-01 00:00:00")
                time_var.setncattr('time_origin',"1950-01-01 00:00:00")
                time_var.setncattr('bounds',"time_bounds")

                # add time, and the time bounds.
                if pername == 'diff':
                    time_var[:] = (other_dat[seas]['time_fut_mean'][:] + other_dat[seas]['time_pres_mean'][:])/2
                else:
                    time_var[:] = other_dat[seas]['time_%s_mean'%pername][:]

                if pername == 'diff':
                    timebnds_var[0,0] = other_dat[seas]['time_bounds_pres_min'].min()
                    timebnds_var[0,1] = other_dat[seas]['time_bounds_fut_max'].max()
                else:
                    timebnds_var[0,0] = other_dat[seas]['time_bounds_%s_min'%pername][:].min()
                    timebnds_var[0,1] = other_dat[seas]['time_bounds_%s_max'%pername][:].max()


                #Add statistic variables.
                nc_2d_var_dict = {}
                for var in var_mat:
                    for ens_stat in ens_stat_lst:
                        nc_2d_var_dict[var+ens_stat ] = rootgrp_out.createVariable( var + ens_stat ,'f4',(time_dim,lat_dim,lon_dim),fill_value = 1.e+20, zlib=compress_variables)

                for var in var_dict[grid_val]:
                    for ens_stat in ens_stat_lst:

                        tmplongname = ens_stat_long_name_format_dict[ens_stat]% (var.upper(),perlab,long_name_dict[var])
                        tmpunit = ens_stat_units_format_dict[ens_stat]%unit_dict[var]

                        nc_2d_var_dict[var + ens_stat].setncattr('long_name',tmplongname)
                        # No Standard_names for Ens Stats #if var in standard_name_dict.keys():                           tmpvardict[var].setncattr('standard_name',standard_name_dict[var])
                        nc_2d_var_dict[var + ens_stat].setncattr('units',tmpunit)
                        nc_2d_var_dict[var + ens_stat].setncattr('online_operation',"average")
                        nc_2d_var_dict[var + ens_stat].setncattr('interval_write',"1 month")
                        nc_2d_var_dict[var + ens_stat].setncattr('interval_operation',"1 time-step")
                        nc_2d_var_dict[var + ens_stat].setncattr('cell_methods',ens_stat_cell_methods_dict[ens_stat])




                # Add data to variables.
                for var in var_dict[grid_val]:
                    for ens_stat in ens_stat_lst:
                        nc_2d_var_dict[var + ens_stat][0,:,:] = proj_dat[seas][var][pername + ens_stat][:,:]


                #Add lat and lon.
                lon_var = rootgrp_out.createVariable('lon', 'f4',(lon_dim))
                lat_var = rootgrp_out.createVariable('lat', 'f4',(lat_dim))

                lon_var.setncattr('long_name',"Longitude")
                lon_var.setncattr('standard_name',"longitude")
                lon_var.setncattr('units',"degrees_east")
                lon_var.setncattr('nav_model',"grid_%s"%grid_val)
                lat_var.setncattr('long_name',"Latitude")
                lat_var.setncattr('standard_name',"latitude")
                lat_var.setncattr('units',"degrees_north")
                lat_var.setncattr('nav_model',"grid_%s"%grid_val)
                lon_var[:] = lonlat_out_nc_dict['lon_' + grid_val]
                lat_var[:] = lonlat_out_nc_dict['lat_' + grid_val]

                #Add global attributes
                out_source = out_source_NWSPPE

                rootgrp_out.setncattr('Conventions',out_Conventions)
                rootgrp_out.setncattr('institution',out_institution)
                rootgrp_out.setncattr('title',out_title)
                rootgrp_out.setncattr('source',out_source)

                rootgrp_out.setncattr('history',out_history)
                rootgrp_out.setncattr('references',out_references)
                rootgrp_out.setncattr('variance_separation_references',out_variance_separation_references)

                # close files.
                rootgrp_out.setncattr('comments',out_comment_EnsStat)
                rootgrp_out.setncattr('comments_grid',out_comment_grid%(grid_val))
                rootgrp_out.close()

                #pdb.set_trace()


def main():

    run_CEDA_regional_means() #20min ##16:11 - 1634; 1734-1754
    run_CEDA_monthly(file_output_freq_ann = file_output_freq_ann)# 2:15 ##17min PD, 9min/NWSPPE ens mem
    run_CEDA_ens_climatologies(file_output_freq_ann = file_output_freq_ann) #52 mins
    run_CEDA_ens_stats() # 1 mins

    pdb.set_trace()

if __name__ == "__main__":
    main()
