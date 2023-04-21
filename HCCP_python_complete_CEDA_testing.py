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

#sys.path.append('/net/home/h01/hadjt/workspace/python3')
#from climatological_tools import create_mon_lst, create_ann_seas_lst

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
    HCCP_output_dir='/project/shelf_sea_projection/HCCP_UKCP_PPE/proc_files/CEDA_python_ann_Ag/'
    HCCP_output_dir='/project/shelf_sea_projection/HCCP_UKCP_PPE/proc_files/CEDA/'
else:
    HCCP_output_dir='/project/shelf_sea_projection/HCCP_UKCP_PPE/proc_files/CEDA_python_mon/'

if not os.path.exists(HCCP_output_dir): os.makedirs(HCCP_output_dir)



#eg_file_lst = ['/scratch/hadjt/HCCP_UKCP_PPE/Results/amm7_ensemble/ap977_ar095_au084/HCCP_CO9_ap977_ar095_au084_r001i1p00000_02/19900101_Monthly2D_grid_%s.nc'%ss for ss in grid_list]


##################################################################################
##################################################################################



HCCP_ens_mat_12 = ['r001i1p00000', 'r001i1p00605', 'r001i1p00834', 'r001i1p01113', 'r001i1p01554','r001i1p01649', 'r001i1p01843', 'r001i1p01935', 'r001i1p02123', 'r001i1p02242','r001i1p02491', 'r001i1p02868']  #


monstr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
monmxx = ['m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12']



HCCP_subdirs = HCCP_ens_mat_12 + ['PDCtrl','EnsStats']

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
long_name_dict['MLD'] = 'Mixed Layer Depth using Kara approach'
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
#long_name_dict['RegAveMLD'] = 'Regional Mean Mixed Layer Depth using Kara approach'
long_name_dict['RegAvePEA'] = 'Regional Mean Potential Energy Anomaly'


long_name_dict['reg_id'] = 'Regional region-id'
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
# % (var.upper(),perlab.replace('to',' to '),attrib_dict[ncvar]['long_name'])
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
#%unit_dict[var]
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


out_Conventions = "CF-1.6"
out_Conventions = "CF-1.8"
out_institution = "Met Office Hadley Centre, Exeter, UK."
out_title = "Marine climate projections for the North West European Shelf Seas.\nA GCM Perturbed Parameter Ensemble was downscaled with a shelf seas model.\n"

#out_production = "A NEMO consortia model"

out_history = 'Model output post processed by NWS_simulations_postproc.py'
out_references = 'Tinker et al. 2023, A set of climate projections for the NW European Shelf Seas, in prep.'
out_comment_monthly_mean = 'These are monthly mean values calculated from every timestep, as output by NEMO shelf version %s.'

out_comment_regmean = "Regional mean time series are output by NEMO shelf version %s, on the Wakelin et al. (2012) region mask (see region_refererence). NEMO reads in the region mask, and every time step it averages the variable within each mask region. The regional mean is the averaged over the month, and output by NEMO. This file includes the regional mean time-series for a selection of variables (as RegAveSST, RegAveSSS etc. for the regional mean SST and SSS respectively), for each month, and each of the 14 regions. There are variables for the region id (reg_id) and the number of grid boxes within each region (cnt). The region mask included (mask, with the associated longitude and latitude variables). This methodology is described by Tinker et al. (2019; see region_methology_reference)."

out_comment_grid = 'Note that NEMO uses the Arakawa "C" grid, where the T, U and V grids are offset. Most variables are on the T grid, apart from the Eastward and Northward components of the ocean barotropic current, which are on the U and V grids respectively. We also provide the barotropic current speed on T grid, where we transform the U and V velocity components onto the T grid, before calculating their magnitude. We therefore separate the variables on the T, U and V grids into separate files. This file contains variables on the %s grid.'


out_comment_EnsStat = 'These are ensemble statistics are calculated from the climatological mean and standard deviation for each NWSPPE ensemble member, which in turn are based on the monthly mean values calculated from every timestep, as output by NEMO shelf version 4.0.4. For a given climatological period (e.g. 2000-2019, 2079-2098), and month, season or year, the climatological mean and standard deviation of all 12 NWSPPE ensmeble members are loaded. The mean of these climatological means gives the Ensemble Mean (ensmean). The variance of the climatological means gives the Ensemble Variance (ensvar), and its square-root gives the Ensemble Standard Deviation (ensstd). The mean of climatological standard deviation squared (i.e. the climatological variance) gives the Interannual Variance (intvar). Each variable (SST, SSS, NBT, etc.) has 4 statistic associated with it, and the NetCDF variable names join the variable name (e.g. SST) with the statistic name (e.g. ensmean), separated by a hypen. For example the ensmeble mean SST is SST_ensmean and the NBT interannual variability is NBT_intvar. This is all captured in the "cell_methods", where "realization" denotes the ensemble members. The Ensemble Variance and Interannual Variance can be combined into the Total Variance following Tinker et al. 2016; see variance_separation_references).'





# compress variables? smaller, but file size varies, so spot incomplete files
compress_variables = False


def test_file_cf_compliance():
    #pdb.set_trace()
    '''
    rootgrp_out = Dataset('/home/h01/hadjt/workspace/python3/HCCP_UKCP_PPE/tmp.nc','w', format='NETCDF4')
    ens_var = rootgrp_out.createVariable('ensemble', 'f4',fill_value = 1.e+20)
    ens_var.setncattr('long_name',"PPE ensemble member number")
    ens_var.setncattr('units',"1")
    ens_var[0] = -1
    rootgrp_out.region_refererence = 'Adapted from Wakelin, S. L., Holt, J., Blackford, J., Allen, I., Butenschön, M., and Artioli, Y.: Modeling the carbon fluxes of the northwest European continental shelf: Validation and budgets, 117, C05020, https://doi.org/10.1029/2011JC007402, 2012.'.encode('utf8')  #
    rootgrp_out.close()
    pdb.set_trace()






    '''







    # open file, and add dimensions, and variables.
    rootgrp_out = Dataset('/home/h01/hadjt/workspace/python3/HCCP_UKCP_PPE/tmp.nc','w', format='NETCDF4')


    #Add dimensions
    lat_dim_nc = rootgrp_out.createDimension(lat_dim, 375)
    lon_dim_nc = rootgrp_out.createDimension(lon_dim, 297)
    time_counter_dim_nc = rootgrp_out.createDimension(time_dim, None)
    axis_nbounds_dim_nc = rootgrp_out.createDimension(bnd_dim, 2)

    #ens_nc = rootgrp_out.createDimension('ensemble', 1)



    time_var = rootgrp_out.createVariable('time', 'f8',(time_dim))
    timebnds_var = rootgrp_out.createVariable('time_bounds', 'f8',(time_dim,bnd_dim))

    time_var.setncattr('long_name',"Time axis")
    time_var.setncattr('calendar',"360_day")
    time_var.setncattr('units',"seconds since 1950-01-01 00:00:00")
    time_var.setncattr('time_origin',"1950-01-01 00:00:00")
    time_var.setncattr('bounds',"time_bounds")

    time_var[0:1] = 60*86400*360.
    timebnds_var[0:1,:] = time_var[0]  + np.array([-86400,86400])

    #ens_var = rootgrp_out.createVariable('ensemble', 'f8',fill_value = 1.e+20
    #ens_var = rootgrp_out.createVariable('ensemble', 'f8')
    #ens_var.setncattr('long_name',"PPE ensemble member number")
    #ens_var.setncattr('units',"1")
    #ens_var[0] = 1.
    perlab='DJF'
    nc_2d_var_dict = {}
    for var in ['SST']:
        #nc_2d_var_dict[var+'_ensmean' ] = rootgrp_out.createVariable( var + '_ensmean' ,'f4',('ensemble',time_dim,lat_dim,lon_dim),fill_value = 1.e+20, zlib=compress_variables)

        nc_2d_var_dict[var+'_ensmean' ] = rootgrp_out.createVariable( var + '_ensmean' ,'f4',(time_dim,lat_dim,lon_dim),fill_value = 1.e+20, zlib=compress_variables)

        nc_2d_var_dict[var + '_ensmean'].long_name = 'Ensemble mean %s for %s: %s'%(var.upper(),perlab.replace('to',' to '),'SST')

        nc_2d_var_dict[var + '_ensmean'].cell_methods = 'time: mean within years time: mean over years realization: mean'
        #nc_2d_var_dict[var + '_ensmean'].cell_methods = 'ensemble: point'
        nc_2d_var_dict[var + '_ensmean'].units = unit_dict[var]
        #nc_2d_var_dict[var + '_ensmean'].setncattr('coordinates',"ensemble")
        nc_2d_var_dict[var + '_ensmean'][:,:,:] = 0.




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



    rootgrp_out.setncattr('Conventions',out_Conventions)
    rootgrp_out.setncattr('institution',out_institution)
    rootgrp_out.setncattr('title',out_title)
    rootgrp_out.setncattr('source','Dummy File')

    rootgrp_out.close()



    rootgrp_out = Dataset('/home/h01/hadjt/workspace/python3/HCCP_UKCP_PPE/tmp2.nc','w', format='NETCDF4')


    #Add dimensions as per example
    lat_dim_nc = rootgrp_out.createDimension('lat', 180)
    lon_dim_nc = rootgrp_out.createDimension('lon', 360)
    time_counter_dim_nc = rootgrp_out.createDimension('time', None)

    #Add variables as per example
    atime_var = rootgrp_out.createVariable('atime', 'f8')
    time_var = rootgrp_out.createVariable('time', 'f8',('time'))
    lon_var = rootgrp_out.createVariable('lon', 'f8',('lon'))
    lat_var = rootgrp_out.createVariable('lat', 'f8',('lat'))
    p500_var = rootgrp_out.createVariable('p500', 'f8')
    height_var = rootgrp_out.createVariable('height', 'f8',('time','lat','lon'))


    #Add attributes as per example
    atime_var.standard_name = "forecast_reference_time" ;
    atime_var.units = "hours since 1999-01-01 00:00" ;
    time_var.standard_name = "time" ;
    time_var.units = "hours since 1999-01-01 00:00" ;
    lon_var.long_name = "station longitude";
    lon_var.units = "degrees_east";
    lat_var.long_name = "station latitude" ;
    lat_var.units = "degrees_north" ;
    p500_var.long_name = "pressure" ;
    p500_var.units = "hPa" ;
    p500_var.positive = "down" ;
    height_var.long_name = "geopotential height" ;
    height_var.standard_name = "geopotential_height" ;
    height_var.units = "m" ;
    height_var.coordinates = "atime p500" ;

    #add data as per example
    time_var[:4] = [6., 12., 18., 24.]
    atime_var[:1] = 0.
    p500_var[:1] = 500.



    # Add additional variables/attributes etc to pass cfchecker

    # WARN: (4.4.1): Use of the calendar and/or month_lengths attributes is recommended for time coordinate variables
    time_var.setncattr('calendar',"360_day")

    #WARN: (2.6.1): No 'Conventions' attribute present
    rootgrp_out.setncattr('Conventions',out_Conventions)

    # add lon/lat data to pass "ERROR: (5): co-ordinate variable not monotonic"
    lon_var[:] = np.linspace(-20,20,360)
    lat_var[:] = np.linspace(40,650,180)

    #Test use of cell methods from an scalar coordinate variables
    height_var.cell_methods = 'p500: mean'

    rootgrp_out.close()

    """

    dimensions:
      lat = 180 ;
      lon = 360 ;
      time = UNLIMITED ;
    variables:
      double atime
        atime:standard_name = "forecast_reference_time" ;
        atime:units = "hours since 1999-01-01 00:00" ;
      double time(time);
        time:standard_name = "time" ;
        time:units = "hours since 1999-01-01 00:00" ;
      double lon(lon) ;
        lon:long_name = "station longitude";
        lon:units = "degrees_east";
      double lat(lat) ;
        lat:long_name = "station latitude" ;
        lat:units = "degrees_north" ;
      double p500
        p500:long_name = "pressure" ;
        p500:units = "hPa" ;
        p500:positive = "down" ;
      float height(time,lat,lon);
        height:long_name = "geopotential height" ;
        height:standard_name = "geopotential_height" ;
        height:units = "m" ;
        height:coordinates = "atime p500" ;
    data:
      time = 6., 12., 18., 24. ;
      atime = 0. ;
      p500 = 500. ;

    """

    exit()

def test_ens_stats(grid_val, var, iind, jind, yrmat,mn):

    from math import isclose
    nyr,nens = yrmat.size ,12
    tmp_var_mat = np.zeros((nyr,nens))

    for ei,ens in enumerate(HCCP_ens_mat_12):
        print(ens, datetime.now())

        fname_ens = ens.capitalize()

        path_out = HCCP_output_dir + '%s/'%ens
        for yi,yr in enumerate(yrmat):


            file_out = path_out +'%04i%02i%sMonthly2DGrid%s.nc'%(yr,mn,fname_ens,grid_val)


            rootgrp_in = Dataset(file_out, 'r', format='NETCDF4')

            tmp_var_mat[yi,ei] = rootgrp_in.variables[var][0,jind,iind].astype('double')
            rootgrp_in.close()



    ensstatfname = '%s/EnsStats/clim%s%04i-%04iGrid%sstats.nc'%(HCCP_output_dir,monstr[mn-1].lower(),yrmat[0],yrmat[-1],grid_val)
    rootgrp_in = Dataset(ensstatfname, 'r', format='NETCDF4')
    tmp_intvar = rootgrp_in.variables[var + '_intvar'][0,jind,iind].astype('double')
    tmp_ensvar = rootgrp_in.variables[var + '_ensvar'][0,jind,iind].astype('double')
    tmp_ensmean = rootgrp_in.variables[var + '_ensmean'][0,jind,iind].astype('double')
    tmp_ensstd = rootgrp_in.variables[var + '_ensstd'][0,jind,iind].astype('double')
    rootgrp_in.close()

    tmp_totvar = tmp_intvar + tmp_ensvar
    #tmp_ensstd = tmp_ensvar**2

    '''
    print ('ensmean:',tmp_var_mat.mean() ,tmp_ensmean ,tmp_var_mat.mean()/tmp_ensmean ,tmp_var_mat.mean()-tmp_ensmean )

    print ('ensvar:',tmp_var_mat.mean(axis = 0).var() ,tmp_ensvar ,tmp_var_mat.mean(axis = 0).var()/tmp_ensvar ,tmp_var_mat.mean(axis = 0).var()-tmp_ensvar )

    print ('intvar:',tmp_var_mat.var(axis = 0).mean() ,tmp_intvar ,tmp_var_mat.var(axis = 0).mean()/tmp_intvar ,tmp_var_mat.var(axis = 0).mean()-tmp_intvar )

    print ('totvar:',tmp_var_mat.ravel().var() ,tmp_totvar ,tmp_var_mat.ravel().var()/tmp_totvar ,tmp_var_mat.ravel().var()-tmp_totvar )

    print ('ensstd:',tmp_var_mat.mean(axis = 0).var()**2 ,tmp_ensstd ,tmp_var_mat.mean(axis = 0).var()**2/tmp_ensstd ,tmp_var_mat.mean(axis = 0).var()**2-tmp_ensstd )
    '''

    print ('ensmean:',isclose(tmp_var_mat.mean() ,tmp_ensmean,abs_tol=1e-6 ))
    print ('ensvar:',isclose(tmp_var_mat.mean(axis = 0).var() ,tmp_ensvar,abs_tol=1e-6 ))
    print ('intvar:',isclose(tmp_var_mat.var(axis = 0).mean() ,tmp_intvar,abs_tol=1e-6  ))
    print ('totvar:',isclose(tmp_var_mat.ravel().var() ,tmp_totvar,abs_tol=1e-6 ))
    print ('ensstd:',isclose(np.sqrt(tmp_var_mat.mean(axis = 0).var() ),tmp_ensstd,abs_tol=1e-6  ))




    pdb.set_trace()


def main():
    pdb.set_trace()
    test_ens_stats()
    test_file_cf_compliance()

main()
