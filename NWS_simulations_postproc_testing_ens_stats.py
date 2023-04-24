# -*- coding: utf-8 -*-
import os
import glob
import pdb
from netCDF4 import Dataset
import numpy as np
from datetime import datetime

from math import isclose

# Include file with NWS configuration info
import NWS_simulations_postproc_config as NWS_config

NWSPPE_output_dir = NWS_config.NWSPPE_output_dir


NWSPPE_ens_mat_12 = [
    "r001i1p00000", "r001i1p00605", "r001i1p00834",
    "r001i1p01113", "r001i1p01554", "r001i1p01649",
    "r001i1p01843", "r001i1p01935", "r001i1p02123",
    "r001i1p02242", "r001i1p02491", "r001i1p02868",
]

monmxx = [
    "m01", "m02", "m03",
    "m04", "m05", "m06",
    "m07", "m08", "m09",
    "m10", "m11", "m12",
]

monstr = [
    "Jan", "Feb", "Mar",
    "Apr", "May", "Jun",
    "Jul", "Aug", "Sep",
    "Oct", "Nov", "Dec",
]


def test_ens_stats_diff(var, grid_val,  iind, jind, seas):
    """
    Test the ensemble statistic calculations

    For a given variable (var) on a given grid (grid_val), at a given (i,j),
    location, a given month (seas, in the form of "m01" - "m12", not yet set
    up for seasons or annual means), and a given climatological period (per
    "2000-2019", or "2079-2098", not yet set up for "2079-2098minus2000-2019"),
    it loads all the data from the annual files and the data from the
    climatologies, and calculates the ensemble statistic, and compares them to
    those in the ensemble statistic files.
    """



    yrmat_1=np.arange(2000, 2019 + 1)
    yrmat_2=np.arange(2079, 2098 + 1)



    # Labeling output date strings for filesnames
    yrmat_1_str = "%04i-%04i" % (yrmat_1[0], yrmat_1[-1])
    yrmat_2_str = "%04i-%04i" % (yrmat_2[0], yrmat_2[-1])

    diff_yrmat_str = "%sminus%s" % (yrmat_2_str, yrmat_1_str)


    if seas in monmxx:
        mn = int(seas[1:])
    else:
        print('not set up for seasonal or annual means')

    nyr,nens = yrmat_1.size ,12
    ann_var_mat = np.zeros((2,nyr,nens))

    for pi,(pername, perlab,yrmat) in enumerate(zip(["pres", "fut"], [yrmat_1_str, yrmat_2_str],[yrmat_1,yrmat_2])):

        for ei,ens in enumerate(NWSPPE_ens_mat_12):
            print(ens, datetime.now())

            fname_ens = ens


            path = NWSPPE_output_dir + "NWSClim/NWSPPE/%s/annual/" % ens

            for yi,yr in enumerate(yrmat):

                fname_date = "%04i" % (yr)

                file_out = (
                    "NWSClim_NWSPPE_%s_%s_grid%s.nc"
                    % (fname_ens,fname_date,grid_val,)
                )


                rootgrp_in = Dataset(path + file_out, 'r', format='NETCDF4')

                ann_var_mat[pi,yi,ei] = rootgrp_in.variables[var][mn-1,jind,iind].astype('double')
                rootgrp_in.close()


    calc_ensmean = ann_var_mat.mean(axis = 1).mean(axis = 1)
    calc_ensvar = ann_var_mat.mean(axis = 1).var(axis = 1)
    calc_ensstd = np.sqrt(ann_var_mat.mean(axis = 1).var(axis = 1) )
    calc_intvar = ann_var_mat.var(axis = 1).mean(axis = 1)
    calc_totvar = ann_var_mat.reshape(2,-1).var(axis = 1)




    clim_mean_mat = np.zeros((2,nens))
    clim_std_mat = np.zeros((2,nens))



    for pi,(pername, perlab,yrmat) in enumerate(zip(["pres", "fut"], [yrmat_1_str, yrmat_2_str],[yrmat_1,yrmat_2])):
        for ei,ens in enumerate(NWSPPE_ens_mat_12):
            print(ens, datetime.now())

            fname_ens = ens


            clim_path = NWSPPE_output_dir + "NWSClim/NWSPPE/%s/clim/" % ens


            tmpfname_clim = (
                "%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_mean.nc"
                % (clim_path,ens,seas,perlab,grid_val)
            )
            tmpfname_clim_stdev = (
                "%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_stddev.nc"
                % (clim_path, ens, seas, perlab, grid_val)
            )
            rootgrp_in_mean = Dataset(tmpfname_clim, 'r', format='NETCDF4')
            rootgrp_in_std = Dataset(tmpfname_clim_stdev, 'r', format='NETCDF4')

            clim_mean_mat[pi,ei] = rootgrp_in_mean.variables[var][0,jind,iind].astype('double')
            clim_std_mat[pi,ei] = rootgrp_in_std.variables[var][0,jind,iind].astype('double')
            rootgrp_in_mean.close()
            rootgrp_in_std.close()




    clim_ensmean = clim_mean_mat.mean(axis = 1)
    clim_ensvar = clim_mean_mat.var(axis = 1)
    clim_ensstd = np.sqrt(clim_mean_mat.var(axis = 1) )
    clim_intvar = ((clim_std_mat)**2 ).mean(axis = 1)
    clim_totvar = clim_ensvar + clim_intvar



    ensstat_path = "%s/NWSClim/EnsStats/" % NWSPPE_output_dir


    ensstat_intvar = np.ma.zeros((3))
    ensstat_ensvar = np.ma.zeros((3))
    ensstat_ensmean = np.ma.zeros((3))
    ensstat_ensstd = np.ma.zeros((3))
    ensstat_projstd = np.ma.zeros((1))

    for pi,(pername, perlab) in enumerate(zip(["pres", "fut", "diff"], [yrmat_1_str, yrmat_2_str, diff_yrmat_str])):

        ensstat_file_out = (
            "%sNWSClim_NWSPPE_EnsStats_clim_%s_%s_grid%s_stats.nc"
            % (ensstat_path, seas, perlab, grid_val)
        )


        #ensstatfname = tmpfname_out
        rootgrp_in = Dataset(ensstat_file_out, 'r', format='NETCDF4')
        ensstat_intvar[pi] = rootgrp_in.variables[var + '_intvar'][0,jind,iind].astype('double')
        ensstat_ensvar[pi] = rootgrp_in.variables[var + '_ensvar'][0,jind,iind].astype('double')
        ensstat_ensmean[pi] = rootgrp_in.variables[var + '_ensmean'][0,jind,iind].astype('double')
        ensstat_ensstd[pi] = rootgrp_in.variables[var + '_ensstd'][0,jind,iind].astype('double')
        if pername == "diff":
            ensstat_projstd = rootgrp_in.variables[var + '_projstd'][0,jind,iind].astype('double')
        rootgrp_in.close()

    ensstat_totvar = ensstat_intvar + ensstat_ensvar

    #pdb.set_trace()
    print ('Tests of climatologies and ensemble statistics for %s (%s %s)'%(var, monstr[mn-1],perlab))
    print ('====================================')
    print ('do the climatologies calculated from the annual files agree with the climatology files?')

    for pi,(pername, perlab,yrmat) in enumerate(zip(["pres", "fut"], [yrmat_1_str, yrmat_2_str],[yrmat_1,yrmat_2])):
        print()
        print(perlab)
        print()
        for ei,ens in enumerate(NWSPPE_ens_mat_12):
            ann_clim_mean = ann_var_mat[pi,:,ei].mean()
            ann_clim_std = ann_var_mat[pi,:,ei].std()

            print (ens,'climmean:',isclose(ann_clim_mean,clim_mean_mat[pi,ei],abs_tol=1e-6 ),'climstd:',isclose(ann_clim_std,clim_std_mat[pi,ei],abs_tol=1e-6 ))
        print()



    for pi,(pername, perlab,yrmat) in enumerate(zip(["pres", "fut"], [yrmat_1_str, yrmat_2_str],[yrmat_1,yrmat_2])):
        print()
        print(perlab)
        print()
        print('Do the ens stat calculated from the Ens Clim agree with the Ens Stats File?')
        print ('ensmean:',isclose(clim_ensmean[pi],ensstat_ensmean[pi],abs_tol=1e-6 ))
        print ('ensvar:',isclose(clim_ensvar[pi],ensstat_ensvar[pi],abs_tol=1e-6 ))
        print ('ensstd:',isclose(clim_ensstd[pi],ensstat_ensstd[pi],abs_tol=1e-6  ))
        print ('intvar:',isclose(clim_intvar[pi],ensstat_intvar[pi],abs_tol=1e-6  ))
        print ('totvar:',isclose(clim_totvar[pi],ensstat_totvar[pi],abs_tol=1e-6 ))
        print()


        print('Do the ens stat calculated from the annual file agree with the Ens Stats File?')
        print ('ensmean:',isclose(calc_ensmean[pi],ensstat_ensmean[pi],abs_tol=1e-6 ))
        print ('ensvar:',isclose(calc_ensvar[pi],ensstat_ensvar[pi],abs_tol=1e-6 ))
        print ('ensstd:',isclose(calc_ensstd[pi],ensstat_ensstd[pi],abs_tol=1e-6  ))
        print ('intvar:',isclose(calc_intvar[pi],ensstat_intvar[pi],abs_tol=1e-6  ))
        print ('totvar:',isclose(calc_totvar[pi],ensstat_totvar[pi],abs_tol=1e-6 ))
        print()


        print('Do the ens stat calculated from the annual file agree with those calculated from the clims?')
        print ('ensmean:',isclose(calc_ensmean[pi],clim_ensmean[pi],abs_tol=1e-6 ))
        print ('ensvar:',isclose(calc_ensvar[pi],clim_ensvar[pi],abs_tol=1e-6 ))
        print ('ensstd:',isclose(calc_ensstd[pi],clim_ensstd[pi],abs_tol=1e-6  ))
        print ('intvar:',isclose(calc_intvar[pi],clim_intvar[pi],abs_tol=1e-6  ))
        print ('totvar:',isclose(calc_totvar[pi],clim_totvar[pi],abs_tol=1e-6 ))
        print()

    print()
    print('Difference between periods')
    print()



    # calc diff ens stats from annual files
    ############################################
    '''
    calc_ensmean_rembl = np.diff(ann_var_mat.mean(axis = 1),axis = 0).mean() # 20 yr mean, diff between periods, mean over ens

    calc_ensvar_rembl = np.diff(ann_var_mat.mean(axis = 1),axis = 0).var() # 20 yr mean, diff between periods, var over ens
    calc_ensstd_rembl = np.diff(ann_var_mat.mean(axis = 1),axis = 0).std() # 20 yr mean, diff between periods, std over ens

    calc_intvar_rembl = np.diff(ann_var_mat.var(axis = 1),axis = 0).mean() # 20 yr var, diff between periods, mean over ens
    calc_totvar_rembl =  calc_ensvar_rembl + calc_intvar_rembl
    '''
    calc_projstd = np.diff(ann_var_mat.mean(axis = 1),axis = 0).std() # 20 yr mean, diff between periods, std over ens

    # calc diff ens stats from climatologies
    clim_mean_mat_rembl = np.diff(clim_mean_mat, axis= 0)
    clim_var_mat_rembl = np.diff(clim_std_mat**2, axis= 0)
    clim_projstd = np.sqrt(clim_mean_mat_rembl.var() )



    clim_ensmean_rembl = clim_mean_mat_rembl.mean()
    clim_ensvar_rembl = clim_mean_mat_rembl.var()
    clim_ensstd_rembl = np.sqrt(clim_mean_mat_rembl.var() )
    clim_intvar_rembl = ((clim_var_mat_rembl) ).mean()
    clim_totvar_rembl = clim_ensvar_rembl + clim_intvar_rembl







    ensstat_ensmean_diff = ensstat_ensmean[1]-ensstat_ensmean[0]
    ensstat_ensvar_diff = ensstat_ensvar[1]-ensstat_ensvar[0]
    ensstat_ensstd_diff = ensstat_ensstd[1]-ensstat_ensstd[0]
    ensstat_intvar_diff = ensstat_intvar[1]-ensstat_intvar[0]
    ensstat_totvar_diff = ensstat_totvar[1]-ensstat_totvar[0]



    print('Which ens stat differ when removing the baseline before calculating the stats vs the differnce of the stats between the period?')
    print ('ensmean:',isclose(ensstat_ensmean_diff, clim_ensmean_rembl,abs_tol=1e-6 ),ensstat_ensmean_diff, clim_ensmean_rembl )
    print ('ensmean at 1e-5:',isclose(ensstat_ensmean_diff, clim_ensmean_rembl,abs_tol=1e-5 ),ensstat_ensmean_diff, clim_ensmean_rembl )
    print ('ensvar:',isclose(ensstat_ensvar_diff, clim_ensvar_rembl,abs_tol=1e-6 ),ensstat_ensvar_diff, clim_ensvar_rembl )
    print ('ensstd:',isclose(ensstat_ensstd_diff, clim_ensstd_rembl,abs_tol=1e-6 ),ensstat_ensstd_diff, clim_ensstd_rembl )
    print ('intvar:',isclose(ensstat_intvar_diff, clim_intvar_rembl,abs_tol=1e-6 ),ensstat_intvar_diff, clim_intvar_rembl )
    print ('totvar:',isclose(ensstat_totvar_diff, clim_totvar_rembl,abs_tol=1e-6 ), ensstat_totvar_diff, clim_totvar_rembl)


    print()
    print()

    print('Do the ens stat diff equal the difference between the future and present ens stats?')

    print ('ensmean:',isclose(ensstat_ensmean[1]-ensstat_ensmean[0],ensstat_ensmean[2],abs_tol=1e-6 ))
    print ('ensmean at 1e-5:',isclose(ensstat_ensmean[1]-ensstat_ensmean[0],ensstat_ensmean[2],abs_tol=1e-5 ))
    print ('ensvar:',isclose(ensstat_ensvar[1]-ensstat_ensvar[0],ensstat_ensvar[2],abs_tol=1e-6 ))
    print ('ensstd:',isclose(ensstat_ensstd[1]-ensstat_ensstd[0],ensstat_ensstd[2],abs_tol=1e-6  ))
    print ('intvar:',isclose(ensstat_intvar[1]-ensstat_intvar[0],ensstat_intvar[2],abs_tol=1e-6  ))
    print ('totvar:',isclose(ensstat_totvar[1]-ensstat_totvar[0],ensstat_totvar[2],abs_tol=1e-6 ))


    print()
    print()

    print('Do the ens difference stat calculated from the Ens Clim agree with the Ens Stats File?')
    print ('ensmean:',isclose(np.diff(clim_ensmean),ensstat_ensmean[2],abs_tol=1e-6 ))
    print ('ensvar:',isclose(np.diff(clim_ensvar),ensstat_ensvar[2],abs_tol=1e-6 ))
    print ('ensstd:',isclose(np.diff(clim_ensstd),ensstat_ensstd[2],abs_tol=1e-6  ))
    print ('intvar:',isclose(np.diff(clim_intvar),ensstat_intvar[2],abs_tol=1e-6  ))
    print ('totvar:',isclose(np.diff(clim_totvar),ensstat_totvar[2],abs_tol=1e-6 ))
    print ('projstd:',isclose(clim_projstd,ensstat_projstd,abs_tol=1e-6 ))


    print()
    print()
    print('Do the ens difference stat calculated from the annual files agree with the Ens Stats File?')
    print ('ensmean:',isclose(np.diff(calc_ensmean),ensstat_ensmean[2],abs_tol=1e-6 ))
    print ('ensvar:',isclose(np.diff(calc_ensvar),ensstat_ensvar[2],abs_tol=1e-6 ))
    print ('ensstd:',isclose(np.diff(calc_ensstd),ensstat_ensstd[2],abs_tol=1e-6  ))
    print ('intvar:',isclose(np.diff(calc_intvar),ensstat_intvar[2],abs_tol=1e-6  ))
    print ('totvar:',isclose(np.diff(calc_totvar),ensstat_totvar[2],abs_tol=1e-6 ))
    print ('projstd:',isclose(calc_projstd,ensstat_projstd,abs_tol=1e-6 ))


    print()
    print()
    print('Do the ens difference stat calculated from the annual files agree with those calculated from the Ens Clim File?')
    print ('ensmean:',isclose(np.diff(calc_ensmean),np.diff(clim_ensmean),abs_tol=1e-6 ))
    print ('ensvar:',isclose(np.diff(calc_ensvar),np.diff(clim_ensvar),abs_tol=1e-6 ))
    print ('ensstd:',isclose(np.diff(calc_ensstd),np.diff(clim_ensstd),abs_tol=1e-6  ))
    print ('intvar:',isclose(np.diff(calc_intvar),np.diff(clim_intvar),abs_tol=1e-6  ))
    print ('totvar:',isclose(np.diff(calc_totvar),np.diff(clim_totvar),abs_tol=1e-6 ))
    print ('projstd:',isclose(calc_projstd,clim_projstd,abs_tol=1e-6 ))


    print()
    print()
    pdb.set_trace()


def main():
    var='SST'
    grid_val='T'
    iind = 120
    jind= 120
    seas = 'm07'# 'm01' - 'm12', not set up for djf, mam, jja, son or ann

    test_ens_stats_diff(var, grid_val,  iind, jind, seas)
    pdb.set_trace()



if __name__ == "__main__":
    main()
