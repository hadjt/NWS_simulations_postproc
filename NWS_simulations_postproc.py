# -*- coding: utf-8 -*-
import os
import glob
import pdb
from netCDF4 import Dataset
import numpy as np
from datetime import datetime

# Include file with NWS configuration info
import NWS_simulations_postproc_config as NWS_config

# Include file with NWS dictionaries
import NWS_simulations_postproc_dict as NWS_dict

# Include file with NWS global attributes
import NWS_simulations_postproc_globattrib as NWS_globattrib

# Include file with NWS standardised lon lat
import NWS_simulations_postproc_lon_lat as NWS_lon_lat

import warnings

warnings.filterwarnings("ignore")

# Use info from NWS configuration info

PDCtrl_path_in = NWS_config.PDCtrl_path_in
NWSPPE_path_in_pattern = NWS_config.NWSPPE_path_in_pattern
reg_mask_mat = NWS_config.reg_mask_mat
NWSPPE_output_dir = NWS_config.NWSPPE_output_dir

# Use info from NWS dictionaries

long_name_dict = NWS_dict.long_name_dict
standard_name_dict = NWS_dict.standard_name_dict
unit_dict = NWS_dict.unit_dict
orig_var_name_dict = NWS_dict.orig_var_name_dict
ens_stat_long_name_format_dict = NWS_dict.ens_stat_long_name_format_dict
ens_stat_cell_methods_dict = NWS_dict.ens_stat_cell_methods_dict
ens_stat_units_format_dict = NWS_dict.ens_stat_units_format_dict
var_dict = NWS_dict.var_dict

# Use info from NWS global attributes

out_Conventions = NWS_globattrib.out_Conventions
out_institution = NWS_globattrib.out_institution
out_title = NWS_globattrib.out_title
out_history = NWS_globattrib.out_history
out_references = NWS_globattrib.out_references
out_comment_monthly_mean = NWS_globattrib.out_comment_monthly_mean
out_comment_regmean = NWS_globattrib.out_comment_regmean
out_region_names = NWS_globattrib.out_region_names
out_region_notes = NWS_globattrib.out_region_notes
out_region_refererence = NWS_globattrib.out_region_refererence
out_region_methology_reference = NWS_globattrib.out_region_methology_reference
out_source_PDCtrl = NWS_globattrib.out_source_PDCtrl
out_source_NWSPPE = NWS_globattrib.out_source_NWSPPE
out_comment_grid = NWS_globattrib.out_comment_grid
out_comment_EnsStat = NWS_globattrib.out_comment_EnsStat
out_variance_separation_references = NWS_globattrib.out_variance_separation_references

# Use info from NWS standardised lon lat

lonlat_out_nc_dict = NWS_lon_lat.lonlat_out_nc_dict

# Create root dir if doesn't exist

if not os.path.exists(NWSPPE_output_dir):
    os.makedirs(NWSPPE_output_dir)

# List of grids, ensemble members, months

grid_list = ["T", "U", "V"]

NWSPPE_ens_mat_12 = [
    "r001i1p00000", "r001i1p00605", "r001i1p00834",
    "r001i1p01113", "r001i1p01554", "r001i1p01649",
    "r001i1p01843", "r001i1p01935", "r001i1p02123",
    "r001i1p02242", "r001i1p02491", "r001i1p02868",
]

monstr = [
    "Jan", "Feb", "Mar",
    "Apr", "May", "Jun",
    "Jul", "Aug", "Sep",
    "Oct", "Nov", "Dec",
]

monmxx = [
    "m01", "m02", "m03",
    "m04", "m05", "m06",
    "m07", "m08", "m09",
    "m10", "m11", "m12",
]

mnmat = np.arange(1, 12 + 1)

# Dimension names
reg_dim = "region"
lat_dim = "lat"
lon_dim = "lon"
time_dim = "time"
bnd_dim = "bnds"

# compress variables? smaller, but file size varies, so harder to
# spot incomplete files
compress_variables = False


def create_ann_seas_lst(yrmat):
    """
    When we are creating seasonal or annual mean climatologies for a given
    period (i.e. 2000-2019, if yrmat = np.arange(2000,2019+1)), we need to
    average a number of months into a seasonal or annual mean first. For
    example, if we're making a spring climatology for 2000-2019, we first
    need to average March, April and May for 2000 into spring 2000, and then
    March, April and May for 2001 into spring 2001.

    This function give a list of years and months that need to be averaged
    within the climatological period (yrmat).

    create_ann_seas_lst returns a list with 5 elements (for annual, winter,
    spring, summer and autumn means), with each element being another list of
    the year and months (in the string format YYYYMM) needed for the averaging.

    For Example:

    ann_lst, djf_lst, mam_lst, jja_lst, son_lst =
        create_ann_seas_lst(np.arange(2000,2019+1))

    # the months needed to average into springs, between 2000 and 2019
    mam_lst
    [['200003', '200004', '200005'], ['200103', '200104', '200105'], ['200203', '200204', '200205'], ['200303', '200304', '200305'], ['200403', '200404', '200405'], ['200503', '200504', '200505'], ['200603', '200604', '200605'], ['200703', '200704', '200705'], ['200803', '200804', '200805'], ['200903', '200904', '200905'], ['201003', '201004', '201005'], ['201103', '201104', '201105'], ['201203', '201204', '201205'], ['201303', '201304', '201305'], ['201403', '201404', '201405'], ['201503', '201504', '201505'], ['201603', '201604', '201605'], ['201703', '201704', '201705'], ['201803', '201804', '201805'], ['201903', '201904', '201905']]

    Jonathan Tinker 29/03/2023
    """

    # Predefine lists for annual and the 4 seasonal means
    ann_lst = []
    djf_lst = []
    jja_lst = []
    mam_lst = []
    son_lst = []

    # cycle through the years of the climatological period
    for yr in yrmat:

        # Predefine temporary lists
        tmp_ann_lst = []
        tmp_djf_lst = []
        tmp_jja_lst = []
        tmp_mam_lst = []
        tmp_son_lst = []

        # cycle through and append the year and month in YYYYMM string format

        # cycle through 12 months for annual means
        for mi in range(12):
            tmp_ann_lst.append("%04i%02i" % (yr, mi + 1))

        # cycle through Mar, Apr and May for Spring (MAM)
        for mi in [2, 3, 4]:
            tmp_mam_lst.append("%04i%02i" % (yr, mi + 1))

        # cycle through Jun, Jul and Aug for Summer (JJA)
        for mi in [5, 6, 7]:
            tmp_jja_lst.append("%04i%02i" % (yr, mi + 1))

        # cycle through Sep, Oct and Nov for Autumn (SON)
        for mi in [8, 9, 10]:
            tmp_son_lst.append("%04i%02i" % (yr, mi + 1))

        # For Winter (DJF) append Dec from previous year
        tmp_djf_lst.append("%04i12" % (yr - 1))

        # and then cycle through Jan and Feb
        for mi in [0, 1]:
            tmp_djf_lst.append("%04i%02i" % (yr, mi + 1))

        # Append current year list to main lists
        ann_lst.append(tmp_ann_lst)
        djf_lst.append(tmp_djf_lst)
        mam_lst.append(tmp_mam_lst)
        jja_lst.append(tmp_jja_lst)
        son_lst.append(tmp_son_lst)

    # Return as
    return ann_lst, djf_lst, mam_lst, jja_lst, son_lst


def create_mon_lst(yrmat):
    """
    The provides a list with 12 elements (one for each month of year). Each
    element has a list of years within the climatological period (i.e.
    2000-2019, if yrmat = np.arange(2000,2019+1)). This is a little redundant,
    but is consistent with create_ann_seas_lst (see details above), and allows
    the same code to produce annual, seasonal and monthly climatologies


    For Example:

    tmp_mon = create_mon_lst(np.arange(2000,2019+1))

    # the years needed for the February climatology for between 2000 and 2019
    tmp_mon[1]
    [['200002'], ['200102'], ['200202'], ['200302'], ['200402'], ['200502'], ['200602'], ['200702'], ['200802'], ['200902'], ['201002'], ['201102'], ['201202'], ['201302'], ['201402'], ['201502'], ['201602'], ['201702'], ['201802'], ['201902']]

    Jonathan Tinker 29/03/2023
    """

    # Predefine list
    mon_lst = []

    # cycle through months
    for mi in range(12):
        # Predefine temporary lists
        tmp_mon_lst = []
        # cycle through the years of the climatological period
        # append the year and month in YYYYMM string format
        for yr in yrmat:
            tmp_mon_lst.append(["%04i%02i" % (yr, mi + 1)])
        # Append current year list to main list
        mon_lst.append(tmp_mon_lst)

    return mon_lst


def run_CEDA_regional_means():
    """
    Read in the regional mean files output from NEMO COx and create a single regional mean file for CEDA, with standardised and corrected meta data.

    This is unlikely to be run outside the MO, as the raw NEMO model output
    will not be made available.

    This code cycles through the PDCtrl and the NWSPPE ensemble members.
    For each, it
        opens an output file,
        adds the dimensions, variables and attributes.
        It then cycles through the years and months, and takes the data from
        the shelf region, and the Wakelin et al. region mask, and adds them
        to the output file.

    NEMO 3.6 and NEMO 4.0.4 (PDCtrl and NWSPPE respectively) both have a
    different time origin (1950 and 1900 respectively) so this is corrected
    for.
    Furthermore, NEMO 3.6 (PDCtrl) had an XIOS bug which meant the time
    (in seconds since 1950) would overflow with values >2**32. We correct
    for this.

    There are occasional files that were not archived, and so have been
    recreated (mostly by averaging daily mean files, see paper for details).
    This is done externally to this program, using netcdf operators.
    Files that have been recreated have been given an extension to the file
    name of the raw NEMO output. To account for this, we use glob.glob with
    the input file names.

    Jonathan Tinker 29/03/2023
    """

    grid_val = "R"

    # Create the directory if missing.
    if not os.path.exists(NWSPPE_output_dir):
        os.makedirs(NWSPPE_output_dir)

    # cycle through the ensemble members
    for ens in ["PDCtrl"] + NWSPPE_ens_mat_12:
        print(ens, datetime.now())

        # input path and years for NWSPPE and PDCtrl differ.
        if ens == "PDCtrl":
            yrmat = np.arange(2050, 2250 + 1)

            path_in = PDCtrl_path_in

            fname_ens = ens
            file_out = "NWSClim_PDCtrl_2050-2250_regmean.nc"
            path_out = NWSPPE_output_dir + "NWSClim/%s/regmean/" % ens

        else:
            yrmat = np.arange(1990, 2098 + 1)

            path_in = NWSPPE_path_in_pattern % (ens)

            fname_ens = "NWSPPE_%s" % ens
            file_out = "NWSClim_%s_1990-2098_regmean.nc" % fname_ens
            path_out = NWSPPE_output_dir + "NWSClim/NWSPPE/%s/regmean/" % ens

        if not os.path.exists(path_out):
            os.makedirs(path_out)

        # Open output file, create variables and attributes.
        rootgrp_out = Dataset(path_out + file_out, "w", format="NETCDF4")

        rootgrp_out.createDimension(reg_dim, 14)
        rootgrp_out.createDimension(lat_dim, 375)
        rootgrp_out.createDimension(lon_dim, 297)
        rootgrp_out.createDimension(time_dim, None)
        rootgrp_out.createDimension(bnd_dim, 2)

        # create time variables, and add attributes
        time_var = rootgrp_out.createVariable("time", "f8", (time_dim))
        timebnds_var = rootgrp_out.createVariable(
            "time_bounds", "f8", (time_dim, bnd_dim)
        )
        time_var.setncattr("long_name", "Time axis")
        time_var.setncattr("calendar", "360_day")
        time_var.setncattr("units", "seconds since 1950-01-01 00:00:00")
        time_var.setncattr("time_origin", "1950-01-01 00:00:00")
        time_var.setncattr("bounds", "time_bounds")

        # Create data variables
        tmpvardict = {}
        for var in var_dict[grid_val]:
            tmpvardict[var] = rootgrp_out.createVariable(
                var,
                "f4",
                (time_dim, reg_dim),
                fill_value=1.0e20,
                zlib=compress_variables,
            )

        # Create region id and count variables
        for var in ["reg_id", "cnt"]:
            tmpvardict[var] = rootgrp_out.createVariable(
                var, "f4", (reg_dim), fill_value=1.0e20, zlib=compress_variables
            )

        # Create mask variables
        for var in ["mask"]:
            tmpvardict[var] = rootgrp_out.createVariable(
                var,
                "f4",
                (lat_dim, lon_dim),
                fill_value=1.0e20,
                zlib=compress_variables,
            )

        # Add variable attributes
        for var in var_dict[grid_val] + ["reg_id", "cnt"]:
            tmpvardict[var].setncattr("long_name", long_name_dict[var])
            if var in standard_name_dict.keys():
                tmpvardict[var].setncattr("standard_name", standard_name_dict[var])
            tmpvardict[var].setncattr("units", unit_dict[var])

        for var in var_dict[grid_val]:
            tmpvardict[var].setncattr("online_operation", "average")
            tmpvardict[var].setncattr("interval_write", "1 month")
            tmpvardict[var].setncattr("interval_operation", "1 time-step")
            tmpvardict[var].setncattr("cell_methods", "Area: mean time: mean")

        for var in ["mask"]:
            tmpvardict[var].setncattr("long_name", long_name_dict[var])
            tmpvardict[var].setncattr("units", unit_dict[var])

        # Create lon and lats variables
        lon_var = rootgrp_out.createVariable("lon", "f4", (lon_dim))
        lat_var = rootgrp_out.createVariable("lat", "f4", (lat_dim))

        # Add lon and lats attributes
        lon_var.setncattr("long_name", "Longitude")
        lon_var.setncattr("standard_name", "longitude")
        lon_var.setncattr("units", "degrees_east")
        lon_var.setncattr("nav_model", "grid_T")
        lat_var.setncattr("long_name", "Latitude")
        lat_var.setncattr("standard_name", "latitude")
        lat_var.setncattr("units", "degrees_north")
        lat_var.setncattr("nav_model", "grid_T")

        # Add lon and lats data
        lon_var[:] = lonlat_out_nc_dict["lon_T"]
        lat_var[:] = lonlat_out_nc_dict["lat_T"]

        # Add mask data
        tmpvardict["mask"][:, :] = reg_mask_mat

        # Add Global Attributes
        if ens == "PDCtrl":
            out_source = out_source_PDCtrl
            tmp_nemo_version = "3.6"
        else:
            out_source = out_source_NWSPPE
            tmp_nemo_version = "4.0.4"

        rootgrp_out.setncattr("Conventions", out_Conventions)
        rootgrp_out.setncattr("institution", out_institution)
        rootgrp_out.setncattr("title", out_title)
        rootgrp_out.setncattr("source", out_source)

        rootgrp_out.setncattr("region_names", out_region_names)
        rootgrp_out.setncattr("region_notes", out_region_notes)
        rootgrp_out.setncattr("region_refererence", out_region_refererence)
        rootgrp_out.setncattr(
            "region_methology_reference", out_region_methology_reference
        )

        rootgrp_out.setncattr("history", out_history)
        rootgrp_out.setncattr("references", out_references)
        rootgrp_out.setncattr("comments", out_comment_regmean % tmp_nemo_version)

        # loop through year and month with a counter.
        counter = 0
        for yr in yrmat:
            # print out date every decade
            if (yr % 10) == 0:
                print(yr, datetime.now())
            for mn in mnmat:

                # Open files, with glob to catch recreated files
                file_in = glob.glob(
                    "%s%04i%02i01_Region_mean_timeseries_Monthly*.nc"
                    % (path_in, yr, mn)
                )
                # stop if more than one file founds by glob.glob
                if len(file_in) != 0:
                    print("incorrect number of files found, stopping")
                    print(file_in)
                    pdb.set_trace()
                rootgrp_in = Dataset(file_in[0], "r", format="NETCDF4")

                # time: Check origin, 1950 or 1900
                orig_time_counter_origin = rootgrp_in.variables[
                    "time_counter"
                ].time_origin

                # find offset between 1900 and 1950
                if orig_time_counter_origin == "1900-01-01 00:00:00":
                    calendar_time_offset = 50.0 * 360.0 * 86400.0
                elif orig_time_counter_origin == "1950-01-01 00:00:00":
                    calendar_time_offset = 0.0
                else:
                    print("Check time_counter origin")
                    pdb.set_trace()

                # PDCtrl overflows when second counter gets over 2**32, so work around.
                if ens == "PDCtrl":
                    # seconds since 1950, for the middle of the month
                    tmp_calc_time_counter = (
                        (yr - 1950.0) * 360.0 * 86400.0
                        + (mn - 1.0) * 30.0 * 86400
                        + 15.0 * 86400
                    )

                    # seconds since 1950, for the start/end of the month
                    tmp_calc_time_counter_bounds = np.zeros((2))
                    tmp_calc_time_counter_bounds[0] = (
                        yr - 1950.0
                    ) * 360.0 * 86400.0 + (mn - 1.0) * 30.0 * 86400.0
                    tmp_calc_time_counter_bounds[1] = (
                        (yr - 1950.0) * 360.0 * 86400.0
                        + (mn - 1.0) * 30.0 * 86400.0
                        + (30.0 * 86400.0)
                    )

                    # test work around.
                    # confirm that the time calculation matches that from NEMO, when accounting from the 2**32 overflow
                    # only check for time bounds as they are instantaneous
                    if not (
                        (tmp_calc_time_counter_bounds % (2**32))
                        == (
                            (
                                rootgrp_in.variables["time_counter_bounds"][:, :]
                                - calendar_time_offset
                            )
                            % (2**32)
                        )
                    ).all():
                        print(
                            "Calculated time bounds (instantaneous values) not offset by 2**32 - check"
                        )
                        pdb.set_trace()

                    # Append time and time_bounds to output file
                    time_var[counter : counter + 1] = tmp_calc_time_counter
                    timebnds_var[
                        counter : counter + 1, :
                    ] = tmp_calc_time_counter_bounds
                else:

                    # if NEMO4.0.4 (NWSPPE), add model time but correct
                    # for 1900-1950 change in origin
                    time_var[counter : counter + 1] = (
                        rootgrp_in.variables["time_counter"][:] - calendar_time_offset
                    )
                    timebnds_var[counter : counter + 1, :] = (
                        rootgrp_in.variables["time_counter_bounds"][:, :]
                        - calendar_time_offset
                    )

                # check which regions are the Shelf, and Wakelin regions.
                tmp_cnt = rootgrp_in.variables["cnt"][:].ravel()
                tmp_reg_id = rootgrp_in.variables["reg_id"][:].ravel()
                tmp_mask_id = rootgrp_in.variables["mask_id"][:].ravel()

                # Pull out the shelf region (mask_id = 2, reg_id = 1), and the wakelin regions (mask_id = 3, reg_id > 1)
                reg_ind = ((tmp_mask_id == 2) & (tmp_reg_id == 1)) | (
                    (tmp_mask_id == 3) & (tmp_reg_id > 0)
                )

                # Add count, and reg_id data
                tmpvardict["cnt"][:] = tmp_cnt[reg_ind]
                tmpvardict["reg_id"][:] = tmp_reg_id[reg_ind]

                # fill the variables with data.
                for var in var_dict[grid_val]:
                    # find the NEMO output nc variable name
                    tmp_orig_var_name_dict = orig_var_name_dict[var]

                    # read from NEMO output file, and write to output file.
                    tmpvardict[var][counter : counter + 1, :] = rootgrp_in.variables[
                        tmp_orig_var_name_dict
                    ][:].ravel()[reg_ind]

                # increment monthly counter
                counter += 1
                # close the input file for that month
                rootgrp_in.close()
        # After cycling through all years and months, close output file.
        rootgrp_out.close()


def run_CEDA_monthly():
    """
    Read in the monthly mean files output from NEMO COx and recreate them for CEDA, with standardised and corrected meta data.

    This is unlikely to be run outside the MO, as the raw NEMO model output
    will not be made available.

    NEMO runs in monthly cycles, so produces a separate file for each month.
    Here we combine 12 monthly means into a single annual file (of 12 months).
    To do this, we write an output file every January, add the dimensions,
    variables, attributes, and then close the file. We then open this file
    with append, to add the data every month.

    This code cycles through the PDCtrl and the NWSPPE ensemble members, and through the T, U and V grids.

    For each,
        every January it
            opens an output file,
            adds the dimensions, variables and attributes.
            closes the file
        It then cycles through the years and months, and adds the
        data to the output file.

    NEMO 3.6 and NEMO 4.0.4 (PDCtrl and NWSPPE respectively) both have a
    different time origin (1950 and 1900 respectively) so this is corrected
    for.
    Furthermore, NEMO 3.6 (PDCtrl) had an XIOS bug which meant the time
    (in seconds since 1950) would overflow with values >2**32. We correct
    for this.

    There are occasional files that were not archived, and so have been
    recreated (mostly by averaging daily mean files, see paper for details).
    This is done externally to this program, using netcdf operators.
    Files that have been recreated have been given an extension to the file
    name of the raw NEMO output. To account for this, we use glob.glob with
    the input file names.

    We have provided barotropic (depth mean) current speed (DMUV). This is the
    magnitude of the U and V velocity components (DMU and DMV). However, these
    are on the U and V grid respectively. We therefore have to interpolate the
    U and V components from the U and V grid, to the T grid, before calculating
    the magnitude. Therefore, for the T grid files, we still have to load the
    U and V files - we therefore use a temporary file name dictionary for each
    month.

    Jonathan Tinker 29/03/2023
    """

    # Create the directory if missing.
    if not os.path.exists(NWSPPE_output_dir):
        os.makedirs(NWSPPE_output_dir)

    # cycle through the ensemble members
    for ens in ["PDCtrl"] + NWSPPE_ens_mat_12:
        print(ens, datetime.now())

        # input path and years for NWSPPE and PDCtrl differ.
        if ens == "PDCtrl":
            path_in = PDCtrl_path_in

            yrmat = np.arange(2050, 2250 + 1)
            fname_ens = ens
            path_out = NWSPPE_output_dir + "NWSClim/%s/annual/" % ens
        else:
            path_in = NWSPPE_path_in_pattern % (ens)
            fname_ens = "NWSPPE_%s" % ens

            yrmat = np.arange(1990, 2098 + 1)
            path_out = NWSPPE_output_dir + "NWSClim/NWSPPE/%s/annual/" % ens

        # output path, create it missing.
        if not os.path.exists(path_out):
            os.makedirs(path_out)

        # Cycle through years, months and grid (T, U, V)
        for yr in yrmat:
            if (yr % 10) == 0:
                print(yr, datetime.now())
            for mi, mn in enumerate(mnmat):

                # input file names
                tmp_file_in_dict = {}
                for grid_val in grid_list:

                    tmp_file_in = glob.glob(
                        "%s%04i%02i01_Monthly2D_grid_%s*.nc"
                        % (path_in, yr, mn, grid_val)
                    )

                    if len(tmp_file_in) != 0:
                        print("incorrect number of files found, stopping")
                        print(tmp_file_in)
                        pdb.set_trace()

                    tmp_file_in_dict[grid_val] = tmp_file_in[0]

                for grid_val in grid_list:

                    #  output and input file names
                    fname_date = "%04i" % (yr)
                    mi_ind = mi

                    file_out = "NWSClim_%s_%s_grid%s.nc" % (
                        fname_ens,
                        fname_date,
                        grid_val,
                    )
                    file_in = tmp_file_in_dict[grid_val]

                    # Open input and output file names, and then open them
                    rootgrp_in = Dataset(file_in, "r", format="NETCDF4")

                    # Write file with basic dimensions, variables, and attributes, then close.
                    # then reopen with append, and add data - this allows annual files of 12 months to be easily written

                    # if first month of the year, or if monthly output write file.

                    if mi == 0:
                        rootgrp_out = Dataset(
                            path_out + file_out, "w", format="NETCDF4"
                        )

                        # Add dimensions, variables and attributes to output files.
                        rootgrp_out.createDimension(lat_dim, 375)
                        rootgrp_out.createDimension(lon_dim, 297)
                        rootgrp_out.createDimension(time_dim, None)
                        rootgrp_out.createDimension(bnd_dim, 2)

                        time_var = rootgrp_out.createVariable("time", "f8", (time_dim))
                        timebnds_var = rootgrp_out.createVariable(
                            "time_bounds", "f8", (time_dim, bnd_dim)
                        )

                        # add output variable, and use a dictionary as the handle.
                        tmpvardict = {}
                        for var in var_dict[grid_val]:
                            tmpvardict[var] = rootgrp_out.createVariable(
                                var,
                                "f4",
                                (time_dim, lat_dim, lon_dim),
                                fill_value=1.0e20,
                                zlib=compress_variables,
                            )

                        # add time variables
                        time_var.setncattr("long_name", "Time axis")
                        time_var.setncattr("calendar", "360_day")
                        time_var.setncattr("units", "seconds since 1950-01-01 00:00:00")
                        time_var.setncattr("time_origin", "1950-01-01 00:00:00")
                        time_var.setncattr("bounds", "time_bounds")

                        # add lon/lat variables
                        lon_var = rootgrp_out.createVariable("lon", "f4", (lon_dim))
                        lat_var = rootgrp_out.createVariable("lat", "f4", (lat_dim))

                        # add attributes to lon/lat
                        lon_var.setncattr("long_name", "Longitude")
                        lon_var.setncattr("standard_name", "longitude")
                        lon_var.setncattr("units", "degrees_east")
                        lon_var.setncattr("nav_model", "grid_%s" % grid_val)
                        lat_var.setncattr("long_name", "Latitude")
                        lat_var.setncattr("standard_name", "latitude")
                        lat_var.setncattr("units", "degrees_north")
                        lat_var.setncattr("nav_model", "grid_%s" % grid_val)

                        # Add attributes to variables
                        for var in var_dict[grid_val]:
                            tmpvardict[var].setncattr("long_name", long_name_dict[var])
                            if var in standard_name_dict.keys():
                                tmpvardict[var].setncattr(
                                    "standard_name", standard_name_dict[var]
                                )
                            tmpvardict[var].setncattr("units", unit_dict[var])
                            tmpvardict[var].setncattr("online_operation", "average")
                            tmpvardict[var].setncattr("interval_write", "1 month")
                            tmpvardict[var].setncattr(
                                "interval_operation", "1 time-step"
                            )
                            tmpvardict[var].setncattr("cell_methods", "time: mean")

                        if ens == "PDCtrl":
                            out_source = out_source_PDCtrl
                            tmp_nemo_version = "3.6"
                        else:
                            out_source = out_source_NWSPPE
                            tmp_nemo_version = "4.0.4"

                        rootgrp_out.setncattr("Conventions", out_Conventions)
                        rootgrp_out.setncattr("institution", out_institution)
                        rootgrp_out.setncattr("title", out_title)
                        rootgrp_out.setncattr("source", out_source)

                        rootgrp_out.setncattr("history", out_history)
                        rootgrp_out.setncattr("references", out_references)
                        rootgrp_out.setncattr(
                            "comments", out_comment_monthly_mean % (tmp_nemo_version)
                        )
                        rootgrp_out.setncattr(
                            "comments_grid", out_comment_grid % (grid_val)
                        )

                        # Add lon/lat
                        lon_var[:] = lonlat_out_nc_dict["lon_" + grid_val]
                        lat_var[:] = lonlat_out_nc_dict["lat_" + grid_val]

                        rootgrp_out.close()

                    if grid_val == "T":
                        rootgrp_U_in = Dataset(
                            tmp_file_in_dict["U"], "r", format="NETCDF4"
                        )
                        rootgrp_V_in = Dataset(
                            tmp_file_in_dict["V"], "r", format="NETCDF4"
                        )

                        tmp_DMU = rootgrp_U_in.variables[orig_var_name_dict["DMU"]][:]
                        tmp_DMV = rootgrp_V_in.variables[orig_var_name_dict["DMV"]][:]
                        DMU_T = np.ma.zeros((1, 375, 297)) * np.ma.masked
                        DMV_T = np.ma.zeros((1, 375, 297)) * np.ma.masked
                        DMUV_T = np.ma.zeros((1, 375, 297)) * np.ma.masked

                        # Interpolate DMU and DMV from U and V grid, to T grid.
                        # Confirmed with Jeff Polton, National Oceanography Centre
                        DMU_T[0, 1:, 1:] = (
                            tmp_DMU[0, 1:, :-1] + tmp_DMU[0, 1:, 1:]
                        ) / 2.0
                        DMV_T[0, 1:, 1:] = (
                            tmp_DMV[0, :-1, 1:] + tmp_DMV[0, 1:, 1:]
                        ) / 2.0

                        # Calculate the barotropic current speed
                        DMUV_T = np.sqrt(DMU_T**2 + DMV_T**2)

                        rootgrp_U_in.close()
                        rootgrp_V_in.close()

                    rootgrp_out = Dataset(path_out + file_out, "a", format="NETCDF4")

                    # Add dimensions, variables and attributes to output files.
                    time_var = rootgrp_out.variables["time"]
                    timebnds_var = rootgrp_out.variables["time_bounds"]

                    # add output variable, and use a dictionary as the handle.
                    tmpvardict = {}
                    for var in var_dict[grid_val]:
                        tmpvardict[var] = rootgrp_out.variables[var]

                    # add attributes to time.
                    # PDCtrl and NWSPPE have different time origins (1900 and 1950)
                    # Also PDCtrl has a time overflow problem, when seconds are > 2**32.
                    orig_time_counter_origin = rootgrp_in.variables[
                        "time_counter"
                    ].time_origin
                    if orig_time_counter_origin == "1900-01-01 00:00:00":
                        calendar_time_offset = 50.0 * 360.0 * 86400.0
                    elif orig_time_counter_origin == "1950-01-01 00:00:00":
                        calendar_time_offset = 0.0
                    else:
                        print("time_counter origin needs checking")

                    # Implement a time overflow work around for when PDCtrl seconds are > 2**32.
                    if ens == "PDCtrl":
                        tmp_calc_time_counter = (
                            (yr - 1950.0) * 360.0 * 86400.0
                            + (mn - 1.0) * 30.0 * 86400
                            + 15.0 * 86400
                        )
                        tmp_calc_time_counter_bounds = np.zeros((2))
                        tmp_calc_time_counter_bounds[0] = (
                            yr - 1950.0
                        ) * 360.0 * 86400.0 + (mn - 1.0) * 30.0 * 86400.0
                        tmp_calc_time_counter_bounds[1] = (
                            (yr - 1950.0) * 360.0 * 86400.0
                            + (mn - 1.0) * 30.0 * 86400.0
                            + (30.0 * 86400.0)
                        )

                        if not (
                            (tmp_calc_time_counter_bounds % (2**32))
                            == (
                                (
                                    rootgrp_in.variables["time_counter_bounds"][:, :]
                                    - calendar_time_offset
                                )
                                % (2**32)
                            )
                        ).all():
                            print(
                                "Calculated time bounds (instantaneous values) not offset by 2**32 - check"
                            )
                            pdb.set_trace()()

                        time_var[mi_ind] = tmp_calc_time_counter
                        timebnds_var[mi_ind, :] = tmp_calc_time_counter_bounds

                    else:
                        time_var[mi_ind] = (
                            rootgrp_in.variables["time_counter"][:]
                            - calendar_time_offset
                        )
                        timebnds_var[mi_ind, :] = (
                            rootgrp_in.variables["time_counter_bounds"][:, :]
                            - calendar_time_offset
                        )

                    # Loop through the variables, finding the NEMO original nc var name, and outputing with the new standardised names
                    for var in var_dict[grid_val]:

                        if var == "DMUV":
                            continue
                        tmp_orig_var_name_dict = orig_var_name_dict[var]
                        # PDCtrl and NWSPPE have different nc var names for MLD, so work around.
                        if ens == "PDCtrl":
                            if var == "MLD":
                                tmp_orig_var_name_dict = (
                                    "Mixed_Layer_Depth_Kara_et_al_definition"
                                )

                        tmpvardict[var][mi_ind, :, :] = rootgrp_in.variables[
                            tmp_orig_var_name_dict
                        ][:]

                    if grid_val == "T":
                        # add SST mask, so setting 15 additional sea points to land
                        DMUV_T.mask = tmpvardict["SST"][:].mask | DMUV_T.mask
                        tmpvardict["DMUV"][mi_ind, :, :] = DMUV_T

                    rootgrp_in.close()
                    rootgrp_out.close()


def run_CEDA_ens_climatologies(
    skip_existing=False,
    Test=False,
    Testvar=None,
    yrmat_1=np.arange(2000, 2019 + 1),
    yrmat_2=np.arange(2079, 2098 + 1),
):
    """
    Read in the CEDA reprocessed monthly mean files, and create climatological
    means and standard deviations,for a near present day period (2000-2019)
    and a end of century period (2079-2098). This is repeated for all
    ensemble member for the NWSPPE,and is produced for monthly, seasonal
    and annual means.

    This may be run outside the MO, it only relies on the processed annual
    files which are available on CEDA. Re-running it may allow the user to
    customise the climatological periods.

    We do not create climatologies for PDCtrl, only NWSPPE.
    The climatological means and standard deviations are in
    different files, and the variable names are unchanged (although the
    processing and statisitic is recorded in the variable attributes,
    cell_methods). As these files have the same structure, they are created at
    the same time, and closed, and the reopened (with append) to add the data.

    For the seasonal and annual mean climatologies, we cycle through the years,
    and average the approprate months, by summing them up, and dividing by the
    number of months.

    When caluculating the climatological mean and standard deviations, we sum
    up the variables for the months/seasons/years, and their square. After
    cycling through the years of the climatological period, we divide the sum
    by the number of years to give the climatological mean. For the
    climatological standard deviation we take the take the difference between
    the mean of the squares, and the square of the mean.

    This code cycles through the climatological periods, NWSPPE ensemble
    members, the monthly, annual and seaonal files, and finally through
    the T, U and V grids.

    For each,
        If opens the mean and standard deviation file
        adds the dimensions, variables and attributes.
        closes the file

        Average the require months for that season or year (or month)
        cycle over the yeras of the climatology and
            increment the sum of the variables
            increment the sum of the square of the variables

        Then it converts these sums into the climatological mean and standard
        deviations and reopens the output files, adds the variables, and
        updates the cell methods.

    Care must be taken to average the time, but to take the min and max of the
    time bounds.

    Jonathan Tinker 29/03/2023
    """

    for yrmat in [yrmat_1, yrmat_2]:
        yrstr = "%04i-%04i" % (yrmat[0], yrmat[-1])

        ann_lst, djf_lst, mam_lst, jja_lst, son_lst = create_ann_seas_lst(yrmat)
        mon_lst = create_mon_lst(yrmat)

        date_name_a_s_mat = ["ann", "djf", "mam", "jja", "son"]
        date_name_m_mat = [ii.lower() for ii in monmxx]
        date_name_mat = date_name_m_mat + date_name_a_s_mat

        for ens in NWSPPE_ens_mat_12:

            input_dir = NWSPPE_output_dir + "NWSClim/NWSPPE/%s/annual/" % ens
            path_out = NWSPPE_output_dir + "NWSClim/NWSPPE/%s/clim/" % ens
            if not os.path.exists(path_out):
                os.makedirs(path_out)

            # cycle through monthly, annual and seasonal climatologies
            for date_lst, date_name in zip(
                mon_lst + [ann_lst, djf_lst, mam_lst, jja_lst, son_lst], date_name_mat
            ):

                # How many months (files) in each period.
                if date_name.capitalize() in monstr:
                    tmp_nfiles = 1
                    monthnamestr = date_name
                elif date_name in monmxx:
                    tmp_nfiles = 1
                    monthnamestr = np.array(monstr)[date_name == np.array(monmxx)][0]
                elif date_name.upper() in ["DJF", "MAM", "JJA", "SON"]:
                    tmp_nfiles = 3
                elif date_name.upper() in ["ANN"]:
                    tmp_nfiles = 12

                # Cycle through the grids
                for grid_val in grid_list:

                    output_file_pref = "NWSClim_NWSPPE_%s_clim_%s_%04i-%04i_grid%s" % (
                        ens,
                        date_name,
                        yrmat[0],
                        yrmat[-1],
                        grid_val,
                    )

                    tmpfname_clim = path_out + output_file_pref + "_mean.nc"
                    tmpfname_clim_stdev = path_out + output_file_pref + "_stddev.nc"

                    print(datetime.now(), tmpfname_clim)
                    file_pref = "NWSClim_NWSPPE_%s_" % (ens)
                    file_suff = "_grid%s.nc" % (grid_val)

                    if skip_existing:
                        if os.path.isfile(tmpfname_clim) & os.path.isfile(
                            tmpfname_clim_stdev
                        ):
                            return

                    time_bnd = {}
                    if Test:
                        test_lst = []

                    # year within climatological period
                    # e.g. for presnet day annual means,  cycling through years, 2000,2001,2002, 2003...
                    for yi, tmp_date_lst in enumerate(date_lst):

                        # cycle through the months within period for the climatological period
                        # e.g. for presnet day annual means, for the frrst year cycling through 200001,200002,200003,200004..
                        int_cnt = 0
                        for td, tmp_date in enumerate(tmp_date_lst):  #
                            tmp_fname = input_dir + file_pref + tmp_date[:4] + file_suff
                            mi = int(tmp_date[4:]) - 1

                            if not os.path.isfile(tmp_fname):
                                print("missing file")
                                pdb.set_trace()

                            rootgrp = Dataset(tmp_fname, "r", format="NETCDF4")
                            var = rootgrp.variables
                            if td == 0:
                                tmp_key_mat = list(var.keys())
                                tmp_key_mat.remove("lon")
                                tmp_key_mat.remove("lat")

                                int_var_mat = {}
                                for tmp_key in tmp_key_mat:
                                    int_var_mat[tmp_key] = var[tmp_key][
                                        mi : mi + 1
                                    ].astype("double")
                                int_cnt = +1
                            else:
                                for tmp_key in tmp_key_mat:
                                    int_var_mat[tmp_key][:] += var[tmp_key][
                                        mi : mi + 1
                                    ].astype("double")
                                int_cnt += 1

                            # identify time and time bound var names, and take a copy of the time bounds.
                            # later iterations, take the min and max of the time.
                            if (td == 0) & (yi == 0):
                                tbnd_key_lst = [
                                    ss
                                    for ss in tmp_key_mat
                                    if (ss[-7:] == "_bounds") & (ss[:4] == "time")
                                ]
                                time_key_lst = [
                                    ss
                                    for ss in tmp_key_mat
                                    if (ss[-7:] != "_bounds") & (ss[:4] == "time")
                                ]
                                for tmp_tbnd_key in tbnd_key_lst:
                                    time_bnd[tmp_tbnd_key] = var[tmp_tbnd_key][
                                        mi : mi + 1, :
                                    ]
                            else:

                                for tmp_tbnd_key in tbnd_key_lst:
                                    time_bnd[tmp_tbnd_key][0, 0] = np.minimum(
                                        time_bnd[tmp_tbnd_key][0, 0],
                                        var[tmp_tbnd_key][mi, 0],
                                    )
                                    time_bnd[tmp_tbnd_key][0, 1] = np.maximum(
                                        time_bnd[tmp_tbnd_key][0, 1],
                                        var[tmp_tbnd_key][mi, 1],
                                    )

                            rootgrp.close()

                        # if didn't count any files, stop
                        if int_cnt == 0:
                            print("no files this season")
                            pdb.set_trace()

                        # if didn't count any files, stop
                        if int_cnt != tmp_nfiles:
                            print(
                                "not enough files for ",
                                date_name,
                                "Expected ",
                                tmp_nfiles,
                                "got",
                                int_cnt,
                            )
                            pdb.set_trace()

                        # internal variable average matrix
                        # monthly, seasonal or annual mean.

                        int_var_ave_mat = {}
                        for tmp_key in tmp_key_mat:
                            int_var_ave_mat[tmp_key] = int_var_mat[tmp_key][:] / int_cnt
                        if Test:
                            test_lst.append(int_var_ave_mat[Testvar][0, 120, 120])

                        # if the first year, take a copy, and a copy of the square
                        if yi == 0:
                            var_mat = {}
                            var_mat_stdev = {}
                            for tmp_key in tmp_key_mat:
                                var_mat[tmp_key] = int_var_ave_mat[tmp_key][:]
                                var_mat_stdev[tmp_key] = (
                                    int_var_ave_mat[tmp_key][:] ** 2
                                )
                            tot_cnt = 1

                        # otherwise increment the sum, and the sum of the square, and the counter
                        else:
                            for tmp_key in tmp_key_mat:
                                var_mat[tmp_key][:] += int_var_ave_mat[tmp_key][:]
                                var_mat_stdev[tmp_key][:] += (
                                    int_var_ave_mat[tmp_key][:] ** 2
                                )
                            tot_cnt += 1

                    if Test:
                        for tmp_key in tmp_key_mat:
                            tmp_mean = var_mat[Testvar][:] / tot_cnt
                            tmp_std = np.sqrt(
                                var_mat_stdev[Testvar][:] / tot_cnt - tmp_mean**2
                            )
                            print(tmp_mean[0, 120, 120], np.array(test_lst).mean())
                            print(tmp_std[0, 120, 120], np.array(test_lst).std())

                            pdb.set_trace()

                    # create mean and standard deviation files at the same time
                    for tmp_fill_fname in [tmpfname_clim, tmpfname_clim_stdev]:

                        rootgrp_out = Dataset(
                            tmp_fill_fname, "w", format="NETCDF4_CLASSIC"
                        )
                        rootgrp_out.createDimension(lat_dim, 375)
                        rootgrp_out.createDimension(lon_dim, 297)
                        rootgrp_out.createDimension(time_dim, None)
                        rootgrp_out.createDimension(bnd_dim, 2)

                        time_var = rootgrp_out.createVariable("time", "f8", (time_dim))
                        timebnds_var = rootgrp_out.createVariable(
                            "time_bounds", "f8", (time_dim, bnd_dim)
                        )

                        # add output variable, and use a dictionary as the handle.
                        tmpvardict = {}
                        for var in var_dict[grid_val]:
                            tmpvardict[var] = rootgrp_out.createVariable(
                                var,
                                "f4",
                                (time_dim, lat_dim, lon_dim),
                                fill_value=1.0e20,
                                zlib=compress_variables,
                            )

                        # add time variables
                        time_var.setncattr("long_name", "Time axis")
                        time_var.setncattr("calendar", "360_day")
                        time_var.setncattr("units", "seconds since 1950-01-01 00:00:00")
                        time_var.setncattr("time_origin", "1950-01-01 00:00:00")
                        time_var.setncattr("bounds", "time_bounds")

                        # add lon/lat variables
                        lon_var = rootgrp_out.createVariable("lon", "f4", (lon_dim))
                        lat_var = rootgrp_out.createVariable("lat", "f4", (lat_dim))

                        # add attributes to lon/lat
                        lon_var.setncattr("long_name", "Longitude")
                        lon_var.setncattr("standard_name", "longitude")
                        lon_var.setncattr("units", "degrees_east")
                        lon_var.setncattr("nav_model", "grid_%s" % grid_val)
                        lat_var.setncattr("long_name", "Latitude")
                        lat_var.setncattr("standard_name", "latitude")
                        lat_var.setncattr("units", "degrees_north")
                        lat_var.setncattr("nav_model", "grid_%s" % grid_val)

                        # Add attributes to variables
                        for var in var_dict[grid_val]:
                            tmpvardict[var].setncattr("long_name", long_name_dict[var])
                            if var in standard_name_dict.keys():
                                tmpvardict[var].setncattr(
                                    "standard_name", standard_name_dict[var]
                                )
                            tmpvardict[var].setncattr("units", unit_dict[var])
                            tmpvardict[var].setncattr("online_operation", "average")
                            tmpvardict[var].setncattr("interval_write", "1 month")
                            tmpvardict[var].setncattr(
                                "interval_operation", "1 time-step"
                            )
                            tmpvardict[var].setncattr("cell_methods", "time: mean")

                        out_source = out_source_NWSPPE

                        rootgrp_out.setncattr("Conventions", out_Conventions)
                        rootgrp_out.setncattr("institution", out_institution)
                        rootgrp_out.setncattr("title", out_title)
                        rootgrp_out.setncattr("source", out_source)

                        rootgrp_out.setncattr("history", out_history)
                        rootgrp_out.setncattr("references", out_references)
                        # Add lon/lat
                        lon_var[:] = lonlat_out_nc_dict["lon_" + grid_val]
                        lat_var[:] = lonlat_out_nc_dict["lat_" + grid_val]

                        rootgrp_out.close()

                    rootgrp_out = Dataset(tmpfname_clim, "a", format="NETCDF4_CLASSIC")
                    rootgrp_stdev_out = Dataset(
                        tmpfname_clim_stdev, "a", format="NETCDF4_CLASSIC"
                    )

                    # cycle through the variables, calcuate the mean and standard devation, and write them out.
                    for tmp_key in tmp_key_mat:
                        tmp_mean = var_mat[tmp_key][:] / tot_cnt
                        tmp_std = np.sqrt(
                            var_mat_stdev[tmp_key][:] / tot_cnt - tmp_mean**2
                        )
                        # tmp_variance=(var_mat_stdev[tmp_key][:]/tot_cnt - tmp_mean**2)
                        out_var = rootgrp_out.variables[tmp_key]
                        out_var[:] = tmp_mean
                        out_stdev_clim = rootgrp_stdev_out.variables[tmp_key]
                        out_stdev_clim[:] = tmp_std

                        # update the cell_methods
                        if tmp_key not in ["time", "time_bounds", "lon", "lat"]:
                            out_var.setncattr(
                                "cell_methods",
                                "time: mean within years time: mean over years",
                            )
                            out_stdev_clim.setncattr(
                                "cell_methods",
                                "time: mean within years time: standard_deviation over years",
                            )

                    # update the time bounds, and the time to reflect the climatologies
                    for tmp_tbnd_key in tbnd_key_lst:
                        rootgrp_out.variables[tmp_tbnd_key][:] = time_bnd[tmp_tbnd_key][
                            :
                        ]
                        rootgrp_stdev_out.variables[tmp_tbnd_key][:] = time_bnd[
                            tmp_tbnd_key
                        ][:]

                    for tmp_time_key in time_key_lst:
                        rootgrp_stdev_out.variables[tmp_time_key][
                            :
                        ] = rootgrp_out.variables[tmp_time_key][:]

                    # update the comment global attribute
                    if (date_name.capitalize() in monstr) | (date_name in monmxx):
                        out_comment_mean = (
                            "These are climatological means for each ensemble member. The %s monthly means (as output by NEMO shelf version 4.0.4.) within the current climatological period (%s) are selected, and their mean is taken to give these climatological means."
                            % (monthnamestr, yrstr)
                        )

                        out_comment_stdev = (
                            "These are climatological standard deviations for each ensemble member. The %s monthly means (as output by NEMO shelf version 4.0.4.) within the current climatological period (%s) are selected, and their standard deviation is taken to give these climatological standard deviation."
                            % (monthnamestr, yrstr)
                        )

                    elif date_name.upper() in ["DJF", "MAM", "JJA", "SON"]:

                        out_comment_mean = (
                            "These are climatological means for each ensemble member. The monthly means (as output by NEMO shelf version 4.0.4.) are averaged into seasonal means for %s. The seasonal means within the current climatological period (%s) are selected, and their mean is taken to give these climatological means."
                            % (date_name.upper(), yrstr)
                        )

                        out_comment_stdev = (
                            "These are climatological standard deviations for each ensemble member. The monthly means (as output by NEMO shelf version 4.0.4.) are averaged into seasonal means for %s. The seasonal means within the current climatological period (%s) are selected, and their standard deviation is taken to give these climatological standard deviation."
                            % (date_name.upper(), yrstr)
                        )

                    elif date_name.upper() in ["ANN"]:

                        out_comment_mean = (
                            "These are climatological means for each ensemble member. The monthly means (as output by NEMO shelf version 4.0.4.) are averaged into annual means. The annual means within the current climatological period (%s) are selected, and their mean is taken to give these climatological means."
                            % (yrstr)
                        )

                        out_comment_stdev = (
                            "These are climatological standard deviations for each ensemble member. The monthly means (as output by NEMO shelf version 4.0.4.) are averaged into annual means. The annual means within the current climatological period (%s) are selected, and their standard deviation is taken to give these climatological standard deviation."
                            % (yrstr)
                        )
                    else:
                        print("check date_name", date_name)

                    rootgrp_out.setncattr("comments", out_comment_mean)
                    rootgrp_stdev_out.setncattr("comments", out_comment_stdev)

                    rootgrp_out.setncattr(
                        "comments_grid", out_comment_grid % (grid_val)
                    )
                    rootgrp_stdev_out.setncattr(
                        "comments_grid", out_comment_grid % (grid_val)
                    )

                    rootgrp_out.close()
                    rootgrp_stdev_out.close()


def run_CEDA_ens_stats(
    Test=False,
    Testvar=None,
    yrmat_1=np.arange(2000, 2019 + 1),
    yrmat_2=np.arange(2079, 2098 + 1),
):
    """
    Read in the CEDA climatological means and standard deviations, for ensemble member for the NWSPPE,
    for the near present day period (2000-2019) and end of century period (2079-2098), for the monthly, seasonal and annual means,
    and calculate ensemble statistics. Ensemble statistics include ensmeble mean,  ensmeble standard deviations,  ensmeble variance and interannual variance,
    for the present day, end of century, and the difference between them. All statistics are in the same file, with the variable name modified with a suffix to indentify the statistic (and the processing and statistic is recorded in the variable attributes, cell_methods)

    Jonathan Tinker 29/03/2023
    """

    print("start: ", datetime.now())

    path_out = "%s/NWSClim/EnsStats/" % NWSPPE_output_dir

    # if the output dir is missing, create it
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    date_name_mat = ["ann", "djf", "mam", "jja", "son"] + [ii.lower() for ii in monmxx]

    ens_stat_lst = ["_ensmean", "_ensvar", "_intvar", "_ensstd"]

    # Labeling output date strings for filesnames
    yrmat_1_str = "%04i-%04i" % (yrmat_1[0], yrmat_1[-1])
    yrmat_2_str = "%04i-%04i" % (yrmat_2[0], yrmat_2[-1])
    diff_yrmat_str = "%sminus%s" % (yrmat_2_str, yrmat_1_str)

    for grid_val in grid_list:

        var_mat = var_dict[grid_val]

        # create processing and projections dictionaries
        # proc dict to increment the sum and the sum of the squares
        proc_dat = {}

        proj_dat = {}
        proj_dat["yrmat_1"] = yrmat_1
        proj_dat["yrmat_2"] = yrmat_2

        # other dictionary for lon, lat, time etc.
        other_dat = {}

        for seas in date_name_mat:
            proc_dat[seas] = {}
            proj_dat[seas] = {}
            other_dat[seas] = {}

            for var in var_mat:
                proc_dat[seas][var] = {}
                proj_dat[seas][var] = {}

            if Test:
                Test_mean = {}
                Test_std = {}

        # loop through date types (months, seasons, annuals)
        #   open files for present day, and future climatologies (means and std devs)
        #   sum and sum squares of the variables.
        for di, seas in enumerate(date_name_mat):

            # counter for looping over ensmeble members
            proc_cnt = 0

            if Test:
                Test_mean[seas] = np.ma.zeros(12)
                Test_std[seas] = np.ma.zeros(12)

            print(
                datetime.now(),
                " Loading data clim_%s_%s_grid%s" % (seas, yrmat_1_str, grid_val),
            )
            for ei, ens in enumerate(NWSPPE_ens_mat_12):
                input_dir = NWSPPE_output_dir + "NWSClim/NWSPPE/%s/clim/" % ens

                # file names
                tmpfname_clim_1 = "%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_mean.nc" % (
                    input_dir,
                    ens,
                    seas,
                    yrmat_1_str,
                    grid_val,
                )
                tmpfname_clim_2 = "%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_mean.nc" % (
                    input_dir,
                    ens,
                    seas,
                    yrmat_2_str,
                    grid_val,
                )
                tmpfname_clim_stdev_1 = (
                    "%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_stddev.nc"
                    % (input_dir, ens, seas, yrmat_1_str, grid_val)
                )
                tmpfname_clim_stdev_2 = (
                    "%sNWSClim_NWSPPE_%s_clim_%s_%s_grid%s_stddev.nc"
                    % (input_dir, ens, seas, yrmat_2_str, grid_val)
                )

                # open files
                rootgrp_clim_1 = Dataset(tmpfname_clim_1, "r", format="NETCDF4")
                rootgrp_clim_2 = Dataset(tmpfname_clim_2, "r", format="NETCDF4")
                rootgrp_clim_stdev_1 = Dataset(
                    tmpfname_clim_stdev_1, "r", format="NETCDF4"
                )
                rootgrp_clim_stdev_2 = Dataset(
                    tmpfname_clim_stdev_2, "r", format="NETCDF4"
                )

                # find variables
                var_clim_1 = rootgrp_clim_1.variables
                var_clim_stdev_1 = rootgrp_clim_stdev_1.variables
                var_clim_2 = rootgrp_clim_2.variables
                var_clim_stdev_2 = rootgrp_clim_stdev_2.variables

                # Time variables (and lon/lat)
                #   note mean, min and max of the variables,
                v_name_mat = rootgrp_clim_1.variables.keys()
                for v_name in v_name_mat:
                    if v_name in var_mat:
                        continue

                    # For the first ensemble member, take the value, and subsequently, increment the sum, and take the new min/max
                    if ei == 0:
                        other_dat[seas][v_name + "_pres_sum"] = var_clim_1[v_name][
                            :
                        ].astype("double")
                        other_dat[seas][v_name + "_fut_sum"] = var_clim_2[v_name][
                            :
                        ].astype("double")
                        other_dat[seas][v_name + "_pres_min"] = var_clim_1[v_name][
                            :
                        ].astype("double")
                        other_dat[seas][v_name + "_fut_min"] = var_clim_2[v_name][
                            :
                        ].astype("double")
                        other_dat[seas][v_name + "_pres_max"] = var_clim_1[v_name][
                            :
                        ].astype("double")
                        other_dat[seas][v_name + "_fut_max"] = var_clim_2[v_name][
                            :
                        ].astype("double")
                    else:
                        other_dat[seas][v_name + "_pres_sum"] += var_clim_1[v_name][
                            :
                        ].astype("double")
                        other_dat[seas][v_name + "_fut_sum"] += var_clim_2[v_name][
                            :
                        ].astype("double")
                        other_dat[seas][v_name + "_pres_min"] = np.minimum(
                            other_dat[seas][v_name + "_pres_min"],
                            var_clim_1[v_name][:].astype("double"),
                        )
                        other_dat[seas][v_name + "_fut_min"] = np.minimum(
                            other_dat[seas][v_name + "_fut_min"],
                            var_clim_2[v_name][:].astype("double"),
                        )
                        other_dat[seas][v_name + "_pres_max"] = np.minimum(
                            other_dat[seas][v_name + "_pres_max"],
                            var_clim_1[v_name][:].astype("double"),
                        )
                        other_dat[seas][v_name + "_fut_max"] = np.minimum(
                            other_dat[seas][v_name + "_fut_max"],
                            var_clim_2[v_name][:].astype("double"),
                        )

                # Variables
                # sum the values and their squares.
                #   First iteration copy them, and following iterations add
                #   Note that we are working in variance, rather than standard deviation, so need to convert the climatology std files into variances... we do this a this stage, but squaring var_clim_stdev_1  etc.
                for var in var_mat:
                    if ei == 0:
                        proc_dat[seas][var]["pres_mean_sum"] = var_clim_1[var][
                            0, :, :
                        ].astype("double")
                        proc_dat[seas][var]["fut_mean_sum"] = var_clim_2[var][
                            0, :, :
                        ].astype("double")
                        proc_dat[seas][var]["diff_mean_sum"] = var_clim_2[var][
                            0, :, :
                        ].astype("double") - var_clim_1[var][0, :, :].astype("double")
                        proc_dat[seas][var]["pres_mean_ss"] = (
                            var_clim_1[var][0, :, :].astype("double")
                        ) ** 2
                        proc_dat[seas][var]["fut_mean_ss"] = (
                            var_clim_2[var][0, :, :].astype("double")
                        ) ** 2
                        proc_dat[seas][var]["diff_mean_ss"] = (
                            var_clim_2[var][0, :, :].astype("double")
                            - var_clim_1[var][0, :, :].astype("double")
                        ) ** 2
                        # we need the interannual variance of each ensmeble member, so at this step need to convert from Std Dev by squaring
                        proc_dat[seas][var]["pres_var_sum"] = (
                            var_clim_stdev_1[var][0, :, :].astype("double") ** 2
                        )
                        proc_dat[seas][var]["fut_var_sum"] = (
                            var_clim_stdev_2[var][0, :, :].astype("double") ** 2
                        )
                        proc_dat[seas][var]["diff_var_sum"] = (
                            var_clim_stdev_2[var][0, :, :].astype("double") ** 2
                        ) - (var_clim_stdev_1[var][0, :, :].astype("double") ** 2)
                    else:
                        proc_dat[seas][var]["pres_mean_sum"] += var_clim_1[var][
                            0, :, :
                        ].astype("double")
                        proc_dat[seas][var]["fut_mean_sum"] += var_clim_2[var][
                            0, :, :
                        ].astype("double")
                        proc_dat[seas][var]["diff_mean_sum"] += var_clim_2[var][
                            0, :, :
                        ].astype("double") - var_clim_1[var][0, :, :].astype("double")
                        proc_dat[seas][var]["pres_mean_ss"] += (
                            var_clim_1[var][0, :, :].astype("double")
                        ) ** 2
                        proc_dat[seas][var]["fut_mean_ss"] += (
                            var_clim_2[var][0, :, :].astype("double")
                        ) ** 2
                        proc_dat[seas][var]["diff_mean_ss"] += (
                            var_clim_2[var][0, :, :].astype("double")
                            - var_clim_1[var][0, :, :].astype("double")
                        ) ** 2
                        # we need the interannual variance of each ensmeble member, so at this step need to convert from Std Dev by squaring
                        proc_dat[seas][var]["pres_var_sum"] += (
                            var_clim_stdev_1[var][0, :, :].astype("double") ** 2
                        )
                        proc_dat[seas][var]["fut_var_sum"] += (
                            var_clim_stdev_2[var][0, :, :].astype("double") ** 2
                        )
                        proc_dat[seas][var]["diff_var_sum"] += (
                            var_clim_stdev_2[var][0, :, :].astype("double") ** 2
                        ) - (var_clim_stdev_1[var][0, :, :].astype("double") ** 2)
                # increment ensemble counter
                proc_cnt += 1

                if Test:
                    Test_mean[seas][ei] = var_clim_1[Testvar][0, 120, 120]
                    Test_std[seas][ei] = var_clim_stdev_1[Testvar][0, 120, 120]

                # close input files
                rootgrp_clim_1.close()
                rootgrp_clim_stdev_1.close()
                rootgrp_clim_2.close()
                rootgrp_clim_stdev_2.close()

        # process ensemble statistics, from sum, and sum of square
        for seas in date_name_mat:
            # calc the mean for other variables (time, lon, lat (even time_bounds))
            for v_name in v_name_mat:
                if v_name in var_mat:
                    continue
                other_dat[seas][v_name + "_pres_mean"] = (
                    other_dat[seas][v_name + "_pres_sum"] / proc_cnt
                )
                other_dat[seas][v_name + "_fut_mean"] = (
                    other_dat[seas][v_name + "_fut_sum"] / proc_cnt
                )
            # Variables
            #   Calc the ens_mean, int_var, ens_var (ens_std inferred later)
            # for var, ncvar in zip(var_mat,ncvar_mat):
            for var in var_mat:
                for pername, perlab in zip(
                    ["pres", "fut", "diff"], [yrmat_1_str, yrmat_2_str, diff_yrmat_str]
                ):

                    # Ens_mean -
                    proj_dat[seas][var]["%s_ensmean" % pername] = (
                        proc_dat[seas][var]["%s_mean_sum" % pername] / proc_cnt
                    )

                    # int_var
                    #   note we have already converted the clim std into variance when we read the climatology files in.
                    proj_dat[seas][var]["%s_intvar" % pername] = (
                        proc_dat[seas][var]["%s_var_sum" % pername] / proc_cnt
                    )

                    # ens_var
                    proj_dat[seas][var]["%s_ensvar" % pername] = (
                        proc_dat[seas][var]["%s_mean_ss" % pername] / proc_cnt
                        - proj_dat[seas][var]["%s_ensmean" % pername] ** 2
                    )

                    # ens_stddev
                    proj_dat[seas][var]["%s_ensstd" % pername] = np.sqrt(
                        proj_dat[seas][var]["%s_ensvar" % pername]
                    )

        if Test:
            pdb.set_trace()
            Test_mean["dec"].mean(), proj_dat["dec"][Testvar]["pres_ensmean"][120, 120]

            Test_mean["dec"].var(), proj_dat["dec"][Testvar]["pres_ensvar"][120, 120]
            ((Test_std["dec"]) ** 2).mean(), proj_dat["dec"][Testvar]["pres_intvar"][
                120, 120
            ]

        # Write output files.
        # loop through seasons, and then pres, fut and diff.

        for seas in date_name_mat:
            for pername, perlab in zip(
                ["pres", "fut", "diff"], [yrmat_1_str, yrmat_2_str, diff_yrmat_str]
            ):
                print(datetime.now(), " Writing files: ", seas, perlab, grid_val)

                # output file name.
                tmpfname_out = (
                    "%sNWSClim_NWSPPE_EnsStats_clim_%s_%s_grid%s_stats.nc"
                    % (path_out, seas, perlab, grid_val)
                )

                # open file, and add dimensions, and variables.
                rootgrp_out = Dataset(tmpfname_out, "w", format="NETCDF4")

                # Add dimensions
                rootgrp_out.createDimension(lat_dim, 375)
                rootgrp_out.createDimension(lon_dim, 297)
                rootgrp_out.createDimension(time_dim, None)
                rootgrp_out.createDimension(bnd_dim, 2)

                # Add time variables and attributes
                time_var = rootgrp_out.createVariable("time", "f8", (time_dim))
                timebnds_var = rootgrp_out.createVariable(
                    "time_bounds", "f8", (time_dim, bnd_dim)
                )

                time_var.setncattr("long_name", "Time axis")
                time_var.setncattr("calendar", "360_day")
                time_var.setncattr("units", "seconds since 1950-01-01 00:00:00")
                time_var.setncattr("time_origin", "1950-01-01 00:00:00")
                time_var.setncattr("bounds", "time_bounds")

                # add time, and the time bounds.
                if pername == "diff":
                    time_var[:] = (
                        other_dat[seas]["time_fut_mean"][:]
                        + other_dat[seas]["time_pres_mean"][:]
                    ) / 2
                else:
                    time_var[:] = other_dat[seas]["time_%s_mean" % pername][:]

                if pername == "diff":
                    timebnds_var[0, 0] = other_dat[seas]["time_bounds_pres_min"].min()
                    timebnds_var[0, 1] = other_dat[seas]["time_bounds_fut_max"].max()
                else:
                    timebnds_var[0, 0] = other_dat[seas][
                        "time_bounds_%s_min" % pername
                    ][:].min()
                    timebnds_var[0, 1] = other_dat[seas][
                        "time_bounds_%s_max" % pername
                    ][:].max()

                # Add statistic variables.
                nc_2d_var_dict = {}
                for var in var_mat:
                    for ens_stat in ens_stat_lst:
                        nc_2d_var_dict[var + ens_stat] = rootgrp_out.createVariable(
                            var + ens_stat,
                            "f4",
                            (time_dim, lat_dim, lon_dim),
                            fill_value=1.0e20,
                            zlib=compress_variables,
                        )

                for var in var_dict[grid_val]:
                    for ens_stat in ens_stat_lst:

                        tmplongname = ens_stat_long_name_format_dict[ens_stat] % (
                            var.upper(),
                            perlab,
                            long_name_dict[var],
                        )
                        tmpunit = ens_stat_units_format_dict[ens_stat] % unit_dict[var]

                        nc_2d_var_dict[var + ens_stat].setncattr(
                            "long_name", tmplongname
                        )
                        # No Standard_names for Ens Stats
                        nc_2d_var_dict[var + ens_stat].setncattr("units", tmpunit)
                        nc_2d_var_dict[var + ens_stat].setncattr(
                            "online_operation", "average"
                        )
                        nc_2d_var_dict[var + ens_stat].setncattr(
                            "interval_write", "1 month"
                        )
                        nc_2d_var_dict[var + ens_stat].setncattr(
                            "interval_operation", "1 time-step"
                        )
                        nc_2d_var_dict[var + ens_stat].setncattr(
                            "cell_methods", ens_stat_cell_methods_dict[ens_stat]
                        )

                # Add data to variables.
                for var in var_dict[grid_val]:
                    for ens_stat in ens_stat_lst:
                        nc_2d_var_dict[var + ens_stat][0, :, :] = proj_dat[seas][var][
                            pername + ens_stat
                        ][:, :]

                # Add lat and lon.
                lon_var = rootgrp_out.createVariable("lon", "f4", (lon_dim))
                lat_var = rootgrp_out.createVariable("lat", "f4", (lat_dim))

                lon_var.setncattr("long_name", "Longitude")
                lon_var.setncattr("standard_name", "longitude")
                lon_var.setncattr("units", "degrees_east")
                lon_var.setncattr("nav_model", "grid_%s" % grid_val)
                lat_var.setncattr("long_name", "Latitude")
                lat_var.setncattr("standard_name", "latitude")
                lat_var.setncattr("units", "degrees_north")
                lat_var.setncattr("nav_model", "grid_%s" % grid_val)
                lon_var[:] = lonlat_out_nc_dict["lon_" + grid_val]
                lat_var[:] = lonlat_out_nc_dict["lat_" + grid_val]

                # Add global attributes
                out_source = out_source_NWSPPE

                rootgrp_out.setncattr("Conventions", out_Conventions)
                rootgrp_out.setncattr("institution", out_institution)
                rootgrp_out.setncattr("title", out_title)
                rootgrp_out.setncattr("source", out_source)

                rootgrp_out.setncattr("history", out_history)
                rootgrp_out.setncattr("references", out_references)
                rootgrp_out.setncattr(
                    "variance_separation_references", out_variance_separation_references
                )

                # close files.
                rootgrp_out.setncattr("comments", out_comment_EnsStat)
                rootgrp_out.setncattr("comments_grid", out_comment_grid % (grid_val))
                rootgrp_out.close()


def main():

    run_CEDA_monthly()
    run_CEDA_ens_climatologies()
    run_CEDA_ens_stats()
    run_CEDA_regional_means()

    pdb.set_trace()


if __name__ == "__main__":
    main()
