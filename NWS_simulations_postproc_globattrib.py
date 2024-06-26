# -*- coding: utf-8 -*-
out_Conventions = "CF-1.8"

out_institution = "Met Office Hadley Centre, Exeter, UK."

out_title = "Marine climate projections for the North West European Shelf Seas.\nA GCM Perturbed Parameter Ensemble and present-day control simulation downscaled with a shelf seas model.\n"

out_history = "Model output post processed by NWS_simulations_postproc.py: https://github.com/hadjt/NWS_simulations_postproc"

out_references = "Tinker et al. 2023, 21st century marine climate projections for the NW European Shelf Seas based on a perturbed parameter ensemble, in prep."
out_comment_monthly_mean = "These are monthly mean values calculated from every timestep, as output by NEMO shelf version %s."

out_comment_regmean = "Regional mean time series are output by NEMO shelf version %s, on the Wakelin et al. (2012) region mask (see region_refererence). NEMO reads in the region mask, and every time step it averages the variable within each mask region. The regional mean is then averaged over the month, and output by NEMO. This file includes the regional mean time-series for a selection of variables (as RegAveSST, RegAveSSS etc. for the regional mean SST and SSS respectively), for each month, and each of the 14 regions. There are variables for the region id (reg_id) and the number of grid boxes within each region (cnt). The region mask included (mask, with the associated longitude and latitude variables). This methodology is described by Tinker et al. (2019; see region_methology_reference)."


out_region_names = "['Shelf','Southern North Sea','Central North Sea','Northern North Sea','English Channel','Skagerrak/Kattegat','Norwegian Trench','Shetland Shelf','Irish Shelf','Irish Sea','Celtic Sea','Armorican Shelf','NE Atlantic (S)','NE Atlantic (N)']"

out_region_notes = "The region Shelf is the combination of Southern North Sea; Central North Sea; Northern North Sea; English Channel; Shetland Shelf; Irish Shelf; Irish Sea; and Celtic Sea. It does not include the Norwegian Trench; Skagerrak/Kattegat; Armorican Shelf; NE Atlantic (S); NE Atlantic (N)"

out_region_refererence = "Adapted from Wakelin, S. L., Holt, J., Blackford, J., Allen, I., Butenschön, M., and Artioli, Y.: Modeling the carbon fluxes of the northwest European continental shelf: Validation and budgets, 117, C05020, https://doi.org/10.1029/2011JC007402, 2012.".encode(
    "utf8"
)

out_region_methology_reference = "Tinker, J., Renshaw, R., Barciela, R., and Wood, R.: Regional mean time series for the Northwest European Shelf seas. In: Copernicus Marine Service Ocean State Report, Issue 3, J. Oper. Oceanogr., 12, s26-s30, https://doi.org/10.1080/1755876X.2019.1633075, 2019."

out_source_PDCtrl = "Underlying GCM: HadGEM3 GC3.05\nDownscaling shelf seas model: NEMO shelf version 3.6.\nPresent-Day Control Simulation representing the year 2000.\nModel simulations developed and run by Dr. Jonathan Tinker. The shelf climate version of NEMO 3.6 used here is available on GitHub: https://github.com/hadjt/NEMO_3.6_CO6_shelf_climate"

out_source_NWSPPE = "Underlying GCM: HadGEM3 GC3.05\nDownscaling shelf seas model: NEMO shelf version 4.0.4.\nRepresentative Concentration Pathway: RCP8.5.\nModel simulations developed and run by Dr. Jonathan Tinker. The shelf climate version of NEMO 4.0.4 used here is available on GitHub: https://github.com/hadjt/NEMO_4.0.4_CO9_shelf_climate"


out_comment_grid = 'Note that NEMO uses the Arakawa "C" grid, where the T, U and V grids are offset. Most variables are on the T grid, apart from the Eastward and Northward components of the ocean barotropic current, which are on the U and V grids respectively. We also provide the barotropic current speed on T grid, where we transform the U and V velocity components onto the T grid, before calculating their magnitude. We therefore separate the variables on the T, U and V grids into separate files. This file contains variables on the %s grid.'


out_comment_EnsStat = 'These ensemble statistics are calculated from the climatological mean and standard deviation for each NWSPPE ensemble member, which in turn are based on the monthly mean values calculated from every timestep, as output by NEMO shelf version 4.0.4. For a given climatological period (e.g. 2000-2019, 2079-2098), and month, season or year, the climatological mean and standard deviation of all 12 NWSPPE ensemble members are loaded. The mean of these climatological means gives the Ensemble Mean (ensmean). The variance of the climatological means gives the Ensemble Variance (ensvar), and its square-root gives the Ensemble Standard Deviation (ensstd). The mean of climatological standard deviation squared (i.e. the climatological variance) gives the Interannual Variance (intvar). Each variable (SST, SSS, NBT, etc.) has 4 statistic associated with it, and the NetCDF variable names join the variable name (e.g. SST) with the statistic name (e.g. ensmean), separated by an underscore. For example the ensemble mean SST is SST_ensmean and the NBT interannual variability is NBT_intvar. This is all captured in the "cell_methods", where "realization" denotes the ensemble members. The Ensemble Variance and Interannual Variance can be combined into the Total Variance following Tinker et al. 2016; see variance_separation_references).'

out_comment_EnsStat_diff = (out_comment_EnsStat + '\nTo calculate the absolute difference between the two periods, we subtract the future ensemble statistics from those of the present day - these are stored as ensmean, ensvar, ensstd and intvar. This is useful to describe how the behaviour of the ensemble changes. To use these statistics for projections, we also provide the likely projected change (projensmean, based on the anomaly ensemble mean), and the associated estimate of the uncertainty (projensstd, based on the anomaly ensemble standard deviation. To calculate the projensmean and projensstd, we calculate the anomaly ensemble, relative to the present-day period: for each ensemble member, we remove the future climatological mean from the present-day climatological mean before calculating ensemble mean and standard deviation across the ensemble. Projensmean is effectively ensmean (any differences are within the noise of the data, and are related to truncation errors), whereas projensstd is a different statistic to ensstd: projensstd cannot be negative, whereas ensstd may be negative. When we discuss the southern North Sea projected change, in terms of x ± 2 std deviation, we take the regional mean of projensmean ± 2 x the regional mean of projensstd.').encode("utf8")

out_variance_separation_references = "Tinker, J., Lowe, J., Pardaens, A., Holt, J., and Barciela, R.: Uncertainty in climate projections for the 21st century northwest European shelf seas, Prog. Oceanogr., 148, 56-73, https://doi.org/10.1016/j.pocean.2016.09.003, 2016."
