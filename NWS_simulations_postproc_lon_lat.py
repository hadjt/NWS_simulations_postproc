import numpy as np

nlat = 375
nlon = 297

lat_m_coef = {}
lat_c_coef = {}
lon_m_coef = {}
lon_c_coef = {}

lat_m_coef['T'],lat_c_coef['T']=0.0666699998, 40.0666700263
lat_m_coef['U'],lat_c_coef['U']=0.0666699998, 40.0666700263
lat_m_coef['V'],lat_c_coef['V']=0.0666699973, 40.1000069691
lon_m_coef['T'],lon_c_coef['T'] =0.1111100001,-19.8888900075
lon_m_coef['U'],lon_c_coef['U'] =0.1111099998,-19.8333346811
lon_m_coef['V'],lon_c_coef['V'] =0.1111100001,-19.8888900075


lonlat_out_nc_dict = {}
for gi,grid_val in enumerate(["T","U","V"]):
    lonlat_out_nc_dict['lon_' + grid_val] = np.arange(nlon)*lon_m_coef[grid_val] + lon_c_coef[grid_val]

for gi,grid_val in enumerate(["T","U","V"]):
    lonlat_out_nc_dict['lat_' + grid_val] = np.arange(nlat)*lat_m_coef[grid_val] + lat_c_coef[grid_val]
