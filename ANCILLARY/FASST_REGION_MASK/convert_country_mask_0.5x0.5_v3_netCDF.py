#===================================================================
# Python script to convert country_mask_0.5x0.5_v3.sav into a csv 
# file
# 
# Author: Jerome Hilaire
# Email : hilaire@pik-potsdam.de
# Date  : 10-05-2016
#===================================================================
# Load required libraries and functions
from scipy.io.idl import readsav
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import csv
import time
import netCDF4

# Readin IDL .sav file
data = readsav("country_mask_0.5x0.5_v3.sav")

# Get keys
keys = data.keys()

#TODO: Add empty entries in iso_code and iso_country_name
country_iso      = pd.Series([x.decode('utf-8') for x in data["country_iso"]])
country_reg      = pd.Series([x.decode('utf-8') for x in data["country_region"]])
country          = pd.concat([country_iso, country_reg], axis=1, keys=['iso_code', 'fasst_reg_id'])
iso_code         = pd.Series([x.decode('utf-8') for x in data["iso_code"][0:234]])
iso_code         = pd.concat([iso_code, pd.Series(['','',''], index=[234,235,236])])
iso_country_name = pd.Series([x.decode('utf-8') for x in data["iso_country_name"][0:234]])
iso_country_name = pd.concat([iso_country_name, pd.Series(['','',''], index=[234,235,236])])
iso              = pd.concat([iso_code, iso_country_name], axis=1, keys=['iso_code', 'iso_country_name'])
df_mapping       = pd.merge(country, iso, on='iso_code', how='left')
df_mapping['fasst_reg_id'] = [x.strip() for x in df_mapping['fasst_reg_id']]

country_poles    = pd.concat([pd.Series([str(x) for x in np.arange(56)+1]), pd.Series([x.decode('utf-8') for x in data["country_poles"]])], axis=1, keys=['fasst_reg_id', 'fasst_reg_name'])
df_mapping       = pd.merge(df_mapping, country_poles, on='fasst_reg_id', how='left')

# Save dataframe in csv file
df_mapping.to_csv("country_mask_0.5x0.5_v3_mapping.csv", index=0)

# Process country mask array
# Write NC file
# Open a new NetCDF file to write the data to. For format, you can choose from
# 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
w_nc_fid = Dataset('country_mask_0.5x0.5_v3.nc', 'w', format='NETCDF4_CLASSIC')
# Global Attributes
w_nc_fid.description = 'Country mask'
w_nc_fid.history = 'Created ' + time.ctime(time.time())
w_nc_fid.source = 'JRC. Obtained from Rita van Dingenen.'
                       
# Create dimensions
lon     = w_nc_fid.createDimension('lon', 720)
lat     = w_nc_fid.createDimension('lat', 360)
country = w_nc_fid.createDimension('country', 57)

# Create coordinat  e variables for 4-dimensions
countries  = w_nc_fid.createVariable('country', 'i4', ('country',))
latitudes  = w_nc_fid.createVariable('lat',     'f4', ('lat',))
longitudes = w_nc_fid.createVariable('lon',     'f4', ('lon',)) 

# Create country_mask variable
country_mask = w_nc_fid.createVariable('country_mask', np.int32, ('country','lat','lon')) 

# Variable Attributes
latitudes.units    = 'degree_north'
longitudes.units   = 'degree_east'
country_mask.units = 'none'

# Write values
longitudes[:] = np.arange(-180,180,0.5)
latitudes[:]  = np.arange(-90,90,0.5)
countries[:]  = np.arange(0,57,1)

country_mask[:,:,:] = data["country_mask"]

# Close the nc file
w_nc_fid.close()  
