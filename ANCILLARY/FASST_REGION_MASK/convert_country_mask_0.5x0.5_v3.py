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
import numpy as np
import pandas as pd
import csv

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
arr = data["country_mask"]

# Case: 0
df0 = pd.DataFrame(arr[0].squeeze().byteswap().newbyteorder())
df0['lat']  = np.arange(360)
df0 = pd.melt(df0, id_vars=['lat'], var_name="lon", value_name="1")

for kc in np.arange(56)+1:
  tmp = pd.DataFrame(arr[1].squeeze().byteswap().newbyteorder())
  tmp['lat']  = np.arange(360)
  tmp = pd.melt(tmp, id_vars=['lat'], var_name="lon", value_name=str(kc))
  tmp.pop('lat')
  tmp.pop('lon')
  df0 = pd.concat([df0, tmp], axis=1)
  
df0 = pd.melt(df0, id_vars=['lat', 'lon'], var_name="fasst_region", value_name='value'))

# Save dataframe in csv file
df0.to_csv("country_mask_0.5x0.5_v3.csv", index=0)

