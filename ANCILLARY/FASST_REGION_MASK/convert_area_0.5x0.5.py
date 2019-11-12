#===================================================================
# Python script to convert embase.sav into a csv file
# Author: Jerome Hilaire
# Email : hilaire@pik-potsdam.de
# Date  : 03-03-2016
#===================================================================
# Load required libraries and functions
from scipy.io.idl import readsav
import pandas as pd
import csv

# Readin IDL .sav file
data = readsav("EMBASE.SAV")

# Get names (header name + species)
names   = list(data.embase.dtype.names)
nspecs  = len(names)
hdname  = names[0]
species = names[1:nspecs]
nspecs  = nspecs -1

# Get sectors
sectors = data.embase[hdname][0].tolist()
nsect   = len(sectors)
sectors = [x.decode('utf-8') for x in sectors]  # transform from binary to string

# Reshape data into a dictionary first ...
d={}
d[hdname] = sectors
for sp in species:
  d[sp]=data.embase[sp][0].tolist()

# ... and then a dataframe
df = pd.DataFrame(d, columns=names)

# Save dataframe in csv file
df.to_csv("embase.csv", index=0)

