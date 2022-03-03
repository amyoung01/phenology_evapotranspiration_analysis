#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np

# Set Working Directory
wdir = "/Volumes/GoogleDrive/My Drive/W/projects/Young_evapotranspiration_phenology_analysis"

# Load in metadata for each study site
os.chdir(wdir + "/data/ancillary_data")
phenocam_flux_metadata_table = pd.read_csv("pheno_flux_sites_to_use.csv")

sites = phenocam_flux_metadata_table.fluxsite

var_names = pd.read_csv("variables_to_import_for_fluxsites.csv")
vars_of_interest = var_names.columns[1:len(var_names.columns)]

AMF_file_names = sorted(os.listdir(wdir + "/data/raw_data/ameriflux/BASE"),key = str.casefold)

for i in range(0,len(sites)):
    
    os.chdir(wdir + "/data/raw_data/ameriflux/BASE")
    
    fluxdat = pd.read_csv(AMF_file_names[i],skiprows=[0,1])
    fluxdat = fluxdat.replace(to_replace=-9999.0,value=np.nan)
    
    datetime_start = pd.to_datetime(fluxdat.TIMESTAMP_START,format='%Y%m%d%H%M')
    datetime_end   = pd.to_datetime(fluxdat.TIMESTAMP_END,format='%Y%m%d%H%M')
    
    dt = {'datetime_start':datetime_start,
          'datetime_end':datetime_end}
    dt_df = pd.DataFrame(data = dt)    
    
    raw_data_var_names = fluxdat.columns    
     
    for p in range(0,len(vars_of_interest)):
        
        if vars_of_interest[p] == "precip":
            continue            
        
        column_id = p + 1
        varname_to_import = var_names.iloc[i,column_id]
        
        if np.isreal(varname_to_import):
            
            fluxdat[vars_of_interest[p]] = np.nan
                
        else:
            
            fluxdat = fluxdat.rename(columns = {varname_to_import:vars_of_interest[p]},\
                                     copy = False)
    
    fluxdat_subset_df = fluxdat[vars_of_interest[vars_of_interest != "precip"]]
    fluxdat_subset_df = pd.concat([dt_df,fluxdat_subset_df],axis=1)
    fluxdat_subset_df = fluxdat_subset_df.replace(to_replace=np.nan,value=-9999.0)
        
    os.chdir(wdir + "/results/flux_data/halfhour")       
    export_fn = sites[i] + ".csv"
    
    pd.DataFrame.to_csv(fluxdat_subset_df,\
                        index=False,\
                        path_or_buf=export_fn)

# End of script ----------------------------------------------------------------