#!/usr/bin/env python3

"""
01_standardized_halfhour_flux_data.py

This script reads in the BASE AmeriFlux datasets for each site and standardizes
so that the variable names are the same across all sites. It does this for all
variables except precip which is handled in a different script.
"""

import os
import pandas as pd

# Set Working Directory
wdir = "/Volumes/GoogleDrive/My Drive/W/projects/phenology_evapotranspiration_analysis"

# Load in metadata and site info for each study site
os.chdir(os.path.join(wdir,"data/ancillary_data"))
phenocam_flux_metadata_table = pd.read_csv("pheno_flux_sites_to_use.csv")

# Get list of site names (in alphabetical order)
sites = phenocam_flux_metadata_table.fluxsite

# Read in table that stores variables names
var_names = pd.read_csv("variables_to_import_for_fluxsites.csv",
                    index_col=0)
vars_of_interest = var_names.columns

# Get list of file names to read/import for each site
fname = phenocam_flux_metadata_table.AMF_filename

for i in range(0,len(sites)):
    
    # Read in BASE AmeriFlux data, skipping first two rows
    os.chdir(os.path.join(wdir,"data/raw_data/ameriflux/BASE"))
    fluxdat = pd.read_csv(fname[i],comment = '#')
    fluxdat = fluxdat.replace(to_replace=-9999.0,value=pd.NA)
    
    # Reformat start and end datetimes
    datetime_start = pd.to_datetime(fluxdat.TIMESTAMP_START,
                                format='%Y%m%d%H%M')
    datetime_end   = pd.to_datetime(fluxdat.TIMESTAMP_END,
                                format='%Y%m%d%H%M')
    
    # Create separate dataframe of just datetime values
    dt = {'datetime_start':datetime_start,
          'datetime_end':datetime_end}
    dt_df = pd.DataFrame(data = dt)    
    
    # Get original list of variable names for BASE file
    raw_data_var_names = fluxdat.columns    
    
    # For each variable ...
    for p in range(0,len(vars_of_interest)):
        
        # Skip precip for now
        if vars_of_interest[p] == "precip":
            continue            
        
        # Variable name to standardize
        varname_to_import = var_names.iloc[i,p]
        
        # If data is not available in BASE file still provide column
        # in output dataframe but make all values NA
        if pd.isna(varname_to_import):
            
            fluxdat[vars_of_interest[p]] = pd.NA
                
        else:
            
            fluxdat = fluxdat.rename(columns = {varname_to_import:vars_of_interest[p]},
                                    copy = False)
    
    # Prepare dataframe for export by adding datetime values back in using concat
    # and removing precip column. Set all NA values to -9999.0
    fluxdat_subset_df = fluxdat[vars_of_interest[vars_of_interest != "precip"]]
    fluxdat_subset_df = pd.concat([dt_df,fluxdat_subset_df],axis=1)
    fluxdat_subset_df = fluxdat_subset_df.replace(to_replace=pd.NA,
                                                value=-9999.0)
        
    # Set export filename with directory
    export_fname = os.path.join(wdir,
                                "results/flux_data/halfhour",sites[i] + ".csv")
    
    # Export as CSV file
    pd.DataFrame.to_csv(fluxdat_subset_df,
                        index=False,
                        path_or_buf=export_fname)

# End of script ----------------------------------------------------------------