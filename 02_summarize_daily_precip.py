#!/usr/bin/env python3

"""
02_summarize_daily_precip.py

This script summarizes daily precip for each site for the following data 
sources:
    1. AmeriFlux halfhour data
    2. Daymet data
    3. And for US-KFS only NOAA weather station data

There appeared to be significant data gaps for US-KFS from AmeriFlux so we uesed 
NOAA weather station data for only this site. All other sites used AmeriFlux precip
data. Daymet was collected and used to evaluate relative completeness of AmeriFlux
precip time series.

"""

import os
import datetime as dt
import pandas as pd

# Set Working Directory
wdir = "/Volumes/GoogleDrive/My Drive/W/projects/phenology_evapotranspiration_analysis"

# Load in metadata for each study site
os.chdir(os.path.join(wdir,"data/ancillary_data"))
phenocam_flux_metadata_table = pd.read_csv("pheno_flux_sites_to_use.csv")

sites = phenocam_flux_metadata_table.fluxsite

var_names = pd.read_csv("variables_to_import_for_fluxsites.csv")
precip_var_names = var_names.iloc[0:,var_names.columns == "precip"]

fname = phenocam_flux_metadata_table.AMF_filename

for i in range(0,len(sites)):
    
    # -------------------------------------------------------------------------
    # Read in AmeriFlux and Daymet datasets
    # -------------------------------------------------------------------------

    # Read in AmeriFlux BASE file for a given site
    os.chdir(os.path.join(wdir,"data/raw_data/ameriflux/BASE"))
    fluxdat = pd.read_csv(fname[i],comment = '#')
    fluxdat = fluxdat.replace(to_replace = -9999,value = pd.NA)
    
    # Read in Daymet climate data for the grid cell where a given site 
    # is located.
    os.chdir(os.path.join(wdir,"data/raw_data/daymet"))
    precip_daymet = pd.read_csv(sites[i] + '.csv',skiprows=7)

    # -------------------------------------------------------------------------
    # Restructure dateframe with precip data to summarize at daily time step
    # and to have the same column names
    # -------------------------------------------------------------------------

    # Restructure AmeriFlux dataframe -----------------------------------------
    precip = fluxdat[precip_var_names.precip[i]]
    
    datetime_start = pd.to_datetime(fluxdat.TIMESTAMP_START,format='%Y%m%d%H%M')
    datetime_end   = pd.to_datetime(fluxdat.TIMESTAMP_END,format='%Y%m%d%H%M')
    
    # Get datenumber values to groupby
    dn = datetime_start.map(dt.datetime.toordinal)
    
    # Identify missing observations
    nan_bool = pd.isna(precip)
    
    # Create new dataframe to summarize precip data
    new_precip_df = pd.DataFrame({'dt':datetime_start,
                                  'dn':dn,
                                  'precip':precip,
                                  'nan_bool':nan_bool})
    
    # Sum daily precip data
    precip_daily = new_precip_df.groupby('dn')['precip'].sum()
    precip_daily.loc[pd.isna(precip_daily)] = -9999.0
    
    # Get number of missing precip values for each day
    nan_bool_daily = new_precip_df.groupby('dn')['nan_bool'].sum()
    nan_bool_daily[pd.isna(nan_bool_daily)] = -9999.0
    
    # Get datetime format values for date in final dataframe
    date_values = precip_daily.index.map(dt.datetime.fromordinal)
    
    precip_flux_df = pd.DataFrame({'date':date_values, \
                                   'precip':precip_daily.values, \
                                   'n_missing':nan_bool_daily.values})

    # Round precip values to 2 decimal places
    precip_flux_df["precip"] = precip_flux_df["precip"].astype("float").round(2)

    # Restructure daymet dataframe --------------------------------------------

    # Use year and doy values to get datetime values.
    yr_str = [str(c) for c in list(precip_daymet.year)]
    doy_str = ["%03d" % c for c in list(precip_daymet.yday)]
    date_str = list(map(lambda a,b: a+b, yr_str, doy_str))

    precip_daymet["date"] = pd.to_datetime(date_str,
                                           format = '%Y%j')
    precip_daymet = precip_daymet.rename(columns={'prcp (mm/day)':'precip'},\
                                         copy=False)
    precip_daymet = precip_daymet[['date','precip']]

    # Round precip values to 2 decimal places
    precip_daymet["precip"] = precip_daymet["precip"].astype("float").round(2)

    # -------------------------------------------------------------------------
    # Create dataframe that has both AmeriFlux and Daymet
    # -------------------------------------------------------------------------

    precip_df = precip_flux_df.merge(precip_daymet,\
                                     on='date',\
                                     suffixes=('_amf','_daymet'))
    
    # If site is US-KFS include NOAA weather station data as well.
    if sites[i] == 'US-KFS':
        
        os.chdir(wdir + '/data/raw_data/noaa_ncei')
        kfs_precip = pd.read_csv('kfs_precip_data.csv',parse_dates = [0])
        precip_df = precip_df.merge(kfs_precip,\
                                    on = 'date')
        export_df = precip_df[['date','precip_amf','precip_daymet','precip_noaa']]

    else:
        
        export_df = precip_df[['date','precip_amf','precip_daymet']]
        
    # Export/write daily precip file
    export_fname = os.path.join(wdir,
                                "results/flux_data/precip_data",
                                sites[i] + "_precip.csv")
    export_df.to_csv(export_fname,index=False)
  
# End of script --------------------------------------------------------------------