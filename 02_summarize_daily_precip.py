#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 08:26:30 2021

@author: ay394
"""

import os

import datetime as dt

import pandas as pd
import numpy as np

# Set Working Directory
wdir = "/Volumes/GoogleDrive/My Drive/Young_evapotranspiration_phenology_analysis"

# Load in metadata for each study site
os.chdir(wdir + "/data/ancillary_data")
phenocam_flux_metadata_table = pd.read_csv("pheno_flux_sites_to_use.csv")

MAP = phenocam_flux_metadata_table.MAP
sites = phenocam_flux_metadata_table.fluxsite
start_date = pd.to_datetime(phenocam_flux_metadata_table.start_date)
end_date = pd.to_datetime(phenocam_flux_metadata_table.end_date)

var_names = pd.read_csv("variables_to_import_for_fluxsites.csv")
precip_var_names = var_names.iloc[0:,var_names.columns == "precip"]

AMF_file_names = sorted(os.listdir(wdir + "/data/raw_data/ameriflux/BASE"),
                        key = str.casefold)

for i in range(0,len(sites)):
    
    os.chdir(wdir + "/data/raw_data/ameriflux/BASE")

    fluxdat = pd.read_csv(AMF_file_names[i],comment = '#')
    fluxdat = fluxdat.replace(to_replace = -9999,value = np.nan)
    
    os.chdir(wdir + '/data/raw_data/daymet')
    precip_daymet = pd.read_csv(sites[i] + '.csv',skiprows=7)
    precip_daymet['date'] = pd.to_datetime(precip_daymet['year']*1000 + precip_daymet['yday'],\
                                           format = '%Y%j')
    precip_daymet = precip_daymet.rename(columns={'prcp (mm/day)':'precip'},\
                                         copy=False)
    precip_daymet = precip_daymet[['date','precip']]
    
    precip = fluxdat[precip_var_names.precip[i]]
    
    datetime_start = pd.to_datetime(fluxdat.TIMESTAMP_START,format='%Y%m%d%H%M')
    datetime_end   = pd.to_datetime(fluxdat.TIMESTAMP_END,format='%Y%m%d%H%M')
    
    dn = datetime_start.map(dt.datetime.toordinal)
    
    nan_bool = np.isnan(precip)
    
    new_precip_df = pd.DataFrame({'dt':datetime_start,
                                  'dn':dn,
                                  'precip':precip,
                                  'nan_bool':nan_bool})
    
    
    
        
    precip_daily = new_precip_df.groupby('dn')['precip'].sum()
    precip_daily.loc[np.isnan(precip_daily)] = -9999.0
    
    nan_bool_daily = new_precip_df.groupby('dn')['nan_bool'].sum()
    nan_bool_daily[np.isnan(nan_bool_daily)] = -9999.0
    
    date_values = precip_daily.index.map(dt.datetime.fromordinal)
    
    precip_flux_df = pd.DataFrame({'date':date_values, \
                                   'precip':precip_daily.values, \
                                   'n_missing':nan_bool_daily.values})
        
    precip_df = precip_flux_df.merge(precip_daymet,\
                                     on='date',\
                                     suffixes=('_amf','_daymet'))
        
    if sites[i] == 'US-KFS':
        
        os.chdir(wdir + '/data/raw_data/noaa_ncei')
        kfs_precip = pd.read_csv('kfs_precip_data.csv',parse_dates = [0])
        precip_df = precip_df.merge(kfs_precip,\
                                    on = 'date')
        export_df = precip_df[['date','precip_amf','precip_daymet','precip_noaa']]

    else:
        
        export_df = precip_df[['date','precip_amf','precip_daymet']]
        
    
    os.chdir(wdir + "/results/flux_data/precip_data")
    export_df_fn = sites[i] + "_precip.csv"
    export_df.to_csv(export_df_fn,index=False)
  
# End of script --------------------------------------------------------------------