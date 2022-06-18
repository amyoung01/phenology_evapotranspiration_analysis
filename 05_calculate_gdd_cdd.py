#!/usr/bin/env python3

import os, sys
import numpy as np
import pandas as pd
import datetime as dt
from scipy.interpolate import interp1d

wdir = "/Volumes/GoogleDrive/My Drive/W/projects/phenology_evapotranspiration_analysis"

# Import custom lowess smoother function
sys.path.append(wdir + '/code/z_functions')
from lowess import lowess

os.chdir(os.path.join(wdir,"data/ancillary_data"))
phenoflux_metadata = pd.read_csv("pheno_flux_sites_to_use.csv")

sites = phenoflux_metadata.fluxsite
phenos = phenoflux_metadata.phenosite
vegtypes = phenoflux_metadata.vegtype

# Pre-determined parameters for GDD-CDD for each PFT
db_params = np.array([[12,24],[50,240],[150,500]]) # DBF        
en_params = np.array([[5,18],[60,250],[200,400]]) # ENF  
gr_params = np.array([[20,30],[50,180],[150,500]]) # GRA    
sh_params = db_params # OSH
sa_params = db_params # WSA

unique_veg = np.array(['DB','EN','GR','SA','SH'])
gdd_cdd_params = np.stack((db_params,en_params,gr_params,sa_params,sh_params))

for i in range(0,len(sites)):
    
    # Find which PFT corresponds to site[i]
    veg_id = unique_veg == vegtypes[i]
    
    # Load in halfhour fluxdata
    os.chdir(wdir + "/results/flux_data/halfhour")
    fluxdat_hh = pd.read_csv(sites[i] + '.csv',\
                             parse_dates = [0])
    fluxdat_hh = fluxdat_hh.replace(to_replace = -9999.0,value = np.nan)
    
    # Look at midday temperature data
    hr = fluxdat_hh.datetime_start.dt.hour
    hr_id = (hr >= 10) & (hr < 14) 
    fluxdat_hh = fluxdat_hh.loc[hr_id] 
    t_air_hh = fluxdat_hh[['datetime_start','t_air']]
    t_air_hh = t_air_hh.reset_index(drop = True)
    
    dn = t_air_hh.datetime_start.map(dt.datetime.toordinal)
    t_air_hh['dn'] = dn
    
    # Get midday average temperature data
    t_air_daily = t_air_hh.groupby('dn')['t_air'].mean()
    dates = t_air_daily.index.map(dt.datetime.fromordinal)
    
    t_air_daily = pd.DataFrame({'date':dates,\
                                't_air':t_air_daily.values})
    
    # Calculate GDD and CDD for each year
    gdd_cdd_i = pd.DataFrame()
        
    gdd_cdd_i['date'] = t_air_daily.date
    gdd_cdd_i['doy'] = dates.dayofyear
    gdd_cdd_i['t_air'] = t_air_daily.t_air
    
    x = np.asarray(range(0,gdd_cdd_i.shape[0]))
    y = gdd_cdd_i.t_air
    
    x0 = x
    x1 = x
    y0 = y
    y1 = y
    
    x0 = x0[~np.isnan(y0)]
    y0 = y0[~np.isnan(y0)]
    
    x1 = x1[np.isnan(y1)]
    x1 = x1[(x1 > np.min(x0)) & (x1 < np.max(x0))]
    
    smooth_y = lowess(x0,y0,0.01)
    f = interp1d(x0,smooth_y,kind = 'nearest')
    y_pred = f(x1)
    
    gdd_cdd_i.loc[x1,'t_air'] = y_pred
    
        
    os.chdir(wdir + "/results/flux_data/daily")
    fluxdat = pd.read_csv(sites[i] + '_daily_values.csv',\
                          parse_dates=[0])
    fluxdat = fluxdat.replace(to_replace=-9999.0,value=np.nan)
    dates = fluxdat.date
    
    fluxdat = fluxdat.merge(gdd_cdd_i,
                            how = 'left',
                            on = 'date',                           
                            copy = False)
    
    yr = np.unique(dates.dt.year)
    complete_yr_id = np.ones((fluxdat.shape[0],1),dtype = 'bool')    

    for y in range(0,len(yr)):
        
        yr_id = np.where(dates.dt.year == yr[y])[0]
        
        if len(yr_id) < 360:
          
           complete_yr_id[yr_id] = False
           
    fluxdat = fluxdat.loc[complete_yr_id]
    fluxdat = fluxdat.reset_index(drop = True)
            
    gdd = np.zeros((fluxdat.shape[0],1))
    cdd = np.zeros((fluxdat.shape[0],1))
    
    gdd_base = gdd_cdd_params[veg_id,0,0]
    cdd_base = gdd_cdd_params[veg_id,0,1]
    gdd_start = gdd_cdd_params[veg_id,1,0]
    cdd_start = gdd_cdd_params[veg_id,1,1]
    gdd_max = gdd_cdd_params[veg_id,2,0]
    cdd_max = gdd_cdd_params[veg_id,2,1]
    
    gdd_0 = np.maximum(fluxdat.t_air_y - gdd_base,0)
    cdd_0 = np.maximum(cdd_base - fluxdat.t_air_y,0)
    
    for k in range(1,len(gdd_0)):
        
        if fluxdat.doy[k] >= gdd_start: 
            gdd[k] = np.minimum(gdd[k-1] + gdd_0[k],gdd_max)
            
        if fluxdat.doy[k] >= cdd_start:
            cdd[k] = np.minimum(cdd[k-1] + cdd_0[k],cdd_max)
           
    
    fluxdat['gdd'] = gdd
    fluxdat['cdd'] = cdd
    fluxdat['gdd_cdd'] = gdd / gdd_max - cdd / cdd_max
    
    fluxdat = fluxdat.replace(to_replace=np.nan,value=-9999)
    fluxdat = fluxdat.drop(['doy','t_air_y'],axis=1)
    fluxdat = fluxdat.rename(columns={'t_air_x': 't_air'})
    fluxdat = fluxdat[[col for col in fluxdat.columns if col != 'to_remove_id' ] + ['to_remove_id'] ]

    os.chdir(wdir + '/results/flux_data/daily')
    pd.DataFrame.to_csv(fluxdat,sites[i] + '_daily_values.csv',\
                        index = False)
    
# End of script ---------------------------------------------------------------------