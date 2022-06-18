# ############################################################################
# Calculate daily summaries of flux data needed for path analysis
# ############################################################################

# Initialize Workspace 
rm(list = ls()) # Clear all data from workspace
graphics.off() # close all current figures and plots
cat("\14") # clear command console

# Load required packages
require(bigleaf) # Physical and Physiological Ecosystem Properties from Eddy
                 # Covariance (version 0.7.1)
require(lubridate) # Make Dealing with Dates a Little Easier (version 1.7.10)
require(zoo) # S3 Infrastructure for Regular and Irregular Time Series (version 1.8-8)
require(data.table) # Extension of 'data.frame' (version 1.14.0)

# Set working directory path
wdir <- "/Volumes/GoogleDrive/My Drive/W/projects/phenology_evapotranspiration_analysis"

# Load in metadata table with site-specific information
setwd(paste0(wdir,"/data/ancillary_data"))
phenoflux_metadata <- read.csv("pheno_flux_sites_to_use.csv")

sites <- phenoflux_metadata$fluxsite
phenos <- phenoflux_metadata$phenosite
emissivity <- phenoflux_metadata$emissivity # PFT specific emsissivity. 
                                            # Needed for radiometric surface
                                            # temp calculations

for (i in 1:length(sites)){
  
  # ############################################################################
  # Load in required datasets for site [i] 
  # ############################################################################
  
  # Load in half-hour (or hour) flux data
  setwd(paste0(wdir,"/results/flux_data/halfhour"))
  fluxdat <- read.csv(sprintf("%s.csv",sites[i]))
  fluxdat[fluxdat == -9999] <- NA
  
  # Load in daily precip summaries
  setwd(paste0(wdir,"/results/flux_data/precip_data"))
  precip <- read.csv(sprintf("%s_precip.csv",sites[i]))
  
  # Load in Gcc data from PhenoCam located at the flux tower site
  setwd(paste0(wdir,"/results/processed_phenocam_data/time_series"))
  phenodat <- read.csv(sprintf("%s_gcc_time_series.csv",phenos[i]))
  phenodat[phenodat == -9999] <- NA
  
  # ############################################################################
  # Summarize to daily averages from half-hour 
  # ############################################################################
  
  # First, get hour values for each 30 minute (or hour) observation using 
  # lubridate package to identify midday observations.
  dt <- strftime(fluxdat$datetime_start,format="%Y-%m-%d %H:%M:%S")
  
  hr <- lubridate::hour(dt)
  midday_bool <- hr >=  10 & hr < 14
  
  fluxdat_midday <- fluxdat[midday_bool,]
  
  # Next, get Tair > 0 and ustar > 0.2 filters
  ustar_bool <- fluxdat_midday$ustar > 0.2
  t_air_bool <- fluxdat_midday$t_air > 0.0
  
  # Use above filters to fill in rows that don't meet the above criteria with 
  # NAs
  to_keep_bool_HH <- ustar_bool & t_air_bool
  to_keep_bool_HH[is.na(to_keep_bool_HH)] <- FALSE
  fluxdat_midday[!to_keep_bool_HH,3:ncol(fluxdat_midday)] <- NA
  
  # Get date without times (hours and minutes) and add to data frame
  fluxdat_midday$date <- as.Date(fluxdat_midday$datetime_start)
  
  # Get daily averages for each variable using the aggregate function by date
  fluxdat_daily <- aggregate.data.frame(fluxdat_midday[,3:(ncol(fluxdat_midday)-1)],
                                        by = list(date = fluxdat_midday$date),
                                        FUN = mean,
                                        na.rm = TRUE,
                                        na.action = NULL)
  
  # Convert NaNs to NAs for consistency. R doesn't make this easy for some reason
  nan_id <- as.data.frame(lapply(fluxdat_daily[,2:ncol(fluxdat_daily)],is.nan))
  fluxdat_daily[,2:ncol(fluxdat_daily)][as.matrix(nan_id)] <- NA
  
  # ############################################################################
  # Precipitation calculations 
  # ############################################################################
  
  # Now calculate cumulative precipitation for 2 and 10 days. 2 days is to
  # remove these observations with recent rainfall events. We use the 10-day
  # summary in the path analysis.
  
  if (sites[i] != "US-KFS"){
    
    precip_2day  <- zoo::rollsum(precip$precip_amf,2,fill = NA,align = "right")
    precip_10day <- zoo::rollsum(precip$precip_amf,10,fill = NA,align = "right")
    
  } else {
    
    precip_2day  <- zoo::rollsum(precip$precip_noaa,2,fill = NA,align = "right")
    precip_10day <- zoo::rollsum(precip$precip_noaa,10,fill = NA,align = "right")    
    
  }
  
  # Roll sum produces very small values when adding zeros (< 1e-16). This 
  # just converts them to 0. Kind of a nuisance, but whatever.
  precip_2day[precip_2day < 0.0001] <- 0
  precip_10day[precip_10day < 0.0001] <- 0
  
  precip_df_to_merge <- cbind(precip,data.frame(precip_2day = precip_2day,precip_10day = precip_10day))
  precip_df_to_merge <- precip_df_to_merge[,c("date","precip_2day","precip_10day")]
  precip_df_to_merge$date <- as.Date(precip_df_to_merge$date)
  
  fluxdat_daily <- data.table::merge.data.table(x = fluxdat_daily,
                                                y = precip_df_to_merge,
                                                by = "date",
                                                all.x = TRUE,all.y = FALSE)
  
  # ############################################################################
  # Add PhenoCam Gcc data 
  # ############################################################################
  
  phenodat$date <- as.Date(phenodat$date)
  
  fluxdat_daily <- data.table::merge.data.table(fluxdat_daily,
                                                phenodat,
                                                by="date",
                                                all.x=TRUE,all.y=FALSE)
  
  # ############################################################################
  # Calculate relevant ET statistics using the 'bigleaf' package 
  # ############################################################################
  
  # First calculate VPD using RH if it is not available from the downloaded 
  # fluxdata. 
  if (sum(is.na(fluxdat$VPD)) == nrow(fluxdat)){
    
    fluxdat_daily$VPD <- bigleaf::rH.to.VPD(fluxdat_daily$RH/100,fluxdat_daily$t_air)
  
  }

  # Calculate surface temperature using LW_IN and LW_OUT. Emissivity values are 
  # PFT specific and available in the metadata table.
  fluxdat_daily$t_surf <- bigleaf::radiometric.surface.temp(fluxdat_daily,
                                                            LW_up = "lw_out",
                                                            LW_down = "lw_in",
                                                            emissivity = emissivity[i])[,"Trad_degC"]
  
  # Obtain/back out aerodynamic conductance(Ga) from H 
  H <- fluxdat_daily$H
  rho <- bigleaf::air.density(fluxdat_daily$t_air,fluxdat_daily$pressure)
  cp <- bigleaf::bigleaf.constants()$cp
  delta_t <- fluxdat_daily$t_surf - fluxdat_daily$t_air
  fluxdat_daily$Ga <- H / (rho * cp * delta_t)
  
  # Calculate available energy (A) from H + LE
  fluxdat_daily$A <- fluxdat_daily$H + fluxdat_daily$LE
  
  # Calculate surface conductance by inverting Penman-Monteith
  fluxdat_daily$Gs <- bigleaf::surface.conductance(fluxdat_daily,
                                                   Tair = "t_air",
                                                   pressure = "pressure",
                                                   Rn = "A",
                                                   VPD = "VPD",
                                                   LE = "LE",
                                                   Ga = "Ga",
                                                   formulation = "Penman-Monteith")[,"Gs_ms"]
  
  # Calculate evaporative fraction
  fluxdat_daily$EF <- fluxdat_daily$LE / fluxdat_daily$A
  
  # Calculate ET using Penman-Monteith
  fluxdat_daily$ET <- bigleaf::LE.to.ET(fluxdat_daily$LE,Tair = fluxdat_daily$t_air)
  
  # ############################################################################
  # Identify days to remove from analysis 
  # ############################################################################
  # (1) Average midday H or LE < 0 
  # (2) Precip recorded in past 2 days
  # (3) Low VPD values causing very high Gs estimates (e.g., Novick et al. 2016)
  # (4) Negative average Tsurf - Tair difference
  # (5) Surface conductance less than 0
  # (6) Negative Ga or anomalously high Ga
  
  to_remove_bool_1 <- fluxdat_daily$H < 0 | fluxdat_daily$LE < 0
  to_remove_bool_2 <- fluxdat_daily$precip_2day > 0
  to_remove_bool_3 <- fluxdat_daily$t_surf - fluxdat_daily$t_air < 0
  to_remove_bool_4 <- fluxdat_daily$VPD < 0.6
  to_remove_bool_5 <- fluxdat_daily$Gs <= 0
  to_remove_bool_6 <- fluxdat_daily$Ga <= 0 | fluxdat_daily$Ga > 0.5
  
  to_remove_bool_all <- to_remove_bool_1 | to_remove_bool_2 | to_remove_bool_3 |
                        to_remove_bool_4 | to_remove_bool_5 | to_remove_bool_6
  to_remove_bool_all[is.na(to_remove_bool_all)] <- TRUE
  fluxdat_daily$to_remove_id <- as.numeric(to_remove_bool_all)

  # ############################################################################
  # Export data frame to csv file 
  # ############################################################################
  
  fluxdat_daily[is.na(fluxdat_daily)] <- -9999 # Replace NA with -9999
  
  # Write out daily flux data file
  setwd(paste0(wdir,"/results/flux_data/daily"))
  write.csv(fluxdat_daily,
            sprintf("%s_daily_values.csv",sites[i]),
            row.names = FALSE)
  
}