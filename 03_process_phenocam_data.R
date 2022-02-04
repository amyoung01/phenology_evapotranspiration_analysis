rm(list=ls())
graphics.off()

library(phenocamr)
library(anytime)

wdir <- "/Volumes/GoogleDrive/My Drive/Young_evapotranspiration_phenology_analysis"

setwd(paste0(wdir,"/data/ancillary_data"))
phenocam_flux_metadata_table <- read.csv("pheno_flux_sites_to_use.csv")

phenos <- phenocam_flux_metadata_table$phenosite
vegtype <- phenocam_flux_metadata_table$vegtype

rmse_stats <- c("mean","50","75","90")

for (i in 1:length(phenos)){
    
    setwd(sprintf('%s//data//raw_data//phenocam//%s',wdir,as.character(phenos[i])))
    pheno_dir <- getwd()
    versions <- files <- list.files()
    
    for (v in 1:length(versions)){
        
        setwd(sprintf('%s//%s',pheno_dir,versions[v]));
        file_to_read = paste0(sprintf("%s_%s_1day_transition_dates.csv",
                                      phenos[i],versions[v]))
        
        fid <- file(description = file_to_read,open = "rt")
        header_lines <- readLines(con = fid,n = 15)
        close(fid)
        
        header_lines <- header_lines[12:15]
        rmse <- NA * numeric(length = length(header_lines))
        
        for (r in 1:4) {
          
          str_r <- unlist(strsplit(header_lines[r]," "))
          rmse[r] <- as.numeric(str_r[5])
          
        }
        
        min_rmse <- which(rmse == min(rmse))
        min_rmse <- ifelse(length(min_rmse) > 1,min_rmse[1],min_rmse)
        
        fn <- sprintf("%s_%s_1day.csv",phenos[i],versions[v])
        
        # pheno_td <- smooth_ts(fn,out_dir = getwd())
        pheno_td <- read_phenocam(fn)
        
        new_dates_fwd <- transition_dates(pheno_td,
                                      lower_thresh = 0.10,
                                      middle_thresh = 0.50,
                                      upper_thresh = 0.90,
                                      percentile = rmse_stats[min_rmse],
                                      reverse = FALSE,
                                      plot = FALSE)
        
        new_dates_rev <- transition_dates(pheno_td,
                                          lower_thresh = 0.10,
                                          middle_thresh = 0.50,
                                          upper_thresh = 0.90,
                                          percentile = rmse_stats[min_rmse],
                                          reverse = TRUE,
                                          plot = FALSE)
        
        new_dates_fwd$direction <- "rising"
        new_dates_rev$direction <- "falling"
        
        new_dates <- rbind(new_dates_fwd,new_dates_rev)
        
        new_dates <- new_dates[,c(ncol(new_dates),1:3)]
        
        new_dates[,2] <- anydate(new_dates[,2])
        new_dates[,3] <- anydate(new_dates[,3])
        new_dates[,4] <- anydate(new_dates[,4])

        if (v == 1){
            
            gcc_trans_dates <- new_dates

        } else {

            gcc_trans_dates <- rbind(gcc_trans_dates, new_dates)
            
        }
        
        pheno_ts_v <- pheno_td$data #read.csv(pheno_ts_name,header = TRUE,sep = ",",skip = 24)
        
        col_id_1 <- which(colnames(pheno_ts_v) == paste0("gcc_",rmse_stats[min_rmse]))
        col_id_2 <- which(colnames(pheno_ts_v) == paste0("smooth_gcc_",rmse_stats[min_rmse]))
        
        new_ts <- data.frame(date = pheno_ts_v$date,
                             gcc = pheno_ts_v[,col_id_1],
                             smooth_gcc = pheno_ts_v[,col_id_2])
        
        if (v == 1){
            
            pheno_ts <- new_ts
            
        } else {
            
            if (as.Date(new_ts$date[1]) < as.Date(pheno_ts$date[nrow(pheno_ts)])){
                
                id <- which(as.Date(new_ts$date) == as.Date(pheno_ts$date[nrow(pheno_ts)]))
                new_ts <- new_ts[(id+1):nrow(new_ts),]
                
            }
            
            pheno_ts <- rbind(pheno_ts, new_ts)
            
        }
        
    }
    
    pheno_ts[is.na(pheno_ts)] <- -9999
    
    pheno_ts_file_to_write <- sprintf("%s_gcc_time_series.csv",phenos[i])
    write.csv(pheno_ts,paste0(wdir,
                              "/results/processed_phenocam_data/time_series/",
                              pheno_ts_file_to_write),
              row.names = FALSE)
    
    pheno_td_file_to_write <- sprintf("%s_gcc_transition_dates.csv",phenos[i])
    
    sort_id <- sort.int(gcc_trans_dates$transition_50,index.return = TRUE)$ix
    gcc_trans_dates <- gcc_trans_dates[sort_id,]
    
    write.csv(gcc_trans_dates,
              paste0(wdir,
                     "/results/processed_phenocam_data/transition_dates/",
                     pheno_td_file_to_write),
              row.names = FALSE) 
    
}