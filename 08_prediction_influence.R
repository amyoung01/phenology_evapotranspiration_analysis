rm(list=ls())
cat("\14")
graphics.off()

require(lavaan) # Latent Variable Analysis (version 0.6-6)
require(car) # Companion to Applied Regression (version 3.0-10)
require(lubridate)

wdir <- "/Volumes/GoogleDrive/My Drive/Young_evapotranspiration_phenology_analysis"

setwd(paste0(wdir,"/data/ancillary_data"))

phenoflux_metadata <- read.csv("pheno_flux_sites_to_use.csv")

sites <- phenoflux_metadata$fluxsite
vegtypes <- phenoflux_metadata$vegtype
vegtypes[sites == "US-Ro4"] <- "GR"
vegtypes[sites == "US-Mpj" | sites == "US-Ton"] <- "SA"

transformations = data.frame(var=c("VPD","precip_10day","Gs"),
                             transformation=c("log","sqrt","log"))

resp_vars <- c("Gs","EF")
pred_vars <- c("VPD","precip_10day","gcc")

r2 <- matrix(NA,nrow = length(sites),ncol = length(pred_vars) + 1)

################################################################################
for (i in 1:length(sites)){
  
  setwd(paste0(wdir,"/results/flux_data/model_matrices"))
  data = read.csv(sprintf("%s_model_matrix.csv",sites[i]))
  data[data == -9999] = NA
  
  model_matrix = data[,unique(c(resp_vars,pred_vars))]
  
  for (e in 1:(length(pred_vars) + 1)){
    
    if (e == 1){
      
      lm_e <- lm(Gs ~ precip_10day + VPD + gcc,data = model_matrix)
      r2[i,e] <- summary(lm_e)$adj.r.squared
      
    } else {
      
      model_matrix_e <- model_matrix
      rsamp <- sample(nrow(model_matrix_e),size = nrow(model_matrix_e),replace = FALSE)
      model_matrix_e[,pred_vars[e-1]] <- model_matrix_e[,pred_vars[e-1]][rsamp]
      
      lm_e <- lm(Gs ~ precip_10day + VPD + gcc,data = model_matrix_e)
      r2[i,e] <- summary(lm_e)$adj.r.squared
      
    }
    
  }
  
}

r2 <- round(r2,2)

relChange_VPD = r2[,1]/r2[,2]
relChange_precip_10day = r2[,1]/r2[,3]
relChange_gcc = r2[,1]/r2[,4]

AI = phenoflux_metadata$AI
PFT = as.factor(phenoflux_metadata$vegtype)

# plot(log(AI),relChange_gcc)

# setwd(paste0(wdir,"/results/path_analysis_results"))
r2 <- as.data.frame(r2)
rownames(r2) <- sites
colnames(r2) <- c("all",pred_vars)
r2$PFT <- PFT
r2$AI <- AI

r2 <- r2[,c(5,6,1:4)]

write.csv(r2,paste0(wdir,"/results/prediction_influence/rsquared_values.csv"))
