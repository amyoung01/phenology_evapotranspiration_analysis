rm(list=ls())
cat("\14")
graphics.off()

require(lavaan) # Latent Variable Analysis (version 0.6-6)
require(car) # Companion to Applied Regression (version 3.0-10)
require(lubridate) # Make Dealing with Dates a Little easier (version 1.7.10)

wdir <- "/Volumes/GoogleDrive/My Drive/W/projects/Young_evapotranspiration_phenology_analysis"

setwd(paste0(wdir,"/data/ancillary_data"))

phenoflux_metadata <- read.csv("pheno_flux_sites_to_use.csv")

sites <- phenoflux_metadata$fluxsite
vegtypes <- phenoflux_metadata$vegtype

##########################
mdl_basic =
'
VPD~t_air+precip_10day+gcc
gcc~gdd+cdd+precip_10day+VPD
Gs~gcc+VPD+precip_10day
EF~Gs+VPD
'

mdl_full =
'
VPD~a1*t_air+a2*precip_10day+a3*gcc
gcc~b1*gdd+b2*cdd+b3*precip_10day+b4*VPD
Gs~c1*gcc+c2*VPD+c3*precip_10day
EF~d1*Gs+d2*VPD

de_gcc := c1
de_vpd := c2
de_p10 := c3

ie_gcc := a3*c2
ie_vpd := b4*c1
ie_p10 := a2*c2 + b3*c1

te_gcc := de_gcc + ie_gcc
te_vpd := de_vpd + ie_vpd
te_p10 := de_p10 + ie_p10
'

resp_vars <- c()
pred_vars <- c()

mdl_parts <- strsplit(mdl_basic,"\n")[[1]]
mdl_parts <- mdl_parts[2:length(mdl_parts)]

for (m in 1:length(mdl_parts)){

  mdl_m = mdl_parts[m]
  mdl_m_split = strsplit(mdl_m,"\\~")[[1]]

  resp_vars <- c(resp_vars,mdl_m_split[1])

  pred_vars_m <- strsplit(mdl_m_split[2],"\\+")[[1]]

  for (p in 1:length(pred_vars_m)){

    pred_vars_mp <- pred_vars_m[p]

    if (any(pred_vars == pred_vars_mp)){

      next

    }

    pred_vars <- c(pred_vars,pred_vars_mp)

  }

}

col_names <- c("site","response","vegtype",pred_vars)

gof_stats <- matrix(NA,nrow = length(sites),2)
r2 <- matrix(NA,nrow = length(sites),length(resp_vars))

rownames(gof_stats) <- sites; colnames(gof_stats) <- c("CFI","SRMR")
colnames(r2) <- resp_vars

# Empty matrices to store parameter estimates and standard errors
# ------------------------------------------------------------------
params <- matrix(-9999,nrow=length(sites)*length(resp_vars),length(col_names))
colnames(params) <- col_names
params <- as.data.frame(params)
params[,1] <- as.vector(t(matrix(rep(sites,length(resp_vars)),ncol=length(resp_vars),byrow=F)))
params[,2] <- rep(resp_vars,length(sites))
params[,3] <- as.vector(t(matrix(rep(vegtypes,length(resp_vars)),ncol=length(resp_vars),byrow=F)))

params_se <- matrix(-9999,nrow=length(sites)*length(resp_vars),length(col_names))
colnames(params_se) <- col_names
params_se <- as.data.frame(params_se)
params_se[,1] <- as.vector(t(matrix(rep(sites,length(resp_vars)),ncol=length(resp_vars),byrow=F)))
params_se[,2] <- rep(resp_vars,length(sites))
params_se[,3] <- as.vector(t(matrix(rep(vegtypes,length(resp_vars)),ncol=length(resp_vars),byrow=F)))

Gs_eff <- matrix(-9999,nrow=length(sites),ncol=9)

nyr = matrix(NA,nrow=length(sites),ncol=1)

##########################################################
for (i in 1:length(sites)){

  setwd(paste0(wdir,"/test/test_linearity/results/flux_data/model_matrices"))
  data = read.csv(sprintf("%s_model_matrix.csv",sites[i]))
  data[data == -9999] = NA

  model_matrix = data[,unique(c(resp_vars,pred_vars))]

  sem_mdl = lavaan::sem(model = mdl_full,data=model_matrix)
  sem_results = summary(sem_mdl,fit=TRUE,rsquare=TRUE)

  mlr_models <- list()

  for (k in 1:length(resp_vars)){

    setwd(paste0(wdir,"/test/test_linearity/results/path_analysis_results/partial_regression_results"))

    lm_i = eval(parse(text = sprintf("lm(%s,data=as.data.frame(model_matrix))",mdl_parts[k])))
    mlr_models[[k]] <- lm_i
    avp <- avPlots(lm_i,main=sites[i])

    resp_var_k <- resp_vars[k]
    pred_vars_k <- strsplit(mdl_parts[k],paste0(resp_var_k,"~"))[[1]][2]
    pred_vars_k <- strsplit(pred_vars_k,"\\+")[[1]]

    mtx_row_id <- which(params$site == sites[i] & params[,2] == resp_var_k)

    for (v in 1:length(pred_vars_k)){

      id_v <- which(sem_results$PE$lhs == resp_var_k &
                      sem_results$PE$rhs == pred_vars_k[v])

      mtx_col_id <- which(colnames(params) == pred_vars_k[v])

      params[mtx_row_id,mtx_col_id] <- round(sem_results$PE$est[id_v],3)
      params_se[mtx_row_id,mtx_col_id] <- round(sem_results$PE$se[id_v],3)

      pp_plot_data <- data.frame(avp[[v]])
      pp_plot_data[,paste0(resp_var_k,"_pred")] <- as.vector(predict(lm(pp_plot_data[,2] ~ pp_plot_data[,1])))

      write.csv(pp_plot_data,file = sprintf("%s_partReg_%s-%s.csv",sites[i],resp_var_k,pred_vars_k[v]),
                row.names = FALSE)

    }

    r2_id <- which(sem_results$PE$op == "r2" &
                     sem_results$PE$lhs == resp_var_k)
    r2[i,k] <- round(sem_results$PE$est[r2_id],2)

  }

  Gs_eff[i,1] <- sem_results$PE$est[sem_results$PE$label == "de_p10"]
  Gs_eff[i,2] <- sem_results$PE$est[sem_results$PE$label == "ie_p10"]
  Gs_eff[i,3] <- sem_results$PE$est[sem_results$PE$label == "te_p10"]

  Gs_eff[i,4] <- sem_results$PE$est[sem_results$PE$label == "de_vpd"]
  Gs_eff[i,5] <- sem_results$PE$est[sem_results$PE$label == "ie_vpd"]
  Gs_eff[i,6] <- sem_results$PE$est[sem_results$PE$label == "te_vpd"]

  Gs_eff[i,7] <- sem_results$PE$est[sem_results$PE$label == "de_gcc"]
  Gs_eff[i,8] <- sem_results$PE$est[sem_results$PE$label == "ie_gcc"]
  Gs_eff[i,9] <- sem_results$PE$est[sem_results$PE$label == "te_gcc"]


  gof_stats[i,1] <- sem_results$FIT[9]
  gof_stats[i,2] <- sem_results$FIT[21]

  nyr[i] = nrow(data)/365

}

gof_stats <- data.frame(gof_stats)
r2 <- data.frame(r2)
params <- data.frame(params)
params_se <- data.frame(params_se)

rownames(Gs_eff) <- sites
colnames(Gs_eff) <- c("de_precip_10day","ie_precip_10day","te_precip_10day",
                      "de_vpd","ie_vpd","te_vpd",
                      "de_gcc","ie_gcc","te_gcc")

rownames(r2) <- sites

rownames(nyr) <- sites
colnames(nyr) <- "site_yr"

setwd(paste0(wdir,"/test/test_linearity/results/path_analysis_results"))
write.csv(params,"sem_param_est.csv",row.names = FALSE)
write.csv(params_se,"sem_param_se.csv",row.names = FALSE)

write.csv(Gs_eff,"total_effect_on_Gs.csv")
write.csv(gof_stats,"sem_gof_stats.csv")
write.csv(r2,"rsquared_stats.csv")
write.csv(nyr,"n_sites_years.csv")
