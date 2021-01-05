# ==== Bayesian meta-analysis workflow ====  
# Appendix for the manuscript; model is described in the paper

rm(list = ls())

# !!! IMPORTANT !!! JAGS needs to be installed prior to running this script: http://mcmc-jags.sourceforge.net/
# Load packages; install R packages as needed:
meta_packs <- c('rjags', 'coda', 'mcmcplots', 'R2jags')
new.packages <- meta_packs[!(meta_packs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(meta_packs, require, character.only = T)

# Model inputs from the Dieleman et al. 2010 dataset; read-in data:
record_input <- read.csv('./data/record_level_input.csv')
attach(record_input)

study_input <- read.csv('./data/study_level_input.csv')
attach(study_input)

Nv <- 3 # number of variable is 3: AA, AB, and SCE 
Nr <- dim(record_input)[1] # number of records
Nk <- max(record_input$Study) # number of studies

# Caculate average MAT and MAT across all studies:
MATave <- mean(MAT, na.rm = T)
MAPave <- mean(MAP, na.rm = T)

# Calculate average study-level duration:
dur_study <- aggregate(Duration ~ Study, data=cbind(Study, Duration), FUN=mean)
aveDUR <- mean(dur_study$Duration)
# Study-level average duration:
Duration_Study <- rep(NA, Nk)
Duration_Study[dur_study$Study] <- dur_study$Duration
Duration_Study[-dur_study$Study] <- aveDUR

# MAT and MAP values at which log-response-ratios will be predicted:
MATvals <- seq(-10, 30, 1)
MAPvals <- seq(-0, 2000, 100)

# Input matrix for specifying Wishart priors for precision matrices:
RO <- diag(c(0.1, 0.1, 0.1)) 

# Collect the list of model inputs:
meta_input_data <- list('Nv' = Nv, 'Nr' = Nr, 'Nk' = Nk, 'RO' = RO, 
                        'R' = R, 'TrCode' = TrCode, 'varR' = varR, 
                        'Study' = Study, 'VarID' = varID, 
                        'MAT' = MAT, 'MATave' = MATave, 'MATvals' = MATvals, 'Nt' = length(MATvals),
                        'MAP' = MAP, 'MAPave' = MAPave, 'MAPvals' = MAPvals, 'Np' = length(MAPvals),
                        'Duration' = Duration, 'Duration_Study' = Duration_Study, 'aveDUR' = aveDUR)

# Source the Bayesian model; the model is described in detail in file Bayes_meta_model_multivar.r
source('./Bayes_meta_model_multivar.r')

# Variables to monitor: 
# Line 1: Site-level treatment effects; 
# Line 2: "global" treatment effects and thier Bayesian p-values; 
# Line 3: coefficients of the covariates and their Bayesian p-values; 
# Line 4: predicted MATXMAP interaction effects and site-level effects; 
# Line 5: other quantities to monitor.
meta_monitor_vars <- c('tr_co2', 'tr_warm', 'tr_cw', 'tr_add',
                       'tr_co2_star', 'tr_warm_star', 'tr_cw_star','tr_add_star', 'Itr_co2_star', 'Itr_warm_star', 'Itr_cw_star','Itr_add_star',
                       'beta_dur', 'beta_MAT', 'beta_MAP', 'beta_MAPMAT', 'Ibeta_dur', 'Ibeta_MAT', 'Ibeta_MAP', 'Ibeta_MAPMAT',
                       'pred_LRR', 'pred_LRR_dur0', 'LRR_co2', 'LRR_warm', 'LRR_cw', 'LRR_add',
                       'R_rep', 'mu_t', 'Sigma', 'rho') 

# Run the Bayesian meta-analysis model and save the results:
load.module("dic")
# Initialize model and run past convergence.
meta_jags <- jags(model.file = meta_model_multivar, data = meta_input_data, parameters.to.save = meta_monitor_vars,
                  n.chains=3, n.iter = 100000, n.thin = 20) 
# Update jags model to obtain posterior samples:
meta_out_jags <- update(meta_jags, 100000, n.thin = 100) # thin every 500 to obtain independent samples
meta_out <- as.mcmc(meta_out_jags)
dir.create(file.path('./output/'), showWarnings = FALSE)
save(meta_out, file = paste('output/Meta_analysis_output_', Sys.Date(), '.rda', sep = ''))

# Create variables to collect posterior samples:
pred_LRR_monitor <- paste('pred_LRR[', outer(1:21, outer(1:36, outer(1:3, 1:3, paste, sep = ','), 
                                                            paste, sep = ','), paste, sep = ','),']', sep = '')
tr_co2_monitor <- paste('tr_co2[', outer(1:Nk, 1:Nv, paste, sep = ','), ']', sep = '')
tr_warm_monitor <- paste('tr_warm[', outer(1:Nk, 1:Nv, paste, sep = ','), ']', sep = '')
tr_cw_monitor <- paste('tr_cw[', outer(1:Nk, 1:Nv, paste, sep = ','), ']', sep = '')
tr_co2_star_monitor <- paste('tr_co2_star[', 1:Nv, ']', sep = '')
tr_warm_star_monitor <- paste('tr_warm_star[', 1:Nv, ']', sep = '')
tr_cw_star_monitor <- paste('tr_cw_star[', 1:Nv, ']', sep = '')
beta_MAP_monitor <- paste('beta_MAP[', outer(1:Nv, 1:3, paste, sep = ','),']', sep = '')
beta_MAT_monitor <- paste('beta_MAT[', outer(1:Nv, 1:3, paste, sep = ','),']', sep = '')
beta_MAPMAT_monitor <- paste('beta_MAPMAT[', outer(1:Nv, 1:3, paste, sep = ','),']', sep = '')
beta_dur_monitor <- paste('beta_dur[', outer(1:Nv, 1:3, paste, sep = ','),']', sep = '')
rho_monitor <- paste('rho[', outer(1:3, outer(1:3, 1:3, paste, sep = ','), paste, sep= ','), ']', sep = '')

# Produce MCMC plots to check for convergence:
mcmcplot(meta_out[, c('deviance', pred_LRR_monitor, rho_monitor, tr_co2_monitor, tr_warm_monitor, tr_cw_monitor, 
                      tr_co2_star_monitor, tr_warm_star_monitor, tr_cw_star_monitor, 
                      beta_MAT_monitor, beta_MAP_monitor, beta_MAPMAT_monitor, beta_dur_monitor)]) #

# Now output can be used to compute posterior statistics and plots (code not provided).

detach(record_input)
detach(study_input)

gc()
