# Multivariate Bayesian model for AB, BB, and SCE data. 
# The "general" version of the model that works with more treatment types starts at Line #222
meta_model_multivar <- function(){
  for(r in 1:Nr){
    # Record-level data are organized with a single data variable (vector) for response ratio (R)
    # and the associated pooled variance for R (varR), along with additional data vectors
    # for variable ID (VarID: 1 = aboveground biomass (AB), 2 = belowground biomass (BB), and 3 = soil CO2
    # exchange (SCE)), study ID (Study), and treatment code (TrCode: 1 = CO2-only, 2 = Warming-only, 
    # 3 = both CO2 and warming). Covariates MAP and MAT are assumed to vary at the study level; whereas
    # Duration varies at the record level.
    
    # Likelihood of log-response ratio "data":
    R[r] ~ dnorm(mu_t[r], tau_y[r])
    tau_y[r] <- 1/varR[r]
    
    # Replicated log-response ratio "data" (e.g., for evaluating model fit)
    R_rep[r] ~ dnorm(mu_t[r], tau_y[r])
    
    # Mean model whereby the study-level main effect of CO2-only (TrCode = 1), warming-only (TrCode = 2), 
    # and both (TrCode = 3) are explicilty accounted for within the same model (i.e., tr_x terms); 
    # include effects (beta terms) of continuous covariates (Duration, MAT, MAP) by modeling the beta's as 
    # matrices indexed by VarID and TrCode. Center covariates about the dataset mean.
    mu_t[r] <- tr_co2[Study[r],VarID[r]]*(TrCode[r] == 1) + tr_warm[Study[r],VarID[r]]*(TrCode[r] == 2) + 
      tr_cw[Study[r],VarID[r]] * (TrCode[r] == 3) + beta_dur[VarID[r],TrCode[r]]*Duration[r] + 
      beta_MAP[VarID[r],TrCode[r]]*(MAP[Study[r]] - MAPave) + beta_MAT[VarID[r],TrCode[r]]*(MAT[Study[r]] - MATave) +
      beta_MAPMAT[VarID[r],TrCode[r]]*(MAP[Study[r]] - MAPave)*(MAT[Study[r]] - MATave)
  }
  
  # Calculate the predicted, study-level log-response ratio (LRR) for each Study/Variable combination.
  # Duration_Study is read in as data and represents the average duration for study s. (This code is not 
  # necessary for implementing the meta-analysis.) The workflow file creates the Duration_Study data.
  # Study loop (s)
  for(s in 1:Nk){ 
    # Variable loop (v). 1 = AB, 2 = BB, and 3 = SCE
    for(v in 1:Nv){
      LRR_co2[s, v] <- tr_co2[s,v] +  beta_dur[v,1]*Duration_Study[s] + 
        beta_MAP[v,1]*(MAP[s] - MAPave) + beta_MAT[v,1]*(MAT[s] - MATave) +
        beta_MAPMAT[v,1]*(MAP[s] - MAPave)*(MAT[s] - MATave)
      
      LRR_warm[s, v] <- tr_warm[s,v] +  beta_dur[v,2]*Duration_Study[s] + 
        beta_MAP[v,2]*(MAP[s] - MAPave) + beta_MAT[v,2]*(MAT[s] - MATave) +
        beta_MAPMAT[v,2]*(MAP[s] - MAPave)*(MAT[s] - MATave)
      
      LRR_cw[s, v] <- tr_cw[s,v] +  beta_dur[v,3]*Duration_Study[s] + 
        beta_MAP[v,3]*(MAP[s] - MAPave) + beta_MAT[v,3]*(MAT[s] - MATave) +
        beta_MAPMAT[v,3]*(MAP[s] - MAPave)*(MAT[s] - MATave)
      
      LRR_add[s, v] <- tr_co2[s,v] + tr_warm[s,v] +  beta_dur[v,3]*Duration_Study[s] + 
        beta_MAP[v,3]*(MAP[s] - MAPave) + beta_MAT[v,3]*(MAT[s] - MATave) +
        beta_MAPMAT[v,3]*(MAP[s] - MAPave)*(MAT[s] - MATave)
    }
  }
  
  # Predict the LRR of AB, BB, and SCE under different MAP and MAT levels, for each treatment type.
  # Code used to generate contour plots of predicted LLR of AB, BB, and SCE across the range of MAT and MAP values; 
  # Code not necessary for implementing the meta-analysis. The workflow file creats the MAPvals and MATvals data.
  # MAT levels: -5, -4, -3, ..., 29, and 30 degree-celcius (36 values); 
  for(t in 1:Nt){
    # MAP: 0, 100, 200, 300, ..., 2000 mm/yr (21 values)
    for(p in 1:Np){
      for(v in 1:Nv){
        # Predicted log-response ratios at the average duration reported across al studies:
        # CO2-only treatment
          pred_LRR[p,t, v, 1] <- tr_co2_star[v] + beta_dur[v,1]*aveDUR +
            beta_MAP[v,1]*(MAPvals[p] - MAPave) + beta_MAT[v,1]*(MATvals[t] - MATave) +
            beta_MAPMAT[v,1]*(MAPvals[p] - MAPave)*(MATvals[t] - MATave)
        
        # Warming-only treatment  
          pred_LRR[p,t, v, 2] <- tr_warm_star[v] + beta_dur[v,2]*aveDUR +
            beta_MAP[v,2]*(MAPvals[p] - MAPave) + beta_MAT[v,2]*(MATvals[t] - MATave) +
            beta_MAPMAT[v,2]*(MAPvals[p] - MAPave)*(MATvals[t] - MATave)
        
        # CO2 and warming treatment
          pred_LRR[p,t, v, 3] <- tr_cw_star[v] + beta_dur[v,3]*aveDUR +
            beta_MAP[v,3]*(MAPvals[p] - MAPave) + beta_MAT[v,3]*(MATvals[t] - MATave) +
            beta_MAPMAT[v,3]*(MAPvals[p] - MAPave)*(MATvals[t] - MATave)
        
        # Predicted log-response ratios at the beginning of a study (Duration = 0 years)
        # CO2-only treatment
          pred_LRR_dur0[p,t, v, 1] <- tr_co2_star[v] + beta_dur[v,1]*0 +
            beta_MAP[v,1]*(MAPvals[p] - MAPave) + beta_MAT[v,1]*(MATvals[t] - MATave) +
            beta_MAPMAT[v,1]*(MAPvals[p] - MAPave)*(MATvals[t] - MATave)
        
        # Warming-only treatment
          pred_LRR_dur0[p,t, v, 2] <- tr_warm_star[v] + beta_dur[v,2]*0 +
            beta_MAP[v,2]*(MAPvals[p] - MAPave) + beta_MAT[v,2]*(MATvals[t] - MATave) +
            beta_MAPMAT[v,2]*(MAPvals[p] - MAPave)*(MATvals[t] - MATave)
        
        # CO2 and warming treatment
          pred_LRR_dur0[p,t, v, 3] <- tr_cw_star[v] + beta_dur[v,3]*0 +
            beta_MAP[v,3]*(MAPvals[p] - MAPave) + beta_MAT[v,3]*(MATvals[t] - MATave) +
            beta_MAPMAT[v,3]*(MAPvals[p] - MAPave)*(MATvals[t] - MATave)
      }
    }
  }
  
  # Mulitvariate, hierarhical priors for study-level main effects:
  # Study loop
  for(k in 1:Nk){
    tr_co2[k,1:Nv] ~ dmnorm(tr_co2_star[1:Nv], Omega_co2[1:Nv,1:Nv])
    tr_warm[k,1:Nv] ~ dmnorm(tr_warm_star[1:Nv], Omega_warm[1:Nv,1:Nv])
    tr_cw[k,1:Nv] ~ dmnorm(tr_cw_star[1:Nv], Omega_cw[1:Nv,1:Nv])
    # Variable loop
    for(v in 1:Nv){
      # Compute predicted additive effects for each response variable:
      tr_add[k,v] <- tr_co2[k,v] + tr_warm[k,v]
    }
  }
  
  # Hyperpriors for global treatment-effect parameters. Loop through first
  # 2 response variables (AB and BB), then define priors for last response 
  # variable (SCE) to account for sparse data associated with SCE, recognizing
  # that some parameters cannot be informed by the SCE data. Otherwise, could
  # have looped through all 3 responses variable.
  for(v in 1:(Nv-1)){
    # Priors for means:
    tr_co2_star[v] ~ dnorm(0, 0.0001)
    tr_warm_star[v] ~ dnorm(0, 0.0001)
    tr_cw_star[v] ~ dnorm(0, 0.0001)
    # Global-level additive effect of CO2 and warming
    tr_add_star[v] <-  tr_co2_star[v] + tr_warm_star[v]
    
    # Priors for covariate effects:
    for(t in 1:3){
      beta_dur[v,t] ~ dnorm(0,0.0001)
      beta_MAT[v,t] ~ dnorm(0,0.0001)
      beta_MAP[v,t] ~ dnorm(0,0.0001)
      beta_MAPMAT[v,t] ~ dnorm(0,0.0001)
    }
  }
  for(v in 3){
    # Priors for means:
    tr_co2_star[v] ~ dnorm(0, 0.0001)
    tr_warm_star[v] ~ dnorm(0, 0.0001)
    tr_cw_star[v] ~ dnorm(0, 0.0001)
    # Global-level additive effect of CO2 and warming
    tr_add_star[v] <-  tr_co2_star[v] + tr_warm_star[v]
    
    # Priors for covariate effects:
    for(t in 1:2){
      beta_dur[v,t] ~ dnorm(0,0.0001)
      beta_MAT[v,t] ~ dnorm(0,0.0001)
      beta_MAP[v,t] ~ dnorm(0,0.0001)
      beta_MAPMAT[v,t] ~ dnorm(0,0.0001)
    }
    # Since there is limited climate data and/or little variation in the
    # climate data for SCE in the CO2+warming studies, cannot estimate
    # effects of MAT and MAP, so set coefficients = 0.
    beta_dur[v,3] ~ dnorm(0,0.0001)
    beta_MAT[v,3] <- 0
    beta_MAP[v,3] <- 0
    beta_MAPMAT[v,3] <- 0
  }
  
  # For computing Bayesian p-values (which are the
  # posterior means of these indicator quantities)
  for(v in 1:Nv){
    for(t in 1:3){
      # Ibeta = 1 if beta > 0, and Ibeta = 0 otherwise.
      # Here, ingore Ibeta and p-value for beta_X[v,3] where
      # X = MAT, MAP or MAPMAT
      Ibeta_dur[v,t] <- (beta_dur[v,t] > 0) * 1
      Ibeta_MAT[v,t] <- (beta_MAT[v,t] > 0) * 1
      Ibeta_MAP[v,t] <- (beta_MAP[v,t] > 0) * 1
      Ibeta_MAPMAT[v,t] <- (beta_MAPMAT[v,t] > 0) * 1
    }
    # Similarly for the global log-response ratios:
    Itr_co2_star[v] <- (tr_co2_star[v] > 0) * 1
    Itr_warm_star[v] <- (tr_warm_star[v] > 0) * 1
    Itr_cw_star[v] <- (tr_cw_star[v] > 0) * 1
    Itr_add_star[v] <- (tr_add_star[v] > 0) * 1
  }
  
  # Compute covariance matrix:
  Sigma[1, 1:Nv,1:Nv] <- inverse(Omega_co2[1:Nv,1:Nv])
  Sigma[2, 1:Nv,1:Nv] <- inverse(Omega_warm[1:Nv,1:Nv])
  Sigma[3, 1:Nv,1:Nv] <- inverse(Omega_cw[1:Nv,1:Nv])
  
  # Compute correlation matrix:
  for(t in 1:3){
    for(r in 1:(Nv-1)){
      # Among study standard deviations for each response variable:
      rho[t,r,r] <- sqrt(Sigma[t,r,r])
      for(c in (r+1):Nv){
        # Pairwise correlations amoung pairs of response variables:
        rho[t,r,c] <- Sigma[t,r,c]/sqrt(Sigma[t,r,r] * Sigma[t,c,c])
        rho[t,c,r] <- rho[t,r,c]
      }
    }
    # Among study std dev for SCE response variable:
    rho[t,Nv,Nv] <- sqrt(Sigma[t,Nv,Nv])
  }
  
  # Fairly non-informative priors for precision matrix:
  # RO <- diag(c(0.1,0.1,0.1)) and is created and read in as data in the workflow file.
  Omega_co2[1:Nv,1:Nv] ~ dwish(RO[1:Nv, 1:Nv], Nv)
  Omega_warm[1:Nv,1:Nv] ~ dwish(RO[1:Nv, 1:Nv], Nv)
  Omega_cw[1:Nv,1:Nv] ~ dwish(RO[1:Nv, 1:Nv], Nv)
  
  # Missing data models for Duration, MAP, and MAT; can be deleted if there is no missing data:
  for(r in 1:Nr){
    # Likelihood of record-level duration data
    Duration[r] ~ dnorm(mu.dur, tau.dur)
  }
  for(k in 1:Nk){
    # Likelihood of study-level MAT and MAP data
    MAT[k] ~ dnorm(mu.mat, tau.mat)
    MAP[k] ~ dnorm(mu.map, tau.map)
  }
  # priors for means and precisions (or std deviations) in above likelihoods
  mu.dur ~ dnorm(0,0.0001)
  mu.map ~ dnorm(0,0.0001)
  mu.mat ~ dnorm(0,0.0001)
  tau.dur <- pow(sig.dur, -2)
  sig.dur ~ dunif(0,1000)
  tau.map <- pow(sig.map, -2)
  sig.map ~ dunif(0,1000)
  tau.mat <- pow(sig.mat, -2)
  sig.mat ~ dunif(0,1000)
}

# The "general" version of the model that works with more treatment types 
general_meta_model_multivar <- function(){
  for(r in 1:Nr){
    # Record-level data are organized with a single data variable (vector) for response ratio (R)
    # and the associated pooled variance for R (varR), along with additional data vectors
    # for variable ID (VarID: 1 = aboveground biomass (AB), 2 = belowground biomass (BB), and 3 = soil CO2
    # exchange (SCE)), study ID (Study), and treatment code (TrCode: 1 = CO2-only, 2 = Warming-only, 
    # 3 = both CO2 and warming). Covariates MAP and MAT are assumed to vary at the study level; whereas
    # Duration varies at the record level.
    
    # Likelihood of log-response ratio "data":
    R[r] ~ dnorm(mu_t[r], tau_y[r])
    tau_y[r] <- 1/varR[r]
    
    # YL modified: now varR is calculated based on Se2C, Se2T, NSeC, NSeT, YC, YT.
    varR[r] <- SE2[r]*(1/(NSeC[r]*YC[r]*YC[r])+1/(NSeT[r]*YT[r]*YT[r]))
    SE2[r] <- 1/(NSeC[r]+NSeT[r]-2)*((NSeC[r]-1)*NSeC[r]*Se2C[r] + (NSeT[r]-1)*NSeT[r]*Se2T[r])
    
    
    NSeC[r] <- NSeC.minus2[r]
    NSeT[r] <- NSeT.minus2[r]
    NSeC.minus2[r] ~ dpois(lamNC[VarID[r]])
    NSeT.minus2[r] ~ dpois(lamNT[VarID[r]])
    
    Se2C[r] ~ dgamma(alphaC[r], betaC[r])
    Se2T[r] ~ dgamma(alphaT[r], betaT[r])
    alphaC[r] <- 0.5*(NSeC[r] - 1)
    alphaT[r] <- 0.5*(NSeT[r] - 1)
    betaC[r] <- 0.5*(NSeC[r] - 1)*NSeC[r]*pow(sigXC[VarID[r]],-2)
    betaT[r] <- 0.5*(NSeT[r] - 1)*NSeT[r]*pow(sigXT[VarID[r]],-2)
    
    # Replicated log-response ratio "data" (e.g., for evaluating model fit)
    R_rep[r] ~ dnorm(mu_t[r], tau_y[r])
    
    # Mean model whereby the study-level main effect of CO2-only (TrCode = 1), warming-only (TrCode = 2), 
    # and both (TrCode = 3) are explicilty accounted for within the same model (i.e., tr_x terms); 
    # include effects (beta terms) of continuous covariates (Duration, MAT, MAP) by modeling the beta's as 
    # matrices indexed by VarID and TrCode. Center covariates about the dataset mean.
    mu_t[r] <- treat[Study[r],TrCode[r],VarID[r]] + # YL modified
      beta_dur[TrCode[r],VarID[r]]*Duration[r] + # YL modified also swaped the dimensions of the beta parameters for consistency
      beta_MAP[TrCode[r],VarID[r]]*(MAP[Study[r]] - MAPave) + 
      beta_MAT[TrCode[r],VarID[r]]*(MAT[Study[r]] - MATave) +
      beta_MAPMAT[TrCode[r],VarID[r]]*(MAP[Study[r]] - MAPave)*(MAT[Study[r]] - MATave)
  }
  
  # Calculate the predicted, study-level log-response ratio (LRR) for each Study/Variable combination.
  # Duration_Study is read in as data and represents the average duration for study s. (This code is not 
  # necessary for implementing the meta-analysis.) The workflow file creates the Duration_Study data.
  # Study loop (s)
  for(k in 1:Nk){ 
    # Variable loop (v). 1 = AB, 2 = BB, and 3 = SCE # YL modified: Eline, you can note your definition of varID here
    for(t in 1:Nt){
      for(v in 1:Nv){
        # YL modified: now calculating LLR for all k
        LRR[k,t,v] <- treat[k,t,v] +  
          beta_dur[t,v]*Duration_Study[k] + 
          beta_MAP[t,v]*(MAP[k] - MAPave) + 
          beta_MAT[t,v]*(MAT[k] - MATave) +
          beta_MAPMAT[t,v]*(MAP[k] - MAPave)*(MAT[k] - MATave)        
      }
    }
  }
  
  # Predict the LRR of AB, BB, and SCE under different MAP and MAT levels, for each treatment type.
  # Code used to generate contour plots of predicted LLR of AB, BB, and SCE across the range of MAT and MAP values; 
  # Code not necessary for implementing the meta-analysis. The workflow file creats the MAPvals and MATvals data.
  # MAT levels: -5, -4, -3, ..., 29, and 30 degree-celcius (36 values); 
  
  # YL modified: REMOVED the pred_LRR section for now. Can add it back in later when we identify some treatment/variable effects of interest.
  
  # Mulitvariate, hierarhical priors for study-level main effects:
  # Study loop
  for(k in 1:Nk){
    for(t in 1:Nt){
      treat[k,t,1:Nv] ~ dmnorm(treat_star[t,1:Nv], Omega[t,1:Nv,1:Nv]) # YL modified
    }
  }
  
  # Hyperpriors for global treatment-effect parameters. Loop through first
  # 2 response variables (AB and BB), then define priors for last response 
  # variable (SCE) to account for sparse data associated with SCE, recognizing
  # that some parameters cannot be informed by the SCE data. Otherwise, could
  # have looped through all 3 responses variable.
  for(v in 1:Nv){
    for(t in 1: Nt){
      # YL modified
      # Priors for means:
      treat_star[t,v] ~ dnorm(0, 0.0001)
      # Priors for covariate effects:
      beta_dur[t,v] ~ dnorm(0,0.0001)
      beta_MAT[t,v] ~ dnorm(0,0.0001)
      beta_MAP[t,v] ~ dnorm(0,0.0001)
      beta_MAPMAT[t,v] ~ dnorm(0,0.0001)
    }
  }
  
  # For computing Bayesian p-values (which are the
  # posterior means of these indicator quantities)
  for(v in 1:Nv){
    for(t in 1:Nt){
      # Ibeta = 1 if beta > 0, and Ibeta = 0 otherwise.
      # Here, ingore Ibeta and p-value for beta_X[v,3] where
      # X = MAT, MAP or MAPMAT
      # YL modified:
      Ibeta_dur[t,v] <- (beta_dur[t,v] > 0) * 1
      Ibeta_MAT[t,v] <- (beta_MAT[t,v] > 0) * 1
      Ibeta_MAP[t,v] <- (beta_MAP[t,v] > 0) * 1
      Ibeta_MAPMAT[t,v] <- (beta_MAPMAT[t,v] > 0) * 1
      # Similarly for the global log-response ratios:
      Itr_star[t,v] <- (treat_star[t,v] > 0) * 1
    }
  }
  
  # Compute correlation matrix:
  # YL modified:
  for(t in 1:Nt){
    Sigma[t, 1:Nv,1:Nv] <- inverse(Omega[t, 1:Nv,1:Nv])
    # Fairly non-informative priors for precision matrix:
    # RO <- diag(c(0.1, ...,0.1)) and is created and read in as data in the workflow file.
    Omega[t, 1:Nv,1:Nv] ~ dwish(RO[1:Nv, 1:Nv], Nv)
    for(r in 1:(Nv-1)){
      # Among study standard deviations for each response variable:
      rho[t,r,r] <- sqrt(Sigma[t,r,r])
      for(c in (r+1):Nv){
        # Pairwise correlations amoung pairs of response variables:
        rho[t,r,c] <- Sigma[t,r,c]/sqrt(Sigma[t,r,r] * Sigma[t,c,c])
        rho[t,c,r] <- rho[t,r,c]
      }
    }
    # Among study std dev for SCE response variable:
    rho[t,Nv,Nv] <- sqrt(Sigma[t,Nv,Nv])
  }
  
  # Missing data models for Duration, MAP, and MAT; can be deleted if there is no missing data:
  for(r in 1:Nr){
    # Likelihood of record-level duration data
    Duration[r] ~ dnorm(mu.dur, tau.dur)
  }
  for(k in 1:Nk){
    # #YL test
    # alpha[k] <- df[k]
    # beta[k] <- df[k] * varR.study[k]
    # df[k] ~ dpois(df.all)
    # varR.study[k] ~ dgamma(0.01, 0.01) # not hier...
    
    # Likelihood of study-level MAT and MAP data
    MAT[k] ~ dnorm(mu.mat, tau.mat)
    MAP[k] ~ dnorm(mu.map, tau.map)
  }
  # priors for means and precisions (or std deviations) in above likelihoods
  df.all ~ dunif(0, 100)
  
  mu.varR ~ dlnorm(0,0.0001)
  tau.varR <- pow(sig.varR, -2)
  sig.varR ~ dunif(0,100)
  mu.dur ~ dnorm(0,0.0001)
  mu.map ~ dnorm(0,0.0001)
  mu.mat ~ dnorm(0,0.0001)
  tau.dur <- pow(sig.dur, -2)
  sig.dur ~ dunif(0,1000)
  tau.map <- pow(sig.map, -2)
  sig.map ~ dunif(0,1000)
  tau.mat <- pow(sig.mat, -2)
  sig.mat ~ dunif(0,1000)
  
  for(t in 1:Nt){
    lamNC[t] ~ dunif(0, 100)
    lamNT[t] ~ dunif(0, 100)
    sigXC[t] ~ dunif(0, 100)
    sigXT[t] ~ dunif(0, 100)    
  }
}
