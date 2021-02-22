# Bayesian meta-analysis
Description:

0. JAGS must installed prior to running the R script: http://mcmc-jags.sourceforge.net/
1. Clone or download the repository to your computer.
2. Open the Rproject file ("meta-analysis.Rproj") with Rstudio.
3. Open and run the workflow file "workflow_meta_model_multivar.r". This will read the data files and source the Bayesian model ("Bayes_meta_model_multivar.r").
5. Output of the model run will be saved in the "./output/" folder.

Note that the Data folder contains two files: 
  - **Meta-analysis_record_level_input.csv** This is record-level data, which includes the variable ID (varID: 1 = AB, 2 = BB, and 3 = SCE), Study indicator (unique integers use to identify studies), the log-response ratio data (R), the assocaited variance of the log-response ratios (varR), the treatment type code (TrCode: 1 = CO2-only, 2 = warming-only, 3 = CO2  + warming), the System.Type (IGBP classification; not used in the analysis), the study duration (Duration, in years). In addition, "Site.Name" and "Source.Reference" entries are defined in Dieleman et al. (2012) GCB 18:2681-2693, in Tables S1 and S6, respectively, where full citations for the records can be found.
  - **Meta-analysis_study_level_input.csv** This file contains study level information, such that the row number aligns with the study number in the record-level data (i.e., study #1 gets study-level data in row 1); this files contains the climate (MAT and MAP) data.
