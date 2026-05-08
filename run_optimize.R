##******************************************************************************
##* Script to run and/or optimize the model wrt a set of optimization parameters. 
##* Optimization is typically performed to align the model output to exogenous
##* data such as observed incidence of cancer and/or CVD.
##* 
##******************************************************************************
rm(list = ls())
# The model code is sourced. Defines simulate_model().
source("ncdsim.R")

start <- Sys.time()
outp <- simulate_model(startyear = 2022, endyear = 2026,cfact_startyear = 2025,nsteps = 50,
                       cfact_endyear = 2026,
                       is_baseline = T)
print(Sys.time() - start)

# The validation script is sourced. Defines validate_ncdsim().
source("validation.R")

# Validation of demographics, stocks and flows
validate_ncdsim(
  # path to NCDSim project
  projectroot = getwd(),
  # timestamp for simulation to validate
  timestamp = NCDSim_timestamp,
  use_scb_demographics=TRUE) 

