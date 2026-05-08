### Loading packages
library(data.table)
library(deSolve)
library(jsonlite)

source("calibrate.R")

##*****************************************************************************
## Miscellaneous functions
##*****************************************************************************

# Return indexed costs based on a KPP and a costs change rate for health care
health_costs <- function(index_year, endyear, costs_change_rate, kpp, stock) {
  n_years <- endyear - index_year
  if (n_years == 0) { 
    return(stock * kpp)
  } else {
    return(
      stock * kpp * (1 + sign(n_years) * costs_change_rate) ** (abs(n_years)))
  }
}


# piecewise function in 3 segments. Used to scale cFact
# Breakpoints at x0 and x1. Horizontal segment (y0) for x<= x0 and y1 for x > X1. 
piecewise <- function(x, x0, y0, x1, y1){
  a <- (y1 - y0) / (x1 - x0)
  b <- y0 - x0 * a
  f <- (x < x0) * y0 + (x >= x0 & x < x1) * (a * x + b) + (x >= x1) * y1
  return(f)
}

# Population Attributable fraction given a relative risk rr and a prevalence p
paf <- function(p, rr) {
  return((p * (rr - 1)) / (p * (rr - 1) + 1))
}

# Cumulative Population Attributable Fraction given an array of PAFs
cpaf <- function(pafs) {
  return(1 - apply(1 - pafs, 1, prod))
}

# read_data_parameters returns a list of data.tables containing the data and
# parameters of the model
read_data_parameters <- function(filepath = "") {
  
  # Mortality rates
  par_rate_dead <- fread(file = paste0(filepath, "r_dead.csv"), sep = ";")
  
  # Mortality rates
  par_rate_dead_all <- fread(file = paste0(filepath, "r_dead_all.csv"),
                             sep = ";")
  
  # Immigration rates
  par_rate_immig <- fread(file = paste0(filepath, "r_immigration.csv"),
                          sep = ";")
  
  # Emigration rates
  par_rate_emig <- fread(file = paste0(filepath, "r_emigration.csv"),
                         sep = ";")
  
  # Birth rates
  par_rate_born <- fread(file = paste0(filepath, "r_born.csv"), sep = ";")
  
  getrates <- function(filename) {
    tmp <- fread(file = paste0(filepath, filename), sep = ";")
    tmp2 <- melt(tmp, id.vars = c("sex", "age"))
    setnames(tmp2, 
             old = c("variable", "value"),
             new = c("year", "rate"))
    return(tmp2[, .(year, sex, age, rate)])
  }
  excess_mortality <- fread(file = paste0(filepath, "excess_mortality_disease.csv"), sep = ";")
  
  get_prevalence <- function(filename) {
    tmp <- fread(file = paste0(filepath, filename), sep = ";")
    cnames <- setdiff(colnames(tmp), c("sex", "age"))
    tmp[, (cnames) := lapply(.SD, as.double), .SDcols = cnames]
    tmp2 <- melt(tmp, id.vars = c("sex", "age"))
    tmp2[, ":="(
      year = variable, 
      prev = as.numeric(value))]
    return(tmp2[, .(year, sex, age, prev)])
  }
  
  ## Risk factor prevalence: smoking
  par_prev_smoking <- get_prevalence("prev_smoking.csv")
  
  ## Risk factor prevalence: alcohol
  par_prev_alcohol <- get_prevalence("prev_alcohol.csv")
  
  ## Risk factor prevalence: physical inactivity
  par_prev_inactivity <- get_prevalence("prev_inactivity.csv")
  
  ## Risk factor prevalence: BMI
  par_prev_bmi  <- get_prevalence("prev_bmi.csv")
  
  # Population counts (2000 - 2120), observed & predicted
  dat_popcounts <- fread(file = paste0(filepath, "pop_counts_scb.csv"), 
                         sep = ";")
  
  # Read prevalences and relative risks related to diet
  diet_prev <- fread(file = paste0(filepath, "prev_diet.csv"), sep = ";")
  
  # Prevalences (for validation and initialization of stock)
  preval_new <- fread(paste0(filepath, "stocks.csv"), sep = ';')
  preval_new <- dcast(preval_new, year+sex+age ~ stock, value.var = 'N')
  colnames(preval_new) <- c("year","sex", "age","prev_cancer","prev_comorb", "prev_cvd" ,"totpop")
  preval <- CJ(year = min(preval_new[, (year)]):max(preval_new[, (year)]),
               sex = 1:2, age = 0:100)
  preval <- merge(preval, 
                  preval_new[, c("year","age", "sex","prev_cancer","prev_comorb", "prev_cvd")], 
                  by=c("year","age", "sex" ), all=TRUE)

  preval <- melt(preval, id.vars = c("sex", "age","year"), 
                 value.vars = c("prev_cancer","prev_comorb", "prev_cvd"))

  colnames(preval) <- c("sex", "age","year", "ncd", "prevalence")
  preval <- setnafill(preval, "const", fill=0, cols=c("prevalence"))
  par_prevalence_cancer <- preval[ncd=="prev_cancer", .(year, sex, age, prevalence)]
  par_prevalence_cvd <- preval[ncd=="prev_cvd", .(year, sex, age, prevalence)]
  par_prevalence_comorb  <- preval[ncd=="prev_comorb", .(year, sex, age, prevalence)]
  
  return(list(rate_dead = par_rate_dead,
              rate_dead_all = par_rate_dead_all, 
              rate_immig = par_rate_immig, 
              rate_emig = par_rate_emig, 
              rate_born = par_rate_born, 
              excess_mortality_cancer = excess_mortality[stock=='cancer', c("age", "sex","year", "rate")], 
              excess_mortality_cvd = excess_mortality[stock=='cvd', c("age", "sex","year", "rate")], 
              excess_mortality_comorb = excess_mortality[stock=='comorb', c("age", "sex","year", "rate")], 
              prev_smoking = par_prev_smoking, 
              prev_alcohol = par_prev_alcohol, 
              prev_inactivity = par_prev_inactivity, 
              prev_bmi = par_prev_bmi, 
              popcounts = dat_popcounts,
              prevalence_cancer = par_prevalence_cancer,
              prevalence_cvd = par_prevalence_cvd,
              prevalence_comorb = par_prevalence_comorb,
              prev_diet = diet_prev))
}


##****************************************************************************
## Main loop: year, sex, age
##****************************************************************************

simulate_model <- function(
    # Indicates if baseline scenario is simulated.
  # If true calibrate incidence rate and mortality
  is_baseline = TRUE,
  # The root file path of the project
  PROJECTROOT = getwd(),
  # path to the json file which contains the baseline parameters 
  # such as relative risks etc.
  # if is_baseline = FALSE a .csv file with the same stem name
  # is needed in the same folder.
  # It is saved automatically when a baseline scenario is simulated
  baseline_parameters_path = paste0(PROJECTROOT, 
                                    "/Input/defaults/baseline.json"),
  input_path = paste0(PROJECTROOT, "/Input/"),
  
  # Base year of simulation. NOTE: the year before the first year of simulation
  startyear = 2015,
  
  # The last year of simulation
  endyear = 2030,
  
  
  
  # Indicates whether output data is written to file
  write_data_to_file = TRUE,
  
  # Use when run by Shiny app
  ui = FALSE,
  
  ui_lang = "eng",
  
  ## Use when UI = TRUE to keep track of out data for each scenario run
  scen_name = NULL,
  
  ## Use when UI = TRUE for session specific path for scenario files
  scen_path = NULL,
  
  # Number of (equidistant) steps within each year for ode()
  nsteps = 20,
  
  # F: normal simulation, T: only demographic simulation (no NCD)
  only_demog = FALSE,
  
  # F: normal output, T: output of debug info (ode: stock/flow)
  debug_ode = FALSE,
  
  # Timestamp for identification of simulation output
  timestamp = gsub("[.]", "", format(Sys.time(), "%Y_%m_%d_%H_%M_%OS2")),
  
  # Adjustment of risk factor prevalences to implement scenarios
  # To do: in the GUI the sliders should control the population average of 
  # prevalences given the sex/age-specific prevalences and the population (SCB) 
  # of some year of reference (the current year?). In this version the
  # adjustment factors denotes the proportional change.
  cfact = c(cfact_smoking = 1.0,
            cfact_alcohol = 1.0,
            cfact_inactivity = 1.0,
            cfact_bmi = 1.0),
  
  cfact_food = c(fruit = 1.0,
                 wholegrains = 1.0,
                 greens = 1.0,
                 meat = 1.0, 
                 salt = 1.0),
  
  # First year of applying interventions, i.e. cfact that is not equal to 1.0
  cfact_startyear = 2024,
  
  # Last year of onset for applying interventions, i.e. 
  # cfact now is equal to cfact.
  # linear decrease
  cfact_endyear = 2025
) {
  
  # Assign timestamp to variable in .GlobalEnv to use in validation
  assign("NCDSim_timestamp", timestamp, envir = .GlobalEnv)
  
  ## Write to log when running from Shiny
  if (ui == TRUE) {
    cat(
      switch(ui_lang,
             "eng" = paste0("Starting simulation ", scen_name,
                            " from ", startyear, " to ", endyear, ".\n"), 
             "swe" = paste0("Startar simulering ", scen_name,
                            " från ", startyear, " till ", endyear, ".\n") 
      ),
      file=paste0(scen_path, "log.txt"))
  }
  
  ##*****************************************************************************
  ##* Create start data
  ##*****************************************************************************
  
  ## Path to /Input directory including files with risk factor prevalences etc
  
  
  # Read the input json file and assign parameters
  baseline_parameters <- read_json(baseline_parameters_path)
  calibrate_cpaf_cancer <- as.integer(baseline_parameters$calibrate_cpaf_cancer)
  calibrate_cpaf_cvd <- as.integer(baseline_parameters$calibrate_cpaf_cvd)
  cpaf_cancer <- as.numeric(baseline_parameters$cpaf_cancer)
  cpaf_cvd <- as.numeric(baseline_parameters$cpaf_cvd)
  age_cutoff_cancer <- as.integer(baseline_parameters$age_cutoff_cancer)
  age_cutoff_cvd <- as.integer(baseline_parameters$age_cutoff_cvd)
  age_cutoff_cancer_high <- as.integer(baseline_parameters$age_cutoff_cancer_high)
  age_cutoff_cvd_high <- as.integer(baseline_parameters$age_cutoff_cvd_high)
  dcost_total_cancer <- as.numeric(baseline_parameters$dcost_total_cancer)
  icost_total_cancer <- as.numeric(baseline_parameters$icost_total_cancer)
  dcost_growth_cancer_input <- as.numeric(
    baseline_parameters$dcost_growth_cancer)
  dcost_total_cvd <- as.numeric(baseline_parameters$dcost_total_cvd)
  
  icost_total_cvd <- as.numeric(baseline_parameters$icost_total_cvd)
  dcost_growth_cvd_input <- as.numeric(baseline_parameters$dcost_growth_cvd)
  dcost_total_base_year <- as.integer(baseline_parameters$dcost_total_base_year)
  rr <- unlist(baseline_parameters$rr)
  rr_cancer_diet <- unlist(baseline_parameters$rr_cancer_diet)
  rr_cvd_diet <- unlist(baseline_parameters$rr_cvd_diet)
  communalities <- unlist(baseline_parameters$communalities)
  direct_costs <- c(cancer = dcost_total_cancer, cvd = dcost_total_cvd)
  indirect_costs <- c(cancer = icost_total_cancer, cvd = icost_total_cvd)
  
  # Perform calibration on the go
  # Run calibration and assign variables. Read into calibrate.R for details
  if (is_baseline == TRUE) {
    calibration_results <- calibrate(input_path=input_path,
                                     calibrate_cpaf_cancer = calibrate_cpaf_cancer,
                                     calibrate_cpaf_cvd = calibrate_cpaf_cvd,
                                     prop_cancer=cpaf_cancer,
                                     prop_cvd=cpaf_cvd, 
                                     age_cutoff_cancer=age_cutoff_cancer, 
                                     age_cutoff_cvd=age_cutoff_cvd,
                                     age_cutoff_cancer_high=age_cutoff_cancer_high, 
                                     age_cutoff_cvd_high=age_cutoff_cvd_high,
                                     rr=rr, 
                                     rr_cancer_diet=rr_cancer_diet,
                                     rr_cvd_diet=rr_cvd_diet,
                                     communalities=communalities,
                                     direct_costs=direct_costs,
                                     indirect_costs=indirect_costs,
                                     dcost_total_base_year=dcost_total_base_year
    )
    
    const <- calibration_results$constants
    paf_other_calibration_0 <- calibration_results$paf_other
    heal_rate <- calibration_results$heal_rate
    dcost_unit <- calibration_results$dcost_unit
    # _ to avoid error with casting to the data table
    icost_unit_ <- calibration_results$icost_unit

    # pad paf_other if the simulation period is longer than the calculated values.
    # If the startyear < first year in paf_other: pad years between with paf(first_year)
    # if last year of simulation is > last_year in paf_other: pad years between with paf(last_year)
    first_year <- min(as.integer(paf_other_calibration_0[, year]))
    last_year <- max(as.integer(paf_other_calibration_0[, year]))
    
    if (startyear < first_year){
      tmp <- paf_other_calibration_0[year==first_year]
      for (y in startyear:(first_year -1)){
        tmp[,year:=y]
        paf_other_calibration_0 <- rbindlist(list(tmp,paf_other_calibration_0),
                                             fill=TRUE)
      }
    }
    if (endyear > last_year) {
      tmp <- paf_other_calibration_0[year==last_year]
      for (y in (last_year + 1):endyear){
        tmp[,year:=y]
        paf_other_calibration_0 <- rbindlist(list(paf_other_calibration_0,tmp),
                                             fill=TRUE)
      }
    }


    calibration_mortality_table <- CJ(year = startyear:endyear, sex = 1:2,
                                      age = 0:100, adj = 0.0)
    
  } 
  
  # If not a baseline run, then read calibration parameters
  # from the correspondent baseline run and format accordingly
  
  else{
    cal_param_file <- paste0(strsplit(baseline_parameters_path, ".json"),".csv")
    calibrated_parameters <- fread(cal_param_file, sep=";")
    calibration_mortality_table <- CJ(
      year = startyear:endyear, age = 0:100, sex = 1:2)
    calibration_mortality_table <- merge(calibration_mortality_table, 
                                         calibrated_parameters[
                                           , .(year, age, sex, calibration_mortality)], 
                                         by = c("year", "age", "sex" ), all = FALSE)
    colnames(calibration_mortality_table) <- c("year","age","sex", "adj")
    
    paf_other_calibration_0 <- calibrated_parameters[,c("sex", "age", "year",
                                                        "paf_cancer_other",
                                                        "paf_cvd_other",
                                                        "paf_cancer_comorb_other",
                                                        "paf_cvd_comorb_other",
                                                        "paf_comorb_other")]
    paf_other_calibration_0 <- setnafill(paf_other_calibration_0, "const",
                                         fill = 0, cols = c("paf_cancer_other",
                                                            "paf_cvd_other",
                                                            "paf_cancer_comorb_other",
                                                            "paf_cvd_comorb_other",
                                                            "paf_comorb_other"))
    colnames(paf_other_calibration_0) <- c("sex", "age", "year", "cpaf_cancer",
                                           "cpaf_cvd", "cpaf_cancer_comorb",
                                           "cpaf_cvd_comorb","cpaf_comorb")
    paf_other_calibration_0 <-  melt(paf_other_calibration_0,
                                     id.vars = c("sex", "age", "year"),
                                     variable.name = 'disease',
                                     value.name = 'cpaf' )
    
    
    const <- calibrated_parameters[year == as.integer((startyear + endyear) / 2),
                                   c("sex", "age", "calibration_cancer",
                                     "calibration_cvd",
                                     "calibration_cancer_comorb",
                                     "calibration_cvd_comorb",
                                     "calibration_comorb")]
    const <- setnafill(const, "const", fill = 0, cols = c("calibration_cancer",
                                                          "calibration_cvd",
                                                          "calibration_cancer_comorb",
                                                          "calibration_cvd_comorb",
                                                          "calibration_comorb"))
    colnames(const) <- c("sex", "age", "const_cancer", "const_cvd",
                         "const_cancer_comorb","const_cvd_comorb",
                         "const_comorb")
    const <- melt(const, id.vars = c("sex", "age"),
                  variable.name = 'disease', value.name = 'constant' )
    
    
    heal_rate <- calibrated_parameters[year == as.integer((startyear + endyear) / 2),
                                       c("sex", "age", "heal_cvd_pop",
                                         "heal_cancer_pop", "heal_comorb_pop",
                                         "heal_comorb_cancer", "heal_comorb_cvd")]
    heal_rate <- setnafill(heal_rate, "const", fill = 0.0,
                           cols = c("heal_cvd_pop", "heal_cancer_pop",
                                    "heal_comorb_pop", "heal_comorb_cancer",
                                    "heal_comorb_cvd"))
    colnames(heal_rate) = c("sex", "age", "cvd", "cancer", "comorb_healthy",
                            "comorb_cancer", "comorb_cvd")
    heal_rate <- melt(heal_rate, id.vars = c("sex", "age"),
                      variable.name = 'disease', value.name = 'heal_rate')
    
    dcost_unit <- CJ(sex = 1:2, age = 0:100)
    dcost_unit <- merge(dcost_unit, calibrated_parameters[
      year==as.integer((startyear + endyear) / 2),
      .(age, sex, dcost_unit_cancer, dcost_unit_cvd)], 
      by = c("age", "sex" ), all = TRUE)
    colnames(dcost_unit) <- c("age", "sex", "cancer", "cvd")
    dcost_unit <- melt(dcost_unit, id.vars=c("age", "sex"))
    colnames(dcost_unit) <- c("age", "sex", "ncd", "dcost_unit")
    # _ to avoid error with casting to the data table
    icost_unit_ <- CJ(sex = 1:2, age = 0:100)
    icost_unit_ <- merge(icost_unit_, 
                         calibrated_parameters[
                           year==as.integer((startyear + endyear) / 2), 
                           .(age, sex, icost_unit_cancer, icost_unit_cvd)],
                         by = c("age", "sex" ), all=TRUE)
    colnames(icost_unit_) <- c("age","sex","cancer", "cvd")
    icost_unit_ <- melt(icost_unit_, id.vars=c("age", "sex"))
    colnames(icost_unit_) <- c("age","sex","ncd","icost_unit")
    

  }
  

  
  ## Get data and parameters
  dat_par <- read_data_parameters(input_path)
  max_year_mortality <- max(dat_par$excess_mortality_comorb$year)
  
  # extract basic values for diet
  food_cat <- colnames(dat_par$prev_diet)[3:length(colnames(dat_par$prev_diet))]
  rr_cancer_diet <- data.frame(as.list(rr_cancer_diet))
  rr_cvd_diet <- data.frame(as.list(rr_cvd_diet))
  
  
  # Creating DT for storage of initial and simulated population
  # Note: dat is initialized to contain all data from the simulation. This 
  #       prevents the model from having to concatenate data (and hence allocate
  #       new memory) within the model loop which would be very time consuming.
  dat <- CJ(year = startyear:endyear, sex = 1:2, age = 0:100, s_pop = 0.0)
  
  
  # Merge initial population to dat
  dat[dat_par$popcounts[year == startyear], 
      on = .(year == year, sex == sex, age == age), 
      s_pop_ := i.pop]
  
  # Creating empty variables for flows
  dat[, c("f_dead", "f_born", "f_immig_pop", "f_immig_cvd",
          "f_immig_cancer", "f_immig_comorb", 
          "f_pop_cancer","f_cancer_pop", "f_pop_cvd", "f_cvd_pop", "f_pop_comorb",
          "f_comorb_pop", "f_cancer_comorb", "f_comorb_cancer",
          "f_cvd_comorb", "f_comorb_cvd","f_cancer_dead",
          "f_cvd_dead", "f_comorb_dead") := as.numeric(NA)]
  
  ### CREATE empty variables for costs:
  dat[, c("dcost_unit_cancer", "icost_unit_cancer",
          "dcost_unit_cvd", "icost_unit_cvd",
          "dcost_cancer", "dcost_cvd", "icost_cancer",
          "icost_cvd", "dcost_growth_cancer",
          "dcost_growth_cvd") := as.numeric(NA)]
  
  dat[, (paste("prev", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_cancer", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_cvd", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_cancer_comorb", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_cvd_comorb", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_comorb", food_cat, sep = "_")) := as.numeric(NA)]
  

  
  ### Costs change rates can be further customized per year, sex and disease
  
  dat[, dcost_growth_cancer := dcost_growth_cancer_input]
  dat[, dcost_growth_cvd := dcost_growth_cvd_input]
  
  # Merge prevalences of cancer and CVD (by reference) to dat
  dat[dat_par$prevalence_cancer[year == startyear], 
      on = .(year == year, sex == sex, age == age),
      s_cancer := as.numeric(i.prevalence)]
  
  dat[dat_par$prevalence_cvd[year == startyear], 
      on = .(year == year, sex == sex, age == age),
      s_cvd := as.numeric(i.prevalence)]
  dat[dat_par$prevalence_comorb[year == startyear], 
      on = .(year == year, sex == sex, age == age),
      s_comorb := as.numeric(i.prevalence)]
  
  # Subtracting stocks of NCD from (health) pop
  # Note: if s_pop becomes negative s_pop_cancer and s_pop_cvd are reduced
  #       proportionally.
  dat[year == startyear, s_pop := s_pop_ - s_cancer - s_cvd - s_comorb]
  idx <- which(dat[, s_pop < 0])
  dat[idx, pct := 1 - abs(s_pop) / (s_cancer + s_cvd + s_comorb)]
  dat[idx, ":="(
    s_cancer = pct * s_cancer,
    s_cvd = pct * s_cvd,
    s_comorb = pct * s_comorb,
    s_pop = 0
  )]
  
  dat[, c("s_pop_", "pct") := NULL]
  
  # Add stocks of dead from NCD:s (initially empty)
  dat[, ":="(
    s_dead_cancer = 0,
    s_dead_cvd = 0,
    s_dead_comorb = 0
  )]
  
  # Initializing stock of (accumulated) dead
  dat[, s_dead := 0]
  
  ##****************************************************************************
  
  ## Demography only mode adjustments to relative risks and data
  if (only_demog == TRUE) {
    rr[] <- 0 ## Set all relative risks to zero
    dat[year == startyear, s_pop := s_pop + s_cancer + s_cvd]
    dat[year == startyear, c("s_cancer", "s_cvd") := as.numeric(0)]
  }
  
  model <- function(time, stocks, auxs){
    with(as.list(c(stocks, auxs)), { 
      
      ## Calculate components of demographic flows
      
      dr <- unname(pmin(1.0, r_dead))
      f_dead <- s_pop * dr
      f_born <- n_born
      f_immig_pop <- s_pop * r_immig
      f_immig_cvd <- s_cvd * r_immig
      f_immig_cancer <- s_cancer * r_immig
      f_immig_comorb <- s_comorb * r_immig
      
      f_emig_pop <- s_pop * pmin(1.0, r_emig)
      f_emig_cvd <- s_cvd * pmin(1.0, r_emig)
      f_emig_cancer <- s_cancer * pmin(1.0, r_emig)
      f_emig_comorb <- s_comorb * pmin(1.0, r_emig)
      
      ##** Calculating flows for cancer
      # Getting sick. Use cPAF functions
      paf_cancer_smoking <- paf(prev_smoking, rr_cancer_smoking)
      paf_cancer_inactivity <- paf(prev_inactivity, rr_cancer_inactivity)
      paf_cancer_bmi <- paf(prev_bmi, rr_cancer_bmi)
      paf_cancer_alcohol <- paf(prev_alcohol, rr_cancer_alcohol)
      
      p_pop_cancer_0_nm <- 1 - (1 - comm_smoking * paf_cancer_smoking) *
        (1 - comm_inactivity * paf_cancer_inactivity) *
        (1 - comm_bmi * paf_cancer_bmi) *
        (1 - comm_alcohol * paf_cancer_alcohol)
      
      p_pop_cancer_0 <- (1 - (1 - cpaf_diet_cancer) * (1 - p_pop_cancer_0_nm))
      p_pop_cancer <- p_pop_cancer_0 * calibration_cancer + paf_cancer_other
      
      f_pop_cancer <- s_pop * p_pop_cancer
      
      
      # Getting well (five years with cancer)
      f_cancer_pop <- s_cancer * heal_cancer_pop # change from 1/5
      
      # Dying
      dr_cancer <- unname(pmin(1.0, r_dead_cancer))
      f_cancer_dead <- s_cancer * dr_cancer
      
      ##** Calculating flows for CVD
      # Getting sick
      paf_cvd_smoking <-  paf(prev_smoking, rr_cvd_smoking)
      paf_cvd_inactivity <- paf(prev_inactivity, rr_cvd_inactivity)
      paf_cvd_bmi <-paf(prev_bmi, rr_cvd_bmi)
      paf_cvd_alcohol <- paf(prev_alcohol, rr_cvd_alcohol)
      
      p_pop_cvd_0_nm <- (1 - (1 - comm_smoking * paf_cvd_smoking) *
                           (1 - comm_inactivity * paf_cvd_inactivity) *
                           (1 - comm_bmi * paf_cvd_bmi) *
                           (1 - comm_alcohol * paf_cvd_alcohol))
      p_pop_cvd_0 <- (1 - (1 - cpaf_diet_cvd) * (1 - p_pop_cvd_0_nm))
      
      p_pop_cvd <- p_pop_cvd_0 * calibration_cvd  + paf_cvd_other
      
      f_pop_cvd <- s_pop * p_pop_cvd
      
      
      # Estimated healing time for cvd from an optimization 
      # between observed and calculated over years 2009-2022
      f_cvd_pop <- s_cvd * heal_cvd_pop # change from (1/11).
      
      # Dying
      dr_cvd = unname(pmin(1.0, r_dead_cvd))
      f_cvd_dead <- s_cvd * dr_cvd
      
      ##** Calculating flows for Cancer and CVD # paf_comorb=max(paf_cvd, paf_ccr))
      # Getting sick
      
      p_pop_comorb_0 <- max(p_pop_cancer_0, p_pop_cvd_0)
      p_pop_comorb <- p_pop_comorb_0 * calibration_comorb + paf_comorb_other
      
      f_pop_comorb <- s_pop * p_pop_comorb
      f_cancer_comorb <- s_cancer * (p_pop_cvd_0 * calibration_cancer_comorb +
                                       paf_cancer_comorb_other)
      f_cvd_comorb <- s_cvd * (p_pop_cancer_0 * calibration_cvd_comorb +
                                 paf_cvd_comorb_other)
      
      # Estimated healing time for cancer and cvd 
      
      f_comorb_pop <- s_comorb * heal_comorb_pop  
      
      f_comorb_cancer <- s_comorb * heal_comorb_cancer
      
      f_comorb_cvd <- s_comorb * heal_comorb_cvd
      
      # Dying
      dr_comorb = unname(pmin(1.0, r_dead_comorb))
      f_comorb_dead <- s_comorb * dr_comorb
      
      ##
      ## Calculating differential equations
      ##
      
      # Equation for stock of non-NCD population
      d_s_pop_dt <- f_born - f_dead + f_immig_pop - f_emig_pop + f_cvd_pop - 
        f_pop_cvd + f_cancer_pop - f_pop_cancer + f_comorb_pop - f_pop_comorb
      
      # Equation for stock of population with cancer 
      d_s_cancer_dt <- f_pop_cancer - f_cancer_pop + f_immig_cancer - 
        f_emig_cancer - f_cancer_dead + f_comorb_cancer - f_cancer_comorb
      
      # Equation for stock of population with CVD
      d_s_cvd_dt <- f_pop_cvd - f_cvd_pop + f_immig_cvd - f_emig_cvd - 
        f_cvd_dead + f_comorb_cvd - f_cvd_comorb
      
      # Equation for stock of population with CVD and cancer
      d_s_comorb_dt <- f_pop_comorb - f_comorb_pop + f_immig_comorb - f_emig_comorb-
        f_comorb_dead + f_cancer_comorb - f_comorb_cancer + f_cvd_comorb - f_comorb_cvd 
      
      # Equation for stock of dead from non-NCD
      d_s_dead_dt <- f_dead
      
      # Equation for stock of dead from cancer
      d_s_dead_cancer_dt <- f_cancer_dead
      
      # Equation for stock of dead from CVD
      d_s_dead_cvd_dt <- f_cvd_dead
      
      # Equation for stock of dead from Comorbiditet
      d_s_dead_comorb_dt <- f_comorb_dead
      
      
      return(
        list(c(
          
          # Stock differentials
          d_s_pop_dt, d_s_cancer_dt, d_s_cvd_dt, d_s_comorb_dt, d_s_dead_dt,
          d_s_dead_cancer_dt, d_s_dead_cvd_dt, d_s_dead_comorb_dt),
          
          # Flows
          f_dead = f_dead, f_born = f_born, 
          
          f_immig_pop = f_immig_pop, f_immig_cvd = f_immig_cvd, 
          f_immig_cancer = f_immig_cancer, f_immig_comorb = f_immig_comorb,
          
          f_emig_pop = f_emig_pop, f_emig_cvd = f_emig_cvd, 
          f_emig_cancer = f_emig_cancer, f_emig_comorb = f_emig_comorb,
          
          f_pop_cancer = f_pop_cancer, f_cancer_pop = f_cancer_pop, 
          f_cancer_dead = f_cancer_dead,
          
          f_pop_cvd = f_pop_cvd, f_cvd_pop = f_cvd_pop,
          f_cvd_dead = f_cvd_dead, 
          
          f_pop_comorb = f_pop_comorb,
          f_comorb_pop = f_comorb_pop,
          f_cvd_comorb = f_cvd_comorb,
          f_comorb_cvd =f_comorb_cvd,
          f_cancer_comorb = f_cancer_comorb,
          f_comorb_cancer = f_comorb_cancer,
          f_comorb_dead = f_comorb_dead,
          
          # Assumptions
          dr = dr,
          dr_cancer = dr_cancer,
          dr_cvd = dr_cvd,
          dr_comorb = dr_comorb,
          
          heal_cvd_pop = heal_cvd_pop, 
          heal_cancer_pop = heal_cancer_pop,
          heal_comorb_pop = heal_comorb_pop,
          heal_comorb_cancer = heal_comorb_cancer,
          heal_comorb_cvd = heal_comorb_cvd,
          
          prev_alcohol = prev_alcohol,
          prev_bmi = prev_bmi,
          prev_smoking = prev_smoking,
          prev_inactivity = prev_inactivity,
          
          paf_cvd_smoking = paf_cvd_smoking,
          paf_cvd_inactivity = paf_cvd_inactivity, 
          paf_cvd_bmi = paf_cvd_bmi,
          paf_cvd_alcohol = paf_cvd_alcohol,
          cpaf_diet_cvd = cpaf_diet_cvd,
          p_pop_cvd = p_pop_cvd,
          
          paf_cancer_smoking = paf_cancer_smoking,
          paf_cancer_inactivity = paf_cancer_inactivity,
          paf_cancer_bmi = paf_cancer_bmi,
          paf_cancer_alcohol = paf_cancer_alcohol,
          cpaf_diet_cancer = cpaf_diet_cancer,
          p_pop_cancer = p_pop_cancer,
          
          
          calibration_cancer = calibration_cancer,
          calibration_cvd = calibration_cvd,
          calibration_comorb = calibration_comorb,
          calibration_cancer_comorb = calibration_cancer_comorb,
          calibration_cvd_comorb = calibration_cvd_comorb,
          
          paf_cvd_other = paf_cvd_other,
          paf_cancer_other = paf_cancer_other,
          paf_comorb_other = paf_comorb_other,
          paf_cancer_comorb_other = paf_cancer_comorb_other,
          paf_cvd_comorb_other = paf_cvd_comorb_other
          
        ))
    })
  }
  
  # Loop over years
  for (starttime in (startyear + 1):endyear) {

    
    ## Check whether to stop simulation when running from Shiny
    if (ui == TRUE) {
      t0 <- Sys.time()
      ui_status <- scan(paste0(scen_path, "ui_status.txt"), what="character")
      if (ui_status == "stop") {
        cat(switch(ui_lang,
                   "eng" = "Simulation stopped\n",
                   "swe" = "Simulering stoppad\n"
        ), 
        file=paste0(scen_path, "log.txt"),
        append = TRUE)
        stop("Stopped from ui")
      }
      
      ## Write current simulation year to file
      cat(starttime, file = paste0(scen_path, "current_simyear.txt"))
      
      ## Write to log
      cat(switch(ui_lang,
                 "eng" = paste0("Year ", starttime, " started. "), 
                 "swe" = paste0("År ", starttime, " påbörjat. ")
      ),
      file=paste0(scen_path, "log.txt"), append = TRUE)
    }
    
    print(paste0("year: ", starttime))
    
    # Total population: for calculation of nr. of births
    n_pop <- dat[year == (starttime - 1) & age > 19 & age < 40,
                 sum(s_pop + s_cancer + s_cvd + s_comorb)]
    
    # Loop over sex
    for (sex_ in unique(dat$sex)) {
      
      # Loop over age
      for (age_ in unique(dat$age)) {
        
        ## Check whether to stop simulation when running from Shiny
        if (ui == TRUE) {
          ui_status <- scan(paste0(scen_path, "ui_status.txt"), what="character")
          if (ui_status == "stop") {
            cat("Simulation stopped\n", file=paste0(scen_path, "log.txt"),
                append = TRUE)
            stop("Stopped from ui")
          }
        }
        
        # Create the start time, finish time, and time step
        START <- (starttime - 1)
        FINISH <- START + 1
        STEP <- 1 / nsteps # Time-steps within each year (arbitrary)
        
        # Create time vector
        simtime <- seq(START, FINISH, by = STEP)
        
        # Defining stock vectors for call to ode()
        # Ageing of population by moving forward data from previous year/age
        #       within sex
        if (age_ == 0) {
          yvec <- c(s_pop = 0.0, s_cancer = 0.0, s_cvd = 0.0, s_dead = 0.0, 
                    s_dead_cancer = 0.0, s_dead_cvd = 0.0)
          yvec <- c(s_pop = 0.0, s_cancer = 0.0, s_cvd = 0.0, s_comorb = 0.0, s_dead = 0.0, 
                    s_dead_cancer = 0.0, s_dead_cvd = 0.0, s_dead_comorb = 0.0)
        } else if (age_ == 100) {
          
          # Age 100+ is an open age group which leads to an accumulation of
          #       counts
          yvec <- unlist(dat[year == START & sex == sex_ & age == (age_ - 1), 
                             .(s_pop, s_cancer, s_cvd, s_comorb, s_dead,
                               s_dead_cancer, s_dead_cvd, s_dead_comorb)]) +
            unlist(dat[year == START & sex == sex_ & age == age_, 
                       .(s_pop, s_cancer, s_cvd, s_comorb, s_dead,
                         s_dead_cancer, s_dead_cvd, s_dead_comorb)])
        } else {
          yvec <- unlist(dat[year == START & sex == sex_ & age == (age_ - 1), 
                             .(s_pop, s_cancer, s_cvd, s_comorb, s_dead,
                               s_dead_cancer, s_dead_cvd, s_dead_comorb)])
        }
        
        # Defining vector of auxiliary parameters: avec
        
        # r_dead = mortality rate for the healthy population
        # r_dead_all = mortality in all the population
        # (used only for calibration purposes)
        r_dead <- dat_par$rate_dead[year == fifelse(FINISH > max_year_mortality,
                                                    max_year_mortality, FINISH)
                                    & sex == sex_ & age == age_, 
                                    rate]
        # In a few cases the observed death rates for children are zero. Since
        # this will result in a zero denominator when doing the proportional
        # adjustment of death rates below we replace the death rate with a very
        # small number
        if (r_dead == 0) r_dead <- 1e-10
        
        # All dead
        r_dead_all <- dat_par$rate_dead_all[year == FINISH &
                                              sex == sex_ & age == age_, rate]
        if (r_dead_all == 0) r_dead_all <- 1e-10
        
        # Calculate death rates for cancer and CVD from general death rates
        # and excess mortality
        excess_mort_cancer <- dat_par$excess_mortality_cancer[
          year == fifelse(FINISH > max_year_mortality, max_year_mortality,
                          FINISH) & sex == sex_ & age == age_, rate]
        r_dead_cancer <- pmin(0.999, r_dead * excess_mort_cancer)
        
        excess_mort_cvd <- dat_par$excess_mortality_cvd[
          year == fifelse(FINISH > max_year_mortality, max_year_mortality,
                          FINISH) & sex == sex_ & age == age_, rate]
        r_dead_cvd <- pmin(0.999, r_dead * excess_mort_cvd)
        
        excess_mort_comorb <- dat_par$excess_mortality_comorb[
          year == fifelse(FINISH > max_year_mortality, max_year_mortality, 
                          FINISH) & sex == sex_ & age == age_, rate]
        r_dead_comorb <- pmin(0.999, r_dead * excess_mort_comorb)

        # Calibration constants for other factors:
        calibration_cancer <- const[sex == sex_ & age == age_ & disease=="const_cancer",
                                    constant]
        calibration_cvd <- const[sex == sex_ & age == age_ & disease=="const_cvd",
                                 constant]
        calibration_comorb <- const[sex == sex_ & age == age_ & disease=="const_comorb",
              constant]
        calibration_cancer_comorb <- const[sex == sex_ & age == age_ &
                                             disease=="const_cancer_comorb",
                                           constant]
        calibration_cvd_comorb <- const[sex == sex_ & age == age_ &
                                          disease=="const_cvd_comorb",
                                           constant]
        
        
        paf_cancer_other <- paf_other_calibration_0[year == FINISH & sex == sex_ &
                                                 age == age_ & disease=='cpaf_cancer', cpaf]
        paf_cvd_other <- paf_other_calibration_0[year == FINISH & sex == sex_ &
                                           age == age_ & disease == 'cpaf_cvd', cpaf]
        
        paf_comorb_other <- paf_other_calibration_0[year == FINISH & sex == sex_ &
                                           age == age_ & disease == 'cpaf_comorb', cpaf]
        paf_cancer_comorb_other <- paf_other_calibration_0[year == FINISH & sex == sex_ &
                                                      age == age_ & disease == 'cpaf_cancer_comorb', cpaf]
        paf_cvd_comorb_other <- paf_other_calibration_0[year == FINISH & sex == sex_ &
                                                             age == age_ & disease == 'cpaf_cvd_comorb', cpaf]

        heal_cvd_pop <- heal_rate[sex == sex_ &
                                   age == age_ & disease == "cvd", heal_rate]
        heal_cancer_pop <- heal_rate[sex == sex_ &
                                      age == age_ & disease == "cancer", heal_rate]
        heal_comorb_pop <- heal_rate[sex == sex_ &
                                      age == age_ & disease == "comorb_healthy", heal_rate]
        heal_comorb_cancer <- heal_rate[sex == sex_ &
                                         age == age_ & disease == "comorb_cancer", heal_rate]
        heal_comorb_cvd <- heal_rate[sex == sex_ &
                                      age == age_ & disease == "comorb_cvd", heal_rate]
        
        
        r_immig <- dat_par$rate_immig[year == FINISH & sex == sex_ &
                                        age == age_, rate]
        r_emig <- dat_par$rate_emig[year == FINISH & sex == sex_ &
                                      age == age_, rate]
        
        # Alignment of death rates for the stock of diagnosed and non-diagnosed
        # to ensure total population death rate according to Statistics Sweden 
        # (observed and forecasted) while maintaining the estimated relative 
        # death rate ratio between the diagnosed and non-diagnosed.
        # The initial mortality rates of the diagnosed are calculated using 
        # empirically estimated excess mortality rates, compared to the non-
        # diagnosed population. For non-diagnosed, which is a larger population, 
        # death rates have been directly estimated.
        # The death rates within each stock are adjusted proportionally to 
        # achieve the correct death rates for the total population.
        #
        # Note: when a baseline scenario is simulated the alignment factors are
        #   calculated and written to file. When an alternative scenario is 
        #   simulated the adjustment factors from the baseline scenario are read 
        #   from file and used. Using pre-calculated alignment factors in 
        #   alternative scenaries prevents the alignment procedure from 
        #   cancelling the effects of health interventions on mortality.
        if (is_baseline == TRUE) {
          

          # No adjustment if the population stock is empty (happens at age zero)
          if (yvec["s_pop"] > 0) {
            
            # Nr of deaths according to Statistics Sweden
            w_deaths <- r_dead_all * (yvec["s_pop"] + yvec["s_cvd"] +
                                        yvec["s_cancer"] + yvec["s_comorb"])
            
            # Simulated (unadjusted) nr of deaths
            denom <- (r_dead * yvec["s_pop"] + r_dead_cvd * yvec["s_cvd"] + 
                        r_dead_cancer * yvec["s_cancer"] + r_dead_comorb * yvec["s_comorb"])
            if (denom == 0) {
              writeLines(paste0("Error: zero denominator in proportional",
                                "adjustment of death rates for sex: ", sex_,
                                " and age: ", age_))
            }
            
            
            # Alpha: (proportional) adjustment factor for death rates
            alpha_ <- r_dead_all * (yvec["s_pop"] + yvec["s_cvd"] +
                                      yvec["s_cancer"] + yvec["s_comorb"]) /
              (r_dead * yvec["s_pop"] + r_dead_cvd * yvec["s_cvd"] +
                 r_dead_cancer * yvec["s_cancer"] + r_dead_comorb * yvec["s_comorb"])
            
            # Adjustment (proportional) of death rates in each stock            
            r_dead <- r_dead * alpha_
            r_dead_cancer <- r_dead_cancer * alpha_
            r_dead_cvd <- r_dead_cvd * alpha_
            r_dead_comorb <- r_dead_comorb * alpha_
            
            # If negative death rate for pop reduce death rate for NCD
            # proportionally so that death rate for pop becomes zero
            if (r_dead < 0) {
              writeLines(paste0("Negative deaths in sex: ", sex_, ", age: ",
                                age_, ", r_dead: ", format(round(r_dead, 6),
                                                           scientific = FALSE), 6))
              alpha_ <- w_deaths / (yvec["s_cvd"] * r_dead_cvd +
                                      yvec["s_cancer"] * r_dead_cancer +  r_dead_comorb * yvec["s_comorb"])
              r_dead_cvd <- r_dead_cvd * alpha_
              r_dead_cancer <- r_dead_cancer * alpha_
              r_dead_comorb <- r_dead_comorb * alpha_
              r_dead <- (w_deaths - yvec["s_cvd"] *
                           r_dead_cvd - yvec["s_cancer"] * r_dead_cancer - r_dead_comorb * yvec["s_comorb"]) /
                yvec["s_pop"]
            }
            
            # store adjustment factors in calibration_mortality_table
            calibration_mortality_table[year == FINISH & sex == sex_ &
                                          age == age_,
                                        adj := alpha_]
            
          }
          else {
            alpha_ <- 1
          }
          
        } else {
          
          # If not simulating a baseline scenario: use mortality adjustment 
          # factors from base scenario.
          alpha_ <- calibration_mortality_table[year == FINISH & sex == sex_ &
                                                  age == age_,
                                                adj]
          
          if (yvec["s_pop"] > 0) {
            r_dead <- r_dead * alpha_
            r_dead_cancer <- r_dead_cancer * alpha_
            r_dead_cvd <- r_dead_cvd * alpha_
            r_dead_comorb <- r_dead_comorb * alpha_
            
            # If negative death rate for pop reduce death rate for NCD
            # proportionally so that death rate for pop becomes zero
            if (r_dead < 0) {
              writeLines(paste0("Negative deaths in sex: ", sex_, ", age: ",
                                age_, ", r_dead: ", format(round(r_dead, 6),
                                                           scientific = F), 6))
              alpha_ <- w_deaths / (yvec["s_cvd"] * r_dead_cvd +
                                      yvec["s_cancer"] * r_dead_cancer + r_dead_comorb * yvec["s_comorb"])
              r_dead_cvd <- r_dead_cvd * alpha_
              r_dead_cancer <- r_dead_cancer * alpha_
              r_dead_comorb <- r_dead_comorb * alpha_
              r_dead <- (w_deaths - yvec["s_cvd"] *
                           r_dead_cvd - yvec["s_cancer"] * r_dead_cancer - r_dead_comorb * yvec["s_comorb"]) /
                yvec["s_pop"]
            }
          }
        }

        # Note(TE221209): since the aggregate n_pop denotes the total population
        #    of then previous year the simulated number of births will tend to
        #    be somewhat underestimated in a growing population. This can be
        #    solved by estimating fertility rates for women, by age, to be
        #    included in avec and passed to ode().
        if (age_ == 0) {
          n_born <- dat_par$rate_born[year == FINISH & sex == sex_, rate] *
            n_pop 
        } else {
          n_born <- 0
        }
        
        # Prevalence of risk groups
        prev_smoking <- dat_par$prev_smoking[year == FINISH & sex == sex_ &
                                               age == age_, prev]
        prev_alcohol <- dat_par$prev_alcohol[year == FINISH & sex == sex_ &
                                               age == age_, prev]
        prev_inactivity <- dat_par$prev_inactivity[year == FINISH &
                                                     sex == sex_ & age == age_,
                                                   prev]
        prev_bmi <- dat_par$prev_bmi[year == FINISH & sex == sex_ &
                                               age == age_, prev]
        
        # Calculate cPAF for diet
        prev_diet <- dat_par$prev_diet[sex==sex_ & age == age_ ,..food_cat]
        
        
        prev_smoking <- pmin(prev_smoking *
                               piecewise(x=FINISH, x0=cfact_startyear,
                                         y0=1, x1=cfact_endyear,
                                         y1=unname(cfact["cfact_smoking"])), 1)
        prev_alcohol <- pmin(prev_alcohol *
                               piecewise(x=FINISH, x0=cfact_startyear,
                                         y0=1, x1=cfact_endyear,
                                         y1=unname(cfact["cfact_alcohol"])), 1)
        prev_inactivity <- pmin(prev_inactivity *
                                  piecewise(x=FINISH, x0=cfact_startyear,
                                            y0=1, x1=cfact_endyear,
                                            y1=unname(cfact["cfact_inactivity"])),
                                1)
        prev_bmi <- pmin(prev_bmi *
                               piecewise(x=FINISH, x0=cfact_startyear, y0=1,
                                         x1=cfact_endyear,
                                         y1=unname(cfact["cfact_bmi"])), 1)
        
        cfact_food_0 <- c(fruit  = 1.0,
                          wholegrains = 1.0,
                          greens = 1.0,
                          meat = 1.0, 
                          salt = 1.0)
        
        
        cfact_food_use <- mapply(piecewise, x=FINISH, x0=cfact_startyear,
                                 y0=cfact_food_0, x1=cfact_endyear,
                                 y1=cfact_food)
        prev_diet_use <- pmin(prev_diet[, Map("*", .SD, cfact_food_use)], 1)
        
        
        
        pafs_diet_cancer <- rbindlist(apply(prev_diet_use, 1, paf,
                                            rr=rr_cancer_diet))
        cpaf_diet_cancer <- (1 - apply(1 - pafs_diet_cancer, 1, prod)) 
        pafs_diet_cvd <- rbindlist(apply(prev_diet_use, 1, paf, rr=rr_cvd_diet))
        cpaf_diet_cvd <- (1 - apply(1 - pafs_diet_cvd, 1, prod)) 
        cpaf_diet_comorb <- max(cpaf_diet_cancer,cpaf_diet_cvd)
        
        avec <- c(
          r_dead = r_dead, r_immig = r_immig, r_emig = r_emig, 
          n_born = n_born, 
          
          r_dead_cancer = r_dead_cancer,
          r_dead_cvd = r_dead_cvd,
          r_dead_comorb = r_dead_comorb,
          
          heal_cvd_pop = heal_cvd_pop, 
          heal_cancer_pop = heal_cancer_pop,
          heal_comorb_pop = heal_comorb_pop,
          heal_comorb_cancer = heal_comorb_cancer,
          heal_comorb_cvd = heal_comorb_cvd,
          
          rr_cancer_smoking = unname(rr["rr_cancer_smoking"]),
          rr_cancer_inactivity = unname(rr["rr_cancer_inactivity"]),
          rr_cancer_bmi = unname(rr["rr_cancer_bmi"]),
          rr_cancer_alcohol = unname(rr["rr_cancer_alcohol"]),
          rr_cvd_smoking = unname(rr["rr_cvd_smoking"]),
          rr_cvd_inactivity = unname(rr["rr_cvd_inactivity"]),
          rr_cvd_bmi = unname(rr["rr_cvd_bmi"]),
          rr_cvd_alcohol = unname(rr["rr_cvd_alcohol"]),
          comm_smoking = unname(communalities["smoking"]),
          comm_alcohol = unname(communalities["alcohol"]),
          comm_bmi = unname(communalities["bmi"]),
          comm_inactivity = unname(communalities["inactivity"]),
          prev_smoking = prev_smoking,
          prev_alcohol = prev_alcohol,
          prev_inactivity = prev_inactivity,
          prev_bmi = prev_bmi,
          
          calibration_cancer = calibration_cancer,
          calibration_cvd = calibration_cvd,
          calibration_comorb = calibration_comorb,
          calibration_cancer_comorb = calibration_cancer_comorb,
          calibration_cvd_comorb = calibration_cvd_comorb,
          cpaf_diet_cancer = cpaf_diet_cancer,
          cpaf_diet_cvd = cpaf_diet_cvd,
          paf_cvd_other = paf_cvd_other,
          paf_cancer_other = paf_cancer_other,
          paf_comorb_other = paf_comorb_other,
          paf_cancer_comorb_other = paf_cancer_comorb_other,
          paf_cvd_comorb_other = paf_cvd_comorb_other
          
        )
        
        # Call to solve ODE

        
        o <- ode(y = yvec,
                 times = simtime, 
                 func = model, 
                 parms = avec,
                 method = "euler")
        
        if (debug_ode) {
          # Merge columns of o to debug data
          thecols <- setdiff(colnames(o), c("year", "sex", "age"))
          dbo <- data.table(o)
          dbo[, ":="(
            year = FINISH,
            sex = sex_,
            age = age_,
            cohort = FINISH - age_
          )]
          
          fwrite(dbo, file = paste0(PROJECTROOT, "/Output/debug_ode_",
                                    timestamp, ".csv"),
                 sep = ";", append = TRUE)
        }
        
        # Write stocks and flows to dat. For stocks the last observation contains
        # the updated stocks. For flows the average flow of the within-year time
        # steps corresponds to the flow consistent with the initial and final
        # stocks.
        # Note: the demographic flows have been estimated as the number of events
        #  in relation to the mid-year population
        tmp <- setDT(data.frame(o))
        lastrow <- dim(tmp)[1]
        
        # Extract costs change rate values and calculate KPP for direct
        # and indirect costs
        dcost_cancer_index <- dat[year == FINISH & sex == sex_ & age == age_,
                                  dcost_growth_cancer]
        dcost_cvd_index <- dat[year == FINISH & sex == sex_ & age == age_,
                               dcost_growth_cvd]
        dcost_unit_index_cancer <- dcost_unit[sex == sex_ & age == age_ &
                                                ncd=="cancer", dcost_unit]
        icost_unit_index_cancer <- icost_unit_[sex == sex_ & age == age_ &
                                                 ncd=="cancer",
                                               icost_unit]
        dcost_unit_index_cvd <- dcost_unit[sex == sex_ & age == age_ &
                                             ncd=="cvd", dcost_unit]
        icost_unit_index_cvd <- icost_unit_[sex == sex_ & age == age_ &
                                              ncd=="cvd",
                                            icost_unit]
        
        # Add new data to dat
        dat[year == FINISH & sex == sex_ & age == age_, ":="(
          # Last value of stocks
          s_pop = tmp[lastrow]$s_pop,
          s_cancer = tmp[lastrow]$s_cancer,
          s_cvd = tmp[lastrow]$s_cvd,
          s_comorb = tmp[lastrow]$s_comorb,
          s_dead = tmp[lastrow]$s_dead,
          s_dead_cancer = tmp[lastrow]$s_dead_cancer,
          s_dead_cvd = tmp[lastrow]$s_dead_cvd,
          s_dead_comorb = tmp[lastrow]$s_dead_comorb,
          # Average value of flows
          f_dead = mean(tmp$f_dead),
          f_born = mean(tmp$f_born),
          f_immig_pop = mean(tmp$f_immig_pop),
          f_immig_cvd = mean(tmp$f_immig_cvd),
          f_immig_cancer = mean(tmp$f_immig_cancer),
          f_immig_comorb = mean(tmp$f_immig_comorb),
          f_emig_pop = mean(tmp$f_emig_pop),
          f_emig_cvd = mean(tmp$f_emig_cvd),
          f_emig_cancer = mean(tmp$f_emig_cancer),
          f_emig_comorb = mean(tmp$f_emig_comorb),
          f_pop_cancer = mean(tmp$f_pop_cancer),
          f_cancer_pop = mean(tmp$f_cancer_pop),
          f_cancer_dead = mean(tmp$f_cancer_dead),
          f_pop_cvd = mean(tmp$f_pop_cvd),
          f_cvd_pop = mean(tmp$f_cvd_pop),
          f_cvd_dead = mean(tmp$f_cvd_dead),
          f_pop_comorb = mean(tmp$f_pop_comorb),
          f_comorb_pop = mean(tmp$f_comorb_pop),
          f_comorb_dead = mean(tmp$f_comorb_dead),
          f_cancer_comorb = mean(tmp$f_cancer_comorb),
          f_comorb_cancer = mean(tmp$f_comorb_cancer),
          f_cvd_comorb = mean(tmp$f_cvd_comorb),
          f_comorb_cvd = mean(tmp$f_comorb_cvd),
          # Assumptions (e.g. rates, prevalences, healing rates)
          heal_cvd_pop = tmp[lastrow]$heal_cvd_pop,
          heal_cancer_pop = tmp[lastrow]$heal_cancer_pop,
          heal_comorb_pop = tmp[lastrow]$heal_comorb_pop,
          heal_comorb_cancer = tmp[lastrow]$heal_comorb_cancer,
          heal_comorb_cvd = tmp[lastrow]$heal_comorb_cvd,
          dr = tmp[lastrow]$dr,
          dr_cancer = tmp[lastrow]$dr_cancer,
          dr_cvd = tmp[lastrow]$dr_cvd,
          dr_comorb = tmp[lastrow]$dr_comorb,
          prev_alcohol = tmp[lastrow]$prev_alcohol,
          prev_bmi = tmp[lastrow]$prev_bmi,
          prev_smoking = tmp[lastrow]$prev_smoking,
          prev_inactivity = tmp[lastrow]$prev_inactivity,
          
          paf_cvd_smoking = fifelse(age_ >= age_cutoff_cvd & age_ <= age_cutoff_cvd_high,
                                    tmp[lastrow]$paf_cvd_smoking, 0.0),
          paf_cvd_inactivity = fifelse(age_ >= age_cutoff_cvd & age_ <= age_cutoff_cvd_high,
                                       tmp[lastrow]$paf_cvd_inactivity, 0.0),
          paf_cvd_bmi = fifelse(age_ >= age_cutoff_cvd & age_ <= age_cutoff_cvd_high,
                                    tmp[lastrow]$paf_cvd_bmi, 0.0),
          paf_cvd_alcohol = fifelse(age_ >= age_cutoff_cvd & age_ <= age_cutoff_cvd_high,
                                    tmp[lastrow]$paf_cvd_alcohol, 0.0),
          cpaf_diet_cvd = fifelse(age_ >= age_cutoff_cvd & age_ <= age_cutoff_cvd_high,
                                  tmp[lastrow]$cpaf_diet_cvd, 0.0),
          p_pop_cvd = tmp[lastrow]$p_pop_cvd,
          
          paf_cancer_smoking = fifelse(age_ >= age_cutoff_cancer & age_ <= age_cutoff_cancer_high,
                                       tmp[lastrow]$paf_cancer_smoking, 0.0),
          paf_cancer_inactivity = fifelse(age_ >= age_cutoff_cancer & age_ <= age_cutoff_cancer_high,
                                          tmp[lastrow]$paf_cancer_inactivity,
                                          0.0),
          paf_cancer_bmi = fifelse(age_ >= age_cutoff_cancer & age_ <= age_cutoff_cancer_high,
                                       tmp[lastrow]$paf_cancer_bmi, 0.0),
          paf_cancer_alcohol = fifelse(age_ >= age_cutoff_cancer & age_ <= age_cutoff_cancer_high,
                                       tmp[lastrow]$paf_cancer_alcohol, 0.0),
          cpaf_diet_cancer = fifelse(age_ >= age_cutoff_cancer & age_ <= age_cutoff_cancer_high,
                                     tmp[lastrow]$cpaf_diet_cancer, 0.0),
          p_pop_cancer = tmp[lastrow]$p_pop_cancer,
          paf_cancer_other = tmp[lastrow]$paf_cancer_other,
          paf_cvd_other = tmp[lastrow]$paf_cvd_other,
          paf_comorb_other = tmp[lastrow]$paf_comorb_other,
          paf_cancer_comorb_other = tmp[lastrow]$paf_cancer_comorb_other,
          paf_cvd_comorb_other = tmp[lastrow]$paf_cvd_comorb_other,
          calibration_cancer = tmp[lastrow]$calibration_cancer,
          calibration_cvd = tmp[lastrow]$calibration_cvd,
          calibration_comorb = tmp[lastrow]$calibration_comorb,
          calibration_cancer_comorb = tmp[lastrow]$calibration_cancer_comorb,
          calibration_cvd_comorb = tmp[lastrow]$calibration_cvd_comorb,
          calibration_mortality = alpha_,
          
          dcost_unit_cancer = dcost_unit_index_cancer,
          icost_unit_cancer = fifelse(is.na(icost_unit_index_cancer), 0.0,
                                      icost_unit_index_cancer),
          dcost_unit_cvd = dcost_unit_index_cvd,
          icost_unit_cvd = fifelse(is.na(icost_unit_index_cvd), 0.0,
                                   icost_unit_index_cvd),
          dcost_cancer = health_costs(dcost_total_base_year, 
                                      FINISH,dcost_cancer_index, 
                                      dcost_unit_index_cancer,
                                      tmp[lastrow]$s_cancer + tmp[lastrow]$s_comorb),
          dcost_cvd = health_costs(dcost_total_base_year, FINISH,
                                   dcost_cvd_index, dcost_unit_index_cvd,
                                   tmp[lastrow]$s_cvd + tmp[lastrow]$s_comorb),
          icost_cancer =  fifelse(is.na(icost_unit_index_cancer *
                                          (tmp[lastrow]$s_cancer + tmp[lastrow]$s_comorb)), 0.0,
                                  icost_unit_index_cancer * (tmp[lastrow]$s_cancer + tmp[lastrow]$s_comorb)),
          icost_cvd =  fifelse(is.na(icost_unit_index_cvd * (tmp[lastrow]$s_cvd + tmp[lastrow]$s_comorb)),
                               0.0, icost_unit_index_cvd * (tmp[lastrow]$s_cvd + tmp[lastrow]$s_comorb))
        )]
        # Write cPAF for diet with the right cutoff
        p_cancer <- ifelse(rep(age_ >= age_cutoff_cancer & age_ <= age_cutoff_cancer_high,
                               length(pafs_diet_cancer)),
                           pafs_diet_cancer, pafs_diet_cancer * 0)
        p_cvd <- ifelse(rep(age_ >= age_cutoff_cvd & age_ <= age_cutoff_cvd_high, length(pafs_diet_cvd)),
                        pafs_diet_cvd, pafs_diet_cvd * 0)
        dat[year == FINISH & sex == sex_ & age == age_,
            paste("prev", names(prev_diet_use), sep = "_") :=
              unname(prev_diet_use)]
        dat[year == FINISH & sex == sex_ & age == age_,
            paste("paf_cancer", names(pafs_diet_cancer), sep = "_") :=
              unname(p_cancer)]
        dat[year == FINISH & sex == sex_ & age == age_, 
            paste("paf_cvd", names(pafs_diet_cvd), sep = "_") := unname(p_cvd)]
        
        
      } # End loop age
    } # End loop sex
    
    ## Save data to file and write status and log if run by Shiny
    
    if (ui == TRUE) { 
      
      cat(switch(ui_lang,
                 "eng" = "Writing data. ", 
                 "swe" = "Skriver data. "
      ),
      file = paste0(scen_path, "log.txt"), 
      append = TRUE)
      
      ## Remove data file from previous year if it exists
      if (file.exists(paste0(scen_path, "data/dat_", starttime - 1, ".rda"))) {
        file.remove(paste0(scen_path, "data/dat_", starttime - 1, ".rda"))
      }
      
      ## Write latest data
      save(dat, file = paste0(scen_path, "data/dat_", starttime, ".rda"))
      
      cat(starttime, file = paste0(scen_path, "last_simyear.txt"))
      
      simyear_timer <- Sys.time() - t0
      cat(switch(ui_lang,
                 "eng" = paste0("Completed in ", round(simyear_timer, 0), " seconds . \n"), 
                 "swe" = paste0("Klart efter ", round(simyear_timer, 0), " sekunder . \n") 
      ),
      file = paste0(scen_path, "log.txt"), append = TRUE)
    }
    
  } # End loop year
  
  
  if (write_data_to_file == TRUE) {
    fwrite(dat, file = paste0(PROJECTROOT, "/Output/output_NCDSim_",
                              timestamp, ".csv"),
           sep = ";")
  }
  if (is_baseline == TRUE) {
    fwrite(dat, file = paste0(strsplit(baseline_parameters_path, ".json"),
                              ".csv"),
           sep = ";")
  }
  
  ## Write to log if run by Shiny
  if (ui == TRUE) {
    cat(switch(ui_lang,
               "eng" = "Simulation completed.",
               "eng" = "Simulering klar."
    ), 
    file = paste0(scen_path, "log.txt"), 
    append = TRUE)
  }
  
  return(dat)
}

  