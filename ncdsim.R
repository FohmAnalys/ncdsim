### Loading packages
library(data.table)
library(deSolve)
library(jsonlite)

source("calibrate.R")

##*****************************************************************************
## Miscellaneous helper functions
##*****************************************************************************

# Return indexed costs based on a KPP and a costs change rate for health care
health_costs <- function(index_year, endyear, costs_change_rate, kpp, stock) {
  n_years <- endyear - index_year
  if (n_years == 0) {
    return(stock * kpp)
  } else {
    return(
      stock * kpp * (1 + sign(n_years) * costs_change_rate)**(abs(n_years))
    )
  }
}


# piecewise function in 3 segments. Used to scale cFact
# Breakpoints at x0 and x1. Horizontal segment (y0) for x<= x0 and y1 for x > X1.
piecewise <- function(x, x0, y0, x1, y1) {
  # Use stats::approxfun for interpolation
  f <- approxfun(x = c(x0, x1), c(y0, y1), yleft = y0, yright = y1, rule = 2)
  return(f(x))
}

# Population Attributable fraction given a relative risk rr and a prevalence p
paf <- function(p, rr) {
  x <- p * (rr - 1)
  return(x / (x + 1))
}

# Cumulative Population Attributable Fraction given an array of PAFs
cpaf <- function(pafs) {
  return(1 - apply(1 - pafs, 1, prod))
}

# Effective paf when a cfact is applied directly to the PAF instead
# than to the prevalence (t.ex. for UV exposure)
eff_paf <- function(paf, c) {
  return(pmin((c * paf) / (paf * (c - 1) + 1),1.0))
}

# read_data_parameters returns a list of data.tables containing the data and
# parameters of the model
read_data_parameters <- function(filepath = "") {
  # Mortality rates
  par_rate_dead <- fread(file = paste0(filepath, "r_dead.csv"), sep = ";")

  # Mortality rates
  par_rate_dead_all <- fread(
    file = paste0(filepath, "r_dead_all.csv"),
    sep = ";"
  )

  # Immigration rates
  par_rate_immig <- fread(
    file = paste0(filepath, "r_immigration.csv"),
    sep = ";"
  )

  # Emigration rates
  par_rate_emig <- fread(file = paste0(filepath, "r_emigration.csv"), sep = ";")

  # Birth rates
  par_rate_born <- fread(file = paste0(filepath, "r_born.csv"), sep = ";")

  excess_mortality <- fread(
    file = paste0(filepath, "excess_mortality_disease.csv"),
    sep = ";"
  )

  get_prevalence <- function(filename) {
    tmp <- fread(file = paste0(filepath, filename), sep = ";")
    cnames <- setdiff(colnames(tmp), c("sex", "age"))
    tmp[, (cnames) := lapply(.SD, as.double), .SDcols = cnames]
    tmp2 <- melt(tmp, id.vars = c("sex", "age"))
    tmp2[, ":="(
      year = as.integer(as.character(variable)),
      prev = as.numeric(value)
    )]
    return(tmp2[, .(year, sex, age, prev)])
  }

  ## Risk factor prevalence: smoking
  par_prev_smoking <- get_prevalence("prev_smoking.csv")

  ## Risk factor prevalence: alcohol
  par_prev_alcohol <- get_prevalence("prev_alcohol.csv")

  ## Risk factor prevalence: physical inactivity
  par_prev_inactivity <- get_prevalence("prev_inactivity.csv")

  ## Risk factor prevalence: BMI
  par_prev_bmi <- get_prevalence("prev_bmi.csv")

  # Population counts (2000 - 2120), observed & predicted
  dat_popcounts <- fread(
    file = paste0(filepath, "pop_counts_scb.csv"),
    sep = ";"
  )

  # Read prevalences and relative risks related to diet
  diet_prev <- fread(file = paste0(filepath, "prev_diet.csv"), sep = ";")

  # Prevalences (for validation and initialization of stock)
  preval_new <- fread(paste0(filepath, "stocks.csv"), sep = ';')
  preval_new <- dcast(preval_new, year + sex + age ~ stock, value.var = 'N')
  colnames(preval_new) <- c(
    "year",
    "sex",
    "age",
    "prev_cancer",
    "prev_comorb",
    "prev_cvd",
    "totpop"
  )
  preval <- CJ(
    year = min(preval_new[, (year)]):max(preval_new[, (year)]),
    sex = 1:2,
    age = 0:100
  )
  preval <- merge(
    preval,
    preval_new[, c(
      "year",
      "age",
      "sex",
      "prev_cancer",
      "prev_comorb",
      "prev_cvd"
    )],
    by = c("year", "age", "sex"),
    all = TRUE
  )

  preval <- melt(
    preval,
    id.vars = c("sex", "age", "year"),
    value.vars = c("prev_cancer", "prev_comorb", "prev_cvd")
  )

  colnames(preval) <- c("sex", "age", "year", "ncd", "prevalence")
  preval <- setnafill(preval, "const", fill = 0, cols = c("prevalence"))
  par_prevalence_cancer <- preval[
    ncd == "prev_cancer",
    .(year, sex, age, prevalence)
  ]
  par_prevalence_cvd <- preval[ncd == "prev_cvd", .(year, sex, age, prevalence)]
  par_prevalence_comorb <- preval[
    ncd == "prev_comorb",
    .(year, sex, age, prevalence)
  ]

  return(list(
    rate_dead = par_rate_dead,
    rate_dead_all = par_rate_dead_all,
    rate_immig = par_rate_immig,
    rate_emig = par_rate_emig,
    rate_born = par_rate_born,
    excess_mortality_cancer = excess_mortality[
      stock == 'cancer',
      c("age", "sex", "year", "rate")
    ],
    excess_mortality_cvd = excess_mortality[
      stock == 'cvd',
      c("age", "sex", "year", "rate")
    ],
    excess_mortality_comorb = excess_mortality[
      stock == 'comorb',
      c("age", "sex", "year", "rate")
    ],
    prev_smoking = par_prev_smoking,
    prev_alcohol = par_prev_alcohol,
    prev_inactivity = par_prev_inactivity,
    prev_bmi = par_prev_bmi,
    popcounts = dat_popcounts,
    prevalence_cancer = par_prevalence_cancer,
    prevalence_cvd = par_prevalence_cvd,
    prevalence_comorb = par_prevalence_comorb,
    prev_diet = diet_prev
  ))
}

# Function to calculate stock vectors
calculate_stock_vectors <- function(last_year_stocks, stock_cols) {
  this_year <- copy(last_year_stocks)[, age := age + 1]
  age_0_rows <- data.table(sex = unique(this_year$sex), age = 0)
  age_0_rows <- age_0_rows[, (stock_cols) := 0]
  this_year <- rbindlist(list(age_0_rows, this_year))
  this_year[age >= 100, age := 100]
  this_year <- this_year[,
    lapply(.SD, sum, na.rm = TRUE),
    by = .(sex, age),
    .SDcols = stock_cols
  ]
  setorder(this_year, sex, age)
}

model <- function(time, stocks, auxs) {
  # Extract stocks
  s_pop <- stocks[["s_pop"]]
  s_cancer <- stocks[["s_cancer"]]
  s_cvd <- stocks[["s_cvd"]]
  s_comorb <- stocks[["s_comorb"]]
  s_dead <- stocks[["s_dead"]]
  s_dead_cancer <- stocks[["s_dead_cancer"]]
  s_dead_cvd <- stocks[["s_dead_cvd"]]
  s_dead_comorb <- stocks[["s_dead_comorb"]]

  # Extract auxs (list) - use [[ to get numeric values
  r_dead <- auxs[["r_dead"]]
  n_born <- auxs[["n_born"]]
  r_immig <- auxs[["r_immig"]]
  r_emig <- auxs[["r_emig"]]

  prev_smoking <- auxs[["prev_smoking"]]
  prev_inactivity <- auxs[["prev_inactivity"]]
  prev_bmi <- auxs[["prev_bmi"]]
  prev_alcohol <- auxs[["prev_alcohol"]]

  rr_cancer_smoking <- auxs[["rr_cancer_smoking"]]
  rr_cancer_inactivity <- auxs[["rr_cancer_inactivity"]]
  rr_cancer_bmi <- auxs[["rr_cancer_bmi"]]
  rr_cancer_alcohol <- auxs[["rr_cancer_alcohol"]]

  comm_smoking <- auxs[["comm_smoking"]]
  comm_inactivity <- auxs[["comm_inactivity"]]
  comm_bmi <- auxs[["comm_bmi"]]
  comm_alcohol <- auxs[["comm_alcohol"]]

  cpaf_diet_cancer <- auxs[["cpaf_diet_cancer"]]
  paf_cancer_uv <- auxs[["paf_cancer_uv"]]
  calibration_cancer <- auxs[["calibration_cancer"]]
  paf_cancer_other <- auxs[["paf_cancer_other"]]
  heal_cancer_pop <- auxs[["heal_cancer_pop"]]
  r_dead_cancer <- auxs[["r_dead_cancer"]]

  rr_cvd_smoking <- auxs[["rr_cvd_smoking"]]
  rr_cvd_inactivity <- auxs[["rr_cvd_inactivity"]]
  rr_cvd_bmi <- auxs[["rr_cvd_bmi"]]
  rr_cvd_alcohol <- auxs[["rr_cvd_alcohol"]]

  cpaf_diet_cvd <- auxs[["cpaf_diet_cvd"]]
  paf_cvd_uv <- auxs[["paf_cvd_uv"]]
  calibration_cvd <- auxs[["calibration_cvd"]]
  paf_cvd_other <- auxs[["paf_cvd_other"]]
  heal_cvd_pop <- auxs[["heal_cvd_pop"]]
  r_dead_cvd <- auxs[["r_dead_cvd"]]

  calibration_comorb <- auxs[["calibration_comorb"]]
  paf_comorb_other <- auxs[["paf_comorb_other"]]
  heal_comorb_pop <- auxs[["heal_comorb_pop"]]
  calibration_cancer_comorb <- auxs[["calibration_cancer_comorb"]]
  paf_cancer_comorb_other <- auxs[["paf_cancer_comorb_other"]]
  calibration_cvd_comorb <- auxs[["calibration_cvd_comorb"]]
  paf_cvd_comorb_other <- auxs[["paf_cvd_comorb_other"]]
  heal_comorb_cancer <- auxs[["heal_comorb_cancer"]]
  heal_comorb_cvd <- auxs[["heal_comorb_cvd"]]
  r_dead_comorb <- auxs[["r_dead_comorb"]]

  # Helper function for PAF calculations
  calculate_paf <- function(
    paf_uv,
    prev_smoking,
    prev_inactivity,
    prev_bmi,
    prev_alcohol,
    rr_smoking,
    rr_inactivity,
    rr_bmi,
    rr_alcohol,
    comm_smoking,
    comm_inactivity,
    comm_bmi,
    comm_alcohol,
    cpaf_diet,
    calibration,
    paf_other
  ) {
    paf_smoking <- paf(prev_smoking, rr_smoking)
    paf_inactivity <- paf(prev_inactivity, rr_inactivity)
    paf_bmi <- paf(prev_bmi, rr_bmi)
    paf_alcohol <- paf(prev_alcohol, rr_alcohol)

    p_pop_0_nm <- 1 -
      (1 - comm_smoking * paf_smoking) *
        (1 - comm_inactivity * paf_inactivity) *
        (1 - comm_bmi * paf_bmi) *
        (1 - comm_alcohol * paf_alcohol) *
      (1 - paf_uv)

    p_pop_0 <- 1 - (1 - cpaf_diet) * (1 - p_pop_0_nm)
    p_pop <- p_pop_0 * calibration + paf_other

    return(list(
      paf_smoking = paf_smoking,
      paf_inactivity = paf_inactivity,
      paf_bmi = paf_bmi,
      paf_alcohol = paf_alcohol,
      paf_uv = paf_uv,
      p_pop_0 = p_pop_0,
      p_pop = p_pop
    ))
  }

  # Demographic flows
  dr <- pmin(1.0, r_dead)
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

  # Cancer flows
  cancer_paf <- calculate_paf(
    paf_cancer_uv,
    prev_smoking,
    prev_inactivity,
    prev_bmi,
    prev_alcohol,
    rr_cancer_smoking,
    rr_cancer_inactivity,
    rr_cancer_bmi,
    rr_cancer_alcohol,
    comm_smoking,
    comm_inactivity,
    comm_bmi,
    comm_alcohol,
    cpaf_diet_cancer,
    calibration_cancer,
    paf_cancer_other
  )
  # Extract PAFs for cancer
  paf_cancer_smoking <- cancer_paf$paf_smoking
  paf_cancer_inactivity <- cancer_paf$paf_inactivity
  paf_cancer_bmi <- cancer_paf$paf_bmi
  paf_cancer_alcohol <- cancer_paf$paf_alcohol
  paf_cancer_uv <- cancer_paf$paf_uv
  p_pop_cancer_0 <- cancer_paf$p_pop_0
  p_pop_cancer <- cancer_paf$p_pop
  f_pop_cancer <- s_pop * p_pop_cancer
  f_cancer_pop <- s_cancer * heal_cancer_pop
  dr_cancer <- pmin(1.0, r_dead_cancer)
  f_cancer_dead <- s_cancer * dr_cancer

  # CVD flows
  cvd_paf <- calculate_paf(
    paf_cvd_uv,
    prev_smoking,
    prev_inactivity,
    prev_bmi,
    prev_alcohol,
    rr_cvd_smoking,
    rr_cvd_inactivity,
    rr_cvd_bmi,
    rr_cvd_alcohol,
    comm_smoking,
    comm_inactivity,
    comm_bmi,
    comm_alcohol,
    cpaf_diet_cvd,
    calibration_cvd,
    paf_cvd_other
  )

  # Extract PAFs for CVD
  paf_cvd_smoking <- cvd_paf$paf_smoking
  paf_cvd_inactivity <- cvd_paf$paf_inactivity
  paf_cvd_bmi <- cvd_paf$paf_bmi
  paf_cvd_alcohol <- cvd_paf$paf_alcohol
  paf_cvd_uv <- cvd_paf$paf_uv
  p_pop_cvd_0 <- cvd_paf$p_pop_0
  p_pop_cvd <- cvd_paf$p_pop
  f_pop_cvd <- s_pop * p_pop_cvd
  f_cvd_pop <- s_cvd * heal_cvd_pop
  dr_cvd <- pmin(1.0, r_dead_cvd)
  f_cvd_dead <- s_cvd * dr_cvd

  # Comorbidity flows
  #p_pop_comorb_0 <- max(p_pop_cancer_0, p_pop_cvd_0)
  p_pop_comorb_0 <- min(p_pop_cancer_0, p_pop_cvd_0)
  p_pop_comorb <- p_pop_comorb_0 * calibration_comorb + paf_comorb_other
  f_pop_comorb <- s_pop * p_pop_comorb
  f_cancer_comorb <- s_cancer *
    (p_pop_cvd_0 * calibration_cancer_comorb + paf_cancer_comorb_other)
  f_cvd_comorb <- s_cvd *
    (p_pop_cancer_0 * calibration_cvd_comorb + paf_cvd_comorb_other)
  f_comorb_pop <- s_comorb * heal_comorb_pop
  f_comorb_cancer <- s_comorb * heal_comorb_cancer
  f_comorb_cvd <- s_comorb * heal_comorb_cvd
  dr_comorb <- pmin(1.0, r_dead_comorb)
  f_comorb_dead <- s_comorb * dr_comorb

  # Differential equations
  d_s_pop_dt <- f_born -
    f_dead +
    f_immig_pop -
    f_emig_pop +
    f_cvd_pop -
    f_pop_cvd +
    f_cancer_pop -
    f_pop_cancer +
    f_comorb_pop -
    f_pop_comorb
  d_s_cancer_dt <- f_pop_cancer -
    f_cancer_pop +
    f_immig_cancer -
    f_emig_cancer -
    f_cancer_dead +
    f_comorb_cancer -
    f_cancer_comorb
  d_s_cvd_dt <- f_pop_cvd -
    f_cvd_pop +
    f_immig_cvd -
    f_emig_cvd -
    f_cvd_dead +
    f_comorb_cvd -
    f_cvd_comorb
  d_s_comorb_dt <- f_pop_comorb -
    f_comorb_pop +
    f_immig_comorb -
    f_emig_comorb -
    f_comorb_dead +
    f_cancer_comorb -
    f_comorb_cancer +
    f_cvd_comorb -
    f_comorb_cvd
  d_s_dead_dt <- f_dead
  d_s_dead_cancer_dt <- f_cancer_dead
  d_s_dead_cvd_dt <- f_cvd_dead
  d_s_dead_comorb_dt <- f_comorb_dead

  # Return results
  return(
    list(
      # Differential equations (unnamed vector)
      c(
        d_s_pop_dt,
        d_s_cancer_dt,
        d_s_cvd_dt,
        d_s_comorb_dt,
        d_s_dead_dt,
        d_s_dead_cancer_dt,
        d_s_dead_cvd_dt,
        d_s_dead_comorb_dt
      ),
      # Named elements
      f_dead = f_dead,
      f_born = f_born,
      f_immig_pop = f_immig_pop,
      f_immig_cvd = f_immig_cvd,
      f_immig_cancer = f_immig_cancer,
      f_immig_comorb = f_immig_comorb,
      f_emig_pop = f_emig_pop,
      f_emig_cvd = f_emig_cvd,
      f_emig_cancer = f_emig_cancer,
      f_emig_comorb = f_emig_comorb,
      f_pop_cancer = f_pop_cancer,
      f_cancer_pop = f_cancer_pop,
      f_cancer_dead = f_cancer_dead,
      f_pop_cvd = f_pop_cvd,
      f_cvd_pop = f_cvd_pop,
      f_cvd_dead = f_cvd_dead,
      f_pop_comorb = f_pop_comorb,
      f_comorb_pop = f_comorb_pop,
      f_cvd_comorb = f_cvd_comorb,
      f_comorb_cvd = f_comorb_cvd,
      f_cancer_comorb = f_cancer_comorb,
      f_comorb_cancer = f_comorb_cancer,
      f_comorb_dead = f_comorb_dead,
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
      paf_cvd_uv = paf_cvd_uv,
      cpaf_diet_cvd = cpaf_diet_cvd,
      p_pop_cvd = p_pop_cvd,
      paf_cancer_smoking = paf_cancer_smoking,
      paf_cancer_inactivity = paf_cancer_inactivity,
      paf_cancer_bmi = paf_cancer_bmi,
      paf_cancer_alcohol = paf_cancer_alcohol,
      paf_cancer_uv = paf_cancer_uv,
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
    )
  )
}


#============== Helper functions ===============#
format_baseline_parameters <- function(baseline_parameters) {
  # Bundle all calibration-related parameters into a named list
  list(
    calibrate_cpaf_cancer = as.integer(
      baseline_parameters$calibrate_cpaf_cancer
    ),
    calibrate_cpaf_cvd = as.integer(baseline_parameters$calibrate_cpaf_cvd),
    cpaf_cancer = as.numeric(baseline_parameters$cpaf_cancer),
    cpaf_cvd = as.numeric(baseline_parameters$cpaf_cvd),
    age_cutoff_cancer = as.integer(baseline_parameters$age_cutoff_cancer),
    age_cutoff_cvd = as.integer(baseline_parameters$age_cutoff_cvd),
    age_cutoff_cancer_high = as.integer(
      baseline_parameters$age_cutoff_cancer_high
    ),
    age_cutoff_cvd_high = as.integer(baseline_parameters$age_cutoff_cvd_high),
    dcost_total_cancer = as.numeric(baseline_parameters$dcost_total_cancer),
    icost_total_cancer = as.numeric(baseline_parameters$icost_total_cancer),
    dcost_growth_cancer_input = as.numeric(
      baseline_parameters$dcost_growth_cancer
    ),
    dcost_total_cvd = as.numeric(baseline_parameters$dcost_total_cvd),
    icost_total_cvd = as.numeric(baseline_parameters$icost_total_cvd),
    dcost_growth_cvd_input = as.numeric(baseline_parameters$dcost_growth_cvd),
    dcost_total_base_year = as.integer(
      baseline_parameters$dcost_total_base_year
    ),
    rr = unlist(baseline_parameters$rr),
    rr_cancer_diet = unlist(baseline_parameters$rr_cancer_diet),
    rr_cvd_diet = unlist(baseline_parameters$rr_cvd_diet),
    communalities = unlist(baseline_parameters$communalities),
    direct_costs = c(
      cancer = baseline_parameters$dcost_total_cancer,
      cvd = baseline_parameters$dcost_total_cvd
    ),
    indirect_costs = c(
      cancer = baseline_parameters$icost_total_cancer,
      cvd = baseline_parameters$icost_total_cvd
    )
  )
}

pad_paf_other <- function(paf_other, startyear, endyear) {
  # If the startyear < first year in paf_other: pad years between with paf(first_year)
  # if last year of simulation is > last_year in paf_other: pad years between with paf(last_year)
  first_year <- min(as.integer(paf_other[, year]))
  last_year <- max(as.integer(paf_other[, year]))

  if (startyear < first_year) {
    tmp <- paf_other[year == first_year]
    years_to_add <- startyear:(first_year - 1)
    padded <- lapply(years_to_add, function(y) {
      tmp_copy <- copy(tmp)
      tmp_copy[, year := y]
      return(tmp_copy)
    })
    paf_other <- rbindlist(c(padded, list(paf_other)), fill = TRUE)
  }

  if (endyear > last_year) {
    tmp <- paf_other[year == last_year]
    years_to_add <- (last_year + 1):endyear
    padded <- lapply(years_to_add, function(y) {
      tmp_copy <- copy(tmp)
      tmp_copy[, year := y]
      return(tmp_copy)
    })
    paf_other <- rbindlist(c(list(paf_other), padded), fill = TRUE)
  }

  return(paf_other)
}

handle_calibration <- function(
  is_baseline,
  input_path,
  baseline_parameters_path,
  baseline_parameters,
  startyear,
  endyear
) {
  # Perform calibration on the go
  # Run calibration and assign variables. Read into calibrate.R for details
  if (is_baseline) {
    calibration_results <- calibrate(
      input_path = input_path,
      params = baseline_parameters,
      base_year = startyear
    )

    calibration_results$icost_unit_ <- calibration_results$icost_unit
    calibration_results$icost_unit <- NULL
    # pad paf_other if the simulation period is longer than the calculated values.
    first_year <- min(as.integer(calibration_results$paf_other[, year]))
    last_year <- max(as.integer(calibration_results$paf_other[, year]))
    calibration_results$paf_other_calibration_0 <- pad_paf_other(
      calibration_results$paf_other,
      startyear,
      endyear
    )
    calibration_results$paf_other <- NULL

    calibration_results$calibration_mortality_table <- CJ(
      year = startyear:endyear,
      sex = 1:2,
      age = 0:100,
      adj = 0.0
    )
    # Set keys to allow faster subsetting
    setkeyv(calibration_results$dcost_unit, c("sex", "age", "ncd"))
    setkeyv(calibration_results$icost_unit_, c("sex", "age", "ncd"))
    setkeyv(calibration_results$const, c("sex", "age", "disease"))
    calibration_results$paf_other_calibration_0 <- calibration_results$paf_other_calibration_0[,
      year := as.numeric(year)
    ]
    setkeyv(
      calibration_results$paf_other_calibration_0,
      c("sex", "age", "year", "disease")
    )
    setkeyv(calibration_results$heal_rate, c("sex", "age", "disease"))
    setkeyv(
      calibration_results$calibration_mortality_table,
      c("year", "sex", "age")
    )
    return(calibration_results)
  } else {
    # If not a baseline run, then read calibration parameters
    # from the correspondent baseline run and format accordingly
    cal_param_file <- paste0(
      strsplit(baseline_parameters_path, ".json"),
      ".csv"
    )
    calibrated_parameters <- fread(cal_param_file, sep = ";")
    calibration_mortality_table <- CJ(
      year = startyear:endyear,
      age = 0:100,
      sex = 1:2
    )
    calibration_mortality_table <- merge(
      calibration_mortality_table,
      calibrated_parameters[,
        .(year, age, sex, calibration_mortality)
      ],
      by = c("year", "age", "sex"),
      all = FALSE
    )
    colnames(calibration_mortality_table) <- c(
      "year",
      "age",
      "sex",
      "adj"
    )

    paf_other_calibration_0 <- calibrated_parameters[, c(
      "sex",
      "age",
      "year",
      "paf_cancer_other",
      "paf_cvd_other",
      "paf_cancer_comorb_other",
      "paf_cvd_comorb_other",
      "paf_comorb_other"
    )]
    paf_other_calibration_0 <- setnafill(
      paf_other_calibration_0,
      "const",
      fill = 0,
      cols = c(
        "paf_cancer_other",
        "paf_cvd_other",
        "paf_cancer_comorb_other",
        "paf_cvd_comorb_other",
        "paf_comorb_other"
      )
    )
    colnames(paf_other_calibration_0) <- c(
      "sex",
      "age",
      "year",
      "cpaf_cancer",
      "cpaf_cvd",
      "cpaf_cancer_comorb",
      "cpaf_cvd_comorb",
      "cpaf_comorb"
    )
    paf_other_calibration_0 <- melt(
      paf_other_calibration_0,
      id.vars = c("sex", "age", "year"),
      variable.name = 'disease',
      value.name = 'cpaf'
    )

    const <- calibrated_parameters[
      year == as.integer((startyear + endyear) / 2),
      c(
        "sex",
        "age",
        "calibration_cancer",
        "calibration_cvd",
        "calibration_cancer_comorb",
        "calibration_cvd_comorb",
        "calibration_comorb"
      )
    ]
    const <- setnafill(
      const,
      "const",
      fill = 0,
      cols = c(
        "calibration_cancer",
        "calibration_cvd",
        "calibration_cancer_comorb",
        "calibration_cvd_comorb",
        "calibration_comorb"
      )
    )
    colnames(const) <- c(
      "sex",
      "age",
      "const_cancer",
      "const_cvd",
      "const_cancer_comorb",
      "const_cvd_comorb",
      "const_comorb"
    )
    const <- melt(
      const,
      id.vars = c("sex", "age"),
      variable.name = 'disease',
      value.name = 'constant'
    )

    heal_rate <- calibrated_parameters[
      year == as.integer((startyear + endyear) / 2),
      c(
        "sex",
        "age",
        "heal_cvd_pop",
        "heal_cancer_pop",
        "heal_comorb_pop",
        "heal_comorb_cancer",
        "heal_comorb_cvd"
      )
    ]
    heal_rate <- setnafill(
      heal_rate,
      "const",
      fill = 0.0,
      cols = c(
        "heal_cvd_pop",
        "heal_cancer_pop",
        "heal_comorb_pop",
        "heal_comorb_cancer",
        "heal_comorb_cvd"
      )
    )
    colnames(heal_rate) <- c(
      "sex",
      "age",
      "cvd",
      "cancer",
      "comorb_healthy",
      "comorb_cancer",
      "comorb_cvd"
    )
    heal_rate <- melt(
      heal_rate,
      id.vars = c("sex", "age"),
      variable.name = 'disease',
      value.name = 'heal_rate'
    )

    dcost_unit <- CJ(sex = 1:2, age = 0:100)
    dcost_unit <- merge(
      dcost_unit,
      calibrated_parameters[
        year == as.integer((startyear + endyear) / 2),
        .(age, sex, dcost_unit_cancer, dcost_unit_cvd)
      ],
      by = c("age", "sex"),
      all = TRUE
    )
    colnames(dcost_unit) <- c("age", "sex", "cancer", "cvd")
    dcost_unit <- melt(dcost_unit, id.vars = c("age", "sex"))
    colnames(dcost_unit) <- c("age", "sex", "ncd", "dcost_unit")
    # _ to avoid error with casting to the data table
    icost_unit_ <- CJ(sex = 1:2, age = 0:100)
    icost_unit_ <- merge(
      icost_unit_,
      calibrated_parameters[
        year == as.integer((startyear + endyear) / 2),
        .(age, sex, icost_unit_cancer, icost_unit_cvd)
      ],
      by = c("age", "sex"),
      all = TRUE
    )
    colnames(icost_unit_) <- c("age", "sex", "cancer", "cvd")
    icost_unit_ <- melt(icost_unit_, id.vars = c("age", "sex"))
    colnames(icost_unit_) <- c("age", "sex", "ncd", "icost_unit")
  }
  # Calibration DONE
  # Set keys to allow faster subsetting
  setkeyv(dcost_unit, c("sex", "age", "ncd"))
  setkeyv(icost_unit_, c("sex", "age", "ncd"))
  setkeyv(const, c("sex", "age", "disease"))
  paf_other_calibration_0 <- paf_other_calibration_0[,
    year := as.numeric(year)
  ]
  setkeyv(
    paf_other_calibration_0,
    c("sex", "age", "year", "disease")
  )
  setkeyv(heal_rate, c("sex", "age", "disease"))
  setkeyv(
    calibration_mortality_table,
    c("year", "sex", "age")
  )
  return(
    list(
      constants = const,
      paf_other_calibration_0 = paf_other_calibration_0,
      heal_rate = heal_rate,
      dcost_unit = dcost_unit,
      icost_unit_ = icost_unit_,
      calibration_mortality_table = calibration_mortality_table
    )
  )
}

initialize_simulation_data <- function(
  startyear,
  endyear,
  dat_par,
  baseline_parameters,
  food_cat
) {
  dat <- CJ(year = startyear:endyear, sex = 1:2, age = 0:100, s_pop = 0.0)

  dat[
    dat_par$popcounts[year == startyear],
    on = .(year == year, sex == sex, age == age),
    s_pop_ := i.pop
  ]

  # Empty flow variables
  dat[,
    c(
      "f_dead",
      "f_born",
      "f_immig_pop",
      "f_immig_cvd",
      "f_immig_cancer",
      "f_immig_comorb",
      "f_pop_cancer",
      "f_cancer_pop",
      "f_pop_cvd",
      "f_cvd_pop",
      "f_pop_comorb",
      "f_comorb_pop",
      "f_cancer_comorb",
      "f_comorb_cancer",
      "f_cvd_comorb",
      "f_comorb_cvd",
      "f_cancer_dead",
      "f_cvd_dead",
      "f_comorb_dead"
    ) := as.numeric(NA)
  ]

  # Empty cost variables
  dat[,
    c(
      "dcost_unit_cancer",
      "icost_unit_cancer",
      "dcost_unit_cvd",
      "icost_unit_cvd",
      "dcost_cancer",
      "dcost_cvd",
      "icost_cancer",
      "icost_cvd",
      "dcost_growth_cancer",
      "dcost_growth_cvd"
    ) := as.numeric(NA)
  ]

  # Empty prevalence and PAF variables by food category
  dat[, (paste("prev", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_cancer", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_cvd", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_cancer_comorb", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_cvd_comorb", food_cat, sep = "_")) := as.numeric(NA)]
  dat[, (paste("paf_comorb", food_cat, sep = "_")) := as.numeric(NA)]

  # Set cost growth rates
  dat[, dcost_growth_cancer := baseline_parameters$dcost_growth_cancer_input]
  dat[, dcost_growth_cvd := baseline_parameters$dcost_growth_cvd_input]

  # Merge prevalences
  dat[
    dat_par$prevalence_cancer[year == startyear],
    on = .(year == year, sex == sex, age == age),
    s_cancer := as.numeric(i.prevalence)
  ]
  dat[
    dat_par$prevalence_cvd[year == startyear],
    on = .(year == year, sex == sex, age == age),
    s_cvd := as.numeric(i.prevalence)
  ]
  dat[
    dat_par$prevalence_comorb[year == startyear],
    on = .(year == year, sex == sex, age == age),
    s_comorb := as.numeric(i.prevalence)
  ]

  # Set keys and adjust stocks
  setkeyv(dat, c("year", "sex", "age"))
  dat[.(startyear), s_pop := s_pop_ - s_cancer - s_cvd - s_comorb]

  idx <- which(dat[, s_pop < 0])
  dat[idx, pct := 1 - abs(s_pop) / (s_cancer + s_cvd + s_comorb)]
  dat[
    idx,
    ":="(
      s_cancer = pct * s_cancer,
      s_cvd = pct * s_cvd,
      s_comorb = pct * s_comorb,
      s_pop = 0
    )
  ]
  dat[, c("s_pop_", "pct") := NULL]

  # Initialize dead stocks
  dat[, ":="(
    s_dead_cancer = 0,
    s_dead_cvd = 0,
    s_dead_comorb = 0,
    s_dead = 0
  )]

  # Empty emigration rates
  dat[,
    c(
      "f_emig_pop",
      "f_emig_cvd",
      "f_emig_cancer",
      "f_emig_comorb"
    ) := as.numeric(NA)
  ]

  # Empty heal rates
  dat[,
    c(
      "heal_cvd_pop",
      "heal_cancer_pop",
      "heal_comorb_pop",
      "heal_comorb_cancer",
      "heal_comorb_cvd"
    ) := as.numeric(NA)
  ]

  # Empty dr
  dat[, c("dr", "dr_cancer", "dr_cvd", "dr_comorb") := as.numeric(NA)]

  # Empty prevalences
  dat[,
    c(
      "prev_alcohol",
      "prev_bmi",
      "prev_smoking",
      "prev_inactivity"
    ) := as.numeric(NA)
  ]

  # Empty PAF variables
  dat[,
    c(
      "paf_cvd_smoking",
      "paf_cvd_inactivity",
      "paf_cvd_bmi",
      "paf_cvd_alcohol",
      "paf_cvd_uv",
      "cpaf_diet_cvd",
      "p_pop_cvd",
      "paf_cancer_smoking",
      "paf_cancer_inactivity",
      "paf_cancer_bmi",
      "paf_cancer_alcohol",
      "paf_cancer_uv",
      "cpaf_diet_cancer",
      "p_pop_cancer",
      "paf_cancer_other",
      "paf_cvd_other",
      "paf_comorb_other",
      "paf_cancer_comorb_other",
      "paf_cvd_comorb_other"
    ) := as.numeric(NA)
  ]

  # Empty calibration
  dat[,
    c(
      "calibration_cancer",
      "calibration_cvd",
      "calibration_comorb",
      "calibration_cancer_comorb",
      "calibration_cvd_comorb",
      "calibration_mortality"
    ) := as.numeric(NA)
  ]

  return(dat)
}

setup_year_context <- function(
  starttime,
  dat,
  dat_par,
  nsteps,
  max_year_mortality
) {
  START <- starttime - 1
  FINISH <- START + 1
  STEP <- 1 / nsteps
  simtime <- seq(START, FINISH, by = STEP)

  last_year <- dat[
    year == START,
    .(
      sex,
      age,
      s_pop,
      s_cancer,
      s_cvd,
      s_comorb,
      s_dead,
      s_dead_cancer,
      s_dead_cvd,
      s_dead_comorb
    )
  ]
  stock_cols <- setdiff(names(last_year), c("sex", "age"))
  stocks <- calculate_stock_vectors(last_year, stock_cols)
  setkeyv(stocks, c("sex", "age"))

  n_pop <- last_year[
    age > 19 & age < 40,
    sum(s_pop + s_cancer + s_cvd + s_comorb)
  ]

  target_year <- if (FINISH > max_year_mortality) max_year_mortality else FINISH
  death_and_mortality_rates <- CJ(sex = unique(dat$sex), age = unique(dat$age))

  death_and_mortality_rates[
    dat_par$rate_dead[year == target_year],
    on = .(sex, age),
    r_dead := i.rate
  ]
  death_and_mortality_rates[r_dead == 0, r_dead := 1e-10]

  # NOTE THAT rate_dead_all does not use
  # target_year as per original implementation
  death_and_mortality_rates[
    dat_par$rate_dead_all[year == FINISH],
    on = .(sex, age),
    r_dead_all := i.rate
  ]
  death_and_mortality_rates[r_dead_all == 0, r_dead_all := 1e-10]

  death_and_mortality_rates[
    dat_par$excess_mortality_cancer[year == target_year],
    on = .(sex, age),
    excess_mort_cancer := i.rate
  ]
  death_and_mortality_rates[, r_dead_cancer := r_dead * excess_mort_cancer]
  death_and_mortality_rates[r_dead_cancer > 0.999, r_dead_cancer := 0.999]

  death_and_mortality_rates[
    dat_par$excess_mortality_cvd[year == target_year],
    on = .(sex, age),
    excess_mort_cvd := i.rate
  ]
  death_and_mortality_rates[, r_dead_cvd := r_dead * excess_mort_cvd]
  death_and_mortality_rates[r_dead_cvd > 0.999, r_dead_cvd := 0.999]

  death_and_mortality_rates[
    dat_par$excess_mortality_comorb[year == target_year],
    on = .(sex, age),
    excess_mort_comorb := i.rate
  ]
  death_and_mortality_rates[, r_dead_comorb := r_dead * excess_mort_comorb]
  death_and_mortality_rates[r_dead_comorb > 0.999, r_dead_comorb := 0.999]

  return(list(
    START = START,
    FINISH = FINISH,
    STEP = STEP,
    simtime = simtime,
    stocks = stocks,
    stock_cols = stock_cols,
    n_pop = n_pop,
    death_and_mortality_rates = death_and_mortality_rates
  ))
}

create_yvec <- function(stocks, sex_, age_, stock_cols) {
  unlist(stocks[.(sex_, age_), .SD, .SDcols = stock_cols])
}

create_avec <- function(
  yvec,
  sex_,
  age_,
  FINISH,
  death_and_mortality_rates,
  rate_immig,
  rate_emig,
  rate_born,
  prev_smoking,
  prev_alcohol,
  prev_inactivity,
  prev_bmi,
  prev_diet,
  heal_rate,
  const,
  paf_other_calibration_0,
  global_params, # Precomputed
  n_pop,
  rr_cancer_diet,
  rr_cvd_diet,
  is_baseline,
  calibration_mortality_table,
  cfact_startyear,
  cfact_endyear,
  cfact,
  cfact_food_use,
  food_cat
) {
  # --- Death rates ---
  r_dead <- death_and_mortality_rates[.(sex_, age_), r_dead]
  r_dead_all <- death_and_mortality_rates[.(sex_, age_), r_dead_all]
  r_dead_cancer <- death_and_mortality_rates[.(sex_, age_), r_dead_cancer]
  r_dead_cvd <- death_and_mortality_rates[.(sex_, age_), r_dead_cvd]
  r_dead_comorb <- death_and_mortality_rates[.(sex_, age_), r_dead_comorb]

  # --- Migration rates ---
  r_immig <- rate_immig[
    .(FINISH, sex_, age_),
    rate
  ]
  r_emig <- rate_emig[.(FINISH, sex_, age_), rate]

  # --- Birth rate ---
  # Note(TE221209): since the aggregate n_pop denotes the total population
  #    of then previous year the simulated number of births will tend to
  #    be somewhat underestimated in a growing population. This can be
  #    solved by estimating fertility rates for women, by age, to be
  #    included in avec and passed to ode().
  n_born <- if (age_ == 0) {
    rate_born[.(FINISH, sex_), rate] * n_pop
  } else {
    0
  }

  # --- Healing rates ---
  heal_cvd_pop <- heal_rate[.(sex_, age_, "cvd"), heal_rate]
  heal_cancer_pop <- heal_rate[.(sex_, age_, "cancer"), heal_rate]
  heal_comorb_pop <- heal_rate[.(sex_, age_, "comorb_healthy"), heal_rate]
  heal_comorb_cancer <- heal_rate[.(sex_, age_, "comorb_cancer"), heal_rate]
  heal_comorb_cvd <- heal_rate[.(sex_, age_, "comorb_cvd"), heal_rate]

  # --- PAF other values ---
  paf_cvd_other <- paf_other_calibration_0[
    .(sex_, age_, FINISH, "cpaf_cvd"),
    cpaf
  ]
  paf_cancer_other <- paf_other_calibration_0[
    .(sex_, age_, FINISH, "cpaf_cancer"),
    cpaf
  ]
  paf_comorb_other <- paf_other_calibration_0[
    .(sex_, age_, FINISH, "cpaf_comorb"),
    cpaf
  ]
  paf_cancer_comorb_other <- paf_other_calibration_0[
    .(sex_, age_, FINISH, "cpaf_cancer_comorb"),
    cpaf
  ]
  paf_cvd_comorb_other <- paf_other_calibration_0[
    .(sex_, age_, FINISH, "cpaf_cvd_comorb"),
    cpaf
  ]

  # --- Calibration factors ---
  calibration_cancer <- const[.(sex_, age_, "const_cancer"), constant]
  calibration_cvd <- const[.(sex_, age_, "const_cvd"), constant]
  calibration_comorb <- const[.(sex_, age_, "const_comorb"), constant]
  calibration_cancer_comorb <- const[
    .(sex_, age_, "const_cancer_comorb"),
    constant
  ]
  calibration_cvd_comorb <- const[.(sex_, age_, "const_cvd_comorb"), constant]

  # --- Prevalences (with adjustment factors) ---
  prev_smoking <- pmin(
    prev_smoking[.(FINISH, sex_, age_), prev] *
      piecewise(
        FINISH,
        cfact_startyear,
        1,
        cfact_endyear,
        unname(cfact["cfact_smoking"])
      ),
    1
  )
  prev_alcohol <- pmin(
    prev_alcohol[.(FINISH, sex_, age_), prev] *
      piecewise(
        FINISH,
        cfact_startyear,
        1,
        cfact_endyear,
        unname(cfact["cfact_alcohol"])
      ),
    1
  )
  prev_inactivity <- pmin(
    prev_inactivity[.(FINISH, sex_, age_), prev] *
      piecewise(
        FINISH,
        cfact_startyear,
        1,
        cfact_endyear,
        unname(cfact["cfact_inactivity"])
      ),
    1
  )
  prev_bmi <- pmin(
    prev_bmi[.(FINISH, sex_, age_), prev] *
      piecewise(
        FINISH,
        cfact_startyear,
        1,
        cfact_endyear,
        unname(cfact["cfact_bmi"])
      ),
    1
  )
  
  paf_cancer_uv <- pmin(
    eff_paf(global_params[["paf_cancer_uv"]],
            piecewise(
              FINISH,
              cfact_startyear,
              1,
              cfact_endyear,
              unname(cfact["cfact_uv"]))
    ),
    1
  )
  
  paf_cvd_uv <- pmin(
    eff_paf(global_params[["paf_cvd_uv"]],
         piecewise(
           FINISH,
           cfact_startyear,
           1,
           cfact_endyear,
           unname(cfact["cfact_uv"]))
        ),
    1
  )

  # --- Diet prevalences and PAFs ---
  prev_diet <- prev_diet[sex == sex_ & age == age_, ..food_cat]
  prev_diet_use <- pmin(prev_diet[, Map("*", .SD, cfact_food_use)], 1)
  pafs_diet_cancer <- rbindlist(apply(
    prev_diet_use,
    1,
    paf,
    rr = rr_cancer_diet
  ))
  cpaf_diet_cancer <- cpaf(pafs_diet_cancer)

  pafs_diet_cvd <- rbindlist(apply(
    prev_diet_use,
    1,
    paf,
    rr = rr_cvd_diet
  ))
  cpaf_diet_cvd <- cpaf(pafs_diet_cvd)

  # --- Death rate alignment (uses yvec) ---
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
  alpha_ <- 1
  if (yvec[["s_pop"]] > 0) {
    if (is_baseline) {
      w_deaths <- r_dead_all *
        (yvec[["s_pop"]] +
          yvec[["s_cvd"]] +
          yvec[["s_cancer"]] +
          yvec[["s_comorb"]])
      denom <- (r_dead *
        yvec[["s_pop"]] +
        r_dead_cvd * yvec[["s_cvd"]] +
        r_dead_cancer * yvec[["s_cancer"]] +
        r_dead_comorb * yvec[["s_comorb"]])
      if (denom == 0) {
        stop(
          "Zero denominator in death rate adjustment for sex: ",
          sex_,
          ", age: ",
          age_
        )
      }
      alpha_ <- w_deaths / denom

      r_dead <- r_dead * alpha_
      r_dead_cancer <- r_dead_cancer * alpha_
      r_dead_cvd <- r_dead_cvd * alpha_
      r_dead_comorb <- r_dead_comorb * alpha_

      if (r_dead < 0) {
        alpha_ <- w_deaths /
          (yvec[["s_cvd"]] *
            r_dead_cvd +
            yvec[["s_cancer"]] * r_dead_cancer +
            r_dead_comorb * yvec[["s_comorb"]])
        r_dead_cvd <- r_dead_cvd * alpha_
        r_dead_cancer <- r_dead_cancer * alpha_
        r_dead_comorb <- r_dead_comorb * alpha_
        r_dead <- (w_deaths -
          yvec[["s_cvd"]] * r_dead_cvd -
          yvec[["s_cancer"]] * r_dead_cancer -
          r_dead_comorb * yvec[["s_comorb"]]) /
          yvec[["s_pop"]]
      }
      calibration_mortality_table[.(FINISH, sex_, age_), adj := alpha_]
    } else {
      alpha_ <- calibration_mortality_table[.(FINISH, sex_, age_), adj]
      if (yvec[["s_pop"]] > 0) {
        r_dead <- r_dead * alpha_
        r_dead_cancer <- r_dead_cancer * alpha_
        r_dead_cvd <- r_dead_cvd * alpha_
        r_dead_comorb <- r_dead_comorb * alpha_
        if (r_dead < 0) {
          alpha_ <- w_deaths /
            (yvec[["s_cvd"]] *
              r_dead_cvd +
              yvec[["s_cancer"]] * r_dead_cancer +
              r_dead_comorb * yvec[["s_comorb"]])
          r_dead_cvd <- r_dead_cvd * alpha_
          r_dead_cancer <- r_dead_cancer * alpha_
          r_dead_comorb <- r_dead_comorb * alpha_
          r_dead <- (w_deaths -
            yvec[["s_cvd"]] * r_dead_cvd -
            yvec[["s_cancer"]] * r_dead_cancer -
            r_dead_comorb * yvec[["s_comorb"]]) /
            yvec[["s_pop"]]
        }
      }
    }
  }

  # --- Combine all parameters ---
  c(
    r_dead = r_dead,
    r_immig = r_immig,
    r_emig = r_emig,
    n_born = n_born,
    r_dead_cancer = r_dead_cancer,
    r_dead_cvd = r_dead_cvd,
    r_dead_comorb = r_dead_comorb,
    heal_cvd_pop = heal_cvd_pop,
    heal_cancer_pop = heal_cancer_pop,
    heal_comorb_pop = heal_comorb_pop,
    heal_comorb_cancer = heal_comorb_cancer,
    heal_comorb_cvd = heal_comorb_cvd,
    prev_smoking = prev_smoking,
    prev_alcohol = prev_alcohol,
    prev_inactivity = prev_inactivity,
    prev_bmi = prev_bmi,
    prev_diet_use = prev_diet_use,
    cpaf_diet_cancer = cpaf_diet_cancer,
    cpaf_diet_cvd = cpaf_diet_cvd,
    paf_cancer_uv = paf_cancer_uv,
    paf_cvd_uv = paf_cvd_uv,
    paf_cvd_other = paf_cvd_other,
    paf_cancer_other = paf_cancer_other,
    paf_comorb_other = paf_comorb_other,
    paf_cancer_comorb_other = paf_cancer_comorb_other,
    paf_cvd_comorb_other = paf_cvd_comorb_other,
    pafs_diet_cancer = pafs_diet_cancer,
    pafs_diet_cvd = pafs_diet_cvd,
    calibration_mortality = alpha_,
    calibration_cancer = calibration_cancer,
    calibration_cvd = calibration_cvd,
    calibration_comorb = calibration_comorb,
    calibration_cancer_comorb = calibration_cancer_comorb,
    calibration_cvd_comorb = calibration_cvd_comorb,
    global_params
  )
}

compute_age_cutoff_pafs <- function(
  tmp_last,
  age_,
  baseline_parameters,
  disease
) {
  age_low <- baseline_parameters[[paste0("age_cutoff_", disease)]]
  age_high <- baseline_parameters[[paste0("age_cutoff_", disease, "_high")]]

  paf_cols <- grep(paste0("paf_", disease), names(tmp_last), value = TRUE)
  pafs <- as.list(tmp_last[1, ])[paf_cols]

  cpaf_col <- paste0("cpaf_diet_", disease)

  # Identify paf_disease_other and paf_disease_comorb_other columns
  other_paf_cols <- grep(
    paste0("paf_", disease, "_other$|paf_", disease, "_comorb_other$"),
    paf_cols,
    value = TRUE
  )
  non_other_paf_cols <- setdiff(names(pafs), other_paf_cols)

  # Set only non-other pafs to 0 outside age interval
  if (age_ < age_low | age_ > age_high) {
    pafs[non_other_paf_cols] <- 0.0
  }

  c(
    pafs,
    cpaf_diet = fifelse(
      age_ >= age_low & age_ <= age_high,
      tmp_last[[cpaf_col]],
      0.0
    )
  )
}

compute_costs <- function(
  tmp_last,
  FINISH,
  dcost_total_base_year,
  dcost_index,
  dcost_unit_index,
  icost_unit_index,
  stock_sum
) {
  dcost <- health_costs(
    dcost_total_base_year,
    FINISH,
    dcost_index,
    dcost_unit_index,
    stock_sum
  )
  icost <- fifelse(
    is.na(icost_unit_index * stock_sum),
    0.0,
    icost_unit_index * stock_sum
  )
  list(dcost = dcost, icost = icost)
}

create_update_row <- function(
  o, # ODE output (data.frame)
  sex_,
  age_,
  FINISH,
  avec,
  dat,
  calibration_results,
  baseline_parameters
) {
  tmp <- setDT(as.data.frame(o))
  tmp_last <- tmp[nrow(tmp), ] # Last row (final state)

  # --- Extract context-specific parameters ---
  dcost_cancer_index <- dat[.(FINISH, sex_, age_), dcost_growth_cancer]
  dcost_cvd_index <- dat[.(FINISH, sex_, age_), dcost_growth_cvd]
  dcost_unit_index_cancer <- calibration_results$dcost_unit[
    .(sex_, age_, "cancer"),
    dcost_unit
  ]
  icost_unit_index_cancer <- calibration_results$icost_unit_[
    .(sex_, age_, "cancer"),
    icost_unit
  ]
  dcost_unit_index_cvd <- calibration_results$dcost_unit[
    .(sex_, age_, "cvd"),
    dcost_unit
  ]
  icost_unit_index_cvd <- calibration_results$icost_unit_[
    .(sex_, age_, "cvd"),
    icost_unit
  ]

  # --- Extract from avec ---
  alpha_ <- avec[["calibration_mortality"]]

  # --- Stocks (last value) ---
  stocks <- tmp_last[,
    .SD,
    .SDcols = c(
      "s_pop",
      "s_cancer",
      "s_cvd",
      "s_comorb",
      "s_dead",
      "s_dead_cancer",
      "s_dead_cvd",
      "s_dead_comorb"
    )
  ]

  # --- Flows (mean) ---
  flows <- tmp[, lapply(.SD, mean), .SDcols = patterns("^f_")]

  # --- Assumptions (last value) ---
  assumptions <- tmp_last[,
    .SD,
    .SDcols = c(
      "heal_cvd_pop",
      "heal_cancer_pop",
      "heal_comorb_pop",
      "heal_comorb_cancer",
      "heal_comorb_cvd",
      "dr",
      "dr_cancer",
      "dr_cvd",
      "dr_comorb",
      "prev_alcohol",
      "prev_bmi",
      "prev_smoking",
      "prev_inactivity"
    )
  ]

  # --- PAFs (with age cutoffs) ---
  pafs_cvd <- compute_age_cutoff_pafs(
    tmp_last,
    age_,
    baseline_parameters,
    "cvd"
  )
  pafs_cvd$cpaf_diet_cvd <- pafs_cvd$cpaf_diet
  pafs_cvd$cpaf_diet <- NULL

  pafs_cancer <- compute_age_cutoff_pafs(
    tmp_last,
    age_,
    baseline_parameters,
    "cancer"
  )
  pafs_cancer$cpaf_diet_cancer <- pafs_cancer$cpaf_diet
  pafs_cancer$cpaf_diet <- NULL

  pafs <- c(
    pafs_cvd,
    pafs_cancer,
    cpaf_diet_cvd = tmp_last[1, cpaf_diet_cvd],
    p_pop_cvd = tmp_last[1, p_pop_cvd],
    p_pop_cancer = tmp_last[1, p_pop_cancer],
    paf_cancer_other = tmp_last[1, paf_cancer_other],
    paf_cvd_other = tmp_last[1, paf_cvd_other],
    paf_comorb_other = tmp_last[1, paf_comorb_other],
    paf_cancer_comorb_other = tmp_last[1, paf_cancer_comorb_other],
    paf_cvd_comorb_other = tmp_last[1, paf_cvd_comorb_other]
  )

  # --- Calibration ---
  calibration <- tmp_last[, .SD, .SDcols = patterns("^calibration_")]

  # --- Costs ---
  cancer_stock_sum <- tmp_last[1, s_cancer] + tmp_last[1, s_comorb]
  cvd_stock_sum <- tmp_last[1, s_cvd] + tmp_last[1, s_comorb]
  costs_cancer <- compute_costs(
    tmp_last,
    FINISH,
    baseline_parameters$dcost_total_base_year,
    dcost_cancer_index,
    dcost_unit_index_cancer,
    icost_unit_index_cancer,
    cancer_stock_sum
  )
  costs_cvd <- compute_costs(
    tmp_last,
    FINISH,
    baseline_parameters$dcost_total_base_year,
    dcost_cvd_index,
    dcost_unit_index_cvd,
    icost_unit_index_cvd,
    cvd_stock_sum
  )

  # Add diet-specific columns (with age cutoffs)
  prev_diet_cols <- grep("prev_diet", names(avec), value = TRUE)
  prev_diet_use <- avec[prev_diet_cols]
  pafs_diet_cancer_cols <- grep("pafs_diet_cancer", names(avec), value = TRUE)
  pafs_diet_cancer <- avec[pafs_diet_cancer_cols]
  pafs_diet_cvd_cols <- grep("pafs_diet_cvd", names(avec), value = TRUE)
  pafs_diet_cvd <- avec[pafs_diet_cvd_cols]

  # Apply age cutoffs
  p_cancer <- lapply(pafs_diet_cancer, function(x) {
    ifelse(
      age_ >= baseline_parameters$age_cutoff_cancer &
        age_ <= baseline_parameters$age_cutoff_cancer_high,
      x,
      0
    )
  })
  p_cvd <- lapply(pafs_diet_cvd, function(x) {
    ifelse(
      age_ >= baseline_parameters$age_cutoff_cvd &
        age_ <= baseline_parameters$age_cutoff_cvd_high,
      x,
      0
    )
  })

  # Add to result
  p_cancer <- setNames(
    as.list(p_cancer),
    sub("pafs_diet_cancer\\.", "paf_cancer_", names(p_cancer))
  )
  p_cvd <- setNames(
    as.list(p_cvd),
    sub("pafs_diet_cvd\\.", "paf_cvd_", names(p_cvd))
  )
  prev_diet_use <- setNames(
    as.list(prev_diet_use),
    sub("prev_diet_use\\.", "prev_", names(prev_diet_use))
  )

  # --- Combine all ---
  unlist(
    c(
      list(
        year = FINISH,
        sex = sex_,
        age = age_,
        calibration_mortality = alpha_
      ),
      stocks,
      flows,
      assumptions,
      pafs,
      calibration,
      list(
        dcost_cancer = costs_cancer$dcost,
        icost_cancer = costs_cancer$icost,
        dcost_cvd = costs_cvd$dcost,
        icost_cvd = costs_cvd$icost
      ),
      dcost_unit_cancer = dcost_unit_index_cancer,
      icost_unit_cancer = fifelse(
        is.na(icost_unit_index_cancer),
        0.0,
        icost_unit_index_cancer
      ),
      dcost_unit_cvd = dcost_unit_index_cvd,
      icost_unit_cvd = fifelse(
        is.na(icost_unit_index_cvd),
        0.0,
        icost_unit_index_cvd
      ),
      p_cancer,
      p_cvd,
      prev_diet_use
    ),
    recursive = FALSE
  )
}

log_debug_ode <- function(
  o_list,
  sex_age_combos,
  FINISH,
  simtime,
  PROJECTROOT
) {
  combined_output <- rbindlist(
    lapply(seq_along(o_list), function(i) {
      ode_output <- as.data.table(o_list[[i]])
      # Add sex and age columns, repeated for each row in the ode output
      ode_output[, sex := rep(sex_age_combos$sex[i], .N)]
      ode_output[, age := rep(sex_age_combos$age[i], .N)]
      ode_output[,
        cohort := rep(FINISH - sex_age_combos$age[i], .N)
      ]
      # Add the simulation time column
      ode_output[, time := simtime]
    }),
    fill = TRUE
  )

  # Reorder columns: time, sex, age, cohort, then the rest
  setcolorder(
    combined_output,
    neworder = c(
      "time",
      "sex",
      "age",
      "cohort",
      setdiff(names(combined_output), c("time", "sex", "age", "cohort"))
    )
  )
  fwrite(
    combined_output,
    file = paste0(
      PROJECTROOT,
      "/Output/debug_ode_",
      timestamp,
      ".csv"
    ),
    sep = ";",
    append = TRUE
  )
}

#================== UI related function ==============#
log_start <- function(scen_path, scen_name, ui_lang, startyear, endyear) {
  log_file <- paste0(scen_path, "log.txt")
  cat(
    switch(
      ui_lang,
      "eng" = paste0(
        "Starting simulation ",
        scen_name,
        " from ",
        startyear,
        " to ",
        endyear,
        ".\n"
      ),
      "swe" = paste0(
        "Startar simulering ",
        scen_name,
        " från ",
        startyear,
        " till ",
        endyear,
        ".\n"
      )
    ),
    file = log_file
  )
}
check_ui_status <- function(scen_path, ui_lang) {
  status_file <- paste0(scen_path, "ui_status.txt")
  ui_status <- tryCatch(
    scan(
      status_file,
      what = "character",
      quiet = TRUE
    ),
    error = function(e) NULL
  )

  if (length(ui_status) > 0 && ui_status == "stop") {
    cat(
      switch(
        ui_lang,
        "eng" = "Simulation stopped\n",
        "swe" = "Simulering stoppad\n"
      ),
      file = paste0(scen_path, "log.txt"),
      append = TRUE
    )
    stop("Stopped from ui")
  }
}

log_year_start <- function(scen_path, ui_lang, starttime) {
  ## Write current simulation year to file
  cat(starttime, file = paste0(scen_path, "current_simyear.txt"))

  ## Write to log
  cat(
    switch(
      ui_lang,
      "eng" = paste0("Year ", starttime, " started. "),
      "swe" = paste0("År ", starttime, " påbörjat. ")
    ),
    file = paste0(scen_path, "log.txt"),
    append = TRUE
  )
}
write_simulation_output <- function(ui_lang, scen_path, dat, starttime, t0) {
  # Log writing data message
  cat(
    switch(
      ui_lang,
      "eng" = "Writing data. ",
      "swe" = "Skriver data. "
    ),
    file = paste0(scen_path, "log.txt"),
    append = TRUE
  )

  # Remove previous year's data file if it exists
  prev_file <- paste0(scen_path, "data/dat_", starttime - 1, ".rda")
  if (file.exists(prev_file)) {
    file.remove(prev_file)
  }

  # Write latest data
  save(dat, file = paste0(scen_path, "data/dat_", starttime, ".rda"))

  # Update last_simyear.txt
  cat(starttime, file = paste0(scen_path, "last_simyear.txt"))

  # Log completion time
  simyear_timer <- Sys.time() - t0
  cat(
    switch(
      ui_lang,
      "eng" = paste0("Completed in ", round(simyear_timer, 0), " seconds . \n"),
      "swe" = paste0("Klart efter ", round(simyear_timer, 0), " sekunder . \n")
    ),
    file = paste0(scen_path, "log.txt"),
    append = TRUE
  )
}

write_simulation_results <- function(
  dat,
  write_data_to_file,
  is_baseline,
  ui,
  ui_lang,
  scen_path,
  PROJECTROOT,
  timestamp,
  baseline_parameters_path
) {
  if (write_data_to_file) {
    fwrite(
      dat,
      file = paste0(PROJECTROOT, "/Output/output_NCDSim_", timestamp, ".csv"),
      sep = ";"
    )
  }

  if (is_baseline) {
    baseline_csv_path <- paste0(
      strsplit(baseline_parameters_path, ".json")[[1]],
      ".csv"
    )
    fwrite(dat, file = baseline_csv_path, sep = ";")
  }

  if (ui) {
    cat(
      switch(
        ui_lang,
        "eng" = "Simulation completed.",
        "swe" = "Simulering klar."
      ),
      file = paste0(scen_path, "log.txt"),
      append = TRUE
    )
  }
}

#===================== Main function =================#
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
  baseline_parameters_path = paste0(
    PROJECTROOT,
    "/Input/defaults/baseline.json"
  ),
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
  cfact = c(
    cfact_smoking = 1.0,
    cfact_alcohol = 1.0,
    cfact_inactivity = 1.0,
    cfact_bmi = 1.0,
    cfact_uv = 1.0
  ),

  cfact_food = c(
    fruit = 1.0,
    wholegrains = 1.0,
    greens = 1.0,
    meat = 1.0,
    salt = 1.0
  ),

  # First year of applying interventions, i.e. cfact that is not equal to 1.0
  cfact_startyear = 2024,

  # Last year of onset for applying interventions, i.e.
  # cfact now is equal to cfact.
  # linear decrease
  cfact_endyear = 2025
) {
  # cat("In simulate_model...\n")
  # Assign timestamp to variable in .GlobalEnv to use in validation
  assign("NCDSim_timestamp", timestamp, envir = .GlobalEnv)

  ## Write to log when running from Shiny
  if (ui) {
    log_start(scen_path, scen_name, ui_lang, startyear, endyear)
  }
  ##*****************************************************************************
  ##* Create start data
  ##*****************************************************************************

  ## Path to /Input directory including files with risk factor prevalences etc
  # Read the input json file and assign parameters
  baseline_parameters <- read_json(baseline_parameters_path)
  baseline_parameters <- format_baseline_parameters(baseline_parameters)

  calibration_results <- handle_calibration(
    is_baseline,
    input_path,
    baseline_parameters_path,
    baseline_parameters,
    startyear,
    endyear
  )
  ## Get data and parameters
  dat_par <- read_data_parameters(input_path)
  max_year_mortality <- max(dat_par$excess_mortality_comorb$year)

  # extract basic values for diet
  food_cat <- colnames(dat_par$prev_diet)[
    !colnames(dat_par$prev_diet) %in% c("sex", "age")
  ]
  baseline_parameters$rr_cancer_diet <- data.frame(as.list(
    baseline_parameters$rr_cancer_diet
  ))
  baseline_parameters$rr_cvd_diet <- data.frame(as.list(
    baseline_parameters$rr_cvd_diet
  ))

  # Creating DT for storage of initial and simulated population
  # Note: dat is initialized to contain all data from the simulation. This
  #       prevents the model from having to concatenate data (and hence allocate
  #       new memory) within the model loop.
  dat <- initialize_simulation_data(
    startyear,
    endyear,
    dat_par,
    baseline_parameters,
    food_cat
  )

  cfact_food_0 <- c(
    fruit = 1.0,
    wholegrains = 1.0,
    greens = 1.0,
    meat = 1.0,
    salt = 1.0
  )
  ##****************************************************************************

  ## Demography only mode adjustments to relative risks and data
  if (only_demog) {
    baseline_parameters$rr[] <- 0 ## Set all relative risks to zero
    dat[.(startyear), s_pop := s_pop + s_cancer + s_cvd]
    dat[.(startyear), c("s_cancer", "s_cvd") := as.numeric(0)]
  }

  # Create minimal lookup tables from dat_par that
  # can be used to compute avec

  rate_emig_dt <- dat_par$rate_emig[
    year %between% c(startyear, endyear),
    .(rate),
    keyby = .(year, sex, age)
  ]
  rate_immig_dt <- dat_par$rate_immig[
    year %between% c(startyear, endyear),
    .(rate),
    keyby = .(year, sex, age)
  ]
  rate_born_dt <- dat_par$rate_born[
    year %between% c(startyear, endyear),
    .(rate),
    keyby = .(year, sex)
  ]
  prev_smoking_dt <- dat_par$prev_smoking[
    year %between% c(startyear, endyear),
    .(prev),
    keyby = .(year, sex, age)
  ]
  prev_alcohol_dt <- dat_par$prev_alcohol[
    year %between% c(startyear, endyear),
    .(prev),
    keyby = .(year, sex, age)
  ]
  prev_inactivity_dt <- dat_par$prev_inactivity[
    year %between% c(startyear, endyear),
    .(prev),
    keyby = .(year, sex, age)
  ]
  prev_bmi_dt <- dat_par$prev_bmi[
    year %between% c(startyear, endyear),
    .(prev),
    keyby = .(year, sex, age)
  ]
  prev_diet_dt <- dat_par$prev_diet
  # Loop over years
  for (starttime in (startyear + 1):endyear) {
    ## Check whether to stop simulation when running from Shiny
    if (ui) {
      t0 <- Sys.time()
      check_ui_status(scen_path, ui_lang)
      log_year_start(scen_path, ui_lang, starttime)
    }

    print(paste0("year: ", starttime))

    year_context <- setup_year_context(
      starttime,
      dat,
      dat_par,
      nsteps,
      max_year_mortality
    )
    #TODO: Possible to vectorize over sex and age?
    sex_age_combos <- CJ(
      sex = unique(year_context$stocks$sex),
      age = unique(year_context$stocks$age)
    )
    # Global parameters (static for all sex/age combinations)
    global_params <- c(
      rr_cancer_smoking = unname(baseline_parameters$rr["rr_cancer_smoking"]),
      rr_cancer_inactivity = unname(baseline_parameters$rr[
        "rr_cancer_inactivity"
      ]),
      rr_cancer_bmi = unname(baseline_parameters$rr["rr_cancer_bmi"]),
      rr_cancer_alcohol = unname(baseline_parameters$rr["rr_cancer_alcohol"]),
      rr_cvd_smoking = unname(baseline_parameters$rr["rr_cvd_smoking"]),
      rr_cvd_inactivity = unname(baseline_parameters$rr["rr_cvd_inactivity"]),
      rr_cvd_bmi = unname(baseline_parameters$rr["rr_cvd_bmi"]),
      rr_cvd_alcohol = unname(baseline_parameters$rr["rr_cvd_alcohol"]),
      paf_cancer_uv = unname(baseline_parameters$rr["paf_cancer_uv"]),
      paf_cvd_uv = unname(baseline_parameters$rr["paf_cvd_uv"]),
      comm_smoking = unname(baseline_parameters$communalities["smoking"]),
      comm_alcohol = unname(baseline_parameters$communalities["alcohol"]),
      comm_bmi = unname(baseline_parameters$communalities["bmi"]),
      comm_inactivity = unname(baseline_parameters$communalities["inactivity"])
    )
    # Does not depend on age/sex
    #cat("Calculating cfact_food_use\n")
    cfact_food_use <- mapply(
      piecewise,
      x = year_context$FINISH,
      x0 = cfact_startyear,
      y0 = cfact_food_0, # Baseline (static)
      x1 = cfact_endyear,
      y1 = cfact_food # Target (from input)
    )

    # TODO: Bottleneck is create_avec.
    # Precompute and create lookup tables for as
    # much as possible.

    ode_inputs <- lapply(
      1:nrow(sex_age_combos),
      function(i) {
        sex_ <- sex_age_combos$sex[i]
        age_ <- sex_age_combos$age[i]

        yvec <- create_yvec(
          year_context$stocks,
          sex_,
          age_,
          year_context$stock_cols
        )
        avec <- create_avec(
          yvec,
          sex_,
          age_,
          year_context$FINISH,
          year_context$death_and_mortality_rates,
          rate_immig_dt,
          rate_emig_dt,
          rate_born_dt,
          prev_smoking_dt,
          prev_alcohol_dt,
          prev_inactivity_dt,
          prev_bmi_dt,
          prev_diet_dt,
          calibration_results$heal_rate,
          calibration_results$constants,
          calibration_results$paf_other_calibration_0,
          global_params,
          year_context$n_pop,
          baseline_parameters$rr_cancer_diet,
          baseline_parameters$rr_cvd_diet,
          is_baseline,
          calibration_results$calibration_mortality_table,
          cfact_startyear,
          cfact_endyear,
          cfact,
          cfact_food_use,
          food_cat
        )

        list(
          yvec = yvec,
          avec = avec
        )
      }
    )

    # Extract results
    yvec_matrix <- do.call(cbind, lapply(ode_inputs, `[[`, "yvec"))
    avec_list <- lapply(ode_inputs, `[[`, "avec")

    o_list <- lapply(1:nrow(sex_age_combos), function(i) {
      ode(
        y = yvec_matrix[, i],
        times = year_context$simtime,
        func = model,
        parms = avec_list[[i]],
        method = "euler"
      )
    })

    if (debug_ode) {
      log_debug_ode(
        o_list,
        sex_age_combos,
        year_context$FINISH,
        year_context$simtime,
        PROJECTROOT
      )
    }
    # Process all ODE outputs
    update_list <- lapply(1:length(o_list), function(i) {
      create_update_row(
        o = o_list[[i]],
        sex_ = sex_age_combos$sex[i],
        age_ = sex_age_combos$age[i],
        FINISH = year_context$FINISH,
        avec = ode_inputs[[i]]$avec,
        dat = dat,
        calibration_results = calibration_results,
        baseline_parameters = baseline_parameters
      )
    })

    # Merge into dat
    key_cols <- c("year", "sex", "age")
    update_dt <- as.data.table(do.call(rbind, update_list))

    # Ensure both data.tables have the same columns (only the common ones)
    common_cols <- intersect(names(dat), names(update_dt))
    update_dt[, (key_cols) := lapply(.SD, as.integer), .SDcols = key_cols]
    # Perform an update join: modify dat in-place for matching rows
    setkeyv(update_dt, key_cols)
    # Use i. prefix to reference columns from update_dt
    dat[
      update_dt,
      (common_cols) := mget(paste0("i.", common_cols))
    ]
    setkeyv(dat, key_cols)
    ## Save data to file and write status and log if run by Shiny

    if (ui) {
      write_simulation_output(ui_lang, scen_path, dat, starttime, t0)
    }
  } # End loop year
  write_simulation_results(
    dat,
    write_data_to_file,
    is_baseline,
    ui,
    ui_lang,
    scen_path,
    PROJECTROOT,
    timestamp,
    baseline_parameters_path
  )

  return(dat)
}
