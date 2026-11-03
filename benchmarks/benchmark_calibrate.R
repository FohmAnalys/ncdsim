library(bench)
library(data.table)
library(jsonlite)
###########################################################################
# Setup for benchmarking
source("function_params.R") # reads the "default" parameters of "simulate_model"
source("ncdsim.R") # reads the paf function
source("calibrate.R")

baseline_parameters <- read_json(baseline_parameters_path)

calibration_params <- list(
  calibrate_cpaf_cancer = as.integer(baseline_parameters$calibrate_cpaf_cancer),
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
  dcost_total_base_year = as.integer(baseline_parameters$dcost_total_base_year),
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
###########################################################################
# Benchmarking base input

read_base_input <- function(input_path) {
  # --- Note ---:
  # Mostly avoids renaming columns multiple times
  # not necessarily faster but cleaner.
  # --- Helper: Read CSV and rename columns ---
  read_csv <- function(file, new_names = NULL, sep = ";") {
    dt <- fread(paste0(input_path, file), sep = sep)
    if (!is.null(new_names)) {
      setnames(dt, new_names)
    }
    dt
  }

  # --- 1. Read raw files ---
  costs <- read_csv("costs_direct.csv", c("age", 1, 2, "ncd"))
  ind_costs <- read_csv("costs_indirect.csv", c("age", 1, 2, "ncd"))
  population <- read_csv("pop_counts_scb.csv")
  preval_new <- read_csv("stocks.csv")
  incid_new <- read_csv("flows.csv")

  # --- 2. Reshape and clean ---
  preval_new <- dcast(preval_new, year + sex + age ~ stock, value.var = "N")
  incid_new <- incid_new[!flow %like% "NA_"]
  incid_new <- dcast(incid_new, year + sex + age ~ flow, value.var = "N")

  # Extract heal_new and drop unused flows
  heal_new <- incid_new[, .(
    year,
    age,
    sex,
    cancer_healthy,
    comorb_cancer,
    comorb_cvd,
    comorb_healthy,
    cvd_healthy
  )]
  incid_new[,
    c(
      "cancer_cancer",
      "cancer_healthy",
      "comorb_cancer",
      "comorb_comorb",
      "comorb_cvd",
      "comorb_healthy",
      "cvd_cvd",
      "cvd_healthy",
      "healthy_healthy"
    ) := NULL
  ]

  # --- 3. Create full tables with CJ and merge ---
  incid <- merge(
    CJ(year = min(incid_new$year):max(incid_new$year), age = 0:100, sex = 1:2),
    incid_new,
    by = c("year", "age", "sex"),
    all = TRUE
  )[,
    .(
      year,
      age,
      sex,
      cancer_comorb,
      cancer_cvd,
      cvd_cancer,
      cvd_comorb,
      healthy_cancer,
      healthy_comorb,
      healthy_cvd
    )
  ]
  heal <- merge(
    CJ(year = min(heal_new$year):max(heal_new$year), sex = 1:2, age = 0:100),
    heal_new,
    by = c("year", "age", "sex"),
    all = TRUE
  )
  preval <- merge(
    CJ(
      year = min(preval_new$year):max(preval_new$year),
      sex = 1:2,
      age = 0:100
    ),
    preval_new[, .(year, age, sex, cvd, cancer, comorb)],
    by = c("year", "age", "sex"),
    all = TRUE
  )

  # --- 4. Rename columns and fill NAs ---
  incid[, c("cancer_cvd", "cvd_cancer") := NULL]
  setnames(
    incid,
    c(
      "cancer_comorb",
      "cvd_comorb",
      "healthy_cancer",
      "healthy_comorb",
      "healthy_cvd"
    ),
    c("inc_cancer_cvd", "inc_cvd_cancer", "inc_cancer", "inc_comorb", "inc_cvd")
  )
  setnames(
    preval,
    c("cvd", "cancer", "comorb"),
    c("prev_cvd", "prev_cancer", "prev_comorb"),
  )
  setnames(
    heal,
    c(
      "cancer_healthy",
      "comorb_cancer",
      "comorb_cvd",
      "comorb_healthy",
      "cvd_healthy"
    ),
    c(
      "heal_cancer_healthy",
      "heal_comorb_cancer",
      "heal_comorb_cvd",
      "heal_comorb_healthy",
      "heal_cvd_healthy"
    )
  )

  setnafill(
    incid,
    "const",
    fill = 0,
    cols = c(
      "inc_cvd",
      "inc_cancer",
      "inc_cvd_cancer",
      "inc_cancer_cvd",
      "inc_comorb"
    )
  )
  setnafill(
    preval,
    "const",
    fill = 0,
    cols = c("prev_cvd", "prev_cancer", "prev_comorb")
  )
  setnafill(
    heal,
    "const",
    fill = 0,
    cols = c(
      "heal_cancer_healthy",
      "heal_comorb_cancer",
      "heal_comorb_cvd",
      "heal_comorb_healthy",
      "heal_cvd_healthy"
    )
  )

  # --- 5. Filter population and finalize preval ---
  list_of_years <- intersect(unique(incid$year), unique(preval$year))
  population <- population[
    year >= list_of_years[1] & year <= tail(list_of_years, 1)
  ]

  preval <- merge(
    population[, .(year, age, sex, pop)],
    preval,
    by = c("year", "age", "sex"),
    all = TRUE
  )[,
    frisk := pop - prev_cancer - prev_cvd - prev_comorb
  ][,
    (c("prev_cvd", "pop", "prev_cancer", "prev_comorb", "frisk")) := lapply(
      .SD,
      as.double
    ),
    .SDcols = c("prev_cvd", "pop", "prev_cancer", "prev_comorb", "frisk")
  ]

  # --- Return all processed tables ---
  list(
    costs = costs,
    ind_costs = ind_costs,
    population = population,
    incid = incid,
    heal = heal,
    preval = preval,
    list_of_years = list_of_years
  )
}

read_base_input_original <- function(input_path) {
  # read base input for costs, population, incidences and prevalences
  # extract and format data tables for incidences, prevalences, healing rates, KPP and population
  costs <- fread(paste0(input_path, "costs_direct.csv"), sep = ";")
  colnames(costs) <- c("age", 1, 2, "ncd")

  ind_costs <- fread(paste0(input_path, "costs_indirect.csv"), sep = ";")
  colnames(ind_costs) <- c("age", 1, 2, "ncd")

  population <- fread(paste0(input_path, "pop_counts_scb.csv"), sep = ";")

  preval_new <- fread(paste0(input_path, "stocks.csv"), sep = ';')
  incid_new <- fread(paste0(input_path, "flows.csv"), sep = ';')

  preval_new <- dcast(preval_new, year + sex + age ~ stock, value.var = 'N')
  incid_new <- incid_new[!flow %like% "NA_"]
  incid_new <- dcast(incid_new, year + sex + age ~ flow, value.var = 'N')
  colnames(preval_new) <- c(
    "AR",
    "Kon",
    "alder",
    "prev_cancer",
    "prev_comorb",
    "prev_cvd",
    "totpop"
  )
  heal_new <- incid_new[, c(
    "year",
    "sex",
    "age",
    "cancer_healthy",
    "comorb_cancer",
    "comorb_cvd",
    "comorb_healthy",
    "cvd_healthy"
  )]
  incid_new[,
    c(
      "cancer_cancer",
      "cancer_healthy",
      "comorb_cancer",
      "comorb_comorb",
      "comorb_cvd",
      "comorb_healthy",
      "cvd_cvd",
      "cvd_healthy",
      "healthy_healthy"
    ) := NULL
  ]

  colnames(incid_new) <- c(
    "AR",
    "Kon",
    "alder",
    "inc_cancer_cvd",
    "inc_cancer_to_cvd",
    "inc_cvd_to_cancer",
    "inc_cvd_cancer",
    "inc_cancer",
    "inc_comorb",
    "inc_cvd"
  )

  incid <- CJ(
    AR = min(incid_new[, (AR)]):max(incid_new[, (AR)]),
    Kon = 1:2,
    alder = 0:100
  )
  heal <- CJ(
    year = min(heal_new[, (year)]):max(heal_new[, (year)]),
    sex = 1:2,
    age = 0:100
  )
  preval <- CJ(
    AR = min(preval_new[, (AR)]):max(preval_new[, (AR)]),
    Kon = 1:2,
    alder = 0:100
  )
  incid <- merge(incid, incid_new, by = c("AR", "alder", "Kon"), all = TRUE)
  heal <- merge(heal, heal_new, by = c("year", "age", "sex"), all = TRUE)
  preval <- merge(
    preval,
    preval_new[, c(
      "AR",
      "alder",
      "Kon",
      "prev_cvd",
      "prev_cancer",
      "prev_comorb"
    )],
    by = c("AR", "alder", "Kon"),
    all = TRUE
  )

  incid[, c('inc_cancer_to_cvd', 'inc_cvd_to_cancer') := NULL]
  colnames(incid) <- c(
    "year",
    "age",
    "sex",
    "inc_cancer_cvd",
    "inc_cvd_cancer",
    "inc_cancer",
    "inc_comorb",
    "inc_cvd"
  )
  colnames(preval) <- c(
    "year",
    "age",
    "sex",
    "prev_cvd",
    "prev_cancer",
    "prev_comorb"
  )
  colnames(heal) <- c(
    "year",
    "age",
    "sex",
    "heal_cancer_healthy",
    "heal_comorb_cancer",
    "heal_comorb_cvd",
    "heal_comorb_healthy",
    "heal_cvd_healthy"
  )

  incid <- setnafill(
    incid,
    "const",
    fill = 0,
    cols = c(
      "inc_cvd",
      "inc_cancer",
      "inc_cvd_cancer",
      "inc_cancer_cvd",
      "inc_comorb"
    )
  )
  preval <- setnafill(
    preval,
    "const",
    fill = 0,
    cols = c("prev_cvd", "prev_cancer", "prev_comorb")
  )
  heal <- setnafill(
    heal,
    "const",
    fill = 0,
    cols = c(
      "heal_cancer_healthy",
      "heal_comorb_cancer",
      "heal_comorb_cvd",
      "heal_comorb_healthy",
      "heal_cvd_healthy"
    )
  )
  list_of_years <- intersect(unique(incid[, year]), unique(preval[, year]))
  population <- population[
    year >= list_of_years[1] &
      year <= list_of_years[length(list_of_years)]
  ]

  preval <- merge(
    population[, .(year, age, sex, pop)],
    preval,
    by = c("year", "age", "sex"),
    all = TRUE
  )
  preval[, frisk := pop - prev_cancer - prev_cvd - prev_comorb]
  cols_to_convert <- c("prev_cvd", "pop", "prev_cancer", "prev_comorb", "frisk")
  preval[,
    (cols_to_convert) := lapply(.SD, as.double),
    .SDcols = cols_to_convert
  ]
  # --- Return all processed tables ---
  list(
    costs = costs,
    ind_costs = ind_costs,
    population = population,
    incid = incid,
    heal = heal,
    preval = preval,
    list_of_years = list_of_years
  )
}

bm_base_input <- bench::mark(
  original <- read_base_input_original(input_path),
  optimized <- read_base_input(input_path),
  relative = TRUE
)

print(bm_base_input)

###########################################################################
# Benchmark calculate_riskgroup
calculate_riskgroup <- function(preval) {
  xpreval <- melt(preval, id.vars = c("year", "sex", "age"))[,
    cohort := year - age # Birth year
  ]

  # Set key for efficient shifting and calculate mean prevalence between years
  setkey(xpreval, variable, sex, cohort, year)
  xpreval[, value_prev_year := shift(value)] # Previous year's value

  # Set NA where cohort/sex/year are not consecutive
  xpreval[
    !(cohort == shift(cohort) & sex == shift(sex) & year == shift(year) + 1),
    value_prev_year := NA
  ]

  # Calculate mean prevalence for each stock
  xpreval[, meanvalue := (value + value_prev_year) / 2]

  # Reshape back to wide format (one row per year/sex/age)
  riskgroup <- dcast(
    xpreval[, .(year, sex, age, variable, rg = meanvalue)],
    year + sex + age ~ variable,
    value.var = "rg"
  )

  # Fill NAs with 0 for risk group prevalences
  riskgroup <- setnafill(
    riskgroup,
    "const",
    fill = 0,
    cols = c("prev_cvd", "prev_cancer", "prev_comorb", "frisk")
  )
  return(riskgroup)
}

calculate_riskgroup_original <- function(preval) {
  xpreval <- melt(preval, id.vars = c("year", "sex", "age"))
  xpreval[, cohort := year - age]
  setkey(xpreval, variable, sex, cohort, year)
  xpreval[, value_prev_year := shift(value)]
  xpreval[
    !(cohort == shift(cohort) & sex == shift(sex) & year == (shift(year) + 1)),
    value_prev_year := NA
  ]
  xpreval[, meanvalue := (value + value_prev_year) / 2]
  riskgroup <- dcast(
    xpreval[, .(year, sex, age, variable, rg = meanvalue)],
    year + sex + age ~ variable,
    value.var = "rg"
  )

  riskgroup <- setnafill(
    riskgroup,
    "const",
    0,
    cols = c("prev_cvd", "prev_cancer", "prev_comorb", "frisk")
  )

  return(riskgroup)
}


base_input <- read_base_input(input_path)

preval <- base_input$preval

bm_riskgroup <- bench::mark(
  original <- calculate_riskgroup_original(preval),
  refactored <- calculate_riskgroup(preval),
  relative = TRUE
)
print(bm_riskgroup)

###########################################################################
# Benchmark calculate_incidence_rates

calculate_incidence_rates <- function(riskgroup, incid) {
  # Merge incid and heal with riskgroup
  incid_rate <- incid[riskgroup, on = c("year", "age", "sex"), nomatch = NULL]

  # Calculate incidence rates (vectorized)
  incid_rate <- incid_rate[,
    ":="(
      ir_cvd = inc_cvd / frisk,
      ir_cancer = inc_cancer / frisk,
      ir_cvd_cancer = inc_cvd_cancer / prev_cvd,
      ir_cancer_cvd = inc_cancer_cvd / prev_cancer,
      ir_comorb = inc_comorb / frisk
    )
  ][
    # Drop unused columns
    ,
    c(
      "pop",
      "inc_cvd",
      "inc_cancer",
      "inc_cvd_cancer",
      "inc_cancer_cvd",
      "inc_comorb",
      "prev_cvd",
      "prev_cancer",
      "prev_comorb",
      "frisk"
    ) := NULL
  ][
    year > min(year)
  ]
  incid_rate <- setnafill(
    incid_rate,
    "const",
    fill = 0,
    cols = colnames(incid_rate)
  )

  return(incid_rate)
}

calculate_incidence_rates_original <- function(riskgroup, incid) {
  incid_rate <- incid[riskgroup, on = c("year", "age", "sex"), nomatch = NULL]

  incid_rate[, ir_cvd := inc_cvd / frisk]
  incid_rate[, ir_cancer := inc_cancer / frisk]
  incid_rate[, ir_cvd_cancer := inc_cvd_cancer / prev_cvd]
  incid_rate[, ir_cancer_cvd := inc_cancer_cvd / prev_cancer]
  incid_rate[, ir_comorb := inc_comorb / frisk]
  incid_rate[,
    c(
      "pop",
      "inc_cvd",
      "inc_cancer",
      "inc_cvd_cancer",
      "inc_cancer_cvd",
      "inc_comorb",
      "prev_cvd",
      "prev_cancer",
      "prev_comorb",
      "frisk"
    ) := NULL
  ]
  incid_rate <- setnafill(
    incid_rate,
    "const",
    fill = 0,
    cols = colnames(incid_rate)
  )
  incid_rate <- incid_rate[year > min(incid_rate$year)]
  return(incid_rate)
}

base_input <- read_base_input(input_path)
preval <- base_input$preval
riskgroup <- calculate_riskgroup_original(preval)
bm_incid_rate <- bench::mark(
  original <- calculate_incidence_rates_original(riskgroup, base_input$incid),
  refactored <- calculate_incidence_rates(riskgroup, base_input$incid),
  relative = TRUE
)
print(bm_incid_rate)

###########################################################################
# Benchmarking of calculate_heal_rates
calculate_heal_rates <- function(heal) {
  heal <- heal[riskgroup, on = c("year", "age", "sex"), nomatch = NULL]
  # Calculate healing rates
  heal <- heal[,
    `:=`(
      hr_cvd = heal_cvd_healthy / prev_cvd,
      hr_cancer = heal_cancer_healthy / prev_cancer,
      hr_comorb_cancer = heal_comorb_cancer / prev_comorb,
      hr_comorb_cvd = heal_comorb_cvd / prev_comorb,
      hr_comorb_healthy = heal_comorb_healthy / prev_comorb
    )
  ]

  # Aggregate healing rates by age/sex
  heal_rate <- heal[,
    .(
      cvd = mean(hr_cvd),
      cancer = mean(hr_cancer),
      comorb_cancer = mean(hr_comorb_cancer),
      comorb_cvd = mean(hr_comorb_cvd),
      comorb_healthy = mean(hr_comorb_healthy)
    ),
    by = c("age", "sex")
  ]

  # Fill NAs with 0
  heal_rate <- setnafill(
    heal_rate,
    "const",
    fill = 0,
    cols = c("cvd", "cancer", "comorb_cancer", "comorb_cvd", "comorb_healthy")
  )
}

calculate_heal_rates_original <- function(heal) {
  heal <- heal[riskgroup, on = c("year", "age", "sex"), nomatch = NULL]
  # calculate healing rates
  heal[, hr_cvd := heal_cvd_healthy / prev_cvd]
  heal[, hr_cancer := heal_cancer_healthy / prev_cancer]
  heal[, hr_comorb_cancer := heal_comorb_cancer / prev_comorb]
  heal[, hr_comorb_cvd := heal_comorb_cvd / prev_comorb]
  heal[, heal_comorb_healthy := heal_comorb_healthy / prev_comorb]
  heal_rate <- heal[,
    .(
      cvd = mean(hr_cvd),
      cancer = mean(hr_cancer),
      comorb_cancer = mean(hr_comorb_cancer),
      comorb_cvd = mean(hr_comorb_cvd),
      comorb_healthy = mean(heal_comorb_healthy)
    ),
    by = c('age', 'sex')
  ]
  heal_rate <- setnafill(
    heal_rate,
    "const",
    0,
    cols = c("cvd", "cancer", "comorb_cancer", "comorb_cvd", "comorb_healthy")
  )
  return(heal_rate)
}
base_input <- read_base_input(input_path)
preval <- base_input$preval
riskgroup <- calculate_riskgroup_original(preval)
bm_heal_rate <- bench::mark(
  original <- calculate_heal_rates_original(base_input$heal),
  refactored <- calculate_heal_rates(base_input$heal),
  relative = TRUE
)
print(bm_heal_rate)

###########################################################################
# Benchmarking of cpaf_function
calculate_cpaf_diet <- function(food_prevalences, rr_cancer_diet, rr_cvd_diet) {
  # Reshape food_prevalences to long format
  food_prevalences_long <- melt(
    food_prevalences,
    id.vars = c("sex", "age"),
    variable.name = "food_category",
    value.name = "prevalence"
  )

  # Prepare rr data
  rr_cancer_diet_dt <- data.table(
    food_category = names(rr_cancer_diet),
    rr = as.numeric(rr_cancer_diet),
    ncd = "cancer"
  )
  rr_cvd_diet_dt <- data.table(
    food_category = names(rr_cvd_diet),
    rr = as.numeric(rr_cvd_diet),
    ncd = "cvd"
  )

  rr_dt <- rbindlist(list(rr_cancer_diet_dt, rr_cvd_diet_dt))

  # Merge with food_prevalences_long
  pafs_diet <- merge(
    food_prevalences_long,
    rr_dt,
    by = "food_category",
    allow.cartesian = TRUE
  )

  pafs_diet[, paf := paf(prevalence, rr)]

  # Aggregate to compute cpaf for each sex, age, and ncd
  cpaf_diet <- pafs_diet[,
    .(cpaf = 1 - prod(1 - paf)),
    by = .(sex, age, ncd)
  ]
  setkeyv(cpaf_diet, c("sex", "age", "ncd"))
}

calculate_cpaf_diet_original <- function(
  food_prevalences,
  rr_cancer_diet,
  rr_cvd_diet
) {
  # read prevalences for risk factors and relative risks for food
  food_categories <- colnames(food_prevalences)[3:length(food_prevalences)]
  rr_cancer_diet <- data.frame(as.list(rr_cancer_diet))
  rr_cvd_diet <- data.frame(as.list(rr_cvd_diet))

  # calculate cPAF for food
  cpaf_diet <- CJ(sex = 1:2, age = 0:100, ncd = c("cancer", "cvd"))
  for (sex_ in 1:2) {
    prev_diet <- food_prevalences[sex == sex_, ..food_categories]
    pafs_diet_cancer <- rbindlist(apply(prev_diet, 1, paf, rr = rr_cancer_diet))
    cpaf_diet_cancer <- (1 - apply(1 - pafs_diet_cancer, 1, prod))
    pafs_diet_cvd <- rbindlist(apply(prev_diet, 1, paf, rr = rr_cvd_diet))
    cpaf_diet_cvd <- (1 - apply(1 - pafs_diet_cvd, 1, prod))
    cpaf_diet[ncd == "cancer" & sex == sex_, cpaf := cpaf_diet_cancer]
    cpaf_diet[ncd == "cvd" & sex == sex_, cpaf := cpaf_diet_cvd]
  }
  return(cpaf_diet)
}

base_input <- read_base_input(input_path)
preval <- base_input$preval
riskgroup <- calculate_riskgroup_original(preval)
prevalences <- read_prevalences(input_path, base_input$list_of_years)

bm_cpaf <- bench::mark(
  original <- calculate_cpaf_diet_original(
    prevalences$food_prevalences,
    calibration_params$rr_cancer_diet,
    calibration_params$rr_cvd_diet
  ),
  optimized <- calculate_cpaf_diet(
    prevalences$food_prevalences,
    calibration_params$rr_cancer_diet,
    calibration_params$rr_cvd_diet
  ),
  check = \(x, y) {
    all.equal(x, y, check.attributes = FALSE)
  },
  min_iterations = 10,
  relative = TRUE
)

print(bm_cpaf)

###########################################################################
# Benchmarking of cpafs_per_year functions

calculate_cpafs_per_year <- function(
  list_of_years,
  communalities,
  cpaf_diet,
  prevalences,
  rr
) {
  # --- 1. Stack and reshape all prevalence tables ---
  prev_stacked <- rbindlist(
    lapply(c("alcohol", "smoking", "bmi", "inactivity"), function(factor) {
      dt <- prevalences[[(paste0("prev_", factor))]]
      dt[, factor := factor]
    }),
    idcol = FALSE
  )

  prev_long_all <- melt(
    prev_stacked,
    id.vars = c("sex", "age", "factor"),
    variable.name = "year",
    value.name = "prevalence"
  )[, year := as.character(year)]

  # --- 2. Prepare risk factor parameters (cancer + CVD) ---
  risk_factor_params <- data.table(
    factor = rep(c("alcohol", "smoking", "bmi", "inactivity"), each = 2),
    ncd = rep(c("cancer", "cvd"), times = 4),
    rr = c(
      rr["rr_cancer_alcohol"],
      rr["rr_cvd_alcohol"],
      rr["rr_cancer_smoking"],
      rr["rr_cvd_smoking"],
      rr["rr_cancer_bmi"],
      rr["rr_cvd_bmi"],
      rr["rr_cancer_inactivity"],
      rr["rr_cvd_inactivity"]
    ),
    communality = rep(
      as.numeric(communalities[c("alcohol", "smoking", "bmi", "inactivity")]),
      each = 2
    )
  )

  # --- 3. Merge and calculate PAF for lifestyle factors ---
  prev_long_all <- merge(
    prev_long_all,
    risk_factor_params,
    by = "factor",
    allow.cartesian = TRUE
  )
  prev_long_all[, cpaf := paf(prevalence, rr) * communality]

  # --- 4. Add diet PAF (from cpaf_diet) ---
  n_years <- length(list_of_years)
  cpaf_diet_year <- cpaf_diet[,
    .(sex, age, ncd, cpaf) # Select columns
  ][
    rep(1:.N, each = n_years) # Repeat each row for every year
  ][,
    year := rep(as.character(list_of_years), times = .N / n_years) # Assign years
  ][,
    factor := "diet" # Add factor column
  ]

  # --- 5. Combine all PAFs (diet + lifestyle) ---
  all_pafs <- rbindlist(list(
    cpaf_diet_year[, .(sex, age, year, ncd, cpaf, factor)],
    prev_long_all[, .(sex, age, year, ncd, cpaf, factor)]
  ))

  # --- 6. Calculate combined cpaf for each group ---
  cpaf_combined <- all_pafs[,
    .(cpaf = 1 - prod(1 - cpaf, na.rm = TRUE)),
    by = .(sex, age, year, ncd)
  ]

  # --- 7. Split into paf_cancer, paf_cvd, and calculate means ---
  paf_cancer <- cpaf_combined[ncd == "cancer", .(sex, age, year, cpaf)]
  paf_cvd <- cpaf_combined[ncd == "cvd", .(sex, age, year, cpaf)]

  paf_comorb <- merge(
    paf_cancer[, .(sex, age, year, cpaf_cancer = cpaf)],
    paf_cvd[, .(sex, age, year, cpaf_cvd = cpaf)],
    by = c("sex", "age", "year")
  )[, cpaf := pmax(cpaf_cancer, cpaf_cvd)]
  paf_comorb <- paf_comorb[, .(sex, age, year, cpaf)]

  cpaf_cvd_mean <- paf_cvd[, .(cpaf = mean(cpaf)), by = .(age, sex)]
  cpaf_cancer_mean <- paf_cancer[, .(cpaf = mean(cpaf)), by = .(age, sex)]
  cpaf_cancer_cvd_mean <- paf_comorb[, .(cpaf = mean(cpaf)), by = .(age, sex)]

  return(
    list(
      paf_cvd = paf_cvd,
      paf_cancer = paf_cancer,
      paf_comorb = paf_comorb,
      cpaf_cvd_mean = cpaf_cvd_mean,
      cpaf_cancer_mean = cpaf_cancer_mean,
      cpaf_cancer_cvd_mean = cpaf_cancer_cvd_mean
    )
  )
}

calculate_cpafs_per_year_original <- function(
  list_of_years,
  communalities,
  cpaf_diet,
  prevalences,
  rr
) {
  # calculate cPAF per year for lifestyle factors
  paf_cancer <- CJ(sex = 1:2, age = 0:100, year = as.character(list_of_years))
  paf_cvd <- CJ(sex = 1:2, age = 0:100, year = as.character(list_of_years))
  for (sex_ in 1:2) {
    for (year_ in as.character(list_of_years)) {
      paf_diet_c_ <- cpaf_diet[sex == sex_ & ncd == "cancer", cpaf]
      paf_alcohol_c <- unname((communalities["alcohol"])) *
        paf(
          prevalences$prev_alcohol[sex == sex_, ..year_],
          unname(rr["rr_cancer_alcohol"])
        )
      paf_smoking_c <- unname((communalities["smoking"])) *
        paf(
          prevalences$prev_smoking[sex == sex_, ..year_],
          unname(rr["rr_cancer_smoking"])
        )
      paf_bmi_c <- unname((communalities["smoking"])) *
        paf(
          prevalences$prev_bmi[sex == sex_, ..year_],
          unname(rr["rr_cancer_bmi"])
        )
      paf_inactivity_c <- unname((communalities["inactivity"])) *
        paf(
          prevalences$prev_inactivity[sex == sex_, ..year_],
          unname(rr["rr_cancer_inactivity"])
        )
      cpaf_c <- 1 -
        (1 - paf_diet_c_) *
          (1 - paf_alcohol_c) *
          (1 - paf_smoking_c) *
          (1 - paf_bmi_c) *
          (1 - paf_inactivity_c)
      paf_cancer[sex == sex_ & year == year_, cpaf := cpaf_c]

      paf_diet_cvd_ <- cpaf_diet[sex == sex_ & ncd == "cvd", cpaf]
      paf_alcohol_cvd <- unname((communalities["alcohol"])) *
        paf(
          prevalences$prev_alcohol[sex == sex_, ..year_],
          unname(rr["rr_cvd_alcohol"])
        )
      paf_smoking_cvd <- unname((communalities["smoking"])) *
        paf(
          prevalences$prev_smoking[sex == sex_, ..year_],
          unname(rr["rr_cvd_smoking"])
        )
      paf_bmi_cvd <- unname((communalities["smoking"])) *
        paf(
          prevalences$prev_bmi[sex == sex_, ..year_],
          unname(rr["rr_cvd_bmi"])
        )
      paf_inactivity_cvd <- unname((communalities["inactivity"])) *
        paf(
          prevalences$prev_inactivity[sex == sex_, ..year_],
          unname(rr["rr_cvd_inactivity"])
        )
      cpaf_cvd <- 1 -
        (1 - paf_diet_cvd_) *
          (1 - paf_alcohol_cvd) *
          (1 - paf_smoking_cvd) *
          (1 - paf_bmi_cvd) *
          (1 - paf_inactivity_cvd)
      paf_cvd[sex == sex_ & year == year_, cpaf := cpaf_cvd]
    }
  }
  paf_comorb <- paf_cvd[, cpaf := pmax(paf_cancer$cpaf, paf_cvd$cpaf)]
  # calculate the total yearly mean cPAF (dvs food + non food)
  cpaf_cvd_mean <- paf_cvd[, .(cpaf = mean(cpaf)), by = .(age, sex)]
  cpaf_cancer_mean <- paf_cancer[, .(cpaf = mean(cpaf)), by = .(age, sex)]
  cpaf_cancer_cvd_mean <- paf_comorb[, .(cpaf = mean(cpaf)), by = .(age, sex)]
  return(
    list(
      paf_cvd = paf_cvd,
      paf_cancer = paf_cancer,
      paf_comorb = paf_comorb,
      cpaf_cvd_mean = cpaf_cvd_mean,
      cpaf_cancer_mean = cpaf_cancer_mean,
      cpaf_cancer_cvd_mean = cpaf_cancer_cvd_mean
    )
  )
}

cpaf_diet <- calculate_cpaf_diet(
  prevalences$food_prevalences,
  calibration_params$rr_cancer_diet,
  calibration_params$rr_cvd_diet
)

bm_cpaf <- bench::mark(
  original <- calculate_cpafs_per_year_original(
    base_input$list_of_years,
    calibration_params$communalities,
    cpaf_diet,
    prevalences,
    calibration_params$rr
  ),
  optimized <- calculate_cpafs_per_year(
    base_input$list_of_years,
    calibration_params$communalities,
    cpaf_diet,
    prevalences,
    calibration_params$rr
  ),
  check = \(x, y) {
    all.equal(x, y, check.attributes = FALSE)
  },
  relative = TRUE
)

print(bm_cpaf)

###########################################################################
# Benchmark calculate_constants
calculate_constants <- function(params, cpafs, ir_mean) {
  constants_df <- CJ(age = 0:100, sex = 1:2)

  # Pre-calculate age ranges
  age_min_cancer <- params$age_cutoff_cancer
  age_max_cancer <- params$age_cutoff_cancer_high
  age_min_cvd <- params$age_cutoff_cvd
  age_max_cvd <- params$age_cutoff_cvd_high
  age_min_comorb <- max(age_min_cancer, age_min_cvd)
  age_max_comorb <- max(age_max_cancer, age_max_cvd)

  # Loop over sexes (only 2 iterations, negligible overhead)
  for (sex_ in 1:2) {
    # Filter data for this sex and age ranges
    pf_c <- cpafs$cpaf_cancer_mean[
      sex == sex_ & age >= age_min_cancer & age <= age_max_cancer
    ]
    pf_cvd <- cpafs$cpaf_cvd_mean[
      sex == sex_ & age >= age_min_cvd & age <= age_max_cvd
    ]
    pf_c_cvd <- cpafs$cpaf_cancer_cvd_mean[
      sex == sex_ & age >= age_min_comorb & age <= age_max_comorb
    ]

    # Cancer constants
    if (params$calibrate_cpaf_cancer == 0) {
      c_c <- ir_mean[
        sex == sex_ & age >= age_min_cancer & age <= age_max_cancer,
        ir_cancer
      ]
      c_cvd_c <- ir_mean[
        sex == sex_ & age >= age_min_cancer & age <= age_max_cancer,
        ir_cvd_cancer
      ]
    } else {
      mean_cpaf_c <- mean(pf_c$cpaf, na.rm = TRUE)
      c_c <- (params$prop_cancer *
        ir_mean[
          sex == sex_ & age >= age_min_cancer & age <= age_max_cancer,
          ir_cancer
        ]) /
        mean_cpaf_c
      c_cvd_c <- (params$prop_cancer *
        ir_mean[
          sex == sex_ & age >= age_min_cancer & age <= age_max_cancer,
          ir_cvd_cancer
        ]) /
        mean_cpaf_c
    }

    # CVD constants
    if (params$calibrate_cpaf_cvd == 0) {
      c_cvd <- ir_mean[
        sex == sex_ & age >= age_min_cvd & age <= age_max_cvd,
        ir_cvd
      ]
      c_c_cvd <- ir_mean[
        sex == sex_ & age >= age_min_cvd & age <= age_max_cvd,
        ir_cancer_cvd
      ]
      c_c_o_cvd <- ir_mean[
        sex == sex_ & age >= age_min_comorb & age <= age_max_comorb,
        ir_comorb
      ]
    } else {
      mean_cpaf_cvd <- mean(pf_cvd$cpaf, na.rm = TRUE)
      mean_cpaf_c_cvd <- mean(pf_c_cvd$cpaf, na.rm = TRUE)
      c_cvd <- (params$prop_cvd *
        ir_mean[
          sex == sex_ & age >= age_min_cvd & age <= age_max_cvd,
          ir_cvd
        ]) /
        mean_cpaf_cvd
      c_c_cvd <- (params$prop_cvd *
        ir_mean[
          sex == sex_ & age >= age_min_cvd & age <= age_max_cvd,
          ir_cancer_cvd
        ]) /
        mean_cpaf_cvd
      c_c_o_cvd <- (mean(params$prop_cvd, params$prop_cancer) *
        ir_mean[
          sex == sex_ & age >= age_min_comorb & age <= age_max_comorb,
          ir_comorb
        ]) /
        mean_cpaf_c_cvd
    }

    # Update constants_df (direct assignment, no merge overhead)
    constants_df[
      sex == sex_ & age >= age_min_cancer & age <= age_max_cancer,
      `:=`(const_cancer = c_c, const_cvd_comorb = c_cvd_c)
    ]
    constants_df[
      sex == sex_ & age >= age_min_cvd & age <= age_max_cvd,
      `:=`(const_cvd = c_cvd, const_cancer_comorb = c_c_cvd)
    ]
    constants_df[
      sex == sex_ & age >= age_min_comorb & age <= age_max_comorb,
      const_comorb := c_c_o_cvd
    ]
  }

  # Fill NAs with 0
  setnafill(
    constants_df,
    "const",
    fill = 0.0,
    cols = c(
      "const_cancer",
      "const_cvd",
      "const_cancer_comorb",
      "const_cvd_comorb",
      "const_comorb"
    )
  )

  constants_df
}

calculate_constants_original <- function(params, cpafs, ir_mean) {
  constants_df <- CJ(age = 0:100, sex = 1:2)

  for (sex_ in 1:2) {
    pf_c <- cpafs$cpaf_cancer_mean[
      sex == sex_ &
        age >= params$age_cutoff_cancer &
        age <= params$age_cutoff_cancer_high
    ]
    pf_cvd <- cpafs$cpaf_cvd_mean[
      sex == sex_ &
        age >= params$age_cutoff_cvd &
        age <= params$age_cutoff_cvd_high
    ]
    pf_c_cvd <- cpafs$cpaf_cancer_cvd_mean[
      sex == sex_ &
        age >= max(params$age_cutoff_cvd, params$age_cutoff_cancer) &
        age <= max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high)
    ]

    if (params$calibrate_cpaf_cancer == 0) {
      c_c <- ir_mean[
        sex == sex_ &
          age >= params$age_cutoff_cancer &
          age <= params$age_cutoff_cancer_high,
        ir_cancer
      ]
      c_cvd_c <- ir_mean[
        sex == sex_ &
          age >= params$age_cutoff_cancer &
          age <= params$age_cutoff_cancer_high,
        ir_cvd_cancer
      ]
    } else {
      c_c <- (prop_cancer *
        ir_mean[
          sex == sex_ &
            age >= params$age_cutoff_cancer &
            age <= params$age_cutoff_cancer_high,
          ir_cancer
        ]) /
        mean(pf_c$cpaf)
      c_cvd_c <- (prop_cancer *
        ir_mean[
          sex == sex_ &
            age >= params$age_cutoff_cancer &
            age <= params$age_cutoff_cancer_high,
          ir_cvd_cancer
        ]) /
        mean(pf_c$cpaf)
    }
    constants_df[
      sex == sex_ &
        age >= params$age_cutoff_cancer &
        age <= params$age_cutoff_cancer_high,
      const_cancer := c_c
    ]
    constants_df[
      sex == sex_ &
        age >= params$age_cutoff_cancer &
        age <= params$age_cutoff_cancer_high,
      const_cvd_comorb := c_cvd_c
    ]

    if (params$calibrate_cpaf_cvd == 0) {
      c_cvd <- ir_mean[
        sex == sex_ &
          age >= params$age_cutoff_cvd &
          age <= params$age_cutoff_cvd_high,
        ir_cvd
      ]
      c_c_cvd <- ir_mean[
        sex == sex_ &
          age >= params$age_cutoff_cvd &
          age <= params$age_cutoff_cvd_high,
        ir_cancer_cvd
      ]
      c_c_o_cvd <- ir_mean[
        sex == sex_ &
          age >= max(params$age_cutoff_cvd, params$age_cutoff_cancer) &
          age <= max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high),
        ir_comorb
      ]
    } else {
      c_cvd <- (prop_cvd *
        params$ir_mean[
          sex == sex_ &
            age >= params$age_cutoff_cvd &
            age <= params$age_cutoff_cvd_high,
          params$ir_cvd
        ]) /
        mean(pf_cvd$cpaf)
      c_c_cvd <- (prop_cvd *
        params$ir_mean[
          sex == sex_ &
            age >= params$age_cutoff_cvd &
            age <= params$age_cutoff_cvd_high,
          params$ir_cancer_cvd
        ]) /
        mean(pf_cvd$cpaf)
      c_c_o_cvd <- (mean(params$prop_cvd, params$prop_cancer) *
        params$ir_mean[
          sex == sex_ &
            age >= max(params$age_cutoff_cvd, params$age_cutoff_cancer) &
            age <=
              max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high),
          params$ir_comorb
        ]) /
        mean(pf_c_cvd$cpaf)
    }
    constants_df[
      sex == sex_ &
        age >= params$age_cutoff_cvd &
        age <= params$age_cutoff_cvd_high,
      const_cvd := c_cvd
    ]
    constants_df[
      sex == sex_ &
        age >= params$age_cutoff_cvd &
        age <= params$age_cutoff_cvd_high,
      const_cancer_comorb := c_c_cvd
    ]

    constants_df[
      sex == sex_ &
        age >= max(params$age_cutoff_cvd, params$age_cutoff_cancer) &
        age <= max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high),
      const_comorb := c_c_o_cvd
    ]
  }
  setnafill(
    constants_df,
    "const",
    fill = 0.0,
    cols = c(
      "const_cancer",
      "const_cvd",
      "const_cancer_comorb",
      "const_cvd_comorb",
      "const_comorb"
    )
  )
}


# Calculate riskgroup
riskgroup <- calculate_riskgroup(preval)
# Calculate incidence rate for flows
incid_rate <- calculate_incidence_rates(riskgroup, base_input$incid)
# Calculate healing rates
heal_rate <- calculate_heal_rates(base_input$heal)
# Read prevalences
prevalences <- read_prevalences(input_path, base_input$list_of_years)
# calculate cPAF for food
cpaf_diet <- calculate_cpaf_diet(
  prevalences$food_prevalences,
  calibration_params$rr_cancer_diet,
  calibration_params$rr_cvd_diet
)
cpafs_per_year <- calculate_cpafs_per_year(
  base_input$list_of_years,
  calibration_params$communalities,
  cpaf_diet,
  prevalences,
  calibration_params$rr
)

ir_mean <- incid_rate[,
  .(
    ir_cvd = mean(ir_cvd),
    ir_cancer = mean(ir_cancer),
    ir_cvd_cancer = mean(ir_cvd_cancer),
    ir_cancer_cvd = mean(ir_cancer_cvd),
    ir_comorb = mean(ir_comorb)
  ),
  by = c('age', 'sex')
]
bm_constants <- bench::mark(
  original <- calculate_constants_original(
    calibration_params,
    cpafs_per_year,
    ir_mean
  ),
  refactored <- calculate_constants(
    calibration_params,
    cpafs_per_year,
    ir_mean
  ),
  relative = TRUE
)
print(bm_constants)

###########################################################################
# Benchmark calculate_paf_other

calculate_paf_other <- function(
  list_of_years,
  params,
  constants_df,
  ir_mean,
  cpafs
) {
  paf_other <- CJ(sex = 1:2, age = 0:100, year = as.character(list_of_years))
  # Pre-calculate age ranges
  age_min_cancer <- params$age_cutoff_cancer
  age_max_cancer <- params$age_cutoff_cancer_high
  age_min_cvd <- params$age_cutoff_cvd
  age_max_cvd <- params$age_cutoff_cvd_high
  age_min_comorb <- max(age_min_cancer, age_min_cvd)
  age_max_comorb <- max(age_max_cancer, age_max_cvd)

  # Rename PAF columns to avoid conflicts during joins
  setnames(cpafs$paf_cancer, "cpaf", "paf_cancer")
  setnames(cpafs$paf_cvd, "cpaf", "paf_cvd")
  setnames(cpafs$paf_comorb, "cpaf", "paf_comorb")

  # Join all static and dynamic data into paf_other
  paf_other <- paf_other[
    constants_df,
    on = c("sex", "age") # Join constants (static)
  ][
    ir_mean,
    on = c("sex", "age") # Join incidence rates (static)
  ][
    cpafs$paf_cancer,
    on = c("year", "sex", "age") # Join PAF cancer (dynamic)
  ][
    cpafs$paf_cvd,
    on = c("year", "sex", "age") # Join PAF CVD (dynamic)
  ][
    cpafs$paf_comorb,
    on = c("year", "sex", "age") # Join PAF comorb (dynamic)
  ]

  # --- Cancer ---
  paf_other[
    age >= age_min_cancer & age <= age_max_cancer,
    cpaf_cancer := ir_cancer - const_cancer * paf_cancer
  ]
  paf_other[
    age < age_min_cancer | age > age_max_cancer,
    cpaf_cancer := ir_cancer
  ]

  # --- CVD ---
  paf_other[
    age >= age_min_cvd & age <= age_max_cvd,
    cpaf_cvd := ir_cvd - const_cvd * paf_cvd
  ]
  paf_other[
    age < age_min_cvd | age > age_max_cvd,
    cpaf_cvd := ir_cvd
  ]

  # --- Cancer + CVD (comorbidity) ---
  paf_other[
    age >= age_min_cvd & age <= age_max_cvd,
    cpaf_cancer_comorb := ir_cancer_cvd - const_cancer_comorb * paf_cvd
  ]
  paf_other[
    age < age_min_cvd | age > age_max_cvd,
    cpaf_cancer_comorb := ir_cancer_cvd
  ]

  # --- CVD + Cancer (comorbidity) ---
  paf_other[
    age >= age_min_cancer & age <= age_max_cancer,
    cpaf_cvd_comorb := ir_cvd_cancer - const_cvd_comorb * paf_cancer
  ]
  paf_other[
    age < age_min_cancer | age > age_max_cancer,
    cpaf_cvd_comorb := ir_cvd_cancer
  ]

  # --- Comorbidity (all 3) ---
  paf_other[
    age >= age_min_comorb & age <= age_max_comorb,
    cpaf_comorb := ir_comorb - const_comorb * paf_comorb
  ]
  paf_other[
    age < max(age_min_cvd, age_min_cancer) |
      age > max(age_max_cvd, age_max_cancer),
    cpaf_comorb := ir_comorb
  ]
  paf_other[,
    c(
      "const_cancer",
      "const_cvd_comorb",
      "const_cvd",
      "const_cancer_comorb",
      "const_comorb",
      "ir_cvd",
      "ir_cancer",
      "ir_cvd_cancer",
      "ir_cancer_cvd",
      "ir_comorb",
      "paf_cancer",
      "paf_cvd",
      "paf_comorb"
    ) := NULL
  ]

  # This is a weird fix but the original code loops over list_of_years[2:length(list_of_years)],
  # leaving all NAs where year == list_of_years[1]
  paf_other[
    year == list_of_years[1],
    c(
      "cpaf_cancer",
      "cpaf_cvd",
      "cpaf_cancer_comorb",
      "cpaf_cvd_comorb",
      "cpaf_comorb"
    ) := NA
  ]
  melt(
    paf_other,
    id.vars = c("sex", "age", "year"),
    variable.name = "disease",
    value.name = "cpaf"
  )
}

calculate_paf_other_original <- function(
  list_of_years,
  params,
  constants_df,
  ir_mean,
  cpafs
) {
  # Calculate paf_other:
  # if age >= age_cutoff then paf_other is the complement to paf from lifestyle
  # if age < age_cutoff then paf_other = incidence rate
  paf_other <- CJ(sex = 1:2, age = 0:100, year = as.character(list_of_years))

  for (sex_ in 1:2) {
    for (y in list_of_years[2:length(list_of_years)]) {
      c_c <- constants_df[
        sex == sex_ &
          age >= params$age_cutoff_cancer &
          age <= params$age_cutoff_cancer_high,
        const_cancer
      ]
      ir_c <- ir_mean[
        sex == sex_ &
          age >= params$age_cutoff_cancer &
          age <= params$age_cutoff_cancer_high,
        ir_cancer
      ]
      c_cvd <- constants_df[
        sex == sex_ &
          age >= params$age_cutoff_cvd &
          age <= params$age_cutoff_cvd_high,
        const_cvd
      ]
      ir_cvd <- ir_mean[
        sex == sex_ &
          age >= params$age_cutoff_cvd &
          age <= params$age_cutoff_cvd_high,
        ir_cvd
      ]
      c_cvd_cancer <- constants_df[
        sex == sex_ &
          age >= params$age_cutoff_cancer &
          age <= params$age_cutoff_cancer_high,
        const_cvd_comorb
      ]
      ir_cvd_cancer <- ir_mean[
        sex == sex_ &
          age >= params$age_cutoff_cancer &
          age <= params$age_cutoff_cancer_high,
        ir_cvd_cancer
      ]
      c_c_cvd <- constants_df[
        sex == sex_ &
          age >= params$age_cutoff_cvd &
          age <= params$age_cutoff_cvd_high,
        const_cancer_comorb
      ]
      ir_c_cvd <- ir_mean[
        sex == sex_ &
          age >= params$age_cutoff_cvd &
          age <= params$age_cutoff_cvd_high,
        ir_cancer_cvd
      ]
      c_c_o_cvd <- constants_df[
        sex == sex_ &
          age >= max(params$ge_cutoff_cvd, params$age_cutoff_cancer) &
          age <= max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high),
        const_comorb
      ]
      ir_c_o_cvd <- ir_mean[
        sex == sex_ &
          age >= max(params$age_cutoff_cvd, params$ge_cutoff_cancer) &
          age <= max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high),
        ir_comorb
      ]

      cut_c <- cpafs$paf_cancer[year == y & sex == sex_, ]
      paf_slice_c <- cut_c[
        age >= params$age_cutoff_cancer & age <= params$age_cutoff_cancer_high,
        cpaf
      ]
      cut_cvd <- cpafs$paf_cvd[year == y & sex == sex_, ]
      paf_slice_cvd <- cut_cvd[
        age >= params$age_cutoff_cvd & age <= params$age_cutoff_cvd_high,
        cpaf
      ]
      cut_c_o_cvd <- cpafs$paf_comorb[year == y & sex == sex_, ]
      paf_slice_c_o_cvd <- cut_c_o_cvd[
        age >= max(params$age_cutoff_cvd, params$age_cutoff_cancer) &
          age <= max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high),
        cpaf
      ]

      lb_c <- ir_c - c_c * paf_slice_c
      res_c <- ir_mean[
        sex == sex_ &
          (age < params$age_cutoff_cancer |
            age > params$age_cutoff_cancer_high),
        ir_cancer
      ]
      paf_other[
        year == y &
          sex == sex_ &
          age >= params$age_cutoff_cancer &
          age <= params$age_cutoff_cancer_high,
        cpaf_cancer := lb_c
      ]
      paf_other[
        year == y &
          sex == sex_ &
          (age < params$age_cutoff_cancer |
            age > params$age_cutoff_cancer_high),
        cpaf_cancer := res_c
      ]
      lb_cvd <- ir_cvd - c_cvd * paf_slice_cvd
      res_cvd <- ir_mean[
        sex == sex_ &
          (age < params$age_cutoff_cvd | age > params$age_cutoff_cvd_high),
        ir_cvd
      ]
      paf_other[
        year == y &
          sex == sex_ &
          age >= params$age_cutoff_cvd &
          age <= params$age_cutoff_cvd_high,
        cpaf_cvd := lb_cvd
      ]
      paf_other[
        year == y &
          sex == sex_ &
          (age < params$age_cutoff_cvd | age > params$age_cutoff_cvd_high),
        cpaf_cvd := res_cvd
      ]
      lb_cancer_cvd <- ir_c_cvd - c_c_cvd * paf_slice_cvd
      res_cancer_cvd <- ir_mean[
        sex == sex_ &
          (age < params$age_cutoff_cvd | age > params$age_cutoff_cvd_high),
        ir_cancer_cvd
      ]
      paf_other[
        year == y &
          sex == sex_ &
          age >= params$age_cutoff_cvd &
          age <= params$age_cutoff_cvd_high,
        cpaf_cancer_comorb := lb_cancer_cvd
      ]
      paf_other[
        year == y &
          sex == sex_ &
          (age < params$age_cutoff_cvd | age > params$age_cutoff_cvd_high),
        cpaf_cancer_comorb := res_cancer_cvd
      ]
      lb_cvd_cancer <- ir_cvd_cancer - c_cvd_cancer * paf_slice_c
      res_cvd_cancer <- ir_mean[
        sex == sex_ &
          (age < params$age_cutoff_cancer |
            age > params$age_cutoff_cancer_high),
        ir_cvd_cancer
      ]
      paf_other[
        year == y &
          sex == sex_ &
          age >= params$age_cutoff_cancer &
          age <= params$age_cutoff_cancer_high,
        cpaf_cvd_comorb := lb_cvd_cancer
      ]
      paf_other[
        year == y &
          sex == sex_ &
          (age < params$age_cutoff_cancer |
            age > params$age_cutoff_cancer_high),
        cpaf_cvd_comorb := res_cvd_cancer
      ]

      lb_comorb <- ir_c_o_cvd - c_c_o_cvd * paf_slice_c_o_cvd
      res_comorb <- ir_mean[
        sex == sex_ &
          (age < max(params$age_cutoff_cvd, params$age_cutoff_cancer) |
            age >
              max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high)),
        ir_comorb
      ]
      paf_other[
        year == y &
          sex == sex_ &
          age >= max(params$age_cutoff_cvd, params$age_cutoff_cancer) &
          age <= max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high),
        cpaf_comorb := lb_comorb
      ]
      paf_other[
        year == y &
          sex == sex_ &
          (age < max(params$age_cutoff_cvd, params$age_cutoff_cancer) |
            age >
              max(params$age_cutoff_cvd_high, params$age_cutoff_cancer_high)),
        cpaf_comorb := res_comorb
      ]
    }
  }
  melt(
    paf_other,
    id.vars = c("sex", "age", "year"),
    variable.name = "disease",
    value.name = "cpaf"
  )
}

constants_df <- calculate_constants(calibration_params, cpafs_per_year, ir_mean)
bm_paf_other <- bench::mark(
  original <- calculate_paf_other_original(
    base_input$list_of_years,
    calibration_params,
    constants_df,
    ir_mean,
    copy(cpafs_per_year)
  ),
  optimized <- calculate_paf_other(
    base_input$list_of_years,
    calibration_params,
    constants_df,
    ir_mean,
    copy(cpafs_per_year)
  ),
  relative = TRUE
)
print(bm_paf_other)

###########################################################################
# Benchmark calculate_costs
calculate_costs <- function(
  preval,
  costs,
  ind_costs,
  direct_costs,
  indirect_costs,
  dcost_total_base_year
) {
  # Helper: Direct calculation (replaces uniroot)
  calc_kpp <- function(kpp_1, prev, del_cost) {
    b <- (del_cost - sum(kpp_1 * prev, na.rm = TRUE)) / sum(prev, na.rm = TRUE)
    kpp_1 + b
  }

  # --- Direct Costs ---
  # Process both NCDs in one go
  ncds <- c("cancer", "cvd")
  out_kpp <- rbindlist(
    lapply(ncds, function(ncd_) {
      ncd_cols <- names(preval)[
        grepl(ncd_, names(preval)) | grepl("comorb", names(preval))
      ]
      tot_sick <- sum(preval[year == dcost_total_base_year, ..ncd_cols])

      # Melt and add sex column
      d_costs <- melt(costs[ncd == ncd_, -"ncd"], id.vars = "age")
      setnames(d_costs, "variable", "sex")
      d_costs[, "sex" := as.integer(sex)]

      # Include age in prev_dt
      prev_dt <- preval[
        year == dcost_total_base_year,
        c("age", "sex", ..ncd_cols)
      ]

      # Merge and calculate
      merged <- merge(d_costs, prev_dt, by = c("age", "sex"), all = TRUE)
      merged[, tot_prev := rowSums(.SD, na.rm = TRUE), .SDcols = ncd_cols]

      # Vectorized over sex
      merged[, tot_p := sum(tot_prev, na.rm = TRUE), by = sex]
      merged[, sum_c_m := sum(value, na.rm = TRUE), by = sex]
      merged[, del_cost := direct_costs[[ncd_]] * tot_p / tot_sick, by = sex]
      merged[, kpp_1 := del_cost * value / (tot_p * sum_c_m), by = sex]
      merged[, dcost_unit := calc_kpp(kpp_1, tot_prev, del_cost), by = sex]

      unique(
        merged[, .(age, sex, ncd = ncd_, dcost_unit)],
        by = c("age", "sex")
      )
    })
  )

  # --- Indirect Costs ---
  out_ik <- rbindlist(
    lapply(ncds, function(ncd_) {
      ncd_cols <- names(preval)[
        grepl(ncd_, names(preval)) | grepl("comorb", names(preval))
      ]
      tot_sick <- sum(preval[year == dcost_total_base_year, ..ncd_cols])

      # Melt and add sex column
      i_costs <- melt(ind_costs[ncd == ncd_, -"ncd"], id.vars = "age")
      setnames(i_costs, "variable", "sex")
      i_costs[, "sex" := as.integer(sex)]

      # Include age in prev_dt
      prev_dt <- preval[
        year == dcost_total_base_year,
        c("age", "sex", ..ncd_cols)
      ]

      # Merge and calculate
      merged <- merge(i_costs, prev_dt, by = c("age", "sex"), all = TRUE)
      merged[, tot_prev := rowSums(.SD, na.rm = TRUE), .SDcols = ncd_cols]

      merged[, tot_p := sum(tot_prev, na.rm = TRUE), by = sex]
      merged[,
        del_cost := (indirect_costs[[ncd_]] / tot_sick) * tot_p,
        by = sex
      ]
      merged[, icost_unit := value * del_cost / tot_prev]
      merged[is.na(icost_unit), icost_unit := 0]

      unique(
        merged[, .(age, sex, ncd = ncd_, icost_unit)],
        by = c("age", "sex")
      )
    })
  )

  return(list(out_kpp = out_kpp, out_ik = out_ik))
}

calculate_costs_original <- function(
  preval,
  costs,
  ind_costs,
  direct_costs,
  indirect_costs,
  dcost_total_base_year
) {
  # Calibrate costs

  # Direct costs, adjust KPP with a and b so that (a + b*KPP)*stock = total costs
  out_kpp <- data.table()
  out_ik <- data.table()

  # function to optimize
  f <- function(b) {
    ret <- sum((kpp_1 + b) * prev) - del_cost
    return(ret)
  }
  f1 <- Vectorize(f)

  for (ncd_ in c("cancer", "cvd")) {
    tot_sick <- sum(preval[
      year == dcost_total_base_year,
      ((names(preval) %like% ncd_) | (names(preval) %like% "comorb")),
      with = FALSE
    ])
    d_costs <- costs[ncd == ncd_]
    d_costs[, ncd := NULL]
    d_costs <- melt(d_costs, id.vars = "age")

    # Adjust indirect costs i_k with a constant c
    # so that c*i_k*stock = total indirect costs
    i_costs <- ind_costs[ncd == ncd_]
    i_costs[, ncd := NULL]
    i_costs <- melt(i_costs, id.vars = "age")

    out_kpp_0 <- CJ(age = 0:100, sex = 1:2, ncd = ncd_)
    # optimization step per sex
    for (s in 1:2) {
      c_m <- d_costs[variable == s, value]
      prev <- preval[
        year == dcost_total_base_year & sex == s,
        ((names(preval) %like% ncd_) | (names(preval) %like% "comorb")),
        with = FALSE
      ]
      tot_p <- sum(prev)
      del_cost <- direct_costs[[ncd_]] * tot_p / tot_sick
      kpp_1 <- del_cost * c_m / (tot_p * sum(c_m))
      rr <- uniroot(f1, c(-max(c_m) * 10, max(c_m) * 10))
      kpp_1 <- kpp_1 + rr$root
      out_kpp_0[sex == s, dcost_unit := kpp_1]
    }
    # indirect costs. Optimize so the unit_costs*NCDSim_stock = total_costs
    out_ik_0 <- CJ(age = 0:100, sex = 1:2, ncd = ncd_)
    for (s in 1:2) {
      c_m <- i_costs[variable == s, value]
      prev <- preval[
        year == dcost_total_base_year & sex == s,
        ((names(preval) %like% ncd_) | (names(preval) %like% "comorb")),
        with = FALSE
      ]
      tot_p <- sum(prev)
      del_cost <- (indirect_costs[[ncd_]] / tot_sick) * tot_p
      ik_i <- c_m * del_cost / rowSums(prev)
      ik_i[is.na(ik_i)] <- 0
      out_ik_0[sex == s, icost_unit := ik_i]
    }

    out_kpp <- rbind(out_kpp, out_kpp_0)
    out_ik <- rbind(out_ik, out_ik_0)
  }
  return(list(out_kpp = out_kpp, out_ik = out_ik))
}

bm_costs <- bench::mark(
  original <- calculate_costs_original(
    base_input$preval,
    base_input$costs,
    base_input$ind_costs,
    calibration_params$direct_costs,
    calibration_params$indirect_costs,
    calibration_params$dcost_total_base_year
  ),
  optimized <- calculate_costs(
    base_input$preval,
    base_input$costs,
    base_input$ind_costs,
    calibration_params$direct_costs,
    calibration_params$indirect_costs,
    calibration_params$dcost_total_base_year
  ),
  relative = TRUE
)

print(bm_costs)
###########################################################################
# Benchmark full calibrate implementations
baseline_parameters <- read_json(baseline_parameters_path)

bm_calibrate <- bench::mark(
  original <- {
    source("calibrate_original.R")
    # Read the input json file and assign parameters

    # Bundle all calibration-related parameters into a named list
    calibrate_cpaf_cancer <- as.integer(
      baseline_parameters$calibrate_cpaf_cancer
    )
    calibrate_cpaf_cvd <- as.integer(baseline_parameters$calibrate_cpaf_cvd)
    cpaf_cancer <- as.numeric(baseline_parameters$cpaf_cancer)
    cpaf_cvd <- as.numeric(baseline_parameters$cpaf_cvd)
    age_cutoff_cancer <- as.integer(baseline_parameters$age_cutoff_cancer)
    age_cutoff_cvd <- as.integer(baseline_parameters$age_cutoff_cvd)
    age_cutoff_cancer_high <- as.integer(
      baseline_parameters$age_cutoff_cancer_high
    )
    age_cutoff_cvd_high <- as.integer(baseline_parameters$age_cutoff_cvd_high)
    dcost_total_cancer <- as.numeric(baseline_parameters$dcost_total_cancer)
    icost_total_cancer <- as.numeric(baseline_parameters$icost_total_cancer)
    dcost_growth_cancer_input <- as.numeric(
      baseline_parameters$dcost_growth_cancer
    )
    dcost_total_cvd <- as.numeric(baseline_parameters$dcost_total_cvd)
    icost_total_cvd <- as.numeric(baseline_parameters$icost_total_cvd)
    dcost_growth_cvd_input <- as.numeric(baseline_parameters$dcost_growth_cvd)
    dcost_total_base_year <- as.integer(
      baseline_parameters$dcost_total_base_year
    )
    rr <- unlist(baseline_parameters$rr)
    rr_cancer_diet <- unlist(baseline_parameters$rr_cancer_diet)
    rr_cvd_diet <- unlist(baseline_parameters$rr_cvd_diet)
    communalities <- unlist(baseline_parameters$communalities)
    direct_costs <- c(
      cancer = baseline_parameters$dcost_total_cancer,
      cvd = baseline_parameters$dcost_total_cvd
    )
    indirect_costs <- c(
      cancer = baseline_parameters$icost_total_cancer,
      cvd = baseline_parameters$icost_total_cvd
    )

    calibration_results <- calibrate(
      input_path,
      calibrate_cpaf_cancer,
      calibrate_cpaf_cvd,
      prop_cancer,
      prop_cvd,
      age_cutoff_cancer,
      age_cutoff_cvd,
      age_cutoff_cancer_high,
      age_cutoff_cvd_high,
      rr,
      rr_cancer_diet,
      rr_cvd_diet,
      communalities,
      direct_costs,
      indirect_costs,
      dcost_total_base_year
    )
    calibration_results
  },
  refactored <- {
    source("calibrate.R")
    calibration_params <- list(
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
    calibration_results <- calibrate(input_path, calibration_params)
    calibration_results
  },
  relative = TRUE
)
print(bm_calibrate)
