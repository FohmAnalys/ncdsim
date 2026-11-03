library(data.table)
library(ggplot2)
library(gridExtra)
library(scales)
library(stringr)
library(jsonlite)

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


calculate_heal_rates <- function(heal, riskgroup, base_year) {
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
  
  setnafill(heal,
            "const",
            fill = 0,
            cols = c("hr_cvd", "hr_cancer", "hr_comorb_cancer", "hr_comorb_cvd", "hr_comorb_healthy"))
  
  heal_rate <- heal[, {
    fit <- lm(cbind(hr_cvd, hr_cancer, hr_comorb_cancer, hr_comorb_cvd, hr_comorb_healthy) ~ year)
    pred <- predict(fit, newdata = list(year = base_year))
    as.list(pred)
  }, by = .(sex, age)]
  setnames(heal_rate, old = paste0("V", 1:5), new = c("cvd", "cancer", "comorb_cancer", "comorb_cvd", "comorb_healthy"))

  return(heal_rate)
}

read_prevalences <- function(input_path, list_of_years) {
  # read prevalences for risk factors and relative risks for food
  food_prevalences <- fread(paste0(input_path, "prev_diet.csv"), sep = ";")
  prev_smoking <- fread(paste0(input_path, "prev_smoking.csv"), sep = ";")
  prev_bmi <- fread(paste0(input_path, "prev_bmi.csv"), sep = ";")
  prev_alcohol <- fread(paste0(input_path, "prev_alcohol.csv"), sep = ";")
  prev_inactivity <- fread(paste0(input_path, "prev_inactivity.csv"), sep = ";")

  filter_cols <- c(c("sex", "age"), as.character(list_of_years))
  prev_bmi <- prev_bmi[, ..filter_cols]
  prev_alcohol <- prev_alcohol[, ..filter_cols]
  prev_inactivity <- prev_inactivity[, ..filter_cols]
  prev_smoking <- prev_smoking[, ..filter_cols]
  return(
    list(
      food_prevalences = food_prevalences,
      prev_alcohol = prev_alcohol,
      prev_inactivity = prev_inactivity,
      prev_smoking = prev_smoking,
      prev_bmi = prev_bmi
    )
  )
}

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

calculate_cpafs_per_year <- function(
  list_of_years,
  communalities,
  cpaf_diet,
  prevalences,
  rr,
  base_year
) {
  base_year <- as.character(base_year)
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
  
  # ---- 4b. Add UV from paf
  paf_uv <- CJ(sex = 1:2,
               age = 0:100,
               year = list_of_years,
               ncd = c("cancer", "cvd"),
               factor = "uv",
               cpaf = as.numeric(NA)
  )

  paf_uv[ncd=="cancer", cpaf := rr["paf_cancer_uv"]]
  paf_uv[ncd=="cvd", cpaf := rr["paf_cvd_uv"]]
  

  # --- 5. Combine all PAFs (diet + lifestyle) ---
  all_pafs <- rbindlist(list(
    cpaf_diet_year[, .(sex, age, year, ncd, cpaf, factor)],
    prev_long_all[, .(sex, age, year, ncd, cpaf, factor)],
    paf_uv[, .(sex, age, year, ncd, cpaf, factor)]
  ))
  

  # --- 6. Calculate combined cpaf for each group ---
  cpaf_combined <- all_pafs[,
    .(cpaf = 1 - prod(1 - cpaf, na.rm = TRUE)),
    by = .(sex, age, year, ncd)
  ]
  #cpaf_combined[, year := as.integer(year)]
  
  # --- 7. Split into paf_cancer, paf_cvd, and calculate means ---
  paf_cancer <- cpaf_combined[ncd == "cancer", .(sex, age, year, cpaf)]
  paf_cvd <- cpaf_combined[ncd == "cvd", .(sex, age, year, cpaf)]

  paf_comorb <- merge(
    paf_cancer[, .(sex, age, year, cpaf_cancer = cpaf)],
    paf_cvd[, .(sex, age, year, cpaf_cvd = cpaf)],
    by = c("sex", "age", "year")
#  )[, cpaf := pmax(cpaf_cancer, cpaf_cvd)]
  )[, cpaf := pmin(cpaf_cancer, cpaf_cvd)]
  paf_comorb <- paf_comorb[, .(sex, age, year, cpaf)]
  
  cpaf_cvd_mean <- paf_cvd[, {
    fit <- lm(cpaf ~ year)
    pred <- predict(fit, newdata = list(year = base_year))
    as.list(pred)
  }, by = .(sex, age)]
  setnames(cpaf_cvd_mean, old= "1", new='cpaf')
  
  cpaf_cancer_mean <- paf_cancer[, {
    fit <- lm(cpaf ~ year)
    pred <- predict(fit, newdata = list(year = base_year))
    as.list(pred)
  }, by = .(sex, age)]
  setnames(cpaf_cancer_mean, old= "1", new='cpaf')
  
  cpaf_cancer_cvd_mean <- paf_comorb[, {
    fit <- lm(cpaf ~ year)
    pred <- predict(fit, newdata = list(year = base_year))
    as.list(pred)
  }, by = .(sex, age)]
  setnames(cpaf_cancer_cvd_mean, old= "1", new='cpaf')

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
  setkeyv(paf_other, c("sex", "age", "year"))

  return(paf_other)
}

calculate_costs <- function(
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


calibrate <- function(input_path, params, base_year) {
  # Read base input
  base_input <- read_base_input(input_path)
  list_of_years <- base_input$list_of_years
  # Calculate riskgroup
  riskgroup <- calculate_riskgroup(base_input$preval)
  # Calculate incidence rate for flows
  incid_rate <- calculate_incidence_rates(base_input$incid, riskgroup)
  # Calculate healing rates
  heal_rate <- calculate_heal_rates(base_input$heal, riskgroup, base_year)
  # Read prevalences
  prevalences <- read_prevalences(input_path, base_input$list_of_years)
  food_categories <- colnames(prevalences$food_prevalences)[
    3:length(prevalences$food_prevalences)
  ]
  # calculate cPAF for food
  cpaf_diet <- calculate_cpaf_diet(
    prevalences$food_prevalences,
    params$rr_cancer_diet,
    params$rr_cvd_diet
  )
  cpafs_per_year <- calculate_cpafs_per_year(
    base_input$list_of_years,
    params$communalities,
    cpaf_diet,
    prevalences,
    params$rr,
    base_year
  )
  
  ir_mean <- incid_rate[, {
    fit <- lm(cbind(ir_cvd, ir_cancer, ir_cvd_cancer, ir_cancer_cvd, ir_comorb) ~ year)
    pred <- predict(fit, newdata = list(year = base_year))
    as.list(pred)
  }, by = .(sex, age)]
  setnames(ir_mean, old = paste0("V", 1:5), new = colnames(incid_rate)[4:length(colnames(incid_rate))])

  constants_df <- calculate_constants(params, cpafs_per_year, ir_mean)
  paf_other <- calculate_paf_other(
    base_input$list_of_years,
    params,
    constants_df,
    ir_mean,
    cpafs_per_year
  )
  
  paf_other <- setnafill(paf_other, "const", fill=0.0, cols = c(
    "cpaf_cancer",
    "cpaf_cvd",
    "cpaf_cancer_comorb",
    "cpaf_cvd_comorb" ,
    "cpaf_comorb" 
  ))

  # Calibrate costs
  costs <- calculate_costs(
    base_input$preval,
    base_input$costs,
    base_input$ind_costs,
    params$direct_costs,
    params$indirect_costs,
    params$dcost_total_base_year
  )

  return(list(
    constants = melt(
      constants_df,
      id.vars = c("sex", "age"),
      variable.name = 'disease',
      value.name = 'constant'
    ),
    heal_rate = melt(
      heal_rate,
      id.vars = c("sex", "age"),
      variable.name = 'disease',
      value.name = 'heal_rate'
    ),
    paf_other = melt(
      paf_other,
      id.vars = c("sex", "age", "year"),
      variable.name = 'disease',
      value.name = 'cpaf'
    ),
    dcost_unit = costs$out_kpp,
    icost_unit = costs$out_ik
  ))
}
