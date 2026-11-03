library(bench)

calculate_death_and_mortality_rates <- function(dat_par) {
    r_dead <- dat_par$rate_dead[
        year ==
            fifelse(FINISH > max_year_mortality, max_year_mortality, FINISH) &
            sex == sex_ &
            age == age_,
        rate
    ]
    # In a few cases the observed death rates for children are zero. Since
    # this will result in a zero denominator when doing the proportional
    # adjustment of death rates below we replace the death rate with a very
    # small number
    if (r_dead == 0) {
        r_dead <- 1e-10
    }

    # All dead
    r_dead_all <- dat_par$rate_dead_all[
        year == FINISH &
            sex == sex_ &
            age == age_,
        rate
    ]
    if (r_dead_all == 0) r_dead_all <- 1e-10
}

calculate_death_and_mortality_rates_vectorized <- function(dat_par) {
    # Calculate all death and mortality rates vectorized
    target_year_mortality <- fifelse(
        FINISH > max_year_mortality,
        max_year_mortality,
        FINISH
    )
    # (used only for calibration purposes)
    death_and_mortality_rates <- CJ(
        sex = unique(dat$sex),
        age = unique(dat$age)
    )

    death_and_mortality_rates[
        dat_par$rate_dead[year == target_year_mortality],
        on = .(sex, age),
        r_dead := i.rate
    ]
    # In a few cases the observed death rates for children are zero. Since
    # this will result in a zero denominator when doing the proportional
    # adjustment of death rates below we replace the death rate with a very
    # small number
    death_and_mortality_rates[r_dead == 0, r_dead := 1e-10]

    # All dead
    death_and_mortality_rates[
        dat_par$rate_dead_all[year == target_year_mortality],
        on = .(sex, age),
        r_dead_all := i.rate
    ]
    death_and_mortality_rates[r_dead_all == 0, r_dead_all := 1e-10]

    # Calculate death rates for cancer and CVD from general death rates
    # and excess mortality
    # Cancer
    death_and_mortality_rates[
        dat_par$excess_mortality_cancer[year == target_year_mortality],
        on = .(sex, age),
        excess_mort_cancer := i.rate
    ]
    death_and_mortality_rates[, r_dead_cancer := r_dead * excess_mort_cancer]
    death_and_mortality_rates[r_dead_cancer > 0.999, r_dead_cancer := 0.999]
    # CVD
    death_and_mortality_rates[
        dat_par$excess_mortality_cvd[year == target_year_mortality],
        on = .(sex, age),
        excess_mort_cvd := i.rate
    ]
    death_and_mortality_rates[, r_dead_cvd := r_dead * excess_mort_cvd]
    death_and_mortality_rates[r_dead_cancer > 0.999, r_dead_cancer := 0.999]
    # Comorbidity
    death_and_mortality_rates[
        dat_par$excess_mortality_comorb[year == target_year_mortality],
        on = .(sex, age),
        excess_mort_comorb := i.rate
    ]
    death_and_mortality_rates[, r_dead_comorb := r_dead * excess_mort_comorb]
    death_and_mortality_rates[r_dead_comorb > 0.999, r_dead_comorb := 0.999]
    return(death_and_mortality_rates)
}
