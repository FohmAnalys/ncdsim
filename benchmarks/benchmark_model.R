library(bench)
source("ncdsim.R")

model_original <- function(time, stocks, auxs) {
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

		p_pop_cancer_0_nm <- 1 -
			(1 - comm_smoking * paf_cancer_smoking) *
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
		paf_cvd_smoking <- paf(prev_smoking, rr_cvd_smoking)
		paf_cvd_inactivity <- paf(prev_inactivity, rr_cvd_inactivity)
		paf_cvd_bmi <- paf(prev_bmi, rr_cvd_bmi)
		paf_cvd_alcohol <- paf(prev_alcohol, rr_cvd_alcohol)

		p_pop_cvd_0_nm <- (1 -
			(1 - comm_smoking * paf_cvd_smoking) *
				(1 - comm_inactivity * paf_cvd_inactivity) *
				(1 - comm_bmi * paf_cvd_bmi) *
				(1 - comm_alcohol * paf_cvd_alcohol))
		p_pop_cvd_0 <- (1 - (1 - cpaf_diet_cvd) * (1 - p_pop_cvd_0_nm))

		p_pop_cvd <- p_pop_cvd_0 * calibration_cvd + paf_cvd_other

		f_pop_cvd <- s_pop * p_pop_cvd

		# Estimated healing time for cvd from an optimization
		# between observed and calculated over years 2009-2022
		f_cvd_pop <- s_cvd * heal_cvd_pop # change from (1/11).

		# Dying
		dr_cvd <- unname(pmin(1.0, r_dead_cvd))
		f_cvd_dead <- s_cvd * dr_cvd

		##** Calculating flows for Cancer and CVD # paf_comorb=max(paf_cvd, paf_ccr))
		# Getting sick

		p_pop_comorb_0 <- max(p_pop_cancer_0, p_pop_cvd_0)
		p_pop_comorb <- p_pop_comorb_0 * calibration_comorb + paf_comorb_other

		f_pop_comorb <- s_pop * p_pop_comorb
		f_cancer_comorb <- s_cancer *
			(p_pop_cvd_0 * calibration_cancer_comorb + paf_cancer_comorb_other)
		f_cvd_comorb <- s_cvd *
			(p_pop_cancer_0 * calibration_cvd_comorb + paf_cvd_comorb_other)

		# Estimated healing time for cancer and cvd

		f_comorb_pop <- s_comorb * heal_comorb_pop

		f_comorb_cancer <- s_comorb * heal_comorb_cancer

		f_comorb_cvd <- s_comorb * heal_comorb_cvd

		# Dying
		dr_comorb <- unname(pmin(1.0, r_dead_comorb))
		f_comorb_dead <- s_comorb * dr_comorb

		##
		## Calculating differential equations
		##

		# Equation for stock of non-NCD population
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

		# Equation for stock of population with cancer
		d_s_cancer_dt <- f_pop_cancer -
			f_cancer_pop +
			f_immig_cancer -
			f_emig_cancer -
			f_cancer_dead +
			f_comorb_cancer -
			f_cancer_comorb

		# Equation for stock of population with CVD
		d_s_cvd_dt <- f_pop_cvd -
			f_cvd_pop +
			f_immig_cvd -
			f_emig_cvd -
			f_cvd_dead +
			f_comorb_cvd -
			f_cvd_comorb

		# Equation for stock of population with CVD and cancer
		d_s_comorb_dt <- f_pop_comorb -
			f_comorb_pop +
			f_immig_comorb -
			f_emig_comorb -
			f_comorb_dead +
			f_cancer_comorb -
			f_comorb_cancer +
			f_cvd_comorb -
			f_comorb_cvd

		# Equation for stock of dead from non-NCD
		d_s_dead_dt <- f_dead

		# Equation for stock of dead from cancer
		d_s_dead_cancer_dt <- f_cancer_dead

		# Equation for stock of dead from CVD
		d_s_dead_cvd_dt <- f_cvd_dead

		# Equation for stock of dead from Comorbiditet
		d_s_dead_comorb_dt <- f_comorb_dead

		return(
			list(
				c(
					# Stock differentials
					d_s_pop_dt,
					d_s_cancer_dt,
					d_s_cvd_dt,
					d_s_comorb_dt,
					d_s_dead_dt,
					d_s_dead_cancer_dt,
					d_s_dead_cvd_dt,
					d_s_dead_comorb_dt
				),

				# Flows
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
			)
		)
	})
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
	calibration_cancer <- auxs[["calibration_cancer"]]
	paf_cancer_other <- auxs[["paf_cancer_other"]]
	heal_cancer_pop <- auxs[["heal_cancer_pop"]]
	r_dead_cancer <- auxs[["r_dead_cancer"]]

	rr_cvd_smoking <- auxs[["rr_cvd_smoking"]]
	rr_cvd_inactivity <- auxs[["rr_cvd_inactivity"]]
	rr_cvd_bmi <- auxs[["rr_cvd_bmi"]]
	rr_cvd_alcohol <- auxs[["rr_cvd_alcohol"]]

	cpaf_diet_cvd <- auxs[["cpaf_diet_cvd"]]
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
				(1 - comm_alcohol * paf_alcohol)

		p_pop_0 <- 1 - (1 - cpaf_diet) * (1 - p_pop_0_nm)
		p_pop <- p_pop_0 * calibration + paf_other

		return(list(
			paf_smoking = paf_smoking,
			paf_inactivity = paf_inactivity,
			paf_bmi = paf_bmi,
			paf_alcohol = paf_alcohol,
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
	p_pop_cancer_0 <- cancer_paf$p_pop_0
	p_pop_cancer <- cancer_paf$p_pop
	f_pop_cancer <- s_pop * p_pop_cancer
	f_cancer_pop <- s_cancer * heal_cancer_pop
	dr_cancer <- pmin(1.0, r_dead_cancer)
	f_cancer_dead <- s_cancer * dr_cancer

	# CVD flows
	cvd_paf <- calculate_paf(
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
	p_pop_cvd_0 <- cvd_paf$p_pop_0
	p_pop_cvd <- cvd_paf$p_pop
	f_pop_cvd <- s_pop * p_pop_cvd
	f_cvd_pop <- s_cvd * heal_cvd_pop
	dr_cvd <- pmin(1.0, r_dead_cvd)
	f_cvd_dead <- s_cvd * dr_cvd

	# Comorbidity flows
	p_pop_comorb_0 <- max(p_pop_cancer_0, p_pop_cvd_0)
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
			# Named elements (match the original exactly)
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
		)
	)
}


# Generate mock data for stocks
set.seed(123) # For reproducibility
stocks <- c(
	s_pop = 1e6,
	s_cancer = 5e4,
	s_cvd = 1e5,
	s_comorb = 2e4,
	s_dead = 1e4,
	s_dead_cancer = 5e3,
	s_dead_cvd = 8e3,
	s_dead_comorb = 2e3
)

# Generate mock data for auxs
auxs <- list(
	# Demographic parameters
	r_dead = 0.01,
	n_born = 1e4,
	r_immig = 0.005,
	r_emig = 0.003,

	# Prevalence of risk factors
	prev_smoking = 0.2,
	prev_inactivity = 0.3,
	prev_bmi = 0.25,
	prev_alcohol = 0.15,

	# Relative risks for cancer
	rr_cancer_smoking = 2.5,
	rr_cancer_inactivity = 1.8,
	rr_cancer_bmi = 1.5,
	rr_cancer_alcohol = 1.2,

	# Relative risks for CVD
	rr_cvd_smoking = 2.0,
	rr_cvd_inactivity = 1.5,
	rr_cvd_bmi = 1.3,
	rr_cvd_alcohol = 1.1,

	# Community factors
	comm_smoking = 0.8,
	comm_inactivity = 0.7,
	comm_bmi = 0.6,
	comm_alcohol = 0.5,

	# Counterfactual PAFs
	cpaf_diet_cancer = 0.1,
	cpaf_diet_cvd = 0.08,

	# Calibration factors
	calibration_cancer = 1.0,
	calibration_cvd = 1.0,
	calibration_comorb = 1.0,
	calibration_cancer_comorb = 0.9,
	calibration_cvd_comorb = 0.85,

	# Other PAFs
	paf_cancer_other = 0.05,
	paf_cvd_other = 0.04,
	paf_comorb_other = 0.03,
	paf_cancer_comorb_other = 0.02,
	paf_cvd_comorb_other = 0.015,

	# Healing rates
	heal_cancer_pop = 0.2,
	heal_cvd_pop = 0.1,
	heal_comorb_pop = 0.05,
	heal_comorb_cancer = 0.03,
	heal_comorb_cvd = 0.02,

	# Disease-specific death rates
	r_dead_cancer = 0.05,
	r_dead_cvd = 0.04,
	r_dead_comorb = 0.08
)

# Mock time value
time <- seq(0, 1, 0.05)

benchmark_results <- bench::mark(
	original <- model_original(time, stocks, auxs),
	optimized <- model(time, stocks, auxs),
	relative = TRUE
)

print(benchmark_results)
