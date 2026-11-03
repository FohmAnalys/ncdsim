library(bench)
library(data.table)

calculate_stock_vectors <- function(last_year_stocks) {
	stock_cols <- setdiff(names(last_year_stocks), c("sex", "age"))
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

	return(this_year)
}

calculate_stock_vectors_new <- function(dat) {
	START <- 2015
	last_year <- dat[
		.(START),
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
	stock <- calculate_stock_vectors(last_year)
	setkeyv(stock, c("sex", "age"))
	return(stock)
}

calculate_stock_vectors_loop <- function(dat) {
	START <- 2015
	this_year <- CJ(sex = 1:2, age = 0:100)
	for (sex_ in unique(dat$sex)) {
		for (age_ in unique(dat$age)) {
			if (age_ == 0) {
				yvec <- c(
					s_pop = 0.0,
					s_cancer = 0.0,
					s_cvd = 0.0,
					s_comorb = 0.0,
					s_dead = 0.0,
					s_dead_cancer = 0.0,
					s_dead_cvd = 0.0,
					s_dead_comorb = 0.0
				)
			} else if (age_ == 100) {
				# Age 100+ is an open age group which leads to an accumulation of
				#       counts
				yvec <- unlist(dat[
					.(START, sex_, age_ - 1),
					.(
						s_pop,
						s_cancer,
						s_cvd,
						s_comorb,
						s_dead,
						s_dead_cancer,
						s_dead_cvd,
						s_dead_comorb
					)
				]) +
					unlist(dat[
						.(START, sex_, age_),
						.(
							s_pop,
							s_cancer,
							s_cvd,
							s_comorb,
							s_dead,
							s_dead_cancer,
							s_dead_cvd,
							s_dead_comorb
						)
					])
			} else {
				yvec <- unlist(dat[
					.(START, sex_, age_ - 1),
					.(
						s_pop,
						s_cancer,
						s_cvd,
						s_comorb,
						s_dead,
						s_dead_cancer,
						s_dead_cvd,
						s_dead_comorb
					)
				])
			}
			this_year[
				.(sex_, age_),
				c(
					"s_pop",
					"s_cancer",
					"s_cvd",
					"s_comorb",
					"s_dead",
					"s_dead_cancer",
					"s_dead_cvd",
					"s_dead_comorb"
				) := as.list(yvec)
			]
		}
	}

	return(this_year)
}

dat <- CJ(year = 2015:2030, sex = 1:2, age = 0:100)
stock_cols <- c("s_pop", "s_cancer", "s_cvd", "s_comorb")
max_values = c(6000, 200, 1000, 200)
for (i in seq_along(stock_cols)) {
	dat[, (stock_cols[i]) := sample(0:max_values[i], .N, replace = TRUE)]
}
dat[, c("s_dead", "s_dead_cancer", "s_dead_cvd", "s_dead_comorb") := 0]

bm <- bench::mark(
	original = calculate_stock_vectors_loop(dat),
	optimized = calculate_stock_vectors_new(dat),
	min_iterations = 100
)

print(bm)
