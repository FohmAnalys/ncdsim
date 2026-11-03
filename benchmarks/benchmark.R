library(bench)

startyear <- 2020
endyear <- 2040
is_baseline <- TRUE

bm <- bench::mark(
  optimized <- {
    source("ncdsim.R")
    optimized <- simulate_model(
      startyear = startyear,
      endyear = endyear,
      is_baseline = is_baseline
    )
  },
  original <- {
    source("ncdsim_original.R")
    original <- simulate_model(
      startyear = startyear,
      endyear = endyear,
      is_baseline = is_baseline
    )
  },
  relative = FALSE
)
print(bm)
