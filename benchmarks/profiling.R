library(profvis)
source("ncdsim.R")
startyear <- 2020
endyear <- 2040

profile <- profvis(
    simulate_model(startyear = startyear, endyear = endyear, is_baseline = TRUE)
)

profile
