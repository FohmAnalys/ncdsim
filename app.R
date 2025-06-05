################################################################################
## UI for running NCD-SIM simulations and visualizing the results
##
## The main model script NCDSim_R.R is invoked for new simulations. In order not
## to freeze the UI while waiting for a simulation, this is done in a background
## R process using futureCall() from the future package. When run from the UI,
## NCDSim_R.R outputs data year-by-year, so that the UI can load and visualize 
## the data continuously. 
##
## Active simulations can be stopped from the UI, by the UI writing to a status 
## file on disk, ui_status.txt, which is checked within each simulation loop by 
## NCDSim_R.R. This is done for each scenario, in separate sub-directories in
## which the simulation data is temporarily stored during a run as well. The
## directory for the baseline scenario is denoted scen_1.
##
## Note that the reading of the simulation output .rda-files into the UI can
## glitch a bit, possibly because R tries to read an .rda file which is being
## written to. This will print an error message to the console, but will 
## normally work fine from there, as the file will be read on the next attempt
## instead.
##
## The following abbreviations are used in the code, referring to the three 
## different plot types used:
## 
## ts = time series
## cs = cross-section
## acc = acccumulated
##
## There is an English and a Swedish version of the GUI. It must be hard-coded
## in this script (under Setup below) by specifying lang <- "swe" / "eng
##
################################################################################

################################################################################
## Setup
################################################################################

library(data.table)
library(future)
library(shiny)
library(shinyjs)
library(shinyFiles)
library(shinyWidgets)
library(shinyBS)
library(DT)
library(ggplot2)
library(plotly)

## Initiate multisession plan, enabling model to be run in separate session
plan(multisession)

PROJECTROOT <- getwd()

## Toggle language in the UI. It can only be set here, not in the GUI.
# lang <- "eng"
lang <- "swe"

## Load model code
source(paste0(PROJECTROOT, "/ncdsim.R"))

################################################################################
## Graphics and formatting settings and help functions
################################################################################

## Load and initialize ggplot fohm template
source(paste0(PROJECTROOT, "/theme_fohm.R"))
theme_set(theme_fohm(text_size = 8))

## Help functions for pretty text formatting of numbers
pretty_print <- function(x, digits = 3) {
  formatC(x, digits = digits, big.mark = " ", format = "f", 
          decimal.mark=switch(lang, "eng" = ".", "swe" = ",")
  )
}

pretty_print_1 <- function(x) {
  pretty_print(x, digits = 1)
}

int_print <- function(x) pretty_print(x, digits=0)

year_print <- function(x) formatC(x, digits = 0, format = "f")

## Help function for scaling line widths in plotly
## Note: should exclude non-line objects?
scale_line_width <- function(plotly_obj, scale_factor = 0.75) {
  n_layers <- length(plotly_obj$x$data)
  
  for (j in 1:n_layers) {
    if (substr(plotly_obj$x$data[[j]]$mode, 1, 4) == "line")  {
      plotly_obj$x$data[[j]]$line$width <- 
        plotly_obj$x$data[[j]]$line$width * scale_factor
    }
  } # for
  plotly_obj
} # function scale_line_width

## Help function for adding legends in plotly
change_legend <- function(plotly_obj, labels=NULL, title) {
  
  n_layers <- length(plotly_obj$x$data)
  
  if (!is.null(labels)) {
    for (j in 1:n_layers) {
      plotly_obj$x$data[[j]]$name <- labels[j]
    } # for
  } # if
  
  plotly_obj$x$layout$legend$title$text <- title
  
  plotly_obj
} # function change_legend

## Convenience function for creating an info icon with hover text
info_icon <- function(x) {
  icon("info-sign", lib = "glyphicon", style = "color:#009F80;", title = x)
}

ok_icon <- function(x) {
  icon("ok", lib = "glyphicon", style = "color:#009F80;", title = x)
}

## Custom slider input for risk factors
risk_factor_slider <- function(inputId, label, info_text) {
  sliderInput(
    inputId = inputId,
    label = span(label, info_icon(info_text)),
    min = 0,
    max = 2.0,
    value = 1.0,
    step = 0.01)
}

################################################################################
## Pre-definitions of various labels and variable lists
################################################################################

risk_factor_info_text <- switch(lang, 
                                "eng" = "Click to expand and use sliders to adjust the size of the risk groups with an adjustment factor ranging from 0 (turned off) to 2 (doubled). The adjustments are phased in linearly with the 'Phase-in period' slider below.", 
                                "swe"= "Klicka för att expandera och använd reglagen för att justera storleken på riskgrupperna med en justeringsfaktor som från 0 (avslagen) till 2 (dubblerad). Justeringarna fasas in linjärt med reglaget för 'Infasningsperiod' nedan.")

label_slider_smoking <- switch(lang, "eng" = "Smoking", "swe" = "Rökning")
info_text_smoking <- switch(lang, 
                            "eng" = "Daily smokers",
                            "swe" = "Dagligrökare")

label_slider_alcohol <- switch(lang, "eng" = "Alcohol", "swe" = "Alkohol")
info_text_alcohol <- switch(lang, 
                        "eng" = "Audit-c score >= 6 for men and >= 5 for women",
                        "swe" = "Audit-c >= 6 för män och >= 5 för kvinnor")

label_slider_inactivity <- switch(lang,
                                  "eng" = "Insufficient physical activity", 
                                  "swe" = "Otillräcklig fysisk aktivitet")

info_text_inactivity <- switch(lang, 
                    "eng" = "< 150 minutes moderate physical activity per week",
                    "swe" = "< 150 minuter måttlig fysisk aktivitet per vecka")

label_slider_obesity <- switch(lang, "eng" = "Obesity", "swe" = "Obesitas")

info_text_obesity <- switch(lang, 
                            "eng" = "BMI > 30",
                            "swe" = "BMI > 30")

label_slider_fruit <- switch(lang, "eng" = "Fruit", "swe" = "Frukt")
info_text_fruit <- switch(lang, 
                          "eng" = "Eats < 200 g fruit per day",
                          "swe" = "Äter < 200 g frukt per dag ")

label_slider_wholegrains <- switch(lang, 
                                   "eng" = "Wholegrains",
                                   "swe" = "Fullkorn")
info_text_wholegrains <- switch(lang, 
                                "eng" = "Eats < 100 g wholegrains per day",
                                "swe" = "Äter < 100 g wholegrains per dag ")

label_slider_greens <- switch(lang, "eng" = "Vegetables", "swe" = "Grönsaker")
info_text_greens <- switch(lang, 
                           "eng" = "Eats < 200 g vegetables per day",
                           "swe" = "Äter < 200 g grönsaker per dag ")

label_slider_meat <- switch(lang, "eng" = "Meat", "swe" = "Kött")
info_text_meat <- switch(lang, 
                         "eng" = "Eats > 50 g red and processed meat per day",
                         "swe" = "Äter > 50 g rött och processat kött per dag ")

label_slider_salt <- switch(lang, "eng" = "Salt", "swe" = "Salt")
info_text_salt <- switch(lang, 
                         "eng" = "Eats > 6 g salt per day",
                         "swe" = "Äter > 6 g salt per dag ")

baseline_slider_names <-
  c("slider_simyears",
    "slider_intervention_years_1",
    "cfact_smoking_1",
    "cfact_alcohol_1",
    "cfact_inactivity_1",
    "cfact_obesity_1",
    "cfact_fruit_1",
    "cfact_wholegrains_1",
    "cfact_greens_1",
    "cfact_meat_1",
    "cfact_salt_1",
    "cpaf_cancer",
    "cpaf_cvd",
    "age_cutoff_cancer",
    "age_cutoff_cvd")

scenario_slider_names <- 
  c("slider_intervention_years_",
    "cfact_smoking_",
    "cfact_alcohol_",
    "cfact_inactivity_",
    "cfact_obesity_",
    "cfact_fruit_",
    "cfact_wholegrains_",
    "cfact_greens_",
    "cfact_meat_",
    "cfact_salt_")

ts_group_list <- switch(lang,
                        "eng" = 
                          list("Total" = "total", 
                               "Sex" = "sex", 
                               "Age group" = c(
                                 "70+ years" = "age_70p",
                                 "Ten-year groups" = "age",
                                 "16-84 years" = "age_16_84")),
                        "swe" =
                          list("Total" = "total", 
                               "Kön" = "sex", 
                               "Åldersgrupp" = c(
                                 "70+ år" = "age_70p",
                                 "Tioårsgrupper" = "age",
                                 "16-84 år" = "age_16_84"))
) # end switch

ts_measure_list <- switch(lang,
  "eng" = 
    list("Absolute" = "abs",
         "Population share" = "pop_share",
         "Difference vz baseline" = c(
           "Absolute difference" = "abs_diff",
           "Relative difference" = "rel_diff",
           "Relative difference of population shares" = "rel_diff_pop_share")),
  "swe" = 
    list("Absoluta tal" = "abs",
         "Befolkningsandel" = "pop_share",
         "Differens gentemot baslinje " = c(
           "Absolut differens" = "abs_diff",
           "Relativ differens" = "rel_diff",
           "Relativ differens av befolkningsandelar" = "rel_diff_pop_share"))
) # end switch

cs_group_list <- switch(lang,
                        "eng" = 
                          list("Total" = "total", 
                               "Sex" = "sex"),
                        "swe" = 
                          list("Total" = "total", 
                               "Kön" = "sex")
) # end switch

cs_measure_list <- ts_measure_list

acc_group_list <- switch(lang,
                         "eng" = 
                           list("Total" = "total", 
                                "Sex" = "sex", 
                                "Age group" = c(
                                  "70+ years" = "age_70p",
                                  "16-84 years" = "age_16_84")),
                         "swe"= list("Total" = "total", 
                                     "Kön" = "sex",
                                     "Åldersgrupp" = c(
                                       "70+ år" = "age_70p",
                                       "16-84 år" = "age_16_84"))
) # end switch

acc_measure_list <- switch(lang,
                           "eng" = 
                             list("Absolute" = "abs", 
                                  "Difference vz baseline" = "abs_diff"),
                           "swe" = 
                             list("Absoluta tal" = "abs", 
                                  "Differens gentemot baslinje" = "abs_diff")
) # end switch

################################################################################

## Read default baseline parameter values
default_parameters_path <- paste0(PROJECTROOT, "/Input/defaults/baseline.json")
default_parameters <- read_json(default_parameters_path)

# browser()

## Min and max years for simulation
start_year_min <- 2020
end_year_max <- 2100

## Max number of scenarios (including baseline)
n_scen_max <- 8

## Create scenarios folder if not existing already
if (file.exists("UI/scenarios") == FALSE) dir.create("UI/scenarios")

## Create folders for all simulation runs
for (scen in 1:n_scen_max) {
  tmp_path <- paste0(PROJECTROOT, "/UI/runs/scen_", scen)
  if (file.exists(tmp_path) == FALSE) {
    dir.create(tmp_path)
  }
  if (file.exists(paste0(tmp_path, "/data")) == FALSE) {
    dir.create(paste0(tmp_path, "/data"))
  }
}

## Reset log files
for (scen in 1:n_scen_max) {
  cat(switch(lang,
             "eng" = "Waiting to start or load", 
             "swe" = "Väntar på att starta eller ladda"),
      file = paste0(PROJECTROOT, "/UI/runs/scen_", scen, "/log.txt"))
}

slider_width <- "100%"

## All variable names (except the ones derived in the UI), used for validation
all_var_names <- c(
  "year",
  "sex",
  "age",
  "s_pop",
  "f_dead",
  "f_born",
  "f_immig_pop",
  "f_immig_cvd",
  "f_immig_cancer",
  "f_cancer_pop",
  "f_pop_cvd",
  "f_cvd_pop",
  "f_cancer_dead",
  "f_cvd_dead",
  "dcost_unit_cancer",
  "icost_unit_cancer",
  "dcost_unit_cvd",
  "icost_unit_cvd",
  "dcost_cancer", 
  "dcost_cvd",
  "icost_cancer",
  "icost_cvd",
  "dcost_growth_cancer",
  "dcost_growth_cvd",
  "prev_fruit",
  "prev_wholegrains",
  "prev_greens",
  "prev_meat",
  "prev_salt",
  "paf_cancer_fruit",
  "paf_cancer_wholegrains",
  "paf_cancer_greens",
  "paf_cancer_meat",
  "paf_cancer_salt",
  "paf_cvd_fruit",
  "paf_cvd_wholegrains",
  "paf_cvd_greens",
  "paf_cvd_meat",
  "paf_cvd_salt",
  "s_cancer",
  "s_cvd",
  "s_dead_cancer",
  "s_dead_cvd",
  "s_dead",
  "f_emig_pop",
  "f_emig_cvd",
  "f_emig_cancer",
  "f_pop_cancer",
  "dr",
  "dr_cancer",
  "dr_cvd",
  "prev_alcohol",
  "prev_obesity",
  "prev_smoking",
  "prev_inactivity",
  "paf_cvd_other",
  "paf_cvd_smoking",
  "paf_cvd_inactivity",
  "paf_cvd_obesity",
  "paf_cvd_alcohol",
  "cpaf_diet_cvd",
  "p_pop_cvd",
  "paf_cancer_other",
  "paf_cancer_smoking",
  "paf_cancer_inactivity",
  "paf_cancer_obesity",
  "paf_cancer_alcohol",
  "cpaf_diet_cancer",
  "p_pop_cancer",
  "calibration_cancer",
  "calibration_cvd",
  "calibration_mortality"
)

## Baseline variables inherited by alternative scenarios, for validation
baseline_vars_to_check <- c(
  "paf_cancer_other", 
  "paf_cvd_other", 
  "icost_unit_cancer",
  "dcost_unit_cancer",
  "icost_unit_cvd",
  "dcost_unit_cvd",
  "calibration_cancer", 
  "calibration_cvd",
  "calibration_mortality"
)

## All parameter names
all_par_names <- c(
  "cpaf_cancer",
  "cpaf_cvd",
  "age_cutoff_cancer",
  "age_cutoff_cvd",
  "dcost_total_cancer",
  "icost_total_cancer",
  "dcost_growth_cancer",
  "dcost_total_cvd",
  "icost_total_cvd",
  "dcost_growth_cvd",
  "dcost_total_base_year",
  "rr",
  "rr_cancer_diet",
  "rr_cvd_diet",
  "communalities")

## Variables to plot (and keys) in time series
var_list  <- switch(lang,
                    "eng" =  
                      list("Cancer and CVD combined" = c(
                        "Ongoing cases" = "s_ncd", 
                        "New cases" = "f_ncd", 
                        "Deaths" = "f_ncd_dead",
                        "Healthcare costs" = "dcost_ncd_msek",
                        "Indirect costs" = "icost_ncd_msek"),
                        "Non-diet risk factors" = c(
                          "Alcohol" = "s_alcohol",
                          "Obesity" = "s_obesity",
                          "Smoking" = "s_smoking",
                          "Insufficient physical inactivity" = "s_inactivity"),
                        "Dietary risk factors" = c(
                          "Fruit" = "s_fruit",
                          "Wholegrains" = "s_wholegrains",
                          "Vegetables" = "s_greens",
                          "Meat" = "s_meat",
                          "Salt" = "s_salt"),
                        "Demographics" = c(
                          "Population size" = "total_pop",
                          "Births" = "f_born",
                          "Deaths" = "f_total_dead"),
                        "Cancer" = c(
                          "Ongoing cases" = "s_cancer",
                          "New cases" = "f_cancer",
                          "Deaths" = "f_cancer_dead",
                          "Healthcare costs" = "dcost_cancer_msek",
                          "Indirect costs" = "icost_cancer_msek"),
                        "CVD" = c(
                          "Ongoing cases" = "s_cvd",
                          "New cases" = "f_cvd",
                          "Deaths" = "f_cvd_dead",
                          "Healthcare costs" = "dcost_cvd_msek",
                          "Indirect costs" = "icost_cvd_msek")
                      ),
                    "swe" =
                      list("Cancer plus hjärt-kärl" = c(
                        "Pågående fall" = "s_ncd",
                        "Nya fall" = "f_ncd",
                        "Dödsfall" = "f_ncd_dead",
                        "Sjukvårdskostnader" = "dcost_ncd_msek",
                        "Indirekta kostnader" = "icost_ncd_msek"),
                        "Riskfaktorer, icke-kost" = c(
                          "Alkohol" = "s_alcohol",
                          "Obesitas" = "s_obesity",
                          "Rökning" = "s_smoking",
                          "Otillräcklig fysisk aktivitet" = "s_inactivity"),
                        "Riskfaktorer, kost" = c(
                          "Frukt" = "s_fruit",
                          "Fullkorn" = "s_wholegrains",
                          "Grönsaker" = "s_greens",
                          "Kött" = "s_meat",
                          "Salt" = "s_salt"),
                        "Demografi" = c(
                          "Befolkningsstorlek" = "total_pop",
                          "Födslar" = "f_born",
                          "Dödsfall" = "f_total_dead"),
                        "Cancer" = c(
                          "Pågående fall" = "s_cancer",
                          "Nya fall" = "f_cancer",
                          "Dödsfall" = "f_cancer_dead",
                          "Sjukvårdskostnader" = "dcost_cancer_msek",
                          "Indirekta kostnader" = "icost_cancer_msek"),
                        "Hjärt-kärl" = c(
                          "Pågående fall" = "s_cvd",
                          "Nya fall" = "f_cvd",
                          "Dödsfall" = "f_cvd_dead",
                          "Sjukvårdskostnader" = "dcost_cvd_msek",
                          "Indirekta kostnader" = "icost_cvd_msek")
                      )
) # end switch

## Labels to be used on y-axis for time series graphs
var_labels  <- switch(lang,
                      "eng" = 
                        list(
                          "f_ncd" = "New cases",
                          "s_ncd" = "Ongoing cases",
                          "f_ncd_dead" = "NCD deaths",
                          "dcost_ncd_msek" = "Costs (million SEK)",
                          "icost_ncd_msek" = "Costs (million SEK)",
                          "s_alcohol" = "Persons in risk group",
                          "s_obesity" = "Persons in risk group",
                          "s_smoking" = "Persons in risk group",
                          "s_inactivity" = "Persons in risk group",
                          "s_fruit" = "Persons in risk group",
                          "s_wholegrains" = "Persons in risk group",
                          "s_greens" = "Persons in risk group",
                          "s_meat" = "Persons in risk group",
                          "s_salt" = "Persons in risk group",
                          "total_pop" = "Population size",
                          "f_born" = "Births",
                          "f_total_dead" = "Deaths (all causes)",
                          "s_cancer" = "Ongoing cases",
                          "f_cancer" = "New cases",
                          "f_cancer_dead" = "Cancer deaths",
                          "dcost_cancer_msek" = "Costs (million SEK)",
                          "icost_cancer_msek" = "Costs (million SEK)",
                          "s_cvd" = "Ongoing cases",
                          "f_cvd" = "New cases",
                          "f_cvd_dead" = "CVD deaths",
                          "dcost_cvd_msek" = "Costs (million SEK)",
                          "icost_cvd_msek" = "Costs (million SEK)"
                        ),
                      "swe" = 
                        list(
                          "f_ncd" = "Nya fall",
                          "s_ncd" = "Pågående fall",
                          "f_ncd_dead" = "Dödsfall",
                          "dcost_ncd_msek" = "Kostnader (MSEK)",
                          "icost_ncd_msek" = "Kostnader (MSEK)",
                          "s_alcohol" = "Personer i riskgrupp",
                          "s_obesity" = "Personer i riskgrupp",
                          "s_smoking" = "Personer i riskgrupp",
                          "s_inactivity" = "Personer i riskgrupp",
                          "s_fruit" = "Personer i riskgrupp",
                          "s_wholegrains" = "Personer i riskgrupp",
                          "s_greens" = "Personer i riskgrupp",
                          "s_meat" = "Personer i riskgrupp",
                          "s_salt" = "Personer i riskgrupp",
                          "total_pop" = "Befolkningsstorlek",
                          "f_born" = "Födslar",
                          "f_total_dead" = "Dödsfall (alla orsaker)",
                          "f_cancer" = "Nya fall",
                          "s_cancer" = "Pågående fall",
                          "f_cancer_dead" = "Dödsfall",
                          "dcost_cancer_msek" = "Kostnader (MSEK)",
                          "icost_cancer_msek" = "Kostnader (MSEK)",
                          "f_cvd" = "Nya fall",
                          "s_cvd" = "Pågående fall",
                          "f_cvd_dead" = "Dödsfall",
                          "dcost_cvd_msek" = "Kostnader (MSEK)",
                          "icost_cvd_msek" = "Kostnader (MSEK)"
                        )
) # end switch                    

################################################################################
## Read text for How-to and About sections from file
################################################################################

how_to_text <- scan(switch(lang,
                           "eng" = "UI/how_to_eng.txt",
                           "swe" = "UI/how_to_swe.txt"
), what = "character"
)

about_text <- scan(switch(lang,
                          "eng" = "UI/about_eng.txt",
                          "swe" = "UI/about_swe.txt"
), what = "character"
)

################################################################################

## Function for running simulations in a separate R session. It's used both for
## baseline simulations (scen == 1) and scenarios
run_scen <- function(scen, input, sim_status, baseline_active, 
                     root=PROJECTROOT) {  
  
  ## Require that a baseline is loaded if scenario run
  if (scen!=1 & baseline_active()!="yes") {
    showNotification(switch(lang,
        "eng" = "You need an active baseline to run a scenario simulation",
        "swe" = "Det krävs en aktiv baslinje för att köra en scenariosimulering"
    ), duration=3
    )
    return(NULL)
  }
  
  ## Extract simulation and intervention years
  simyears <- input$slider_simyears
  intervention_years <- input[[paste0("slider_intervention_years_", scen)]]
  
  ## Require that start year is max 2022
  if (simyears[1] > 2022) {
    showNotification(switch(lang,
                  "eng" = "The simulation must start from year 2022 or earlier",
                  "swe" = "Simuleringen måste påbörjas år 2022 eller tidigare"
    ), duration=3
    )
    return(NULL)
  }
  
  ## Require at least three simulation years
  if (simyears[2] - simyears[1] < 2) {
    showNotification(switch(lang,
                    "eng" = "The simulation period must be minimum three years", 
                    "swe" = "Simuleringsperioden måste vara minst tre år lång"
    ), duration=3)
    return(NULL)
  }
  
  ## Require that intervention period is at least one year
  if (intervention_years[2] - intervention_years[1] < 1) {
    showNotification(switch(lang,
                        "eng" = "The phase-in period must be minimum one year", 
                        "swe" = "Infasningperioden måste vara minst ett år lång"
    ), duration=3)
    return(NULL)
  }
  
  ## Extract paths
  scen_name <- paste0("scen_", scen)
  scen_path <-  paste0(root, "/UI/runs/", scen_name, "/")
  baseline_path <- paste0(root, "/UI/scenarios/", input$scen_1_name, "/")
  
  ## Separate controls if baseline simulation, when scen == 1
  
  if (scen == 1) {
    
    ## Create a separate baseline directory if it doesn't exist already
    ## else give a confirmation prompt
    if (file.exists(baseline_path)) {
      showNotification(switch(lang, 
                      "eng" = "A baseline folder with that name already exists", 
                      "swe" = "En baslinjemapp med det namnet existerar redan"
      ), duration = 5
      )
      return(NULL)
    }
    
    dir.create(baseline_path)
    
    baseline_active("running")
    
    ## Write a json file with the selected parameters
    ## Note: make scenario names flexible later
    write_json(
      path = paste0(baseline_path, input$scen_1_name, ".json"),
      
      x = list(
        cpaf_cancer = input$cpaf_cancer,
        cpaf_cvd = input$cpaf_cvd,
        age_cutoff_cancer = input$age_cutoff_cancer,
        age_cutoff_cvd = input$age_cutoff_cvd,
        dcost_total_cancer = default_parameters$dcost_total_cancer,
        icost_total_cancer = default_parameters$icost_total_cancer,
        dcost_growth_cancer = default_parameters$dcost_growth_cancer,
        dcost_total_cvd = default_parameters$dcost_total_cvd,
        icost_total_cvd = default_parameters$icost_total_cvd,
        dcost_growth_cvd = default_parameters$dcost_growth_cvd,
        dcost_total_base_year = default_parameters$dcost_total_base_year,
        rr = default_parameters$rr,
        rr_cancer_diet = default_parameters$rr_cancer_diet,
        rr_cvd_diet = default_parameters$rr_cvd_diet,
        communalities = default_parameters$communalities
      ), auto_unbox = TRUE, pretty = TRUE
    ) # end write_json()
  } # if (scen == 1) {
  
  sim_status(replace(sim_status(), scen, "run"))
  showNotification(switch(lang,
                          "eng" = "Starting simulation ...",
                          "swe" = "Startar simulering ..."
  ), duration=3
  )
  cat("run", file = paste0(scen_path, "ui_status.txt"))
  cat(0, file=paste0(scen_path, "current_simyear.txt"))
  cat(0, file=paste0(scen_path, "last_simyear.txt"))
  futureCall(
    simulate_model, args = list(
      PROJECTROOT = root,
      is_baseline = (scen == 1),
      baseline_parameters = paste0(baseline_path, input$scen_1_name, ".json"),
      startyear = simyears[1],
      endyear = simyears[2],
      write_data_to_file = FALSE,
      ui = TRUE,
      ui_lang = lang,
      scen_name = scen_name,
      cfact = c(
        cfact_smoking = input[[paste0("cfact_smoking_", scen)]],
        cfact_alcohol = input[[paste0("cfact_alcohol_", scen)]],
        cfact_inactivity = input[[paste0("cfact_inactivity_", scen)]],
        cfact_obesity = input[[paste0("cfact_obesity_", scen)]]
      ),
      cfact_food = c(
        fruit = input[[paste0("cfact_fruit_", scen)]],
        wholegrains = input[[paste0("cfact_wholegrains_", scen)]],
        greens = input[[paste0("cfact_greens_", scen)]],
        meat = input[[paste0("cfact_meat_", scen)]],
        salt = input[[paste0("cfact_salt_", scen)]]
      ),
      cfact_startyear = intervention_years[1],
      cfact_endyear = intervention_years[2]
    )
  )
  NULL
} # end run_scen()

## Function for clearing a scenario
clear_scen <- function(scen, sim_status, dat_base, baseline_active, 
                       root=PROJECTROOT) {
  showNotification(paste0(switch(lang,
                                 "eng" = "Stopping/clearing scenario ", 
                                 "swe" = "Stoppar/rensar scenario "
  ),  
  scen, "..."), duration=3)
  scen_path <-  paste0(root, "/UI/runs/scen_", scen, "/")
  cat("stop", file = paste0(scen_path, "ui_status.txt"))
  sim_status(replace(sim_status(), scen, "stop"))
  ## Remove any .rda files from simulations, ignore warning if there was no file
  system(paste0("rm ", scen_path, "data/*.rda"), ignore.stderr = TRUE)
  dat_base[[paste0("scen_", scen)]] <- NULL
  
  ## Reset baseline controls and clear all scenarios if baseline is cleared
  if (scen == 1) {
    reset("scen_1_name")
    for (i in 1:length(baseline_slider_names)) {
      reset(baseline_slider_names[i])
    }
    sapply(2:n_scen_max, clear_scen, sim_status = sim_status, 
           dat_base = dat_base, baseline_active = baseline_active, root = root)
    baseline_active("no")
  } else if (scen %in% 2:n_scen_max) {
    reset(paste0("scen_", scen, "_name"))
    for (i in 1:length(scenario_slider_names)) {
      reset(paste0(scenario_slider_names[i], scen))
    }
  }
  NULL
} # end clear_scen()

## Function for loading a (non-baseline) scenario from a csv-file
load_scen <- function(scen, input, sim_status, dat_base) {
  
  ## Check that there is not an ongoing simulation 
  if (sim_status()[[scen]] %in% c("run", "active")) {
    showNotification(switch(lang,
                "eng" = "Stop/clear current scenario before loading from file",
                "swe" = "Stoppa/rensa nuvarande scenario innan du laddar från fil"
    )
    )
    return(NULL)
  }
  
  file <- input[[paste0("button_load_scen_", scen)]]
  
  dat <- tryCatch(
    fread(file$datapath, sep = ";", header=TRUE),
    error = function(e) {
      showNotification(paste0(switch(lang,
                       "eng" = "Problem reading data, error message: ",
                       "swe" = "Problem med att läsa in data, felmeddelande: "
      ),
      e$message)) }
  ) # end tryCatch
  
  if (is.null(dat)) return(NULL)
  
  if (all(colnames(dat) == all_var_names) == FALSE) {
    showNotification(switch(lang,
              "eng" = "Variable names of data not matching NCD-SIM output data",
              "swe" = "Variabelnamn i data matchar inte utdata från NCD-SIM"
    )
    )
    return(NULL)
  }
  
  ## Check that scenario comes from the loaded baseline
  if (!isTRUE(all.equal(dat_base[["scen_1"]][, ..baseline_vars_to_check],
                        dat[, ..baseline_vars_to_check]))) {
    showNotification(switch(lang,
                            "eng" = "Scenario does not match baseline",
                            "swe" = "Scenario matchar inte baslinje"
    )
    )
    return(NULL)
  }
  
  ## Update scenario name
  scen_name <- strsplit(file$name, ".", fixed = TRUE)[[1]][1]
  updateTextInput(inputId = paste0("scen_", scen, "_name"), value = scen_name)
  
  ## Create scenario id variable
  dat[, scen_id := paste0("scen_", scen)]
  
  ## Update data and sim_status()
  dat_base[[paste0("scen_", scen)]] <- dat
  sim_status(replace(sim_status(), scen, "active"))
  
  ## Disable run button. Note that this perhaps should be done at the top of
  ## this function, but then it should be reverted if check vz baseline fails
  disable(paste0("button_run_scen_", scen))
} # end load_scen

## Function for scanning a simulation log file, returning its contents
scan_log <- function(scen, root = PROJECTROOT) {
  paste(scan(file =
               paste0(root, "/UI/runs/scen_", scen, "/log.txt"),
             what="character", sep="\n", quiet = TRUE), collapse="\n")
}

## Function for creating plots by aggregating and modifying the raw
## simulation output collected from reactiveVals dat_base(). Some other objects 
## for plotting are also created

create_plot <- function(dat_list, plot_type, var_to_plot, group, measure, years, 
                        var_list, var_labels, input, data_only = FALSE) {
  
  
  ## Combine data from all scenarios
  dat_wide <- rbindlist(dat_list, fill = TRUE)
  
  ## Extract scenario ids and number of scenarios
  scen_ids <- dat_wide[, unique(scen_id)]
  n_scen <- length(scen_ids)
  
  if (measure %in% c("abs_diff", "rel_diff", "rel_diff_pop_share") & n_scen<2) {
    showNotification(switch(lang,
          "eng" = "Please load a comparison scenario", 
          "swe" = "Du behöver ladda ett jämförelsescenario för denna diagramtyp"
    ), duration = 3
    )
    measure <- "abs"
  }
  
  ## Derived variables
  dat_wide[, s_ncd := s_cvd + s_cancer]
  dat_wide[, f_cancer := f_pop_cancer + f_immig_cancer]
  dat_wide[, f_cvd := f_pop_cvd + f_immig_cvd]
  dat_wide[, f_ncd := f_cancer + f_cvd]
  dat_wide[, f_ncd_dead := f_cancer_dead + f_cvd_dead]
  dat_wide[, f_total_dead := f_dead + f_ncd_dead]
  dat_wide[, total_pop := s_pop + s_cancer + s_cvd]
  dat_wide[, s_alcohol := prev_alcohol * total_pop]
  dat_wide[, s_obesity := prev_obesity * total_pop]
  dat_wide[, s_smoking := prev_smoking * total_pop]
  dat_wide[, s_inactivity := prev_inactivity * total_pop]
  
  dat_wide[, s_fruit := prev_fruit * total_pop]
  dat_wide[, s_wholegrains := prev_wholegrains * total_pop]
  dat_wide[, s_greens := prev_greens * total_pop]
  dat_wide[, s_meat := prev_meat * total_pop]
  dat_wide[, s_salt := prev_salt * total_pop]
  
  ## Compute costs in million SEK
  dat_wide[, dcost_ncd_msek := 
             (dcost_cancer + dcost_cvd) / 10^6]
  dat_wide[, dcost_cancer_msek := dcost_cancer / 10^6]
  dat_wide[, dcost_cvd_msek := dcost_cvd / 10^6]
  dat_wide[, icost_ncd_msek := 
             (icost_cancer + icost_cvd) / 10^6]
  dat_wide[, icost_cancer_msek := icost_cancer / 10^6]
  dat_wide[, icost_cvd_msek := icost_cvd / 10^6]
  
  ## Transform into long format. Warning due to different data types of the 
  ## value variables is suppressed here
  dat_long <- suppressWarnings(
    melt(dat_wide, id.vars = c("scen_id", "year", "sex", "age"))
  )
  
  
  ## Extract scenario names
  scen_names <- character()
  
  ## Set default formatting function for axes and hover label
  if (measure %in% c("pop_share", "rel_diff", "rel_diff_pop_share")) {
    format_fun <- pretty_print
  } else {
    format_fun <- int_print
  }
  
  for (scen in 1:n_scen) {
    scen_names[scen] <- input[[paste0(scen_ids[scen], "_name")]]
  }
  
  dat_long <- merge(dat_long,
                    data.table(scen_id = scen_ids, scen_name = scen_names),
                    by = "scen_id")
  
  dat_long[, group := ""]
  
  ## Customize group variable
  if (group == "total") {
    dat_long[, group := 
               fifelse(lang == "eng", "Total pop.", "Hela befolkningen")]
  } else if (group == "sex") {
    if (lang == "eng") {
      dat_long[, group := fifelse(sex == 1, "Men", "Women")]
    } else if (lang == "swe") {
      dat_long[, group := fifelse(sex == 1, "Män", "Kvinnor")]
    }
  } else if (group == "age") {
    dat_long[age %in% 0:9, group := "0-9 years"]
    dat_long[age %in% 10:19, group := "10-19 years"]
    dat_long[age %in% 20:29, group := "20-29 years"]
    dat_long[age %in% 30:39, group := "30-39 years"]
    dat_long[age %in% 40:49, group := "40-49 years"]
    dat_long[age %in% 50:59, group := "50-59 years"]
    dat_long[age %in% 60:69, group := "60-69 years"]
    dat_long[age %in% 70:79, group := "70-79 years"]
    dat_long[age %in% 80:89, group := "80-89 years"]
    dat_long[age >= 90, group := "90+ years"]
    if (lang == "swe") {
      dat_long[, group := sub("years", "år", group)]
    }
  } else if (group == "age_70p") {
    if (lang == "eng") {
      dat_long[, group := fifelse(age %in% 0:69, "0-69 years", "70+ years")]
    } else if (lang == "swe") {
      dat_long[, group := fifelse(age %in% 0:69, "0-69 år", "70+ år")]
    }
  } else if (group == "age_16_84") {
    dat_long <- dat_long[age >= 16, ]
    if (lang == "eng") {
      #dat_long[age <= 15, group := "0-15 years"]
      dat_long[age %in% 16:84, group := "16-84 years"]
      dat_long[age >= 85, group := "85+ years"]
    } else if (lang == "swe") {
      #dat_long[age <= 15, group := "0-15 år"]
      dat_long[age %in% 16:84, group := "16-84 år"]
      dat_long[age >= 85, group := "85+ år"]
    }
  }
  
  ## Aggregate by the chosen group, and on year or age, depending on plot type
  if (plot_type == "ts") {
    dat_by_group <- dat_long[, .(value = sum(value)), 
                             by = c("scen_id", "scen_name", "year",
                                    "group", "variable")]
    aes_x <- aes(x = year)
    
  } else if (plot_type == "cs") {
    
    dat_by_group <- dat_long[year == years, .(value = sum(value)), 
                             by = c("scen_id", "scen_name", "age",
                                    "group", "variable")]
    aes_x <- aes(x = age)
    
  } else if (plot_type == "acc") {
    
    dat_by_group <- dat_long[year %in% years[1]:years[2], 
                            .(value = sum(value, na.rm = TRUE)), 
                            by = c("scen_id", "scen_name", "group", "variable")]
  }
  
  ## Customize plot objects to group
  aes_group <- aes(group = interaction(scen_id, group), lty = group)
  group_labels <- dat_by_group[, unique(group)]
  if (group == "total") {
    legend_labels <- scen_names
  } else {
    legend_labels <- c(sapply(scen_names,
                              function(x) paste(x, group_labels, sep = ": ")))
  }
  
  ## Plot formatting. Some of these will be overridden below.
  group_names <- switch(lang,
                        "eng" = 
                          list("total" = "Total population.",
                               "sex" = "By sex.",
                               "age" = "By age group.", 
                               "age_70p" = "By age group.",
                               "age_16_84" = "By age group."),
                        "swe" = 
                          list("total" = "Hela befolkningen.",
                               "sex" = "Per kön.",
                               "age" = "Per åldersgrupp.", 
                               "age_70p" = "Per åldersgrupp.",
                               "age_16_84" = "Per åldersgrupp.")
  ) # end switch                     
  
  ylim_custom <- c(0, NA)
  line_colours <- fohm_colours()[!sapply(dat_list, is.null)]
  scale_y_custom <- scale_y_continuous(labels = int_print)
  ylab_custom <- var_labels[[var_to_plot]]
  title_custom <- paste0(gsub(".", ": ", names(unlist(var_list)), 
                              fixed = TRUE)[unlist(var_list) == var_to_plot],
                         ". ", 
                         group_names[group]
  )
  
  ## Create graph with population shares (per capita)
  if (measure == "pop_share") {
    
    dat_by_group[variable == var_to_plot, value := value / 
                   dat_by_group[variable == "total_pop", value]]
    if (var_to_plot %in% c("dcost_ncd_msek", 
                           "dcost_cvd_msek", 
                           "dcost_cancer_msek",
                           "icost_ncd_msek", 
                           "icost_cvd_msek", 
                           "icost_cancer_msek")) {
      dat_by_group[variable == var_to_plot, value := value * 10^6]
      ylab_custom <- switch(lang,
                            "eng" = "Cost (SEK), per capita",
                            "swe" = "Kostnad (SEK), per capita")
      format_fun <- int_print
    } else if (var_to_plot %in% c("s_ncd",
                                  "f_ncd",
                                  "s_cancer",
                                  "f_cancer",
                                  "s_cvd",
                                  "f_cvd")
    ) {
      dat_by_group[variable == var_to_plot, value := value * 10^3]
      ylab_custom <- switch(lang,
                            "eng" = "Cases per 1 000 population",
                            "swe" = "Fall per tusen invånare")
      format_fun <- int_print
    } else if (var_to_plot %in% c("f_ncd_dead",
                                  "f_cancer_dead",
                                  "f_cvd_dead",
                                  "f_total_dead")
    ) {
      dat_by_group[variable == var_to_plot, value := value * 10^3]
      ylab_custom <- switch(lang,
                            "eng" = "Deaths per 1 000 population",
                            "swe" = "Dödsfall per tusen invånare")
      format_fun <- pretty_print_1
    } else if (var_to_plot == "f_born") {
      dat_by_group[variable == var_to_plot, value := value * 10^3]
      ylab_custom <- switch(lang,
                            "eng" = "Births per 1 000 population",
                            "swe" = "Födslar per tusen invånare")
      format_fun <- int_print
    }
    else {
      ylab_custom <- switch(lang,
                            "eng" = paste0(ylab_custom, ", population share"),
                            "swe" = paste0(ylab_custom, ", befolkningsandel")
      )
    }
  }
  
  ## Create graph diffed in absolutes vz baseline
  if (measure == "abs_diff" & n_scen > 1) {
    
    for (scen in n_scen:2) {
      dat_by_group[scen_id == scen_ids[scen], value := value -
                     dat_by_group[scen_id == scen_ids[1], value]]
    }
    
    dat_by_group <- dat_by_group[scen_id!="scen_1", ]
    legend_labels <- setdiff(legend_labels, 
                             grep(scen_names[1], legend_labels, value = TRUE))
    ylim_custom <- c(NA, NA)
    line_colours <- line_colours[-1]
    ylab_custom <- paste0(ylab_custom, switch(lang,
                                          "eng" = ", scenario minus baseline",
                                          "swe" = ", scenario minus baslinje"
    )
    )
  }
  
  ## Create graph diffed as a proportion vz baseline
  if (measure == "rel_diff" & n_scen > 1) {
    
    for (scen in n_scen:2) {
      dat_by_group[scen_id == scen_ids[scen], 
                   value := value / dat_by_group[scen_id == scen_ids[1], value]]
    }
    dat_by_group <- dat_by_group[scen_id!="scen_1", ]
    
    legend_labels <- setdiff(legend_labels, 
                             grep(scen_names[1], legend_labels, value=TRUE))
    ylim_custom <- c(NA, NA)
    line_colours <- line_colours[-1]
    ylab_custom <- paste0(ylab_custom, switch(lang,
                                        "eng" = ", scenario to baseline ratio",
                                        "swe" = ", kvot scenario/baslinje"
    )
    )
  }
  
  ## Create graph as a proportion of baseline, of population shares
  if (measure == "rel_diff_pop_share" & n_scen > 1) {
    
    dat_by_group[variable == var_to_plot, value := value / 
                   dat_by_group[variable == "total_pop", value]]
    
    for (scen in n_scen:2) {
      dat_by_group[scen_id == scen_ids[scen], 
                   value := value / dat_by_group[scen_id == scen_ids[1], value]]
    }
    dat_by_group <- dat_by_group[scen_id!="scen_1", ]
    
    legend_labels <- setdiff(legend_labels, 
                             grep(scen_names[1], legend_labels, value=TRUE))
    ylim_custom <- c(NA, NA)
    line_colours <- line_colours[-1]
    ylab_custom <- paste0(ylab_custom, switch(lang,
                    "eng" = ", scenario to baseline ratio of population shares",
                    "swe" = ", kvot scenario/baslinje, befolkningsandelar"
    )
    )
  }
  
  dat_to_plot <- dat_by_group[variable == var_to_plot, ]
  
  if (plot_type == "cs") {
    dat_to_plot[, year := years]
  }
  
  if (lang == "eng") {
    hover_pre_group <- "<br>Group: "
    hover_pre_year <- "<br>Year: "
    hover_pre_value <- "<br>Value: "
  } else if (lang == "swe") {
    hover_pre_group <- "<br>Grupp: "
    hover_pre_year <- "<br>År: "
    hover_pre_value <- "<br>Värde: "
  }
  
  if (plot_type %in% c("ts", "cs")) {
    dat_to_plot[, hover_label := paste0(scen_name,
                                      hover_pre_group, group,
                                      hover_pre_year, year,
                                      hover_pre_value, sapply(value, format_fun)
    )]
  } else if (plot_type == c("acc")) {
    dat_to_plot[, hover_label := paste0(scen_name,
                                      hover_pre_group, group,
                                      hover_pre_value, sapply(value, format_fun)
    )]
  }
  
  if (data_only == TRUE) {
    return(dat_to_plot)
  }
  
  ## Create custom x-axis breaks and labels
  
  if (plot_type == "ts") {
    xlab_custom <- switch(lang, "eng" = "Year", "swe" = "År")
    year_min <- years[1]
    year_max <- years[2]
    if (year_max - year_min < 35) {
      x_breaks <- year_min:year_max
    } else {
      x_breaks <- seq(from = year_min, to = year_max, by = 5)
    }
  } else if (plot_type == "cs") {
    xlab_custom <- switch(lang, "eng" = "Age", "swe" = "Ålder")
    x_breaks <- seq(0, 100, 5)
  }
  
  if (plot_type %in% c("ts", "cs")) {
    plot_tmp <- ggplot(dat_to_plot, 
                       aes(y = value, 
                           col = scen_id,
                           text = hover_label)) +
      geom_line() + 
      scale_colour_manual(values = line_colours) + 
      theme(legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      scale_x_continuous(breaks = x_breaks) +
      scale_y_continuous(labels = format_fun) +
      coord_cartesian(ylim = ylim_custom) +
      aes_x + 
      aes_group + 
      labs(title = title_custom, y = ylab_custom, x = xlab_custom)
    
  } else if (plot_type == "acc") {
    
    bar_labels <- legend_labels
    names(bar_labels) <- dat_to_plot[,
                                     as.character(interaction(scen_id, group))]
    plot_tmp <- ggplot(dat_to_plot, 
                       aes(x = group,
                           y = value, 
                           group = scen_id,
                           alpha = group,
                           fill = scen_id,
                           col = scen_id,
                           text = hover_label)) +
      geom_col(width=0.5, position = position_dodge()) +
      scale_colour_manual(values = line_colours) + 
      scale_fill_manual(values = line_colours) + 
      scale_alpha_manual(values = c(1, 0.6)) + 
      theme(legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      scale_x_discrete(labels = bar_labels) +
      scale_y_continuous(labels = format_fun) +
      labs(title = title_custom, y = ylab_custom, x = NULL)
    
    
  }
  ## Plotly formatting
  plot_tmp <- layout(ggplotly(plot_tmp, tooltip = "text"),
                     legend = list(font = list(size = 12), orientation = "v"))
  ## Fix Swedish labels for sex here so that they appear in the correct order
  if (lang == "swe" & group == "sex") {
    plot_tmp <- change_legend(plot_tmp, title = "", labels = rev(legend_labels))
  } else {
    plot_tmp <- change_legend(plot_tmp, title = "", labels = legend_labels)
  }
  return(plot_tmp)
} # end function create_plot()

################################################################################
## UI
################################################################################
ui <- fluidPage(
  ## Customize colour of sliders and headers
  useShinyjs(),
  chooseSliderSkin("Shiny", color = fohm_colours()[1]),
  tags$head(tags$style("h2 {color:#009F80;}")),
  tags$head(tags$style("h4 {color:#009F80;}")),
  tags$head(tags$style("h5 {color:#009F80;}")),
  tags$head(tags$style("a {color:#009F80}")),
  tags$style(HTML(".panel-title {font-size: 13px;}")),
  
  # Application title
  titlePanel("NCD-SIM"),
  
  # Sidebar with controls of simulation and for loading data
  sidebarLayout(
    sidebarPanel(width = 4,
                 tabsetPanel(
                   
                   ## Text input for custom scenario name 
                   tabPanel(uiOutput("baseline_header"),
                            
                        tags$style(HTML("
        .label-left .form-group {
          display: flex;              /* Use flexbox for positioning children */
          flex-direction: row;        /* Place children on a row (default) */
          width: 100%;                /* Set width for container */
          max-width: 350px;
        }
      
        .label-left label {
          margin-right: 1rem;         /* Add spacing between label and slider */
          align-self: center;         /* Vertical align in center of row */
          text-align: right;
          width: 75px;          /* Target width for label */
        }
      
        .label-left .irs {
          width: 275px;          /* Target width for slider */
        }
        "
                        )),
                            
  ## File input for loading pre-simulated baseline data
  h5(switch(lang, 
            "eng" = "Load baseline from folder", 
            "swe"= "Ladda baslinje från mapp"),
     info_icon("Load pre-simulated baseline data and parameters from folder")),
  
  shinyDirButton("button_load_scen_1", 
                 switch(lang, 
                    "eng" = "Select a folder", 
                    "swe"= "Välj mapp"), 
                 switch(lang, 
                    "eng" = "Please select a folder", 
                    "swe"= "Välj mapp"), 
                 FALSE),
  
  h5(switch(lang, 
            "eng" = "Baseline name", 
            "swe"= "Baslinjenamn"), 
     info_icon(switch(lang, 
      "eng" = "Choose a custom name for the baseline", 
      "swe"= "Välj namn för baslinjen")
     )),
  
                    textInput("scen_1_name", label = NULL, value = "baseline"),
                            
                            ## Slider input for simulation period
                        h5(switch(lang, 
                                  "eng" = "Simulation period", 
                                  "swe"= "Simuleringsperiod"),
                               info_icon(
                                 switch(lang, 
          "eng" = "Choose start year (<= 2022) and end year for new simulation", 
          "swe"= "Välj startår (<= 2022) och slutår för ny simulering")
                               )),
                            
                            sliderInput(
                              "slider_simyears",
                              label = NULL,
                              min = start_year_min,
                              max = end_year_max,
                              value = c(2020, 2040),
                              step = 1,
                              sep = "",
                              width = slider_width),
                            
                            ## Sets of slider inputs for adjusting risk factors
                            h5(switch(lang, 
                                      "eng" = "Risk group size", 
                                      "swe"= "Storlek på riskgrupp"), 
                               info_icon(risk_factor_info_text)
                            ),
                            
bsCollapse(id = "risk-factors-baseline",
           bsCollapsePanel(switch(lang, 
                                  "eng" = "(+) Non-dietary", 
                                  "swe"= "(+) Icke-kost"),
                           
                           risk_factor_slider("cfact_smoking_1",
                                              label = label_slider_smoking,
                                              info_text = info_text_smoking
                           ),
                           risk_factor_slider("cfact_alcohol_1",
                                              label = label_slider_alcohol,
                                              info_text = info_text_alcohol
                           ),
                           risk_factor_slider("cfact_inactivity_1",
                                              label = label_slider_inactivity, 
                                              info_text = info_text_inactivity 
                           ),
                           risk_factor_slider("cfact_obesity_1",
                                              label = label_slider_obesity,
                                              info_text = info_text_obesity
                           )
           ), # end bsCollapsePanel("Non-dietary"
           bsCollapsePanel(switch(lang, 
                                  "eng" = "(+) Dietary", 
                                  "swe"= "(+) Kost"),
                           
                           risk_factor_slider("cfact_fruit_1",
                                              label = label_slider_fruit,
                                              info_text = info_text_fruit
                           ),
                           risk_factor_slider("cfact_wholegrains_1",
                                              label = label_slider_wholegrains,
                                              info_text = info_text_wholegrains
                           ),
                           risk_factor_slider("cfact_greens_1",
                                              label = label_slider_greens,
                                              info_text = info_text_greens
                           ),
                           risk_factor_slider("cfact_meat_1",
                                              label = label_slider_meat,
                                              info_text = info_text_meat
                           ),
                           risk_factor_slider("cfact_salt_1",
                                              label = label_slider_salt,
                                              info_text = info_text_salt
                           )
                                       ) # end bsCollapsePanel("Dietary"
                            ), # end bsCollapse(id = "risk-factors-baseline"
                            
                            ## Slider input for linear phase-in of risk group adjustment
                            h5(switch(lang, 
                                      "eng" = "Phase-in period", 
                                      "swe"= "Infasningsperiod"),
                               info_icon(switch(lang, 
                                              "eng" = "Changes in the size of the risk groups are phased in linearly during this period and remain thereafter", 
                                              "swe"= "Förändringar i storleken på riskgrupperna fasas in linjärt under denna period och kvarstår därefter"))
                            ),
                            
                            sliderInput(
                              "slider_intervention_years_1",
                              label = NULL,
                              min = start_year_min,
                              max = end_year_max,
                              value = c(2025, 2026),
                              step = 1,
                              sep = "",
                              width = slider_width),
                            
                            ## Sets of sliders for advanced options
                            h5(switch(lang, 
                                      "eng" = "Advanced options", 
                                      "swe"= "Avancerade inställningar"),
         info_icon(switch(lang, 
          "eng" = "Click to expand and use sliders to adjust", 
          "swe"= "Klicka för att expandera och använd reglagen för att justera")
                               )),
                            
bsCollapse(id = "advanced-options",
   bsCollapsePanel(switch(lang, 
                          "eng" = "(+) Advanced options", 
                          "swe"= "(+) Avancerade inställningar"),
                   
                   sliderInput(
                     inputId = "cpaf_cancer",
                     label = switch(lang,
                                    "eng" = span("Cancer explained percentage",
                                                 info_icon("Percentage share of cancer accounted for by the risk factors in the model")),
                                    "swe" = span("Cancer förklaringsandel",
                                                 info_icon("Procentandelen av cancer som tillskrivs modellens riskfaktorer")),
                     ), # end switch
                     min = 0,
                     max = 1.0,
                     value = default_parameters$cpaf_cancer,
                     step = 0.01),
                   
                   sliderInput(
                     inputId = "cpaf_cvd",
                     label = switch(lang,
                                    "eng" = span("CVD explained percentage",
                                                 info_icon("Percentage share of CVD accounted for by the risk factors in the model")),
                                    "swe" = span("Hjärt-kärl förklaringsandel",
                                                 info_icon("Procentandelen av hjärt-kärlsjukdom som tillskrivs modellens riskfaktorer")),
                     ), # end switch
                     min = 0,
                     max = 1.0,
                     value = default_parameters$cpaf_cvd,
                     step = 0.01),
                   
                   sliderInput(
                     inputId = "age_cutoff_cancer",
                     label = switch(lang,
                                    "eng" = span("Age threshold cancer risk",
                                                 info_icon("The minimum age from which cancer is affected by risk factors")),
                                    "swe" = span("Ålderströskel för cancerrisk",
                                                 info_icon("Lägsta ålder från vilken cancer påverkas av riskfaktorer")),
                     ), # end switch
                     min = 0,
                     max = 100,
                     value = default_parameters$age_cutoff_cancer,
                     step = 1),
                   
                   sliderInput(
                     inputId = "age_cutoff_cvd",
                     label = switch(lang,
                                    "eng" = span("Age threshold CVD risk",
                                                 info_icon("The minimum age from which CVD is affected by risk factors")),
                                    "swe" = span("Ålderströskel för hjärt-kärlrisk",
                                                 info_icon("Lägsta ålder från vilken hjärt-kärlsjukdomar påverkas av riskfaktorer")),
                     ), # end switch
                     min = 0,
                     max = 100,
                     value = default_parameters$age_cutoff_cvd,
                     step = 1)
   ) # end bsCollapsePanel("advanced-options"
),
                            
                            ## Buttons for start/stop and pausing simulation
                            h5(switch(lang,
                                      "eng" = "Simulation controls",
                                      "swe" = "Simuleringskontroller")
                            ),
                            actionButton("button_run_scen_1", 
                                         switch(lang, 
                                                "eng" = "Run", 
                                                "swe"= "Starta")),
                            
                            actionButton("button_clear_scen_1", 
                                         switch(lang, 
                                                "eng" = "Stop/clear", 
                                                "swe"= "Stoppa/rensa")),
                            
                            ## Text output showing simulation progress
                          h5(switch(lang, 
                                    "eng" = "Log", 
                                    "swe"= "Logg")),
                          
                          tags$small(verbatimTextOutput("log_text_1")),
                          
                          disabled(downloadButton("button_save_baseline_data", 
                                                  switch(lang, 
                                                         "eng" = "Save data", 
                                                         "swe"= "Spara data"))),
                          
                          downloadButton("button_save_baseline_parameters", 
                                         switch(lang, 
                                                "eng" = "Save parameters", 
                                                "swe"= "Spara parametrar"))
                            
                   ),
                   
                   tabPanel(h5(switch(lang, 
                                      "eng" = "Scenarios", 
                                      "swe"= "Scenarier")),
                            uiOutput("scenarioPanels")
                   ) # end tabPanel
                 ) # end tabsetPanel
    ), # end sidebarPanel
    
    # Main panel
    mainPanel(width = 8,
              tabsetPanel(
                tabPanel(tags$h5(switch(lang, 
                                        "eng" = "Time series", 
                                        "swe"= "Tidsserier")),
                         fluidRow(
                           column(4, selectInput("ts_var_to_plot", 
                                                 switch(lang, 
                                                        "eng" = "Variable", 
                                                        "swe"= "Variabel"), 
                                                 choices = var_list)
                           ),
                           column(4, selectInput("ts_group",
                                                 switch(lang, 
                                                        "eng" = "Group", 
                                                        "swe"= "Grupp"), 
                                                 choices = ts_group_list)
                           ),
                           column(4, selectInput("ts_measure", 
                                                 switch(lang, 
                                                        "eng" = "Measure", 
                                                        "swe"= "Utfallsmått"), 
                                                 choices = ts_measure_list)
                           )
                         ), # end fluidRow
                         plotlyOutput("plot_time_series"),
                         uiOutput("render_button_save_ts_data")
                ), # tabPanel(tags$h5("Time series"),
                tabPanel(tags$h5(switch(lang, 
                                        "eng" = "Cross section", 
                                        "swe"= "Tvärsnitt")),
                         fluidRow(
                           column(4, selectInput("cs_var_to_plot", 
                                                 switch(lang, 
                                                        "eng" = "Variable", 
                                                        "swe"= "Variabel"), 
                                                 choices = var_list)
                           ),
                           column(4, selectInput("cs_group", 
                                                 switch(lang, 
                                                        "eng" = "Group", 
                                                        "swe"= "Grupp"), 
                                                 choices = cs_group_list)
                           ),
                           column(4, selectInput("cs_measure", 
                                                 switch(lang, 
                                                        "eng" = "Measure", 
                                                        "swe"= "Utfallsmått"), 
                                                 choices = cs_measure_list)
                           )
                         ), # end fluidRow
                         plotlyOutput("plot_cross_section"),
                         uiOutput("slider_year_cs_plot"),
                         uiOutput("render_button_save_cs_data")
                ), # tabPanel(tags$h5("Cross section"),
                tabPanel(tags$h5(switch(lang, 
                                        "eng" = "Accumulated", 
                                        "swe"= "Ackumulerat")),
                         fluidRow(
                           column(4, selectInput("acc_var_to_plot", 
                                                 switch(lang, 
                                                        "eng" = "Variable", 
                                                        "swe"= "Variabel"), 
                                                 choices = var_list)
                           ),
                           column(4, selectInput("acc_group", 
                                                 switch(lang, 
                                                        "eng" = "Group", 
                                                        "swe"= "Grupp"), 
                                                 choices = acc_group_list)
                           ),
                           column(4, selectInput("acc_measure", 
                                                 switch(lang, 
                                                        "eng" = "Measure", 
                                                        "swe"= "Utfallsmått"), 
                                                 choices = acc_measure_list)
                           )
                         ), # end fluidRow
                         plotlyOutput("plot_accumulated"),
                         uiOutput("slider_years_acc_plot"),
                         uiOutput("render_button_save_acc_data")
                ), # end tabPanel(tags$h5("Accumulated"),
                tabPanel(tags$h5(switch(lang, 
                                        "eng" = "How-to", 
                                        "swe"= "Instruktioner")
                ),
                lapply(how_to_text, p)
                ), # tabPanel How-to
                tabPanel(tags$h5(switch(lang, 
                                        "eng" = "About", 
                                        "swe"= "Om modellen")),
                         lapply(about_text, p),
                         
                         img(src = switch(lang,
                                          "eng" = "model_flowchart_eng.png", 
                                          "swe" = "model_flowchart_swe.png"),
                             align = "left", width="75%")
                         
                ) # tabPanel(tags$h5("About"),
              ) # tabsetPanel(
    ) # end mainPanel
  ) # end sidebarLayout
) # end UI

################################################################################
## Server
################################################################################
server <- function(input, output, session) {
  
  ## Increase max limit for loading files
  options(shiny.maxRequestSize=32*1024^2)
  
  ## Create reactive values
  baseline_active <- reactiveVal("no")
  
  sim_status <- reactiveVal(rep("stop", n_scen_max))
  
  ## Reactive data object for all scenario data
  dat_base <- reactiveValues() 
  
  for (scen in 1:n_scen_max) {
    dat_base[[paste0("scen_", scen)]] <- NULL
  }
  
  baseline_parameters <- reactiveVal(read_json(default_parameters_path))
  
  baseline_controls_active <- reactiveVal(TRUE)
  
  timer_log <- reactiveTimer(500)
  timer_data <- reactiveTimer(5000)
  
  ## Observer for reading logs for active simulations
  observeEvent(timer_log(), {
    lapply(1:n_scen_max, function(scen) {
      outputId <- paste0("log_text_", scen)
      
      if (sim_status()[[scen]] == "run") {
        log_text <- scan_log(scen)
        output[[outputId]] <- renderText(log_text)
      } else if (sim_status()[[scen]] == "stop") {
        log_text <- switch(lang,
                           "eng" = "Waiting to start or load", 
                           "swe" = "Väntar på att starta eller ladda")
        output[[outputId]] <- renderText(log_text)
      } else if (sim_status()[[scen]] == "active") {
        log_text <- switch(lang,
                           "eng" = "Scenario active",
                           "swe" = "Scenario aktivt"
        )
        output[[outputId]] <- renderText(log_text)
      }
    }) # end lapply
    NULL
  })
  
  ## Observer for reading data for active simulations
  ## sim_status() is also set to "stop" when run completed
  observeEvent(timer_data(), {
    lapply(1:n_scen_max, function(scen) {
      if (sim_status()[[scen]] == "run") {
        
      data_path <- paste0("UI/runs/scen_", scen, "/data/")
      ## Extract path to the .rda file. Normally there is only one file, but if
      ## the simulation script fails to remove the previous one, this will be a
      ## vector and so we should account for that
      data_file <- dir(data_path)
        if (length(data_file) == 1) {
          ## Load data and assign to reactive data list
          
          tmp <- tryCatch(load(file = paste0(data_path, data_file)), 
                          error = function(e) {
                            cat("Problem reading simulation data. If the simulations stops entirely, stop/clear the simulation and restart. Error message: ", 
                                e$message) }
          )
          
          ## Proceed with this scenario only if data loading worked
          if (is.null(tmp) == FALSE) {
            ## Data check
            if (all(colnames(dat) == all_var_names) == FALSE) {
              showNotification(switch(lang,
                                      "eng" = "Variable names of data not matching NCD-SIM output data, perhaps try clearing scenario and restart.",
                                      "swe" = "Variabelnamn i data matchar inte utdata från NCD-SIM, prova kanske att rensa scenario och starta om."
              )
              )
              return(NULL)
            }
            
            ## Create scenario id variable
        dat[, scen_id := paste0("scen_", scen)]
        
        dat_base[[paste0("scen_", scen)]] <- dat
        ## Check if loaded data was the last year and if so update sim_status()
        if(dat[.N, s_pop!=0]) {
          sim_status(replace(sim_status(), scen, "active"))
          ## Activate baseline if a baseline run was completed
          if (scen == 1) baseline_active("yes")
        } # if(dat[.N, s_pop!=0])
      } # if (!is.null(tmp))
    } else if (length(data_file) > 1) {
      file.remove(paste0(data_path, dir(data_path)))
    }
      } # if (sim_status()[[scen]] == "run")
    }) # end lapply()
    NULL
  })
  
 ## Observer for when baseline becomes active/inactive/starts running, to update
 ## the controls accordingly
  observeEvent(baseline_active(), {
    
    if (baseline_active() == "yes") {
      baseline_controls_active(FALSE)
      disable("button_run_scen_1")
      disable("button_load_scen_1")
      disable("scen_1_name")
      sim_status(replace(sim_status(), 1, "active"))
      enable("button_save_baseline_data")
    }
    
    if (baseline_active() == "running") {
      baseline_controls_active(FALSE)
      disable("button_load_scen_1")
      disable("scen_1_name")
    }
    
    if (baseline_active() == "no") {
      baseline_controls_active(TRUE)
      disable("button_save_baseline_data")
      enable("button_load_scen_1")
      enable("scen_1_name")
      enable("button_run_scen_1")
    }
  })
  
  ## Observer for when status of scenarios change, updating controls accordingly
  observeEvent(sim_status(), {
    for (scen in 2:n_scen_max) {
      if (sim_status()[scen] == "active") {
        enable(paste0("button_save_scen_", scen))
        disable(paste0("button_run_scen_", scen))
      } else if (sim_status()[scen] == "stop") {
        disable(paste0("button_save_scen_", scen))
        enable(paste0("button_run_scen_", scen))
      } else if (sim_status()[scen] == "run") {
        disable(paste0("button_save_scen_", scen))
        disable(paste0("button_run_scen_", scen))
      }
    } # end for
  })
  
  ## Setup for baseline folder selection
  volumes <- getVolumes()
  roots <- c(wd="UI/scenarios")
  shinyDirChoose(input, "button_load_scen_1", roots=roots, session=session,
                 allowDirCreate = FALSE)
  
  observeEvent(input$button_load_scen_1, {
    
    baseline_path <- parseDirPath(roots, input$button_load_scen_1)
    
    if (length(baseline_path) > 0) {
      
      tmp <- strsplit(baseline_path, "/")[[1]]
      
      baseline_name <- tmp[length(tmp)]
      
      ## Sanity checks
      baseline_data_path <- paste0(baseline_path, "/", baseline_name, ".csv")
      baseline_parameters_path <- paste0(baseline_path, "/", baseline_name, 
                                         ".json")
      
      if ((file.exists(baseline_data_path) & 
           file.exists(baseline_parameters_path)) == FALSE) {
        showNotification(switch(lang,
            "eng" = "Data or parameter file with name matching folder missing",
            "swe" = "Data eller parameterfil med namn som matchar mapp saknas"
)
        )
        return(NULL)
      }
      
      ## Read baseline data file
      dat <- tryCatch(fread(file = baseline_data_path, sep=";", header=TRUE),
              error = function(e) {
                showNotification(paste0(switch(lang,
                       "eng" = "Problem reading data, error message: ",
                       "swe" = "Problem med att läsa in data, felmeddelande: "
                ), e$message)) }
      ) # end tryCatch
      
      if (is.null(dat)) return(NULL)
      
      if (all(colnames(dat) == all_var_names) == FALSE) {
        showNotification(switch(lang,
              "eng" = "Variable names of data not matching NCD-SIM output data",
              "swe" = "Variabelnamn i data matchar inte utdata fråb NCD-SIM"
        )
        )
        return(NULL)
      }
      
      ## Read parameter file from same directory
      baseline_par_tmp <- 
        tryCatch(read_json(baseline_parameters_path),
                 error = function(e) {
                   showNotification(paste0(switch(lang,
                    "eng" = "Problem reading parameters, error message: ",
                    "swe" = "Fel vid inläsning av parametrar, felmeddelande: "
                   ), e$message)) }
        ) # end tryCatch
      
      if (is.null(baseline_par_tmp)) return(NULL)
      
      if (all(names(baseline_par_tmp) == all_par_names) == FALSE) {
        showNotification(switch(lang,
            "eng" = "Parameter file includes erroneous parameter names",
            "swe" = "Parameterfilen innehåller felaktiga parameternamn"
        )  
        )
        return(NULL)
      }
      
      baseline_parameters(baseline_par_tmp)
      
      ## Update baseline name
      updateTextInput(inputId = "scen_1_name", value = baseline_name)
      
      ## Create scenario id variable
      dat[, scen_id := "scen_1"]
      
      ## Update reactive value data
      dat_base[["scen_1"]] <- dat
      
      ## Update simulation year slider
      updateSliderInput(inputId = "slider_simyears",
                        value = c(dat[, min(year)], dat[, max(year)]))
      
      ## Update baseline active status
      baseline_active("yes")
    }
  }) # end observeEvent(input$button_load_scen_1, {
  
  ## Read uploaded baseline parameters and update sliders
  observeEvent(baseline_parameters(), {
    updateSliderInput(
      inputId = "cpaf_cancer",
      value = baseline_parameters()$cpaf_cancer)
    
    updateSliderInput(
      inputId = "cpaf_cvd",
      value = baseline_parameters()$cpaf_cvd)
    
    updateSliderInput(
      inputId = "age_cutoff_cancer",
      value = baseline_parameters()$age_cutoff_cancer)
    
    updateSliderInput(
      inputId = "age_cutoff_cvd",
      value = baseline_parameters()$age_cutoff_cvd)
    
  })
  
  ## Observer toggling baseline controls
  observeEvent(baseline_controls_active(), {
    if (baseline_controls_active() == FALSE) {
      for (i in 1:length(baseline_slider_names)) {
        disable(baseline_slider_names[i])
      }
    } else if (baseline_controls_active() == TRUE) {
      for (i in 1:length(baseline_slider_names)) {
        enable(baseline_slider_names[i])
      }
    }
  })
  
  
  ## Download controls
  output$button_save_baseline_data <- downloadHandler(
    filename = function() {
      paste0(input$scen_1_name, ".csv")
    },
    content = function(file) {
      fwrite(dat_base[["scen_1"]][, ..all_var_names], file, row.names = FALSE, 
             sep = ";", encoding = "UTF-8")
    }
  )
  
  output$button_save_baseline_parameters <- downloadHandler(
    filename = function() {
      paste0(input$scen_1_name, ".json")
    },
    content = function(path) {
      write_json(
        path = path,
        x = list(
          cpaf_cancer = input$cpaf_cancer,
          cpaf_cvd = input$cpaf_cvd,
          age_cutoff_cancer = input$age_cutoff_cancer,
          age_cutoff_cvd = input$age_cutoff_cvd,
          dcost_total_cancer = default_parameters$dcost_total_cancer,
          icost_total_cancer = default_parameters$icost_total_cancer,
          dcost_growth_cancer = default_parameters$dcost_growth_cancer,
          dcost_total_cvd = default_parameters$dcost_total_cvd,
          icost_total_cvd = default_parameters$icost_total_cvd,
          dcost_growth_cvd = default_parameters$dcost_growth_cvd,
          dcost_total_base_year = default_parameters$dcost_total_base_year,
          rr = default_parameters$rr,
          rr_cancer_diet = default_parameters$rr_cancer_diet,
          rr_cvd_diet = default_parameters$rr_cvd_diet,
          communalities = default_parameters$communalities
        ), auto_unbox = TRUE, pretty = TRUE
      )
    }
  )
  
  lapply(2:n_scen_max, function(scen) {
    output[[paste0("button_save_scen_", scen)]] <- downloadHandler(
      filename = function() {
        paste0(input[[paste0("scen_", scen, "_name")]], ".csv")
      },
      content = function(file) {
        fwrite(dat_base[[paste0("scen_", scen)]][, ..all_var_names], file, 
               row.names = FALSE, sep = ";", encoding = "UTF-8")
      }
    )
  })
  
  ################################################################################
  ## Dynamic UI elements
  ################################################################################
  
  ## Header of Baseline tab in the control panel. Adds checkmark when active.
  output$baseline_header <- renderUI({
    if (baseline_active() != "yes") {
      h5(switch(lang,
                "eng" = "Baseline ",
                "swe" = "Baslinje "
      )
      )
    } else if (baseline_active() == "yes") {
      switch(lang,
 "eng" = h5("Baseline ",
            ok_icon("The checkmark indicates that a baseline is active")),
 "swe" = h5("Baslinje ",
            ok_icon("Bocken indikerar aktiv baslinje"))
      )
    }
  })
  
  ## Create all scenario controls. To change number of scenarios, set n_scen_max
  ## at top of this script
  output$scenarioPanels <- renderUI({
    
    ui_scen_list <- tagList()
    
    for (scen in 2:n_scen_max) {
      ui_scen_list[[scen]] <-  
        bsCollapse(id = paste0("panel_scen_", scen),
bsCollapsePanel(paste0("(+) Scenario ", scen),
               
               h5(switch(lang, 
                         "eng" = "Load scenario data",
                         "swe" = "Ladda scenariodata"),
                  info_icon(switch(lang,
                       "eng" = "Load pre-simulated scenario data from csv-file",
                       "swe" = "Ladda försimulerade scenariodata från csv-fil")
                  )
               ),
               if (baseline_active() == "yes") {
                 fileInput(paste0("button_load_scen_", scen), label=NULL, 
                           accept=".csv", placeholder = switch(lang,
                                                   "eng" = "Choose a csv-file",
                                                   "swe" = "Välj en csv-fil")
                 )
               } else {
                 disabled(fileInput(paste0("button_load_scen_", scen), label=NULL, 
                                    accept=".csv", placeholder = switch(lang,
                                                    "eng" = "Baseline needed",
                                                    "swe" = "Baslinje behövs")
                 ))
               },
               
               h5(switch(lang, 
                         "eng" = "Scenario name",
                         "swe" = "Scenarionamn"),
                  info_icon(switch(lang,
                               "eng" = "Choose a custom name for the scenario",
                               "swe" = "Välj ett anpassat scenarionamn")
                  )
               ),
               textInput(paste0("scen_", scen, "_name"), label = NULL,
                         value = paste0("scen_", scen)),
               
               h5(switch(lang, 
                         "eng" = "Risk group size", 
                         "swe"= "Storlek på riskgrupp"), 
                  info_icon(risk_factor_info_text)
               ),
               
               bsCollapse(id = "risk-factors",
                          bsCollapsePanel(switch(lang, 
                                                 "eng" = "(+) Non-dietary", 
                                                 "swe"= "(+) Icke-kost"),
                                          
                                          risk_factor_slider(
                                            paste0("cfact_smoking_", scen),
                                            label = label_slider_smoking,
                                            info_text = info_text_smoking
                                          ),
                                          risk_factor_slider(
                                            paste0("cfact_alcohol_", scen),
                                            label = label_slider_alcohol,
                                            info_text = info_text_alcohol
                                          ),
                                          risk_factor_slider(
                                            paste0("cfact_inactivity_", scen),
                                            label = label_slider_inactivity, 
                                            info_text = info_text_inactivity 
                                          ),
                                          risk_factor_slider(
                                            paste0("cfact_obesity_", scen),
                                            label = label_slider_obesity,
                                            info_text = info_text_obesity
                                          )
                          ), # end bsCollapsePanel("Non-dietary"
                          bsCollapsePanel(switch(lang, 
                                                 "eng" = "(+) Dietary", 
                                                 "swe"= "(+) Kost"),
                                          
                                          risk_factor_slider(
                                            paste0("cfact_fruit_", scen),
                                            label = label_slider_fruit,
                                            info_text = info_text_fruit
                                          ),
                                          risk_factor_slider(
                                            paste0("cfact_wholegrains_", scen),
                                            label = label_slider_wholegrains,
                                            info_text = info_text_wholegrains
                                          ),
                                          risk_factor_slider(
                                            paste0("cfact_greens_", scen),
                                            label = label_slider_greens,
                                            info_text = info_text_greens
                                          ),
                                          risk_factor_slider(
                                            paste0("cfact_meat_", scen),
                                            label = label_slider_meat,
                                            info_text = info_text_meat
                                          ),
                                          risk_factor_slider(
                                            paste0("cfact_salt_", scen),
                                            label = label_slider_salt,
                                            info_text = info_text_salt
                                          )
                          ) # end bsCollapsePanel("Dietary"
               ), # end bsCollapse(id = "risk-factors"
               
               h5(switch(lang, 
                         "eng" = "Phase-in period", 
                         "swe"= "Infasningsperiod"),
                  info_icon(switch(lang, 
                                   "eng" = "Changes in the size of the risk groups are phased in linearly during this period and remain thereafter", 
                                   "swe"= "Förändringar i storleken på riskgrupperna fasas in linjärt under denna period och kvarstår därefter"))
               ),
                     
                     sliderInput(
                       paste0("slider_intervention_years_", scen),
                       label = NULL,
                       min = start_year_min,
                       max = end_year_max,
                       value = c(2025, 2026),
                       step = 1,
                       sep = "",
                       width = slider_width),
                     
                     h5(switch(lang,
                               "eng" = "Simulation controls",
                               "swe" = "Simuleringskontroller")
                     ),
                     
                     actionButton(paste0("button_run_scen_", scen), 
                                  switch(lang, 
                                         "eng" = "Run", 
                                         "swe"= "Starta")),
                     
                     actionButton(paste0("button_clear_scen_", scen), 
                                  switch(lang, 
                                         "eng" = "Stop/clear", 
                                         "swe"= "Stoppa/rensa")),
                     
                     ## Text output showing simulation progress
                     h5(switch(lang, 
                               "eng" = "Log", 
                               "swe"= "Logg")
                     ),
                     
                     tags$small(verbatimTextOutput(paste0("log_text_", scen))),
                     disabled(downloadButton(paste0("button_save_scen_", scen),
                                             switch(lang, 
                                                    "eng" = "Save data", 
                                                    "swe"= "Spara data")))
     )  # end bsCollapse
        ) # end bsCollapsePanel
    }
    
    ui_scen_list
    
  }) # end renderUI
  
  ## Slider for changing year of the cross section graph
  output$slider_year_cs_plot <- renderUI({
    dat_base_1 <- dat_base[["scen_1"]]
    if (is.null(dat_base_1)) return(NULL)
    tagList(
      h5(switch(lang, 
                "eng" = "Year to plot",
                "swe" = "År i diagram")
      ),
      
      sliderInput("slider_year_cs_plot",
                  label = NULL,
                  min = input$slider_simyears[1],
                  max = max(input$slider_simyears[1],
                            dat_base_1[!is.na(s_cancer), unique(year)]),
                  value = max(input$slider_simyears[1],
                              dat_base_1[!is.na(s_cancer), unique(year)]),
                  step = 1,
                  sep = "")
    ) # tagList
  })
  
  ## Slider for changing years to aggregate in accumulated graph
  output$slider_years_acc_plot <- renderUI({
    dat_base_1 <- dat_base[["scen_1"]]
    if (is.null(dat_base_1)) return(NULL)
    tagList(
      h5(switch(lang, 
                "eng" = "Years to aggregate",
                "swe" = "År att aggregera")
      ),
      
      sliderInput("slider_years_acc_plot",
                  label = NULL,
                  min = input$slider_simyears[1],
                  max = input$slider_simyears[2],
                  value = c(input$slider_simyears[1], input$slider_simyears[2]),
                  step = 1,
                  sep = ""),
    ) # tagList
  })
  
  ## Buttons for saving graph data
  
  output$render_button_save_ts_data <- renderUI({
    if (is.null(dat_base[["scen_1"]])) return(NULL)
    downloadButton("button_save_ts_data", switch(lang, 
                                                 "eng" = "Save plot data",
                                                 "swe" = "Spara diagramdata"))
  })
  
  output$render_button_save_cs_data <- renderUI({
    if (is.null(dat_base[["scen_1"]])) return(NULL)
    downloadButton("button_save_cs_data", switch(lang, 
                                                 "eng" = "Save plot data",
                                                 "swe" = "Spara diagramdata"))
  })
  
  output$render_button_save_acc_data <- renderUI({
    if (is.null(dat_base[["scen_1"]])) return(NULL)
    downloadButton("button_save_acc_data", switch(lang, 
                                                  "eng" = "Save plot data",
                                                  "swe" = "Spara diagramdata"))
  })
  
  ##############################################################################
  ## Simulation controls
  ##############################################################################
  
  ## Observers for running, clearing and loading scenarios
  
  lapply(1:n_scen_max, function(scen) {
    observeEvent(input[[paste0("button_run_scen_", scen)]], {
      run_scen(scen, input = input, sim_status = sim_status, 
               baseline_active = baseline_active)
    })
  })
  
  lapply(1:n_scen_max, function(scen) {
    observeEvent(input[[paste0("button_clear_scen_", scen)]], {
      clear_scen(scen, sim_status = sim_status, dat_base = dat_base, 
                 baseline_active = baseline_active)
    })
  })
  
  lapply(2:n_scen_max, function(scen) {
    observeEvent(input[[paste0("button_load_scen_", scen)]], {
    load_scen(scen, input = input, sim_status = sim_status, dat_base = dat_base)
    })
  })
  
  ## Stop all runs and delete data on exit
  onStop(function() {
    for (scen in 1:n_scen_max) {
      scen_path <-  paste0(PROJECTROOT, "/UI/runs/scen_", scen, "/")
      cat("stop", file = paste0(scen_path, "ui_status.txt"))
      if (length(grep(".rda", dir(paste0(scen_path, "data"))))>0) {
        system(paste0("rm ", scen_path, "data/*.rda"))
      } # end if
    } # end for
  })
  
  ##############################################################################
  ## Plots
  ##############################################################################
  
  output$plot_time_series <- renderPlotly({
    
    ## Extract reactive data object as a normal list
    dat_list <- reactiveValuesToList(dat_base)
    
    ## Plot only if there is some baseline data
    if (is.null(dat_list[["scen_1"]])) return(NULL)
    
    create_plot(dat_list = dat_list,
                plot_type = "ts",
                var_to_plot = input$ts_var_to_plot,
                group = input$ts_group, 
                measure = input$ts_measure,
                years = input$slider_simyears,
                var_list = var_list,
                var_labels = var_labels,
                input = input)
  })
  
  output$plot_cross_section <- renderPlotly({
    
    ## Extract reactive data object as a normal list
    dat_list <- reactiveValuesToList(dat_base)
    
    ## Plot only if there is some baseline data and years slider loaded
    if (is.null(dat_list[["scen_1"]])) return(NULL)
    if (is.null(input$slider_year_cs_plot)) return(NULL)
    
    create_plot(dat_list = dat_list,
                plot_type = "cs",
                var_to_plot = input$cs_var_to_plot,
                group = input$cs_group, 
                measure = input$cs_measure,
                years = input$slider_year_cs_plot,
                var_list = var_list,
                var_labels = var_labels,
                input = input)
  })
  
  output$plot_accumulated <- renderPlotly({
    
    ## Extract reactive data object as a normal list
    dat_list <- reactiveValuesToList(dat_base)
    
    ## Plot only if there is some baseline data and years slider loaded
    if (is.null(dat_list[["scen_1"]])) return(NULL)
    if (is.null(input$slider_years_acc_plot)) return(NULL)
    
    create_plot(dat_list = dat_list,
                plot_type = "acc",
                var_to_plot = input$acc_var_to_plot,
                group = input$acc_group, 
                measure = input$acc_measure,
                years = input$slider_years_acc_plot,
                var_list = var_list,
                var_labels = var_labels,
                input = input)
  })
  
  ## Handlers for downloading plot data
  
  output$button_save_ts_data <- downloadHandler(
    filename = function() {
      paste0("plot_data.csv")
    },
    content = function(file) {
      fwrite(create_plot(dat_list = reactiveValuesToList(dat_base),
                         plot_type = "ts",
                         var_to_plot = input$ts_var_to_plot,
                         group = input$ts_group,
                         measure = input$ts_measure,
                         years = input$slider_simyears,
                         var_list = var_list,
                         var_labels = var_labels,
                         input = input,
                         data_only = TRUE),
             file, row.names = FALSE, sep = ";", encoding = "UTF-8")
    }
  )
  
  output$button_save_cs_data <- downloadHandler(
    filename = function() {
      paste0("plot_data.csv")
    },
    content = function(file) {
      fwrite(create_plot(dat_list = reactiveValuesToList(dat_base),
                         plot_type = "cs",
                         var_to_plot = input$cs_var_to_plot,
                         group = input$cs_group,
                         measure = input$cs_measure,
                         years = input$slider_year_cs_plot,
                         var_list = var_list,
                         var_labels = var_labels,
                         input = input,
                         data_only = TRUE),
             file, row.names = FALSE, sep = ";", encoding = "UTF-8")
    }
  )
  
  output$button_save_acc_data <- downloadHandler(
    filename = function() {
      paste0("plot_data.csv")
    },
    content = function(file) {
      fwrite(create_plot(dat_list = reactiveValuesToList(dat_base),
                         plot_type = "acc",
                         var_to_plot = input$acc_var_to_plot,
                         group = input$acc_group,
                         measure = input$acc_measure,
                         years = input$slider_years_acc_plot,
                         var_list = var_list,
                         var_labels = var_labels,
                         input = input,
                         data_only = TRUE),
             file, row.names = FALSE, sep = ";", encoding = "UTF-8")
    }
  )
  
} # End server

## Build app
shinyApp(ui = ui, server = server)
