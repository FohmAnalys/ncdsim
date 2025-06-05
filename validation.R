##*****************************************************************************
##* Script for validation of model output:
##*  - Demographic validation using data from Statistics Sweden (observed
##*    outcomes + forecast)
##*  - Validation of stocks and flows
##*    - Graphical presentation
##*    - Consistency control between stocks and flows
##* - If debug information from ode() exists a cohort based validation is done 
##*   to evaluate how many time steps that are needed.
##*****************************************************************************

# Packages are loaded
library(data.table)
library(ggplot2)
library(gridExtra)
library(scales)
library(stringr)


##* FUNCTION DEFINITIONS ******************************************************

fagegrp10 <- function(age) {
  agegrp = cut(age, 
               breaks = c(seq(0, 80, by = 10), Inf), 
               right = FALSE)
  return(agegrp)
}
# table(age = 0:100, agegrp = fagegrp10(0:100))

fsex <- function(sex) {
  vret <- rep("", length = length(sex))
  vret[sex == 1] <- "Men"
  vret[sex == 2] <- "Women"
  return(vret)
}
# table(sex = 1:2, fsex = fsex(1:2))

##* END FUNCTION DEFINITIONS **************************************************


##*****************************************************************************
##* validate_ncdsim(): general validation of demographics, stocks and flows.
##*****************************************************************************

validate_ncdsim <- function(
  projectroot = getwd(), # path to NCDSim project
  timestamp = NA,  # timestamp for simulation to validate
  use_scb_demographics = TRUE) {     

  ##****************************************************************************
  ## Reading simulation data
  ##****************************************************************************
  simdat <- fread(file = paste0(projectroot, "/Output/output_NCDSim_", 
                                timestamp,".csv"))

  firstyear <- min(unique(simdat$year))
  lastyear <- max(unique(simdat$year))

  ##****************************************************************************
  ## Reading validation data
  ##****************************************************************************
  if(use_scb_demographics){
    scbpop <- fread(file = paste0(projectroot, "/Input/pop_counts_scb.csv"))
  }
  

  ##****************************************************************************
  ## Reading alignment data
  ##****************************************************************************
  
  input_path <- paste0(projectroot, "/Input/")
  
  # Incidences and prevalences
  incid_new <- fread(paste0(input_path, "incidenser.csv"), sep = ";")
  preval_new <- fread(paste0(input_path, "prevalenser.csv"), sep = ";")
  colnames(incid_new) <- c("year","age", "sex", "cvd", "cancer")
  colnames(preval_new) <- c("year","age", "sex", "cvd", "cancer", "totpop")
  
  incid <- CJ(year = min(incid_new[,(year)]):max(incid_new[, (year)]), 
              sex = 1:2, age = 0:100)
  preval <- CJ(year = min(preval_new[,(year)]):max(preval_new[, (year)]), 
               sex = 1:2, age = 0:100)
  
  incid <- merge(incid, incid_new, by = c("year", "age", "sex"), all = TRUE)
  preval <- merge(preval, preval_new[,c("year", "age", "sex", "cvd", "cancer")], 
                  by = c("year","age", "sex"), all = TRUE)
  
  incid <- melt(incid, id.vars = c("sex", "age", "year"), 
                value.vars = c("cvd", "cancer"))
  preval <- melt(preval, id.vars = c("sex", "age", "year"), 
                 value.vars = c("cvd", "cancer"))
  colnames(incid) <- c("sex", "age", "year", "NCD", "incidence")
  colnames(preval) <- c("sex", "age", "year", "NCD", "prevalence")
  incid <- setnafill(incid, "const", fill = 0, cols = c("incidence")) 
  preval <- setnafill(preval, "const", fill = 0, cols = c("prevalence")) 
  
  prevalence_cancer <- preval[NCD == "cancer", .(year, sex, age, prevalence)]
  prevalence_cvd <- preval[NCD == "cvd", .(year, sex, age, prevalence)]
  incidence_cancer <- incid[NCD == "cancer", .(year, sex, age, incidence)]
  incidence_cvd <- incid[NCD == "cvd" & year > 2009, 
                         .(year, sex, age, incidence)]
  
  ##****************************************************************************
  ## Graphical presentation of simulation results: stocks and flows
  ##****************************************************************************
  
  #  Initialize pdf output
  pdf(file = paste0(projectroot, "/output/validation_ncdsim_", timestamp, 
                    ".pdf"), 
      paper = "a4r", width = 10, height = 6, pointsize = 8)          
  devpdf <- dev.cur()
  
  ##********************************************
  ## Population frequencies by year, sex and age
  ##********************************************
  
  t1 <- simdat[, .(pop = sum(s_pop + s_cancer + s_cvd), lbl = "NCDSim"), 
               by = .(year, sex, age)]
  if (use_scb_demographics) {
    t2 <- scbpop[, .(pop = sum(pop), lbl = "SCB"), by = .(year, sex, age)]
    totdat <- rbindlist(list(t1, t2), use.names = TRUE)
  } else {
    totdat <- t1
  }

  for (y in firstyear:lastyear){
    pdat <- totdat[year == y]
    p <- ggplot(data = pdat, aes(y = pop, x = age, color = lbl)) +
      geom_line(size = 1) +
      facet_wrap(~fsex(sex)) +
      ggtitle(paste0("Simulated and observed/assumed population in ", y)) +
      ylim(0, NA) +
      theme(legend.position = "bottom")
    print(p)
  }
  
  ##********************************************
  ## Total population by year
  ##********************************************
  
  t1 <- simdat[, .(pop = sum(s_pop + s_cancer + s_cvd), lbl = "NCDSim"), 
               by = .(year)]
  if (use_scb_demographics) {
    t2 <- scbpop[, .(pop = sum(pop), lbl = "SCB"), by = .(year)]
    totdat <- rbindlist(list(t1, t2), use.names = TRUE)
  } else {
    totdat <- t1
  }
  
  p <- ggplot(data = totdat[year < lastyear], 
              aes(y = pop, x = year, color = lbl)) +
    geom_line(size = 1) +
    ggtitle("Total population per year ") +
    ylim(0, NA) +
    theme(legend.position = "bottom")
  print(p)
  
  
  ##********************************************
  ## Demographic components by year
  ##********************************************

  t1 <- simdat[, .(ndead = sum(f_dead + f_cancer_dead + f_cvd_dead),
                   nborn = sum(f_born),
                   nimmig = sum(f_immig_pop + f_immig_cvd + f_immig_cancer),
                   nemig = sum(f_emig_pop + f_emig_cvd + f_emig_cancer),
                   lbl = "NCDSim"), 
               by = year]
  if (use_scb_demographics) {
    t2 <- scbpop[, .(ndead = sum(dead), nborn = sum(born), nimmig = sum(immig),
                     nemig = sum(emig), lbl = "SCB"), 
                 by = year]
    totdat <- rbindlist(list(t1, t2), use.names = TRUE)
  } else {
    totdat <- t1
  }
  
  pdat <- melt(totdat[year >= (firstyear + 1) & year <= lastyear], 
               id.vars = c("year", "lbl"))
  
  p <- ggplot(data = pdat, aes(y = value, x = year, color = lbl)) +
    geom_line(size = 1) +
    facet_wrap(~variable) +
    ggtitle("Simulated and observed births, deaths and migration") +
    ylim(0, NA) +
    theme(legend.position = "bottom")
  print(p)
  
  
  ##*********************************************
  ## Demographic components by year and age group
  ##*********************************************
  
  t1 <- simdat[, .(ndead = sum(f_dead + f_cancer_dead + f_cvd_dead),
                   nborn = sum(f_born),
                   nimmig = sum(f_immig_pop + f_immig_cvd + f_immig_cancer),
                   nemig = sum(f_emig_pop + f_emig_cvd + f_emig_cancer),
                   lbl = "NCDSim"), 
               by = .(year, agegrp = fagegrp10(age))]
  if (use_scb_demographics) {
    t2 <- scbpop[, .(ndead = sum(dead), nborn = sum(born), nimmig = sum(immig),
                     nemig = sum(emig), lbl = "SCB"), 
                 by = .(year, agegrp = fagegrp10(age))]
    totdat <- rbindlist(list(t1, t2), use.names = TRUE)
  } else {
    totdat <- t1
  }
  
  pdat <- melt(totdat[year >= (firstyear + 1) & year <= lastyear], 
               id.vars = c("year", "lbl", "agegrp"))
  
  for (v in setdiff(unique(pdat$variable), "nborn")) {
    p <- ggplot(data = pdat[variable == v], 
                aes(y = value, x = year, color = lbl)) +
      geom_line(size = 1) +
      facet_wrap(~agegrp) +
      ggtitle(paste0("Simulated and observed values by agegoup, variable: ", 
                     v)) +
      ylim(0, NA) +
      theme(legend.position = "bottom")
    print(p)
  }
  

  ##*********************************************
  ## Flows by year, sex and age
  ##*********************************************

  # pdat <- melt(simdat[year >= (firstyear + 1), 
  #                     .(year, sex, age, f_dead, f_cvd_dead, f_cancer_dead,
  #                       f_immig_pop, f_immig_cvd,f_immig_cancer, f_emig_pop,
  #                       f_emig_cvd, f_emig_cancer)],
  #              id.vars = c("year", "sex", "age"))
  # for (y in unique(pdat$year)) {
  #   p <- ggplot(data = pdat[year == y], 
  #               aes(y = value, x = age, color = variable)) +
  #     geom_line(size = 1) +
  #     facet_wrap(~sex) +
  #     ggtitle(paste0("Population flows, year: ", y))
  #   print(p)
  # }

  
  ##*********************************************
  ## Demographic flows by year
  ##*********************************************
    
  tmp <- melt(simdat[year >= (firstyear + 1), 
                     .(year, f_born, f_dead, f_cvd_dead, f_cancer_dead, 
                       f_immig_pop, f_immig_cvd, f_immig_cancer, f_emig_pop, 
                       f_emig_cvd, f_emig_cancer)], 
              id.vars = "year")
  pdat <- tmp[, .(sumval = sum(value)), by = .(year, variable)]
  p <- ggplot(data = pdat, aes(y = sumval, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0, NA) +
    ggtitle("Number of births, deaths, emigrations and immigrations, by year")
  print(p)  

  
  ##*********************************************
  ## NCD stocks by year
  ##*********************************************
  
  tmp <- simdat[, .(ncancer = sum(s_cancer), ncvd = sum(s_cvd)), 
             by = year]
  pdat <- melt(tmp, id.vars = "year")
  tmp <- prevalence_cancer[, .(variable = "ncancer_sos",
                               value = sum(prevalence)), by = year]
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  tmp <- prevalence_cvd[, .(variable = "ncvd_sos",
                            value = sum(prevalence)), by = year]
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0, NA) +
    ggtitle("Stocks of cancer and CVD (and validation data from SoS), by year")
  print(p)  

  
  ##*********************************************
  ## Total NCD stocks by year
  ##*********************************************
  
  tmp <- simdat[, .(nstock = sum(s_cancer) + sum(s_cvd)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  tmp <- prevalence_cancer[, .(variable = "nstock_sos",
                               value = sum(prevalence)), by = year]
  
  tmp1 <- prevalence_cvd[, .(variable = "nstock_sos",
                            value1 = sum(prevalence)), by = year]
  sosd <- merge(tmp,tmp1[,c("year","value1")], by = c("year"), all = TRUE)
  sosd[,value := value + value1]
  sosd[,value1 := NULL]
  
  pdat <- rbindlist(list(pdat, sosd), use.names = TRUE)
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0,NA) + 
    ggtitle("Total stock of diseases (and validation data from SoS), by year")
  print(p)  
  
  
  ##*********************************************
  ## NCD stocks by year, sex and age group
  ##*********************************************

  tmp <- simdat[, .(ncancer = sum(s_cancer), ncvd = sum(s_cvd)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  tmp <- prevalence_cancer[, .(variable = "ncancer_sos",
                               value = sum(prevalence)), 
                           by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  tmp <- prevalence_cvd[, .(variable = "ncvd_sos",
                            value = sum(prevalence)), 
                        by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  p <- ggplot(data = pdat[sex == 1], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Stocks of cancer and CVD (and validation data from SoS), ", 
                   "by year and agegroup, men"))
  print(p)  

  p <- ggplot(data = pdat[sex == 2], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Stocks of cancer and CVD (and validation data from SoS), ", 
                   "by year and agegroup, women"))
  print(p)  
  
  
  ##*********************************************
  ## Total stocks by year, sex and age group
  ##*********************************************

  tmp <- simdat[, .(nstock = sum(s_cancer) + sum(s_cvd)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  tmp <- prevalence_cancer[, .(variable = "nstock_sos",
                               value = sum(prevalence)), 
                           by = .(year, sex, agegrp = fagegrp10(age))]
  tmp1 <- prevalence_cvd[, .(variable = "ncvd_sos",
                            value1 = sum(prevalence)), 
                        by = .(year, sex, agegrp = fagegrp10(age))]
  sosd <- merge(tmp, tmp1[,c("year","sex", "agegrp", "value1")], 
                by = c("year","sex", "agegrp"), all = TRUE)
  sosd[, value := value + value1]
  sosd[, value1 := NULL]
  
  pdat <- rbindlist(list(pdat, sosd), use.names = TRUE)
  p <- ggplot(data = pdat[sex == 1], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Total stock of diseases (and validation data from SoS), ", 
                   "by year and agegroup, men"))
  print(p)  
  
  p <- ggplot(data = pdat[sex == 2], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Total stock of diseases (and validation data from SoS), ", 
                   "by year and agegroup, women"))
  print(p)  
  

  ##*********************************************
  ## NCD flows by year
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(f_pop_cancer = sum(f_pop_cancer),
                  f_cancer_pop = sum(f_cancer_pop),
                  f_cancer_dead = sum(f_cancer_dead),
                  f_pop_cvd = sum(f_pop_cvd),
                  f_cvd_pop = sum(f_cvd_pop),
                  f_cvd_dead = sum(f_cvd_dead)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  tmp <- incidence_cancer[, .(variable = "f_pop_cancer_sos",
                              value = sum(incidence)), by = year]
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  tmp <- incidence_cvd[, .(variable = "f_pop_cvd_sos",
                           value = sum(incidence)), by = year]
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  
  p <- ggplot(data = pdat[variable %like% "cancer"], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0, NA) +
    ggtitle("Flows related to cancer (and validation data from SoS), by year")
  print(p)  

  p <- ggplot(data = pdat[variable %like% "cvd"], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0, NA) +
    ggtitle("Flows related to cvd (and validation data from SoS), by year")
  print(p)  
  
  
  ##*********************************************
  ## Total NCD flows by year
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(flow_disease = sum(f_pop_cancer) + sum(f_pop_cvd),
                  flow_healty = sum(f_cancer_pop) + sum(f_cvd_pop),
                  flow_dead = sum(f_cancer_dead) + sum(f_cvd_dead)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")

  tmp <- incidence_cancer[, .(variable = "flow_disease_sos",
                              value = sum(incidence)), 
                          by = year]
  tmp1 <- incidence_cvd[, .(variable = "f_pop_cvd_sos",
                            value1 = sum(incidence)), 
                        by = year]
  sosd <- merge(tmp, tmp1[, c("year", "value1")], by = c("year"), all = FALSE)
  sosd[, value := value + value1]
  sosd[, value1 := NULL]
  
  pdat <- rbindlist(list(pdat, sosd), use.names = TRUE)
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ggtitle("Total flows (and validation data from SoS), by year")
  print(p)  
  
  
  ##*********************************************
  ## NCD flows by year, sex and age group
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(f_pop_cancer = sum(f_pop_cancer),
                  f_cancer_pop = sum(f_cancer_pop),
                  f_cancer_dead = sum(f_cancer_dead),
                  f_pop_cvd = sum(f_pop_cvd),
                  f_cvd_pop = sum(f_cvd_pop),
                  f_cvd_dead = sum(f_cvd_dead)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  tmp <- incidence_cancer[, .(variable = "f_pop_cancer_sos",
                              value = sum(incidence)), 
                          by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  tmp <- incidence_cvd[, .(variable = "f_pop_cvd_sos",
                           value = sum(incidence)), 
                       by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  
  p <- ggplot(data = pdat[sex == 1 & variable %like% "cancer"], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to cancer (and validation data from SoS), ",
                   "by year, men"))
  print(p)  

  p <- ggplot(data = pdat[sex == 2 & variable %like% "cancer"], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to cancer (and validation data from SoS), ",
                   "by year, women"))
  print(p)  

  p <- ggplot(data = pdat[sex == 1 & variable %like% "cvd"], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to CVD (and validation data from SoS), ", 
                   "by year, men"))
  print(p)  
  
  p <- ggplot(data = pdat[sex == 2 & variable %like% "cvd"], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to CVD (and validation data from SoS), ", 
                   "by year, women"))
  print(p)  
  
  
  ##*********************************************
  ## total flows by year, sex and age group
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1),
                .(flow_disease = sum(f_pop_cancer) + sum(f_pop_cvd),
                  flow_healty = sum(f_cancer_pop) + sum(f_cvd_pop),
                  flow_dead = sum(f_cancer_dead) + sum(f_cvd_dead)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  tmp <- incidence_cancer[, .(variable = "f_pop_cancer_sos",
                              value = sum(incidence)), 
                          by = .(year, sex, agegrp = fagegrp10(age))]
  tmp1 <- incidence_cvd[, .(variable = "f_pop_cvd_sos",
                           value1 = sum(incidence)), 
                       by = .(year, sex, agegrp = fagegrp10(age))]
  sosd <- merge(tmp, tmp1[, c("year", "sex", "agegrp", "value1")], 
                by = c("year", "sex", "agegrp"), all = TRUE)
  sosd[, value := value + value1]
  sosd[, value1 := NULL]
  
  pdat <- rbindlist(list(pdat, sosd), use.names = TRUE)
  
  p <- ggplot(data = pdat[sex == 1], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle("Total flows (and validation data from SoS), by year, men")
  print(p)  
  
  p <- ggplot(data = pdat[sex == 2], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle("Total flows (and validation data from SoS), by year, women")
  print(p)  

  
  ##***********************************************
  ## NCD Costs by year
  ##***********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost_cancer = sum(dcost_cancer), 
                  cost_cvd = sum(dcost_cvd)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0, NA) +
    ggtitle("Cost of cancer and CVD in SEK, by year")
  print(p)  
  
  
  ##***********************************************
  ## NCD Costs by year, indexed to 2024
  ##***********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost_cancer = sum(dcost_cancer), 
                  cost_cvd = sum(dcost_cvd)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  nd_cancer = pdat[year == 2024 & variable == "cost_cancer", value]
  nd_cvd = pdat[year == 2024 & variable == "cost_cvd", value]
  pdat[variable == "cost_cancer", value := (value / nd_cancer) * 100]
  pdat[variable == "cost_cvd", value := (value / nd_cvd) * 100]

  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ggtitle("Cost of cancer and CVD compared to 2024 in percent, by year")
  print(p) 
  
  
  ##***********************************************
  ## Total Costs by year
  ##***********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost = sum(dcost_cancer) + sum(dcost_cvd)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0, NA) +
    ggtitle("Total cost in SEK, by year")
  print(p) 
  
  
  ##***********************************************
  ## Total Costs by year, indexed to 2024
  ##***********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost = sum(dcost_cancer) + sum(dcost_cvd)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  nd = pdat[year == 2024,value]
  pdat[,value := (value / nd) * 100]
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ggtitle("Total cost compared to 2024 in percent, by year")
  print(p) 
  
  
  ##*********************************************
  ## NCD costs by year, sex and age group
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost_cancer = sum(dcost_cancer), 
                  cost_cvd = sum(dcost_cvd)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))

  p <- ggplot(data = pdat[sex == 1], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp) +
    ylim(0, NA) +
    ggtitle("Costs of cancer and CVD, by year and agegroup, men")
  print(p)  
  
  p <- ggplot(data = pdat[sex == 2], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp) +
    ylim(0, NA) +
    ggtitle("costs of cancer and CVD in SEK, by year and agegroup, women")
  print(p)  
  
  
  ##*********************************************
  ## NCD costs by year, sex and age group, indexed to 2024
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost_cancer = sum(dcost_cancer), 
                  cost_cvd = sum(dcost_cvd)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  
  nd_cancer_m = pdat[year == 2024 & variable == "cost_cancer" & sex == 1, value]
  nd_cancer_w = pdat[year == 2024 & variable == "cost_cancer" & sex == 2, value]
  nd_cvd_m = pdat[year == 2024 & variable == "cost_cvd" & sex == 1, value]
  nd_cvd_w = pdat[year == 2024 & variable == "cost_cvd" & sex == 2, value]
  pdat[variable == "cost_cancer" & sex == 1, 
       value := (value / nd_cancer_m) * 100]
  pdat[variable == "cost_cancer" & sex == 2, 
       value := (value / nd_cancer_w) * 100]
  pdat[variable == "cost_cvd" & sex == 1, value := (value / nd_cvd_m) * 100]
  pdat[variable == "cost_cvd" & sex == 2, value := (value / nd_cvd_w) * 100]
  
  p <- ggplot(data = pdat[sex == 1], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Costs of cancer and CVD compared to 2024, by year and ", 
                   "agegroup, men"))
  print(p)  
  
  p <- ggplot(data = pdat[sex == 2], 
              aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("costs of cancer and CVD compared to 2024, by year and ", 
                   "agegroup, women"))
  print(p)  
  
  
  ##*********************************************
  ## total costs by year, sex and age group
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost = sum(dcost_cancer) + sum(dcost_cvd)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  
  p <- ggplot(data = pdat, aes(y = value, x = year, color = fsex(sex))) +
    geom_line(size = 1) +
    facet_wrap(~agegrp) +
    ylim(0, NA) +
    ggtitle("Total costs in SEK, by year and agegroup")
  print(p)  
  
  
  ##*********************************************
  ## total costs by year, sex and age group, indexed to 2024
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost = sum(dcost_cancer) + sum(dcost_cvd)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  
  nd_m = pdat[year == 2024 & sex == 1, value]
  nd_w = pdat[year == 2024 & sex == 2, value]
  pdat[sex == 1, value := (value/nd_m) * 100]
  pdat[sex == 2, value := (value/nd_w) * 100]
  
  p <- ggplot(data = pdat, aes(y = value, x = year, color = fsex(sex))) +
    geom_line(size = 1) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle("Total costs compared to 2024, by year and agegroup")
  print(p)  
  
  
  ##*********************************************
  ## Death rates by year, sex, age and group
  ##*********************************************
  
  tmp <- simdat[year %% 2 == 0, .(year, sex, age, dr, dr_cancer, dr_cvd)]
  pdat <- melt(tmp, id.vars = c("year", "sex", "age"))
  
  p <- ggplot(data = pdat[sex == 1], 
              aes(y = value, x = age, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~year) +
    ylim(0, NA) +
    ggtitle("Death rates by year, age and group, men")
  print(p)  
  
  p <- ggplot(data = pdat[sex == 2], 
              aes(y = value, x = age, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~year) +
    ylim(0, NA) +
    ggtitle("Death rates by year, age and group, women")
  print(p)  
  
  
  ##*********************************************
  ## Log of Death rates by year, sex, age and group
  ##*********************************************
  
  tmp <- simdat[year %% 2 == 0, .(year, sex, age, dr, dr_cancer, dr_cvd)]
  pdat <- melt(tmp, id.vars = c("year", "sex", "age"))
  
  p <- ggplot(data = pdat[sex == 1], 
              aes(y = log10(value), x = age, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~year) +
    ylim(NA, 0) +
    ggtitle("Death rates by year, age and group, men, log scale")
  print(p)  
  
  p <- ggplot(data = pdat[sex == 2], 
              aes(y = log10(value), x = age, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~year) +
    ylim(NA, 0) +
    ggtitle("Death rates by year, age and group, women, log scale")
  print(p) 
  
  
  ##********************************************
  ## Risk factor prevalence by year, sex and age 
  ##********************************************
  
  tmp <- simdat[, .(year, sex, age, prev_alcohol, prev_obesity, prev_smoking, 
                    prev_inactivity)]
  pdat <- melt(tmp, id.vars = c("year", "sex", "age"))
  
  for (y in seq(min(simdat$year), max(simdat$year), by = 2)) {
    p <- ggplot(data = pdat[year == y], 
                aes(y = value, x = age, color = variable)) +
      geom_line(size = 1) +
      facet_wrap(~fsex(sex)) +
      ylim(0, NA) +
      ggtitle(paste0("Risk factor prevalences by sex, age and group, year: ", 
                     y))
    print(p)  
  }
  
  cl <- c("year", "sex", "age", "prev_fruit", "prev_wholegrains", "prev_greens", 
          "prev_meat", "prev_salt")
  tmp <- simdat[, ..cl]
  pdat <- melt(tmp, id.vars = c("year", "sex", "age"))
  
  p <- ggplot(data = pdat[year == min(simdat$year)+2], 
              aes(y = value, x = age, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~fsex(sex)) +
    ylim(0, NA) +
    ggtitle(paste0("Risk factor prevalences by sex, age and group for food"))
  print(p)  
    
  
  ##********************************************
  ## PAF by year, sex and age 
  ##********************************************
  
  tmp <- simdat[, .(year, sex, age, paf_cvd_smoking, paf_cvd_inactivity,
                    paf_cvd_obesity, paf_cvd_alcohol, 
                    paf_cancer_smoking, paf_cancer_inactivity, 
                    paf_cancer_obesity, paf_cancer_alcohol, cpaf_diet_cancer,
                    cpaf_diet_cvd)]
  pdat <- melt(tmp, id.vars = c("year", "sex", "age"))
  pdat[, grp := fifelse(length(grep("cancer", variable)) > 0, "cancer", "CVD"),
       by = .(year, sex, age, variable)]
  
  for (y in seq(min(simdat$year), max(simdat$year), by = 2)) {
    p <- ggplot(data = pdat[year == y], 
                aes(y = value, x = age, color = variable)) +
      geom_line(size = 1) +
      facet_wrap(~fsex(sex) + grp) +
      ylim(0, NA) +
      ggtitle(paste0("PAF by sex, age and group, year: ", y))
    print(p)  
  }
  
  
  ##********************************************
  ## Food related PAF by year, sex and age 
  ##********************************************
  
  cl <- c("year", "sex", "age", "paf_cancer_fruit", "paf_cancer_wholegrains", 
          "paf_cancer_greens", "paf_cancer_meat", "paf_cancer_salt", 
          "paf_cvd_fruit", "paf_cvd_wholegrains", "paf_cvd_greens", 
          "paf_cvd_meat", "paf_cvd_salt" )
  tmp <- simdat[, ..cl]
  pdat <- melt(tmp, id.vars = c("year", "sex", "age"))
  pdat[, grp := fifelse(length(grep("cancer", variable)) > 0, "cancer", "CVD"),
       by = .(year, sex, age, variable)]
  p <- ggplot(data = pdat[year == min(simdat$year)+2], 
              aes(y = value, x = age, color = variable)) +
    geom_line(size = 1) +
    facet_wrap(~fsex(sex) + grp) +
    ylim(0, NA) +
    ggtitle(paste0("Food related PAF by age and group"))
  print(p)  
 
  
  ##********************************************
  ## PAF from other riskfactors by year, sex and age 
  ##********************************************
  
  tmp <- simdat[, .(year, sex, age,
                    paf_cancer_other, paf_cvd_other)]
  pdat <- melt(tmp, id.vars = c("year", "sex", "age"))
  
  for (y in seq(min(simdat$year), max(simdat$year), by = 2)) {
    p <- ggplot(data = pdat[year == y], 
                aes(y = value, x = age, color = variable)) +
      geom_line(size = 1) +
      facet_wrap(~fsex(sex)) +
      ylim(0, NA) +
      ggtitle(paste0("PAF from other riskfactors by sex and age, year: ", y))
    print(p)  
  } 
  
  
  ##********************************************
  ## Incidence rates
  ##********************************************
  
  tmp <- simdat[, .(year, sex, age, p_pop_cancer, p_pop_cvd)]
  pdat <- melt(tmp, id.vars = c("year", "sex", "age"))

  for (y in seq(min(simdat$year), max(simdat$year), by = 2)) {
    p <- ggplot(data = pdat[year == y], 
                aes(y = value, x = age, color = variable)) +
      geom_line(size = 1) +
      facet_wrap(~fsex(sex)) +
      ylim(0, NA) +
      ggtitle(paste0("Incidence rate by sex, age and group, year: ", y))
    print(p)  
  }
  
  
  ##****************************************************************************
  ## Checking the consistency between stocks and flows. Comparison between the 
  ## stocks returned from ode() and stocks calculated from lagged values of 
  ## stocks and current values of flows.
  ## Also: checking for negative values of "control stocks".
  ##****************************************************************************
  
  simdat[, cohort := year - age]
  
  # Sort before using lagged values
  setkey(simdat, cohort, sex, age)
  simdat[age < 100, ":="(
    s_pop_ = shift(s_pop) + f_born + f_immig_pop - f_dead - f_emig_pop + 
      f_cancer_pop + f_cvd_pop - f_pop_cancer - f_pop_cvd,
    s_cancer_ = shift(s_cancer) + f_immig_cancer - f_emig_cancer - 
      f_cancer_dead + f_pop_cancer - f_cancer_pop,
    s_cvd_ = shift(s_cvd) + f_immig_cvd - f_emig_cvd - f_cvd_dead + f_pop_cvd - 
      f_cvd_pop
  )]
  simdat[year == firstyear | cohort != shift(cohort) | sex != shift(sex) |
           age != (shift(age) + 1) | age == 100, ":="(
             # NA for non-valid lags.
             s_pop_ = NA,
             s_cancer_ = NA,
             s_cvd_ = NA
           )]
  
  simdat[, ":="(
    diff_s_pop = s_pop_ - s_pop,
    rdiff_s_pop = s_pop_ / s_pop,
    diff_s_cancer = s_cancer_ - s_cancer,
    rdiff_s_cancer = s_cancer_ / s_cancer,
    diff_s_cvd = s_cvd_ - s_cvd,
    rdiff_s_cvd = s_cvd_ / s_cvd
  )]
  
  p <- ggplot(data = simdat[year > firstyear & cohort %% 5 == 0 & age <= 100], 
              aes(y = diff_s_pop, x = age, color = factor(cohort))) +
    geom_line(size = 1) + 
    facet_wrap(~fsex(sex)) +
    labs(title = paste0("Difference between healthy population calculated ", 
                        "from flows or directly simulated"))
  print(p)
  
  p <- ggplot(data = simdat[year > firstyear & cohort %% 5 == 0 & age <= 100], 
              aes(y = rdiff_s_pop, x = age, color = factor(cohort))) +
    geom_line(size = 1) + 
    facet_wrap(~fsex(sex)) +
    labs(title = paste0("Relative difference between healthy population ", 
                        "calculated from flows or directly simulated"))
  print(p)
  
  if (sum(simdat$s_cancer) > 0) {
    
    p <- ggplot(data = simdat[year > firstyear & cohort %% 5 == 0 & age <= 100], 
                aes(y = diff_s_cancer, x = age, color = factor(cohort))) +
      geom_line(size = 1) + 
      facet_wrap(~fsex(sex)) +
      labs(title = paste0("Difference between population with cancer ", 
                          "calculated from flows or directly simulated"))
    print(p)
    
    p <- ggplot(data = simdat[year > firstyear & cohort %% 5 == 0 & age <= 100], 
                aes(y = rdiff_s_cancer, x = age, color = factor(cohort))) +
      geom_line(size = 1) + 
      facet_wrap(~fsex(sex)) +
      labs(title = paste0("Relative difference between population with cancer ", 
                          "calculated from flows or directly simulated"))
    print(p)
  }
  
  if (sum(simdat$s_cvd) > 0) {
    
    p <- ggplot(data = simdat[year > firstyear & cohort %% 5 == 0 & age <= 100], 
                aes(y = diff_s_cvd, x = age, color = factor(cohort))) +
      geom_line(size = 1) + 
      facet_wrap(~fsex(sex)) +
      labs(title = paste0("Difference between population with CVD calculated ", 
                          "from flows or directly simulated"))
    print(p)
    
    p <- ggplot(data = simdat[year > firstyear & cohort %% 5 == 0 & age <= 100], 
                aes(y = rdiff_s_cvd, x = age, color = factor(cohort))) +
      geom_line(size = 1) + 
      facet_wrap(~fsex(sex)) +
      labs(title = paste0("Relative difference between population with CVD ", 
                          "calculated from flows or directly simulated"))
    print(p)
    
  }
   
  # Close pdf output
  dev.off(devpdf)
    
  
  ##****************************************************************************
  ##* If file for output from ode() exists validation a cohort based validation 
  ##* of stocks/flows using time step data from ode() is created.
  ##****************************************************************************
  
  filename <- paste0(projectroot, "/output/debug_ode_", timestamp, ".csv")
  if (file.exists(filename)) {

    pdffilename <- paste0(projectroot, "/output/debug_ode_", timestamp, ".pdf")
    
    pdf(file = pdffilename, paper = "a4r", width = 10, height = 6, 
        pointsize = 8)          
    devpdf <- dev.cur()

    ## Read csv file with debug data
    simdat <- fread(file = paste0(projectroot, "/output/debug_ode_", 
                                  timestamp, ".csv"))
    simdat[, cohort := year - age]

    # Determine time step within year
    minyear <- min(simdat$year)
    minsex <- min(simdat[year == minyear]$sex)
    minage <- min(simdat[year == minyear]$age)
    timestep <- length(unique(simdat[year == minyear & sex == minsex & 
                                       age == minage, 
                                     (steps = year - time)])) - 1
    
    # Select cohorts
    cohorts <- unique(simdat$cohort)
    cohorts_ <- sort(cohorts[cohorts %% 5 == 0])
    
    
    ##********************************************
    ## Stocks by cohort
    ##********************************************
    
    for (c in cohorts_) {
      pdat <- simdat[cohort == c, .(cohort, sex, time, s_pop, s_cvd, s_cancer)]
      pdat_ <- melt(pdat, id.vars = c("cohort", "sex", "time"))
      p <- ggplot(data = pdat_, aes(y = value, x = time, color = variable)) +
        geom_line(size = 1) +
        facet_wrap(~fsex(sex)) +
        labs(title = paste0("Simulated stocks for cohort: ", c),
             subtitle = paste0("Nr. of timesteps within year: ", timestep)) +
        ylim(0, NA) +
        theme(legend.position = "bottom")
      print(p)
    }

        
    ##********************************************
    ## Flows by cohort
    ##********************************************
    
    for (c in cohorts_) {
      pdat <- simdat[cohort == c, 
                     .(cohort, sex, time, f_dead, f_born, f_immig_pop,
                       f_immig_cvd, f_immig_cancer, f_emig_pop, f_emig_cvd, 
                       f_emig_cancer, f_pop_cancer, f_pop_cvd, f_cvd_pop, 
                       f_cancer_pop, f_cvd_dead, f_cancer_dead)]
      pdat_ <- melt(pdat, id.vars = c("cohort", "sex", "time"))
      p <- ggplot(data = pdat_, aes(y = value, x = time, color = variable)) +
        geom_line(size = 1) +
        facet_wrap(~fsex(sex)) +
        labs(title = paste0("Simulated flows for cohort: ", c),
             subtitle = paste0("Nr. of timesteps within year: ", timestep)) +
        ylim(0, NA) +
        theme(legend.position = "bottom")
      print(p)
    }
    
    dev.off(devpdf)
  } else {
    cat("No ode() output to validate!\n")
  }

}


# # Validation of demographics, stocks and flows
#validate_ncdsim(
#  # path to NCDSim project
#  projectroot = getwd(),
  # timestamp for simulation to validate
#  timestamp = "2025_05_27_14_22_1840")
