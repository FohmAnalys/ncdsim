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

validate_ncdsim <- function(
    projectroot = getwd(), # path to NCDSim project
    timestamp = NA,  # timestamp for simulation to validate
    use_scb_demographics = TRUE) { 
  
  ##****************************************************************************
  ## Reading simulation data
  ##****************************************************************************
  ##*

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
  ##*
  preval_new <- fread(paste0(projectroot, "/Input/stocks.csv"), sep = ';')
  incid_new <- fread(paste0(projectroot, "/Input/flows.csv"), sep = ';')

  
  preval_new <- dcast(preval_new, year+sex+age ~ stock, value.var = 'N')
  incid_new <- incid_new[!flow %like% "NA_"]
  incid_new <- dcast(incid_new, year+sex+age ~ flow, value.var = 'N' )
  heal_new <- incid_new[, c("year", "sex", "age", "cancer_healthy", "comorb_cancer", "comorb_cvd", "comorb_healthy","cvd_healthy")]
  colnames(preval_new) <- c("AR","Kon","alder","prev_cancer","prev_comorb", "prev_cvd" ,"totpop")
  incid_new[, c("cancer_cancer","cancer_healthy", "comorb_cancer", "comorb_comorb",
                "comorb_cvd", "comorb_healthy", "cvd_cvd", "cvd_healthy", "healthy_healthy") := NULL]
  
  colnames(incid_new) <- c("AR","Kon","alder","inc_cancer_cvd","inc_cancer_to_cvd","inc_cvd_to_cancer","inc_cvd_cancer","inc_cancer","inc_comorb","inc_cvd")
  
  incid <- CJ(AR = min(incid_new[, (AR)]):max(incid_new[,(AR)]),
              Kon = 1:2, alder = 0:100)
  heal <- CJ(year = min(heal_new[, (year)]):max(heal_new[,(year)]),
             sex = 1:2, age = 0:100)
  preval <- CJ(AR = min(preval_new[, (AR)]):max(preval_new[,(AR)]),
               Kon = 1:2, alder = 0:100)
  incid <- merge(incid, incid_new, by = c("AR", "alder", "Kon" ), all = TRUE)
  heal <- merge(heal, heal_new, by = c("year", "age", "sex" ), all = TRUE)
  preval <- merge(preval, preval_new[, c("AR", "alder", "Kon" , "prev_cvd", "prev_cancer", "prev_comorb")],
                  by = c("AR", "alder", "Kon"  ), all=TRUE)
  
  incid[, c('inc_cancer_to_cvd', 'inc_cvd_to_cancer') := NULL]
  colnames(incid) <- c("year","age","sex","inc_cancer_cvd","inc_cvd_cancer","inc_cancer","inc_comorb","inc_cvd")
  colnames(preval) <- c("year","age","sex","prev_cvd", "prev_cancer", "prev_comorb")
  colnames(heal) <- c("year", "age", "sex", "heal_cancer_healthy", "heal_comorb_cancer",  "heal_comorb_cvd", "heal_comorb_healthy", "heal_cvd_healthy")
  
  incid <- setnafill(incid, "const", fill = 0, cols = c("inc_cvd","inc_cancer","inc_cvd_cancer","inc_cancer_cvd","inc_comorb")) 
  heal <- setnafill(heal, "const", fill = 0, cols = c("heal_cancer_healthy", "heal_comorb_cancer",  "heal_comorb_cvd", "heal_comorb_healthy", "heal_cvd_healthy"))  
  
  
  preval <- setnafill(preval, "const", fill = 0, cols = c("prev_cvd","prev_cancer","prev_comorb"))
  
  
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
  
  t1 <- simdat[, .(pop = sum(s_pop + s_cancer + s_cvd + s_comorb), lbl = "NCDSim"), 
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
  
  t1 <- simdat[, .(pop = sum(s_pop + s_cancer + s_cvd + s_comorb), lbl = "NCDSim"), 
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
  
  t1 <- simdat[, .(ndead = sum(f_dead + f_cancer_dead + f_cvd_dead + f_comorb_dead),
                   nborn = sum(f_born),
                   nimmig = sum(f_immig_pop + f_immig_cvd + f_immig_cancer + f_immig_comorb),
                   nemig = sum(f_emig_pop + f_emig_cvd + f_emig_cancer + f_emig_comorb),
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
  
  t1 <- simdat[, .(ndead = sum(f_dead + f_cancer_dead + f_cvd_dead + f_comorb_dead),
                   nborn = sum(f_born),
                   nimmig = sum(f_immig_pop + f_immig_cvd + f_immig_cancer + f_immig_comorb),
                   nemig = sum(f_emig_pop + f_emig_cvd + f_emig_cancer + f_emig_comorb),
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
  ## Demographic flows by year
  ##*********************************************
  
  tmp <- melt(simdat[year >= (firstyear + 1), 
                     .(year, f_born, f_dead, f_cvd_dead, f_cancer_dead, f_comorb_dead,  
                       f_immig_pop, f_immig_cvd, f_immig_cancer, f_immig_comorb, f_emig_pop, 
                       f_emig_cvd, f_emig_cancer, f_emig_comorb)], 
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
  
  tmp <- simdat[, .(ncancer = sum(s_cancer), ncvd = sum(s_cvd), ncomorb = sum(s_comorb)), 
              by = year]
  pdat <- melt(tmp, id.vars = "year")
  tmp = preval[, .(ncancer_sos = sum(prev_cancer), ncvd_sos = sum(prev_cvd), ncomorb_sos = sum(prev_comorb)),
               by= year]
  tmp <- melt(tmp, id.vars = "year")
  
  pdat <- rbindlist(list(pdat, tmp), use.names = TRUE)
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0, NA) +
    ggtitle("Stocks of cancer, CVD and comorbidity (and validation data from SoS), by year")
  print(p)  
  
  
  
  ##*********************************************
  ## Total NCD stocks by year
  ##*********************************************
  
  tmp <- simdat[, .(nstock = sum(s_cancer) + sum(s_cvd) + sum(s_comorb)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  tmp = preval[, .(nstock_sos = sum(prev_cancer)+ sum(prev_cvd)+ sum(prev_comorb)),
               by= year]
  sosd <- melt(tmp, id.vars = "year")
  
  pdat <- rbindlist(list(pdat, sosd), use.names = TRUE)
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0,NA) + 
    ggtitle("Total stock of diseases (and validation data from SoS), by year")
  print(p)  
  
  
  ##*********************************************
  ## NCD stocks by year, sex and age group
  ##*********************************************
  
  tmp <- simdat[, .(ncancer = sum(s_cancer), ncvd = sum(s_cvd), ncomorb = sum(s_comorb)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  tmp = preval[, .(ncancer_sos = sum(prev_cancer), ncvd_sos = sum(prev_cvd), ncomorb_sos = sum(prev_comorb)),
               by= .(year, sex, agegrp = fagegrp10(age))]
  tmp <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
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
  ## NCD flows by year
  ##*********************************************
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(f_pop_cancer = sum(f_pop_cancer),
                  f_cancer_pop = sum(f_cancer_pop),
                  f_cancer_dead = sum(f_cancer_dead),
                  f_pop_cvd = sum(f_pop_cvd),
                  f_cvd_pop = sum(f_cvd_pop),
                  f_cvd_dead = sum(f_cvd_dead),
                  f_pop_comorb = sum(f_pop_comorb),
                  f_comorb_pop = sum(f_comorb_pop),
                  f_cancer_comorb = sum(f_cancer_comorb),
                  f_comorb_cancer = sum(f_comorb_cancer),
                  f_cvd_comorb = sum(f_cvd_comorb),
                  f_comorb_cvd = sum(f_comorb_cvd),
                  f_comorb_dead = sum(f_comorb_dead)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  pdat$source <- "NCDSim"

  sos_incid <- incid[,
                     .(f_pop_cancer = sum(inc_cancer),
                       f_pop_cvd = sum(inc_cvd),
                       f_cvd_comorb = sum(inc_cvd_cancer),
                       f_cancer_comorb = sum(inc_cancer_cvd),
                       f_pop_comorb = sum(inc_comorb)),
                     by = year]
  sos_incid = melt(sos_incid, id.vars = "year")
  sos_incid$source <- "SoS"
  heal_cvd <- heal[,
                   .(f_cancer_pop = sum(heal_cancer_healthy),
                   f_cvd_pop = sum(heal_cvd_healthy),
                   f_comorb_cvd = sum(heal_comorb_cvd),
                   f_comorb_cancer = sum(heal_comorb_cancer),
                   f_comorb_pop = sum(heal_comorb_healthy)),
                  by = year]
  
  heal_cvd <- melt(heal_cvd, id.vars = "year")
  heal_cvd$source <- "SoS"
  pdat <- rbindlist(list(pdat, heal_cvd, sos_incid), use.names = TRUE)

  p <- ggplot(data = pdat[variable %like% "cancer" & !variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    ylim(0, NA) +
    ggtitle("Flows related to cancer (and validation data from SoS), by year")
  print(p)  
  
  p <- ggplot(data = pdat[variable %like% "cvd" & !variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    ylim(0, NA) +
    ggtitle("Flows related to cvd (and validation data from SoS), by year")
  print(p)  
  
  p <- ggplot(data = pdat[variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    ylim(0, NA) +
    ggtitle("Flows related to comorbidity (and validation data from SoS), by year")
  print(p) 
  
  
  ##*********************************************
  ## NCD flows by year, sex and age group
  ##*********************************************
  ##*
  tmp <- simdat[year >= (firstyear + 1), 
                .(f_pop_cancer = sum(f_pop_cancer),
                  f_cancer_pop = sum(f_cancer_pop),
                  f_cancer_dead = sum(f_cancer_dead),
                  f_pop_cvd = sum(f_pop_cvd),
                  f_cvd_pop = sum(f_cvd_pop),
                  f_cvd_dead = sum(f_cvd_dead),
                  f_pop_comorb = sum(f_pop_comorb),
                  f_comorb_pop = sum(f_comorb_pop),
                  f_cancer_comorb = sum(f_cancer_comorb),
                  f_comorb_cancer = sum(f_comorb_cancer),
                  f_cvd_comorb = sum(f_cvd_comorb),
                  f_comorb_cvd = sum(f_comorb_cvd),
                  f_comorb_dead = sum(f_comorb_dead)), 
                by = .(year, sex, agegrp = fagegrp10(age))]
  pdat <- melt(tmp, id.vars = c("year", "sex", "agegrp"))
  pdat$source <- "NCDSim"

  sos_incid <- incid[,
                     .(f_pop_cancer = sum(inc_cancer),
                       f_pop_cvd = sum(inc_cvd),
                       f_cvd_comorb = sum(inc_cvd_cancer),
                       f_cancer_comorb = sum(inc_cancer_cvd),
                       f_pop_comorb = sum(inc_comorb)),
                     by = .(year, sex, agegrp = fagegrp10(age))]
  sos_incid = melt(sos_incid, c("year", "sex", "agegrp"))
  sos_incid$source <- "SoS"
  heal_cvd <- heal[,
                   .(f_cancer_pop = sum(heal_cancer_healthy),
                     f_cvd_pop = sum(heal_cvd_healthy),
                     f_comorb_cvd = sum(heal_comorb_cvd),
                     f_comorb_cancer = sum(heal_comorb_cancer),
                     f_comorb_pop = sum(heal_comorb_healthy)),
                   by = .(year, sex, agegrp = fagegrp10(age))]
  
  heal_cvd <- melt(heal_cvd, id.vars = c("year", "sex", "agegrp"))
  heal_cvd$source <- "SoS"
  pdat <- rbindlist(list(pdat, heal_cvd, sos_incid), use.names = TRUE)
  #pdat[ !variable %like% "sos", source := 'NCDSim']
  #pdat[ variable %like% "sos", source := 'SoS']
  
  p <- ggplot(data = pdat[sex == 1 & variable %like% "cancer" & !variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to cancer (and validation data from SoS), ",
                   "by year, men"))
  print(p)  
  p <- ggplot(data = pdat[sex == 2 & variable %like% "cancer" & !variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to cancer (and validation data from SoS), ",
                   "by year, women"))
  print(p)  
  
  p <- ggplot(data = pdat[sex == 1 & variable %like% "cvd" & !variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to cvd (and validation data from SoS), ",
                   "by year, men"))
  print(p) 
  
  p <- ggplot(data = pdat[sex == 2 & variable %like% "cvd" & !variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to cvd (and validation data from SoS), ",
                   "by year, women"))
  print(p) 
  
  p <- ggplot(data = pdat[sex == 1 & variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to comorbidity (and validation data from SoS), ",
                   "by year, men"))
  print(p) 
  
  p <- ggplot(data = pdat[sex == 2 & variable %like% "comorb"], 
              aes(y = value, x = year, color = variable, linetype = source)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = c("NCDSim" = "solid", "SoS" = "dashed")) +
    facet_wrap(~agegrp, scales = "free") +
    ylim(0, NA) +
    ggtitle(paste0("Flows related to comorbidity (and validation data from SoS), ",
                   "by year, women"))
  print(p) 
  
  
  
  ##*********************************************
  ## Death rates by year, sex, age and group
  ##*********************************************
  
  tmp <- simdat[year %% 2 == 0, .(year, sex, age, dr, dr_cancer, dr_cvd, dr_comorb)]
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
  
  tmp <- simdat[year %% 2 == 0, .(year, sex, age, dr, dr_cancer, dr_cvd, dr_comorb)]
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
    ggtitle("Direct cost of cancer and CVD in SEK, by year")
  print(p)  
  
  tmp <- simdat[year >= (firstyear + 1), 
                .(cost_cancer = sum(icost_cancer), 
                  cost_cvd = sum(icost_cvd)), 
                by = year]
  pdat <- melt(tmp, id.vars = "year")
  p <- ggplot(data = pdat, aes(y = value, x = year, color = variable)) +
    geom_line(size = 1) +
    ylim(0, NA) +
    ggtitle("Indirect cost of cancer and CVD in SEK, by year")
  print(p)  
  
  ##*********************************************
  ## calibration, raw and normalized för IR
  ##*********************************************  

  tmp = simdat[year >= (firstyear + 1) & age > 40, .(calibration_cancer = mean(calibration_cancer / (f_pop_cancer/s_pop)),
                                calibration_cvd = mean(calibration_cvd / (f_pop_cvd/s_pop)),
                                calibration_comorb = mean(calibration_comorb / (f_pop_comorb/s_pop)),
                                calibration_cancer_comorb = mean(calibration_cancer_comorb /(f_cancer_comorb/s_cancer) ),
                                calibration_cvd_comorb = mean(calibration_cvd_comorb / (f_cvd_comorb/s_cvd))),
                                by = .(sex, age)]
  
  tmp2 = simdat[year >= (firstyear + 1) & age > 40, .(calibration_cancer = mean(calibration_cancer),
                                                     calibration_cvd = mean(calibration_cvd),
                                                     calibration_comorb = mean(calibration_comorb),
                                                     calibration_cancer_comorb = mean(calibration_cancer_comorb),
                                                     calibration_cvd_comorb = mean(calibration_cvd_comorb)),
                                                     by = .(sex, age)]
  
  pdat = melt(tmp2, id.vars=c("age", "sex"))
  p <- ggplot(data = pdat, 
              aes(y = value, x = age, color = variable)) +
    geom_line(size = 1) + 
    facet_wrap(~fsex(sex)) +
    ggtitle("Calibration constants, per sex")
  print(p)
  
  
  pdat = melt(tmp, id.vars=c("age", "sex"))
  p <- ggplot(data = pdat, 
              aes(y = value, x = age, color = variable)) +
    geom_line(size = 1) + 
    facet_wrap(~fsex(sex)) +
    ggtitle("Constants normalized by IR, per sex")
  print(p)
  
  ###***************************************************************************
  ###* mortality calibration
  t <- simdat[year >= (firstyear + 1), .(mortcal = mean(calibration_mortality)),
              by=c("year", "sex")]
  p = ggplot(data = t, aes(x = year, y= mortcal, color = as.factor(sex))) +
    geom_line(size = 1) + 
    ggtitle("Mortality adjustment per year and sex")
  print(p)
  
  
  for (y in (firstyear+1):lastyear){
    pdat <- simdat[year == y]
    p <- ggplot(data = pdat, aes(y = calibration_mortality, x = age)) +
      geom_line(size = 1) +
      facet_wrap(~fsex(sex)) +
      ggtitle(paste0("Mortality adjustment in year", y)) +
      ylim(0, NA) +
      theme(legend.position = "bottom")
    print(p)
  }
  
  ##********************************************
  ## PAF by year, sex and age 
  ##********************************************
  
  tmp <- simdat[, .(year, sex, age, paf_cvd_smoking, paf_cvd_inactivity,
                    paf_cvd_bmi, paf_cvd_alcohol, 
                    paf_cancer_smoking, paf_cancer_inactivity, 
                    paf_cancer_bmi, paf_cancer_alcohol, cpaf_diet_cancer,
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
      f_cancer_pop + f_cvd_pop - f_pop_cancer - f_pop_cvd + f_comorb_pop - f_pop_comorb,
    s_cancer_ = shift(s_cancer) + f_pop_cancer - f_cancer_pop + f_immig_cancer - 
      f_emig_cancer - f_cancer_dead + f_comorb_cancer - f_cancer_comorb,
    s_cvd_ = shift(s_cvd) + f_pop_cvd - f_cvd_pop + f_immig_cvd - f_emig_cvd - 
      f_cvd_dead + f_comorb_cvd - f_cvd_comorb,
    s_comorb_ = shift(s_comorb) +  f_pop_comorb - f_comorb_pop + f_immig_comorb - f_emig_comorb-
      f_comorb_dead + f_cancer_comorb - f_comorb_cancer + f_cvd_comorb - f_comorb_cvd 
  )]
  simdat[year == firstyear | cohort != shift(cohort) | sex != shift(sex) |
           age != (shift(age) + 1) | age == 100, ":="(
             # NA for non-valid lags.
             s_pop_ = NA,
             s_cancer_ = NA,
             s_cvd_ = NA,
             s_comorb_ = NA
           )]
  
  simdat[, ":="(
    diff_s_pop = s_pop_ - s_pop,
    rdiff_s_pop = s_pop_ / s_pop,
    diff_s_cancer = s_cancer_ - s_cancer,
    rdiff_s_cancer = s_cancer_ / s_cancer,
    diff_s_cvd = s_cvd_ - s_cvd,
    rdiff_s_cvd = s_cvd_ / s_cvd,
    diff_s_comorb = s_comorb_ - s_comorb,
    rdiff_s_comorb = s_comorb_ / s_comorb
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
  
  if (sum(simdat$s_comorb) > 0) {
    
    p <- ggplot(data = simdat[year > firstyear & cohort %% 5 == 0 & age <= 100], 
                aes(y = diff_s_comorb, x = age, color = factor(cohort))) +
      geom_line(size = 1) + 
      facet_wrap(~fsex(sex)) +
      labs(title = paste0("Difference between population with comorbidity calculated ", 
                          "from flows or directly simulated"))
    print(p)
    
    p <- ggplot(data = simdat[year > firstyear & cohort %% 5 == 0 & age <= 100], 
                aes(y = rdiff_s_comorb, x = age, color = factor(cohort))) +
      geom_line(size = 1) + 
      facet_wrap(~fsex(sex)) +
      labs(title = paste0("Relative difference between population with comorbidity ", 
                          "calculated from flows or directly simulated"))
    print(p)
    
  }
  
  # Close pdf output
  dev.off(devpdf)
  
}




