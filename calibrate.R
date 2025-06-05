##### Calibration of NCDSim
# calibrate incidence rate, direct costs and indirect costs to match the observed.
# Returns:
# -- constants: calibration constants (per age and sex):
#       constant factors so that the cpaf from ALL the lifestyle factors equals
#       the total observed cpaf from IHE reports.
# -- paf_other_cancer: paf from other riskfactors for cancer (per year, age and sex):
#       the contribution from other riskfactors (such as genetics and chance)
#       to match the total incidence rate, for cancer.
# -- paf_other_cvd: paf from other riskfactors for cvd (per year, age and sex):
#       the contribution from other riskfactors (such as genetics and chance)
#       to match the total incidence rate, for cvd.
# -- dcost_unit: unit (i.e. per patient) direct cost (per age and sex):
#       calibrated unit cost so that the total cost matches the total direct cost
#       for the given base year
# -- icost_unit: unit (i.e. per patient) indirect cost (per age and sex):
#       calibrated unit cost so that the total cost matches the total indirect cost
#       for the given base year



calibrate <- function(input_path, prop_cancer, prop_cvd, age_cutoff_cancer, 
                      age_cutoff_cvd,  rr, rr_cancer_diet,
                      rr_cvd_diet, communalities, direct_costs, indirect_costs,
                      dcost_total_base_year){
  
  # Read input files
  incid_new <- fread(paste0(input_path, "incidenser.csv"), sep = ";")
  preval_new <- fread(paste0(input_path, "prevalenser.csv"), sep = ";")
  population <- fread(paste0(input_path, "pop_counts_scb.csv"), sep = ";")
  colnames(incid_new) <- c("year", "age", "sex", "cvd", "cancer" )
  colnames(preval_new) <- c("year", "age", "sex", "cvd", "cancer", "totpop")
  incid <- CJ(year = min(incid_new[, (year)]):max(incid_new[,(year)]), 
              sex = 1:2, age = 0:100)
  preval <- CJ(year = min(preval_new[, (year)]):max(preval_new[,(year)]),
               sex = 1:2, age = 0:100)
  incid <- merge(incid, incid_new, by = c("year", "age", "sex" ), all = TRUE)
  preval <- merge(preval, preval_new[, c("year", "age", "sex", "cvd", 
                                         "cancer")],
                  by = c("year","age", "sex" ), all=TRUE)
  incid <- melt(incid, id.vars = c("sex", "age","year"),
                value.vars = c("cvd", "cancer"))
  preval <- melt(preval, id.vars = c("sex", "age","year"),
                 value.vars = c("cvd", "cancer"))
  colnames(incid) <- c("sex", "age", "year", "ncd", "incidence")
  colnames(preval) <- c("sex", "age", "year", "ncd", "prevalence")
  incid <- setnafill(incid, "const", fill = 0, cols = c("incidence")) 
  preval <- setnafill(preval, "const", fill = 0, cols = c("prevalence")) 
  
  costs <- fread(paste0(input_path, "costs_direct.csv"), sep = ";")
  colnames(costs) <- c("age", 1, 2, "ncd")
  
  ind_costs <- fread(paste0(input_path,"costs_indirect.csv"), sep = ";")
  colnames(ind_costs) <- c("age", 1, 2, "ncd")
  
  # Compute incidence rate (dvs incidence/healthy_pop) from observed values,
  # only for the years where we have both incidences and prevalences.
  list_of_years <- intersect(unique(incid[,year]), unique(preval[, year]))
  incid_rate <- CJ(sex = 1:2, age = 0:100, ncd = c("cancer", "cvd"))
  for (sex_ in c(1, 2)){
    for (ncd_ in c("cvd", "cancer")){
      const_df <- data.table()
      for (year_ in lapply(list_of_years,sort, decreasing = TRUE)){
        pop <- population[year == year_ & sex == sex_, pop]
        prev_cancer <- preval[year == year_ & sex == sex_ & ncd == "cancer",
                              prevalence]
        prev_cvd <- preval[year == year_ & sex == sex_ & ncd == "cvd",
                           prevalence]
        incidence <- incid[year == year_ & sex == sex_ & ncd == ncd_, incidence]
        risk <- incidence / (pop - prev_cvd - prev_cancer)
        const_df[, (as.character(year_)) := risk]
      }
      incid_rate[ncd == ncd_ & sex == sex_, rr := apply(const_df, 1, mean)]
    }
  }

  # read prevalences for risk factors and relative risks for food
  filter_cols <- c(c("sex", "age"), as.character(list_of_years))
  food_prevalences <- fread(paste0(input_path, "prev_diet.csv"), sep = ";")
  prev_smoking <- fread(paste0(input_path, "prev_smoking.csv"), sep = ";")
  prev_obesity <- fread(paste0(input_path, "prev_obesity.csv"), sep = ";")
  prev_alcohol <- fread(paste0(input_path, "prev_alcohol.csv"), sep = ";")
  prev_inactivity <- fread(paste0(input_path, "prev_inactivity.csv"), sep = ";")
  prev_smoking <- prev_smoking[, ..filter_cols]
  prev_obesity <- prev_obesity[, ..filter_cols]
  prev_alcohol <- prev_alcohol[, ..filter_cols]
  prev_inactivity <- prev_inactivity[, ..filter_cols]
  
  food_categories <- colnames(food_prevalences)[3:length(food_prevalences)]
  rr_cancer_diet <- data.frame(as.list(rr_cancer_diet))
  rr_cvd_diet <- data.frame(as.list(rr_cvd_diet))

  
  # calculate cPAF for food
  cpaf_diet <- CJ(sex = 1:2, age= 0:100, ncd = c("cancer","cvd"))
  for (sex_ in 1:2) {
    prev_diet <- food_prevalences[sex == sex_, ..food_categories]
    pafs_diet_cancer <- rbindlist(apply(prev_diet, 1, paf, rr = rr_cancer_diet))
    cpaf_diet_cancer <- (1 - apply(1-pafs_diet_cancer, 1, prod)) 
    pafs_diet_cvd <- rbindlist(apply(prev_diet, 1, paf, rr = rr_cvd_diet))
    cpaf_diet_cvd <-( 1 - apply(1-pafs_diet_cvd, 1, prod))
    cpaf_diet[ncd == "cancer" & sex == sex_, cpaf := cpaf_diet_cancer]
    cpaf_diet[ncd=="cvd" & sex == sex_,cpaf := cpaf_diet_cvd]
  }
  
  # calculate cPAF per year for lifestyle factors
  paf_cancer <- CJ(sex = 1:2, age= 0:100, year = as.character(list_of_years))
  paf_cvd <- CJ(sex = 1:2, age= 0:100, year = as.character(list_of_years))
  for (sex_ in 1:2) {
    for (year_ in as.character(list_of_years)) {
      paf_diet_c_ <- cpaf_diet[sex == sex_ & ncd == "cancer", cpaf]
      paf_alcohol_c <- unname((communalities["alcohol"])) *
        paf(prev_alcohol[sex == sex_, ..year_], unname(rr["rr_cancer_alcohol"]))
      paf_smoking_c <- unname((communalities["smoking"])) *
        paf(prev_smoking[sex == sex_, ..year_], unname(rr["rr_cancer_smoking"]))
      paf_obesity_c <- unname((communalities["smoking"])) *
        paf(prev_obesity[sex == sex_, ..year_], unname(rr["rr_cancer_obesity"]))
      paf_inactivity_c <- unname((communalities["inactivity"])) *
        paf(prev_inactivity[sex == sex_, ..year_],
            unname(rr["rr_cancer_inactivity"]))
      cpaf_c <- 1 - (1 - paf_diet_c_) *
        (1 - paf_alcohol_c) *
        (1 - paf_smoking_c) *
        (1 - paf_obesity_c) *
        (1 - paf_inactivity_c)
      paf_cancer[sex == sex_ & year == year_, cpaf := cpaf_c]
      
      paf_diet_cvd_ <- cpaf_diet[sex == sex_ & ncd == "cvd", cpaf]
      paf_alcohol_cvd <- unname((communalities["alcohol"])) *
        paf(prev_alcohol[sex == sex_, ..year_], unname(rr["rr_cvd_alcohol"]))
      paf_smoking_cvd <- unname((communalities["smoking"])) *
        paf(prev_smoking[sex == sex_, ..year_], unname(rr["rr_cvd_smoking"]))
      paf_obesity_cvd <- unname((communalities["smoking"])) *
        paf(prev_obesity[sex == sex_, ..year_], unname(rr["rr_cvd_obesity"]))
      paf_inactivity_cvd <- unname((communalities["inactivity"])) *
        paf(prev_inactivity[sex == sex_, ..year_],
            unname(rr["rr_cvd_inactivity"]))
      cpaf_cvd <- 1 - (1 - paf_diet_cvd_) *
        (1 - paf_alcohol_cvd) *
        (1 - paf_smoking_cvd) *
        (1 - paf_obesity_cvd) *
        (1 - paf_inactivity_cvd)
      paf_cvd[sex == sex_ & year == year_, cpaf := cpaf_cvd]
      
    } 
  }
  
  # calculate the total yearly mean cPAF (dvs food + non food) 
  cpaf_cvd_mean <- paf_cvd[, .(cpaf = mean(cpaf)), 
                                 by = .(age, sex)]
  cpaf_cancer_mean <- paf_cancer[, .(cpaf = mean(cpaf)), 
                                 by = .(age, sex)]
  
  # compute the calibration constants
  # If the age >= age_cutoff then calibrate to the observed cpaf
  # If the age < age_cutof, then the calibration is 0 and all cases come from cpaf_other
  constants_df <- CJ(age = 0:100, sex = 1:2, ncd=c("cancer","cvd"))
  for (sex_ in 1:2){
    pf_c <- cpaf_cancer_mean[sex == sex_ & age >= age_cutoff_cancer]
    ir_c <- incid_rate[sex == sex_ & ncd == "cancer" & age >= age_cutoff_cancer]
    c_c <- (prop_cancer * ir_c$rr) / pf_c$cpaf
    constants_df[sex == sex_ & ncd == "cancer" &
                   age >= age_cutoff_cancer, const := c_c]
    
    pf_cvd <- cpaf_cvd_mean[sex == sex_ & age >= age_cutoff_cvd]
    ir_cvd <- incid_rate[sex == sex_ & ncd == "cvd" & age >= age_cutoff_cvd]
    c_cvd <- (prop_cvd * ir_cvd$rr) / pf_cvd$cpaf
    constants_df[sex == sex_ & ncd == "cvd" &
                   age >= age_cutoff_cvd, const := c_cvd]    
    
  }
  setnafill(constants_df, "const", fill=0.0, cols=c("const"))
  
  # Calculate paf_other:
  # if age >= age_cutoff then paf_other is the complement to paf from lifestyle
  # if age < age_cutoff then paf_other = incidence rate
  paf_other_cancer <- CJ(sex = 1:2, age= 0:100,
                         year = as.character(list_of_years))
  paf_other_cvd <- CJ(sex = 1:2, age= 0:100,
                      year = as.character(list_of_years))
  for (sex_ in 1:2){
    for (y in list_of_years){
      c_c <- constants_df[sex == sex_ & ncd == "cancer" &
                            age >= age_cutoff_cancer, const]
      ir_c <- incid_rate[sex == sex_ & ncd == "cancer" &
                           age >= age_cutoff_cancer, rr]
      cut_c <- paf_cancer[year == y & sex == sex_,]
      
      paf_slice_c <- cut_c[age >= age_cutoff_cancer, cpaf]
      lb_c <- (c_c*paf_slice_c - ir_c) / (c_c*paf_slice_c - 1)
      res_c <- incid_rate[sex == sex_ & ncd == "cancer" &
                            age < age_cutoff_cancer, rr]
      paf_other_cancer[year == y & sex == sex_ &
                      age >= age_cutoff_cancer, cpaf := lb_c]
      paf_other_cancer[year == y & sex == sex_ &
                      age < age_cutoff_cancer, cpaf := res_c]
      
      c_cvd <- constants_df[sex == sex_ & ncd == "cvd" &
                              age >= age_cutoff_cvd, const]
      ir_cvd <- incid_rate[sex == sex_ & ncd == "cvd" &
                             age >= age_cutoff_cvd, rr]
      cut_cvd <- paf_cvd[year == y & sex == sex_,]
      paf_slice_cvd <- cut_cvd[age >= age_cutoff_cvd, cpaf]
      lb_cvd <- (c_cvd*paf_slice_cvd - ir_cvd) / (c_cvd*paf_slice_cvd - 1)
      res_cvd <- incid_rate[sex == sex_ & ncd == "cvd" &
                              age < age_cutoff_cvd, rr]
      paf_other_cvd[year == y & sex == sex_ &
                      age >= age_cutoff_cvd, cpaf := lb_cvd]
      paf_other_cvd[year == y & sex == sex_ &
                      age < age_cutoff_cvd, cpaf := res_cvd]
    }
  }
  
 
  # Calibrate costs

  # Direct costs, adjust KPP with a and b so that (a + b*KPP)*stock = total costs
  out_kpp <- data.table()
  out_ik <- data.table()
  
  # function to optimize
  f <- function (b) {
    ret <- sum((kpp_1 + b) * prev) - del_cost
    return(ret)
  }
  f1 <- Vectorize(f)
  
  for( ncd_ in c("cancer","cvd")) {
    tot_sick <- sum(preval[year == dcost_total_base_year & ncd == ncd_, 
                           prevalence])
    tot_sick
    d_costs <- costs[ncd == ncd_]
    d_costs[,ncd := NULL]
    d_costs <- melt(d_costs, id.vars="age")
    
    # Adjust indirect costs i_k with a constant c
    # so that c*i_k*stock = total indirect costs  
    
    i_costs <- ind_costs[ncd == ncd_]
    i_costs[,ncd := NULL]
    i_costs <- melt(i_costs, id.vars="age")
    
    out_kpp_0 <- CJ(age = 0:100, sex = 1:2, ncd = ncd_)
    # optimization step per sex
    for (s in 1:2 ){
      c_m <- d_costs[variable == s, value]
      prev <- preval[year == dcost_total_base_year & sex == s & ncd == ncd_, 
                     prevalence]
      tot_p <- sum(prev)
      del_cost <- direct_costs[[ncd_]] * tot_p / tot_sick
      kpp_1 <- del_cost * c_m / (tot_p * sum(c_m))
      rr <- uniroot(f1, c(-max(c_m)*10, max(c_m) * 10))
      kpp_1 <- kpp_1 + rr$root
      out_kpp_0[sex == s, dcost_unit := kpp_1]
    }
    
    # indirect costs. Optimize so the unit_costs*NCDSim_stock = total_costs
    out_ik_0 <- CJ(age=0:100, sex=1:2, ncd = ncd_)
    for (s in 1:2 ){
      c_m <- i_costs[variable == s, value]
      prev <- preval[year == dcost_total_base_year & sex == s & ncd == ncd_,
                     prevalence]
      tot_p <- sum(prev)
      del_cost <- (indirect_costs[[ncd_]]/tot_sick)  * tot_p 
      ik_i <- c_m * del_cost / prev
      ik_i[is.na(ik_i)] <- 0
      out_ik_0[sex == s, icost_unit := ik_i]
    }
    
    out_kpp <- rbind(out_kpp, out_kpp_0)
    out_ik <- rbind(out_ik, out_ik_0)
  }
  
  
  return(list(constants = constants_df, paf_other_cancer = paf_other_cancer,
              paf_other_cvd = paf_other_cvd, dcost_unit = out_kpp, icost_unit = out_ik))
}

