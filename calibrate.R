library(data.table)
library(ggplot2)
library(gridExtra)
library(scales)
library(stringr)
library(jsonlite)


calibrate <- function(input_path, calibrate_cpaf_cancer, calibrate_cpaf_cvd,
                      prop_cancer, prop_cvd, age_cutoff_cancer, 
                      age_cutoff_cvd,  age_cutoff_cancer_high, 
                      age_cutoff_cvd_high, rr, rr_cancer_diet,
                      rr_cvd_diet, communalities, direct_costs, indirect_costs,
                      dcost_total_base_year){

# read base input for costs, population, incidences and prevalences
# extract and format data tables for incidences, prevalences, healing rates, KPP and population
costs <- fread(paste0(input_path, "costs_direct.csv"), sep = ";")
colnames(costs) <- c("age", 1, 2, "ncd")

ind_costs <- fread(paste0(input_path,"costs_indirect.csv"), sep = ";")
colnames(ind_costs) <- c("age", 1, 2, "ncd")

population <- fread(paste0(input_path, "pop_counts_scb.csv"), sep = ";")

preval_new <- fread(paste0(input_path, "stocks.csv"), sep = ';')
incid_new <- fread(paste0(input_path, "flows.csv"), sep = ';')

preval_new <- dcast(preval_new, year+sex+age ~ stock, value.var = 'N')
incid_new <- incid_new[!flow %like% "NA_"]
incid_new <- dcast(incid_new, year+sex+age ~ flow, value.var = 'N' )
colnames(preval_new) <- c("AR","Kon","alder","prev_cancer","prev_comorb",
                          "prev_cvd" ,"totpop")
heal_new <- incid_new[, c("year", "sex", "age", "cancer_healthy",
                          "comorb_cancer", "comorb_cvd", "comorb_healthy",
                          "cvd_healthy")]
incid_new[, c("cancer_cancer","cancer_healthy", "comorb_cancer", "comorb_comorb",
                "comorb_cvd", "comorb_healthy", "cvd_cvd",
              "cvd_healthy", "healthy_healthy") := NULL]

colnames(incid_new) <- c("AR","Kon","alder","inc_cancer_cvd","inc_cancer_to_cvd",
                         "inc_cvd_to_cancer","inc_cvd_cancer","inc_cancer",
                         "inc_comorb","inc_cvd")

incid <- CJ(AR = min(incid_new[, (AR)]):max(incid_new[,(AR)]),
            Kon = 1:2, alder = 0:100)
heal <- CJ(year = min(heal_new[, (year)]):max(heal_new[,(year)]),
           sex = 1:2, age = 0:100)
preval <- CJ(AR = min(preval_new[, (AR)]):max(preval_new[,(AR)]),
             Kon = 1:2, alder = 0:100)
incid <- merge(incid, incid_new, by = c("AR", "alder", "Kon" ), all = TRUE)
heal <- merge(heal, heal_new, by = c("year", "age", "sex" ), all = TRUE)
preval <- merge(preval, preval_new[, c("AR", "alder", "Kon" , "prev_cvd",
                                       "prev_cancer", "prev_comorb")],
                by = c("AR", "alder", "Kon"  ), all=TRUE)

incid[, c('inc_cancer_to_cvd', 'inc_cvd_to_cancer') := NULL]
colnames(incid) <- c("year","age","sex","inc_cancer_cvd","inc_cvd_cancer",
                     "inc_cancer","inc_comorb","inc_cvd")
colnames(preval) <- c("year","age","sex","prev_cvd", "prev_cancer", "prev_comorb")
colnames(heal) <- c("year", "age", "sex", "heal_cancer_healthy",
                    "heal_comorb_cancer",  "heal_comorb_cvd",
                    "heal_comorb_healthy", "heal_cvd_healthy")

incid <- setnafill(incid, "const", fill = 0, cols = c("inc_cvd","inc_cancer",
                                                      "inc_cvd_cancer",
                                                      "inc_cancer_cvd",
                                                      "inc_comorb")) 
preval <- setnafill(preval, "const", fill = 0, cols = c("prev_cvd",
                                                        "prev_cancer",
                                                        "prev_comorb")) 
heal <- setnafill(heal, "const", fill = 0, cols = c("heal_cancer_healthy",
                                                    "heal_comorb_cancer",
                                                    "heal_comorb_cvd", 
                                                    "heal_comorb_healthy",
                                                    "heal_cvd_healthy")) 
list_of_years <- intersect(unique(incid[,year]), unique(preval[, year]))
population <- population[year>=list_of_years[1] &
                           year<= list_of_years[length(list_of_years)]]

preval <- merge(population[,.(year,age,sex,pop)], preval,
                by = c("year","age", "sex"), all=TRUE)
preval[, frisk := pop - prev_cancer - prev_cvd - prev_comorb]
cols_to_convert = c("prev_cvd", "pop", "prev_cancer", "prev_comorb", "frisk")
preval[, (cols_to_convert) := lapply(.SD, as.double), .SDcols = cols_to_convert]


# calculate riskgroup by taking the mean of the stocks between year x and year x +1
# To account for the yearly resolution (31 december) for the input files

xpreval <- melt(preval, id.vars = c("year", "sex", "age"))
xpreval[, cohort := year - age]
setkey(xpreval, variable, sex, cohort, year)
xpreval[, value1 := shift(value)]
xpreval[!(cohort == shift(cohort) & sex == shift(sex) &
            year == (shift(year) + 1)), value1 := NA]
xpreval[, meanvalue := (value + value1) / 2]
riskgroup <- dcast(xpreval[, .(year, sex, age, variable, rg = meanvalue)],
                   year + sex + age ~ variable, value.var = "rg")


riskgroup <- setnafill(riskgroup, "const", 0, cols=c("prev_cvd","prev_cancer", 
                                                     "prev_comorb", "frisk"))

# calculate incidence rates as incidence/risk group for the different flows
incid_rate <- incid[riskgroup, on = c("year","age", "sex"), nomatch=NULL]
heal <- heal[riskgroup, on = c("year","age", "sex"), nomatch=NULL]

incid_rate[,ir_cvd := inc_cvd / frisk]
incid_rate[,ir_cancer := inc_cancer / frisk]
incid_rate[, ir_cvd_cancer := inc_cvd_cancer / prev_cvd]
incid_rate[, ir_cancer_cvd := inc_cancer_cvd / prev_cancer]
incid_rate[, ir_comorb := inc_comorb / frisk]
incid_rate[,c("pop","inc_cvd","inc_cancer","inc_cvd_cancer","inc_cancer_cvd",
              "inc_comorb",  "prev_cvd", "prev_cancer",
              "prev_comorb", "frisk"):=NULL]
incid_rate <- setnafill(incid_rate, "const", fill=0, cols=colnames(incid_rate))
incid_rate <- incid_rate[year > min(incid_rate$year)]


# calculate healing rates

heal[,hr_cvd := heal_cvd_healthy/prev_cvd]
heal[,hr_cancer := heal_cancer_healthy/prev_cancer]
heal[,hr_comorb_cancer := heal_comorb_cancer/prev_comorb]
heal[,hr_comorb_cvd := heal_comorb_cvd/prev_comorb]
heal[,heal_comorb_healthy := heal_comorb_healthy/prev_comorb]
heal_rate <- heal[, .(cvd = mean(hr_cvd), cancer = mean(hr_cancer),
                      comorb_cancer = mean(hr_comorb_cancer),
                      comorb_cvd = mean(hr_comorb_cvd),
                      comorb_healthy = mean(heal_comorb_healthy)), by=c('age', 'sex')]
heal_rate <- setnafill(heal_rate, "const", 0, cols=c("cvd","cancer",
                                                     "comorb_cancer",
                                                     "comorb_cvd",
                                                     "comorb_healthy"))


# read prevalences for risk factors and relative risks for food
filter_cols <- c(c("sex", "age"), as.character(list_of_years))
food_prevalences <- fread(paste0(input_path, "prev_diet.csv"), sep = ";")
prev_smoking <- fread(paste0(input_path, "prev_smoking.csv"), sep = ";")
prev_bmi <- fread(paste0(input_path, "prev_bmi.csv"), sep = ";")
prev_alcohol <- fread(paste0(input_path, "prev_alcohol.csv"), sep = ";")
prev_inactivity <- fread(paste0(input_path, "prev_inactivity.csv"), sep = ";")
prev_bmi <- prev_bmi[, ..filter_cols]
prev_alcohol <- prev_alcohol[, ..filter_cols]
prev_inactivity <- prev_inactivity[, ..filter_cols]
prev_smoking <- prev_smoking[, ..filter_cols]

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
    paf_bmi_c <- unname((communalities["smoking"])) *
      paf(prev_bmi[sex == sex_, ..year_], unname(rr["rr_cancer_bmi"]))
    paf_inactivity_c <- unname((communalities["inactivity"])) *
      paf(prev_inactivity[sex == sex_, ..year_],
          unname(rr["rr_cancer_inactivity"]))
    cpaf_c <- 1 - (1 - paf_diet_c_) *
      (1 - paf_alcohol_c) *
      (1 - paf_smoking_c) *
      (1 - paf_bmi_c) *
      (1 - paf_inactivity_c)
    paf_cancer[sex == sex_ & year == year_, cpaf := cpaf_c]
    
    paf_diet_cvd_ <- cpaf_diet[sex == sex_ & ncd == "cvd", cpaf]
    paf_alcohol_cvd <- unname((communalities["alcohol"])) *
      paf(prev_alcohol[sex == sex_, ..year_], unname(rr["rr_cvd_alcohol"]))
    paf_smoking_cvd <- unname((communalities["smoking"])) *
      paf(prev_smoking[sex == sex_, ..year_], unname(rr["rr_cvd_smoking"]))
    paf_bmi_cvd <- unname((communalities["smoking"])) *
      paf(prev_bmi[sex == sex_, ..year_], unname(rr["rr_cvd_bmi"]))
    paf_inactivity_cvd <- unname((communalities["inactivity"])) *
      paf(prev_inactivity[sex == sex_, ..year_],
          unname(rr["rr_cvd_inactivity"]))
    cpaf_cvd <- 1 - (1 - paf_diet_cvd_) *
      (1 - paf_alcohol_cvd) *
      (1 - paf_smoking_cvd) *
      (1 - paf_bmi_cvd) *
      (1 - paf_inactivity_cvd)
    paf_cvd[sex == sex_ & year == year_, cpaf := cpaf_cvd]
    
  } 
}
paf_comorb <- paf_cvd[, cpaf := pmax(paf_cancer$cpaf, paf_cvd$cpaf)]
# calculate the total yearly mean cPAF (dvs food + non food) 
cpaf_cvd_mean <- paf_cvd[, .(cpaf = mean(cpaf)), 
                         by = .(age, sex)]
cpaf_cancer_mean <- paf_cancer[, .(cpaf = mean(cpaf)), 
                               by = .(age, sex)]
cpaf_cancer_cvd_mean <- paf_comorb[, .(cpaf = mean(cpaf)), 
                                         by = .(age, sex)]


ir_mean <- incid_rate[, .(ir_cvd = mean(ir_cvd), ir_cancer = mean(ir_cancer),
                    ir_cvd_cancer = mean(ir_cvd_cancer), ir_cancer_cvd = mean(ir_cancer_cvd),
                    ir_comorb = mean(ir_comorb)), by=c('age', 'sex')]



# compute the calibration constants
# If the age >= age_cutoff then calibrate to the observed cpaf
# If the age < age_cutof, then the calibration is 0 and all cases come from cpaf_other
constants_df <- ir_mean
constants_df <- CJ(age = 0:100, sex = 1:2)

for (sex_ in 1:2){
  pf_c <- cpaf_cancer_mean[sex == sex_ & age >= age_cutoff_cancer & age <= age_cutoff_cancer_high]
  pf_cvd <- cpaf_cvd_mean[sex == sex_ & age >= age_cutoff_cvd & age <= age_cutoff_cvd_high]
  pf_c_cvd <- cpaf_cancer_cvd_mean[sex == sex_ & age >= max(age_cutoff_cvd, age_cutoff_cancer) & age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high)]
  
  if(calibrate_cpaf_cancer == 0) {
  c_c <- ir_mean[sex == sex_ & age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, ir_cancer]
  c_cvd_c <- ir_mean[sex == sex_ & age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, ir_cvd_cancer]
  }
  else{
  c_c <- (prop_cancer * ir_mean[sex == sex_ & age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, ir_cancer]) / mean(pf_c$cpaf)
  c_cvd_c <- (prop_cancer * ir_mean[sex == sex_ & age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, ir_cvd_cancer]) / mean(pf_c$cpaf)

  }
  constants_df[sex == sex_ &
                 age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, const_cancer := c_c]
  constants_df[sex == sex_ &
                 age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, const_cvd_comorb := c_cvd_c]
  
  if(calibrate_cpaf_cvd == 0) {
    c_cvd <- ir_mean[sex == sex_ & age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, ir_cvd]
    c_c_cvd <- ir_mean[sex == sex_ & age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, ir_cancer_cvd]
    c_c_o_cvd <- ir_mean[sex == sex_ & age >= max(age_cutoff_cvd, age_cutoff_cancer) & age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high), ir_comorb]
  }
  else{
    c_cvd <- (prop_cvd * ir_mean[sex == sex_ & age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, ir_cvd]) / mean(pf_cvd$cpaf)
    c_c_cvd <- (prop_cvd * ir_mean[sex == sex_ & age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, ir_cancer_cvd]) / mean(pf_cvd$cpaf)
    c_c_o_cvd <- (mean(prop_cvd, prop_cancer)* ir_mean[sex == sex_ & age >= max(age_cutoff_cvd, age_cutoff_cancer) & age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high), ir_comorb] ) /
      mean(pf_c_cvd$cpaf)
  }
  constants_df[sex == sex_ &
                 age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, const_cvd := c_cvd]
  constants_df[sex == sex_ &
                 age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, const_cancer_comorb := c_c_cvd]
  
  constants_df[sex == sex_ &
                 age >= max(age_cutoff_cvd, age_cutoff_cancer) & age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high), const_comorb := c_c_o_cvd]
  
}
setnafill(constants_df, "const", fill=0.0, cols=c("const_cancer", "const_cvd",
                                                  "const_cancer_comorb", "const_cvd_comorb",
                                                  "const_comorb"))


# Calculate paf_other:
# if age >= age_cutoff then paf_other is the complement to paf from lifestyle
# if age < age_cutoff then paf_other = incidence rate
paf_other <- CJ(sex = 1:2, age= 0:100,
                       year = as.character(list_of_years))

for (sex_ in 1:2){
  for (y in list_of_years[2:length(list_of_years)]){
    c_c <- constants_df[sex == sex_ &
                          age >= age_cutoff_cancer  & age <= age_cutoff_cancer_high,
                        const_cancer]
    ir_c <- ir_mean[sex == sex_ &
                         age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, ir_cancer]
    c_cvd <- constants_df[sex == sex_ &
                          age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, const_cvd]
    ir_cvd <- ir_mean[sex == sex_ &
                         age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, ir_cvd]
    c_cvd_cancer <- constants_df[sex == sex_ &
                            age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, const_cvd_comorb]
    ir_cvd_cancer <- ir_mean[sex == sex_ &
                           age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, ir_cvd_cancer]
    c_c_cvd <- constants_df[sex == sex_ & 
                          age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, const_cancer_comorb]
    ir_c_cvd <- ir_mean[sex == sex_ &
                         age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, ir_cancer_cvd]
    c_c_o_cvd <- constants_df[sex == sex_ &
                                age >= max(age_cutoff_cvd, age_cutoff_cancer) & age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high), const_comorb]
    ir_c_o_cvd <- ir_mean[sex == sex_ &
                            age >= max(age_cutoff_cvd, age_cutoff_cancer) & age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high), ir_comorb]
    
    
    cut_c <- paf_cancer[year == y & sex == sex_,]
    paf_slice_c <- cut_c[age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, cpaf]
    cut_cvd <- paf_cvd[year == y & sex == sex_,]
    paf_slice_cvd <- cut_cvd[age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, cpaf]
    cut_c_o_cvd <- paf_comorb[year == y & sex == sex_,]
    paf_slice_c_o_cvd <- cut_c_o_cvd[age >= max(age_cutoff_cvd, age_cutoff_cancer) & age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high), cpaf]
    
    lb_c <- ir_c - c_c*paf_slice_c
    res_c <- ir_mean[sex == sex_ & 
                          (age < age_cutoff_cancer | age > age_cutoff_cancer_high), ir_cancer]
    paf_other[year == y & sex == sex_ &
                       age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, cpaf_cancer := lb_c]
    paf_other[year == y & sex == sex_ &
                       (age < age_cutoff_cancer | age > age_cutoff_cancer_high), cpaf_cancer := res_c]
    lb_cvd <- ir_cvd - c_cvd*paf_slice_cvd
    res_cvd <- ir_mean[sex == sex_ &
                       (age < age_cutoff_cvd | age > age_cutoff_cvd_high), ir_cvd]
    paf_other[year == y & sex == sex_ &
                age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, cpaf_cvd := lb_cvd]
    paf_other[year == y & sex == sex_ &
                (age < age_cutoff_cvd | age > age_cutoff_cvd_high), cpaf_cvd := res_cvd]
    lb_cancer_cvd <- ir_c_cvd - c_c_cvd*paf_slice_cvd
    res_cancer_cvd <- ir_mean[sex == sex_ &
                         (age < age_cutoff_cvd | age > age_cutoff_cvd_high), ir_cancer_cvd]
    paf_other[year == y & sex == sex_ &
                age >= age_cutoff_cvd & age <= age_cutoff_cvd_high, cpaf_cancer_comorb := lb_cancer_cvd]
    paf_other[year == y & sex == sex_ &
                (age < age_cutoff_cvd | age > age_cutoff_cvd_high), cpaf_cancer_comorb := res_cancer_cvd]
    lb_cvd_cancer <- ir_cvd_cancer - c_cvd_c*paf_slice_c
    res_cvd_cancer <- ir_mean[sex == sex_ &
                             (age < age_cutoff_cancer | age > age_cutoff_cancer_high), ir_cvd_cancer]
    paf_other[year == y & sex == sex_ &
                age >= age_cutoff_cancer & age <= age_cutoff_cancer_high, cpaf_cvd_comorb := lb_cvd_cancer]
    paf_other[year == y & sex == sex_ &
                (age < age_cutoff_cancer | age > age_cutoff_cancer_high), cpaf_cvd_comorb := res_cvd_cancer]
    
    lb_comorb <- ir_c_o_cvd - c_c_o_cvd*paf_slice_c_o_cvd
    res_comorb <- ir_mean[sex == sex_ &
                             (age < max(age_cutoff_cvd, age_cutoff_cancer) | age > max(age_cutoff_cvd_high, age_cutoff_cancer_high)), ir_comorb]
    paf_other[year == y & sex == sex_ &
                age >= max(age_cutoff_cvd, age_cutoff_cancer) & age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high), cpaf_comorb := lb_comorb]
    paf_other[year == y & sex == sex_ &
                (age < max(age_cutoff_cvd, age_cutoff_cancer) | age > max(age_cutoff_cvd_high, age_cutoff_cancer_high)), cpaf_comorb := res_comorb]
    
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
  tot_sick <- sum(preval[year == dcost_total_base_year, ((names(preval) %like% ncd_) | (names(preval) %like% "comorb")), with=FALSE])
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
    prev <- preval[year == dcost_total_base_year & sex == s, ((names(preval) %like% ncd_) | (names(preval) %like% "comorb")), with=FALSE]
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
    prev <- preval[year == dcost_total_base_year & sex == s, ((names(preval) %like% ncd_) | (names(preval) %like% "comorb")), with=FALSE]
    tot_p <- sum(prev)
    del_cost <- (indirect_costs[[ncd_]]/tot_sick)  * tot_p 
    ik_i <- c_m * del_cost / rowSums(prev)
    ik_i[is.na(ik_i)] <- 0
    out_ik_0[sex == s, icost_unit := ik_i]
  }
  
  out_kpp <- rbind(out_kpp, out_kpp_0)
  out_ik <- rbind(out_ik, out_ik_0)

}

return(list(constants = melt(constants_df, id.vars = c("sex", "age"), variable.name = 'disease', value.name = 'constant' ),
            heal_rate = melt(heal_rate, id.vars = c("sex", "age"), variable.name = 'disease', value.name = 'heal_rate' ), 
            paf_other = melt(paf_other, id.vars = c("sex", "age", "year"), variable.name = 'disease', value.name = 'cpaf'),
            dcost_unit = out_kpp, icost_unit = out_ik))

}
