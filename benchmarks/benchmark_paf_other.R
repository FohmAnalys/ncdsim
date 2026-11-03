library(bench)

calculate_paf_other_original <- function() {
    paf_other <- CJ(sex = 1:2, age = 0:100, year = as.character(list_of_years))

    for (sex_ in 1:2) {
        for (y in list_of_years[2:length(list_of_years)]) {
            c_c <- constants_df[
                sex == sex_ &
                    age >= age_cutoff_cancer &
                    age <= age_cutoff_cancer_high,
                const_cancer
            ]
            ir_c <- ir_mean[
                sex == sex_ &
                    age >= age_cutoff_cancer &
                    age <= age_cutoff_cancer_high,
                ir_cancer
            ]
            c_cvd <- constants_df[
                sex == sex_ &
                    age >= age_cutoff_cvd &
                    age <= age_cutoff_cvd_high,
                const_cvd
            ]
            ir_cvd <- ir_mean[
                sex == sex_ &
                    age >= age_cutoff_cvd &
                    age <= age_cutoff_cvd_high,
                ir_cvd
            ]
            c_cvd_cancer <- constants_df[
                sex == sex_ &
                    age >= age_cutoff_cancer &
                    age <= age_cutoff_cancer_high,
                const_cvd_comorb
            ]
            ir_cvd_cancer <- ir_mean[
                sex == sex_ &
                    age >= age_cutoff_cancer &
                    age <= age_cutoff_cancer_high,
                ir_cvd_cancer
            ]
            c_c_cvd <- constants_df[
                sex == sex_ &
                    age >= age_cutoff_cvd &
                    age <= age_cutoff_cvd_high,
                const_cancer_comorb
            ]
            ir_c_cvd <- ir_mean[
                sex == sex_ &
                    age >= age_cutoff_cvd &
                    age <= age_cutoff_cvd_high,
                ir_cancer_cvd
            ]
            c_c_o_cvd <- constants_df[
                sex == sex_ &
                    age >= max(age_cutoff_cvd, age_cutoff_cancer) &
                    age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high),
                const_comorb
            ]
            ir_c_o_cvd <- ir_mean[
                sex == sex_ &
                    age >= max(age_cutoff_cvd, age_cutoff_cancer) &
                    age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high),
                ir_comorb
            ]

            cut_c <- paf_cancer[year == y & sex == sex_, ]
            paf_slice_c <- cut_c[
                age >= age_cutoff_cancer & age <= age_cutoff_cancer_high,
                cpaf
            ]
            cut_cvd <- paf_cvd[year == y & sex == sex_, ]
            paf_slice_cvd <- cut_cvd[
                age >= age_cutoff_cvd & age <= age_cutoff_cvd_high,
                cpaf
            ]
            cut_c_o_cvd <- paf_comorb[year == y & sex == sex_, ]
            paf_slice_c_o_cvd <- cut_c_o_cvd[
                age >= max(age_cutoff_cvd, age_cutoff_cancer) &
                    age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high),
                cpaf
            ]

            lb_c <- ir_c - c_c * paf_slice_c
            res_c <- ir_mean[
                sex == sex_ &
                    (age < age_cutoff_cancer | age > age_cutoff_cancer_high),
                ir_cancer
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    age >= age_cutoff_cancer &
                    age <= age_cutoff_cancer_high,
                cpaf_cancer := lb_c
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    (age < age_cutoff_cancer | age > age_cutoff_cancer_high),
                cpaf_cancer := res_c
            ]
            lb_cvd <- ir_cvd - c_cvd * paf_slice_cvd
            res_cvd <- ir_mean[
                sex == sex_ &
                    (age < age_cutoff_cvd | age > age_cutoff_cvd_high),
                ir_cvd
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    age >= age_cutoff_cvd &
                    age <= age_cutoff_cvd_high,
                cpaf_cvd := lb_cvd
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    (age < age_cutoff_cvd | age > age_cutoff_cvd_high),
                cpaf_cvd := res_cvd
            ]
            lb_cancer_cvd <- ir_c_cvd - c_c_cvd * paf_slice_cvd
            res_cancer_cvd <- ir_mean[
                sex == sex_ &
                    (age < age_cutoff_cvd | age > age_cutoff_cvd_high),
                ir_cancer_cvd
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    age >= age_cutoff_cvd &
                    age <= age_cutoff_cvd_high,
                cpaf_cancer_comorb := lb_cancer_cvd
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    (age < age_cutoff_cvd | age > age_cutoff_cvd_high),
                cpaf_cancer_comorb := res_cancer_cvd
            ]
            lb_cvd_cancer <- ir_cvd_cancer - c_cvd_cancer * paf_slice_c
            res_cvd_cancer <- ir_mean[
                sex == sex_ &
                    (age < age_cutoff_cancer | age > age_cutoff_cancer_high),
                ir_cvd_cancer
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    age >= age_cutoff_cancer &
                    age <= age_cutoff_cancer_high,
                cpaf_cvd_comorb := lb_cvd_cancer
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    (age < age_cutoff_cancer | age > age_cutoff_cancer_high),
                cpaf_cvd_comorb := res_cvd_cancer
            ]

            lb_comorb <- ir_c_o_cvd - c_c_o_cvd * paf_slice_c_o_cvd
            res_comorb <- ir_mean[
                sex == sex_ &
                    (age < max(age_cutoff_cvd, age_cutoff_cancer) |
                        age > max(age_cutoff_cvd_high, age_cutoff_cancer_high)),
                ir_comorb
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    age >= max(age_cutoff_cvd, age_cutoff_cancer) &
                    age <= max(age_cutoff_cvd_high, age_cutoff_cancer_high),
                cpaf_comorb := lb_comorb
            ]
            paf_other[
                year == y &
                    sex == sex_ &
                    (age < max(age_cutoff_cvd, age_cutoff_cancer) |
                        age > max(age_cutoff_cvd_high, age_cutoff_cancer_high)),
                cpaf_comorb := res_comorb
            ]
        }
    }
}
