library(data.table)
library(bench)

dt <- CJ(year = 2015:2115, sex = 1:2, age = 0:100)
dt[, value := runif(.N, min = 0, max = 1000)]

setkeyv(dt, c("year", "sex", "age"))

bm_select <- bench::mark(
  matrix <- dt[, value],
  dollar <- dt$value
)

print(bm_select)

bm_subset <- bench::mark(
  matrix <- dt[.(2020, 1, 50)],
  boolean <- dt[year == 2020 & sex == 1 & age == 50]
)

print(bm_subset)
