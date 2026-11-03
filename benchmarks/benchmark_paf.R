library(bench)

paf_optimized <- function(p, rr) {
  x <- p * (rr - 1)
  return(x / (x + 1))
}

paf_original <- function(p, rr) {
  return((p * (rr - 1)) / (p * (rr - 1) + 1))
}

p <- runif(1e6)
rr <- runif(1e6)

bm <- bench::mark(
  original = paf_original(p, rr),
  optimized = paf_optimized(p, rr),
  min_iterations = 1e4
)

print(bm)
