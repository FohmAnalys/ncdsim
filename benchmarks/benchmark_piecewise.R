library(bench)

# Define the original piecewise function
piecewise <- function(x, x0, y0, x1, y1) {
  a <- (y1 - y0) / (x1 - x0)
  b <- y0 - x0 * a
  f <- (x < x0) * y0 + (x >= x0) * (x < x1) * (a * x + b) + (x >= x1) * y1
  return(f)
}

# Define the optimized piecewise function using approxfun
piecewise_optimized <- function(x, x0, y0, x1, y1) {
  # Use approxfun for interpolation
  f <- approxfun(x = c(x0, x1), c(y0, y1), yleft = y0, yright = y1, rule = 2)
  return(f(x))
}

x <- seq(-10, 20, length.out = 10000)
x0 <- 0
y0 <- 0
x1 <- 10
y1 <- 10


# Benchmark
bm <- bench::mark(
  original = piecewise(x, x0, y0, x1, y1),
  approxfun = piecewise_optimized(x, x0, y0, x1, y1),
  min_iterations = 1000
)

print(bm)
