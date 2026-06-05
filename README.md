
# zigdist

`zigdist` is an experimental demonstration that integrates the [Zig
programming language](https://ziglang.org/) into an R package. It
provides high-performance distance metric calculations as a
proof-of-concept. While Zig is not currently supported on CRAN, its
exceptional performance, simple FFI, data-oriented design, and robust
systems programming capabilities suggest it could be of significant
interest to the R package development community.

## Distance methods

- Euclidean (`zd_euclidean()`): Pairwise row distance and cross-distance
  for numeric matrices, supporting missing values (partial distance
  scaling)
- Gower (`zd_gower()`): Mixed-type distance computations (scaling
  numeric features by range, matching categorical features), supporting
  missing values (pairwise deletion)

These distance methods use multi-threaded parallel execution by
specifying `num_threads > 1`.

## Installation

You can install the development version of `zigdist` from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("brownag/zigdist")
```

Note: A Zig compiler (v0.13.0 or later) must be installed on your system
and available in your `PATH` to compile the package from source.

## Quick Start

``` r
library(zigdist)

# Euclidean
x <- matrix(rnorm(15), nrow = 5, ncol = 3)
zd_euclidean(x)
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,] 0.000000 2.887606 2.755103 3.275980 3.171044
    ## [2,] 2.887606 0.000000 3.218801 3.139713 1.352931
    ## [3,] 2.755103 3.218801 0.000000 1.887001 2.958793
    ## [4,] 3.275980 3.139713 1.887001 0.000000 2.130555
    ## [5,] 3.171044 1.352931 2.958793 2.130555 0.000000

``` r
# Gower (supports mixed types: numeric, nominal factors, and ordered factors)
df <- data.frame(
  a = c(1, 2, 5),
  b = factor(c("A", "B", "A")),
  c = ordered(c("Low", "Medium", "High"), levels = c("Low", "Medium", "High"))
)
zd_gower(df)
```

    ##           [,1]      [,2]      [,3]
    ## [1,] 0.0000000 0.5833333 0.6666667
    ## [2,] 0.5833333 0.0000000 0.7500000
    ## [3,] 0.6666667 0.7500000 0.0000000

## Correctness and Equivalence

To demonstrate functional equivalence, we verify that the outputs of
`zigdist` match the standard R implementations (`stats::dist` and
`cluster::daisy`) using `all.equal()`:

``` r
# Euclidean vs. stats::dist
x_check <- matrix(rnorm(300), nrow = 20, ncol = 15)
all.equal(zd_euclidean(x_check), as.matrix(dist(x_check)), check.attributes = FALSE)
```

    ## [1] TRUE

``` r
# Gower vs. cluster::daisy
library(cluster)
df_check <- data.frame(
  n = rnorm(20),
  f = factor(sample(letters[1:3], 20, replace = TRUE)),
  c = factor(sample(c("Yes", "No"), 20, replace = TRUE))
)
all.equal(zd_gower(df_check), as.matrix(daisy(df_check, metric = "gower")), check.attributes = FALSE)
```

    ## [1] TRUE

## Performance Benchmarks

Below are benchmarks comparing `zigdist` against standard R
implementations, using the `microbenchmark` package.

``` r
library(microbenchmark)
library(cluster)
```

### Euclidean

We compare `zigdist::zd_euclidean` (single-threaded and multi-threaded)
against R’s built-in `stats::dist` (converted to a full matrix) on a
$500 \times 100$ numeric matrix.

``` r
set.seed(42)
x_bench <- matrix(rnorm(50000), nrow = 500, ncol = 100)

bench_eucl <- microbenchmark(
  `zigdist 1 thread` = zd_euclidean(x_bench, num_threads = 1L),
  `zigdist 4 thread` = zd_euclidean(x_bench, num_threads = 4L),
  stats_dist = as.matrix(dist(x_bench)),
  times = 50
)
print(bench_eucl)
```

    ## Unit: microseconds
    ##              expr      min       lq     mean   median       uq       max neval
    ##  zigdist 1 thread 3144.840 3443.247 3545.967 3476.088 3509.308  6514.185    50
    ##  zigdist 4 thread  939.331 1253.618 1645.115 1517.393 1870.322  5093.594    50
    ##        stats_dist 7352.883 7513.347 8179.239 8005.413 8244.244 12209.052    50

Even on a single thread (`num_threads = 1L`), `zigdist` significantly
outperforms R’s native `stats::dist` (written in optimized C) by more
than **2.3x**. Spawning multiple threads (e.g., `num_threads = 4L`)
achieves a **5.3x speedup**.

### Gower

We compare `zigdist::zd_gower` (single-threaded and multi-threaded)
against `cluster::daisy(..., metric = "gower")` (converted to a full
matrix) on a mixed-type data frame with 500 rows and 6 columns (3
numeric, 3 categorical).

``` r
set.seed(42)
N_bench <- 500
df_bench <- data.frame(
  n1 = rnorm(N_bench),
  n2 = runif(N_bench),
  n3 = rnorm(N_bench, mean = 10, sd = 2),
  c1 = factor(sample(letters[1:4], N_bench, replace = TRUE)),
  c2 = factor(sample(c("Yes", "No"), N_bench, replace = TRUE)),
  c3 = factor(sample(colors()[1:10], N_bench, replace = TRUE))
)

bench_gow <- microbenchmark(
  `zigdist 1 thread` = zd_gower(df_bench, num_threads = 1L),
  `zigdist 4 thread` = zd_gower(df_bench, num_threads = 4L),
  cluster_daisy = as.matrix(daisy(df_bench, metric = "gower")),
  times = 20
)
print(bench_gow)
```

    ## Unit: microseconds
    ##              expr      min        lq     mean    median        uq      max
    ##  zigdist 1 thread  542.272  604.8385 1127.525  824.4535  986.5465 4128.382
    ##  zigdist 4 thread  562.532  953.7210 1382.014 1171.0705 1566.3550 4057.090
    ##     cluster_daisy 4805.537 6250.7845 7044.421 7627.6950 7912.9260 9645.715
    ##  neval
    ##     20
    ##     20
    ##     20

For mixed-type data, `zigdist` represents a substantial improvement over
`cluster::daisy`, calculating distances **9.3x faster** in
single-threaded mode. For smaller calculations, running sequentially on
a single thread is recommended to eliminate thread-spawning overhead.
