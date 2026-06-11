
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

Note: A Zig compiler (v0.13.0 or later) must be in your `PATH` to
compile the package from source.

The recommended way to install and manage Zig versions is using `zvm`
(Zig Version Manager).

**1. Install `zvm`**

``` bash
curl -fsSL https://www.zvm.app/install.sh | bash
```

This script will install `zvm` to `~/.zvm` and add the necessary
environment variables to your shell profile (e.g., `.bashrc`, `.zshrc`).
You might need to restart your terminal or source your shell profile for
the changes to take effect.

**2. Install a Zig version**

Once `zvm` is installed, you can install the current stable Zig version:

``` bash
zvm install stable
zvm use stable
```

You can also set a default version with `zvm alias default stable`.

## Quick Start

``` r
library(zigdist)

# Euclidean
x <- matrix(rnorm(15), nrow = 5, ncol = 3)
zd_euclidean(x)
```

    ##           [,1]     [,2]      [,3]      [,4]      [,5]
    ## [1,] 0.0000000 2.718995 1.0371808 0.8806271 0.7104837
    ## [2,] 2.7189954 0.000000 2.3454602 3.4465129 2.2092641
    ## [3,] 1.0371808 2.345460 0.0000000 1.4490295 0.5355371
    ## [4,] 0.8806271 3.446513 1.4490295 0.0000000 1.2705891
    ## [5,] 0.7104837 2.209264 0.5355371 1.2705891 0.0000000

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
against R’s built-in `stats::dist` (converted to a full matrix) on a 500
by 100 numeric matrix.

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

    ## Unit: milliseconds
    ##              expr       min        lq      mean    median        uq       max
    ##  zigdist 1 thread  6.429038  7.427399 11.090695  8.315266 15.527293 18.913611
    ##  zigdist 4 thread  2.526398  3.016101  4.354398  3.385486  6.104136  9.398731
    ##        stats_dist 13.214829 13.844654 19.579176 14.417614 28.740439 40.654006
    ##  neval
    ##     50
    ##     50
    ##     50

Even on a single thread (`num_threads = 1L`), `zigdist` significantly
outperforms R’s native `stats::dist` (written in optimized C) by more
than **1.7x**. Spawning multiple threads (e.g., `num_threads = 4L`)
achieves a **4.3x speedup**.

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
    ##              expr      min       lq      mean    median        uq       max
    ##  zigdist 1 thread 1124.568 1192.451  2290.963  1729.693  1867.452  8402.135
    ##  zigdist 4 thread  617.543 1050.959  1771.018  1579.399  1892.815  7103.058
    ##     cluster_daisy 7315.448 9224.679 11691.348 12963.284 13780.539 16096.776
    ##  neval
    ##     20
    ##     20
    ##     20

For mixed-type data, `zigdist` represents a substantial improvement over
`cluster::daisy`, calculating distances **7.5x faster** in
single-threaded mode. For smaller calculations, running sequentially on
a single thread is recommended to eliminate thread-spawning overhead.
