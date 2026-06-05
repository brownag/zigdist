# Compute Gower Distance

Computes the Gower (1971) distance matrix for mixed-type datasets.
Supports missing values (NA) by setting their weight contribution to
zero.

## Usage

``` r
zd_gower(x, y = NULL, weight = NULL, num_threads = 1L)
```

## Arguments

- x:

  A data frame or matrix of size `N_x x P`.

- y:

  An optional data frame or matrix of size `N_y x P`.

- weight:

  An optional numeric vector of weights for each column. If `NULL`, all
  columns are weighted equally.

- num_threads:

  An integer specifying the number of threads to use. Defaults to `1L`
  (single-threaded).

## Value

A numeric matrix of Gower distances. If `y` is `NULL`, returns an
`N_x x N_x` matrix. If `y` is provided, returns an `N_x x N_y` matrix.

## Examples

``` r
df_x <- data.frame(
  a = c(1, 2, 5),
  b = factor(c("A", "B", "A")),
  c = c(TRUE, FALSE, TRUE)
)
zd_gower(df_x, num_threads = 1L)
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.7500000 0.3333333
#> [2,] 0.7500000 0.0000000 0.9166667
#> [3,] 0.3333333 0.9166667 0.0000000

df_y <- data.frame(
  a = c(2, 4),
  b = factor(c("A", "C")),
  c = c(FALSE, TRUE)
)
zd_gower(df_x, df_y, weight = c(2.0, 0.5, 0.5), num_threads = 2L)
#>           [,1]      [,2]
#> [1,] 0.3333333 0.6666667
#> [2,] 0.1666667 0.6666667
#> [3,] 0.6666667 0.3333333
```
