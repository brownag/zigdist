# Compute Pairwise Euclidean Distance

Computes the pairwise Euclidean distance matrix for rows of a matrix x,
or the cross-distance matrix between rows of x and rows of y. Supports
missing values (NA) using partial distance scaling.

## Usage

``` r
zd_euclidean(x, y = NULL, num_threads = 1L)
```

## Arguments

- x:

  A numeric matrix or data frame of size `N_x x P`.

- y:

  An optional numeric matrix or data frame of size `N_y x P`.

- num_threads:

  An integer specifying the number of threads to use. Defaults to `1L`
  (single-threaded).

## Value

A numeric matrix of distances. If `y` is `NULL`, returns an `N_x x N_x`
matrix. If `y` is provided, returns an `N_x x N_y` matrix.

## Examples

``` r
x <- matrix(rnorm(15), nrow = 5, ncol = 3)
zd_euclidean(x, num_threads = 1L)
#>          [,1]     [,2]     [,3]     [,4]     [,5]
#> [1,] 0.000000 3.600169 3.143505 2.245990 2.696585
#> [2,] 3.600169 0.000000 3.433837 2.768468 1.586373
#> [3,] 3.143505 3.433837 0.000000 4.424213 3.430477
#> [4,] 2.245990 2.768468 4.424213 0.000000 2.233607
#> [5,] 2.696585 1.586373 3.430477 2.233607 0.000000

y <- matrix(rnorm(9), nrow = 3, ncol = 3)
zd_euclidean(x, y, num_threads = 2L)
#>          [,1]      [,2]     [,3]
#> [1,] 1.192105 2.3640126 1.986519
#> [2,] 3.185974 2.2730507 2.313141
#> [3,] 1.962506 3.9327768 2.821367
#> [4,] 2.836470 0.9066484 2.474007
#> [5,] 2.622439 2.2378004 1.033952
```
