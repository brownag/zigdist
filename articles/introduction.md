# Introduction to zigdist

## Introduction to zigdist

The `zigdist` package provides high-performance, minimal-dependency
implementations of distance metrics in the Zig programming language. It
currently implements:

- **Euclidean Distance** (self-distance and cross-distance for numeric
  matrices/data frames).
- **Gower (1971) Distance** (self-distance and cross-distance for
  mixed-type data frames containing numeric, categorical, and ordinal
  features).

### Euclidean Distance

To calculate the pairwise Euclidean distance matrix of a numeric matrix
`x`:

``` r

library(zigdist)

# Generate a random matrix
set.seed(42)
x <- matrix(rnorm(15), nrow = 5, ncol = 3)

# Calculate pairwise distances (returns 5x5 matrix)
zd_euclidean(x)
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,] 0.000000 2.706923 2.876115 2.750721 1.733427
    ## [2,] 2.706923 0.000000 4.117041 2.876205 3.045218
    ## [3,] 2.876115 4.117041 0.000000 2.402111 1.256619
    ## [4,] 2.750721 2.876205 2.402111 0.000000 2.098702
    ## [5,] 1.733427 3.045218 1.256619 2.098702 0.000000

To calculate the cross-distance matrix between two matrices `x` (size
$`N \times P`$) and `y` (size $`M \times P`$):

``` r

y <- matrix(rnorm(9), nrow = 3, ncol = 3)

# Calculate cross distance (returns 5x3 matrix)
zd_euclidean(x, y)
```

    ##          [,1]     [,2]     [,3]
    ## [1,] 3.938767 2.637191 4.033411
    ## [2,] 5.797243 2.481898 2.971578
    ## [3,] 2.394005 1.975256 3.992646
    ## [4,] 4.705240 1.157653 4.296040
    ## [5,] 2.902283 1.545239 3.353300

### Gower Distance

Gower distance is designed for mixed data types. Numeric and ordinal
columns are scaled by their range:

``` math
d_{ijk} = \frac{|x_{ik} - x_{jk}|}{\text{range}(x_{\cdot k})}
```

Nominal (categorical/logical) columns are compared using simple
matching:

``` math
d_{ijk} = I(x_{ik} \neq x_{jk})
```

The total Gower distance is the weighted average of these components:

``` math
D_{ij} = \frac{\sum_k w_k d_{ijk}}{\sum_k w_k}
```

#### Ordinal Factors (Ordered Categories)

Ordinal factors are natively supported. Following Podani’s (1999) metric
approach, `zd_gower` automatically converts ordered factors to their
integer ranks ($`1, 2, \dots`$) and handles them as range-scaled numeric
columns. This preserves the relative ordering and “closeness” of the
categories (e.g., the distance between `Low` and `Medium` is smaller
than between `Low` and `High`), matching R’s standard
[`cluster::daisy`](https://rdrr.io/pkg/cluster/man/daisy.html) behavior.

#### Example

``` r

# Create a mixed-type data frame with numeric, nominal, and ordinal columns
df_x <- data.frame(
  a = c(1, 2, 5),
  b = factor(c("A", "B", "A")),
  c = ordered(c("Low", "Medium", "High"), levels = c("Low", "Medium", "High"))
)

# Calculate pairwise Gower distances
zd_gower(df_x)
```

    ##           [,1]      [,2]      [,3]
    ## [1,] 0.0000000 0.5833333 0.6666667
    ## [2,] 0.5833333 0.0000000 0.7500000
    ## [3,] 0.6666667 0.7500000 0.0000000

You can also compute cross-distances between two mixed-type datasets,
and provide custom column weights:

``` r

df_y <- data.frame(
  a = c(2, 4),
  b = factor(c("A", "C")),
  c = ordered(c("Medium", "High"), levels = c("Low", "Medium", "High"))
)

# Weighted Gower cross-distance
zd_gower(df_x, df_y, weight = c(2.0, 0.5, 0.5))
```

    ##           [,1]      [,2]
    ## [1,] 0.2500000 0.8333333
    ## [2,] 0.1666667 0.5833333
    ## [3,] 0.5833333 0.3333333

Columns with zero range (i.e. constant values) are automatically ignored
(excluded from the calculations) without raising division-by-zero
errors.
