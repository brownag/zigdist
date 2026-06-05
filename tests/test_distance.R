library(tinytest)
library(zigdist)

# --- 1. Test Euclidean Distance ---

# Pairwise self-distance comparison
set.seed(42)
x <- matrix(rnorm(50), nrow = 10, ncol = 5)
res_self <- zd_euclidean(x)
ref_self <- as.matrix(dist(x))
expect_equal(unname(res_self), unname(ref_self), tolerance = 1e-7)

# Diagonal is exactly zero
expect_equal(diag(res_self), rep(0, 10))

# Cross-distance comparison
y <- matrix(rnorm(40), nrow = 8, ncol = 5)
res_cross <- zd_euclidean(x, y)

# Manual reference calculation for cross-distance
ref_cross <- matrix(0, nrow = 10, ncol = 8)
for (i in 1:10) {
  for (j in 1:8) {
    ref_cross[i, j] <- sqrt(sum((x[i, ] - y[j, ])^2))
  }
}
expect_equal(res_cross, ref_cross, tolerance = 1e-7)

# Dimension check
expect_equal(dim(res_self), c(10, 10))
expect_equal(dim(res_cross), c(10, 8))

# Error on column mismatch
expect_error(zd_euclidean(x, y[, 1:4]))

# Error on non-numeric inputs
expect_error(zd_euclidean(matrix("a", 2, 2)))


# --- 2. Test Euclidean Distance with Missing Values (NA) ---

# Partial distance scaling: d = sqrt( sum(diff^2) * P / P_present )
# Row 1: 1, NA, 3.  Row 2: 4, 5, 6.
# present: 1 vs 4 (diff^2=9), 3 vs 6 (diff^2=9). P_present = 2, P = 3.
# sum = 18. scaled = 18 * 3/2 = 27. dist = sqrt(27) = 5.196152
x_na <- matrix(c(1, NA, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)
res_na <- zd_euclidean(x_na)
expect_equal(res_na[1, 2], sqrt(27), tolerance = 1e-7)

# Completely missing row comparison yields NaN
x_all_na <- matrix(c(1, 2, 3, NA, NA, NA), nrow = 2, ncol = 3, byrow = TRUE)
res_all_na <- zd_euclidean(x_all_na)
expect_true(is.nan(res_all_na[1, 2]))


# --- 3. Test Gower Distance ---

# Setup test dataset with mixed types
df_x <- data.frame(
  a = c(1, 2, 5),
  b = factor(c("A", "B", "A")),
  c = c(TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

# Hand-calculated self Gower distance matrix (equal weights = c(1, 1, 1))
expected_gower_self <- matrix(c(
  0.0,   0.75,      1.0/3.0,
  0.75,  0.0,       2.75/3.0,
  1.0/3.0, 2.75/3.0, 0.0
), nrow = 3, ncol = 3, byrow = TRUE)

res_gower_self <- zd_gower(df_x)
expect_equal(res_gower_self, expected_gower_self, tolerance = 1e-7)

# Gower cross distance setup
df_y <- data.frame(
  a = c(2, 4),
  b = factor(c("A", "C")),
  c = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

# Hand-calculated Gower cross distance matrix (equal weights)
expected_gower_cross <- matrix(c(
  1.25/3.0, 1.75/3.0,
  1.00/3.0, 2.50/3.0,
  1.75/3.0, 1.25/3.0
), nrow = 3, ncol = 2, byrow = TRUE)

res_gower_cross <- zd_gower(df_x, df_y)
expect_equal(res_gower_cross, expected_gower_cross, tolerance = 1e-7)

# Gower distance with custom weights
res_gower_weighted <- zd_gower(df_x, weight = c(2, 0.5, 0.5))
expect_equal(res_gower_weighted[1, 2], 0.5, tolerance = 1e-7)

# Gower distance with zero-range columns
df_const <- data.frame(
  a = c(1, 2, 5),
  b = c(10, 10, 10)
)
res_const <- zd_gower(df_const)
expect_equal(res_const[1, 2], 0.25, tolerance = 1e-7)
expect_equal(res_const[1, 3], 1.00, tolerance = 1e-7)

# Gower input checks
expect_error(zd_gower(df_x, weight = c(1, 1)))
expect_error(zd_gower(df_x, df_y[, 1:2]))
expect_error(zd_gower(df_x, data.frame(a=1:3, b=1:3, c=1:3)))


# --- 4. Test Gower Distance with Ordinal Factors ---

# Ordered factor self-distance
df_ord_self <- data.frame(
  a = ordered(c("Low", "Medium", "High"), levels = c("Low", "Medium", "High")),
  b = c(1, 2, 5)
)
expected_ord_self <- matrix(c(
  0.000, 0.375, 1.000,
  0.375, 0.000, 0.625,
  1.000, 0.625, 0.000
), nrow = 3, ncol = 3, byrow = TRUE)

res_ord_self <- zd_gower(df_ord_self)
expect_equal(res_ord_self, expected_ord_self, tolerance = 1e-7)
expect_equal(unname(res_ord_self), unname(as.matrix(cluster::daisy(df_ord_self, metric = "gower"))), tolerance = 1e-7)

# Ordered factor cross-distance
df_ord_cross_y <- data.frame(
  a = ordered(c("Medium", "High"), levels = c("Low", "Medium", "High")),
  b = c(2, 4)
)
expected_ord_cross <- matrix(c(
  0.375, 0.875,
  0.000, 0.500,
  0.625, 0.125
), nrow = 3, ncol = 2, byrow = TRUE)

res_ord_cross <- zd_gower(df_ord_self, df_ord_cross_y)
expect_equal(res_ord_cross, expected_ord_cross, tolerance = 1e-7)

# Check errors on mismatch of ordered factors vs non-ordered/numeric
expect_error(zd_gower(df_ord_self, data.frame(a = factor(c("Low", "Medium")), b = c(2, 4))))
expect_error(zd_gower(df_ord_self, data.frame(a = c(1, 2), b = c(2, 4))))


# --- 5. Test Gower Distance with Missing Values (NA) ---

# df_na:
# Row 1: (1.0, "A", TRUE)
# Row 2: (NA,  "B", FALSE) -> Col 'a' is NA (ignored)
# Row 3: (5.0, NA,  TRUE)  -> Col 'b' is NA (ignored)
# Range of a: 4
# dist(1, 2): a ignored. b differs (1), c differs (1). dist = 2/2 = 1.0
# dist(1, 3): a diff = 1. b ignored. c same (0). dist = 1/2 = 0.5
# dist(2, 3): a ignored. b ignored. c differs (1). dist = 1/1 = 1.0
df_na <- data.frame(
  a = c(1, NA, 5),
  b = factor(c("A", "B", NA)),
  c = c(TRUE, FALSE, TRUE)
)
res_na_gow <- zd_gower(df_na)
expected_na_gow <- matrix(c(
  0.0, 1.0, 0.5,
  1.0, 0.0, 1.0,
  0.5, 1.0, 0.0
), nrow = 3, ncol = 3, byrow = TRUE)
expect_equal(res_na_gow, expected_na_gow, tolerance = 1e-7)


# --- 6. Test Thread Control (num_threads) ---

# Verify that 1 thread and 2 threads yield identical results
x_thread <- matrix(rnorm(100), nrow = 20, ncol = 5)
res_t1 <- zd_euclidean(x_thread, num_threads = 1L)
res_t2 <- zd_euclidean(x_thread, num_threads = 2L)
expect_equal(res_t1, res_t2, tolerance = 1e-7)

df_thread <- data.frame(
  a = rnorm(20),
  b = factor(sample(c("X", "Y"), 20, replace = TRUE))
)
res_g_t1 <- zd_gower(df_thread, num_threads = 1L)
res_g_t2 <- zd_gower(df_thread, num_threads = 2L)
expect_equal(res_g_t1, res_g_t2, tolerance = 1e-7)

# Verify invalid thread count errors
expect_error(zd_euclidean(x_thread, num_threads = 0L))
expect_error(zd_euclidean(x_thread, num_threads = -5L))
expect_error(zd_euclidean(x_thread, num_threads = "two"))
expect_error(zd_gower(df_thread, num_threads = 0L))

