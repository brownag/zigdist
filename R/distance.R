#' @useDynLib zigdist, .registration = TRUE
NULL

# Helper: coerce to numeric matrix
as_numeric_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix or data frame")
  }
  x
}

# Helper: validate and coerce num_threads
validate_num_threads <- function(num_threads) {
  if (is.null(num_threads)) return(1L)
  num_threads <- suppressWarnings(as.integer(num_threads))
  if (is.na(num_threads) || num_threads < 1L) {
    stop("num_threads must be a positive integer >= 1")
  }
  num_threads
}

#' Compute Pairwise Euclidean Distance
#'
#' Computes the pairwise Euclidean distance matrix for rows of a matrix x,
#' or the cross-distance matrix between rows of x and rows of y.
#' Supports missing values (NA) using partial distance scaling.
#'
#' @param x A numeric matrix or data frame of size \code{N_x x P}.
#' @param y An optional numeric matrix or data frame of size \code{N_y x P}.
#' @param num_threads An integer specifying the number of threads to use.
#'   Defaults to \code{1L} (single-threaded).
#' @return A numeric matrix of distances. If \code{y} is \code{NULL}, returns an
#'   \code{N_x x N_x} matrix. If \code{y} is provided, returns an \code{N_x x N_y} matrix.
#' @export
#' @examples
#' x <- matrix(rnorm(15), nrow = 5, ncol = 3)
#' zd_euclidean(x, num_threads = 1L)
#'
#' y <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' zd_euclidean(x, y, num_threads = 2L)
zd_euclidean <- function(x, y = NULL, num_threads = 1L) {
  x <- as_numeric_matrix(x)
  num_threads <- validate_num_threads(num_threads)
  has_na <- anyNA(x)

  if (is.null(y)) {
    return(.Call(zd_euclidean_dist_self_R, x, num_threads, has_na))
  }

  y <- as_numeric_matrix(y)
  if (ncol(x) != ncol(y)) {
    stop("x and y must have the same number of columns")
  }
  has_na <- has_na || anyNA(y)
  .Call(zd_euclidean_dist_cross_R, x, y, num_threads, has_na)
}

#' Compute Gower Distance
#'
#' Computes the Gower (1971) distance matrix for mixed-type datasets.
#' Supports missing values (NA) by setting their weight contribution to zero.
#'
#' @param x A data frame or matrix of size \code{N_x x P}.
#' @param y An optional data frame or matrix of size \code{N_y x P}.
#' @param weight An optional numeric vector of weights for each column.
#'   If \code{NULL}, all columns are weighted equally.
#' @param num_threads An integer specifying the number of threads to use.
#'   Defaults to \code{1L} (single-threaded).
#' @return A numeric matrix of Gower distances. If \code{y} is \code{NULL}, returns an
#'   \code{N_x x N_x} matrix. If \code{y} is provided, returns an \code{N_x x N_y} matrix.
#' @references
#' Gower, J. C. (1971). A general coefficient of similarity and some of its properties.
#' \emph{Biometrics}, 27(4), 857–871. \doi{10.2307/2528823}
#'
#' Podani, J. (1999). Extending Gower's general coefficient of similarity to ordinal characters.
#' \emph{Taxon}, 48(2), 331–340. \doi{10.2307/1224024}
#' @export
#' @examples
#' df_x <- data.frame(
#'   a = c(1, 2, 5),
#'   b = factor(c("A", "B", "A")),
#'   c = c(TRUE, FALSE, TRUE)
#' )
#' zd_gower(df_x, num_threads = 1L)
#'
#' df_y <- data.frame(
#'   a = c(2, 4),
#'   b = factor(c("A", "C")),
#'   c = c(FALSE, TRUE)
#' )
#' zd_gower(df_x, df_y, weight = c(2.0, 0.5, 0.5), num_threads = 2L)
zd_gower <- function(x, y = NULL, weight = NULL, num_threads = 1L) {
  x <- as.data.frame(x)
  N_x <- nrow(x)
  P <- ncol(x)

  if (is.null(weight)) {
    weight <- rep(1.0, P)
  } else {
    if (length(weight) != P) {
      stop("weight must have length equal to number of columns in x")
    }
    weight <- as.numeric(weight)
    if (any(is.na(weight)) || any(weight < 0)) {
      stop("weights must be non-negative numbers")
    }
  }

  num_threads <- validate_num_threads(num_threads)

  # Determine column types
  is_num <- sapply(x, is.numeric)
  is_ord <- sapply(x, is.ordered)

  # Pre-process numeric variables for x (including converting ordered factors to ranks)
  for (j in which(is_ord)) {
    x[[j]] <- as.numeric(x[[j]])
  }

  is_num_or_ord <- is_num | is_ord
  x_num <- as.matrix(x[, is_num_or_ord, drop = FALSE])
  mode(x_num) <- "double"

  # Pre-process categorical variables
  cat_indices <- which(!is_num_or_ord)
  P_cat <- length(cat_indices)

  if (is.null(y)) {
    # Self distance
    if (P_cat > 0) {
      # Represent categorical variables as an integer matrix of factor levels
      # Convert factors while preserving missing values as NA
      x_cat_list <- lapply(x[, !is_num_or_ord, drop = FALSE], function(col) {
        if (any(is.na(col))) {
          codes <- as.integer(as.factor(col))
          codes[is.na(col)] <- NA_integer_
          return(codes)
        } else {
          return(as.integer(as.factor(col)))
        }
      })
      x_cat <- do.call(cbind, x_cat_list)
      mode(x_cat) <- "integer"
    } else {
      x_cat <- matrix(integer(0), nrow = N_x, ncol = 0)
    }

    # Calculate range of each numeric column of x (ignoring NA values)
    ranges <- numeric(ncol(x_num))
    if (ncol(x_num) > 0) {
      for (j in seq_len(ncol(x_num))) {
        r_vals <- range(x_num[, j], na.rm = TRUE)
        r <- r_vals[2] - r_vals[1]
        ranges[j] <- if (is.na(r) || !is.finite(r)) 0.0 else r
      }
    }

    w_num <- weight[is_num_or_ord]
    w_cat <- weight[!is_num_or_ord]

    has_na <- anyNA(x_num) || (P_cat > 0 && anyNA(x_cat))
    return(.Call(zd_gower_dist_self_R, x_num, x_cat, ranges, w_num, w_cat, num_threads, has_na))

  } else {
    # Cross distance
    y <- as.data.frame(y)
    if (ncol(y) != P) {
      stop("x and y must have the same number of columns")
    }
    N_y <- nrow(y)

    # Check that column types match between x and y
    is_num_y <- sapply(y, is.numeric)
    is_ord_y <- sapply(y, is.ordered)
    if (!all(is_num == is_num_y) || !all(is_ord == is_ord_y)) {
      stop("Columns in x and y must have matching types (numeric/ordered/categorical)")
    }

    # Align levels of ordered factors across x and y before converting to ranks
    for (j in which(is_ord)) {
      if (!identical(levels(x[[j]]), levels(y[[j]]))) {
        levels_combined <- union(levels(x[[j]]), levels(y[[j]]))
        x[[j]] <- ordered(as.character(x[[j]]), levels = levels_combined)
        y[[j]] <- ordered(as.character(y[[j]]), levels = levels_combined)
      }
      x[[j]] <- as.numeric(x[[j]])
      y[[j]] <- as.numeric(y[[j]])
    }

    y_num <- as.matrix(y[, is_num_or_ord, drop = FALSE])
    mode(y_num) <- "double"

    # Align factors/categories between x and y
    if (P_cat > 0) {
      x_cat_list <- list()
      y_cat_list <- list()
      for (idx in cat_indices) {
        # Unify factor levels across x and y, preserving NAs
        col_x <- x[[idx]]
        col_y <- y[[idx]]
        combined <- factor(c(as.character(col_x), as.character(col_y)))
        
        codes_x <- as.integer(combined[1:N_x])
        codes_y <- as.integer(combined[(N_x + 1):(N_x + N_y)])
        
        codes_x[is.na(col_x)] <- NA_integer_
        codes_y[is.na(col_y)] <- NA_integer_
        
        x_cat_list[[length(x_cat_list) + 1]] <- codes_x
        y_cat_list[[length(y_cat_list) + 1]] <- codes_y
      }
      x_cat <- do.call(cbind, x_cat_list)
      y_cat <- do.call(cbind, y_cat_list)
      mode(x_cat) <- "integer"
      mode(y_cat) <- "integer"
    } else {
      x_cat <- matrix(integer(0), nrow = N_x, ncol = 0)
      y_cat <- matrix(integer(0), nrow = N_y, ncol = 0)
    }

    # Ranges are calculated over the combined numeric columns of x and y
    ranges <- numeric(ncol(x_num))
    if (ncol(x_num) > 0) {
      for (j in seq_len(ncol(x_num))) {
        combined_vals <- c(x_num[, j], y_num[, j])
        r_vals <- range(combined_vals, na.rm = TRUE)
        r <- r_vals[2] - r_vals[1]
        ranges[j] <- if (is.na(r) || !is.finite(r)) 0.0 else r
      }
    }

    w_num <- weight[is_num_or_ord]
    w_cat <- weight[!is_num_or_ord]

    has_na <- anyNA(x_num) || (P_cat > 0 && anyNA(x_cat)) || anyNA(y_num) || (P_cat > 0 && anyNA(y_cat))
    return(.Call(zd_gower_dist_cross_R, x_num, x_cat, y_num, y_cat, ranges, w_num, w_cat, num_threads, has_na))
  }
}
