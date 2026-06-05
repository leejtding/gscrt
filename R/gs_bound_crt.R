# gs_upper_crt roxy [sinew] ----
#' @title Boundary derivation for efficacy stopping only.
#'
#' @description  \code{gs_upper_crt} is used to calculate the stopping boundaries
#' for a group sequential trial with 1 or 2-sided efficacy-only stopping. Code
#' adapted from gsDesign package.
#'
#' @param theta Effect size under null hypothesis. Default value is 0.
#' @param I Information levels at interim analyses for which stopping
#' boundaries are calculated. At any given interim analysis this should include
#' the information levels at the current and previous analyses.
#' @param a Lower futility boundary values for the interim analyses at the
#' information levels specified by I.
#' @param false_pos Type I error spent for the interim analyses at the
#' information levels specified by I.
#' @param sides \code{1=} 1-sided test (default) \cr \code{2=} 2-sided test
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#' be changed by the user.  This does not translate directly to number of
#' digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#' Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#' values provide larger number of grid points and greater accuracy.  Normally
#' \code{r} will not be changed by the user.
#' @param printerr Print output for debugging.
#'
#' @return Object containing the following elements: \item{k}{Number of interim
#'  analyses.} \item{theta}{As input.} \item{I}{As input.} \item{a}{As input.}
#'  \item{b}{Computed upper efficacy boundaries at the specified interim
#'  analyses.} \item{r}{As input.} \item{error}{Error flag returned; 0 if
#'  convergence; 1 indicates error.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#' @useDynLib gsDesignCRT gs_upper_1
#' @useDynLib gsDesignCRT gs_upper_2
#'
#' @rdname gs_upper_crt
# gs_upper_crt function [sinew] ----
gs_upper_crt <- function(theta = 0, I, a = NULL, false_pos, sides = 1,
                       tol = 0.000001, r = 18, printerr = 0) {
  # Check input arguments
  k <- as.integer(length(I))
  if (is.null(a)) {
    if (sides == 1) {
      a <- rep(-10, k)
    } else {
      a <- rep(0, k)
    }
  }

  check_scalar(sides, "integer", c(1, 2))
  if (sides == 2) {
    check_scalar(theta, "numeric", c(0, Inf))
    check_vector(a, "numeric", c(0, Inf))
  } else {
    check_scalar(theta, "numeric")
    check_vector(a, "numeric")
  }
  check_vector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  false_pos[false_pos < 1e-15] <- 1e-15
  check_vector(false_pos, "numeric", c(0, 1), c(FALSE, FALSE))
  check_scalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  check_scalar(r, "integer", c(1, 80))
  check_scalar(printerr, "integer")
  check_lengths(a, false_pos, I)

  # Set up numerical calculations
  if (false_pos[k] < 0.) stop("Final spend must be >= 0")
  r <- as.integer(r)
  printerr <- as.integer(printerr)

  storage.mode(theta) <- "double"
  storage.mode(I) <- "double"
  storage.mode(a) <- "double"
  storage.mode(false_pos) <- "double"
  storage.mode(tol) <- "double"

  prob_lo <- a
  b <- a
  retval <- as.integer(0)

  # Solve for stopping boundaries
  if (sides == 1) {
    xx <- .C("gs_upper_1",
             k, theta, I, a, b, prob_lo, false_pos, tol, r, retval, printerr)
  } else {
    xx <- .C("gs_upper_2",
             k, theta, I, a, b, prob_lo, false_pos, tol, r, retval, printerr)
  }

  y <- list(k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]],
            false_neg = xx[[6]], false_pos = xx[[7]], tol = xx[[8]], r = xx[[9]],
            error = xx[[10]])

  # Set last interim analysis boundaries equal if needed
  if (y$error == 0 && min(y$b - y$a) < 0) {
    indx <- (y$b - y$a < 0)
    y$b[indx] <- y$a[indx]
  }
  y
}

# gs_lower_crt roxy [sinew] ----
#' @title Boundary derivation for binding or non-binding futility stopping only.
#'
#' @description  \code{gs_lower_crt} is used to calculate the stopping boundaries
#'  for a group sequential trial with 1 or 2-sided futility-only stopping. Code
#'  adapted from gsDesign package.
#'
#' @param theta Effect size under null hypothesis. Default value is 0.
#' @param I Information levels at interim analyses for which stopping
#'  boundaries are calculated. At any given interim analysis this should include
#'  the information levels at the current and previous analyses.
#' @param false_neg Type II error spent for the interim analyses at the
#'  information levels specified by I.
#' @param b Upper efficacy boundary values for the interim analyses at the
#'  information levels specified by I.
#' @param sides \code{1=} 1-sided test (default) \cr \code{2=} 2-sided test
#' @param binding \code{TRUE=} binding (default) \cr \code{FALSE=} non-binding
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#'  be changed by the user.  This does not translate directly to number of
#'  digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#'  Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#'  values provide larger number of grid points and greater accuracy.  Normally
#'  \code{r} will not be changed by the user.
#' @param printerr Print output for debugging.
#'
#' @return Object containing the following elements: \item{k}{Number of interim
#'  analyses.} \item{theta}{As input.} \item{I}{As input.} \item{a}{Computed
#'  lower futility boundaries at the specified interim analyses.} \item{b}{As
#'  input.} \item{r}{As input.} \item{error}{error flag returned; 0 if
#'  convergence; 1 indicates error.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#'  Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#' @useDynLib gsDesignCRT gs_lower_1
#' @useDynLib gsDesignCRT gs_lower_2
#'
#' @rdname gs_lower_crt
# gs_lower_crt function [sinew] ----
gs_lower_crt <- function(theta = 0, I, false_neg, b = NULL, sides = 1,
                       binding = TRUE, tol = 0.000001, r = 18, printerr = 0) {
  # Check input arguments
  k <- as.integer(length(I))
  if (is.null(b)) {
    b <- rep(10, k)
  }
  check_scalar(sides, "integer", c(1, 2))
  if (sides == 2) {
    check_scalar(theta, "numeric", c(0, Inf))
    check_vector(b, "numeric", c(0, Inf))
  } else {
    check_scalar(theta, "numeric")
    check_vector(b, "numeric")
  }
  check_vector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  false_neg[false_neg < 1e-15] <- 1e-15
  check_vector(false_neg, "numeric", c(0, 1), c(FALSE, FALSE))
  check_scalar(as.numeric(binding), "integer", c(0, 1))
  check_scalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  check_scalar(r, "integer", c(1, 80))
  check_scalar(printerr, "integer")
  check_lengths(b, false_neg, I)

  # Set up numerical calculations
  if (false_neg[k] < 0.) stop("Final spend must be >= 0")
  r <- as.integer(r)
  printerr <- as.integer(printerr)

  storage.mode(theta) <- "double"
  storage.mode(I) <- "double"
  storage.mode(b) <- "double"
  storage.mode(false_neg) <- "double"
  storage.mode(tol) <- "double"

  prob_hi <- b
  a <- b
  retval <- as.integer(0)

  # Solve for stopping boundaries
  if (sides == 1) {
    xx <- .C("gs_lower_1",
             k, theta, I, a, b, false_neg, prob_hi, tol, r, retval, printerr)
  } else {
    xx <- .C("gs_lower_2",
             k, theta, I, a, b, false_neg, prob_hi, tol, r, retval, printerr)
  }

  y <- list(k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]],
            false_neg = xx[[6]], false_pos = xx[[7]], tol = xx[[8]], r = xx[[9]],
            error = xx[[10]])

  # Set last interim analysis boundaries equal if needed
  if (y$error == 0 && min(y$b - y$a) < 0) {
    indx <- (y$b - y$a < 0)
    y$b[indx] <- y$a[indx]
  }
  y
}

# gs_bounds_crt roxy [sinew] ----
#' @title Boundary derivation for efficacy and binding or non-binding futility
#'  stopping.
#'
#' @description \code{gs_bounds_crt} is used to calculate the stopping boundaries
#'  for a group sequential trial with 1 or 2-sided efficacy and binding or
#'  non-binding futility stopping. Code adapted from gsDesign package.
#'
#' @param theta Effect size under null hypothesis. Default value is 0.
#' @param I Information levels at interim analyses for which stopping
#'  boundaries are calculated. At any given interim analysis this should include
#'  the information levels at the current and previous analyses.
#' @param false_neg Type II error spent for the interim analyses at the
#'  information levels specified by I.
#' @param false_pos Type I error spent for the interim analyses at the
#'  information levels specified by I.
#' @param sides \code{1=} 1-sided test (default) \cr \code{2=} 2-sided test
#' @param binding \code{TRUE=} binding (default) \cr \code{FALSE=} non-binding
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#'  be changed by the user.  This does not translate directly to number of
#'  digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#'  Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#'  values provide larger number of grid points and greater accuracy.  Normally
#'  \code{r} will not be changed by the user.
#' @param printerr Print output for debugging.
#'
#' @return Object containing the following elements: \item{k}{Number of interim
#'  analyses.} \item{theta}{As input.} \item{I}{As input.} \item{a}{Computed
#'  lower futility boundaries at the specified interim analyses.}
#'  \item{b}{Computed upper efficacy boundaries at the specified interim
#'  analyses.} \item{r}{As input.} \item{error}{error flag returned; 0 if
#'  convergence; 1 indicates error.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#'  Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#' @useDynLib gsDesignCRT gs_bounds_1
#' @useDynLib gsDesignCRT gs_bounds_nb_1
#' @useDynLib gsDesignCRT gs_bounds_2
#' @useDynLib gsDesignCRT gs_bounds_nb_2
#'
#' @rdname gs_bounds_crt
# gs_bounds_crt function [sinew] ----
gs_bounds_crt <- function(theta = 0, I, false_neg, false_pos, sides = 1,
                        binding = TRUE, tol = 0.000001, r = 18, printerr = 0) {
  # Check input arguments
  check_scalar(theta, "numeric", c(0, Inf))
  check_vector(I, "numeric", c(0, Inf), c(FALSE, TRUE))
  false_neg[false_neg < 1e-15] <- 1e-15
  false_pos[false_pos < 1e-15] <- 1e-15
  check_vector(false_neg, "numeric", c(0, 1), c(FALSE, FALSE))
  check_vector(false_pos, "numeric", c(0, 1), c(FALSE, FALSE))
  check_scalar(sides, "integer", c(1, 2))
  check_scalar(as.numeric(binding), "integer", c(0, 1))
  check_scalar(tol, "numeric", c(0, Inf), c(FALSE, TRUE))
  check_scalar(r, "integer", c(1, 80))
  check_scalar(printerr, "integer")
  check_lengths(false_neg, false_pos, I)

  # Set up numerical calculations
  k <- as.integer(length(I))
  if (false_neg[k] < 0.) stop("Final futility spend must be >= 0")
  if (false_pos[k] < 0.) stop("Final spend must be >= 0")
  r <- as.integer(r)
  printerr <- as.integer(printerr)
  storage.mode(I) <- "double"
  storage.mode(false_neg) <- "double"
  storage.mode(false_pos) <- "double"
  storage.mode(tol) <- "double"
  a <- false_neg
  b <- false_pos
  retval <- as.integer(0)

  # Solve for stopping boundaries
  if (sides == 1) {
    if (binding) {
      xx <- .C("gs_bounds_1",
               k, theta, I, a, b, false_neg, false_pos, tol, r, retval, printerr)
    } else {
      xx <- .C("gs_bounds_nb_1",
               k, theta, I, a, b, false_neg, false_pos, tol, r, retval, printerr)
    }
  } else {
    if (binding) {
      xx <- .C("gs_bounds_2",
               k, theta, I, a, b, false_neg, false_pos, tol, r, retval, printerr)
    } else {
      xx <- .C("gs_bounds_nb_2",
               k, theta, I, a, b, false_neg, false_pos, tol, r, retval, printerr)
    }
  }

  y <- list(
    k = xx[[1]], theta = xx[[2]], I = xx[[3]], a = xx[[4]], b = xx[[5]],
    false_neg = xx[[7]], false_pos = xx[[6]], tol = xx[[8]], r = xx[[9]],
    error = xx[[10]], diff = xx[[5]] - xx[[4]]
  )
  y
}

# gs_probability_crt roxy [sinew] ----
#' @title Compute stopping boundary crossing probabilities.
#'
#' @param theta Effect size. Default value is 0.
#' @param I Information levels at interim analyses for which boundary crossing
#'  probabilities are computed.
#' @param a Lower futility boundaries for the interim analyses at the
#'  information levels specified by I.
#' @param b Upper efficacy boundaries for the interim analyses at the
#'  information levels specified by I.
#' @param sides \code{1=} 1-sided test (default) \cr \code{2=} 2-sided test.
#' @param r Integer value controlling grid for numerical integration as in
#'  Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#'  values provide larger number of grid points and greater accuracy.  Normally
#'  \code{r} will not be changed by the user.
#'
#' @return Object containing the following elements: \item{k}{Number of interim
#'  analyses.} \item{theta}{As input.} \item{I}{As input.} \item{lower}{List
#'  containing the lower futility boundaries (\code{bound}) and the
#'  corresponding probabilities of crossing the lower boundaries given
#'  \code{theta} (\code{prob}).} \item{upper}{List containing the upper efficacy
#'  boundaries (\code{bound}) and the corresponding probabilities of crossing
#'  the upper boundaries given theta (\code{prob}).} \item{power}{Estimated
#'  power for trial.} \item{futile}{Estimated probability of futility for trial}
#'  \item{r}{As input.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#'  Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#' @rdname gs_probability_crt
#' @useDynLib gsDesignCRT prob_rej_1
#' @useDynLib gsDesignCRT prob_rej_2
#'
# gs_probability_crt function [sinew] ----
gs_probability_crt <- function(theta = 0, I = 1, a = 0, b = 1,
                             sides = 1, r = 18) {
  # Check input arguments
  check_vector(theta, "numeric")

  # Check remaining input arguments
  check_vector(I - c(0, I[1:length(I) - 1]),
              "numeric", c(0, Inf), c(FALSE, FALSE))
  check_scalar(sides, "integer", c(1, 2))
  check_scalar(r, "integer", c(1, 80))
  check_lengths(I, a, b)

  # Cast integer scalars
  k <- as.integer(length(I))
  n_theta <- as.integer(length(theta))
  r <- as.integer(r)

  phi <- as.double(c(1:(k * n_theta)))
  plo <- as.double(c(1:(k * n_theta)))

  # Compute crossing probabilities
  if (sides == 1) {
    xx <- .C("prob_rej_1",
             k, n_theta, as.double(theta), as.double(I),
             as.double(a), as.double(b), plo, phi, r)
  } else if (sides == 2) {
    xx <- .C("prob_rej_2",
             k, n_theta, as.double(theta), as.double(I),
             as.double(a), as.double(b), plo, phi, r)
  } else {
    stop("invalid test type")
  }
  plo <- matrix(xx[[7]], k, n_theta)
  phi <- matrix(xx[[8]], k, n_theta)
  power <- as.vector(rep(1, k) %*% phi)
  futile <- rep(1, k) %*% plo

  x <- list(k = xx[[1]], theta = xx[[3]], I = xx[[4]],
            lower = list(bound = xx[[5]], prob = plo),
            upper = list(bound = xx[[6]], prob = phi),
            power = power, futile = futile, r = r)

  x
}


######### Hidden Functions #########

# gs_beta_diff_crt function [sinew] ----
gs_beta_diff_crt <- function(i_max, theta, beta, time, a, b, sides, r = 18) {
  # Compute difference between actual and desired Type II error
  I <- time * i_max
  x <- gs_prob_crt(theta = theta, I = I, a = a, b = b, sides = sides, r = r)
  (sum(x$prob_hi) + sum(x$prob_lo)) - (1 - beta)
}

# gs_alpha_diff_crt function [sinew] ----
gs_alpha_diff_crt <- function(i_max, theta, alpha, false_neg, time, b, sides,
                           binding, tol = 0.000001, r = 18) {
  # Compute difference between actual and desired Type I error
  I <- time * i_max
  lower <- gs_lower_crt(theta = theta, I = I, false_neg = false_neg, b = b,
                      sides = sides, binding = binding, tol = tol, r = r)
  x <- gs_prob_crt(theta = 0, I = I, a = lower$a, b = b, sides = sides, r = r)
  (sum(x$prob_hi) + sum(x$prob_lo)) - alpha
}

# gs_bound_diff_crt function [sinew] ----
gs_bound_diff_crt <- function(i_max, theta, false_neg, false_pos, time, sides,
                           binding, tol = 0.000001, r = 18) {
  # Compute difference between last upper and lower stopping boundaries
  I <- time * i_max
  bounds <- gs_bounds_crt(theta = theta, I = I, false_neg = false_neg,
                        false_pos = false_pos, sides = sides, binding = binding,
                        tol = tol, r = r)
  return(bounds$diff[length(time)])
}

# gs_prob_crt function [sinew] ----
gs_prob_crt <- function(theta, I, a, b, sides = 1, r = 18) {
  n_anal <- as.integer(length(I))
  n_theta <- as.integer(length(theta))
  phi <- as.double(c(1:(n_anal * n_theta)))
  plo <- as.double(c(1:(n_anal * n_theta)))
  if (sides == 1) {
    xx <- .C(
      "prob_rej_1", n_anal, n_theta, as.double(theta), as.double(I),
      as.double(a), as.double(b), plo, phi, as.integer(r)
    )
  } else if (sides == 2) {
    xx <- .C(
      "prob_rej_2", n_anal, n_theta, as.double(theta), as.double(I),
      as.double(a), as.double(b), plo, phi, as.integer(r)
    )
  } else {
    stop("invalid number of sides for test")
  }
  plo <- matrix(xx[[7]], n_anal, n_theta)
  phi <- matrix(xx[[8]], n_anal, n_theta)
  power <- rep(1, n_anal) %*% phi
  futile <- rep(1, n_anal) %*% plo
  list(
    k = xx[[1]], theta = xx[[3]], I = xx[[4]], a = xx[[5]], b = xx[[6]],
    prob_lo = plo, prob_hi = phi, power = power, futile = futile, r = r
  )
}