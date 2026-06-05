# sf_logistic roxy [sinew] ----
#' @title Two-parameter Spending Function Families
#' @description The functions \code{sf_logistic()}, \code{sf_normal()},
#' \code{sf_extreme_value()}, \code{sf_extreme_value2()}, \code{sf_cauchy()}, and
#' \code{sf_beta_dist()} are all 2-parameter spending function families. These
#' provide increased flexibility in some situations where the flexibility of a
#' one-parameter spending function family is not sufficient. These functions
#' all allow fitting of two points on a cumulative spending function curve; in
#' this case, four parameters are specified indicating an x and a y coordinate
#' for each of 2 points.
#'
#' \code{sf_beta_dist(alpha,t,param)} is simply \code{alpha} times the incomplete
#' beta cumulative distribution function with parameters \eqn{a} and \eqn{b}
#' passed in \code{param} evaluated at values passed in \code{t}.
#'
#' The other spending functions take the form \deqn{f(t;\alpha,a,b)=\alpha
#' F(a+b_f^{-1}(t))}{f(t;alpha,a,b)=alpha F(a+b_f^{-1}(t))} where \eqn{F()} is a
#' cumulative distribution function with values \eqn{> 0} on the real line
#' (logistic for \code{sf_logistic()}, normal for \code{sf_normal()}, extreme
#' value for \code{sf_extreme_value()} and Cauchy for \code{sf_cauchy()}) and
#' \eqn{F^{-1}()} is its inverse.
#'
#' For the logistic spending function this simplifies to
#' \deqn{f(t;\alpha,a,b)=\alpha (1-(1+e^a(t/(1-t))^b)^{-1}).}
#'
#' For the extreme value distribution with \deqn{F(x)=\exp(-\exp(-x))} this
#' simplifies to \deqn{f(t;\alpha,a,b)=\alpha \exp(-e^a (-\ln t)^b).} Since the
#' extreme value distribution is not symmetric, there is also a version where
#' the standard distribution is flipped about 0. This is reflected in
#' \code{sf_extreme_value2()} where \deqn{F(x)=1-\exp(-\exp(x)).}
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size or information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#' Values of the proportion of sample size or information for which the
#' spending function will be computed.
#' @param param In the two-parameter specification, \code{sf_beta_dist()}
#' requires 2 positive values, while \code{sf_logistic()}, \code{sf_normal()},
#' \code{sf_extreme_value()},
#'
#' \code{sf_extreme_value2()} and \code{sf_cauchy()} require the first parameter
#' to be any real value and the second to be a positive value.  The four
#' parameter specification is \code{c(t1,t2,u1,u2)} where the objective is that
#' \code{sf(t1)=alpha*u1} and \code{sf(t2)=alpha*u2}.  In this
#' parameterization, all four values must be between 0 and 1 and \code{t1 <
#' t2}, \code{u1 < u2}.
#' @return An object of type \code{spendfn}. See
#' \code{vignette("spending_function_overview")} for further details.
#'
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_crt_package_overview")}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @keywords design
#' @export
#' @rdname sf_distribution
# sf_logistic function [sinew] ----
sf_logistic <- function(alpha, t, param) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  check_vector(param, "numeric")
  len <- length(param)

  if (len == 2) {
    if (!is.numeric(param[1])) {
      stop("Numeric first logistic spending parameter not given")
    }
    if (param[2] <= 0.) {
      stop("Second logistic spending parameter param[2] must be real value > 0")
    }

    a <- param[1]
    b <- param[2]
  }
  else if (len == 4) {
    check_range(param, inclusion = c(FALSE, FALSE))
    t0 <- param[1:2]
    p0 <- param[3:4]

    if (t0[2] <= t0[1] || p0[2] <= p0[1]) {
      stop("4-parameter specification of logistic function incorrect")
    }

    xv <- log(t0 / (1 - t0))
    y <- log(p0 / (1 - p0))
    b <- (y[2] - y[1]) / (xv[2] - xv[1])
    a <- y[2] - b * xv[2]
    param <- c(a, b)
  }
  else {
    stop("Logistic spending function parameter must be of length 2 or 4")
  }

  xv <- log(t / (1 - 1 * (!is.element(t, 1)) * t))
  y <- exp(a + b * xv)
  y <- y / (1 + y)
  t[t > 1] <- 1

  x <- list(
    name = "Logistic", param = param, parname = c("a", "b"), sf = sf_logistic,
    spend = alpha * (1 * (!is.element(t, 1)) * y + 1 * is.element(t, 1)),
    bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  x
}

# sf_exponential roxy [sinew] ----
#' @title Exponential Spending Function
#' @description The function \code{sf_exponential} implements the exponential
#' spending function (Anderson and Clark, 2009).
#'
#' An exponential spending function is defined for any positive \code{nu} and
#' \eqn{0\le t\le 1} as
#' \deqn{f(t;\alpha,\nu)=\alpha(t)=\alpha^{t^{-\nu}}.}{f(t;alpha,nu)=alpha^(t^(-nu)).}
#' A value of \code{nu=0.8} approximates an O'Brien-Fleming spending function
#' well.
#'
#' The general class of spending functions this family is derived from requires
#' a continuously increasing cumulative distribution function defined for
#' \eqn{x>0} and is defined as \deqn{f(t;\alpha,
#' \nu)=1-F\left(F^{-1}(1-\alpha)/ t^\nu\right).}{% f(t; alpha,
#' nu)=1-F(F^(-1)(1-alpha)/ t^nu).} The exponential spending function can be
#' derived by letting \eqn{F(x)=1-\exp(-x)}, the exponential cumulative
#' distribution function. This function was derived as a generalization of the
#' Lan-DeMets (1983) spending function used to approximate an O'Brien-Fleming
#' spending function (\code{sf_ldof()}), \deqn{f(t; \alpha)=2-2\Phi \left(
#' \Phi^{-1}(1-\alpha/2)/ t^{1/2} \right).}{% f(t;
#' alpha)=2-2*Phi(Phi^(-1)(1-alpha/2)/t^(1/2)).}
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size/information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#' Values of the proportion of sample size/information for which the spending
#' function will be computed.
#' @param param A single positive value specifying the nu parameter for which
#' the exponential spending is to be computed; allowable range is (0, 1.5].
#' @return An object of type \code{spendfn}. See
#' \code{vignette("spending_function_overview")} for further details.
#'
#' @note The gsDesign technical manual shows how to use \code{sf_exponential()}
#'   to closely approximate an O'Brien-Fleming design.
#'   The manual is available at <https://keaven.github.io/gsd-tech-manual/>.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_crt_package_overview")}
#' @references Anderson KM and Clark JB (2009), Fitting spending functions.
#' \emph{Statistics in Medicine}; 29:321-327.
#'
#' Jennison C and Turnbull BW (2000), \emph{Group Sequential Methods with
#' Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' Lan, KKG and DeMets, DL (1983), Discrete sequential boundaries for clinical
#' trials. \emph{Biometrika}; 70:659-663.
#' @keywords design
#' @export
#' @rdname sf_exponential
# sf_exponential function [sinew] ----
sf_exponential <- function(alpha, t, param) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  # K. Wills 12/11/08: restrict param range
  # check_scalar(param, "numeric", c(0, 10), c(FALSE, TRUE))
  check_scalar(param, "numeric", c(0, 1.5), c(FALSE, TRUE))

  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1

  x <- list(
    name = "Exponential", param = param, parname = "nu", sf = sf_exponential,
    spend = alpha^(t^(-param)), bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  x
}

# sf_beta_dist roxy [sinew] ----
#' @rdname sf_distribution
#' @export
#' @importFrom stats nlminb pbeta
# sf_beta_dist function [sinew] ----
sf_beta_dist <- function(alpha, t, param) {
  if (length(param) == 4) {
    if (param[1] >= param[2]) stop("sf_beta_dist: param[1] < param[2] required in 4 point parameterization")
    if (param[3] >= param[4]) stop("sf_beta_dist: param[3] < param[4] required in 4 point parameterization")
  }

  x <- list(
    name = "Beta distribution", param = param, parname = c("a", "b"),
    sf = sf_beta_dist, spend = NULL, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  len <- length(param)

  if (len == 2) {
    check_vector(param, "numeric", c(0, Inf), c(FALSE, TRUE))
  } else if (len == 4) {
    check_vector(param, "numeric", c(0, 1), c(FALSE, FALSE))

    tem <- stats::nlminb(c(1, 1), diff_beta_dist, lower = c(0, 0),
                         xval = param[1:2], uval = param[3:4])

    if (tem$convergence != 0) {
      stop("Solution to 4-parameter specification of Beta distribution spending function not found.")
    }

    x$param <- tem$par
  } else {
    stop("Beta distribution spending function parameter must be of length 2 or 4")
  }

  t[t > 1] <- 1

  x$spend <- alpha * stats::pbeta(t, x$param[1], x$param[2])

  x
}

# sf_cauchy roxy [sinew] ----
#' @rdname sf_distribution
#' @export
#' @importFrom stats qcauchy pcauchy
# sf_cauchy function [sinew] ----
sf_cauchy <- function(alpha, t, param) {
  x <- list(
    name = "Cauchy", param = param, parname = c("a", "b"), sf = sf_cauchy,
    spend = NULL, bound = NULL, prob = NULL
  )

  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  class(x) <- "spendfn"

  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  check_vector(param, "numeric")
  len <- length(param)

  if (len == 2) {
    if (param[2] <= 0.) {
      stop("Second Cauchy spending parameter param[2] must be real value > 0")
    }

    a <- param[1]
    b <- param[2]
  } else if (len == 4) {
    t0 <- param[1:2]
    p0 <- param[3:4]

    if (t0[2] <= t0[1] || p0[2] <= p0[1]) {
      stop("4-parameter specification of Cauchy function incorrect")
    }

    xv <- stats::qcauchy(t0)
    y <- stats::qcauchy(p0)
    b <- (y[2] - y[1]) / (xv[2] - xv[1])
    a <- y[2] - b * xv[2]
    x$param <- c(a, b)
  } else {
    stop("Cauchy spending function parameter must be of length 2 or 4")
  }

  t[t > 1] <- 1
  xv <- stats::qcauchy(1 * (!is.element(t, 1)) * t)
  y <- stats::pcauchy(a + b * xv)
  x$spend <- alpha * (1 * (!is.element(t, 1)) * y + 1 * is.element(t, 1))

  x
}

# sf_extreme_value roxy [sinew] ----
#' @rdname sf_distribution
#' @export
# sf_extreme_value function [sinew] ----
sf_extreme_value <- function(alpha, t, param) {
  x <- list(
    name = "Extreme value", param = param, parname = c("a", "b"), sf = sf_extreme_value, spend = NULL,
    bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  check_vector(param, "numeric")
  len <- length(param)

  if (len == 2) {
    if (param[2] <= 0.) {
      stop("Second extreme value spending parameter param[2] must be real value > 0")
    }

    a <- param[1]
    b <- param[2]
  } else if (len == 4) {
    t0 <- param[1:2]
    p0 <- param[3:4]

    if (t0[2] <= t0[1] || p0[2] <= p0[1]) {
      stop("4-parameter specification of extreme value function incorrect")
    }

    xv <- -log(-log(t0))
    y <- -log(-log(p0))
    b <- (y[2] - y[1]) / (xv[2] - xv[1])
    a <- y[2] - b * xv[2]
    x$param <- c(a, b)
  } else {
    stop("Extreme value spending function parameter must be of length 2 or 4")
  }

  t[t > 1] <- 1
  xv <- -log(-log((!is.element(t, 1)) * t))
  y <- exp(-exp(-a - b * xv))
  x$spend <- alpha * (1 * (!is.element(t, 1)) * y + 1 * is.element(t, 1))

  x
}

# sf_extreme_value2 roxy [sinew] ----
#' @rdname sf_distribution
#' @export
# sf_extreme_value2 function [sinew] ----
sf_extreme_value2 <- function(alpha, t, param) {
  x <- list(
    name = "Extreme value 2", param = param, parname = c("a", "b"),
    sf = sf_extreme_value2, spend = NULL, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  check_vector(param, "numeric")
  len <- length(param)

  if (len == 2) {
    if (param[2] <= 0.) {
      stop("Second extreme value (2) spending parameter param[2] must be real value > 0")
    }

    a <- param[1]
    b <- param[2]
  } else if (len == 4) {
    t0 <- param[1:2]
    p0 <- param[3:4]

    if (t0[2] <= t0[1] || p0[2] <= p0[1]) {
      stop("4-parameter specification of extreme value (2) function incorrect")
    }

    xv <- log(-log(1 - t0))
    y <- log(-log(1 - p0))
    b <- (y[2] - y[1]) / (xv[2] - xv[1])
    a <- y[2] - b * xv[2]
    x$param <- c(a, b)
  } else {
    stop("Extreme value (2) spending function parameter must be of length 2 or 4")
  }

  t[t > 1] <- 1
  xv <- log(-log(1 - 1 * (!is.element(t, 1)) * t))
  y <- 1 - exp(-exp(a + b * xv))
  x$spend <- alpha * (1 * (!is.element(t, 1)) * y + 1 * is.element(t, 1))

  x
}

# sf_hsd roxy [sinew] ----
#' @title Hwang-Shih-DeCani Spending Function
#' @description The function \code{sf_hsd} implements a Hwang-Shih-DeCani
#' spending function.
#'
#' A Hwang-Shih-DeCani spending function takes the form \deqn{f(t;\alpha,
#' \gamma)=\alpha(1-e^{-\gamma t})/(1-e^{-\gamma})}{f(t; alpha, gamma) = alpha
#' * (1-exp(-gamma * t))/(1 - exp(-gamma))} where \eqn{\gamma}{gamma} is the
#' value passed in \code{param}. A value of \eqn{\gamma=-4}{gamma=-4} is used
#' to approximate an O'Brien-Fleming design (see \code{\link{sf_exponential}}
#' for a better fit), while a value of \eqn{\gamma=1}{gamma=1} approximates a
#' Pocock design well.
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size/information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#' Values of the proportion of sample size/information for which the spending
#' function will be computed.
#' @param param A single real value specifying the gamma parameter for which
#' Hwang-Shih-DeCani spending is to be computed; allowable range is [-40, 40]
#' @return An object of type \code{spendfn}. See
#' \code{vignette("spending_function_overview")} for further details.
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_crt_package_overview")}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @keywords design
#' @export
#' @rdname sf_hsd
#'
# sf_hsd function [sinew] ----
sf_hsd <- function(alpha, t, param) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_scalar(param, "numeric", c(-40, 40))

  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1

  x <- list(
    name = "Hwang-Shih-DeCani", param = param, parname = "gamma", sf = sf_hsd,
    spend = if (param == 0) t * alpha else alpha * (1. - exp(-t * param)) / (1 - exp(-param)),
    bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  x
}

# sf_ldof roxy [sinew] ----
#' @title Lan-DeMets Spending function overview
#' @description Lan and DeMets (1983) first published the method of using
#' spending functions to set boundaries for group sequential trials. In this
#' publication they proposed two specific spending functions: one to approximate
#' an O'Brien-Fleming design and the other to approximate a Pocock design.
#' The spending function to approximate O'Brien-Fleming has been generalized as
#' proposed by Liu, et al (2012)
#'
#' With \code{param=1=rho}, the Lan-DeMets (1983) spending function to
#' approximate an O'Brien-Fleming bound is implemented in the function
#' (\code{sf_ldof()}): \deqn{f(t; 
#' \alpha)=2-2\Phi\left(\Phi^{-1}(1-\alpha/2)/ t^{\rho/2}\right).}{%
#' f(t; alpha)=2-2*Phi(Phi^(-1)(1-alpha/2)/t^(rho/2)\right)}
#' For \code{rho} otherwise in \code{[.005,2]}, this is the generalized version
#' of Liu et al (2012). For \code{param} outside of \code{[.005,2]}, \code{rho}
#' is set to 1. The Lan-DeMets (1983) spending function to approximate a Pocock
#' design is implemented in the function \code{sf_ld_pocock()}:
#' \deqn{f(t;\alpha)=\alpha ln(1+(e-1)t).}{f(t;alpha)= alpha ln(1+(e-1)t).} As
#' shown in examples below, other spending functions can be used to ge t as good
#' or better approximations to Pocock and O'Brien-Fleming bounds. In particular,
#' O'Brien-Fleming bounds can be closely approximated using
#' \code{\link{sf_exponential}}.
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size/information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#' Values of the proportion of sample size/information for which the spending
#' function will be computed.
#' @param param This parameter is not used for \code{sf_ld_pocock}, not required
#' for \code{sf_ldof} and need not be specified. For \code{sf_ld_pocock} it is here
#' so that the calling sequence conforms to the standard for spending functions
#' used with \code{gs_design_crt()}. For \code{sf_ldof} it will default to 1
#' (Lan-DeMets function to approximate O'Brien-Fleming) if \code{NULL} or if
#' outside of the range \code{[.005,2]}. otherwise, it will be use to set rho
#' from Liu et al (2012).
#' @return An object of type \code{spendfn}. See
#' \code{vignette("spending_function_overview")} for further details.
#'
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_package_overview")}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' Lan, KKG and DeMets, DL (1983), Discrete sequential boundaries for clinical
#' trials. \emph{Biometrika};70: 659-663.
#'
#' Liu, Q, Lim, P, Nuamah, I, and Li, Y (2012), On adaptive error spending
#' approach for group sequential trials with random information levels.
#' \emph{Journal of biopharmaceutical statistics}; 22(4), 687-699.
#' @keywords design
#' @rdname sf_ldof
#' @export
# sf_ldof function [sinew] ----
sf_ldof <- function(alpha, t, param = NULL) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  # Following 2 lines updated 10/11/17; fix needed since default for gs_design_crt
  # is param=-4 is out of range for LDOF
  if (is.null(param) || param < .005 || param > 20) param <- 1
  check_scalar(param, "numeric", c(.005, 20), c(TRUE, TRUE))
  t[t > 1] <- 1
  if (param == 1) {
    rho <- 1
    txt <- "Lan-DeMets O'Brien-Fleming approximation"
    parname <- "none"
  } else {
    rho <- param
    txt <- "Generalized Lan-DeMets O'Brien-Fleming"
    parname <- "rho"
  }
  z <- - qnorm(alpha / 2)

  x <- list(name = txt, param = param, parname = parname, sf = sf_ldof,
            spend = 2 * (1 - pnorm(z / t^(rho / 2))), bound = NULL, prob = NULL)

  class(x) <- "spendfn"

  x
}

# sf_ld_pocock roxy [sinew] ----
#' @rdname sf_ldof
#' @export
# sf_ld_pocock function [sinew] ----
sf_ld_pocock <- function(alpha, t, param) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1

  x <- list(
    name = "Lan-DeMets Pocock approximation", param = NULL, parname = "none",
    sf = sf_ld_pocock, spend = alpha * log(1 + (exp(1) - 1) * t), bound = NULL,
    prob = NULL
  )

  class(x) <- "spendfn"

  x
}

# sf_normal roxy [sinew] ----
#' @export
#' @rdname sf_distribution
#' @importFrom stats qnorm pnorm
# sf_normal function [sinew] ----
sf_normal <- function(alpha, t, param) {
  x <- list(
    name = "Normal", param = param, parname = c("a", "b"), sf = sf_normal,
    spend = NULL, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  check_vector(param, "numeric")
  len <- length(param)

  if (len == 2) {
    if (param[2] <= 0.) {
      stop("Second Normal spending parameter param[2] must be real value > 0")
    }

    a <- param[1]
    b <- param[2]
  } else if (len == 4) {
    t0 <- param[1:2]
    p0 <- param[3:4]

    if (t0[2] <= t0[1] || p0[2] <= p0[1]) {
      stop("4-parameter specification of Normal function incorrect")
    }

    xv <- stats::qnorm(t0)
    y <- stats::qnorm(p0)
    b <- (y[2] - y[1]) / (xv[2] - xv[1])
    a <- y[2] - b * xv[2]
    x$param <- c(a, b)
  } else {
    stop("Normal spending function parameter must be of length 2 or 4")
  }

  t[t > 1] <- 1
  xv <- stats::qnorm(1 * (!is.element(t, 1)) * t)
  y <- stats::pnorm(a + b * xv)
  x$spend <- alpha * (1 * (!is.element(t, 1)) * y + 1 * is.element(t, 1))

  x
}

# sf_linear roxy [sinew] ----
#' @title  Piecewise Linear and Step Function Spending Functions
#'
#' @description The function \code{sf_linear()} allows specification of a
#' piecewise linear spending function. The function \code{sf_step()} specifies
#' a step function spending function. Both functions provide complete
#' flexibility in setting spending at desired timepoints in a group sequential
#' design. Normally these function will be passed to \code{gs_design_crt()} in the
#' parameter \code{sfu} for the upper bound or \code{sfl} for the lower bound to
#' specify a spending function family for a design. When passed to
#' \code{gs_design_crt()}, the value of \code{param} would be passed to
#' \code{sf_linear()} or \code{sf_step()} through the \code{gs_design_crt()}
#' arguments \code{sfupar} for the upper bound and \code{sflpar} for the lower
#' bound.
#'
#' Note that \code{sf_step()} allows setting a particular level of spending when
#' the timing is not strictly known; an example shows how this can inflate Type
#' I error when timing of analyses are changed based on knowing the treatment
#' effect at an interim.
#'
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size or information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#' Values of the proportion of sample size or information for which the
#' spending function will be computed.
#' @param param A vector with a positive, even length. Values must range from 0
#' to 1, inclusive. Letting \code{m <- length(param/2)}, the first \code{m}
#' points in param specify increasing values strictly between 0 and 1
#' corresponding to interim timing (proportion of final total statistical
#' information). The last \code{m} points in \code{param} specify
#' non-decreasing values from 0 to 1, inclusive, with the cumulative proportion
#' of spending at the specified timepoints.
#' @return An object of type \code{spendfn}.  The cumulative spending returned
#' in \code{sf_linear$spend} is 0 for \code{t <= 0} and \code{alpha} for
#' \code{t>=1}.  For \code{t} between specified points, linear interpolation is
#' used to determine \code{sf_linear$spend}.
#'
#' The cumulative spending returned in \code{sf_step$spend} is 0 for
#' \code{t<param[1]} and \code{alpha} for \code{t>=1}.  Letting \code{m <-
#' length(param/2)}, for \code{i=1,2,...m-1} and \code{ param[i]<= t <
#' param[i+1]}, the cumulative spending is set at \code{alpha * param[i+m]}
#' (also for \code{param[m]<=t<1}).
#'
#' Note that if \code{param[2m]} is 1, then the first time an analysis is
#' performed after the last proportion of final planned information
#' (\code{param[m]}) will be the final analysis, using any remaining error that
#' was not previously spent.
#'
#' See \code{vignette("spending_function_overview")} for further details.
#'
#' @aliases sf_linear
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_crt_package_overview")}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @keywords design
#' @rdname sf_linear
#' @export
#'
# sf_linear function [sinew] ----
sf_linear <- function(alpha, t, param) {
  x <- list(
    name = "Piecewise linear", param = param, parname = "line points",
    sf = sf_linear, spend = NULL, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1

  if (!is.numeric(param)) {
    stop("sf_linear parameter param must be numeric")
  }
  j <- length(param)
  if (floor(j / 2) * 2 != j) {
    stop("sf_linear parameter param must have even length")
  }
  k <- j / 2
  if (max(param) > 1 || min(param) < 0) {
    stop("Timepoints and cumulative proportion of spending must be >= 0 and <= 1 in sf_linear")
  }
  if (k > 1) {
    inctime <- x$param[1:k] - c(0, x$param[1:(k - 1)])
    incspend <- x$param[(k + 1):j] - c(0, x$param[(k + 1):(j - 1)])
    if ((j > 2) && (min(inctime) <= 0)) {
      stop("Timepoints must be strictly increasing in sf_linear")
    }
    if ((j > 2) && (min(incspend) < 0)) {
      stop("Spending must be non-decreasing in sf_linear")
    }
  }
  s <- t
  s[t <= 0] <- 0
  s[t >= 1] <- 1
  ind <- (0 < t) & (t <= param[1])
  s[ind] <- param[k + 1] * t[ind] / param[1]
  ind <- (1 > t) & (t >= param[k])
  s[ind] <- param[j] + (t[ind] - param[k]) / (1 - param[k]) * (1 - param[j])
  if (k > 1) {
    for (i in 2:k) {
      ind <- (param[i - 1] < t) & (t <= param[i])
      s[ind] <- param[k + i - 1] + (t[ind] - param[i - 1]) /
        (param[i] - param[i - 1]) *
        (param[k + i] - param[k + i - 1])
    }
  }
  x$spend <- alpha * s
  x
}

#' @export 
#' @rdname sf_linear
# sf_step function [sinew] ----
sf_step <- function(alpha, t, param) {
  x <- list(
    name = "Step ", param = param, parname = "line points", sf = sf_step,
    spend = NULL, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1

  if (!is.numeric(param)) {
    stop("sf_step parameter param must be numeric")
  }

  j <- length(param)
  if (floor(j / 2) * 2 != j) {
    stop("sf_step parameter param must have even length")
  }
  k <- j / 2

  if (max(param) > 1 || min(param) < 0) {
    stop("Timepoints and cumulative proportion of spending must be >= 0 and <= 1 in sf_step")
  }
  inctime <- param[1]
  if (k > 1) inctime <- c(inctime, param[2:k] - param[1:(k - 1)])
  if (min(inctime <= 0)) stop("Timepoints in param must be strictly increasing in sf_step")
  incspend <- param[k + 1]
  if (k > 1) incspend <- c(incspend, param[(k + 2):j] - param[(k + 1):(j - 1)])
  if (min(incspend) < 0) stop("Spending in param must be non-decreasing in sf_step")
  s <- rep(-3, length(t))
  s[t < param[1]] <- 0
  s[t >= param[k]] <- param[j]
  s[t >= 1] <- 1
  if (k > 1) {
    for (i in 1:(k - 1)) {
      ind <- (param[i] <= t) & (t < param[i + 1])
      s[ind] <- param[k + i]
    }
  }
  x$spend <- alpha * s
  x
}

# sf_points roxy [sinew] ----
#' @title  Pointwise Spending Function
#'
#' @description The function \code{sf_points} implements a spending function with
#' values specified for an arbitrary set of specified points. It is now
#' recommended to use sf_linear rather than sf_points. If using \code{sf_points()}
#' in a design, it is recommended to specify how to interpolate between the
#' specified points (e.g,, linear interpolation); also consider fitting smooth
#' spending functions; see \code{vignette("spending_function_overview")}.
#'
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size/information.
#' @param t A vector of points with increasing values from >0 and <=1.  Values
#' of the proportion of sample size/information for which the spending function
#' will be computed.
#' @param param A vector of the same length as \code{t} specifying the
#' cumulative proportion of spending to corresponding to each point in
#' \code{t}; must be >=0 and <=1.
#' @return An object of type \code{spendfn}. See
#' \code{vignette("spending_function_overview")} for further details.
#'
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_crt_package_overview")},
#' \link{sf_logistic}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @keywords design
#' @rdname sf_points
#' @export
# sf_points function [sinew] ----
sf_points <- function(alpha, t, param) {
  x <- list(
    name = "User-specified", param = param, parname = "Points", sf = sf_points,
    spend = NULL, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  k <- length(t)
  j <- length(param)

  if (j == k - 1) {
    x$param <- c(param, 1)
    j <- k
  }

  if (j != k) {
    stop("Cumulative user-specified proportion of spending must be specified for each interim analysis")
  }

  if (!is.numeric(param)) {
    stop("Numeric user-specified spending levels not given")
  }

  incspend <- x$param - c(0, x$param[1:k - 1])

  if (min(incspend) < 0.) {
    stop("Cumulative user-specified spending levels must be non-decreasing with each analysis")
  }

  if (max(x$param) > 1.) {
    stop("Cumulative user-specified spending must be >= 0 and <= 1")
  }

  x$spend <- alpha * x$param

  x
}

# sf_power roxy [sinew] ----
#' @title Kim-DeMets (power) Spending Function
#'
#' @description The function \code{sf_power()} implements a Kim-DeMets (power)
#' spending function. This is a flexible, one-parameter spending function
#' recommended by Jennison and Turnbull (2000).
#'
#' A Kim-DeMets spending function takes the form \deqn{f(t;\alpha,\rho)=\alpha
#' t^\rho}{f(t; alpha, rho) = alpha t^rho} where \eqn{\rho}{rho} is the value
#' passed in \code{param}. See examples below for a range of values of
#' \eqn{\rho}{rho} that may be of interest (\code{param=0.75} to \code{3} are
#' documented there).
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size/information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#' Values of the proportion of sample size/information for which the spending
#' function will be computed.
#' @param param A single, positive value specifying the \eqn{\rho}{rho}
#' parameter for which Kim-DeMets spending is to be computed; allowable range
#' is (0,50]
#' @return An object of type \code{spendfn}.
#' See \code{vignette("spending_function_overview")} for further details.
#'
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_crt_package_overview")}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @keywords design
#' @rdname sf_power
#' @export
# sf_power function [sinew] ----
sf_power <- function(alpha, t, param) {
  # K. Wills 12/11/08: restrict param range to (0, 15]
  # K. Anderson 6/11/23: restrict param range to (0, 50]
  check_scalar(param, "numeric", c(0, 50), c(FALSE, TRUE))

  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1

  x <- list(
    name = "Kim-DeMets (power)", param = param, parname = "rho", sf = sf_power,
    spend = alpha * t^param, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  x
}

# sf_t_dist roxy [sinew] ----
#' @title t-distribution Spending Function
#'
#' @description The function \code{sf_t_dist()} provides perhaps the maximum
#' flexibility among spending functions provided in the \code{gs_design_crt}
#' package. This function allows fitting of three points on a cumulative
#' spending function curve; in this case, six parameters are specified
#' indicating an x and a y coordinate for each of 3 points.
#'
#' The t-distribution spending function takes the form \deqn{f(t;\alpha)=\alpha
#' F(a+b_f^{-1}(t))} where \eqn{F()} is a cumulative t-distribution function
#' with \code{df} degrees of freedom and \eqn{F^{-1}()} is its inverse.
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size/information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#' Values of the proportion of sample size/information for which the spending
#' function will be computed.
#' @param param In the three-parameter specification, the first paramater (a)
#' may be any real value, the second (b) any positive value, and the third
#' parameter (df=degrees of freedom) any real value 1 or greater. When
#' \code{gs_design_crt()} is called with a t-distribution spending function, this
#' is the parameterization printed.  The five parameter specification is
#' \code{c(t1,t2,u1,u2,df)} where the objective is that the resulting
#' cumulative proportion of spending at \code{t} represented by \code{sf(t)}
#' satisfies \code{sf(t1)=alpha*u1}, \code{sf(t2)=alpha*u2}. The t-distribution
#' used has \code{df} degrees of freedom.  In this parameterization, all the
#' first four values must be between 0 and 1 and \code{t1 < t2}, \code{u1 <
#' u2}.  The final parameter is any real value of 1 or more. This
#' parameterization can fit any two points satisfying these requirements.  The
#' six parameter specification attempts to fit 3 points, but does not have
#' flexibility to fit any three points.  In this case, the specification for
#' \code{param} is c(t1,t2,t3,u1,u2,u3) where the objective is that
#' \code{sf(t1)=alpha*u1}, \code{sf(t2)=alpha*u2}, and \code{sf(t3)=alpha*u3}.
#' See examples to see what happens when points are specified that cannot be
#' fit.
#' @return An object of type \code{spendfn}. See
#' \code{vignette("spending_function_overview")} for further details.
#'
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_crt_package_overview")}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @keywords design
#' @rdname sf_t_dist
#' @export
#' @importFrom stats qt uniroot pt
# sf_t_dist function [sinew] ----
sf_t_dist <- function(alpha, t, param) {
  x <- list(
    name = "t-distribution", param = param, parname = c("a", "b", "df"),
    sf = sf_t_dist, spend = NULL, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  check_vector(param, "numeric")
  len <- length(param)

  if (len == 3) {
    if (param[3] < 1.) {
      stop("Final t-distribution spending parameter must be real value at least 1")
    }

    a <- param[1]
    b <- param[2]
    df <- param[3]
  } else if (len == 5) {
    t0 <- param[1:2]
    p0 <- param[3:4]
    df <- param[5]

    if (param[5] < 1.) {
      stop("Final t-distribution spending parameter must be real value at least 1")
    }

    if (t0[2] <= t0[1] || p0[2] <= p0[1]) {
      stop("5-parameter specification of t-distribution spending function incorrect")
    }

    xv <- stats::qt(t0, df)
    y <- stats::qt(p0, df)
    b <- (y[2] - y[1]) / (xv[2] - xv[1])
    a <- y[2] - b * xv[2]
  } else if (len == 6) {
    t0 <- param[1:3]
    p0 <- param[4:6]

    if (t0[2] <= t0[1] || p0[2] <= p0[1]) {
      stop("6-parameter specification of t-distribution spending function incorrect")
    }

    # check Cauchy and normal which must err in opposite directions for a 
    # solution to exist
    unorm <- sf_normal(alpha, t0[3],
                      param = c(t0[1:2], p0[1:2]))$spend / alpha - p0[3]
    ucauchy <- sf_cauchy(alpha, t0[3],
                        param = c(t0[1:2], p0[1:2]))$spend / alpha - p0[3]

    if (unorm * ucauchy >= 0.) {
      stop("6-parameter specification of t-distribution spending function did not produce a solution")
    }

    sol <- stats::uniroot(t_dist_diff, interval = c(1, 200), t0 = t0, p0 = p0)
    df <- sol$root
    xv <- stats::qt(t0, df)
    y <- stats::qt(p0, df)
    b <- (y[2] - y[1]) / (xv[2] - xv[1])
    a <- y[2] - b * xv[2]
  } else {
    stop("t-distribution spending function parameter must be of length 3, 5 or 6")
  }

  x$param <- c(a, b, df)
  t[t > 1] <- 1
  xv <- stats::qt(1 * (!is.element(t, 1)) * t, df)
  y <- stats::pt(a + b * xv, df)
  x$spend <- alpha * (1 * (!is.element(t, 1)) * y + 1 * is.element(t, 1))

  x
}

# sf_truncated roxy [sinew] ----
#' @title Truncated, trimmed and gapped spending functions
#' @description The functions \code{sf_truncated()} and \code{sf_trimmed} apply
#' any other spending function over a restricted range. This allows eliminating
#' spending for early interim analyses when you desire not to stop for the bound
#' being specified; this is usually applied to eliminate early tests for a
#' positive efficacy finding. The truncation can come late in the trial if you
#' desire to stop a trial any time after, say, 90 percent of information is
#' available and an analysis is performed. This allows full Type I error
#' spending if the final analysis occurs early. Both functions set cumulative
#' spending to 0 below a 'spending interval' in the interval [0,1], and set
#' cumulative spending to 1 above this range. \code{sf_trimmed()} otherwise does
#' not change an input spending function that is specified; probably the
#' preferred and more intuitive method in most cases. \code{sf_truncated()}
#' resets the time scale on which the input spending function is computed to the
#' 'spending interval.'
#'
#' \code{sf_gapped()} allows elimination of analyses after some time point in
#' the trial; see details and examples.
#'
#' \code{sf_trimmed} simply computes the value of the input spending function
#' and parameters in the sub-range of [0,1], sets spending to 0 below this
#' range and sets spending to 1 above this range.
#'
#' \code{sf_gapped} spends outside of the range provided in trange. Below
#' trange, the input spending function is used. Above trange, full spending is
#' used; i.e., the first analysis performed above the interval in trange is the
#' final analysis. As long as the input spending function is strictly
#' increasing, this means that the first interim in the interval trange is the
#' final interim analysis for the bound being specified.
#'
#' \code{sf_truncated} compresses spending into a sub-range of [0,1]. The
#' parameter \code{param$trange} specifies the range over which spending is to
#' occur. Within this range, spending is spent according to the spending
#' function specified in \code{param$sf} along with the corresponding spending
#' function parameter(s) in \code{param$param}. See example using
#' \code{sf_linear} that spends uniformly over specified range.
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Normally,
#' \code{alpha=0.025} for one-sided Type I error specification or
#' \code{alpha=0.1} for Type II error specification. However, this could be set
#' to 1 if for descriptive purposes you wish to see the proportion of spending
#' as a function of the proportion of sample size or information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#' Values of the proportion of sample size or information for which the
#' spending function will be computed.
#' @param param a list containing the elements sf (a spendfn object such as
#' sf_hsd), trange (the range over which the spending function increases from 0
#' to 1; 0 <= trange[1]<trange[2] <=1; for sf_gapped, trange[1] must be > 0),
#' and param (null for a spending function with no parameters or a scalar or
#' vector of parameters needed to fully specify the spending function in sf).
#' @return An object of type \code{spendfn}. See
#' \code{vignette("spending_function_overview")} for further details.
#'
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#' @seealso \code{vignette("spending_function_overview")},
#' \code{\link{gs_design_crt}}, \code{vignette("gs_design_crt_package_overview")}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @keywords design
#' @export
#' @rdname sf_special
# sf_truncated function [sinew] ----
sf_truncated <- function(alpha, t, param) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  if (!is.list(param)) stop("param must be a list. See help(sf_truncated)")
  if (!max(names(param) == "trange")) stop("param must include trange, sf, param. See help(sf_truncated)")
  if (!max(names(param) == "sf")) stop("param must include trange, sf, param. See help(sf_truncated)")
  if (!max(names(param) == "param")) stop("param must include trange, sf, param. See help(sf_truncated)")
  if (!is.vector(param$trange)) stop("param$trange must be a vector of length 2 with 0 <= param$trange[1] <param$trange[2]<=1. See help(sf_truncated)")
  if (length(param$trange) != 2) stop("param$trange parameter must be a vector of length 2 with 0 <= param$trange[1] <param$trange[2]<=1. See help(sf_truncated)")
  if (param$trange[1] >= 1. | param$trange[2] <= param$trange[1] | param$trange[2] <= 0) {
    stop("param$trange must be a vector of length 2 with 0 <= param$trange[1] < param$trange[2]<=1. See help(sf_truncated)")
  }
  if (!inherits(param$sf, "function")) stop("param$sf must be a spending function")
  if (!is.numeric(param$param)) stop("param$param must be numeric")
  spend <- as.vector(rep(0, length(t)))
  spend[t >= param$trange[2]] <- alpha
  indx <- param$trange[1] < t & t < param$trange[2]
  if (max(indx)) {
    s <- param$sf(alpha = alpha, t = (t[indx] - param$trange[1]) / (param$trange[2] - param$trange[1]), param$param)
    spend[indx] <- s$spend
  }
  # the following line is awkward, but necessary to get the input spending
  # function name in some cases
  param2 <- param$sf(alpha = alpha, t = .5, param = param$param)
  param$name <- param2$name
  param$parname <- param2$parname
  x <- list(
    name = "Truncated", param = param, parname = "range",
    sf = sf_truncated, spend = spend, bound = NULL, prob = NULL
  )
  class(x) <- "spendfn"
  x
}

# sf_trimmed roxy [sinew] ----
#' @export
#' @rdname sf_special
# sf_trimmed function [sinew] ----
sf_trimmed <- function(alpha, t, param) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  if (!is.list(param)) stop("param must be a list. See help(sf_trimmed)")
  if (!max(names(param) == "trange")) stop("param must include trange, sf, param. See help(sf_trimmed)")
  if (!max(names(param) == "sf")) stop("param must include trange, sf, param. See help(sf_trimmed)")
  if (!max(names(param) == "param")) stop("param must include trange, sf, param. See help(sf_trimmed)")
  if (!is.vector(param$trange)) stop("param$trange must be a vector of length 2 with 0 <= param$trange[1] <param$trange[2]<=1. See help(sf_trimmed)")
  if (length(param$trange) != 2) stop("param$trange parameter must be a vector of length 2 with 0 <= param$trange[1] <param$trange[2]<=1. See help(sf_trimmed)")
  if (param$trange[1] >= 1. | param$trange[2] <= param$trange[1] | param$trange[2] <= 0) {
    stop("param$trange must be a vector of length 2 with 0 <= param$trange[1] < param$trange[2]<=1. See help(sf_trimmed)")
  }
  if (!inherits(param$sf, "function")) stop("param$sf must be a spending function")
  if (!is.numeric(param$param)) stop("param$param must be numeric")
  spend <- as.vector(rep(0, length(t)))
  spend[t >= param$trange[2]] <- alpha
  indx <- param$trange[1] < t & t < param$trange[2]
  if (max(indx)) {
    s <- param$sf(alpha = alpha, t = t[indx], param$param)
    spend[indx] <- s$spend
  }
  # the following line is awkward, but necessary to get the input spending
  # function name in some cases
  param2 <- param$sf(alpha = alpha, t = .5, param = param$param)
  param$name <- param2$name
  param$parname <- param2$parname
  x <- list(
    name = "Trimmed", param = param, parname = "range",
    sf = sf_trimmed, spend = spend, bound = NULL, prob = NULL
  )
  class(x) <- "spendfn"
  x
}

# sf_gapped roxy [sinew] ----
#' @export
#' @rdname sf_special
# sf_gapped function [sinew] ----
sf_gapped <- function(alpha, t, param) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  if (!is.list(param)) stop("param must be a list. See help(sf_trimmed)")
  if (!max(names(param) == "trange")) stop("param must include trange, sf, param. See help(sf_gapped)")
  if (!max(names(param) == "sf")) stop("param must include trange, sf, param. See help(sf_gapped)")
  if (!max(names(param) == "param")) stop("param must include trange, sf, param. See help(sf_gapped)")
  if (!is.vector(param$trange)) stop("param$trange must be a vector of length 2 with 0 < param$trange[1] < param$trange[2]<=1. See help(sf_gapped)")
  if (length(param$trange) != 2) stop("param$trange parameter must be a vector of length 2 with 0 < param$trange[1] <param$trange[2]<=1. See help(sf_gapped)")
  if (param$trange[1] >= 1. | param$trange[2] <= param$trange[1] | param$trange[2] <= 0 | param$trange[1] <= 0) {
    stop("param$trange must be a vector of length 2 with 0 < param$trange[1] < param$trange[2]<=1. See help(sf_trimmed)")
  }
  if (!inherits(param$sf, "function")) stop("param$sf must be a spending function")
  if (!is.numeric(param$param)) stop("param$param must be numeric")
  spend <- as.vector(rep(0, length(t)))
  spend[t >= param$trange[2]] <- alpha
  indx <- param$trange[1] > t
  if (max(indx)) {
    s <- param$sf(alpha = alpha, t = t[indx], param$param)
    spend[indx] <- s$spend
  }
  indx <- (param$trange[1] <= t & param$trange[2] > t)
  if (max(indx)) {
    spend[indx] <- param$sf(alpha = alpha, t = param$trange[1], param$param)$spend
  }
  # the following line is awkward, but necessary to get the input spending
  # function name in some cases
  param2 <- param$sf(alpha = alpha, t = .5, param = param$param)
  param$name <- param2$name
  param$parname <- param2$parname
  x <- list(
    name = "Gapped", param = param, parname = "range",
    sf = sf_gapped, spend = spend, bound = NULL, prob = NULL
  )
  class(x) <- "spendfn"
  x
}

# spending_function roxy [sinew] ----
#' @title Spending Function
#'
#' @param alpha Real value \eqn{> 0} and no more than 1. Defaults in calls to
#'   \code{gs_design_crt()} are \code{alpha=0.025} for one-sided Type I error
#'   specification and \code{alpha=0.1} for Type II error specification.
#'   However, this could be set to 1 if, for descriptive purposes, you wish to
#'   see the proportion of spending as a function of the proportion of sample
#'   size/information.
#' @param t A vector of points with increasing values from 0 to 1, inclusive.
#'   Values of the proportion of sample size/information for which the spending
#'   function will be computed.
#' @param param A single real value or a vector of real values specifying the
#'   spending function parameter(s); this must be appropriately matched to the
#'   spending function specified.
#'
#' @return
#' \code{spending_function} and spending functions in general produce an
#' object of type \code{spendfn}.
#' \item{name}{A character string with the name of the spending function.}
#' \item{param}{any parameters used for the spending function.}
#' \item{parname}{a character string or strings with the name(s) of
#' the parameter(s) in \code{param}.}
#' \item{sf}{the spending function specified.}
#' \item{spend}{a vector of cumulative spending values corresponding to
#' the input values in \code{t}.}
#' \item{bound}{this is null when returned from the spending function,
#' but is set in \code{gs_design_crt()} if the spending function is called
#' from there.  Contains z-values for bounds of a design.}
#' \item{prob}{this is null when returned from the spending function,
#' but is set in \code{gs_design_crt()} if the spending function is called
#' from there.  Contains probabilities of boundary crossing at \code{i}-th
#' analysis for \code{j}-th theta value input to \code{gs_design_crt()} in
#' \code{prob[i,j]}.}
#'
#' @export
#'
#' @rdname spending_function
#'
#' @note The gsDesign technical manual is available at
#'   \url{https://keaven.github.io/gsd-tech-manual/}.
#'
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#'
#' @seealso \code{\link{gs_design_crt}}, \code{\link{sf_hsd}},
#' \code{\link{sf_power}}, \code{\link{sf_logistic}}, \code{\link{sf_exponential}},
#' \code{\link{sf_truncated}}, \code{vignette("gs_design_crt_package_overview")}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @keywords design
#'
# spending_function function [sinew] ----
spending_function <- function(alpha, t, param) {
  check_scalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1

  x <- list(
    name = "Linear", param = param, parname = "none", sf = spending_function,
    spend = alpha * t, bound = NULL, prob = NULL
  )

  class(x) <- "spendfn"

  x
}

# diff_beta_dist roxy [sinew] ----
#' @importFrom stats pbeta
# diff_beta_dist function [sinew] ----
diff_beta_dist <- function(aval, xval, uval) {
  if (min(aval) <= 0.) {
    return(1000)
  }

  diff <- uval - stats::pbeta(xval, aval[1], aval[2])

  sum(diff^2)
}

# t_dist_diff roxy [sinew] ----
#' @importFrom stats qt
# t_dist_diff function [sinew] ----
t_dist_diff <- function(x, t0, p0) {
  xv <- stats::qt(t0, x)
  y <- stats::qt(p0, x)
  b <- (y[2] - y[1]) / (xv[2] - xv[1])
  a <- y[2] - b * xv[2]

  a + b * xv[3] - y[3]
}
