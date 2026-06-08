# gs_design_crt roxy [sinew] ----
#' @title Compute stopping boundaries, maximum sample size, and expected sample
#'  sizes for a group sequential cluster randomized trial.
#'
#' @description  \code{gs_design_crt()} is used to determine the maximum sample
#'  size needed for a specified parallel group sequential cluster randomized
#'  trial to detect a clinically meaningful effect size with some Type I error
#'  rate and power. Code adapted from gsDesign package.
#'
#' @param k Number of analyses planned, including interim and final.
#' @param outcome_type \code{1=}continuous difference of means
#' \cr \code{2=}binary difference of proportions
#' @param test_type \code{1=} early stopping for efficacy only
#' \cr \code{2=} early stopping for binding futility only
#' \cr \code{3=} early stopping for non-binding futility only
#' \cr \code{4=} early stopping for either efficacy or binding futility
#' \cr \code{5=} early stopping for either efficacy or non-binding futility
#' @param test_sides \code{1=} one-sided test \cr \code{2=} two-sided test
#' @param size_type \code{1=}clusters per arm \cr \code{2=}cluster size
#' @param timing_type \code{1=} maximum and expected sample sizes based on
#'  specified information fractions in \code{info_timing}; recruit according to
#'  design specified in \code{recruit_type}
#'  \cr \code{2=} maximum sample size based on specified information fractions
#'  in \code{info_timing} and expected sample sizes based on specified sample
#'  size fractions in \code{m_timing} and \code{n_timing}
#'  \cr \code{3=} maximum and expected sample sizes based on specified
#'  sample size fractions in \code{m_timing} and \code{n_timing}
#'  @param recruit_type \code{1=}recruit clusters with fixed sizes
#'  \cr \code{2=}recruit individuals into fixed number of clusters
#'  \cr \code{3=}recruit at both cluster and individual levels
#' @param delta Effect size for theta under alternative hypothesis. Must be > 0.
#' @param sigma_vec Standard deviations for control and treatment groups
#'  (continuous case).
#' @param p_vec Probabilities of event for control and treatment groups
#'  (binary case).
#' @param rho Intraclass correlation coefficient. Default value is 0.
#' @param alpha Type I error, default value is 0.05.
#' @param beta Type II error, default value is 0.1 (90\% power).
#' @param m Number of clusters for finding maximum mean cluster size. If
#'  \code{m} is a scalar, it is treated as the total number of clusters across
#'  arms. If \code{m} is a vector of length 2, it is treated as the number of
#'  clusters per arm.
#' @param m_alloc Allocation ratio of clusters per arm. Default is
#'  \code{c(0.5, 0.5)}.
#' @param n Mean cluster size for finding maximum number of clusters. If
#'  \code{n} is a scalar, it is treated as the average cluster size for both
#'  arms. If \code{n} is a vector of length 2, it is treated as the average
#'  cluster size per arm.
#' @param n_cv Coefficient of variation for cluster size. If \code{n_cv} is a
#'  scalar, it is treated as the coefficient of variation for both arms. If
#'  \code{n_cv} is a vector of length 2, it is treated as the coefficient of
#'  variation per arm.
#' @param info_timing Sets timing of interim analyses based on information
#'  fractions. Default of 1 produces analyses at equal-spaced increments.
#'  Otherwise, this is a vector of length \code{k} or \code{k-1}. The values
#'  should satisfy \code{0 < info_timing[1] < info_timing[2] < ... <
#'  info_timing[k-1] < info_timing[k]=1}.
#' @param m_timing Sets timing of interim analyses based on fractions of the
#'  number of clusters.
#' @param n_timing Sets timing of interim analyses based on cluster size
#' @param alpha_sf A spending function or a character string indicating an
#'  upper boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#'  \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#'  bounds). The default value is \code{sf_ldof} which is a Lan-DeMets
#'  O'Brien-Fleming spending function. See details,
#'  \code{vignette("spending_function_overview")}, manual and examples.
#' @param alpha_sf_par Real value, default is \eqn{-4} which is an
#'  O'Brien-Fleming-like conservative bound when used with a
#'  Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#'  functions. The parameter \code{alpha_sf_par} specifies any parameters needed
#'  for the spending function specified by \code{alpha_sf}; this will be ignored
#'  for spending functions (\code{sf_ldof}, \code{sf_ld_pocock}) or bound types
#'  (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param beta_sf A spending function or a character string indicating an
#'  lower boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#'  \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#'  bounds). The default value is \code{sf_ldof} which is a Lan-DeMets
#'  O'Brien-Fleming spending function. See details,
#'  \code{vignette("spending_function_overview")}, manual and examples.
#' @param beta_sf_par Real value, default is \eqn{-4} which is an
#'  O'Brien-Fleming-like conservative bound when used with a
#'  Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#'  functions. The parameter \code{beta_sf_par} specifies any parameters needed
#'  for the spending function specified by \code{beta_sf}; this will be ignored
#'  for spending functions (\code{sf_ldof}, \code{sf_ld_pocock}) or bound types
#'  (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#'  be changed by the user.  This does not translate directly to number of
#'  digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#'  Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#'  values provide larger number of grid points and greater accuracy.  Normally
#'  \code{r} will not be changed by the user.
#'
#' @return Object containing the following elements: \item{k}{As
#'  input.} \item{outcome_type}{As input.} \item{test_type}{As input.}
#'  \item{test_sides}{As input.} \item{size_type}{As input.}
#'  \item{recruit_type}{As input.} \item{timing_type}{As input.} \item{delta}{As
#'  input.} \item{sigma_vec}{As input.} \item{p_vec}{As input.} \item{rho}{As
#'  input.} \item{alpha}{As input.} \item{beta}{As input.} \item{info_timing}{As
#'  input.} \item{m_timing}{As input.} \item{n_timing}{As input.}
#'  \item{info_schedule}{Fraction of maximum information at each planned interim
#'  analysis based on \code{timing_type} and \code{info_timing}.}
#'  \item{m_schedule}{Fraction of maximum number of clusters per arm at each
#'  planned interim analysis based on \code{timing_type} and \code{m_timing}.}
#'  \item{n_schedule}{Fraction of maximum cluster size at each planned interim
#'  analysis based on \code{timing_type} and \code{n_timing}.} \item{i}{Fisher
#'  information at each planned interim analysis based on \code{timing_type}.}
#'  \item{max_i}{Maximum information corresponding to design specifications.}
#'  \item{m}{Number of clusters per arm at each planned interim analysis.}
#'  \item{max_m}{Maximum number of clusters per arm.} \item{e_m}{A vector of
#'  length 2 with expected number of clusters per arm under the null and
#'  alternative hypotheses. For simplicity, the expected sizes with non-binding
#'  futility boundaries are calculated assuming the boundaries are binding
#'  futility.} \item{n}{Mean cluster size at each planned interim analysis.}
#'  \item{max_n}{Maximum mean cluster size.} \item{e_n}{A vector of length 2
#'  with expected mean cluster sizes under the null and alternative hypotheses.
#'  For simplicity, the expected sizes with non-binding futility boundaries are
#'  calculated assuming the futility boundaries are binding.}
#'  \item{max_total}{Maximum number of individuals in the trial.}
#'  \item{e_total}{A vector of length 2 with expected total number of
#'  individuals in the trial under the null and alternative hypotheses. For
#'  simplicity, the expected sizes with non-binding futility boundaries are
#'  calculated assuming the futility boundaries are binding.} \item{sufficient}{
#'  Value denoting whether calculated sample size will be sufficient to achieve
#'  specified Type I error rate and power given the trial specifications.}
#'  \item{lower_bound}{Calculated lower futility boundaries under analysis
#'  schedule specified by \code{timing_type}} \item{upper_bound}{Calculated
#'  upper efficacy boundaries under analysis schedule specified by
#'  \code{timing_types}.} \item{tol}{As input.} \item{r}{As input.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#'  Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @references Anderson K (2025). gsDesign: Group Sequential Design.
#'  R package version 3.8.0, https://keaven.github.io/gsDesign/.
#'
#' @export
#'
#' @name gs_design_crt
# gs_design_crt function [sinew] ----
gs_design_crt <- function(k = 3, outcome_type = 1, test_type = 1, test_sides = 1,
                        size_type = 1, timing_type = 2, recruit_type = 1,
                        delta = 1, sigma_vec = c(1, 1), p_vec = c(0.5, 0.5),
                        rho = 0, alpha = 0.05, beta = 0.1,
                        m = 2, m_alloc = c(0.5, 0.5), n = 1, n_cv = c(0, 0),
                        info_timing = 1, m_timing = c(1, 1), n_timing = c(1, 1),
                        alpha_sf = sf_ldof, alpha_sf_par = -4,
                        beta_sf = sf_ldof, beta_sf_par = -4,
                        tol = 1e-6, r = 18) {

  ### Check and pre-process inputs ###

  # Design parameters
  check_scalar(k, "integer", c(1, Inf))
  check_scalar(outcome_type, "integer", c(1, 2))
  check_scalar(test_type, "integer", c(1, 5))
  check_scalar(test_sides, "integer", c(1, 2))
  check_scalar(size_type, "integer", c(1, 2))
  check_scalar(recruit_type, "integer", c(1, 3))
  check_scalar(timing_type, "integer", c(1, 3))

  # Outcome parameters
  if (outcome_type == 1) {
    check_scalar(delta, "numeric", c(0, Inf), c(FALSE, FALSE))
    check_vector(sigma_vec, "numeric", c(0, Inf), c(FALSE, FALSE), length = 2)
    var_vec <- sigma_vec^2
  } else if (outcome_type == 2) {
    check_scalar(delta, "numeric", c(0, 1), c(FALSE, TRUE))
    check_vector(p_vec, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)
    var_vec <- p_vec * (1 - p_vec)
  }

  # Type I error, power, and ICC
  check_scalar(alpha, "numeric", 0:1, c(FALSE, FALSE))
  check_scalar(beta, "numeric", c(0, 1 - alpha), c(FALSE, FALSE))
  rho <- expand_two(rho)
  check_vector(rho, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)

  # Number of clusters
  check_alloc(m_alloc)
  if (length(m) == 1) { # Input is total clusters
    m_total_fix <- m
    m_fix <- round(m * m_alloc)
  } else { # Input is clusters per arm
    m_fix <- m
    m_total_fix <- sum(m)
    if (size_type == 2) {
      m_alloc <- m_fix / sum(m_fix)
    }
  }
  check_vector(m_fix, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_scalar(m_total_fix, "integer", c(0, Inf), c(FALSE, FALSE))

  # Average cluster size
  n_fix  <- expand_two(n)
  check_vector(n_fix, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  n_cv   <- expand_two(n_cv)
  check_vector(n_cv, "numeric", c(0, Inf), c(TRUE, FALSE), length = 2)

  # Error spending
  check_sf(alpha_sf)
  check_sf(beta_sf)

  # Numerical parameters
  check_scalar(tol, "numeric", c(0, 0.1), c(FALSE, TRUE))
  check_scalar(r, "integer", c(1, 80))
  if (length(info_timing) > 0) {
    check_vector(info_timing, "numeric", c(0, 1), c(FALSE, TRUE))
  }

  ### Initialize design object ###
  x <- list(k = k, outcome_type = outcome_type, test_type = test_type,
            test_sides = test_sides, size_type = size_type,
            recruit_type = recruit_type, timing_type = timing_type,
            delta = delta, sigma_vec = sigma_vec, p_vec = p_vec,
            var_vec = var_vec, rho = rho, alpha = alpha, beta = beta,
            info_timing = info_timing, m_timing = m_timing, n_timing = n_timing,
            info_schedule = NULL, m_schedule = NULL, n_schedule = NULL,
            i = 1, max_i = 1, m = m_fix, m_alloc = m_alloc, max_m = m_fix,
            m_total = m_total_fix, max_m_total = m_total_fix, e_m = m_fix,
            n = n_fix, n_cv = n_cv, max_n = n_fix, e_n = n_fix,
            max_total = m_fix * n_fix, e_total = m_fix * n_fix, sufficient = 1,
            alpha_sf = alpha_sf, alpha_sf_par = alpha_sf_par,
            beta_sf = beta_sf, beta_sf_par = beta_sf_par,
            lower_bound = rep(0, k), upper_bound = rep(0, k), tol = tol, r = r)

  ### Compute maximum sample size based on information fractions ###
  if (x$timing_type <= 2) {
    # Compute maximum information
    x_max <- gs_max_info_crt(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par)
    x$i <- x_max$i
    x$max_i <- x_max$i[x$k]

    # Convert maximum information to sample sizes
    if (size_type == 1) { # Clusters per arm
      x$max_m <- m_diff(i = x$max_i,
                       m_alloc = x$m_alloc,
                       n_fix = x$max_n,
                       n_cv = x$n_cv,
                       var_vec = x$var_vec,
                       rho = x$rho) / 0.9
      x$max_m_total <- sum(x$max_m)
    } else if (size_type == 2) { # Mean cluster size
      x$sufficient <- check_suff(m_total = x$max_m_total,
                                m_alloc = x$m_alloc,
                                n_cv = x$n_cv,
                                max_i = x$max_i,
                                var_vec = x$var_vec,
                                rho = x$rho,
                                alpha = x$alpha,
                                beta = x$beta)

      if (x$sufficient == 0) {
        warning("trial not feasible with given specifications")
        x$max_n <- x$n_fix
      } else {
        x$max_n <- rep(n_diff(i = x$max_i,
                             m_fix = x$max_m,
                             n_cv = x$n_cv,
                             var_vec = x$var_vec,
                             rho = x$rho), 2)
      }
    } else {
      stop("invalid sample size variable")
    }
  } else {
    if (size_type == 1) { # Clusters per arm
      # Compute maximum information based on specified sample size fractions
      x_max <- gs_max_size_crt(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par)
      x$max_m <- x_max$max_m
      x$i <- x_max$info
      x$max_i <- x_max$info[x$k]
    } else if (size_type == 2) { # Mean cluster size
      # Compute maximum information based on specified sample size fractions
      x_max <- gs_max_size_crt(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par)
      x$sufficient <- x_max$sufficient

      if (x$sufficient == 0) {
        warning("trial not feasible with given specifications")
      } else {
        x$max_n <- rep(x_max$max_n, 2)
        x$i <- x_max$info
        x$max_i <- x_max$info[x$k]
      }
    } else {
      stop("invalid sample size variable")
    }
  }

  ### Compute information levels, expected sample size, and boundaries according
  # to specified analysis schedule if feasibility criteria is met ###
  if (x$sufficient == 1) {
    x_expect <- gs_expect_crt(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par)

    # Schedule of analyses
    x$info_schedule <- x_expect$info_schedule
    x$m_schedule <- x_expect$m_schedule
    x$n_schedule <- x_expect$n_schedule

    # Information and sample sizes at scheduled analyses
    x$i <- x_expect$info
    x$m <- x_expect$m
    x$n <- x_expect$n

    # Expected sample sizes at scheduled analyses
    x$e_m <- x_expect$e_m
    x$e_n <- x_expect$e_n
    x$e_total <- x_expect$e_total

    # Compute total number of participants in trial
    x$total <- x$m * x$n
    x$max_total <- x$max_m * x$max_n

    # Check spacing of analyses
    x$i_spacing <- x$i[2:x$k] - x$i[1:(x$k - 1)]
    if (any(x$i_spacing < 1)) {
      x$sufficient <- 0
    }

    # Calculate stopping boundaries and stopping probabilities
    x$lower_bound <- x_expect$lower_bound
    x$upper_bound <- x_expect$upper_bound
    x$power <- x_expect$power
    x$futile <- x_expect$futile
  }
  class(x) <- c("gs_design_crt", class(x))
  return(x)
}

# gs_print_crt function [sinew] ----
#' @title Print group sequential CRT design
#'
#' @param x Object of class \code{gs_design_crt} created by \code{gs_design_crt()}.
#'
#' @return Printed output of design parameters and table of stopping boundaries
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @export
#' @name gs_print_crt
# gs_print_crt function [sinew] ----
gs_print_crt <- function(x) {
  cat("Group sequential CRT design\n")
  cat("Outcome: ", if (x$outcome_type == 1) "continuous" else "binary", "\n")
  cat("Total Number of Analyses:", x$k, "\n")
  cat("Max clusters per arm:", x$max_m, "\n")
  cat("Expected clusters per arm under H0:", x$e_m[, 1], "\n")
  cat("Expected clusters per arm under H1:", x$e_m[, 2], "\n")
  cat("Max mean cluster size:", x$max_n, "\n")
  cat("Expected mean cluster size under H0:", x$e_n[, 1], "\n")
  cat("Expected mean cluster size under H1:", x$e_n[, 2], "\n")
  cat("Max total sample size:", x$max_total, "\n")
  cat("Expected total sample size under H0 and H1:", x$e_total, "\n")
  cat("Overall Alpha:", x$alpha, " Overall Beta:", x$beta, "\n\n")
  # Boundary table
  boundary_table <- data.frame(
    Analysis = 1:x$k,
    Information = x$i,
    Clusters_Arm0 = x$m[1, ],
    Clusters_Arm1 = x$m[2, ],
    Mean_Cluster_Size_Arm0 = x$n[1, ],
    Mean_Cluster_Size_Arm1 = x$n[2, ],
    Lower_Bound = x$lower_bound,
    Upper_Bound = x$upper_bound
  )
  print(boundary_table)
}

# gs_plot_crt function [sinew] ----
#' @title Plot group sequential CRT design
#'
#' @param x Object of class \code{gs_design_crt} created by \code{gs_design_crt()}.
#'
#' @return Plot of stopping boundaries for the group sequential CRT design
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @export
#' @name gs_plot_crt
# gs_plot_crt function [sinew] ----
gs_plot_crt <- function(x) {
  # Plot lower and upper boundaries with ggplot
  boundary_data <- data.frame(
    analysis = 1:x$k,
    lower_bound = x$lower_bound,
    upper_bound = x$upper_bound
  )
  boundary_data_melted <- melt(boundary_data, id.vars = "analysis",
                               variable.name = "boundary_type",
                               value.name = "boundary_value")
  ggplot(boundary_data_melted, aes(x = analysis,
                                   y = boundary_value,
                                   color = boundary_type)) +
    geom_line() +
    geom_point() +
    theme_bw() + geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Stopping Boundaries for Group Sequential CRT Design",
         x = "k",
         y = "Boundary Value") +
    scale_color_manual(values = c("lower_bound" = "red",
                                  "upper_bound" = "blue"),
                       labels = c("Futility Bound", "Efficacy Bound"),
                       name = "Boundary Type")
}


######### Main hidden functions for gs_design_crt #########

# gs_max_info_crt function [sinew] ----
gs_max_info_crt <- function(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par) {
  # Timing of clusters per arm - 2 x k matrix
  x$info_schedule <- info_set_schedule(x$info_timing, x$k)

  # Set initial information guess
  i0 <- set_i0(x)

  # Partition error probabilities based on error spending
  part <- partition_error(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par)
  false_pos <- part$false_pos
  false_neg <- part$false_neg
  binding <- part$binding

  # Determine maximum information that achieves a_K = b_K with error spending
  x$i <- stats::uniroot(gs_bound_diff_crt, lower = i0, upper = 10 * i0,
                        extendInt = "yes", theta = x$delta,
                        false_neg = false_neg, false_pos = false_pos,
                        time = x$info_schedule, sides = x$test_sides,
                        binding = binding, tol = x$tol,
                        r = x$r)$root * x$info_schedule
  # x$i <- ((stats::qnorm(false_pos[x$k]) + stats::qnorm(x$beta)) / x$delta)^2 * x$info_schedule

  # Calculate corresponding boundaries based on maximum information
  bounds <- gs_bounds_crt(x$delta, x$i, false_neg, false_pos, x$test_sides,
                        binding, x$tol, x$r)
  x$lower_bound <- bounds$a
  x$upper_bound <- bounds$b

  # Save information for trial with no interim analyses
  x$fix_i <- i0

  return(x)
}


# gs_max_size_crt function [sinew] ----
gs_max_size_crt <- function(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par) {
  # Timing of clusters per arm - 2 x k matrix
  x$m_schedule <- m_set_schedule(x$m_timing, x$k)

  # Timing of average cluster size per arm - 2 x k matrix
  x$n_schedule <- n_set_schedule(x$n_timing, x$k)

  # Set information for initial guess
  x$i0 <- set_i0(x)

  # Solve for maximum sample size
  if (x$size_type == 1) {
    # Calculate initial guess for number of clusters
    x$m0 <- m_diff(i = x$i0, m_alloc = x$m_alloc, n_fix = x$max_n, n_cv = x$n_cv,
                  var_vec = x$var_vec, rho = x$rho)
    x$m0_total <- sum(x$m0)

    # Solve for maximum sample size
    obj_func_m <- function(m_total) {
      gs_size_obj_crt(x = x, m_total = m_total, n = x$max_n,
                   alpha_sf = alpha_sf, alpha_sf_par = alpha_sf_par,
                   beta_sf = beta_sf, beta_sf_par = beta_sf_par)
    }
    x$max_m_total <- stats::uniroot(obj_func_m, lower = x$m0_total,
                                    upper = x$m0_total * 10,
                                    extendInt = "yes")$root
    x$max_m <- x$max_m_total * x$m_alloc
  } else if (x$size_type == 2) {
    # Check feasibility of trial
    x$sufficient <- check_suff(m_total = x$max_m_total, m_alloc = x$m_alloc,
                              n_cv = x$n_cv, max_i = x$i0, var_vec = x$var_vec,
                              rho = x$rho, alpha = x$alpha, beta = x$beta)

    if (x$sufficient == 1) {
      # Average initial guess for average cluster sizes
      x$n0 <- n_diff(i = x$i0, m_fix = x$max_m, n_cv = x$n_cv,
                    var_vec = x$var_vec, rho = x$rho)

      # Solve for maximum sample size
      obj_func_n <- function(n) {
        gs_size_obj_crt(x = x, m_total = x$max_m_total, n = n,
                     alpha_sf = alpha_sf, alpha_sf_par = alpha_sf_par,
                     beta_sf = beta_sf, beta_sf_par = beta_sf_par)
      }
      x$max_n <- stats::uniroot(obj_func_n, lower = x$n0, upper = x$n0 * 10,
                                extendInt = "yes")$root
    }
  } else {
    stop("invalid sample size variable")
  }

  # Get information levels for maximum sample size
  x$m <- x$max_m * x$m_schedule
  x$n <- x$max_n * x$n_schedule
  x$info <- i_diff(x$m, x$n, x$n_cv, x$var_vec, x$rho)
  x$info_schedule <- x$info / x$info[x$k]
  return(x)
}


# gs_size_obj_crt function [sinew] ----
gs_size_obj_crt <- function(x, m_total, n, alpha_sf, alpha_sf_par, beta_sf,
                         beta_sf_par) {
  # Confirm sample sizes to information levels
  m <- m_total * x$m_alloc
  x$m <- m * x$m_schedule
  x$n <- n * x$n_schedule
  if (any(x$m < 0)) {
    return(1e6) # Return large value if negative cluster size
  }
  if (any(x$n < 0)) {
    return(1e6) # Return large value if negative cluster size
  }

  x$info <- i_diff(x$m, x$n, x$n_cv, x$var_vec, x$rho)
  x$info_schedule <- x$info / x$info[x$k]

  # Partition error probabilities based on error spending
  part <- partition_error(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par)
  false_pos <- part$false_pos
  false_neg <- part$false_neg
  binding <- part$binding

  # Solve for stopping boundaries
  bounds <- gs_bounds_crt(x$delta, x$info, false_neg, false_pos, x$test_sides,
                        binding, x$tol, x$r)
  return(bounds$diff[x$k])
}


# gs_expect_crt function [sinew] ----
gs_expect_crt <- function(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par) {
  # Specify timing depending on provided information
  if (x$timing_type == 1) {
    if (x$recruit_type == 1) {
      # Recruit clusters with fixed sizes
      x$n <- matrix(x$max_n, nrow = 2, ncol = x$k)
      x$m <- sapply(x$i, function(i) {
        m_diff(i, x$m_alloc, x$max_n, x$n_cv, x$var_vec, x$rho)
      })
    } else if (x$recruit_type == 2) {
      # Recruit individuals into fixed number of clusters
      x$m <- matrix(x$max_m, nrow = 2, ncol = x$k)
      x$n <- sapply(x$i, function(i) {
        n_diff(i, x$max_m, x$n_cv, x$var_vec, x$rho)
      })
      x$n <- rbind(x$n, x$n)
    } else if (x$recruit_type == 3) {
      # Recruit at both cluster and individual levels
      n_schedule <- n_set_schedule(x$n_timing, x$k)
      x$n <- x$max_n * n_schedule
      x$m <- sapply(1:x$k, function(k_i) {
        m_diff(x$i[k_i], x$m_alloc, x$n[, k_i], x$n_cv, x$var_vec, x$rho)
      })
    }
  } else {
    # Timing of clusters per arm - 2 x k matrix
    x$m_schedule <- m_set_schedule(x$m_timing, x$k)
    x$m <- x$max_m * x$m_schedule

    # Timing of average cluster size per arm - 2 x k matrix
    x$n_schedule <- n_set_schedule(x$n_timing, x$k)
    x$n <- x$max_n * x$n_schedule
  }

  # Round sample sizes to integers and get corresponding information
  x$m[, x$k] <- x$max_m
  x$n[, x$k] <- x$max_n

  x$info <- i_diff(m = x$m, n = x$n, n_cv = x$n_cv, var_vec = x$var_vec,
                  rho = x$rho)
  x$info_schedule <- x$info / x$max_i
  x$info_schedule[x$info_schedule > 1] <- 1

  # Partition error probabilities based on error spending
  part <- partition_error(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par)
  false_pos <- part$false_pos
  false_neg <- part$false_neg
  binding <- part$binding

  # Determine boundary values to achieve Type I error = alpha
  bounds <- gs_bounds_crt(x$delta, x$info, false_neg, false_pos, x$test_sides,
                        binding, x$tol, x$r)
  x$lower_bound <- bounds$a
  x$upper_bound <- bounds$b

  # Compute crossing probabilities
  theta <- as.double(c(0, x$delta)) # For H0 and H1
  y <- gs_prob_crt(theta, x$info, x$lower_bound, x$upper_bound, x$test_sides, x$r)
  x$power <- y$power
  x$futile <- y$futile

  # Compute expected sample sizes from crossing probabilities
  p_nc <- sapply(t(rep(1, length(theta))) - x$power - x$futile,
                 function(x) max(x, 0))

  if (x$k == 1) {
    x$e_m <- cbind(x$m, x$m)
    x$e_n <- cbind(x$n, x$n)
    x$e_total <- cbind(x$m * x$n, x$m * x$n)
  } else {
    x$e_m <- x$m %*% (y$prob_lo + y$prob_hi) +
      (x$m[, x$k] %*% t(p_nc))
    x$e_n <- x$n %*% (y$prob_lo + y$prob_hi) +
      (x$n[, x$k] %*% t(p_nc))
    x$e_total <- (x$m * x$n) %*% (y$prob_lo + y$prob_hi) +
      (x$m[, x$k] * x$n[, x$k] %*% t(p_nc))
  }
  return(x)
}


### Hidden sub-functions for gs_design_crt ###

# set_i0 function [sinew] ----
set_i0 <- function(x) {
  if (x$test_sides == 1) {
    i0 <- ((stats::qnorm(x$alpha) + stats::qnorm(x$beta)) / x$delta)^2
  } else {
    i0 <- ((stats::qnorm(x$alpha / 2) + stats::qnorm(x$beta)) / x$delta)^2
  }
  return(i0)
}

# partition_error function [sinew] ----
partition_error <- function(x, alpha_sf, alpha_sf_par, beta_sf, beta_sf_par) {
  k <- x$k
  test_sides <- x$test_sides
  test_type  <- x$test_type

  if (x$k == 1) {
    # No interim analyses, so all error spent at final analysis
    if (test_sides == 1) {
      false_pos <- c(x$alpha)
    } else {
      false_pos <- c(x$alpha / 2)
    }
    false_neg <- c(x$beta)
    binding <- TRUE
    return(list(false_pos = false_pos, false_neg = false_neg, binding = binding))
  } else if (test_type == 1) {
    # Early stop for efficacy only
    if (test_sides == 1) {
      upper <- alpha_sf(x$alpha, x$info_schedule, alpha_sf_par)
    } else {
      upper <- alpha_sf(x$alpha / 2, x$info_schedule, alpha_sf_par)
    }
    false_pos <- upper$spend
    false_pos <- false_pos - c(0, false_pos[1:(k - 1)])
    false_neg <- c(rep(1e-15, k - 1), x$beta)
    binding  <- TRUE

  } else if (test_type == 2 || test_type == 3) {
    # Futility only
    if (x$test_sides == 1) {
      false_pos <- c(rep(1e-15, x$k - 1), x$alpha)
    } else {
      false_pos <- c(rep(1e-15, x$k - 1), x$alpha / 2)
    }
    lower <- beta_sf(x$beta, x$info_schedule, beta_sf_par)
    false_neg <- lower$spend
    false_neg <- false_neg - c(0, false_neg[1:(k - 1)])
    binding  <- (test_type == 2)

  } else if (test_type == 4 || test_type == 5) {
    # Both upper and lower
    if (x$test_sides == 1) {
      upper <- alpha_sf(x$alpha, x$info_schedule, alpha_sf_par)
    } else {
      upper <- alpha_sf(x$alpha / 2, x$info_schedule, alpha_sf_par)
    }
    false_pos <- upper$spend
    false_pos <- false_pos - c(0, false_pos[1:(k - 1)])

    lower <- beta_sf(x$beta, x$info_schedule, beta_sf_par)
    false_neg <- lower$spend
    false_neg <- false_neg - c(0, false_neg[1:(k - 1)])
    binding  <- (test_type == 4)

  } else {
    stop("Unsupported test_type: ", test_type)
  }

  return(list(false_pos = false_pos, false_neg = false_neg, binding = binding))
}

# i_diff function [sinew] ----
i_diff <- function(m, n, n_cv, var_vec, rho) {
  k <- dim(m)[2]
  d_eff <- (1 + ((1 + ((m - 1) / m) * n_cv^2) * n - 1) * rho)
  i <- colSums(matrix(var_vec, nrow = 2, ncol = k) * d_eff / (m * n))^(-1)
  return(i)
}

# m_diff function [sinew] ----
m_diff <- function(i, m_alloc, n_fix, n_cv, var_vec, rho) {
  d_eff <- (1 + ((1 + n_cv^2) * n_fix - 1) * rho)
  ai <- sum(i * var_vec * d_eff / (m_alloc * n_fix))
  bi <- sum(i * var_vec * n_cv^2 * rho / m_alloc^2)

  if (ai^2 - 4 * bi < 0) {
    stop("trial not feasible with given specifications")
  }
  m_total <- (ai + sqrt(ai^2 - 4 * bi)) / 2
  return(m_total * m_alloc)
}

# n_diff function [sinew] ----
n_diff <- function(i, m_fix, n_cv, var_vec, rho) {
  num <- i * var_vec * (1 - rho) / m_fix
  denom <- 1 - sum(i * var_vec *
                     (1 + ((m_fix - 1) / m_fix) * n_cv^2) * rho / m_fix)
  n <- sum(num / denom)
  return(n)
}

# check_suff function [sinew] ----
check_suff <- function(m_total, m_alloc, n_cv, max_i, var_vec, rho,
                       alpha = 0.05, beta = 0.1) {
  # Conduct feasibility check
  check <- 0
  comp <- sum(max_i * var_vec *
                (1 + ((m_alloc - 1) / m_alloc) * n_cv^2) * rho / m_alloc)
  if (m_total > comp) {
    check <- 1
  }
  return(check)
}

# info_set_schedule function [sinew] ----
info_set_schedule <- function(info_timing, k) {
  # Specify timing depending on provided information
  if (k == 1) {
    return(1)
  } else if (length(info_timing) < 1 ||
               (length(info_timing) == 1 &&
                  (k > 2 || (k == 2 && (info_timing[1] <= 0 || info_timing[1] >= 1))))) {
    return(seq(k) / k)
  } else if (length(info_timing) == k - 1 || length(info_timing) == k) {
    if (length(info_timing) == k - 1) {
      return(c(info_timing, 1))
    } else if (info_timing[k] != 1) {
      stop("if analysis timing for final analysis is input, it must be 1")
    }
    if (length(info_timing) > 1 &&
          min(info_timing - c(0, info_timing[1:(k - 1)])) <= 0) {
      stop("input timing of interim analyses must be increasing strictly
           between 0 and 1")
    }
  } else {
    stop("value input for timing must be length 1, k-1 or k")
  }
}

# m_set_schedule function [sinew] ----
m_set_schedule <- function(m_timing, k) {
  # Timing of clusters per arm - 2 x k matrix
  m_timing <- as.matrix(m_timing)
  if (dim(m_timing)[1] == 2 && dim(m_timing)[2] == 1 &&
        sum(!(m_timing %in% c(0, 1))) > 0) {
    m_timing <- rbind(t(m_timing), t(m_timing))
  } else if (dim(m_timing)[1] != 2) {
    if (dim(m_timing)[1] == 1 ||
          dim(m_timing)[1] == k && dim(m_timing)[2] == 1) {
      m_timing <- rbind(t(m_timing), t(m_timing))
    } else {
      stop("m_timing must be a matrix with 2 rows")
    }
  }
  if (dim(m_timing)[2] == 1) {
    if (m_timing[1] == 0) {
      m_timing1 <- rep(1, k)
    } else if (m_timing[1] == 1) {
      m_timing1 <- seq(k) / k
    } else {
      stop("if analysis timing is one number, it must be 0 or 1")
    }
    if (m_timing[2] == 0) {
      m_timing2 <- rep(1, k)
    } else if (m_timing[2] == 1) {
      m_timing2 <- seq(k) / k
    } else {
      stop("if analysis timing is one number, it must be 0 or 1")
    }
    m_timing <- rbind(unname(m_timing1), unname(m_timing2))
  } else if (dim(m_timing)[2] == k - 1 || dim(m_timing)[2] == k) {
    if (dim(m_timing)[2] == k - 1) {
      m_timing <- cbind(m_timing, 1)
    } else if (m_timing[1, k] != 1 || m_timing[2, k] != 1) {
      stop("if analysis timing for final analysis is input, it must be 1")
    }
    if (dim(m_timing)[2] > 1 &&
          min(m_timing - cbind(0, m_timing[, 1:(k - 1)])) <= 0) {
      stop("input timing of interim analyses must be increasing strictly
           between 0 and 1")
    }
  } else {
    stop("value input for timing must be length 1, k-1 or k")
  }
  return(m_timing)
}

# n_set_schedule function [sinew] ----
n_set_schedule <- function(n_timing, k) {
  n_timing <- as.matrix(n_timing)
  if (dim(n_timing)[1] == 2 && dim(n_timing)[2] == 1 &&
        sum(!(n_timing %in% c(0, 1))) > 0) {
    n_timing <- rbind(t(n_timing), t(n_timing))
  } else if (dim(n_timing)[1] != 2) {
    if (dim(n_timing)[1] == 1 ||
          dim(n_timing)[1] == k && dim(n_timing)[2] == 1) {
      n_timing <- rbind(t(n_timing), t(n_timing))
    } else {
      stop("n_timing must be a matrix with 2 rows")
    }
  }
  if (dim(n_timing)[2] == 1) {
    if (n_timing[1] == 0) {
      n_timing1 <- rep(1, k)
    } else if (n_timing[1] == 1) {
      n_timing1 <- seq(k) / k
    } else {
      stop("if analysis timing is one number, it must be 0 or 1")
    }
    if (n_timing[2] == 0) {
      n_timing2 <- rep(1, k)
    } else if (n_timing[2] == 1) {
      n_timing2 <- seq(k) / k
    } else {
      stop("if analysis timing is one number, it must be 0 or 1")
    }
    n_timing <- rbind(unname(n_timing1), unname(n_timing2))
  } else if (dim(n_timing)[2] == k - 1 || dim(n_timing)[2] == k) {
    if (dim(n_timing)[2] == k - 1) {
      n_timing <- cbind(n_timing, 1)
    } else if (n_timing[1, k] != 1 || n_timing[2, k] != 1) {
      stop("if analysis timing for final analysis is input, it must be 1")
    }
    if (dim(n_timing)[2] > 1 &&
          min(n_timing - cbind(0, n_timing[, 1:(k - 1)])) <= 0) {
      stop("input timing of interim analyses must be increasing strictly
          between 0 and 1")
    }
  } else {
    stop("value input for timing must be length 1, k-1 or k")
  }
  return(n_timing)
}
