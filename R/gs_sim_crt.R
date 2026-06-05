# gs_sim_crt roxy [sinew] ----
#' @title Simulate group sequential cluster-randomized trial.
#'
#' @param design Object of class \code{gs_design_crt} specifying the group
#' sequential cluster-randomized trial design.
#' @param data Simulated outcomes. Should be an n x 4 matrix with the
#' columns encoding the treatment arm, cluster, individual, and response.
#' @param stat Type of test statistic to use. Options are \code{"Z_known"},
#'  \code{"Z_unknown"}, and \code{"t_unknown"}.
#'
#' @return Object containing the following elements: \item{reject}{Whether the
#' null hypothesis was rejected in the simulated trial.} \item{k_i}{Interim
#' analysis at which simulated trial was stopped.} \item{m_i}{Number of clusters
#' per arm when the simulated trial was stopped.} \item{n_i}{Average number of
#' individuals per cluster when the simulated trial was stopped.}
#' \item{total_i}{Total number of individuals per arm when the simulated trial
#' was stopped.} \item{i_frac}{Information fraction when the simulated trial
#' was stopped.}
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @export
#' @name gs_sim_crt
# gs_sim_crt function [sinew] ----
gs_sim_crt <- function(design, data,
                     stat = c("Z_known", "Z_unknown", "t_unknown")) {
  stopifnot(inherits(design, "gs_design_crt"))
  stat <- match.arg(stat)
  if (design$outcome_type == 1) {
    # continuous
    out <- gs_sim_cont_crt(
      k           = design$k,
      data        = data,
      test_type   = design$test_type,
      test_sides  = design$test_sides,
      recruit_type = design$recruit_type,
      stat_type   = match(stat, c("Z_known", "Z_reestimate", "t_reestimate")),
      delta       = design$delta,
      sigma_vec   = design$sigma_vec,
      rho         = design$rho,
      alpha       = design$alpha,
      beta        = design$beta,
      i_max       = design$max_i,
      m_max       = ceiling(design$max_m),
      n_max       = ceiling(design$max_n),
      n_cv        = design$n_cv,
      schedule_m  = design$m_schedule,
      schedule_n  = design$n_schedule,
      alpha_sf    = design$alpha_sf,
      alpha_sf_par = design$alpha_sf_par,
      beta_sf     = design$beta_sf,
      beta_sf_par  = design$beta_sf_par,
      tol         = design$tol,
      r           = design$r
    )
  } else {
    # binary
    out <- gs_sim_bin_crt(
      k           = design$k,
      data        = data,
      test_type   = design$test_type,
      test_sides  = design$test_sides,
      recruit_type = design$recruit_type,
      stat_type   = match(stat, c("Z_known", "Z_reestimate", "t_reestimate")),
      delta       = design$delta,
      p_vec       = design$p_vec,
      rho         = design$rho,
      alpha       = design$alpha,
      beta        = design$beta,
      i_max       = design$max_i,
      m_max       = ceiling(design$max_m),
      n_max       = ceiling(design$max_n),
      n_cv        = design$n_cv,
      schedule_m  = design$m_schedule,
      schedule_n  = design$n_schedule,
      alpha_sf    = design$alpha_sf,
      alpha_sf_par = design$alpha_sf_par,
      beta_sf     = design$beta_sf,
      beta_sf_par  = design$beta_sf_par,
      tol         = design$tol,
      r           = design$r
    )
  }

  class(out) <- c("gs_sim_crt", class(out))
  out
}

# gs_sim_cont_crt roxy [sinew] ----
#' @title Simulate group sequential cluster-randomized trial with continuous
#' outcomes.
#'
#' @param k Number of analyses planned, including interim and final.
#' @param data Simulated continuous outcomes. Should be an n x 4 matrix with the
#' columns encoding the treatment arm, cluster, individual, and response.
#' @param test_type \code{1=} early stopping for efficacy only
#' \cr \code{2=} early stopping for binding futility only
#' \cr \code{3=} early stopping for non-binding futility only
#' \cr \code{4=} early stopping for either efficacy or binding futility
#' \cr \code{5=} early stopping for either efficacy or non-binding futility
#' @param test_sides \code{1=} one-sided test \cr \code{2=} two-sided test
#' @param recruit_type \code{1=}by clusters with all individuals recruited
#' \cr \code{2=}by individuals in recruited cluster
#' \cr \code{3=}by both clusters and individuals in clusters
#' @param stat_type \code{1=} Z-test with known variance and ICC
#' \cr \code{2=} Z-test with re-estimated variance and ICC
#' \cr \code{3=} t-test with re-estimated variance and ICC
#' @param delta Effect size for theta under alternative hypothesis.
#' @param sigma_vec Standard deviations for control and treatment groups.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#' @param alpha Desired Type I error, always one-sided. Default value is 0.05.
#' @param beta Desired Type II error, default value is 0.1 (90\% power).
#' @param i_max Maximum information.
#' @param m_max Number of clusters.
#' @param n_max Mean size of each cluster.
#' @param n_cv Coefficient of variation for cluster size
#' @param schedule_m Number of clusters at each interim look. Interim analyses
#' will be conducted according to the information levels in \code{schedule_m}
#' and \code{schedule_n} if provided. Otherwise, interim analyses will be
#' conducted at equal-sized sample increments according to \code{recruit_type}.
#' @param schedule_n Average cluster size at each interim look. Interim
#' analyses will be conducted according to the information levels in
#' \code{schedule_m} and \code{schedule_n} if provided. Otherwise, interim
#' analyses will be conducted at equal-sized sample increments according to
#' \code{recruit_type}.
#' @param alpha_sf A spending function or a character string indicating an
#' upper boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#' \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#' bounds). The default value is \code{sf_ldof} which is a Lan-DeMets
#' O'Brien-Fleming spending function. See details,
#' \code{vignette("spending_function_overview")}, manual and examples.
#' @param alpha_sf_par Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with a
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions. The parameter \code{alpha_sf_par} specifies any parameters needed
#' for the spending function specified by \code{alpha_sf}; this will be ignored
#' for spending functions (\code{sf_ldof}, \code{sf_ld_pocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param beta_sf A spending function or a character string indicating an
#' lower boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#' \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#' bounds). The default value is \code{sf_ldof} which is a Lan-DeMets
#' O'Brien-Fleming spending function. See details,
#' \code{vignette("spending_function_overview")}, manual and examples.
#' @param beta_sf_par Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with a
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions. The parameter \code{beta_sf_par} specifies any parameters needed
#' for the spending function specified by \code{beta_sf}; this will be ignored
#' for spending functions (\code{sf_ldof}, \code{sf_ld_pocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#' be changed by the user.  This does not translate directly to number of
#' digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#' Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#' values provide larger number of grid points and greater accuracy.  Normally
#' \code{r} will not be changed by the user.
#'
#' @return Object containing the following elements: \item{reject}{Whether the
#' null hypothesis was rejected in the simulated trial.} \item{k_i}{Interim
#' analysis at which simulated trial was stopped.} \item{m_i}{Number of clusters
#' per arm when the simulated trial was stopped.} \item{n_i}{Average number of
#' individuals per cluster when the simulated trial was stopped.}
#' \item{total_i}{Total number of individuals per arm when the simulated trial
#' was stopped.} \item{i_frac}{Information fraction when the simulated trial
#' was stopped.}
#'
#' @importFrom lme4 lmer
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @export
#' @name gs_sim_cont_crt
# gs_sim_cont_crt function [sinew] ----
gs_sim_cont_crt <- function(k, data, test_type, test_sides, recruit_type,
                         stat_type, delta, sigma_vec = c(1, 1), rho,
                         alpha = 0.05, beta = 0.1, i_max = 1,
                         m_max = c(1, 1), n_max = c(1, 1), n_cv = c(0, 0),
                         schedule_m = NULL, schedule_n = NULL,
                         alpha_sf, alpha_sf_par = -4, beta_sf, beta_sf_par = -4,
                         tol = 0.000001, r = 18) {
  # Check inputs
  check_scalar(k, "integer", c(1, Inf))
  check_scalar(test_type, "integer", c(1, 5))
  check_scalar(test_sides, "integer", c(1, 2))
  check_scalar(recruit_type, "integer", c(1, 3))
  check_scalar(stat_type, "integer", c(1, 3))
  check_scalar(delta, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(sigma_vec, "numeric", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(rho, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)
  check_scalar(alpha, "numeric", 0:1, c(FALSE, FALSE))
  check_scalar(beta, "numeric", c(0, 1 - alpha), c(FALSE, FALSE))
  check_scalar(i_max, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(m_max, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(n_max, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(n_cv, "numeric", c(0, Inf), c(TRUE, FALSE), length = 2)
  check_scalar(tol, "numeric", c(0, 0.1), c(FALSE, TRUE))
  check_scalar(r, "integer", c(1, 80))

  # Set interim look parameters
  ka_i <- 0
  kb_i <- 0

  info_vec <- c()
  info_vec_ret <- c()

  ifrac_prev <- 0
  ifrac_vec <- c()
  ifrac_vec_ret <- c()

  false_pos_vec <- c()
  false_neg_vec <- c()
  stop <- FALSE

  # Equal sample size fractions if not provided
  if (is.null(schedule_m) || is.null(schedule_n)) {
    equal_frac <- seq_len(k) / k
    schedule_m <- rbind(equal_frac, equal_frac)
    schedule_n <- rbind(equal_frac, equal_frac)
  }

  # Set schedule of recruited clusters - draw from multinomial distribution
  n_max_gen <- n_max + round(5 * n_cv * n_max)
  n1_max_vec <- gen_cluster_sizes(m_max[1], n_max[1], n_cv[1],
                                rep(1, m_max[1]), rep(n_max_gen[1], m_max[1]))
  n2_max_vec <- gen_cluster_sizes(m_max[2], n_max[2], n_cv[2],
                                rep(1, m_max[2]), rep(n_max_gen[2], m_max[2]))

  while (!stop && ka_i < k) {
    # Update number of looks
    ka_i <- ka_i + 1 # number of interim looks
    kb_i <- kb_i + 1 # number of crossing bounds (usually = ka_i)

    # "Recruit" clusters
    m1 <- round(schedule_m[1, ka_i] * m_max[1])
    m2 <- round(schedule_m[2, ka_i] * m_max[2])

    n1_vec <- round(schedule_n[1, ka_i] * n1_max_vec[1:m1])
    n2_vec <- round(schedule_n[2, ka_i] * n2_max_vec[1:m2])
    data_i <- data.frame()

    size_vec1_i <- c()
    for (m_i in 1:m1) {
      size_vec1_i <- c(size_vec1_i, n1_vec[m_i])
      data1_i <- data[(data$arm == 0 & data$cluster == m_i &
                         data$individual %in% 1:n1_vec[m_i]), ]
      data_i <- rbind(data_i, data1_i)
    }

    size_vec2_i <- c()
    for (m_i in 1:m2) {
      size_vec2_i <- c(size_vec2_i, n2_vec[m_i])
      data2_i <- data[(data$arm == 1 & data$cluster == m_max[1] + m_i &
                         data$individual %in% 1:n2_vec[m_i]), ]
      data_i <- rbind(data_i, data2_i)
    }

    # Get responses
    data1_i <- data_i[data_i$arm == 0, ]
    data2_i <- data_i[data_i$arm == 1, ]

    x1_i <- data1_i$response
    x2_i <- data2_i$response

    # Compute corresponding means
    x_bar1_i <- mean(x1_i)
    x_bar2_i <- mean(x2_i)

    # Compute test statistics
    if (stat_type == 1) {
      rho_i <- rho
      se <- se_cont_diff(x1_i, x2_i, size_vec1_i, size_vec2_i, rho_i, sigma_vec)
    } else {
      fit1 <- lmer(response ~ (1 | cluster), data = data1_i)
      vc1 <- as.data.frame(VarCorr(fit1))
      var1_b <- vc1$vcov[vc1$grp == "cluster"]
      var1 <- var1_b + vc1$vcov[vc1$grp == "Residual"]
      icc_est1_i <- var1_b / var1
      if (length(icc_est1_i) > 1) {
        rho1_i <- icc_est1_i$ICC_adjusted
      } else {
        rho1_i <- 0
      }
      fit2 <- lmer(response ~ (1 | cluster), data = data2_i)
      vc2 <- as.data.frame(VarCorr(fit2))
      var2_b <- vc2$vcov[vc2$grp == "cluster"]
      var2 <- var2_b + vc2$vcov[vc2$grp == "Residual"]
      icc_est2_i <- var2_b / var2
      if (length(icc_est2_i) > 1) {
        rho2_i <- icc_est2_i$ICC_adjusted
      } else {
        rho2_i <- 0
      }
      rho_i <- c(rho1_i, rho2_i)
      se <- se_cont_diff(x1_i, x2_i, size_vec1_i, size_vec2_i, rho_i, NULL)
    }

    z_i <- (x_bar2_i - x_bar1_i) / se
    if (test_sides == 2) {
      z_i <- abs(z_i)
    }

    # Conduct group sequential test
    info <- 1 / (se^2)
    info_vec_ret <- c(info_vec_ret, info)

    ifrac <- min(info / i_max, 1)
    ifrac_vec_ret <- c(ifrac_vec_ret, ifrac)

    if (ifrac <= ifrac_prev) {
      kb_i <- kb_i - 1
    } else {
      if (ifrac == 1) {
        stop <- TRUE
      }
      if (test_type == 1) {
        # Alpha spending
        if (test_sides == 1) {
          spend_i <- alpha_sf(alpha, ifrac, alpha_sf_par)
          spend_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sf_par)
        } else {
          spend_i <- alpha_sf(alpha / 2, ifrac, alpha_sf_par)
          spend_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sf_par)
        }
        false_pos_i <- c(false_pos_vec, spend_i$spend - spend_i_prev$spend)

        # Beta spending
        if (kb_i == k) {
          false_neg_i <- c(false_neg_vec, beta)
        } else {
          false_neg_i <- c(false_neg_vec, 1e-15)
        }

        binding <- TRUE
      } else if (test_type == 2 || test_type == 3) {
        # Alpha spending
        if (kb_i == k) {
          false_pos_i <- c(false_pos_vec, 1e-15)
        } else {
          if (test_sides == 1) {
            false_pos_i <- c(false_pos_vec, alpha)
          } else {
            false_pos_i <- c(false_pos_vec, alpha / 2)
          }
        }

        # Beta spending
        spend_i <- beta_sf(beta, ifrac, beta_sf_par)
        spend_i_prev <- beta_sf(beta, ifrac_prev, beta_sf_par)
        false_neg_i <- c(false_neg_vec, spend_i$spend - spend_i_prev$spend)

        # Get stopping bounds depending on test
        if (test_type == 3) {
          binding <- FALSE
        } else {
          binding <- TRUE
        }
      } else if (test_type == 4 || test_type == 5) {
        # Alpha spending
        if (test_sides == 1) {
          spend1_i <- alpha_sf(alpha, ifrac, alpha_sf_par)
          spend1_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sf_par)
        } else {
          spend1_i <- alpha_sf(alpha / 2, ifrac, alpha_sf_par)
          spend1_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sf_par)
        }
        false_pos_i <- c(false_pos_vec, spend1_i$spend - spend1_i_prev$spend)

        # Beta spending
        spend2_i <- beta_sf(beta, ifrac, beta_sf_par)
        spend2_i_prev <- beta_sf(beta, ifrac_prev, beta_sf_par)
        false_neg_i <- c(false_neg_vec, spend2_i$spend - spend2_i_prev$spend)

        # Get stopping bounds depending on test
        if (test_type == 5) {
          binding <- FALSE
        } else {
          binding <- TRUE
        }
      }

      bounds_i <- gs_bounds_crt(theta = delta, I = c(info_vec, info),
                              false_neg = false_neg_i, false_pos = false_pos_i,
                              sides = test_sides, binding = binding,
                              tol = tol, r = r)
      if (stat_type <= 2) {
        lb_i <- bounds_i$a[kb_i]
        ub_i <- bounds_i$b[kb_i]
      } else {
        df_i <- (m1 + m2) - 2
        lb_i <- qt(p = pnorm(bounds_i$a[kb_i]), df = df_i)
        ub_i <- qt(p = pnorm(bounds_i$b[kb_i]), df = df_i)
      }

      # Conduct group sequential test for interim look
      if (z_i <= lb_i) {
        reject <- FALSE
        stop <- TRUE
      } else if (z_i >= ub_i) {
        reject <- TRUE
        stop <- TRUE
      } else {
        reject <- FALSE
      }

      # Update previous information fractions
      info_vec <- c(info_vec, info)
      ifrac_prev <- ifrac
      ifrac_vec <- c(ifrac_vec, ifrac)
      false_pos_vec <- false_pos_i
      false_neg_vec <- false_neg_i
    }
  }

  # Return results
  m_i <- mean(length(size_vec1_i), length(size_vec2_i))
  n_i <- mean(c(size_vec1_i, size_vec2_i))
  total_i <- sum(size_vec1_i) + sum(size_vec2_i)
  out <- list("reject" = reject,
              "t_i" = ka_i,
              "m_i" = m_i,
              "n_i" = n_i,
              "total_i" = total_i,
              "i_frac" = ifrac_vec)
  return(out)
}

# gs_sim_bin_crt roxy [sinew] ----
#' @title Simulate group sequential cluster-randomized trial with binary
#' outcomes
#'
#' @param k Number of analyses planned, including interim and final.
#' @param data Simulated binary outcomes. Should be an n x 4 matrix with the
#' columns encoding the treatment arm, cluster, individual, and response.
#' @param test_type \code{1=} early stopping for efficacy only
#' \cr \code{2=} early stopping for binding futility only
#' \cr \code{3=} early stopping for non-binding futility only
#' \cr \code{4=} early stopping for either efficacy or binding futility
#' \cr \code{5=} early stopping for either efficacy or non-binding futility
#' @param test_sides \code{1=} one-sided test \cr \code{2=} two-sided test
#' @param recruit_type \code{1=}by clusters with all individuals recruited
#' \cr \code{2=}by individuals in recruited cluster
#' \cr \code{3=}by both clusters and individuals in clusters
#' @param stat_type \code{1=} Z-test with known variance and ICC
#' \cr \code{2=} Z-test with re-estimated variance and ICC
#' \cr \code{3=} t-test with re-estimated variance and ICC
#' \cr \code{2=} randomized sample increments from multinomial distribution
#' according to the scheduled interim analyses.
#' @param delta Effect size for theta under alternative hypothesis.
#' @param p_vec Probabilities of event for control and treatment groups.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#' @param alpha Type I error, always one-sided. Default value is 0.05.
#' @param beta Type II error, default value is 0.1 (90\% power).
#' @param i_max Maximum information.
#' @param m_max Number of clusters.
#' @param n_max Mean size of each cluster.
#' @param schedule_m Number of clusters at each interim look. Interim analyses
#' will be conducted according to the information levels in \code{schedule_m}
#' and \code{schedule_n} if provided. Otherwise, interim analyses will be
#' conducted at equal-sized sample increments according to \code{recruit_type}.
#' @param schedule_n Average cluster size at each interim look. Interim
#' analyses will be conducted according to the information levels in
#' \code{schedule_m} and \code{schedule_n} if provided. Otherwise, interim
#' analyses will be conducted at equal-sized sample increments according to
#' \code{recruit_type}.
#' @param alpha_sf A spending function or a character string indicating an
#' upper boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#' \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#' bounds). The default value is \code{sf_ldof} which is a Lan-DeMets
#' O'Brien-Fleming spending function. See details,
#' \code{vignette("spending_function_overview")}, manual and examples.
#' @param alpha_sf_par Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with a
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions. The parameter \code{alpha_sf_par} specifies any parameters needed
#' for the spending function specified by \code{alpha_sf}; this will be ignored
#' for spending functions (\code{sf_ldof}, \code{sf_ld_pocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param beta_sf A spending function or a character string indicating an
#' lower boundary type (that is, \dQuote{WT} for Wang-Tsiatis bounds,
#' \dQuote{OF} for O'Brien-Fleming bounds, and \dQuote{Pocock} for Pocock
#' bounds). The default value is \code{sf_ldof} which is a Lan-DeMets
#' O'Brien-Fleming spending function. See details,
#' \code{vignette("spending_function_overview")}, manual and examples.
#' @param beta_sf_par Real value, default is \eqn{-4} which is an
#' O'Brien-Fleming-like conservative bound when used with a
#' Hwang-Shih-DeCani spending function. This is a real-vector for many spending
#' functions. The parameter \code{beta_sf_par} specifies any parameters needed
#' for the spending function specified by \code{beta_sf}; this will be ignored
#' for spending functions (\code{sf_ldof}, \code{sf_ld_pocock}) or bound types
#' (\dQuote{OF}, \dQuote{Pocock}) that do not require parameters.
#' @param tol Tolerance for error (default is 0.000001). Normally this will not
#' be changed by the user.  This does not translate directly to number of
#' digits of accuracy, so use extra decimal places.
#' @param r Integer value controlling grid for numerical integration as in
#' Jennison and Turnbull (2000); default is 18, range is 1 to 80.  Larger
#' values provide larger number of grid points and greater accuracy.  Normally
#' \code{r} will not be changed by the user.
#'
#' @return Object containing the following elements: \item{reject}{Whether the
#' null hypothesis was rejected in the simulated trial.} \item{k_i}{Interim
#' analysis at which simulated trial was stopped.} \item{m_i}{Number of clusters
#' per arm when the simulated trial was stopped.} \item{n_i}{Average number of
#' individuals per cluster when the simulated trial was stopped.}
#' \item{total_i}{Total number of individuals per arm when the simulated trial
#' was stopped.} \item{i_frac}{Information fraction when the simulated trial
#' was stopped.}
#'
#' @importFrom lme4 lmer
#'
#' @export
#' @name gs_sim_bin_crt
# gs_sim_bin_crt function [sinew] ----
gs_sim_bin_crt <- function(k, data, test_type, test_sides, recruit_type,
                        stat_type, delta, p_vec, rho, alpha = 0.05, beta = 0.1,
                        i_max = 1, m_max = c(1, 1), n_max = c(1, 1),
                        n_cv = c(0, 0), schedule_m = NULL, schedule_n = NULL,
                        alpha_sf, alpha_sf_par = -4, beta_sf, beta_sf_par = -4,
                        tol = 0.000001, r = 18) {
  # Check inputs
  check_scalar(k, "integer", c(1, Inf))
  check_scalar(test_type, "integer", c(1, 5))
  check_scalar(test_sides, "integer", c(1, 2))
  check_scalar(recruit_type, "integer", c(1, 3))
  check_scalar(stat_type, "integer", c(1, 3))
  check_scalar(delta, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(p_vec, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)
  check_vector(rho, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)
  check_scalar(alpha, "numeric", 0:1, c(FALSE, FALSE))
  check_scalar(beta, "numeric", c(0, 1 - alpha), c(FALSE, FALSE))
  check_scalar(i_max, "numeric", c(0, Inf), c(FALSE, FALSE))
  check_vector(m_max, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(n_max, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(n_cv, "numeric", c(0, Inf), c(TRUE, FALSE), length = 2)
  check_scalar(tol, "numeric", c(0, 0.1), c(FALSE, TRUE))
  check_scalar(r, "integer", c(1, 80))

  # Set interim look parameters
  ka_i <- 0
  kb_i <- 0

  info_vec <- c()
  info_vec_ret <- c()

  ifrac_prev <- 0
  ifrac_vec <- c()
  ifrac_vec_ret <- c()

  false_pos_vec <- c()
  false_neg_vec <- c()
  stop <- FALSE

  # Equal sample size fractions if not provided
  if (is.null(schedule_m) || is.null(schedule_n)) {
    equal_frac <- seq_len(k) / k
    schedule_m <- rbind(equal_frac, equal_frac)
    schedule_n <- rbind(equal_frac, equal_frac)
  }

  # Set schedule of recruited clusters - draw from multinomial distribution
  n_max_gen <- n_max + round(5 * n_cv * n_max)
  n1_max_vec <- gen_cluster_sizes(m_max[1], n_max[1], n_cv[1],
                                rep(1, m_max[1]), rep(n_max_gen[1], m_max[1]))
  n2_max_vec <- gen_cluster_sizes(m_max[2], n_max[2], n_cv[2],
                                rep(1, m_max[2]), rep(n_max_gen[2], m_max[2]))

  while (!stop && ka_i < k) {
    # Update number of looks
    ka_i <- ka_i + 1 # number of interim looks
    kb_i <- kb_i + 1 # number of crossing bounds (usually = ka_i)

    # "Recruit" clusters
    m1 <- round(schedule_m[1, ka_i] * m_max[1])
    m2 <- round(schedule_m[2, ka_i] * m_max[2])

    n1_vec <- round(schedule_n[1, ka_i] * n1_max_vec[1:m1])
    n2_vec <- round(schedule_n[2, ka_i] * n2_max_vec[1:m2])
    data_i <- data.frame()

    size_vec1_i <- c()
    for (m_i in 1:m1) {
      size_vec1_i <- c(size_vec1_i, n1_vec[m_i])
      data1_i <- data[(data$arm == 0 & data$cluster == m_i &
                         data$individual %in% 1:n1_vec[m_i]), ]
      data_i <- rbind(data_i, data1_i)
    }

    size_vec2_i <- c()
    for (m_i in 1:m2) {
      size_vec2_i <- c(size_vec2_i, n2_vec[m_i])
      data2_i <- data[(data$arm == 1 & data$cluster == m_max[1] + m_i &
                         data$individual %in% 1:n2_vec[m_i]), ]
      data_i <- rbind(data_i, data2_i)
    }

    # Get responses
    data1_i <- data_i[data_i$arm == 0, ]
    data2_i <- data_i[data_i$arm == 1, ]

    x1_i <- data1_i$response
    x2_i <- data2_i$response

    # Compute corresponding proportions
    p_hat1_i <- mean(x1_i)
    p_hat2_i <- mean(x2_i)

    # Compute test statistics
    if (stat_type == 1) {
      rho_i <- rho
      se <- se_bin_diff(x1_i, x2_i, size_vec1_i, size_vec2_i, rho_i, p_vec)
    } else {
      fit1 <- lmer(response ~ (1 | cluster), data = data1_i)
      vc1 <- as.data.frame(VarCorr(fit1))
      var1_b <- vc1$vcov[vc1$grp == "cluster"]
      var1 <- var1_b + vc1$vcov[vc1$grp == "Residual"]
      icc_est1_i <- var1_b / var1
      if (length(icc_est1_i) > 1) {
        rho1_i <- icc_est1_i$ICC_adjusted
      } else {
        rho1_i <- 0
      }
      fit2 <- lmer(response ~ (1 | cluster), data = data2_i)
      vc2 <- as.data.frame(VarCorr(fit2))
      var2_b <- vc2$vcov[vc2$grp == "cluster"]
      var2 <- var2_b + vc2$vcov[vc2$grp == "Residual"]
      icc_est2_i <- var2_b / var2
      if (length(icc_est2_i) > 1) {
        rho2_i <- icc_est2_i$ICC_adjusted
      } else {
        rho2_i <- 0
      }
      rho_i <- c(rho1_i, rho2_i)
      se <- se_bin_diff(x1_i, x2_i, size_vec1_i, size_vec2_i, rho_i, NULL)
    }

    z_i <- (p_hat2_i - p_hat1_i) / se
    if (test_sides == 2) {
      z_i <- abs(z_i)
    }

    # Conduct group sequential test
    info <- 1 / (se^2)
    info_vec_ret <- c(info_vec_ret, info)

    ifrac <- min(info / i_max, 1)
    ifrac_vec_ret <- c(ifrac_vec_ret, ifrac)

    if (ifrac <= ifrac_prev) {
      kb_i <- kb_i - 1
    } else {
      if (ifrac == 1) {
        stop <- TRUE
      }
      if (test_type == 1) {
        # Alpha spending
        if (test_sides == 1) {
          spend_i <- alpha_sf(alpha, ifrac, alpha_sf_par)
          spend_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sf_par)
        } else {
          spend_i <- alpha_sf(alpha / 2, ifrac, alpha_sf_par)
          spend_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sf_par)
        }
        false_pos_i <- c(false_pos_vec, spend_i$spend - spend_i_prev$spend)

        # Beta spending
        if (kb_i == k) {
          false_neg_i <- c(false_neg_vec, beta)
        } else {
          false_neg_i <- c(false_neg_vec, 1e-15)
        }

        binding <- TRUE
      } else if (test_type == 2 || test_type == 3) {
        # Alpha spending
        if (kb_i == k) {
          false_pos_i <- c(false_pos_vec, 1e-15)
        } else {
          if (test_sides == 1) {
            false_pos_i <- c(false_pos_vec, alpha)
          } else {
            false_pos_i <- c(false_pos_vec, alpha / 2)
          }
        }

        # Beta spending
        spend_i <- beta_sf(beta, ifrac, beta_sf_par)
        spend_i_prev <- beta_sf(beta, ifrac_prev, beta_sf_par)
        false_neg_i <- c(false_neg_vec, spend_i$spend - spend_i_prev$spend)

        # Get stopping bounds depending on test
        if (test_type == 3) {
          binding <- FALSE
        } else {
          binding <- TRUE
        }
      } else if (test_type == 4 || test_type == 5) {
        # Alpha spending
        if (test_sides == 1) {
          spend1_i <- alpha_sf(alpha, ifrac, alpha_sf_par)
          spend1_i_prev <- alpha_sf(alpha, ifrac_prev, alpha_sf_par)
        } else {
          spend1_i <- alpha_sf(alpha / 2, ifrac, alpha_sf_par)
          spend1_i_prev <- alpha_sf(alpha / 2, ifrac_prev, alpha_sf_par)
        }
        false_pos_i <- c(false_pos_vec, spend1_i$spend - spend1_i_prev$spend)

        # Beta spending
        spend2_i <- beta_sf(beta, ifrac, beta_sf_par)
        spend2_i_prev <- beta_sf(beta, ifrac_prev, beta_sf_par)
        false_neg_i <- c(false_neg_vec, spend2_i$spend - spend2_i_prev$spend)

        # Get stopping bounds depending on test
        if (test_type == 5) {
          binding <- FALSE
        } else {
          binding <- TRUE
        }
      }

      bounds_i <- gs_bounds_crt(theta = delta, I = c(info_vec, info),
                              false_neg = false_neg_i, false_pos = false_pos_i,
                              sides = test_sides, binding = binding,
                              tol = tol, r = r)
      if (stat_type <= 2) {
        lb_i <- bounds_i$a[kb_i]
        ub_i <- bounds_i$b[kb_i]
      } else {
        df_i <- (m1 + m2) - 2
        lb_i <- qt(p = pnorm(bounds_i$a[kb_i]), df = df_i)
        ub_i <- qt(p = pnorm(bounds_i$b[kb_i]), df = df_i)
      }

      # Conduct group sequential test for interim look
      if (z_i <= lb_i) {
        reject <- FALSE
        stop <- TRUE
      } else if (z_i >= ub_i) {
        reject <- TRUE
        stop <- TRUE
      } else {
        reject <- FALSE
      }

      # Update previous information fractions
      info_vec <- c(info_vec, info)
      ifrac_prev <- ifrac
      ifrac_vec <- c(ifrac_vec, ifrac)
      false_pos_vec <- false_pos_i
      false_neg_vec <- false_neg_i
    }
  }

  # Return results
  m_i <- mean(length(size_vec1_i), length(size_vec2_i))
  n_i <- mean(c(size_vec1_i, size_vec2_i))
  total_i <- sum(size_vec1_i) + sum(size_vec2_i)
  out <- list("reject" = reject,
              "t_i" = ka_i,
              "m_i" = m_i,
              "n_i" = n_i,
              "total_i" = total_i,
              "i_frac" = ifrac_prev)
  return(out)
}

# gen_cont_crt roxy [sinew] ----
#' @title Simulate cluster-randomized trial data with continuous outcomes
#'
#' @param m Number of clusters.
#' @param m_alloc Allocation ratio of clusters per arm. Default is
#'  \code{c(0.5, 0.5)}.
#' @param n Mean size of each cluster.
#' @param mu_vec Vector of means for control and treatment groups, respectively.
#' @param sigma_vec Vector of standard deviations for control and treatment
#'  groups, respectively.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#'
#' @return Simulated continuous outcomes represented as an n x 4 matrix with the
#'  columns encoding the treatment arm, cluster, individual, and response.
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @importFrom stats rnorm
#'
#' @export
#' @name gen_cont_crt
# gen_cont_crt function [sinew] ----
gen_cont_crt <- function(m = c(1, 1), m_alloc = c(0.5, 0.5), n = c(1, 1),
                       n_cv = c(0, 0), mu_vec = c(0, 1), sigma_vec = c(1, 1),
                       rho = c(0, 0)) {
  # Check inputs
  if (length(m) == 1) { # Input is total clusters
    m <- round(m * m_alloc)
  }
  if (length(n) == 1) { # Assume same average cluster size
    n <- rep(n, 2)
  }
  if (length(n_cv) == 1) {
    n_cv <- rep(n_cv, 2)
  }
  if (length(rho) == 1) {
    rho <- rep(rho, 2)
  }
  check_vector(m, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(m_alloc, "numeric", c(0, 1), c(FALSE, FALSE), length = 2)
  check_vector(n, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(n_cv, "numeric", c(0, Inf), c(TRUE, FALSE), length = 2)
  check_vector(mu_vec, "numeric", c(-Inf, Inf), c(FALSE, FALSE), length = 2)
  check_vector(sigma_vec, "numeric", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(rho, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)

  # Specify population means
  mu1 <- mu_vec[1]
  mu2 <- mu_vec[2]

  # Specify between and within-cluster variances
  sigma_b_1 <- sqrt(rho[1] * sigma_vec[1]^2)
  sigma_w_1 <- sqrt(sigma_vec[1]^2 - sigma_b_1^2)
  sigma_b_2 <- sqrt(rho[2] * sigma_vec[2]^2)
  sigma_w_2 <- sqrt(sigma_vec[2]^2 - sigma_b_2^2)

  # Generate samples by cluster
  n <- n + round(5 * n_cv * n)

  a1 <- rep(0, m[1] * n[1])
  r1 <- rep(0, m[1] * n[1])
  c1 <- rep(1:m[1], each = n[1])
  i1 <- rep(1:n[1], m[1])

  a2 <- rep(1, m[2] * n[2])
  r2 <- rep(0, m[2] * n[2])
  c2 <- rep(m[1] + 1:m[2], each = n[2])
  i2 <- rep(1:n[2], m[2])

  for (m1i in 1:m[1]) {
    re1 <- rnorm(1, mean = 0, sd = sigma_b_1)
    r1[(m1i - 1) * n[1] + 1:n[1]] <- (mu1 + re1 +
                                        rnorm(n[1], mean = 0, sd = sigma_w_1))
  }
  for (m2i in 1:m[2]) {
    re2 <- rnorm(1, mean = 0, sd = sigma_b_2)
    r2[(m2i - 1) * n[2] + 1:n[2]] <- (mu2 + re2 +
                                        rnorm(n[2], mean = 0, sd = sigma_w_2))
  }

  # Reshape samples into dataframe with cluster assignment
  df1 <- cbind.data.frame(a1, c1, i1, r1)
  colnames(df1) <- c("arm", "cluster", "individual", "response")

  df2 <- cbind.data.frame(a2, c2, i2, r2)
  colnames(df2) <- c("arm", "cluster", "individual", "response")

  # Combine samples
  df <- rbind(df1, df2)
  df <- df[order(df[, 1], df[, 2]), ]
  return(df)
}

# gen_bin_crt roxy [sinew] ----
#' @title Simulate cluster-randomized trial data with binary outcomes
#'
#' @param m Number of clusters.
#' @param m_alloc Allocation ratio of clusters per arm. Default is
#'  \code{c(0.5, 0.5)}.
#' @param n Mean size of each cluster.
#' @param p_vec Probabilities of event for control and treatment groups.
#' @param rho Intraclass correlation coefficient. Default value is 0.
#'
#' @return Simulated binary outcomes represented as an n x 4 matrix with the
#'  columns encoding the treatment arm, cluster, individual, and response.
#'
#' @author Lee Ding \email{lee_ding@g.harvard.edu}
#'
#' @export
#' @name gen_bin_crt
# gen_bin_crt function [sinew] ----
gen_bin_crt <- function(m = c(1, 1), m_alloc = c(0.5, 0.5), n = c(1, 1),
                      n_cv = c(0, 0), p_vec = c(0.5, 0.5), rho = c(0, 0)) {
  # Check inputs
  if (length(m) == 1) { # Input is total clusters
    m <- round(m * m_alloc)
  }
  if (length(n) == 1) { # Assume same average cluster size
    n <- rep(n, 2)
  }
  if (length(n_cv) == 1) {
    n_cv <- rep(n_cv, 2)
  }
  if (length(rho) == 1) {
    rho <- rep(rho, 2)
  }
  check_vector(m, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(m_alloc, "numeric", c(0, 1), c(FALSE, FALSE), length = 2)
  check_vector(n, "integer", c(0, Inf), c(FALSE, FALSE), length = 2)
  check_vector(n_cv, "numeric", c(0, Inf), c(TRUE, FALSE), length = 2)
  check_vector(p_vec, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)
  check_vector(rho, "numeric", c(0, 1), c(TRUE, TRUE), length = 2)

  # Generate samples by cluster
  n <- n + round(5 * n_cv * n)

  a1 <- rep(0, m[1] * n[1])
  r1 <- rep(0, m[1] * n[1])
  c1 <- rep(1:m[1], each = n[1])
  i1 <- rep(1:n[1], m[1])

  a2 <- rep(1, m[2] * n[2])
  r2 <- rep(0, m[2] * n[2])
  c2 <- rep(m[1] + 1:m[2], each = n[2])
  i2 <- rep(1:n[2], m[2])

  # Generate data according to Qaqish 2003 paper
  b1 <- sim_b(p_vec[1], n[1], rho[1])
  b2 <- sim_b(p_vec[2], n[2], rho[2])

  for (m1i in 1:m[1]) {
    r1[(m1i - 1) * n[1] + 1:n[1]] <- sim_response(p_vec[1], n[1], b1)
  }
  for (m2i in 1:m[2]) {
    r2[(m2i - 1) * n[2] + 1:n[2]] <- sim_response(p_vec[2], n[2], b2)
  }

  # Reshape samples into dataframe with cluster assignment
  df1 <- cbind.data.frame(a1, c1, i1, r1)
  colnames(df1) <- c("arm", "cluster", "individual", "response")

  df2 <- cbind.data.frame(a2, c2, i2, r2)
  colnames(df2) <- c("arm", "cluster", "individual", "response")

  # Combine samples
  df <- rbind(df1, df2)
  df <- df[order(df[, 1], df[, 2]), ]
  return(df)
}

# sim_b function [sinew] ----
sim_b <- function(p, n, rho) {
  B <- NULL
  R <- (1 - rho) * diag(n) + rho * matrix(1, nrow = n, ncol = n)
  u <- rep(p, each = n)
  A_half <- diag(sqrt(u * (1 - u)))
  v <- A_half %*% R %*% A_half
  b <- v
  for (f in 2:n) {
    f1 <- f - 1
    gf <- v[1:f1, 1:f1]
    sf <- v[1:f1, f]
    bf <- solve(gf, sf)
    b[1:f1, f] <- bf
  }
  B <- cbind(B, b)
  return(B)
}

# sim_response function [sinew] ----
#' @importFrom stats rbinom
# sim_response function [sinew] ----
sim_response <- function(p, n, B) {
  u <- rep(p, each = n)
  y_out <- rep(-1, n)
  y_out[1] <- rbinom(1, 1, u[1])
  for (l in 2:n) {
    l1 <- l - 1
    res <- y_out[1:l1] - u[1:l1]
    cl <- u[l] + sum(res * B[1:l1, l])
    y_out[l] <- rbinom(1, 1, cl)
  }
  y <- c(y_out)
  return(y)
}

# gen_cluster_sizes roxy [sinew] ----
#' @importFrom stats rnbinom
# gen_cluster_sizes function [sinew] ----
gen_cluster_sizes <- function(m, n, n_cv, n_min, n_max) {
  # Check inputs
  check_scalar(m, "integer", c(0, Inf), c(FALSE, FALSE))
  check_scalar(n, "integer", c(0, Inf), c(FALSE, FALSE))
  check_scalar(n_cv, "numeric", c(0, Inf), c(TRUE, FALSE))
  check_vector(n_min, "integer", c(1, Inf), c(TRUE, FALSE), length = m)
  check_vector(n_max, "integer", c(1, Inf), c(TRUE, FALSE), length = m)

  if (n_cv == 0) {
    return(rep(n, m))
  } else {
    if (n_cv^2 <= 1 / n) {
      stop(
        "In gen_cluster_sizes(): n_cv^2 must be greater than 1 / n ",
        "to obtain a valid negative binomial dispersion."
      )
    }
    r_nb <- 1 / (n_cv^2 - (1 / n))
    p_nb <- r_nb / (r_nb + n)

    cluster_sizes <- rnbinom(m, size = r_nb, prob = p_nb)
    cluster_sizes <- pmax(cluster_sizes, n_min)
    cluster_sizes <- pmin(cluster_sizes, n_max)
    return(cluster_sizes)
  }
}

# seContCRT roxy [sinew] ----
#' @importFrom stats var
# se_cont_diff function [sinew] ----
se_cont_diff <- function(x1, x2, size_vec1, size_vec2, rho, sigma_vec = NULL) {
  n_sq_vec1 <- size_vec1^2
  n_sq_vec2 <- size_vec2^2

  if (is.null(sigma_vec)) {
    pooled_num <- ((sum(size_vec1) - 1) * var(x1) +
                     (sum(size_vec2) - 1) * var(x2))
    pooled_denom <- (sum(size_vec1) - 1) + (sum(size_vec2) - 1)
    pooled <- pooled_num / pooled_denom
    v1 <- pooled / sum(size_vec1)
    v2 <- pooled / sum(size_vec2)
  } else {
    v1 <- sigma_vec[1]^2 / sum(size_vec1)
    v2 <- sigma_vec[2]^2 / sum(size_vec2)
  }

  d_eff1 <- (1 + ((sum(n_sq_vec1) / sum(size_vec1)) - 1) * rho[1])
  d_eff2 <- (1 + ((sum(n_sq_vec2) / sum(size_vec2)) - 1) * rho[2])

  se <- sqrt(v1 * d_eff1 + v2 * d_eff2)
  return(se)
}

# se_bin_diff function [sinew] ----
se_bin_diff <- function(x1, x2, size_vec1, size_vec2, rho, p_vec = NULL) {
  n_sq_vec1 <- size_vec1^2
  n_sq_vec2 <- size_vec2^2

  if (is.null(p_vec)) {
    p <- (sum(x1) + sum(x2)) / (sum(size_vec1) + sum(size_vec2))
  } else {
    p <- mean(p_vec)
  }

  if ((p * (1 - p)) == 0) {
    warning("se_bin_diff(): p is 0 or 1; variance is 0 and SE is undefined.
            Returning NA.")
    return(NA_real_)
  } else {
    v1 <- p * (1 - p) / sum(size_vec1)
    v2 <- p * (1 - p) / sum(size_vec2)

    d_eff1 <- (1 + ((sum(n_sq_vec1) / sum(size_vec1)) - 1) * rho[1])
    d_eff2 <- (1 + ((sum(n_sq_vec2) / sum(size_vec2)) - 1) * rho[2])

    se <- sqrt(v1 * d_eff1 + v2 * d_eff2)
    return(se)
  }
}
