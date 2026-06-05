context("base tests")

### Utility functions for gsDesign stress tests

"alpha.beta.range.util" <- function(alpha, beta, test_type, sf) {
  no.err <- TRUE

  for (a in alpha) {
    for (b in beta) {
      if (b < 1 - a - 0.1) {
        res <- try(gs_design_crt(test_type = test_type,
                               alpha = a, beta = b,
                               alpha_sf = sf, beta_sf = sf))

        if (inherits(res, "try-error")) {
          no.err <- FALSE
        }
      }
    }
  }

  no.err
}

"param.range.util" <- function(param, test_type, sf) {
  no.err <- TRUE

  for (p in param)
  {
    res <- try(gs_design_crt(test_type = test_type,
                           alpha_sf = sf, alpha_sf_par = p,
                           beta_sf = sf, beta_sf_par = p))

    if (inherits(res, "try-error")) {
      no.err <- FALSE
    }
  }

  no.err
}

a1 <- round(seq(from = 0.05, to = 0.95, by = 0.05), 2)
a2 <- round(seq(from = 0.05, to = 0.45, by = 0.05), 2)
b <- round(seq(from = 0.05, to = 0.95, by = 0.05), 2)

# nu: sf_exponential parameter
nu <- round(seq(from = 0.1, to = 1.5, by = 0.1), 1)

# rho: sf_power parameter
rho <- round(seq(from = 1, to = 15, by = 1), 0)

# gamma: sf_hsd parameter
gamma <- round(seq(from = -5, to = 5, by = 1), 0)

######### gs_upper_crt tests #########

test_that("test.gs_upper_crt.theta", {
  expect_error(gs_upper_crt(theta = "abc", I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_upper_crt(theta = rep(1, 2), I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable length"
  )
  expect_error(gs_upper_crt(
    theta = -1, I = 1, false_pos = 0.5, sides = 2
  ), info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_upper_crt.I", {
  expect_error(gs_upper_crt(I = "abc", false_pos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_upper_crt(I = 0, false_pos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_upper_crt(I = rep(1, 2), false_pos = 0.5),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gs_upper_crt.a", {
  expect_error(gs_upper_crt(a = "abc", I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_upper_crt(a = rep(0.5, 2), I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable length"
  )
  expect_error(gs_upper_crt(a = -0.5, I = 1, false_pos = 0.5, sides = 2),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_upper_crt.false_pos", {
  expect_error(gs_upper_crt(false_pos = "abc", I = 1),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_upper_crt(false_pos = 1, I = 1),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_upper_crt(false_pos = rep(0.5, 2), I = 1),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gs_upper_crt.sides", {
  expect_error(gs_upper_crt(sides = "abc", I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_upper_crt(sides = 0, I = 1, false_pos = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_upper_crt.tol", {
  expect_error(gs_upper_crt(tol = "abc", I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_upper_crt(tol = 0, I = 1, false_pos = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_upper_crt.r", {
  expect_error(gs_upper_crt(r = "abc", I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_upper_crt(r = 0, I = 1, false_pos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_upper_crt(r = 81, I = 1, false_pos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_upper_crt(r = rep(1, 2), I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable length"
  )
})

######## gs_lower_crt tests #########

test_that("test.gs_lower_crt.theta", {
  expect_error(gs_lower_crt(theta = "abc", I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_lower_crt(theta = rep(1, 2), I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable length"
  )
  expect_error(gs_lower_crt(
    theta = -1, I = 1, false_neg = 0.5, sides = 2
  ), info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_lower_crt.I", {
  expect_error(gs_lower_crt(I = "abc", false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_lower_crt(I = 0, false_neg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_lower_crt(I = rep(1, 2), false_neg = 0.5),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gs_lower_crt.b", {
  expect_error(gs_lower_crt(b = "abc", I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_lower_crt(b = rep(0.5, 2), I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable length"
  )
  expect_error(gs_lower_crt(b = -0.5, I = 1, false_neg = 0.5, sides = 2),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_lower_crt.false_neg", {
  expect_error(gs_lower_crt(false_neg = "abc", I = 1),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_lower_crt(false_neg = 1, I = 1),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_lower_crt(false_neg = rep(0.5, 2), I = 1),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gs_lower_crt.sides", {
  expect_error(gs_lower_crt(sides = "abc", I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_lower_crt(sides = 0, I = 1, false_neg = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_lower_crt.binding", {
  expect_error(gs_lower_crt(binding = 2, I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_lower_crt(binding = "TRUE", I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
})

test_that("test.gs_lower_crt.tol", {
  expect_error(gs_lower_crt(tol = "abc", I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_lower_crt(tol = 0, I = 1, false_neg = 0.5),
    info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_lower_crt.r", {
  expect_error(gs_lower_crt(r = "abc", I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_lower_crt(r = 0, I = 1, false_neg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_lower_crt(r = 81, I = 1, false_neg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_lower_crt(r = rep(1, 2), I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable length"
  )
})

test_that("test.gs_bounds_crt.theta", {
  expect_error(gs_bounds_crt(
    theta = "abc", I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gs_bounds_crt(
    theta = rep(1, 2), I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable length"
  )
  expect_error(gs_bounds_crt(
    theta = -1, I = 1, false_neg = 0.5, false_pos = 0.5, sides = 2
  ), info = "Checking for out-of-range variable value"
  )
})

######## gs_bounds_crt tests #########

test_that("test.gs_bounds_crt.I", {
  expect_error(gs_bounds_crt(I = "abc", false_neg = 0.5, false_pos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_bounds_crt(I = 0, false_neg = 0.5, false_pos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_bounds_crt(
    I = rep(1, 2), false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable length"
  )
})

test_that("test.gs_bounds_crt.false_pos", {
  expect_error(gs_bounds_crt(false_pos = "abc", I = 1, false_neg = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_bounds_crt(false_pos = 1, I = 1, false_neg = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_bounds_crt(
    false_pos = rep(0.5, 2), I = 1, false_neg = 0.5
  ), info = "Checking for incorrect variable length"
  )
})

test_that("test.gs_bounds_crt.false_neg", {
  expect_error(gs_bounds_crt(false_neg = "abc", I = 1, false_pos = 0.5),
    info = "Checking for incorrect variable type"
  )
  expect_error(gs_bounds_crt(false_neg = 1, I = 1, false_pos = 0.5),
    info = "Checking for out-of-range variable value"
  )
  expect_error(gs_bounds_crt(
    false_neg = rep(0.5, 2), I = 1, false_pos = 0.5
  ), info = "Checking for incorrect variable length"
  )
})

test_that("test.gs_bounds_crt.sides", {
  expect_error(gs_bounds_crt(
    sides = "abc", I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gs_bounds_crt(
    sides = 0, I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_bounds_crt.binding", {
  expect_error(gs_bounds_crt(
    binding = 2, I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gs_bounds_crt(
    binding = "TRUE", I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable type"
  )
})

test_that("test.gs_bounds_crt.tol", {
  expect_error(gs_bounds_crt(
    tol = "abc", I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gs_bounds_crt(
    tol = 0, I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for out-of-range variable value"
  )
})

test_that("test.gs_bounds_crt.r", {
  expect_error(gs_bounds_crt(
    r = "abc", I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable type"
  )
  expect_error(gs_bounds_crt(
    r = 0, I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for out-of-range variable value"
  )
  expect_error(gs_bounds_crt(
    r = 81, I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for out-of-range variable value"
  )
  expect_error(gs_bounds_crt(
    r = rep(1, 2), I = 1, false_neg = 0.5, false_pos = 0.5
  ), info = "Checking for incorrect variable length"
  )
})

######### gs_design_crt tests #########

test_that("test.gs_design_crt.k", {
  expect_error(gs_design_crt(k = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(k = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(k = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(k = seq(2)),
               info = "Checking for incorrect variable length")
  expect_error(gs_design_crt(
    k = 3, alpha_sf = sf_points, alpha_sf_par = c(0.05, 0.1, 0.15, 0.2, 1)
  ), info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.outcome.type", {
  expect_error(gs_design_crt(outcome_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(outcome_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(outcome_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(outcome_type = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(outcome_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.test.type", {
  expect_error(gs_design_crt(test_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(test_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(test_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(test_type = 6),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(test_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.test.sides", {
  expect_error(gs_design_crt(test_sides = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(test_sides = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(test_sides = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(test_sides = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(test_sides = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.size.type", {
  expect_error(gs_design_crt(size_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(size_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(size_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(size_type = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(size_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.timing.type", {
  expect_error(gs_design_crt(timing_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(timing_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(timing_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 4),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.recruit.type", {
  expect_error(gs_design_crt(recruit_type = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(recruit_type = 1.2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(recruit_type = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(recruit_type = 4),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(recruit_type = seq(2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.delta", {
  expect_error(gs_design_crt(delta = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(delta = -1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(delta = 2, outcome_type = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(delta = rep(0.1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.sigma_vec", {
  expect_error(gs_design_crt(outcome_type = 1, sigma_vec = c("abc", 1)),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(outcome_type = 1, sigma_vec = 1),
               info = "Checking for incorrect variable length")
  expect_error(gs_design_crt(outcome_type = 1, sigma_vec = seq(3)),
               info = "Checking for incorrect variable length")
  expect_error(gs_design_crt(outcome_type = 1, sigma_vec = c(1, -1)),
               info = "Checking for out-of-range variable value")
})

test_that("test.gs_design_crt.p_vec", {
  expect_error(gs_design_crt(outcome_type = 2, p_vec = c("abc", 0.5)),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(outcome_type = 2, p_vec = 0.5),
               info = "Checking for incorrect variable length")
  expect_error(gs_design_crt(outcome_type = 2, p_vec = rep(0.5, 3)),
               info = "Checking for incorrect variable length")
  expect_error(gs_design_crt(outcome_type = 2, p_vec = c(0.5, -0.5)),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(outcome_type = 2, p_vec = c(1.5, 0.5)),
               info = "Checking for out-of-range variable value")
})

test_that("test.gs_design_crt.rho", {
  expect_error(gs_design_crt(rho = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(rho = -0.5),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(rho = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(rho = rep(1, 3)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.alpha", {
  expect_error(gs_design_crt(alpha = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(alpha = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(alpha = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(alpha = rep(0.5, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.beta", {
  expect_error(gs_design_crt(beta = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(beta = 0.5, alpha = 0.5),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(beta = 1, alpha = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(beta = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(beta = rep(0.1, 2), alpha = 0.5),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.m", {
  expect_error(gs_design_crt(m = "abc", size_type = 2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(m = -1, size_type = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(m = rep(2, 3), size_type = 2),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.m.alloc", {
  expect_error(gs_design_crt(m_alloc = "abc", size_type = 2),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(m_alloc = c(0.5, 0.6), size_type = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(m_alloc = rep(0.5, 3), size_type = 2),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.n", {
  expect_error(gs_design_crt(n = "abc", size_type = 1),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(n = -1, size_type = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(n = rep(2, 3), size_type = 1),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.n.cv", {
  expect_error(gs_design_crt(n_cv = "abc", size_type = 1),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(n_cv = -1, size_type = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(n_cv = rep(2, 3), size_type = 1),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.info.timing", {
  expect_error(gs_design_crt(timing_type = 1, info_timing = "abc", k = 1),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(timing_type = 1, info_timing = -1, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 1, info_timing = 2, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 1, info_timing = c(0.1, 1.1), k = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 1, info_timing = c(0.5, 0.1), k = 2),
               info = "Checking for incorrect variable specification")
  expect_error(gs_design_crt(timing_type = 1, info_timing = c(0.1, 0.5, 1),
                           k = 2),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.m.timing", {
  expect_error(gs_design_crt(timing_type = 3, m_timing = "abc", k = 1),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(timing_type = 3, m_timing = -1, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 3, m_timing = 2, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 3,
                           m_timing = matrix(c(0.5, 0.5, 1, 2), nrow = 2),
                           k = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 3,
                           m_timing = matrix(c(0.5, 0.5, 0.25, 0.25, 1, 1),
                                             nrow = 2),
                           k = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 3,
                           m_timing = matrix(c(0.25, 0.25, 0.5, 0.5, 1, 2),
                                             nrow = 2),
                           k = 2),
               info = "Checking for incorrect variable specification")
})

test_that("test.gs_design_crt.n.timing", {
  expect_error(gs_design_crt(timing_type = 3, n_timing = "abc", k = 1),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(timing_type = 3, n_timing = -1, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 3, n_timing = 2, k = 1),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 3,
                           n_timing = matrix(c(0.5, 0.5, 1, 2), nrow = 2),
                           k = 2),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 3,
                           n_timing = matrix(c(0.5, 0.5, 0.25, 0.25, 1, 1),
                                             nrow = 2),
                           k = 3),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(timing_type = 3,
                           n_timing = matrix(c(0.25, 0.25, 0.5, 0.5, 1, 2),
                                             nrow = 2),
                           k = 2),
               info = "Checking for incorrect variable specification")
})

test_that("test.gs_design_crt.alpha_sf", {
  expect_error(gs_design_crt(alpha_sf = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(alpha_sf = rep(sf_ldof, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.beta_sf", {
  expect_error(gs_design_crt(beta_sf = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(beta_sf = rep(sf_ldof, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.r", {
  expect_error(gs_design_crt(r = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(r = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(r = 81),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(r = rep(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_design_crt.tol", {
  expect_error(gs_design_crt(tol = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_design_crt(tol = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(tol = 0.10000001),
               info = "Checking for out-of-range variable value")
  expect_error(gs_design_crt(tol = rep(0.1, 2)),
               info = "Checking for incorrect variable length")
})

######### Stress tests for spending functions #########

test_that("test.stress.sf_exp.type1", {
  no_errors <- param.range.util(param = nu, test_type = 1, sf = sf_exponential)
  expect_true(no_errors, info = "Type 1 sf_exponential stress test")
})

test_that("test.stress.sf_exp.type2", {
  no_errors <- param.range.util(param = nu, test_type = 2, sf = sf_exponential)
  expect_true(no_errors, info = "Type 2 sf_exponential stress test")
})

test_that("test.stress.sf_exp.type3", {
  no_errors <- param.range.util(param = nu, test_type = 3, sf = sf_exponential)
  expect_true(no_errors, info = "Type 3 sf_exponential stress test")
})

test_that("test.stress.sf_exp.type4", {
  no_errors <- param.range.util(param = nu, test_type = 4, sf = sf_exponential)
  expect_true(no_errors, info = "Type 4 sf_exponential stress test")
})

test_that("test.stress.sf_exp.type5", {
  no_errors <- param.range.util(param = nu, test_type = 5, sf = sf_exponential)
  expect_true(no_errors, info = "Type 5 sf_exponential stress test")
})

test_that("test.stress.sf_hsd.type1", {
  no_errors <- param.range.util(param = gamma, test_type = 1, sf = sf_hsd)
  expect_true(no_errors, info = "Type 1 sf_hsd stress test")
})

test_that("test.stress.sf_hsd.type2", {
  no_errors <- param.range.util(param = gamma, test_type = 2, sf = sf_hsd)
  expect_true(no_errors, info = "Type 2 sf_hsd stress test")
})

test_that("test.stress.sf_hsd.type3", {
  no_errors <- param.range.util(param = gamma, test_type = 3, sf = sf_hsd)
  expect_true(no_errors, info = "Type 3 sf_hsd stress test")
})

test_that("test.stress.sf_hsd.type4", {
  no_errors <- param.range.util(param = gamma, test_type = 4, sf = sf_hsd)
  expect_true(no_errors, info = "Type 4 sf_hsd stress test")
})

test_that("test.stress.sf_hsd.type5", {
  no_errors <- param.range.util(param = gamma, test_type = 5, sf = sf_hsd)
  expect_true(no_errors, info = "Type 5 sf_hsd stress test")
})

test_that("test.stress.sf_ldof.type1", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 1, sf = sf_ldof
  )
  expect_true(no_errors, info = "Type 1 LDOF stress test")
})

test_that("test.stress.sf_ldof.type2", {
  no_errors <- alpha.beta.range.util(
    alpha = a2, beta = b,
    test_type = 2, sf = sf_ldof
  )
  expect_true(no_errors, info = "Type 2 LDOF stress test")
})

test_that("test.stress.sf_ldof.type3", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 3, sf = sf_ldof
  )
  expect_true(no_errors, info = "Type 3 LDOF stress test")
})

test_that("test.stress.sf_ldof.type4", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 4, sf = sf_ldof
  )
  expect_true(no_errors, info = "Type 4 LDOF stress test")
})

test_that("test.stress.sf_ldof.type5", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 5, sf = sf_ldof
  )
  expect_true(no_errors, info = "Type 5 LDOF stress test")
})

test_that("test.stress.sf_ld_pocock.type1", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 1, sf = sf_ld_pocock
  )
  expect_true(no_errors, info = "Type 1 LDPocock stress test")
})

test_that("test.stress.sf_ld_pocock.type2", {
  no_errors <- alpha.beta.range.util(
    alpha = a2, beta = b,
    test_type = 2, sf = sf_ld_pocock
  )
  expect_true(no_errors, info = "Type 2 LDPocock stress test")
})

test_that("test.stress.sf_ld_pocock.type3", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 3, sf = sf_ld_pocock
  )
  expect_true(no_errors, info = "Type 3 LDPocock stress test")
})

test_that("test.stress.sf_ld_pocock.type4", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 4, sf = sf_ld_pocock
  )
  expect_true(no_errors, info = "Type 4 LDPocock stress test")
})

test_that("test.stress.sf_ld_pocock.type5", {
  no_errors <- alpha.beta.range.util(
    alpha = a1, beta = b,
    test_type = 5, sf = sf_ld_pocock
  )
  expect_true(no_errors, info = "Type 5 LDPocock stress test")
})

test_that("test.stress.sf_power.type1", {
  no_errors <- param.range.util(param = rho, test_type = 1, sf = sf_power)
  expect_true(no_errors, info = "Type 1 sf_power stress test")
})

test_that("test.stress.sf_power.type2", {
  no_errors <- param.range.util(param = rho, test_type = 2, sf = sf_power)
  expect_true(no_errors, info = "Type 2 sf_power stress test")
})

test_that("test.stress.sf_power.type3", {
  no_errors <- param.range.util(param = rho, test_type = 3, sf = sf_power)
  expect_true(no_errors, info = "Type 3 sf_power stress test")
})

test_that("test.stress.sf_power.type4", {
  no_errors <- param.range.util(param = rho, test_type = 4, sf = sf_power)
  expect_true(no_errors, info = "Type 4 sf_power stress test")
})

test_that("test.stress.sf_power.type5", {
  no_errors <- param.range.util(param = rho, test_type = 5, sf = sf_power)
  expect_true(no_errors, info = "Type 5 sf_power stress test")
})

######## gs_probability_crt tests #########

test_that("test.gs_probability_crt.theta", {
  expect_error(gs_probability_crt(theta = "abc"),
               info = "Checking for incorrect variable type")
})

test_that("test.gs_probability_crt.I", {
  expect_error(gs_probability_crt(I = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_probability_crt(I = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_probability_crt(I = c(2, 1), a = rep(0, 2), b = rep(1, 2)),
               info = "Checking for out-of-order input sequence")
  expect_error(gs_probability_crt(I = c(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_probability_crt.a", {
  expect_error(gs_probability_crt(a = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_probability_crt(a = c(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_probability_crt.b", {
  expect_error(gs_probability_crt(b = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_probability_crt(b = c(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.gs_probability_crt.sides", {
  expect_error(gs_probability_crt(sides = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_probability_crt(sides = 0),
               info = "Checking for out-of-range variable value")
})

test_that("test.gs_probability.r", {
  expect_error(gs_probability_crt(r = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(gs_probability_crt(r = 0),
               info = "Checking for out-of-range variable value")
  expect_error(gs_probability_crt(r = 81),
               info = "Checking for out-of-range variable value")
  expect_error(gs_probability_crt(r = rep(1, 2)),
               info = "Checking for incorrect variable length")
})

test_that("test.Deming.gsProb", {
  w <- sum(gs_probability_crt(theta = 0, I = 1:2,
                            a = stats::qnorm(0.025) * c(1, 1),
                            b = stats::qnorm(0.975) * c(1, 1))$upper$prob)
  x <- sum(gs_probability_crt(theta = 0, I = 1:4,
                            a = stats::qnorm(0.025) * rep(1, 4),
                            b = stats::qnorm(0.975) * rep(1, 4))$upper$prob)
  y <- sum(gs_probability_crt(theta = 0, I = 1:10,
                            a = stats::qnorm(0.025) * array(1, 10),
                            b = stats::qnorm(0.975) * array(1, 10))$upper$prob)
  z <- sum(gs_probability_crt(theta = 0, I = 1:20,
                            a = stats::qnorm(0.025) * array(1, 20),
                            b = stats::qnorm(0.975) * array(1, 20))$upper$prob)
  expect_equal(0.042, round(w, 3), info = "Checking Type I error, k = 2")
  expect_equal(0.063, round(x, 3), info = "Checking Type I error, k = 4")
  expect_equal(0.097, round(y, 3), info = "Checking Type I error, k = 10")
  expect_equal(0.124, round(z, 3), info = "Checking Type I error, k = 20")
})

######### Spending function parameter tests #########

test_that("test.sfcauchy.param", {
  expect_error(sfcauchy(param = rep(1, 3)),
               info = "Checking for incorrect variable length")
  expect_error(sfcauchy(
    param = c(0.1, 0.6, 0.2, 0.05), k = 5,
    timing = c(0.1, 0.25, 0.4, 0.6)
  ), info = "Checking for out-of-order input sequence")
})

test_that("test.sfcauchy.param ", {
  expect_error(sfcauchy(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfcauchy(param = c(1, 0)),
               info = "Checking for out-of-range variable value")
})

test_that("test.sfexp.param", {
  expect_error(sfexp(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfexp(param = rep(1, 2)),
               info = "Checking for incorrect variable length")
  expect_error(sfexp(param = 0),
               info = "Checking for out-of-range variable value")
  expect_error(sfexp(param = 11),
               info = "Checking for out-of-range variable value")
})

test_that("test.sf_hsd.param", {
  expect_error(sf_hsd(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sf_hsd(param = rep(1, 2)),
               info = "Checking for incorrect variable length")
  expect_error(sf_hsd(param = -41),
               info = "Checking for out-of-range variable value")
  expect_error(sf_hsd(param = 41),
               info = "Checking for out-of-range variable value")
})

test_that("test.sflogistic.param", {
  expect_error(sflogistic(param = rep(1, 3)),
               info = "Checking for incorrect variable length")
  expect_error(sflogistic(
    param = c(0.1, 0.6, 0.2, 0.05), k = 5,
    timing = c(0.1, 0.25, 0.4, 0.6)
  ), info = "Checking for out-of-order input sequence")
})

test_that("test.sflogistic.param ", {
  expect_error(sflogistic(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sflogistic(param = c(1, 0)),
               info = "Checking for out-of-range variable value")
})


test_that("test.sfnorm.param", {
  expect_error(sfnorm(param = rep(1, 3)),
               info = "Checking for incorrect variable length")
  expect_error(sfnorm(
    param = c(0.1, 0.6, 0.2, 0.05), k = 5,
    timing = c(0.1, 0.25, 0.4, 0.6)
  ), info = "Checking for out-of-order input sequence")
})

test_that("test.sfnorm.param ", {
  expect_error(sfnorm(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfnorm(param = c(1, 0)),
               info = "Checking for out-of-range variable value")
})

test_that("test.sfpower.param", {
  expect_error(sfpower(param = "abc"),
               info = "Checking for incorrect variable type")
  expect_error(sfpower(param = rep(1, 2)),
               info = "Checking for incorrect variable length")
  expect_error(sfpower(param = -1),
               info = "Checking for out-of-range variable value")
})

test_that("test.sf_t_dist.param", {
  expect_error(sf_t_dist(param = rep(1, 4)),
               info = "Checking for incorrect variable length")
  expect_error(sf_t_dist(param = c(1, 0, 1)),
               info = "Checking for out-of-range variable value")
  expect_error(sf_t_dist(param = c(1, 1, 0.5)),
               info = "Checking for out-of-range variable value")
  expect_error(sf_t_dist(param = 1, 1:3 / 4, c(
    0.25, 0.5, 0.75,
    0.1, 0.2, 0.3
  )), info = "Checking for out-of-range variable value")
})

test_that("test.sf_t_dist.param ", {
  expect_error(sf_t_dist(param = "abc"),
               info = "Checking for incorrect variable type")
})

######### analysis schedule setting tests #########
test_that("m_set_schedule handles scalar (0/1) inputs correctly", {
  k <- 3

  # Both arms: pattern 0 -> constant max clusters at each analysis
  mt0 <- m_set_schedule(c(0, 0), k = k)
  expect_equal(dim(mt0), c(2L, k))
  expect_equal(mt0[1, ], rep(1, k))
  expect_equal(mt0[2, ], rep(1, k))

  # First arm constant, second arm linear accrual
  mt01 <- m_set_schedule(c(0, 1), k = k)
  expect_equal(dim(mt01), c(2L, k))
  expect_equal(mt01[1, ], rep(1, k))
  expect_equal(mt01[2, ], seq_len(k) / k)

  # Length-1 input recycled to both arms
  mt1 <- m_set_schedule(1, k = k)
  expect_equal(dim(mt1), c(2L, k))
  expect_equal(mt1[1, ], seq_len(k) / k)
  expect_equal(mt1[2, ], seq_len(k) / k)
})

test_that("n_set_schedule handles scalar (0/1) inputs correctly", {
  k <- 4

  # Both arms: pattern 1 -> linear accrual in cluster size
  nt11 <- n_set_schedule(c(1, 1), k = k)
  expect_equal(dim(nt11), c(2L, k))
  expect_equal(nt11[1, ], seq_len(k) / k)
  expect_equal(nt11[2, ], seq_len(k) / k)

  # First arm linear, second arm constant
  nt10 <- n_set_schedule(c(1, 0), k = k)
  expect_equal(dim(nt10), c(2L, k))
  expect_equal(nt10[1, ], seq_len(k) / k)
  expect_equal(nt10[2, ], rep(1, k))
})

test_that("m_set_schedule and n_set_schedule expand 2 x (k-1) matrices and enforce final timing = 1", {
  k <- 3

  # Valid 2 x (k-1) matrix: last column should be appended as 1s
  m_in <- matrix(c(0.3, 0.3,
                   0.7, 0.7), nrow = 2)  # 2 x 2, k = 3
  mt <- m_set_schedule(m_in, k = k)
  expect_equal(dim(mt), c(2L, k))
  expect_equal(mt[, 1:2], m_in)
  expect_equal(mt[, 3], rep(1, 2))

  n_in <- matrix(c(0.25, 0.25,
                   0.6,  0.6), nrow = 2) # 2 x 2, k = 3
  nt <- n_set_schedule(n_in, k = k)
  expect_equal(dim(nt), c(2L, k))
  expect_equal(nt[, 1:2], n_in)
  expect_equal(nt[, 3], rep(1, 2))

  # Invalid: full 2 x k matrix with final timing != 1 should error
  m_bad_final <- matrix(c(0.5, 0.5,
                          0.9, 0.9), nrow = 2) # k = 2
  expect_error(
    m_set_schedule(m_bad_final, k = 2),
    info = "final analysis is input, it must be 1"
  )

  n_bad_final <- matrix(c(0.4, 0.4,
                          0.8, 0.8), nrow = 2)
  expect_error(
    n_set_schedule(n_bad_final, k = 2),
    info = "final analysis is input, it must be 1"
  )
})

test_that("m_set_schedule and n_set_schedule enforce 2-row matrices and strictly increasing timings", {
  # Wrong number of rows
  m_bad_rows <- matrix(0.5, nrow = 4, ncol = 1)
  expect_error(
    m_set_schedule(m_bad_rows, k = 3),
    info = "if analysis timing for final analysis is input, it must be 1"
  )

  n_bad_rows <- matrix(0.5, nrow = 4, ncol = 1)
  expect_error(
    n_set_schedule(n_bad_rows, k = 3),
    info = "if analysis timing for final analysis is input, it must be 1"
  )

  # Not strictly increasing: flat timing should fail
  m_flat <- matrix(c(0.5, 0.5,
                     0.5, 0.5), nrow = 2)  # 2 x 2, k = 2
  expect_error(
    m_set_schedule(m_flat, k = 2),
    info = "Must be increasing strictly"
  )

  n_flat <- matrix(c(0.3, 0.3,
                     0.3, 0.3), nrow = 2)
  expect_error(
    n_set_schedule(n_flat, k = 2),
    info = "Must be increasing strictly"
  )
})

######## gen_cluster_sizes tests #########

test_that("gen_cluster_sizes with zero CV returns constant sizes", {
  sizes <- gen_cluster_sizes(
    m     = 3,
    n     = 10,
    n_cv  = 0,
    n_min = rep(1, 3),
    n_max = rep(20, 3)
  )

  expect_equal(length(sizes), 3)
  expect_true(all(sizes == 10))
})

test_that("gen_cluster_sizes errors for invalid negative binomial dispersion", {
  # n_cv^2 <= 1 / n should error
  expect_error(
    gen_cluster_sizes(
      m     = 3,
      n     = 10,
      n_cv  = sqrt(1 / 10),  # exactly on the boundary
      n_min = rep(1, 3),
      n_max = rep(20, 3)
    ),
    info = "n_cv^2 must be greater"
  )
})

######## gen_cont_crt tests #########

test_that("gen_cont_crt generates continuous CRT data properly", {
  df <- gen_cont_crt(
    m        = c(2, 2),
    m_alloc  = c(0.5, 0.5),
    n        = c(5, 5),
    n_cv     = c(0, 0),
    mu_vec   = c(0, 1),
    sigma_vec = c(1, 1),
    rho      = c(0.01, 0.01)
  )

  expect_s3_class(df, "data.frame")
  expect_equal(colnames(df), c("arm", "cluster", "individual", "response"))

  expect_true(all(df$arm %in% 0:1))
  expect_true(all(df$cluster >= 1))
  expect_true(is.numeric(df$response))
  expect_true(nrow(df) > 0)
})

######### gen_bin_crt tests #########

test_that("gen_bin_crt generates binary CRT data properly", {
  df <- gen_bin_crt(
    m        = c(2, 2),
    m_alloc  = c(0.5, 0.5),
    n        = c(5, 5),
    n_cv     = c(0, 0),
    p_vec    = c(0.3, 0.5),
    rho      = c(0.01, 0.01)
  )

  expect_s3_class(df, "data.frame")
  expect_equal(colnames(df), c("arm", "cluster", "individual", "response"))

  expect_true(all(df$arm %in% 0:1))
  expect_true(all(df$cluster >= 1))
  expect_true(all(df$response %in% c(0, 1)))
  expect_true(nrow(df) > 0)
})

######### gs_sim_cont_crt tests #########

test_that("gs_sim_cont_crt runs with known variance/ICC and returns expected structure", {
  # Simple continuous data set compatible with m_max / n_max
  m_max <- c(3, 3)
  n_max <- c(8, 8)
  data_cont <- gen_cont_crt(
    m        = m_max,
    n        = n_max,
    n_cv     = c(0, 0),
    mu_vec   = c(0, 0.5),
    sigma_vec = c(1, 1),
    rho      = c(0.01, 0.01)
  )

  sim <- gs_sim_cont_crt(
    k           = 3,
    data        = data_cont,
    test_type   = 1,
    test_sides  = 1,
    recruit_type = 1,
    stat_type   = 1,                 # Z_known
    delta       = 0.5,
    sigma_vec   = c(1, 1),
    rho         = c(0.01, 0.01),
    alpha       = 0.05,
    beta        = 0.1,
    i_max       = 1,                  # any positive value; just affects ifrac
    m_max       = m_max,
    n_max       = n_max,
    n_cv        = c(0, 0),
    schedule_m  = NULL,               # let function choose equal fractions
    schedule_n  = NULL,
    alpha_sf    = sf_ldof,
    beta_sf     = sf_ldof
  )

  expect_type(sim, "list")
  expect_true(all(c("reject", "t_i", "m_i", "n_i", "total_i",
                    "i_frac") %in% names(sim)))

  expect_true(is.logical(sim$reject))
  expect_true(length(sim$reject) == 1L)
  expect_true(is.numeric(sim$t_i))
  expect_true(sim$t_i >= 1)
  expect_true(is.numeric(sim$m_i))
  expect_true(is.numeric(sim$n_i))
  expect_true(is.numeric(sim$total_i))
})

######### gs_sim_bin_crt tests #########

test_that("gs_sim_bin_crt runs with known variance/ICC and returns expected structure", {
  set.seed(123)

  m_max <- c(3L, 3L)
  n_max <- c(8L, 8L)
  data_bin <- gen_bin_crt(
    m        = m_max,
    n        = n_max,
    n_cv     = c(0, 0),
    p_vec    = c(0.3, 0.5),
    rho      = c(0.01, 0.01)
  )

  sim <- gs_sim_bin_crt(
    k           = 3L,
    data        = data_bin,
    test_type   = 1L,
    test_sides  = 1L,
    recruit_type = 1L,
    stat_type   = 1L,                 # Z_known
    delta       = 0.2,
    p_vec       = c(0.3, 0.5),
    rho         = c(0.01, 0.01),
    alpha       = 0.05,
    beta        = 0.1,
    i_max       = 1,
    m_max       = m_max,
    n_max       = n_max,
    n_cv        = c(0, 0),
    schedule_m  = NULL,
    schedule_n  = NULL,
    alpha_sf    = sf_ldof,
    beta_sf     = sf_ldof
  )

  expect_type(sim, "list")
  expect_true(all(c("reject", "t_i", "m_i", "n_i", "total_i",
                    "i_frac") %in% names(sim)))

  expect_true(is.logical(sim$reject))
  expect_true(length(sim$reject) == 1L)
  expect_true(is.numeric(sim$t_i))
  expect_true(sim$t_i >= 1)
  expect_true(is.numeric(sim$m_i))
  expect_true(is.numeric(sim$n_i))
  expect_true(is.numeric(sim$total_i))
})