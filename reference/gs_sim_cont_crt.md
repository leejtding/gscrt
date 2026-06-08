# Simulate group sequential cluster-randomized trial with continuous outcomes.

Simulate group sequential cluster-randomized trial with continuous
outcomes.

## Usage

``` r
gs_sim_cont_crt(
  k,
  data,
  test_type,
  test_sides,
  recruit_type,
  stat_type,
  delta,
  sigma_vec = c(1, 1),
  rho,
  alpha = 0.05,
  beta = 0.1,
  i_max = 1,
  m_max = c(1, 1),
  n_max = c(1, 1),
  n_cv = c(0, 0),
  schedule_m = NULL,
  schedule_n = NULL,
  alpha_sf,
  alpha_sf_par = -4,
  beta_sf,
  beta_sf_par = -4,
  tol = 1e-06,
  r = 18
)
```

## Arguments

- k:

  Number of analyses planned, including interim and final.

- data:

  Simulated continuous outcomes. Should be an n x 4 matrix with the
  columns encoding the treatment arm, cluster, individual, and response.

- test_type:

  `1=` early stopping for efficacy only  
  `2=` early stopping for binding futility only  
  `3=` early stopping for non-binding futility only  
  `4=` early stopping for either efficacy or binding futility  
  `5=` early stopping for either efficacy or non-binding futility

- test_sides:

  `1=` one-sided test  
  `2=` two-sided test

- recruit_type:

  `1=`by clusters with all individuals recruited  
  `2=`by individuals in recruited cluster  
  `3=`by both clusters and individuals in clusters

- stat_type:

  `1=` Z-test with known variance and ICC  
  `2=` Z-test with re-estimated variance and ICC  
  `3=` t-test with re-estimated variance and ICC

- delta:

  Effect size for theta under alternative hypothesis.

- sigma_vec:

  Standard deviations for control and treatment groups.

- rho:

  Intraclass correlation coefficient. Default value is 0.

- alpha:

  Desired Type I error, always one-sided. Default value is 0.05.

- beta:

  Desired Type II error, default value is 0.1 (90% power).

- i_max:

  Maximum information.

- m_max:

  Number of clusters.

- n_max:

  Mean size of each cluster.

- n_cv:

  Coefficient of variation for cluster size

- schedule_m:

  Number of clusters at each interim look. Interim analyses will be
  conducted according to the information levels in `schedule_m` and
  `schedule_n` if provided. Otherwise, interim analyses will be
  conducted at equal-sized sample increments according to
  `recruit_type`.

- schedule_n:

  Average cluster size at each interim look. Interim analyses will be
  conducted according to the information levels in `schedule_m` and
  `schedule_n` if provided. Otherwise, interim analyses will be
  conducted at equal-sized sample increments according to
  `recruit_type`.

- alpha_sf:

  A spending function or a character string indicating an upper boundary
  type (that is, “WT” for Wang-Tsiatis bounds, “OF” for O'Brien-Fleming
  bounds, and “Pocock” for Pocock bounds). The default value is
  `sf_ldof` which is a Lan-DeMets O'Brien-Fleming spending function. See
  details,
  [`vignette("spending_function_overview")`](https://leejtding.github.io/gscrt/articles/spending_function_overview.md),
  manual and examples.

- alpha_sf_par:

  Real value, default is \\-4\\ which is an O'Brien-Fleming-like
  conservative bound when used with a Hwang-Shih-DeCani spending
  function. This is a real-vector for many spending functions. The
  parameter `alpha_sf_par` specifies any parameters needed for the
  spending function specified by `alpha_sf`; this will be ignored for
  spending functions (`sf_ldof`, `sf_ld_pocock`) or bound types (“OF”,
  “Pocock”) that do not require parameters.

- beta_sf:

  A spending function or a character string indicating an lower boundary
  type (that is, “WT” for Wang-Tsiatis bounds, “OF” for O'Brien-Fleming
  bounds, and “Pocock” for Pocock bounds). The default value is
  `sf_ldof` which is a Lan-DeMets O'Brien-Fleming spending function. See
  details,
  [`vignette("spending_function_overview")`](https://leejtding.github.io/gscrt/articles/spending_function_overview.md),
  manual and examples.

- beta_sf_par:

  Real value, default is \\-4\\ which is an O'Brien-Fleming-like
  conservative bound when used with a Hwang-Shih-DeCani spending
  function. This is a real-vector for many spending functions. The
  parameter `beta_sf_par` specifies any parameters needed for the
  spending function specified by `beta_sf`; this will be ignored for
  spending functions (`sf_ldof`, `sf_ld_pocock`) or bound types (“OF”,
  “Pocock”) that do not require parameters.

- tol:

  Tolerance for error (default is 0.000001). Normally this will not be
  changed by the user. This does not translate directly to number of
  digits of accuracy, so use extra decimal places.

- r:

  Integer value controlling grid for numerical integration as in
  Jennison and Turnbull (2000); default is 18, range is 1 to 80. Larger
  values provide larger number of grid points and greater accuracy.
  Normally `r` will not be changed by the user.

## Value

Object containing the following elements:

- reject:

  Whether the null hypothesis was rejected in the simulated trial.

- k_i:

  Interim analysis at which simulated trial was stopped.

- m_i:

  Number of clusters per arm when the simulated trial was stopped.

- n_i:

  Average number of individuals per cluster when the simulated trial was
  stopped.

- total_i:

  Total number of individuals per arm when the simulated trial was
  stopped.

- i_frac:

  Information fraction when the simulated trial was stopped.

## Author

Lee Ding <lee_ding@g.harvard.edu>
