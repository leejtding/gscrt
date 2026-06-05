# Compute stopping boundaries, maximum sample size, and expected sample sizes for a group sequential cluster randomized trial.

`gs_design_crt()` is used to determine the maximum sample size needed
for a specified parallel group sequential cluster randomized trial to
detect a clinically meaningful effect size with some Type I error rate
and power. Code adapted from gsDesign package.

## Usage

``` r
gs_design_crt(
  k = 3,
  outcome_type = 1,
  test_type = 1,
  test_sides = 1,
  size_type = 1,
  timing_type = 2,
  recruit_type = 1,
  delta = 1,
  sigma_vec = c(1, 1),
  p_vec = c(0.5, 0.5),
  rho = 0,
  alpha = 0.05,
  beta = 0.1,
  m = 2,
  m_alloc = c(0.5, 0.5),
  n = 1,
  n_cv = c(0, 0),
  info_timing = 1,
  m_timing = c(1, 1),
  n_timing = c(1, 1),
  alpha_sf = sf_ldof,
  alpha_sf_par = -4,
  beta_sf = sf_ldof,
  beta_sf_par = -4,
  tol = 1e-06,
  r = 18
)
```

## Arguments

- k:

  Number of analyses planned, including interim and final.

- outcome_type:

  `1=`continuous difference of means  
  `2=`binary difference of proportions

- test_type:

  `1=` early stopping for efficacy only  
  `2=` early stopping for binding futility only  
  `3=` early stopping for non-binding futility only  
  `4=` early stopping for either efficacy or binding futility  
  `5=` early stopping for either efficacy or non-binding futility

- test_sides:

  `1=` one-sided test  
  `2=` two-sided test

- size_type:

  `1=`clusters per arm  
  `2=`cluster size

- timing_type:

  `1=` maximum and expected sample sizes based on specified information
  fractions in `info_timing`; recruit according to design specified in
  `recruit_type`  
  `2=` maximum sample size based on specified information fractions in
  `info_timing` and expected sample sizes based on specified sample size
  fractions in `m_timing` and `n_timing`  
  `3=` maximum and expected sample sizes based on specified sample size
  fractions in `m_timing` and `n_timing` @param recruit_type `1=`recruit
  clusters with fixed sizes  
  `2=`recruit individuals into fixed number of clusters  
  `3=`recruit at both cluster and individual levels

- delta:

  Effect size for theta under alternative hypothesis. Must be \> 0.

- sigma_vec:

  Standard deviations for control and treatment groups (continuous
  case).

- p_vec:

  Probabilities of event for control and treatment groups (binary case).

- rho:

  Intraclass correlation coefficient. Default value is 0.

- alpha:

  Type I error, default value is 0.05.

- beta:

  Type II error, default value is 0.1 (90% power).

- m:

  Number of clusters for finding maximum mean cluster size. If `m` is a
  scalar, it is treated as the total number of clusters across arms. If
  `m` is a vector of length 2, it is treated as the number of clusters
  per arm.

- m_alloc:

  Allocation ratio of clusters per arm. Default is `c(0.5, 0.5)`.

- n:

  Mean cluster size for finding maximum number of clusters. If `n` is a
  scalar, it is treated as the average cluster size for both arms. If
  `n` is a vector of length 2, it is treated as the average cluster size
  per arm.

- n_cv:

  Coefficient of variation for cluster size. If `n_cv` is a scalar, it
  is treated as the coefficient of variation for both arms. If `n_cv` is
  a vector of length 2, it is treated as the coefficient of variation
  per arm.

- info_timing:

  Sets timing of interim analyses based on information fractions.
  Default of 1 produces analyses at equal-spaced increments. Otherwise,
  this is a vector of length `k` or `k-1`. The values should satisfy
  `0 < info_timing[1] < info_timing[2] < ... < info_timing[k-1] < info_timing[k]=1`.

- m_timing:

  Sets timing of interim analyses based on fractions of the number of
  clusters.

- n_timing:

  Sets timing of interim analyses based on cluster size

- alpha_sf:

  A spending function or a character string indicating an upper boundary
  type (that is, “WT” for Wang-Tsiatis bounds, “OF” for O'Brien-Fleming
  bounds, and “Pocock” for Pocock bounds). The default value is
  `sf_ldof` which is a Lan-DeMets O'Brien-Fleming spending function. See
  details,
  [`vignette("spending_function_overview")`](https://leejtding.github.io/gsDesignCRT/articles/spending_function_overview.md),
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
  [`vignette("spending_function_overview")`](https://leejtding.github.io/gsDesignCRT/articles/spending_function_overview.md),
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

- k:

  As input.

- outcome_type:

  As input.

- test_type:

  As input.

- test_sides:

  As input.

- size_type:

  As input.

- recruit_type:

  As input.

- timing_type:

  As input.

- delta:

  As input.

- sigma_vec:

  As input.

- p_vec:

  As input.

- rho:

  As input.

- alpha:

  As input.

- beta:

  As input.

- info_timing:

  As input.

- m_timing:

  As input.

- n_timing:

  As input.

- info_schedule:

  Fraction of maximum information at each planned interim analysis based
  on `timing_type` and `info_timing`.

- m_schedule:

  Fraction of maximum number of clusters per arm at each planned interim
  analysis based on `timing_type` and `m_timing`.

- n_schedule:

  Fraction of maximum cluster size at each planned interim analysis
  based on `timing_type` and `n_timing`.

- i:

  Fisher information at each planned interim analysis based on
  `timing_type`.

- max_i:

  Maximum information corresponding to design specifications.

- m:

  Number of clusters per arm at each planned interim analysis.

- max_m:

  Maximum number of clusters per arm.

- e_m:

  A vector of length 2 with expected number of clusters per arm under
  the null and alternative hypotheses. For simplicity, the expected
  sizes with non-binding futility boundaries are calculated assuming the
  boundaries are binding futility.

- n:

  Mean cluster size at each planned interim analysis.

- max_n:

  Maximum mean cluster size.

- e_n:

  A vector of length 2 with expected mean cluster sizes under the null
  and alternative hypotheses. For simplicity, the expected sizes with
  non-binding futility boundaries are calculated assuming the futility
  boundaries are binding.

- max_total:

  Maximum number of individuals in the trial.

- e_total:

  A vector of length 2 with expected total number of individuals in the
  trial under the null and alternative hypotheses. For simplicity, the
  expected sizes with non-binding futility boundaries are calculated
  assuming the futility boundaries are binding.

- sufficient:

  Value denoting whether calculated sample size will be sufficient to
  achieve specified Type I error rate and power given the trial
  specifications.

- lower_bound:

  Calculated lower futility boundaries under analysis schedule specified
  by `timing_type`

- upper_bound:

  Calculated upper efficacy boundaries under analysis schedule specified
  by `timing_types`.

- tol:

  As input.

- r:

  As input.

## References

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

Anderson K (2025). gsDesign: Group Sequential Design. R package version
3.8.0, https://keaven.github.io/gsDesign/.

## Author

Lee Ding <lee_ding@g.harvard.edu>
