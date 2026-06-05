# Simulate cluster-randomized trial data with continuous outcomes

Simulate cluster-randomized trial data with continuous outcomes

## Usage

``` r
gen_cont_crt(
  m = c(1, 1),
  m_alloc = c(0.5, 0.5),
  n = c(1, 1),
  n_cv = c(0, 0),
  mu_vec = c(0, 1),
  sigma_vec = c(1, 1),
  rho = c(0, 0)
)
```

## Arguments

- m:

  Number of clusters.

- m_alloc:

  Allocation ratio of clusters per arm. Default is `c(0.5, 0.5)`.

- n:

  Mean size of each cluster.

- mu_vec:

  Vector of means for control and treatment groups, respectively.

- sigma_vec:

  Vector of standard deviations for control and treatment groups,
  respectively.

- rho:

  Intraclass correlation coefficient. Default value is 0.

## Value

Simulated continuous outcomes represented as an n x 4 matrix with the
columns encoding the treatment arm, cluster, individual, and response.

## Author

Lee Ding <lee_ding@g.harvard.edu>
