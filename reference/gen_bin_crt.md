# Simulate cluster-randomized trial data with binary outcomes

Simulate cluster-randomized trial data with binary outcomes

## Usage

``` r
gen_bin_crt(
  m = c(1, 1),
  m_alloc = c(0.5, 0.5),
  n = c(1, 1),
  n_cv = c(0, 0),
  p_vec = c(0.5, 0.5),
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

- p_vec:

  Probabilities of event for control and treatment groups.

- rho:

  Intraclass correlation coefficient. Default value is 0.

## Value

Simulated binary outcomes represented as an n x 4 matrix with the
columns encoding the treatment arm, cluster, individual, and response.

## Author

Lee Ding <lee_ding@g.harvard.edu>
