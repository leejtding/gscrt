# Simulate group sequential cluster-randomized trial.

Simulate group sequential cluster-randomized trial.

## Usage

``` r
gs_sim_crt(design, data, stat = c("Z_known", "Z_unknown", "t_unknown"))
```

## Arguments

- design:

  Object of class `gs_design_crt` specifying the group sequential
  cluster-randomized trial design.

- data:

  Simulated outcomes. Should be an n x 4 matrix with the columns
  encoding the treatment arm, cluster, individual, and response.

- stat:

  Type of test statistic to use. Options are `"Z_known"`, `"Z_unknown"`,
  and `"t_unknown"`.

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
