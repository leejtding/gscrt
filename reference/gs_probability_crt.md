# Compute stopping boundary crossing probabilities.

Compute stopping boundary crossing probabilities.

## Usage

``` r
gs_probability_crt(theta = 0, I = 1, a = 0, b = 1, sides = 1, r = 18)
```

## Arguments

- theta:

  Effect size. Default value is 0.

- I:

  Information levels at interim analyses for which boundary crossing
  probabilities are computed.

- a:

  Lower futility boundaries for the interim analyses at the information
  levels specified by I.

- b:

  Upper efficacy boundaries for the interim analyses at the information
  levels specified by I.

- sides:

  `1=` 1-sided test (default)  
  `2=` 2-sided test.

- r:

  Integer value controlling grid for numerical integration as in
  Jennison and Turnbull (2000); default is 18, range is 1 to 80. Larger
  values provide larger number of grid points and greater accuracy.
  Normally `r` will not be changed by the user.

## Value

Object containing the following elements:

- k:

  Number of interim analyses.

- theta:

  As input.

- I:

  As input.

- lower:

  List containing the lower futility boundaries (`bound`) and the
  corresponding probabilities of crossing the lower boundaries given
  `theta` (`prob`).

- upper:

  List containing the upper efficacy boundaries (`bound`) and the
  corresponding probabilities of crossing the upper boundaries given
  theta (`prob`).

- power:

  Estimated power for trial.

- futile:

  Estimated probability of futility for trial

- r:

  As input.

## References

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

## Author

Lee Ding <lee_ding@g.harvard.edu>
