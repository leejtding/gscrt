# Boundary derivation for efficacy stopping only.

`gs_upper_crt` is used to calculate the stopping boundaries for a group
sequential trial with 1 or 2-sided efficacy-only stopping. Code adapted
from gsDesign package.

## Usage

``` r
gs_upper_crt(
  theta = 0,
  I,
  a = NULL,
  false_pos,
  sides = 1,
  tol = 1e-06,
  r = 18,
  printerr = 0
)
```

## Arguments

- theta:

  Effect size under null hypothesis. Default value is 0.

- I:

  Information levels at interim analyses for which stopping boundaries
  are calculated. At any given interim analysis this should include the
  information levels at the current and previous analyses.

- a:

  Lower futility boundary values for the interim analyses at the
  information levels specified by I.

- false_pos:

  Type I error spent for the interim analyses at the information levels
  specified by I.

- sides:

  `1=` 1-sided test (default)  
  `2=` 2-sided test

- tol:

  Tolerance for error (default is 0.000001). Normally this will not be
  changed by the user. This does not translate directly to number of
  digits of accuracy, so use extra decimal places.

- r:

  Integer value controlling grid for numerical integration as in
  Jennison and Turnbull (2000); default is 18, range is 1 to 80. Larger
  values provide larger number of grid points and greater accuracy.
  Normally `r` will not be changed by the user.

- printerr:

  Print output for debugging.

## Value

Object containing the following elements:

- k:

  Number of interim analyses.

- theta:

  As input.

- I:

  As input.

- a:

  As input.

- b:

  Computed upper efficacy boundaries at the specified interim analyses.

- r:

  As input.

- error:

  Error flag returned; 0 if convergence; 1 indicates error.

## References

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

## Author

Lee Ding <lee_ding@g.harvard.edu>
