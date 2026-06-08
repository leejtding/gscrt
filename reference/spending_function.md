# Spending Function

Spending Function

## Usage

``` r
spending_function(alpha, t, param)
```

## Arguments

- alpha:

  Real value \\\> 0\\ and no more than 1. Defaults in calls to
  [`gs_design_crt()`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md)
  are `alpha=0.025` for one-sided Type I error specification and
  `alpha=0.1` for Type II error specification. However, this could be
  set to 1 if, for descriptive purposes, you wish to see the proportion
  of spending as a function of the proportion of sample
  size/information.

- t:

  A vector of points with increasing values from 0 to 1, inclusive.
  Values of the proportion of sample size/information for which the
  spending function will be computed.

- param:

  A single real value or a vector of real values specifying the spending
  function parameter(s); this must be appropriately matched to the
  spending function specified.

## Value

`spending_function` and spending functions in general produce an object
of type `spendfn`.

- name:

  A character string with the name of the spending function.

- param:

  any parameters used for the spending function.

- parname:

  a character string or strings with the name(s) of the parameter(s) in
  `param`.

- sf:

  the spending function specified.

- spend:

  a vector of cumulative spending values corresponding to the input
  values in `t`.

- bound:

  this is null when returned from the spending function, but is set in
  [`gs_design_crt()`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md)
  if the spending function is called from there. Contains z-values for
  bounds of a design.

- prob:

  this is null when returned from the spending function, but is set in
  [`gs_design_crt()`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md)
  if the spending function is called from there. Contains probabilities
  of boundary crossing at `i`-th analysis for `j`-th theta value input
  to
  [`gs_design_crt()`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md)
  in `prob[i,j]`.

## Note

The gsDesign technical manual is available at
<https://keaven.github.io/gsd-tech-manual/>.

## References

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

## See also

[`gs_design_crt`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md),
[`sf_hsd`](https://leejtding.github.io/gscrt/reference/sf_hsd.md),
[`sf_power`](https://leejtding.github.io/gscrt/reference/sf_power.md),
[`sf_logistic`](https://leejtding.github.io/gscrt/reference/sf_distribution.md),
[`sf_exponential`](https://leejtding.github.io/gscrt/reference/sf_exponential.md),
[`sf_truncated`](https://leejtding.github.io/gscrt/reference/sf_special.md),
`vignette("gs_design_crt_package_overview")`

## Author

Keaven Anderson <keaven_anderson@merck.com>
