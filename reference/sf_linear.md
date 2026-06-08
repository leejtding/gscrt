# Piecewise Linear and Step Function Spending Functions

The function `sf_linear()` allows specification of a piecewise linear
spending function. The function `sf_step()` specifies a step function
spending function. Both functions provide complete flexibility in
setting spending at desired timepoints in a group sequential design.
Normally these function will be passed to
[`gs_design_crt()`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md)
in the parameter `sfu` for the upper bound or `sfl` for the lower bound
to specify a spending function family for a design. When passed to
[`gs_design_crt()`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md),
the value of `param` would be passed to `sf_linear()` or `sf_step()`
through the
[`gs_design_crt()`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md)
arguments `sfupar` for the upper bound and `sflpar` for the lower bound.

Note that `sf_step()` allows setting a particular level of spending when
the timing is not strictly known; an example shows how this can inflate
Type I error when timing of analyses are changed based on knowing the
treatment effect at an interim.

## Usage

``` r
sf_linear(alpha, t, param)

sf_step(alpha, t, param)
```

## Arguments

- alpha:

  Real value \\\> 0\\ and no more than 1. Normally, `alpha=0.025` for
  one-sided Type I error specification or `alpha=0.1` for Type II error
  specification. However, this could be set to 1 if for descriptive
  purposes you wish to see the proportion of spending as a function of
  the proportion of sample size or information.

- t:

  A vector of points with increasing values from 0 to 1, inclusive.
  Values of the proportion of sample size or information for which the
  spending function will be computed.

- param:

  A vector with a positive, even length. Values must range from 0 to 1,
  inclusive. Letting `m <- length(param/2)`, the first `m` points in
  param specify increasing values strictly between 0 and 1 corresponding
  to interim timing (proportion of final total statistical information).
  The last `m` points in `param` specify non-decreasing values from 0 to
  1, inclusive, with the cumulative proportion of spending at the
  specified timepoints.

## Value

An object of type `spendfn`. The cumulative spending returned in
`sf_linear$spend` is 0 for `t <= 0` and `alpha` for `t>=1`. For `t`
between specified points, linear interpolation is used to determine
`sf_linear$spend`.

The cumulative spending returned in `sf_step$spend` is 0 for
`t<param[1]` and `alpha` for `t>=1`. Letting `m <- length(param/2)`, for
`i=1,2,...m-1` and ` param[i]<= t < param[i+1]`, the cumulative spending
is set at `alpha * param[i+m]` (also for `param[m]<=t<1`).

Note that if `param[2m]` is 1, then the first time an analysis is
performed after the last proportion of final planned information
(`param[m]`) will be the final analysis, using any remaining error that
was not previously spent.

See
[`vignette("spending_function_overview")`](https://leejtding.github.io/gscrt/articles/spending_function_overview.md)
for further details.

## Note

The gsDesign technical manual is available at
<https://keaven.github.io/gsd-tech-manual/>.

## References

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

## See also

[`vignette("spending_function_overview")`](https://leejtding.github.io/gscrt/articles/spending_function_overview.md),
[`gs_design_crt`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md),
`vignette("gs_design_crt_package_overview")`

## Author

Keaven Anderson <keaven_anderson@merck.com>
