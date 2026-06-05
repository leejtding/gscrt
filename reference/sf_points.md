# Pointwise Spending Function

The function `sf_points` implements a spending function with values
specified for an arbitrary set of specified points. It is now
recommended to use sf_linear rather than sf_points. If using
`sf_points()` in a design, it is recommended to specify how to
interpolate between the specified points (e.g,, linear interpolation);
also consider fitting smooth spending functions; see
[`vignette("spending_function_overview")`](https://leejtding.github.io/gsDesignCRT/articles/spending_function_overview.md).

## Usage

``` r
sf_points(alpha, t, param)
```

## Arguments

- alpha:

  Real value \\\> 0\\ and no more than 1. Normally, `alpha=0.025` for
  one-sided Type I error specification or `alpha=0.1` for Type II error
  specification. However, this could be set to 1 if for descriptive
  purposes you wish to see the proportion of spending as a function of
  the proportion of sample size/information.

- t:

  A vector of points with increasing values from \>0 and \<=1. Values of
  the proportion of sample size/information for which the spending
  function will be computed.

- param:

  A vector of the same length as `t` specifying the cumulative
  proportion of spending to corresponding to each point in `t`; must be
  \>=0 and \<=1.

## Value

An object of type `spendfn`. See
[`vignette("spending_function_overview")`](https://leejtding.github.io/gsDesignCRT/articles/spending_function_overview.md)
for further details.

## Note

The gsDesign technical manual is available at
<https://keaven.github.io/gsd-tech-manual/>.

## References

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

## See also

[`vignette("spending_function_overview")`](https://leejtding.github.io/gsDesignCRT/articles/spending_function_overview.md),
[`gs_design_crt`](https://leejtding.github.io/gsDesignCRT/reference/gs_design_crt.md),
`vignette("gs_design_crt_package_overview")`,
[sf_logistic](https://leejtding.github.io/gsDesignCRT/reference/sf_distribution.md)

## Author

Keaven Anderson <keaven_anderson@merck.com>
