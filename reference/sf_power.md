# Kim-DeMets (power) Spending Function

The function `sf_power()` implements a Kim-DeMets (power) spending
function. This is a flexible, one-parameter spending function
recommended by Jennison and Turnbull (2000).

A Kim-DeMets spending function takes the form
\$\$f(t;\alpha,\rho)=\alpha t^\rho\$\$ where \\\rho\\ is the value
passed in `param`. See examples below for a range of values of \\\rho\\
that may be of interest (`param=0.75` to `3` are documented there).

## Usage

``` r
sf_power(alpha, t, param)
```

## Arguments

- alpha:

  Real value \\\> 0\\ and no more than 1. Normally, `alpha=0.025` for
  one-sided Type I error specification or `alpha=0.1` for Type II error
  specification. However, this could be set to 1 if for descriptive
  purposes you wish to see the proportion of spending as a function of
  the proportion of sample size/information.

- t:

  A vector of points with increasing values from 0 to 1, inclusive.
  Values of the proportion of sample size/information for which the
  spending function will be computed.

- param:

  A single, positive value specifying the \\\rho\\ parameter for which
  Kim-DeMets spending is to be computed; allowable range is (0,50\]

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
`vignette("gs_design_crt_package_overview")`

## Author

Keaven Anderson <keaven_anderson@merck.com>
