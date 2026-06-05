# Hwang-Shih-DeCani Spending Function

The function `sf_hsd` implements a Hwang-Shih-DeCani spending function.

A Hwang-Shih-DeCani spending function takes the form \$\$f(t;\alpha,
\gamma)=\alpha(1-e^{-\gamma t})/(1-e^{-\gamma})\$\$ where \\\gamma\\ is
the value passed in `param`. A value of \\\gamma=-4\\ is used to
approximate an O'Brien-Fleming design (see
[`sf_exponential`](https://leejtding.github.io/gsDesignCRT/reference/sf_exponential.md)
for a better fit), while a value of \\\gamma=1\\ approximates a Pocock
design well.

## Usage

``` r
sf_hsd(alpha, t, param)
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

  A single real value specifying the gamma parameter for which
  Hwang-Shih-DeCani spending is to be computed; allowable range is
  \[-40, 40\]

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
