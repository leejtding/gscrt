# Two-parameter Spending Function Families

The functions `sf_logistic()`, `sf_normal()`, `sf_extreme_value()`,
`sf_extreme_value2()`, `sf_cauchy()`, and `sf_beta_dist()` are all
2-parameter spending function families. These provide increased
flexibility in some situations where the flexibility of a one-parameter
spending function family is not sufficient. These functions all allow
fitting of two points on a cumulative spending function curve; in this
case, four parameters are specified indicating an x and a y coordinate
for each of 2 points.

`sf_beta_dist(alpha,t,param)` is simply `alpha` times the incomplete
beta cumulative distribution function with parameters \\a\\ and \\b\\
passed in `param` evaluated at values passed in `t`.

The other spending functions take the form \$\$f(t;\alpha,a,b)=\alpha
F(a+b_f^{-1}(t))\$\$ where \\F()\\ is a cumulative distribution function
with values \\\> 0\\ on the real line (logistic for `sf_logistic()`,
normal for `sf_normal()`, extreme value for `sf_extreme_value()` and
Cauchy for `sf_cauchy()`) and \\F^{-1}()\\ is its inverse.

For the logistic spending function this simplifies to
\$\$f(t;\alpha,a,b)=\alpha (1-(1+e^a(t/(1-t))^b)^{-1}).\$\$

For the extreme value distribution with \$\$F(x)=\exp(-\exp(-x))\$\$
this simplifies to \$\$f(t;\alpha,a,b)=\alpha \exp(-e^a (-\ln t)^b).\$\$
Since the extreme value distribution is not symmetric, there is also a
version where the standard distribution is flipped about 0. This is
reflected in `sf_extreme_value2()` where \$\$F(x)=1-\exp(-\exp(x)).\$\$

## Usage

``` r
sf_logistic(alpha, t, param)

sf_beta_dist(alpha, t, param)

sf_cauchy(alpha, t, param)

sf_extreme_value(alpha, t, param)

sf_extreme_value2(alpha, t, param)

sf_normal(alpha, t, param)
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

  In the two-parameter specification, `sf_beta_dist()` requires 2
  positive values, while `sf_logistic()`, `sf_normal()`,
  `sf_extreme_value()`,

  `sf_extreme_value2()` and `sf_cauchy()` require the first parameter to
  be any real value and the second to be a positive value. The four
  parameter specification is `c(t1,t2,u1,u2)` where the objective is
  that `sf(t1)=alpha*u1` and `sf(t2)=alpha*u2`. In this
  parameterization, all four values must be between 0 and 1 and
  `t1 < t2`, `u1 < u2`.

## Value

An object of type `spendfn`. See
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
