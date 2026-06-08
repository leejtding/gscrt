# t-distribution Spending Function

The function `sf_t_dist()` provides perhaps the maximum flexibility
among spending functions provided in the `gs_design_crt` package. This
function allows fitting of three points on a cumulative spending
function curve; in this case, six parameters are specified indicating an
x and a y coordinate for each of 3 points.

The t-distribution spending function takes the form
\$\$f(t;\alpha)=\alpha F(a+b_f^{-1}(t))\$\$ where \\F()\\ is a
cumulative t-distribution function with `df` degrees of freedom and
\\F^{-1}()\\ is its inverse.

## Usage

``` r
sf_t_dist(alpha, t, param)
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

  In the three-parameter specification, the first paramater (a) may be
  any real value, the second (b) any positive value, and the third
  parameter (df=degrees of freedom) any real value 1 or greater. When
  [`gs_design_crt()`](https://leejtding.github.io/gscrt/reference/gs_design_crt.md)
  is called with a t-distribution spending function, this is the
  parameterization printed. The five parameter specification is
  `c(t1,t2,u1,u2,df)` where the objective is that the resulting
  cumulative proportion of spending at `t` represented by `sf(t)`
  satisfies `sf(t1)=alpha*u1`, `sf(t2)=alpha*u2`. The t-distribution
  used has `df` degrees of freedom. In this parameterization, all the
  first four values must be between 0 and 1 and `t1 < t2`, `u1 < u2`.
  The final parameter is any real value of 1 or more. This
  parameterization can fit any two points satisfying these requirements.
  The six parameter specification attempts to fit 3 points, but does not
  have flexibility to fit any three points. In this case, the
  specification for `param` is c(t1,t2,t3,u1,u2,u3) where the objective
  is that `sf(t1)=alpha*u1`, `sf(t2)=alpha*u2`, and `sf(t3)=alpha*u3`.
  See examples to see what happens when points are specified that cannot
  be fit.

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
