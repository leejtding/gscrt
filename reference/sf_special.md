# Truncated, trimmed and gapped spending functions

The functions `sf_truncated()` and `sf_trimmed` apply any other spending
function over a restricted range. This allows eliminating spending for
early interim analyses when you desire not to stop for the bound being
specified; this is usually applied to eliminate early tests for a
positive efficacy finding. The truncation can come late in the trial if
you desire to stop a trial any time after, say, 90 percent of
information is available and an analysis is performed. This allows full
Type I error spending if the final analysis occurs early. Both functions
set cumulative spending to 0 below a 'spending interval' in the interval
\[0,1\], and set cumulative spending to 1 above this range.
`sf_trimmed()` otherwise does not change an input spending function that
is specified; probably the preferred and more intuitive method in most
cases. `sf_truncated()` resets the time scale on which the input
spending function is computed to the 'spending interval.'

`sf_gapped()` allows elimination of analyses after some time point in
the trial; see details and examples.

`sf_trimmed` simply computes the value of the input spending function
and parameters in the sub-range of \[0,1\], sets spending to 0 below
this range and sets spending to 1 above this range.

`sf_gapped` spends outside of the range provided in trange. Below
trange, the input spending function is used. Above trange, full spending
is used; i.e., the first analysis performed above the interval in trange
is the final analysis. As long as the input spending function is
strictly increasing, this means that the first interim in the interval
trange is the final interim analysis for the bound being specified.

`sf_truncated` compresses spending into a sub-range of \[0,1\]. The
parameter `param$trange` specifies the range over which spending is to
occur. Within this range, spending is spent according to the spending
function specified in `param$sf` along with the corresponding spending
function parameter(s) in `param$param`. See example using `sf_linear`
that spends uniformly over specified range.

## Usage

``` r
sf_truncated(alpha, t, param)

sf_trimmed(alpha, t, param)

sf_gapped(alpha, t, param)
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

  a list containing the elements sf (a spendfn object such as sf_hsd),
  trange (the range over which the spending function increases from 0 to
  1; 0 \<= trange\[1\]\<trange\[2\] \<=1; for sf_gapped, trange\[1\]
  must be \> 0), and param (null for a spending function with no
  parameters or a scalar or vector of parameters needed to fully specify
  the spending function in sf).

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
