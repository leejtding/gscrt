# Exponential Spending Function

The function `sf_exponential` implements the exponential spending
function (Anderson and Clark, 2009).

An exponential spending function is defined for any positive `nu` and
\\0\le t\le 1\\ as \$\$f(t;\alpha,\nu)=\alpha(t)=\alpha^{t^{-\nu}}.\$\$
A value of `nu=0.8` approximates an O'Brien-Fleming spending function
well.

The general class of spending functions this family is derived from
requires a continuously increasing cumulative distribution function
defined for \\x\>0\\ and is defined as \$\$f(t;\alpha,
\nu)=1-F\left(F^{-1}(1-\alpha)/ t^\nu\right).\$\$ The exponential
spending function can be derived by letting \\F(x)=1-\exp(-x)\\, the
exponential cumulative distribution function. This function was derived
as a generalization of the Lan-DeMets (1983) spending function used to
approximate an O'Brien-Fleming spending function
([`sf_ldof()`](https://leejtding.github.io/gsDesignCRT/reference/sf_ldof.md)),
\$\$f(t; \alpha)=2-2\Phi \left( \Phi^{-1}(1-\alpha/2)/ t^{1/2}
\right).\$\$

## Usage

``` r
sf_exponential(alpha, t, param)
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

  A single positive value specifying the nu parameter for which the
  exponential spending is to be computed; allowable range is (0, 1.5\].

## Value

An object of type `spendfn`. See
[`vignette("spending_function_overview")`](https://leejtding.github.io/gsDesignCRT/articles/spending_function_overview.md)
for further details.

## Note

The gsDesign technical manual shows how to use `sf_exponential()` to
closely approximate an O'Brien-Fleming design. The manual is available
at \<https://keaven.github.io/gsd-tech-manual/\>.

## References

Anderson KM and Clark JB (2009), Fitting spending functions. *Statistics
in Medicine*; 29:321-327.

Jennison C and Turnbull BW (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

Lan, KKG and DeMets, DL (1983), Discrete sequential boundaries for
clinical trials. *Biometrika*; 70:659-663.

## See also

[`vignette("spending_function_overview")`](https://leejtding.github.io/gsDesignCRT/articles/spending_function_overview.md),
[`gs_design_crt`](https://leejtding.github.io/gsDesignCRT/reference/gs_design_crt.md),
`vignette("gs_design_crt_package_overview")`

## Author

Keaven Anderson <keaven_anderson@merck.com>
