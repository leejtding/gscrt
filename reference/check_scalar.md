# Utility functions to verify variable properties

Utility functions to verify an objects's properties including whether it
is a scalar or vector, the class, the length, and (if numeric) whether
the range of values is on a specified interval. Additionally, the
`check_lengths` function can be used to ensure that all the supplied
inputs have equal lengths.

`is_integer` is similar to
[`is.integer`](https://rdrr.io/r/base/integer.html) except that
`is_integer(1)` returns `TRUE` whereas `is.integer(1)` returns `FALSE`.

`check_scalar` is used to verify that the input object is a scalar as
well as the other properties specified above.

`check_vector` is used to verify that the input object is an atomic
vector as well as the other properties as defined above.

`check_range` is used to check whether the numeric input object's values
reside on the specified interval. If any of the values are outside the
specified interval, a `FALSE` is returned.

`check_length` is used to check whether all of the supplied inputs have
equal lengths.

## Usage

``` r
check_lengths(..., allow_single = FALSE)

check_range(
  x,
  interval = 0:1,
  inclusion = c(TRUE, TRUE),
  varname = deparse(substitute(x)),
  tol = 0
)

check_scalar(x, is_type = "numeric", ...)

check_vector(x, is_type = "numeric", ..., length = NULL)

is_integer(x)
```

## Arguments

- ...:

  For the `check_scalar` and `check_vector` functions, this input
  represents additional arguments sent directly to the `check_range`
  function. For the

  `check_lengths` function, this input represents the arguments to check
  for equal lengths.

- allow_single:

  logical flag. If `TRUE`, arguments that are vectors comprised of a
  single element are not included in the comparative length test in the
  `check_lengths` function. Partial matching on the name of this
  argument is not performed so you must specify 'allow_single' in its
  entirety in the call.

- x:

  any object.

- interval:

  two-element numeric vector defining the interval over which the input
  object is expected to be contained. Use the `inclusion` argument to
  define the boundary behavior.

- inclusion:

  two-element logical vector defining the boundary behavior of the
  specified interval. A `TRUE` value denotes inclusion of the
  corresponding boundary. For example, if `interval=c(3,6)` and
  `inclusion=c(FALSE,TRUE)`, then all the values of the input object are
  verified to be on the interval (3,6\].

- varname:

  character string defining the name of the input variable as sent into
  the function by the caller. This is used primarily as a mechanism to
  specify the name of the variable being tested when `check_range` is
  being called within a function.

- tol:

  numeric scalar defining the tolerance to use in testing the intervals
  of the

  `check_range` function.

- is_type:

  character string defining the class that the input object is expected
  to be.

- length:

  integer specifying the expected length of the object in the case it is
  a vector. If `length=NULL`, the default, then no length check is
  performed.

## Value

`is_integer`: Boolean value as checking result Other functions have no
return value, called for side effects

## Examples

``` r

# check whether input is an integer
is_integer(1)
#> [1] TRUE
is_integer(1:5)
#> [1] TRUE
try(is_integer("abc")) # expect error
#> [1] FALSE

# check whether input is an integer scalar
check_scalar(3, "integer")

# check whether input is an integer scalar that resides
# on the interval on [3, 6]. Then test for interval (3, 6].
check_scalar(3, "integer", c(3, 6))
try(check_scalar(3, "integer", c(3, 6), c(FALSE, TRUE))) # expect error
#> Error in check_range(x, ..., varname = deparse(substitute(x))) : 
#>   3 not on interval (3, 6]

# check whether the input is an atomic vector of class numeric,
# of length 3, and whose value all reside on the interval [1, 10)
x <- c(3, pi, exp(1))
check_vector(x, "numeric", c(1, 10), c(TRUE, FALSE), length = 3)

# do the same but change the expected length; expect error
try(check_vector(x, "numeric", c(1, 10), c(TRUE, FALSE), length = 2))
#> Error in check_vector(x, "numeric", c(1, 10), c(TRUE, FALSE), length = 2) : 
#>   object 'varstr' not found

# create faux function to check input variable
foo <- function(moo) check_vector(moo, "character")
foo(letters)
try(foo(1:5)) # expect error with function and argument name in message
#> Error in check_vector(moo, "character") : 
#>   In function foo : variable moo  must be vector of class  character

# check for equal lengths of various inputs
check_lengths(1:2, 2:3, 3:4)
try(check_lengths(1, 2, 3, 4:5)) # expect error
#> Error in check_lengths(1, 2, 3, 4:5) : 
#>   In function doTryCatch :lengths of inputs are not all equal

# check for equal length inputs but ignore single element vectors
check_lengths(1, 2, 3, 4:5, 7:8, allow_single = TRUE)

```
