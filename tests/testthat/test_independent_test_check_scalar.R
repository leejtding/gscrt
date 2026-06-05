# Test check_scalar function
#----------------------------------------------

test_that("Test check_scalar for invalid value for is_type", code = {
  x <- 5
  expect_error(check_scalar(x, is_type = 12),
    info = "Test check_scalar for invalid range"
  )
})


test_that("Test check_scalar with valid value", code = {
  x <- 5
  expect_invisible(check_scalar(x, is_type = "numeric"))
})