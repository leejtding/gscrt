# Test check_vector function
#----------------------------------------------

test_that("Test check_vector for invalid value for length", code = {
  x <- c(1, 5, 2, 3)
  expect_error(check_vector(x, length = "e"),
    info = "Test check_vector for invalid value for length"
  )
})


test_that("Test check_vector length not matching vector length", code = {
  x <- c(1, 5, 2, 3)
  expect_error(check_vector(x, length = 2),
    info = "Test check_vector length not matching vector length"
  )
})


test_that("Test check_vector with correct value of the length parameter ", code = {
  x <- c(1, 5, 2, 3)
  expect_invisible(check_vector(x, length = 4))
})