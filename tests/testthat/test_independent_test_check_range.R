# Test check_range function
#----------------------------------------------

test_that("Test check_range for invaild value for variable interval", code = {
  M <- c(3, 6, 2)
  expect_error(check_range(M, 1:10),
    info = "Test check_range for invaild value for variable interval"
  )
})


test_that("Test check_range valid value", code = {
  expect_invisible(check_range(1:5, interval = c(1, 5), inclusion = c(TRUE, TRUE)))
})
