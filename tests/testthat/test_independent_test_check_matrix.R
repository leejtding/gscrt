# --------------------------------------------
# Test check_matrix function
#----------------------------------------------

test_that("Test check_matrix with invalid input", code = {
  M <- 5
  expect_error(check_matrix(M, nrows = 3, ncols = 3),
    info = "Test check_matrix with invalid input"
  )
})


test_that("Test check_matrix function for invalid nrows", code = {
  M <- matrix(c(32, 42, 54, 16, 7, 87, 49, 10, 11, 22, 24, 45),
    nrow = 4,
    byrow = TRUE
  )
  expect_error(check_rue_matrix(M, nrows = 3, ncols = 3),
    info = "Test check_matrix function for invalid nrows"
  )
})

test_that("Test check_matrix function for invalid ncols", code = {
  M <- matrix(c(32, 42, 54, 16, 7, 87, 49, 10, 11, 22, 24, 45),
    nrow = 4,
    byrow = TRUE
  )
  expect_error(check_matrix(M, nrows = 4, ncols = 2),
    info = "Test check_matrix function for invalid ncols"
  )
})


test_that("Test check_matrix function for invalid is_type", code = {
  M <- matrix(c(32, 42, 54, 16, 7, 87, 49, 10, 11, 22, 24, 45),
    nrow = 4,
    byrow = TRUE
  )
  expect_error(check_matrix(M, is_type = "character"),
    info = "Test check_matrix function for invalid is_type"
  )
})

test_that("Test check_matrix function for invalid is_type", code = {
  M <- matrix(c("a", "b", "cc", "aa", "cc", "dd"), nrow = 2, byrow = TRUE)
  expect_error(check_matrix(M, is_type = "numeric"),
    info = "Test check_matrix function for invalid is_type"
  )
})

test_that("Test check_matrix function for mixed data ", code = {
  M <- matrix(c(32, 42, 54, 16, 7, "a"), nrow = 2, byrow = TRUE)
  expect_error(check_matrix(M, is_type = "numeric"),
    info = "Test check_matrix function for mixed data"
  )
})

test_that("Test check_matrix function for data in valid interval",
  code = {
    M <- matrix(c(32, 42, 54, 16, 7, 34), nrow = 2, byrow = TRUE)
    expect_error(check_matrix(M, interval = c(21, 45), is_type = "numeric"),
      info = "Test check_matrix function for data in valid interval"
    )
  }
)


test_that("Test check_matrix function for integer data",
  code = {
    M <- matrix(c(32L, 42.5, 54L, 16L), nrow = 2, byrow = TRUE)
    expect_error(check_matrix(M, is_type = "integer"),
      info = "Test check_matrix function for integer data"
    )
  }
)


test_that("Test check_matrix function for valid data", code = {
  M <- matrix(c(32, 42, 54, 16, 7, 87, 49, 10, 11, 22, 24, 45),
    nrow = 4,
    byrow = TRUE
  )
  expect_invisible(check_matrix(M, nrows = 4, ncols = 3))
})