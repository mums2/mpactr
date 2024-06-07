test_that("summary class constructor works correctly", {
  summary <- summary$new("filter1", c(1, 2, 3), c(4, 5, 6))

  expect_equal(summary$get_filter(), "filter1")
  expect_equal(summary$get_failed_ions(), c(1, 2, 3))
  expect_equal(summary$get_passed_ions(), c(4, 5, 6))
})
