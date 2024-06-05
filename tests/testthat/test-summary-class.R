test_that("summary class constructor works correctly", {
  summary <- summary$new("filter1", c("1", "2", "3"), c("4", "5", "6"))

  expect_equal(summary$filter, "filter1")
  expect_equal(summary$failed_ions, c("1", "2", "3"))
  expect_equal(summary$passed_ions, c("4", "5", "6"))
})
