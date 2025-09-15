test_that("example path returns data paths properly", {
  limit_cores()
  expect_true(length(example_path()) > 0)
  expect_true(length(example_path("coculture_peak_table.csv")) >
                0)
  expect_error(example_path("a"))
})
