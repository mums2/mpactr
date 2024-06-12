test_that("mpactr class setup works properly", {
  mpactr_class <- mpactr$new(test_path("exttestdata","102623 peaktable coculture simple.csv"),
                             test_path("exttestdata", "102623_metadata_correct.csv"))
  mpactr_class$setup()

  expect_true(all(c("kmd", "mz", "rt") %in% colnames(mpactr_class$get_peak_table())))
  expect_equal(nrow(mpactr_class$get_peak_table()), 1303)
})
