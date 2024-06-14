test_that("test that filter_pactr-class constructs properly", {
  mpactr_class <- mpactr$new(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                             test_path("exttestdata", "102623_metadata_correct.csv"))
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  expect_true(all(class(filter_class) == c("filter_pactr", "R6")))
  expect_equal(address(mpactr_class), address(filter_class$mpactr_data))
  expect_true(exists("list_of_summaries", filter_class$logger))
  expect_true(all(sapply(filter_class$logger[["list_of_summaries"]], is.null) == TRUE))
})

test_that("get_log returns an error when an incorrect fitler argument is provided", {
  mpactr_class <- mpactr$new(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                             test_path("exttestdata", "102623_metadata_correct.csv"))
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  
  expect_error(filter_class$get_log(filter = "cv"), "`filter` must be one of mpactR's")
})

test_that("get_log returns an error when the fitler argument provided has not yet been run (e.g., not in the log)", {
  mpactr_class <- mpactr$new(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                             test_path("exttestdata", "102623_metadata_correct.csv"))
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  
  expect_error(filter_class$get_log(filter = "mispicked"), "`filter` mispicked has not yet been applied to the data")
})

test_that("get_log returns the correct fitler summary list", {
  mpactr_class <- mpactr$new(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                             test_path("exttestdata", "102623_metadata_correct.csv"))
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
    TRUE)
  
  mispicked_summary <- filter_class$get_log(filter = "mispicked")
  
  expect_type(mispicked_summary, "list")
  expect_equal(length(mispicked_summary), 2)
  expect_equal(length(mispicked_summary$passed_ions), 1233)
})
