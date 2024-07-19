
test_that("filter_summary returns an error when an incorrect fitler argument is provided", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"), format = "Progenesis")

  expect_error(filter_summary(data, filter = "cv"), "`filter` must be one of mpactR's")
  
})

test_that("filter_summary returns an error when the fitler argument provided has not yet been run (e.g., not in the log)", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"), format = "Progenesis")
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  expect_error(filter_summary(data_mpactr, filter = "replicability"), "`filter` replicability has not yet been applied to the data")
})

test_that(" returns the correct fitler summary list", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"), format = "Progenesis")
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  mispicked_summary <- filter_summary(data_mpactr, filter = "mispicked")
  
  expect_type(mispicked_summary, "list")
  expect_equal(length(mispicked_summary), 2)
  expect_equal(length(mispicked_summary$passed_ions), 1233)
})

test_that("get_similar_ions returns error if check_mismatched_peaks has not been called", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"), format = "Progenesis")
  
  expect_error(get_similar_ions(data), "The mispicked filter has not yet been")
})

test_that("get_similar_ions correctly returns the check_mismatched_peaks list", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"), format = "Progenesis")
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  mispicked_groups <- get_similar_ions(data_mpactr)
  
  expect_equal(class(mispicked_groups), c("data.table", "data.frame"))
  expect_equal(length(mispicked_groups), 2)
  expect_equal(names(mispicked_groups), c("main_ion", "similar_ions"))
})

test_that("get_group_averages calculates a group table", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"), format = "Progenesis")
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  avgs <- get_group_averages(data_mpactr)
  expect_equal(class(avgs), c("data.table", "data.frame"))
  expect_equal(nrow(avgs), (1233 * 6))
})


test_that("get_cv_data returns the cv table if fitler cv has been run", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"), format = "Progenesis")
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE) |>
    filter_group(group_to_remove = "Blanks")
  
  expect_error(get_cv_data(data_mpactr), "The cv filter has not yet")
  
  data_mpactr |>
    filter_cv(cv_threshold = 0.2, cv_param = "median")
  
  cv <- get_cv_data(data_mpactr)

  expect_equal(class(cv), c("data.table", "data.frame"))
})