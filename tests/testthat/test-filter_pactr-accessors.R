
test_that("filter_summary returns an error when an incorrect fitler argument is provided", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))

  expect_error(filter_summary(data, filter = "cv"), "`filter` must be one of mpactR's")
  
})

test_that("filter_summary returns an error when the fitler argument provided has not yet been run (e.g., not in the log)", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  expect_error(filter_summary(data_mpactr, filter = "replicability"), "`filter` replicability has not yet been applied to the data")
})

test_that(" returns the correct fitler summary list", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  mispicked_summary <- filter_summary(data_mpactr, filter = "mispicked")
  
  expect_type(mispicked_summary, "list")
  expect_equal(length(mispicked_summary), 2)
  expect_equal(length(mispicked_summary$passed_ions), 1233)
})

test_that("similar_ions returns error if check_mismatched_peaks has not been called", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  expect_error(similar_ions(data), "The mispicked filter has not yet been")
})

test_that("similar_ions correctly returns the check_mismatched_peaks list", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  mispicked_groups <- similar_ions(data_mpactr)
  
  expect_equal(class(mispicked_groups), c("data.table", "data.frame"))
  expect_equal(length(mispicked_groups), 2)
  expect_equal(names(mispicked_groups), c("main_ion", "similar_ions"))
})

test_that("group_averages calculates a group table if fitler blank hasn't been run", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  avgs <- group_averages(data_mpactr)
  expect_equal(class(avgs), c("data.table", "data.frame"))
})


test_that("cv_values returns the cv table if fitler cv has been run", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE) |>
    filter_group(group_to_remove = "Blanks") 
    
  data_mpactr$get_cv()
  data_mpactr$logger[["cv_values"]]
  
  expect_error(cv_values(data_mpactr), "The CV filter has not yet")
  
  data_mpactr |>
    filter_cv(cv_param = "median")
  
  cv <- cv_values(data_mpactr)

  expect_equal(class(cv), c("data.table", "data.frame"))
})