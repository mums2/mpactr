test_that("test that filter_pactr-class constructs properly", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory, peak_table_name))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  expect_true(all(class(filter_class) == c("filter_pactr", "R6")))
  expect_equal(address(mpactr_class), address(filter_class$mpactr_data))
  expect_true(exists("list_of_summaries", filter_class$logger))
  expect_true(all(sapply(filter_class$logger[["list_of_summaries"]],
                         is.null) == TRUE))
})

test_that("get_log returns an error 
when an incorrect filter argument is provided", {

            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"

            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))

            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)

            expect_error(filter_class$get_log(filter = "cv"),
                         "`filter` must be one of mpactr's")
          })

test_that("get_log returns an error when the 
filter argument provided has not yet been run (e.g., not in the log)", {
            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"

            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))

            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)
            er <- "`filter` mispicked has not yet been applied to the data"
            expect_error(filter_class$get_log(filter = "mispicked"), er)
          })

test_that("get_log returns the correct filter summary list", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"

  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory, peak_table_name))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5,
    isowin = 0.01,
    trwin = 0.005,
    max_iso_shift = 3,
    merge_peaks = TRUE,
    merge_method = "sum"
  )

  mispicked_summary <- filter_class$get_log(filter = "mispicked")

  expect_type(mispicked_summary, "list")
  expect_equal(length(mispicked_summary), 2)
  expect_equal(length(mispicked_summary$passed_ions), 1233)
})

test_that("get_log returns the correct filter
 summary list when group is supplied", {
            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"
            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))

            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)

            filter_class$filter_blank()
            filter_class$parse_ions_by_group(group_threshold = 0.01)
            filter_class$apply_group_filter("Blanks", remove_ions = TRUE)

            group_summary <- filter_class$get_log(filter = "group",
                                                  group = "Blanks")

            expect_type(group_summary, "list")
          })


test_that("get_mispicked_ions returns error if 
check_mismatched_peaks has not been called", {
            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"
            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))
            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)

            expect_error(filter_class$get_mispicked_ions(),
                         "The mispicked filter has not yet been")
          })

test_that("get_mispicked_ions correctly returns
 the check_mismatched_peaks list", {
            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"
            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))

            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)
            filter_class$check_mismatched_peaks(
              ringwin = 0.5,
              isowin = 0.01,
              trwin = 0.005,
              max_iso_shift = 3,
              merge_peaks = TRUE,
              merge_method = "sum"
            )

            mispicked_groups <- filter_class$get_mispicked_ions()

            expect_equal(class(mispicked_groups), c("data.table", "data.frame"))
            expect_equal(length(mispicked_groups), 2)
            expect_equal(names(mispicked_groups), c("main_ion", "similar_ions"))
          })


test_that("get_group_averages calculates a group table", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory,
                                            peak_table_name))


  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)

  filter_class$check_mismatched_peaks(
    ringwin = 0.5,
    isowin = 0.01,
    trwin = 0.005,
    max_iso_shift = 3,
    merge_peaks = TRUE,
    merge_method = "sum"
  )

  avgs <- filter_class$get_group_averages()
  expect_equal(class(avgs), c("data.table", "data.frame"))
  expect_equal(nrow(avgs), (1233 * 6))

  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory,
                                            peak_table_name))



  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)

  filter_class$check_mismatched_peaks(
    ringwin = 0.5,
    isowin = 0.01,
    trwin = 0.005,
    max_iso_shift = 3,
    merge_peaks = TRUE,
    merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)

  avgs <- filter_class$get_group_averages()
  expect_equal(class(avgs), c("data.table", "data.frame"))
  expect_equal(nrow(avgs), (484 * 6))
})

test_that("get_cv returns the cv filter has been applied", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory,
                                            peak_table_name))


  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5,
    isowin = 0.01,
    trwin = 0.005,
    max_iso_shift = 3,
    merge_peaks = TRUE,
    merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)

  expect_error(filter_class$get_cv(), "The cv filter has not yet")

  filter_class$cv_filter(cv_threshold = 0.2)

  cv <- filter_class$get_cv()
  expect_equal(class(cv), c("data.table", "data.frame"))
})

test_that("Test that mpactr can be printed from the filter-pactR class", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory,
                                            peak_table_name))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  expect_output(filter_class$print())
})

test_that("is_filter_run correctly assesses if a filter has been run", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory,
                                            peak_table_name))


  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5,
    isowin = 0.01,
    trwin = 0.005,
    max_iso_shift = 3,
    merge_peaks = TRUE,
    merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)

  expect_true(filter_class$is_filter_run(filter = "group",
                                         group = "Blanks"))

  expect_false(filter_class$is_filter_run(filter = "insource"))
})

test_that("get_log returns an error when 
an incorrect filter argument is provided", {
            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"
            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))
            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)
            expect_error(filter_class$get_log(filter = "cv"),
                         "`filter` must be one of mpactr's")
          })

test_that("get_log returns an error when the filter
 argument provided has not yet been run (e.g., not in the log)", {
            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"
            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))
            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)
            err_msg <- "`filter` mispicked has not yet been applied to the data"
            expect_error(filter_class$get_log(filter = "mispicked"), err_msg)
          })

test_that("get_log returns the correct filter summary list", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory,
                                            peak_table_name))
  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(ringwin = 0.5,
                                      isowin = 0.01,
                                      trwin = 0.005,
                                      max_iso_shift = 3,
                                      merge_peaks = TRUE,
                                      merge_method = "sum")
  mispicked_summary <- filter_class$get_log(filter = "mispicked")
  expect_type(mispicked_summary, "list")
  expect_equal(length(mispicked_summary), 2)
  expect_equal(length(mispicked_summary$passed_ions), 1233)
})

test_that("get_mispicked_ions returns error if 
check_mismatched_peaks has not been called", {
            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"
            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))
            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)

            expect_error(filter_class$get_mispicked_ions(),
                         "The mispicked filter has not yet been")
          })

test_that("get_mispicked_ions correctly 
returns the check_mismatched_peaks list", {
            limit_cores()
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"
            meta <- fread(test_path(directory, meta_data_name))
            pt_list <- progenesis_formatter(test_path(directory,
                                                      peak_table_name))
            mpactr_class <- mpactr$new(
              pt_list,
              meta
            )
            mpactr_class$setup()
            filter_class <- filter_pactr$new(mpactr_class)
            filter_class$check_mismatched_peaks(ringwin = 0.5,
                                                isowin = 0.01,
                                                trwin = 0.005,
                                                max_iso_shift = 3,
                                                merge_peaks = TRUE,
                                                merge_method = "sum")

            mispicked_groups <- filter_class$get_mispicked_ions()

            expect_equal(class(mispicked_groups), c("data.table", "data.frame"))
            expect_equal(length(mispicked_groups), 2)
            expect_equal(names(mispicked_groups), c("main_ion", "similar_ions"))
          })


test_that("get_group_averages calculates a group table", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory,
                                            peak_table_name))
  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(ringwin = 0.5,
                                      isowin = 0.01,
                                      trwin = 0.005,
                                      max_iso_shift = 3,
                                      merge_peaks = TRUE,
                                      merge_method = "sum")
  avgs <- filter_class$get_group_averages()
  expect_equal(class(avgs), c("data.table", "data.frame"))
  expect_equal(nrow(avgs), (1233 * 6))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)

  filter_class$check_mismatched_peaks(ringwin = 0.5,
                                      isowin = 0.01,
                                      trwin = 0.005,
                                      max_iso_shift = 3,
                                      merge_peaks = TRUE,
                                      merge_method = "sum")
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)

  avgs <- filter_class$get_group_averages()
  expect_equal(class(avgs), c("data.table", "data.frame"))
  expect_equal(nrow(avgs), (484 * 6))
})

test_that("get_cv returns the cv filter has been applied", {
  limit_cores()
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta <- fread(test_path(directory, meta_data_name))
  pt_list <- progenesis_formatter(test_path(directory,
                                            peak_table_name))
  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(ringwin = 0.5,
                                      isowin = 0.01,
                                      trwin = 0.005,
                                      max_iso_shift = 3,
                                      merge_peaks = TRUE,
                                      merge_method = "sum")
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)

  expect_error(filter_class$get_cv(), "The cv filter has not yet")

  filter_class$cv_filter(cv_threshold = 0.2)

  cv <- filter_class$get_cv()
  expect_equal(class(cv), c("data.table", "data.frame"))
})
