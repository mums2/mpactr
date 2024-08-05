test_that("import_data creates a proper mpactr and filter-pactr object", {

  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"

  data <- import_data(
    peak_table = test_path(directory, peak_table_name),
    meta_data = test_path(directory,  meta_data_name),
    format = "Progenesis"
  )

  expect_true(all(class(data) == c("filter_pactr", "R6")))
  expect_true(exists("list_of_summaries", data$logger))
  expect_true(all(sapply(data$logger[["list_of_summaries"]], is.null) == TRUE))
  expect_true(nrow(data$mpactr_data$get_peak_table()) > 1)
  expect_true(nrow(data$mpactr_data$get_meta_data()) > 1)
})


test_that("import_data aborts when expected
 metadata columns are not provided", {
            directory <- "exttestdata"
            peak_table_name <- "102623_peaktable_coculture_simple.csv"
            meta_data_name <- "102623_metadata_correct.csv"
            meta_data_abort <- read_csv(test_path(directory,  meta_data_name),
                                        show_col_types = FALSE)
            colnames(meta_data_abort) <- c("Injection",
                                           "Sample", "Biological_Group")

            expect_error(import_data(
              peak_table = test_path(directory, peak_table_name),
              meta_data = meta_data_abort,
              format = "Progenesis"
            ))
          })
