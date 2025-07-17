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

test_that("We can use a data.frame as input our peak table in import_data", {
  directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta_data_path <- test_path(directory,  meta_data_name)
  data <- import_data(
    peak_table = test_path(directory, peak_table_name),
    meta_data = meta_data_path,
    format = "Progenesis"
  )
  peak_table <- get_peak_table(data)

  data_df <- import_data(get_peak_table(data), get_meta_data(data), "none")
  peak_table_df <- get_peak_table(data_df)

  expect_true(all(peak_table == peak_table_df))

  data_df_meta_data_file <- import_data(get_peak_table(data),
                                        meta_data_path, "none")
  peak_table_meta_data_file <- get_peak_table(data_df_meta_data_file)
  expect_true(all(peak_table == peak_table_meta_data_file))
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





test_that("import_data aborts when a data.frame is given as input without a proper metadata file", {
  #Pre-filtered
 directory <- "exttestdata"
  peak_table_name <- "102623_peaktable_coculture_simple.csv"
  meta_data_name <- "102623_metadata_correct.csv"
  meta_data_path <- test_path(directory,  meta_data_name)
  data <- import_data(
    peak_table = test_path(directory, peak_table_name),
    meta_data = meta_data_path,
    format = "Progenesis"
  )
  meta <- get_meta_data(data)[1:15, ]
  expect_error(import_data(get_peak_table(data), meta, "none"), "These columns are not")

  meta <- get_meta_data(data)
  meta <- rbind(meta, data.frame(Injection = "UM143", Sample_Code = "143", Biological_Group = "Blanks"))
  expect_error(import_data(get_peak_table(data), meta, "none"), "You have extra columns")

  # Post-filtered
  data <- import_data(
    peak_table = test_path(directory, peak_table_name),
    meta_data = meta_data_path,
    format = "Progenesis"
  )|> 
    filter_mispicked_ions() |>
    filter_cv(0.1) |>
    filter_group(0.1, "Blanks") |>
    filter_insource_ions()

  meta <- get_meta_data(data)[1:15, ]
  expect_error(import_data(get_peak_table(data), meta, "none"), "These columns are not")

  meta <- get_meta_data(data)
  meta <- rbind(meta, data.frame(Injection = "UM143", Sample_Code = "143", Biological_Group = "Blanks"))
  expect_error(import_data(get_peak_table(data), meta, "none"), "You have extra columns")
})


test_that("unique_compounds annotate duplicates properly", {
  df <- data.frame(Compound = c(1, 2, 3, 1:3, 4:7))
  ls <- list(peak_table = df, raw_table = df)
  uniqued_list <- unique_compounds(ls)

  expect_true(length(unique(uniqued_list$peak_table$Compound)) ==
                length(uniqued_list$peak_table$Compound))
  expect_true(length(unique(uniqued_list$raw_table$Compound)) ==
                length(uniqued_list$raw_table$Compound))

  expect_false(length(unique(df$Compound)) == length(df$Compound))

  df <- data.frame(Compound = c("1", "1", "1_1", "1_1_1"))
  ls <- list(peak_table = df, raw_table = df)
  uniqued_list <- unique_compounds(list(peak_table = df, raw_table = df))

  expect_true(uniqued_list$peak_table$Compound[[3]] == "1_1_1")
  expect_true(uniqued_list$peak_table$Compound[[4]] == "1_1_1_1")

  df <- data.frame(Compound = c("1", "1", "1_1", "1_1_1"))
  ls <- list(peak_table = df, raw_table = df)
  expect_message(unique_compounds(ls))

  df <- data.frame(Compound = c("1", "1", NA, "1_1_1"))
  ls <- list(peak_table = df, raw_table = df)
  expect_error(unique_compounds(ls))
})
