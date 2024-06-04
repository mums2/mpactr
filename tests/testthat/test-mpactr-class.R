test_that("mpactr class initialize works correctly", {
  mpactr_class <- mpactr$new(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"),
  here::here("tests/exttestdata/102623_metadata_correct.csv"))
  
  expect_true(all(class(mpactr_class) == c("mpactr", "R6")))
  expect_equal(length(mpactr_class$meta_data_csv), 3)
  expect_equal(length(mpactr_class$peak_table_csv), 21)
  expect_error(mpactr$new(peak_table_path = 2, meta_data_path = 5), NULL)
})

test_that("mpactr class setup works properly", {
  mpactr_class <- mpactr$new(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"),
      here::here("tests/exttestdata/102623_metadata_correct.csv"))
  mpactr_class$setup()
  
  expect_true("kmd" %in% colnames(mpactr_class$peak_table_csv))
})


