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
