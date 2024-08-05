## code to prepare `cultures_data` dataset goes here

cultures_data <- import_data(
  here::here("inst/extdata/cultures_peak_table.csv"),
  here::here("inst/extdata/cultures_metadata.csv"),
  format = "Progenesis"
)

usethis::use_data(cultures_data, overwrite = TRUE)
