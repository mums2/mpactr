## code to prepare `CulturesData` dataset goes here

CulturesData <- import_data(here::here("inst/extdata/cultures_peak_table.csv"),
                      here::here("inst/extdata/cultures_metadata.csv"))

usethis::use_data(CulturesData, overwrite = TRUE)
