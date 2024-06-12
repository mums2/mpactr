## code to prepare `CulturesData` dataset goes here

CulturesData <- import_data(here::here("../tests/exttestdata/102623_peaktable_coculture_simple.csv"),
                      here::here("tests/exttestdata/102623_metadata_correct.csv"))

usethis::use_data(CulturesData, overwrite = TRUE)
