#' Import data into an mpactr object.
#'
#' @description 
#' `import_data()` takes two file paths, one for the pre-processed feature table and one for sample metadata. Both files should be .csv.
#'
#' @details
#' Expected feature table format.
#'
#' Expected metadata table format.
#'
#' @param peak_table The file path to your feature table file or a or `data.frame`.
#' @param meta_data The file path to your meta_data file or `data.frame`. 
#'
#' @return an `mpactr_object`
#' @export 
#'
#' @examples 
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
#' peak_table <- read.csv(example("coculture_peak_table.csv"), skip = 2, check.names = FALSE)
#' meta_data <- read.csv(example("metadata.csv"))
#' data <- import_data(peak_table,
#'                       meta_data)
#'
import_data <- function(peak_table, meta_data)
{
    mpactr_object <- mpactr$new(peak_table_path = peak_table,
                                meta_data_path = meta_data)
    mpactr_object$setup()
    filter_object <- filter_pactr$new(mpactr_object)
    return(filter_object)
}