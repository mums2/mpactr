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
#' @param peak_table_file_path The file path to your feature table file.
#' @param meta_data_file_path The file path to your meta_data file. 
#'
#' @return an `mpactr_object`
#' @export 
#'
#' @examples 
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
import_data <- function(peak_table_file_path, meta_data_file_path)
{
    mpactr_object <- mpactr$new(peak_table_path = peak_table_file_path, meta_data_path = meta_data_file_path)
    mpactr_object$setup()
    filter_object <- filter_pactr$new(mpactr_object)
    return(filter_object)
}