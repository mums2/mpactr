#' Import data into an mpactr object.
#'
#' @description
#' `import_data()` takes two file paths, one for the pre-processed feature
#' table and one for sample metadata. Both files should be .csv.
#'
#' @details
#' mpactR requires 2 files as imput: a feature table and metadata file. Both
#' are expected to be comma separated files (*.csv*).
#'
#' 1. peak_table: a peak table in Progenesis format is expected. To export a
#' compatable peak table in Progenesis, navigate to the
#' *Review Compounds* tab then File -> Export Compound Measurements. Select
#' the following properties: Compound, m/z,
#' Retention time (min), and Raw abundance and click ok.
#' 2. metadata: a table with sample information. At minimum the following
#' columns are expected: Injection, Sample_Code, Biological_Group. Injection
#' is the sample name and is expected to match sample column names in the
#' peak_table. Sample_Code is the id for technical replicate groups.
#' Biological Group is the id for biological replicate groups. Other sample
# metadata can be added, and is encouraged for downstream analysis following
#' filtering with mpactR.
#'
#' @param peak_table The file path to your feature table file.
#' @param meta_data The file path to your meta_data file or `data.frame`.
#' @param format The expected exported type of your peak table, can be
#' one of Progenesis, MzMine, Metaboscape, or None.
#'
#' @return an `mpactr_object`
#' @export
#'
#' @examples
#' data <- import_data(example("coculture_peak_table.csv"),
#'   example("metadata.csv"),
#'   format = "Progenesis"
#' )
#'
#' meta_data <- read.csv(example("metadata.csv"))
#' data <- import_data(example("coculture_peak_table.csv"),
#'   meta_data,
#'   format = "Progenesis"
#' )
#'
import_data <- function(peak_table, meta_data, format = "none") {
  if (!any(class(meta_data) %in% c("data.table", "data.frame"))) {
    meta_data <- data.table(readr::read_csv(meta_data, show_col_types = FALSE))
  }

  #*** check for Injection, Sample_Code, Biological_Group
  cols <- c("Injection", "Sample_Code", "Biological_Group")
  if (any(cols %in% colnames(meta_data) == FALSE)) {
    cli::cli_abort("{.cls {cols[which(!(cols %in% colnames(meta_data)))]}}
                    are not columns in the provided metadata. Please see
                     function documentation for more details.")
  }

  df <- format_by_type(
    peak_table_path = peak_table,
    type_of_peak_table = format,
    sample_names = meta_data$Injection
  )

  mpactr_object <- mpactr$new(
    peak_table = df,
    meta_data = data.table(meta_data)
  )
  mpactr_object$setup()
  filter_object <- filter_pactr$new(mpactr_object)
  return(filter_object)
}
