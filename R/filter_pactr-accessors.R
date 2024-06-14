#' Return the summary for a single mpactR filter.
#'
#' @description 
#' `filter_summary()` is a wrapper function to return the summary from a single filter within the given mpactr object. 
#'
#' @param mpactr_object The mpactr object that is created by calling the import_data() function.
#'
#' @param filter The name of a filter whose summary is to be extracted. Must be one of: "mispicked", "group", "replicability", or "insource".
#'
#' @return a `list` reporting 1) compound ids for compounds which failed the filter and 2) compound ids for compounds which passed the filter.
#' @export 
#'
#' @examples 
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
#' data_filter <- filter_mispicked_ions(data)
#'
#' mispicked_summary <- filter_summary(data_filter, filter = "mispicked")
#' mispicked_summary
#'
filter_summary <- function(mpactr_object, filter) {
    return(mpactr_object$get_log(filter = filter))
}

