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

#' Get similar ion groups.
#'
#' @description 
#' `similar_ions()` is a wrapper function to return similar ion groups determined with the [filter_mispicked_ions()]. 
#'
#' @param mpactr_object The mpactr object that is created by calling the import_data() function.
#'
#' @return a `data.table` retporting the main ion and those found to be similar with [filter_mispicked_ions()].
#' @export 
#'
#' @examples 
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
#' data_filter <- filter_mispicked_ions(data)
#'
#' mispicked_ion_groups <- similar_ions(data_filter)
#' mispicked_ion_groups
#'
similar_ions <- function(mpactr_object) {
    return(mpactr_object$get_mispicked_ions())
}

#' Get groups averages.
#'
#' @description 
#' `group_averages()` is a wrapper function to return group averages for calculated with [filter_group()]. 
#'
#' @param mpactr_object The mpactr object that is created by calling the import_data() function.
#'
#' @return a `data.table` retporting the average and relative standard deviation across biological groups and techincal replicates within each group].
#' @export 
#'
#' @examples 
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
#' data_filter <- filter_group(data, group_to_remove = "Blanks")
#'
#' group_averages <- group_averages(data_filter)
#' head(group_averages)
#'
group_averages <- function(mpactr_object) {
    return(mpactr_object$get_group_averages())
}
