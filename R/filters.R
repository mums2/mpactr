
#######################
### Mistatched peaks ##
#######################

#' Mispicked ions filter
#'
#' @description
#' `filter_mispicked_ions()` identifies ions that were incorrectly split into separate features during preprocessing. This filter checks the feature table for similar ions in terms of mass and retention time. Peaks found to be similar are merged into a single feature given `merge_peaks` is `TRUE`.
#'
#' @param mpactr_object An `mpactr_object`. See [import_data()].
#' @param ringwin Ringing mass window. Default = 0.5.
#' @param isowin Isotopic mass window. Defualt = 0.01.
#' @param trwin A `numeric` denoting the retention time threhold for assessing if ions should be merged. Defulat = 0.005.
#' @param max_iso_shift A `numeric`. Default = 3.
#' @param merge_peaks A `logical` to determine if peaks found to belong to the same ion should be merged in the feature table.
#' @param copy_object A `boolean` paramter that allows users to return a copied object instead of modifying the object.
#' @return an `mpactr_object`
#' @export 
#'
#' @examples 
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
#' data_filter <- filter_mispicked_ions(data,
#'                               ringwin = 0.5, 
#'                               isowin = 0.01,
#'                               trwin = 0.005,
#'                               max_iso_shift = 3,
#'                               merge_peaks = TRUE)
#'
filter_mispicked_ions <- function(mpactr_object, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3,
                                merge_peaks = TRUE, copy_object = FALSE)
{
  if(copy_object)
  {
    mpactr_object <- clone(mpactr_object)
  }
  mpactr_object$check_mismatched_peaks(ringwin = ringwin, isowin = isowin, trwin = trwin, max_iso_shift = max_iso_shift,
                                merge_peaks = merge_peaks)
  return(mpactr_object)
}

#######################
### Group filter     ##
#######################
#' Filter Ions by Group
#'
#' @details 
#' `filter_group()` removes ions which have a signigifanct presence in the defined group (typically solvent blanks). Significance is determined by relative ion abundance arcoss groups at the user-defined `group_threshold`.
#'
#' Filtering is assessed on function parameter `group_threshold`. The defualt is 0.01, meaning features present in `group_to_remove` at greater than 1% are removed. 
#'
#' @param mpactr_object An `mpactr_object`. See [import_data()].
#' @param group_threshold Relative abundance threshold at which to remove ions. Default = 0.01.
#' @param group_to_remove Biological group name to remove ions from.
#' @param remove_ions A `logical`. If `TRUE` failing ions will be removed from the peak table. Default = TRUE.
#' @param copy_object A `boolean` paramter that allows users to return a copied object instead of modifying the object.

#' @return an `mpactr_object`
#' @export 
#'
#' @examples 
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
#' data_filter <- filter_group(data,
#'                               group_threshold = 0.01,
#'                               group_to_remove = "Blanks",
#'                               remove_ions = TRUE)
#'
filter_group <- function(mpactr_object, group_threshold = 0.01, group_to_remove, remove_ions = TRUE,
                          copy_object = FALSE)
{
   if(copy_object)
  {
    mpactr_object <- clone(mpactr_object)
  }
  mpactr_object$filter_blank()
  mpactr_object$parse_ions_by_group(group_threshold = group_threshold)
  mpactr_object$apply_group_filter(group_to_remove, remove_ions = remove_ions)
  return(mpactr_object)
}

#######################
### CV filter        ##
#######################
#' Filter Non-reproducible ions
#'
#' @description 
#' `filter_cv()` removes feature ions that are found to be non-reproducible between technical inejction replicates. Replicability is assessed via mean or median coefficient of variation (CV) between technical replicates. As such, this fitler is expecting an input dataset with at least two replicate injections per sample. 
#'
#'
#' @param mpactr_object An `mpactr_object`. See [import_data()].
#' @param cv_threshold Coefficient of variation threshold. Default = 0.2.
#' @param cv_param Coefficient of variation (CV) to use for filtering. Options are "mean" or "median", corresponding to mean and median CV, respectively. 
#' @param copy_object A `boolean` paramter that allows users to return a copied object instead of modifying the object.
#'
#' @return an `mpactr_object`
#' @export 
#'
#' @examples
#'
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
#' data_filter <- filter_cv(data,
#'                               cv_threshold = 0.01,
#'                               cv_param = "mean")
#'
#' data_filter <- filter_cv(data,
#'                               cv_threshold = 0.01,
#'                               cv_param = "median")
#'
filter_cv <- function(mpactr_object, cv_threshold = 0.2, cv_param, copy_object = FALSE) {
   if(copy_object)
  {
    mpactr_object <- clone(mpactr_object)
  }
   mpactr_object$cv_filter(cv_threshold = cv_threshold, cv_params = cv_param)
  return(mpactr_object)
}


###########################
### Insource ions filter ##
###########################
#' Filter Insource ions
#'
#' @description 
#' `filter_insource_ions()` determines insource ion fragments deconvolution via hierarchical clustering. Highly correlated ions with the same retention times are identified and removed. 
#'
#'
#' @param mpactr_object An `mpactr_object`. See [import_data()].
#' @param cluster_threshold cluster threshold for ion deconvolution. Default = 0.95.
#' @param copy_object A `boolean` paramter that allows users to return a copied object instead of modifying the object.
#'
#' @return an `mpactr_object`
#' @export 
#'
#' @examples 
#' data <- import_data(example("coculture_peak_table.csv"),
#'                     example("metadata.csv"))
#'
#' data_filter <- filter_insource_ions(data,
#'                               cluster_threshold = 0.95)
#'
filter_insource_ions <- function(mpactr_object, cluster_threshold = 0.95, copy_object = FALSE) {
 if(copy_object)
  {
    mpactr_object <- clone(mpactr_object)
  }
  mpactr_object$filter_insource_ions(cluster_threshold = cluster_threshold)
  return(mpactr_object)
}

