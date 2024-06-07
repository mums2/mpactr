#######################
### Mistatched peaks ##
#######################

filter_mismatch_ions <- function(mpactr_object, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3,
                                merge_peaks = TRUE)
{
  mpactr_object$check_mismatched_peaks(ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3,
                                merge_peaks = merge_peaks)
  return(mpactr_object)
}

#######################
### Group filter     ##
#######################
filter_group <- function(mpactr_object, group_threshold, group_to_remove, remove_ions)
{
  mpactr_object$filter_blank()
  mpactr_object$parse_ions_by_group(group_threshold = group_threshold)
  mpactr_object$apply_group_filter(group_to_remove, remove_ions = remove_ions)
  return(mpactr_object)
}

#######################
### CV filter        ##
#######################
filter_cv <- function(mpactr_object, cv_threshold, cv_params) {
   mpactr_object$cv_filter(cv_threshold = cv_threshold, cv_params = cv_params)
  return(mpactr_object)
}


###########################
### Insource ions filter ##
###########################
filter_insource_ions <- function(mpactr_object, cluster_threshold = 0.95) {

  mpactr_object$filter_insource_ions(cluster_threshold = cluster_threshold)
  return(mpactr_object)
}

