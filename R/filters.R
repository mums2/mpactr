#' @params data_frame with first three columns: Compound, mz, rt
solvent_blank_filter <- function(data_frame, sample_to_filter)
{
  result <- apply(data_frame[ , sample_to_filter, drop = FALSE], 1, function(x) {sum(x) <= 0})
  return(data_frame[result, !names(data_frame) %in% sample_to_filter])
}



# add params: df
check_mismatched_peaks <- function(data_frame, ringwin, isowin, trwin, max_iso_shift, merge_peaks) {

ion_filter_list <- list()
current_rt <- 0
current_mass <- 0
mass_diff <- 0
kmd_diff <- 0
cut_ions <- c() # list
merge_groups <- list() # dictonary

  data_frame <- data_frame[order(data_frame$mz, decreasing = FALSE), ]
  number_of_rows <- nrow(data_frame)
  for(i in seq_along(1:number_of_rows)) {
    current_feature <- data_frame[i, c("Compound", "mz", "rt")]
    current_rt <- current_feature$rt
    current_mass <- current_feature$mz
    current_ion <- current_feature$Compound
    if (!(current_ion %in% cut_ions))
    {
      for(j in (i + 1):number_of_rows) {
        if(j > number_of_rows)
        {
          break
        }
        if(data_frame$Compound[j] %in% cut_ions) {
          next
        }

        mass_diff <- data_frame$mz[j] - current_mass
        kmd_diff <- mass_diff - floor(mass_diff)

        if (abs(mass_diff) > max_iso_shift - 0.4) { # BL  - why 0.4??
          break
        }

        rt_diff <- data_frame$rt[j] - current_rt
        ring_band <- floor(abs(mass_diff) * (1/ringwin)) %% (1/ringwin)
        double_band <- kmd_diff - .5004 -
        (floor(mass_diff) * .007) # BL - why 0.500 and 0.007?

        if((abs(rt_diff) <= trwin) && 
           (mass_diff <= max_iso_shift - 0.4) &&
           (ring_band == 0 ||
            double_band < isowin)) {
            cut_ions <- c(cut_ions, data_frame$Compound[j])
            merge_groups[[as.character(current_ion)]] <- c(merge_groups[[as.character(current_ion)]], data_frame$Compound[j])
        }
      }
    }
  }
  # TODO Look into removing merge groups
  ion_filter_list[["cut_ions"]] <- cut_ions
  ion_filter_list[["merge_groups"]] <- merge_groups
  
  if (isTRUE(merge_peaks)) {
    return(data_frame[which(!(data_frame$Compound %in% cut_ions)), ])
  }
  return(ion_filter_list)

}
# Generate a class-like structure for the tree-plot
merge_ions <- function(data_frame, ion_filter_list) {

  data_frame[which(!(data_frame$Compound %in% cut_ions)), ]
  
}

#