# TODO Generate a class-like structure for the tree-plot



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
  if(current_ion == 1248)
        {
          print("stop")
        }
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

        if(current_ion == 1248)
        {
          print("stop")
           if( data_frame$Compound[j] == 1250)
          {print("stop")}
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

# blank filter
filter_blank <- function(data_frame, metadata) {
  
}

# grpave code: 
  msdata_errprop = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), sep=',', header=[0, 1, 2], index_col=[0, 1, 2])

    # Find technical averages and RSDs
    msdata_errprop_tidy = msdata_errprop.stack([0, 1, 2])
    msdata_errprop_mean = msdata_errprop_tidy.groupby(level=[0, 1, 2, 3, 4]).mean()
    msdata_errprop_mean = msdata_errprop_mean.groupby(level=[0, 1, 2, 3]).mean().to_frame()
    msdata_errprop_mean.columns = ['average']
    msdata_errprop_mean['biolRSD'] = msdata_errprop_tidy.groupby(level=[0, 1, 2, 3]).std().fillna(0) / msdata_errprop_tidy.groupby(level=[0, 1, 2, 3]).mean()
    msdata_errprop_mean['bioln'] = msdata_errprop_tidy.groupby(level=[0, 1, 2, 3]).count()
    msdata_errprop_sd = (msdata_errprop_tidy.groupby(level=[0, 1, 2, 3, 4]).std() / msdata_errprop_tidy.groupby(level=[0, 1, 2, 3, 4]).mean()).to_frame()
    msdata_errprop_mean['techRSD'] = msdata_errprop_sd.groupby(level=[0, 1, 2, 3]).mean()
    msdata_errprop_n = msdata_errprop_tidy.groupby(level=[0, 1, 2, 3, 4]).count()
    msdata_errprop_mean['techn'] = msdata_errprop_n.groupby(level=[0, 1, 2, 3]).mean()

    # Save summary data and group averages
    msdata_errprop_mean.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_summarydata.csv'), header=True, index=True)
    msdata_errprop_grpav = msdata_errprop_mean.loc[:, msdata_errprop_mean.columns.intersection(['average'])]
    msdata_errprop_grpav.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_groupaverages.csv'), header=True, index=True)

if analysis_params.grpave:
        stats.groupave(analysis_params)
        print('Parsing ion lists')
        groupionlists = filter.parsionlists(analysis_params)
if analysis_params.blnkfltr:
        msdata = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), sep=',', header=[0, 1, 2], index_col=None)
        msdata = filter.listfilter(msdata, groupionlists[analysis_params.blnkgrp], False)
        msdata = msdata.drop(analysis_params.blnkgrp, axis=1, level=0)
        msdata.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), header=True, index=False)
        iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=0)
        iondict['pass_blnkfil'] = ~iondict.index.isin(groupionlists[analysis_params.blnkgrp])
        iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=True)