#' @params data_frame with first three columns: Compound, mz, rt
solvent_blank_filter <- function(data_frame, sample_to_filter)
{
  result <- apply(data_frame[ , sample_to_filter, drop = FALSE], 1, function(x) {sum(x) <= 0})
  return(data_frame[result, !names(data_frame) %in% sample_to_filter])
}



# add params: df
check_mismatched_peaks <- function(data_frame, ringwin, isowin, trwin, max_iso_shift, merge_peaks) {

ion_position <- 0
ion_comp_pos <- 0
current_rt <- 0
current_mass <- 0
mass_diff <- 0
kmd_diff <- 0

cut_ions <- c() # list
mergegroups <- list() # dictonary

  data_frame <- data_frame[order("mz", decreasing = FALSE),]
  for(i in seq_along(1:range(data_frame)))
  {
    current_feature <- data_frame[i,]
    current_rt <- current_feature$rt
    current_mass <- current_feature$mz
  }
  


}




# #def relationalfilter(analysis_params, ionfilterlist):
#     """
#     Implements a de-ringing and relational filtering algorithm for MS data.

#     Parameters:
#         analysis_params (object): Analysis parameters.
#         ionfilterlist (dict): Dictionary of ion filters.

#     Returns:
#         dict: Dictionary of ion filters.

#     The algorithm works by first sorting the MS data by m/z, and then iterating over the data one ion at a time, 
#     comparing each ion to all subsequent ions that fall within a certain mass and retention time window. 
#     Ions that meet certain criteria are then merged into a single ion using the ionmerge class. 

#     The filtering criteria includes:
#     - A mass difference of less than a specified value
#     - A retention time difference of less than a specified value
#     - A difference between the mass difference and its floor when multiplied by a factor, indicating the presence of ring artifacts
#     - A difference between the mass difference and 0.5004 when accounting for the doubly-charged dimer, indicating the presence of dimer peaks

#     The function returns a dictionary of ion filters, with the 'relfil' filter containing a list of ions to remove, 
#     and a 'merge' attribute which contains the merge groups (as an ionmerge object) that specify which ions to merge.
#     """
    
#     print('Running relational filter')
#     pdindex = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), sep = ',', header = [2], index_col = None)
#     pdindex.drop(pdindex.columns.difference(['Compound','m/z','Retention time (min)']), 1, inplace=True)
#     pdindex = pdindex.sort_values(['m/z'], ascending=[1])
#     npindex = pdindex.to_numpy()

#     ion_position, ion_comp_pos, current_rt, current_mass, mass_diff, kmd_diff, cut_ions, mergegroups = 0, 0, 0, 0, 0, 0, [], {}
#     for ions in npindex[0:,0]: #scrolls through ions in order and checks if any subsequent ions fall within time and mass filter range, any ions that are are copied to a list
#         current_rt = npindex[ion_position, 2]
#         current_mass = npindex[ion_position, 1]
#         current_ion = npindex[ion_position, 0]
#         ion_comp_pos = ion_position + 1
#         if (current_ion not in cut_ions):
#             for ions in npindex[ion_comp_pos:, 0]: # the following statement checks if any peaks are within 0.05 min and in the mass filtering window.
#                     #The mass filtering window works by checking if the difference of the two ions, when multiplied by two and floored, is even this
#                     #corresponds to filtering out ions within 0.5 daltons of a peak or any of its first five isotopic peaks.
#                 mass_diff = abs((npindex[ion_comp_pos, 1] - current_mass))
#                 kmd_diff = mass_diff - math.floor(mass_diff)
#                 if ((npindex[ion_comp_pos, 1] - current_mass) > analysis_params.maxisowin - .4):
#                     break    # breaks loop if we are already past the maximum consideration window for isotope shift, saves considerable processsing time.
#                 if ((abs(npindex[ion_comp_pos, 2]-current_rt) <= analysis_params.RTwin) # RT band
#                     and ((npindex[ion_comp_pos, 1] - current_mass) <= analysis_params.maxisowin - .4)
#                     and ((math.floor(mass_diff * (1/analysis_params.ringingwin)) % (1/analysis_params.ringingwin) == 0) #band for ringing
#                     or ((kmd_diff - .5004 - (math.floor(mass_diff) * .007))  < analysis_params.dimerpeakwin)) #band for doubly charged dimer removal, CHECK TO SEE IF CAN USE MOD FOR THIS
#                     and (npindex[ion_comp_pos, 0] not in cut_ions)): # check if ion already removed
                    
#                     cut_ions.append(npindex[ion_comp_pos, 0]) #creates or appends merge group to indicate which ions to be merged. Merge not yet implimented
#                     if current_ion in mergegroups:
#                         mergegroups[current_ion].sources.append(npindex[ion_comp_pos, 0])
#                     else:
#                         mergegroups[current_ion] = ionmerge(current_ion, npindex[ion_comp_pos, 0])
#                 ion_comp_pos += 1
#         ion_position += 1
#     ionfilterlist['relfil'] = ionfilter('', cut_ions)
#     ionfilterlist['relfil'].merge = mergegroups

#     iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep = ',', header = [0], index_col = 0)
#     iondict['pass_relfil'] = ~iondict.index.isin(ionfilterlist['relfil'].ions)
#     iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header = True, index = True) #saves formatted backup for later use

#     return(ionfilterlist)
