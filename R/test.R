solvent_blank_filter <- function(sample_to_filter, data_frame)
{
  return(subset(data_frame, data_frame[sample_to_filter] <= 0))
}
