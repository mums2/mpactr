solvent_blank_filter <- function(data_frame, sample_to_filter)
{
  result <- apply(data_frame[ , sample_to_filter, drop = FALSE], 1, function(x) {sum(x) <= 0})
  return(data_frame[result, !names(data_frame) %in% sample_to_filter])
}



# test_ex_df <- data.frame("sample1_1" = c(0, 344, 592, 4, 1),
#                          "sample1_2" = c(1, 350, 600, 8, 3),
#                          "blank_1" = c(0, 0, 0, 3, 0),
#                          "blank_2" = c(0, 0, 0, 0, 1))
# filtered_df <- solvent_blank_filter(c("blank_1", "blank_2"), test_ex_df)

# to_filter <- c("blank_1", "blank_2")
#   # find which rows across all samples in group are less <= 0
# result <- apply(test_ex_df[ , to_filter, drop = FALSE], 1, function(x) {sum(x) <= 0})
  
#   # remove rows with detection in blanks
# test_ex_df[result, !names(test_ex_df) %in% to_filter] 
    