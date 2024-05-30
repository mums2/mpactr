set_kmd <- function(data_frame) {
  data_frame$kmd <- data_frame$mz - floor(data_frame$mz)
  return(data_frame)
}

initialize_data <- function(data_frame) {
  # data_frame$Compound <- as.numeric(data_frame$Compound)
  data_frame <- data_frame[which( rowSums(data_frame[, !(colnames(data_frame) %in% c("Compound", "mz", "rt")) ] ) > 0), ]

  # calc kmd
  data_frame <- set_kmd(data_frame)

  return(data_frame)
}