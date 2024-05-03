set_kmd <- function(data_frame) {
  data_frame$kmd <- data_frame$mz - floor(data_frame$mz)
  return(data_frame)
}

