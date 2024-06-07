import_data <- function(peak_table_flle_path, meta_data_file_path)
{
    mpactr_object <- mpactr$new(peak_table = peak_table_flle_path, meta_data = meta_data_file_path)
    mpactr_object$setup()
    filter_object <- filter_pactr$new(mpactr_object)
    return(filter_object)
}