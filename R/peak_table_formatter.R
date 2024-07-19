format_by_type <- function(peak_table, sample_names, type_of_peak_table)
{
    if(type_of_peak_table == "Progensis"){
        return(progensis_formatter(peak_table))    
    }
    else if(type_of_peak_table == "MzMine"){
        return(mz_mine_formatter(peak_table))
    }   
    else if(type_of_peak_table == "Metaboscape"){
        return(metaboscape_formatter(peak_table, samples_names))
    }
    else{}# default condition = NULL
}

progensis_formatter <- function(peak_table){

    peak_table <- data.table(readr::read_csv(peak_table, skip = 2, show_col_types = FALSE))
    raw_peak_table <- peak_table
    setnames(peak_table, c("m/z", "Retention time (min)"), c("mz", "rt"))
    
    return(list("peak_table" = peak_table,
        "raw_table" = raw_peak_table))
}

mz_mine_formatter <- function(peak_table){
    peak_table <- data.table(readr::read_csv(peak_table_path, skip = 2, show_col_types = FALSE))
    raw_peak_table <- peak_table




    return(list("peak_table" = peak_table,
        "raw_table" = raw_peak_table))
}
metaboscape_formatter <- function(peak_table, sample_names){
    peak_table <- data.table(readr::read_csv(peak_table, show_col_types = FALSE))

    adduct_data <- utils::read.csv(system.file("extdata/ion_masses", "DefinedIons.csv", package = "mpactR") )
    
    peak_table_convert <- peak_table[ , ion := gsub(".*\\[(.+)\\].*", "\\1", ADDUCT)][ 
                                                        , charge_string := gsub(".*\\](.+)", "\\1", ADDUCT)][
                                                            charge_string == "+", charge := 1][
                                                                charge_string == "2+", charge := 2][
                                                                    charge_string == "3+", charge := 3][
                                                                        adduct_data, on = .(ion = IONS)][
                                                                            , mz := (PEPMASS / charge) + MASS]

    setnames(peak_table_convert, c("FEATURE_ID", "RT"), c("Compound", "rt"))
    
    peak_table_mpactr <- peak_table_convert[ , .SD, .SDcols = c("Compound", "mz", "rt", sample_names)]

    return(list("raw_table" = peak_table_convert,
        "peak_table" = peak_table_mpactr))
}

