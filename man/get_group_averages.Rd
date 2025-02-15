% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_pactr-accessors.R
\name{get_group_averages}
\alias{get_group_averages}
\title{Get groups averages.}
\usage{
get_group_averages(mpactr_object)
}
\arguments{
\item{mpactr_object}{The mpactr object that is created by calling the
import_data() function.}
}
\value{
a \code{data.table} reporting the average and relative standard
deviation across biological groups and technical replicates within
each group.
}
\description{
\code{get_group_averages()} is a wrapper function to return group averages
for the filtered peak table.
}
\examples{
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_group(data, group_to_remove = "Blanks")

group_averages <- get_group_averages(data_filter)
head(group_averages)

}
