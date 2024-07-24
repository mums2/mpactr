
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mpactR

<!-- badges: start -->
<!-- badges: end -->

## Overview

mpactR is a collection of filters for the purpose of identifying high
quality MS1 features by correcting peak selection errors introduced
during the pre-processing of tandem mass spectrometry data.

Filters in this package address the following issues:

- `filter_mispicked_ions()`: removal of mispicked peaks, or those
  isotopic patterns that are incorrectly split during preprocessing.
- `filter_group()`: removal of features overrepresented in a specific
  group of samples; for example removal of features present is solvent
  blanks due to carryover between samples.
- `filter_cv()`: removal of non-reproducible features, or those that are
  inconsistent between technical replicates.
- `filter_insource_ions()`: removal of fragment ions created during the
  first ionization in the tandem MS/MS workflow.

All filters are independent, meaning they can be used to create a
project-specific workflow, or you can learn more in **cite vignette**.

## Installation

You can install the development version of mpactR from
[GitHub](https://github.com/SchlossLab/mpactR) with:

``` r
# install.packages("devtools")
devtools::install_github("SchlossLab/mpactR")
```

## Get started

See articles to get started.

## Getting help

If you encounter an issue, please file an issue on
[GitHub](https://github.com/SchlossLab/mpactR/issues). Please include a
minimal reproducible example with your issue.

## Contributing

Is there a feature you’d like to see included, please let us know! Pull
requests are welcome on
[GitHub](https://github.com/SchlossLab/mpactR/pulls).