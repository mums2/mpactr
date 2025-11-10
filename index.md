# mpactr

## Overview

mpactr is a collection of filters for the purpose of identifying high
quality MS1 features by correcting peak selection errors introduced
during the pre-processing of tandem mass spectrometry data.

Filters in this package address the following issues:

- [`filter_mispicked_ions()`](https://www.mums2.org/mpactr/reference/filter_mispicked_ions.md):
  removal of mispicked peaks, or those isotopic patterns that are
  incorrectly split during preprocessing.
- [`filter_group()`](https://www.mums2.org/mpactr/reference/filter_group.md):
  removal of features overrepresented in a specific group of samples;
  for example removal of features present in solvent blanks due to
  carryover between samples.
- [`filter_cv()`](https://www.mums2.org/mpactr/reference/filter_cv.md):
  removal of non-reproducible features, or those that are inconsistent
  between technical replicates.
- [`filter_insource_ions()`](https://www.mums2.org/mpactr/reference/filter_insource_ions.md):
  removal of fragment ions created during the first ionization in the
  tandem MS/MS workflow.

All filters are independent, meaning they can be used to create a
project-specific workflow, or you can learn more in [the Getting Started
page](https://www.mums2.org/mpactr/articles/mpactr.html).

## Installation

You can install the CRAN version with:

``` r
install.packages("mpactr")
```

You can install the development version of mpactr from
[GitHub](https://github.com/mums2/mpactr) with:

``` r
# install.packages("devtools")
devtools::install_github("mums2/mpactr")
```

## Get started

See the [Getting
Started](https://www.mums2.org/mpactr/articles/mpactr.html) page to get
started.

## Getting help

If you encounter an issue, please file an issue on
[GitHub](https://github.com/mums2/mpactr/issues). Please include a
minimal reproducible example with your issue.

## Contributing

Is there a feature you’d like to see included, please let us know! Pull
requests are welcome on [GitHub](https://github.com/mums2/mpactr/pulls).
