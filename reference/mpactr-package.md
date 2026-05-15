# mpactr: Correction of Preprocessed MS Data

A fast tool for the correction of preprocessed tandem MS/MS features.
This is an 'R' implementation of the 'python' program Metabolomics Peak
Analysis Computational Tool ('MPACT') (Robert M. Samples, Sara P.
Puckett, and Marcy J. Balunas (2023)
[doi:10.1021/acs.analchem.2c04632](https://doi.org/10.1021/acs.analchem.2c04632)
). Filters in the package serve to address common errors in MS/MS
preprocessing, including: (1) isotopic patterns that are incorrectly
split during preprocessing, (2) features present in solvent blanks due
to carryover between samples, (3) features whose abundance is greater
than user-defined abundance threshold in a specific group of samples,
for example media blanks, (4) ions that are inconsistent between
technical replicates, and (5) in-source fragment ions created during
ionization before fragmentation in the tandem MS/MS workflow.

## See also

Useful links:

- <https://mums2.github.io/mpactr/>

- <https://github.com/mums2/mpactr>

- Report bugs at <https://github.com/mums2/mpactr/issues>

## Author

**Maintainer**: Patrick Schloss <pschloss@umich.edu>
([ORCID](https://orcid.org/0000-0002-6935-4275)) \[copyright holder\]

Authors:

- Allison Mason <masonar@umich.edu>
  ([ORCID](https://orcid.org/0000-0003-1339-1592))

- Gregory Johnson <grejoh@umich.edu>
