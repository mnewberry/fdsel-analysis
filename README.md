# Frequency-dependent selection in timeseries data, analysis
  
This repository contains code and data to reproduce the results and figures of
Newberry and Plotkin (2022) _Measuring frequency-dependent selection in
culture_.

## Dependencies

The analysis relies on the software `fdsel`, available in a separate repository
and R (3.6.1) with packages optparse and ggplot2, as well as a typical unix
environment.

The binaries `fdsel` and `fdsel-simdemo` (for Fig. S6) are required to be in
the path.

## Name analysis

The name analysis is in the `names` subdirectory. Executing `./run.sh` and
`Rscript plot.R` reproduces the analysis. Please note that bootstrap
confidence intervals required to reproduce main text figures take several
CPU-days to execute on modern CPUs.

## Dog popularity analysis

The dog popularity analysis is in the `akc` subdirectory. Executing `./run.sh`
and `Rscript plot.R` reproduces most of the analysis.

The novelty bias grid search (SI Section 5) is extremely computationally
intensive and takes dozens of CPU-days to execute. To run this analysis,
uncomment the relevant lines from `run.sh`.

## Author

All software was written by Mitchell Newberry <mitchell@localpost.io> and is
(c) Mitchell Newberry 2020-2022.  Bug reports and comments are welcome.

Reproductions of public data included in this repository for archival purposes
and convenience. Their copyright and licensing terms lie with thier respective
authors as described in `LICENSE.md`, `names/inp/PROVENANCE.md` and
`akc/inp/PROVENANCE.md`.
