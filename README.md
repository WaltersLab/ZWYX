
# ZWYX

<!-- badges: start -->
<!-- badges: end -->

# Overview

The ZWYX package aims to provide a streamlined workflow for identifying sex-linked regions of a genome assembly by comparing sequencing depth between males and females.

## Principles: Sex linkage and sequencing coverage

In species with differentiated sex chromosomes, the X or Z chromosomes should show a two-fold difference in sequencing coverage between sexes, while the sex-specific Y or W chromosomes should show an opposing and far more biased pattern of coverage. Autosomes should show no differences.

Using sex-specific sequencing coverage is an established method for identifying scaffolds derived from X, Y, W, or Z chromosomes. ZWYX provides data structures and functions for identifying and visualizing genome scaffolds with sex-biased coverage on average as a whole, but also via windows along each scaffold. 

ZWYX also incorporates "changepoint" algorithms for detecting shifts in sex-specific coverage that occur within a scaffold, which typically indicate a mis-assembly that erroneously joins autosomal and sex-linked regions.

## Installation

You can install the development version of ZWYX as follows:

``` r
remotes::install_github("WaltersLab/ZWYX", build_vignettes = TRUE)
```

Using `build_vignettes = TRUE` is important so that the vignette is rendered.

Then you can access the vignette via:

``` r
browseVignettes("ZWYX")
```

