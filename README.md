
# CorMID

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/CorMID)](https://CRAN.R-project.org/package=CorMID)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/janlisec/CorMID/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janlisec/CorMID/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/janlisec/CorMID/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/janlisec/CorMID/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

**CorMID** is an R-package providing functions to solve problems during 
metabolic flux analysis using HR-APCI-MS.

In metabolic flux experiments tracer molecules (often glucose containing 
labelled carbon) are incorporated in compounds measured using mass spectrometry. 
The mass isotopologue distributions (MIDs) of these compounds needs to be 
corrected for natural abundance of labelled carbon and other effects, which are 
specific on the compound and ionization technique applied. This package provides 
functions to correct such effects in high resolution gas chromatography 
atmospheric pressure chemical ionization mass spectrometry (GC-HR-APCI-MS) 
analyses.

## Installation

You can install the development version of CorMID from 
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("janlisec/CorMID")
```

or install the version from [CRAN](https://cran.r-project.org/package=CorMID) 
instead.

## Quick Example

**CorMID** is supposed to disentangle a complex MID. Complex means that the ion
intensities of the isotopes are influenced by natural abundance, artificial
labeling (e.g. by a ^13^C-Glucose tracer) and mass spectrometry artifacts (i.e.
several potential adducts).

You can create and visualize such a complex mass spectrum by providing a 
chemical formula, the true labeling status and an adduct distribution like 
follows:

``` r
library(CorMID)
fml <- "C9H20O3Si2"
mid <- c(0.9, 0, 0, 0.1)
r <- list("M+H" = 0.8, "M-H" = 0.1, "M+H2O-CH4" = 0.1)
(rMID <- CorMID::recMID(mid = mid, r = r, fml = fml))
plot(rMID)
```
Assuming that you have measured these intensities in your experiment, **CorMID**
could estimate the underlying *MID* and *r* for you:

``` r
out <- CorMID::CorMID(int = rMID, fml=fml, prec=0.001, r=unlist(r))
print(out)
```

## Detailed documentation

You might either read the Vignette describing the package functions in detail

``` r
vignette("CorMID", package = "CorMID")
```

or read the [publication](https://doi.org/10.3390/metabo12050408) which shows a
evaluation of the performance of **CorMID** on real data sets.
