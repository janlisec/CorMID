
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CorMID

<!-- badges: start -->

[![R-CMD-check](https://github.com/janlisec/CorMID/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janlisec/CorMID/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![test-coverage](https://github.com/janlisec/CorMID/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/janlisec/CorMID/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/janlisec/CorMID/branch/main/graph/badge.svg?token=NSY6DITZVH)](https://app.codecov.io/gh/janlisec/CorMID)
[![CRAN
status](https://www.r-pkg.org/badges/version/CorMID)](https://CRAN.R-project.org/package=CorMID)
[![DOI](https://img.shields.io/badge/doi-10.3390/metabo12050408-yellow.svg)](https://doi.org/10.3390/metabo12050408)
<!-- badges: end -->

**CorMID** is an R-package providing functions to solve problems during
metabolic flux analysis using HR-APCI-MS.

In metabolic flux experiments tracer molecules (often glucose containing
labelled carbon) are incorporated in compounds measured using mass
spectrometry. The mass isotopologue distributions (MIDs) of these
compounds needs to be corrected for natural abundance of labelled carbon
and other effects, which are specific on the compound and ionization
technique applied. This package provides functions to correct such
effects in high resolution gas chromatography atmospheric pressure
chemical ionization mass spectrometry (GC-HR-APCI-MS) analyses.

## Installation

You can install the development version of CorMID from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("janlisec/CorMID")
```

or install the version from
[CRAN](https://cran.r-project.org/package=CorMID) instead.

## Quick Example

**CorMID** is supposed to disentangle a complex MID. Complex means that
the ion intensities of the isotopes are influenced by natural abundance,
artificial labeling (e.g. by a <sup>13</sup>C-Glucose tracer) and mass
spectrometry artifacts (i.e. several potential adducts).

You can create and visualize such a complex mass spectrum by providing a
chemical formula, the true labeling status and an adduct distribution
like follows:

``` r
# a chemical formula, here: Lactic acid 2 TMS
fml <- "C9H22O3Si2"

# the true labeling status, here: 10% U13C enriched
mid <- c(0.9, 0, 0, 0.1)

# adduct distribution, here: Three different APCI adducts formed
r <- list("M+H" = 0.8, "M-H" = 0.1, "M+H2O-CH4" = 0.1)

# reconstruct and plot the measured intensity vector
(rMID <- CorMID::recMID(mid = mid, r = r, fml = fml))
#>        M-2        M-1        M+0        M+1        M+2        M+3        M+4 
#> 0.06901842 0.01398425 0.55840285 0.12063392 0.12056965 0.08474875 0.01838686 
#>        M+5 
#> 0.01425529 
#> attr(,"class")
#> [1] "recMID"
```

**CorMID** provides a class specific plotting function for such a
reconstructed MID:

``` r
plot(rMID, ylim=c(0,0.6))
mtext(text = "Reconstructed MID", side = 3, line = -1.25, adj = 0.98, outer = T)
text(x = 3, y = 0.4, labels = "[M+H]+", pos=2)
```

![](man/figures/README-exmpl1_plot-1.png)<!-- -->

Assuming that you have measured these intensities in your experiment,
now **CorMID** could estimate the underlying *MID* and *r* for you:

``` r
# disentangle the adduct ratios and true enrichment from the above test data
out <- CorMID::CorMID(int = rMID, fml = fml)
print(out)
#> [class] 'CorMID'
#> MID [%] (estimated)
#>     M0    M1    M2    M3
#>  89.06 00.00 00.78 10.16
#> [attr] 'r' (estimated)
#>   M+H   M+   M-H   M+H2O-CH4
#>  0.81 0.00  0.10        0.09
#> [attr] 'err'
#> 0.003167
```

Please note: no information regarding the true labeling status and the
adduct distribution was provided in this function call. **CorMID** is
able to *guess* the most likely combination.

This allows you to perform the correction for natural abundance and
technical artifacts in a single step and extract the relevant labelling
status for flux analysis or other statistical evaluations.

## Detailed documentation

You might either read the
[Vignette](https://cran.r-project.org/package=CorMID/vignettes/CorMID.html)
describing the package functions in detail or read the
[publication](https://doi.org/10.3390/metabo12050408) which shows a
evaluation of the performance of **CorMID** on real data sets.
