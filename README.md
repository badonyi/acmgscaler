# `acmgscaler` <a href='https://colab.research.google.com/github/badonyi/acmgscaler/blob/main/acmgscaler.ipynb'><img src='hexlogo.png' align="right" height="198" /></a>

<!-- badges: start -->
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/badonyi/acmgscaler/blob/main/acmgscaler.ipynb)
[![CRANstatus](https://www.r-pkg.org/badges/version/acmgscaler)](https://cran.r-project.org/package=acmgscaler)
[![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/acmgscaler?color=blue)](https://CRAN.R-project.org/package=acmgscaler)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/badonyi/acmgscaler?branch=main&svg=true)](https://ci.appveyor.com/project/badonyi/acmgscaler/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/badonyi/acmgscaler/blob/main/LICENSE)
[![DOI:10.1093/bioinformatics/btaf503](https://img.shields.io/badge/DOI-10.1093/bioinformatics/btaf503-B31B1B.svg)](https://doi.org/10.1093/bioinformatics/btaf503)
<!-- badges: end -->

The goal of the `acmgscaler` R package is to provide a robust approach for gene-level calibration of variant effect scores, such as computational predictions or functional assay results, against ACMG/AMP evidence thresholds.
The package is lightweight and written entirely in base R, without additional dependencies.


## Colab notebook
A plug-and-play [Google Colab](https://colab.research.google.com/github/badonyi/acmgscaler/blob/main/acmgscaler.ipynb) notebook with a simple interface is available for all users.


## Installation
You can install the stable version of `acmgscaler` directly from R:

```r
install.packages('acmgscaler')

# or install the development version from GitHub
devtools::install_github('badonyi/acmgscaler')
```


## Quick start

```r
library(acmgscaler)
data('variant_data', package = 'acmgscaler')

# calibrate the example data
calib <- calibrate(
  df = variant_data,
  value = 'score', 
  prior = 0.1,
  group = 'gene'
)

# likelihood ratios for each variant
calib$BRCA1$likelihood_ratios

# score thresholds for ACMG/AMP evidence levels
calib$BRCA1$score_thresholds
```


## Reference
```
@article{acmgscaler,
  title={acmgscaler: An R package and Colab for standardised gene-level variant effect score calibration within the ACMG/AMP framework},
  author={Badonyi, Mihaly and Marsh, Joseph A},
  journal={Bioinformatics},
  pages={btaf503},
  year={2025},
  publisher={Oxford University Press}
}
```
