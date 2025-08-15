# `acmgscaler` <a href='https://colab.research.google.com/github/badonyi/acmgscaler/blob/main/acmgscaler.ipynb'><img src='hexlogo.png' align="right" height="198" /></a>

<!-- badges: start -->
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/badonyi/acmgscaler/blob/main/acmgscaler.ipynb)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/badonyi/acmgscaler?branch=main&svg=true)](https://ci.appveyor.com/project/badonyi/acmgscaler/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/badonyi/acmgscaler/blob/main/LICENSE)
[![DOI:10.1101/2025.05.16.654507v2](https://img.shields.io/badge/DOI-10.1101/2025.05.16.654507v1-B31B1B.svg)](https://www.biorxiv.org/content/10.1101/2025.05.16.654507v2)
<!-- badges: end -->

The goal of the `acmgscaler` R package is to provide a robust approach for gene-level calibration of variant effect scores, such as computational predictions or functional assay results, against ACMG/AMP evidence thresholds.
The package is lightweight and written entirely in base R, without additional dependencies.


## Colab notebook
A plug-and-play [Google Colab](https://colab.research.google.com/github/badonyi/acmgscaler/blob/main/acmgscaler.ipynb) notebook with a simple interface is available for all users.


## Installation
You can install the stable version of `acmgscaler` directly from GitHub using the `devtools` package:

```r
# install devtools if you haven't already
install.packages('devtools')

# install the acmgscaler package from GitHub
devtools::install_github('badonyi/acmgscaler')
```


## Quick start

```r
library(acmgscaler)
data('variant_data', package = 'acmgscaler')

# calibrate example data
calib <- calibrate(
  df = variant_data,
  value = 'score', 
  prior = 0.1,
  group = 'gene'
)

# likelihood_ratio data for each variant
calib$BRCA1$likelihood_ratios

# score thresholds for ACMG/AMP evidence levels
calib$BRCA1$score_thresholds
```


The Colab uses a non-exported internal function to display the score intervals:

```r
acmgscaler:::prettify_score_thresholds(calib$BRCA1$score_thresholds)
```

```
  ACMG/AMP evidence strength      score
1          Benign-VeryStrong       <NA>
2              Benign-Strong > -0.61243
3            Benign-Moderate > -0.77564
4          Benign-Supporting > -0.84782
5      Pathogenic-Supporting < -1.01955
6        Pathogenic-Moderate < -1.08209
7          Pathogenic-Strong < -1.17951
8      Pathogenic-VeryStrong < -2.14538
```

## How to cite `acmgscaler`
If you find this package useful, refer to Badonyi & Marsh, *bioRxiv* (2025); doi: [10.1101/2025.05.16.654507v2](https://www.biorxiv.org/content/10.1101/2025.05.16.654507v2)
