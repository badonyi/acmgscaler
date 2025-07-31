#' Example variant data
#'
#' Example data included in the package containing MAVE-derived functional
#' scores and class labels for missense variants in BRCA1 and TP53. Functional
#' scores are from high-throughput assays (Findlay et al., 2018; Giacomelli
#' et al., 2018), and class labels are P/LP and B/LB from ClinVar.
#'
#' @format A dataframe with 459 observations and 4 variables:
#' \describe{
#'   \item{gene}{Gene symbol in which the variant occurs (e.g., `BRCA1`).}
#'   \item{variant}{Missense variant, represented in single-letter amino acid
#'   notation (e.g., `L3F`).}
#'   \item{class}{Binary class label indicating pathogenicity:
#'     \itemize{
#'       \item `P`: Pathogenic
#'       \item `B`: Benign
#'     }
#'   }
#'   \item{score}{Functional assay score, typically representing the degree of
#'   functional disruption (numeric).}
#' }
#'
#' @references
#' Giacomelli et al., 2018. Mutational processes shape the landscape of TP53
#' mutations in human cancer.
#' *Nature Genetics*, 50(10), 1381–1387.
#' [DOI: 10.1038/s41588-018-0204-y](https://doi.org/10.1038/s41588-018-0204-y)
#'
#' Findlay et al., 2018. Accurate classification of BRCA1 variants with
#' saturation genome editing.
#' *Nature*, 562(7726), 217–222.
#' [DOI: 10.1038/s41586-018-0461-z](https://doi.org/10.1038/s41586-018-0461-z)
#'
#' Landrum et al., 2020. ClinVar: improvements to accessing data.
#' *Nucleic Acids Research*, 48(D1), D835–D844.
#' [DOI: 10.1093/nar/gkz972](https://doi.org/10.1093/nar/gkz972)
#'
#' @keywords datasets
"variant_data"
