#' GSEexplore Package
#'
#' An interactive Shiny application for exploring Gene Expression Omnibus (GEO)
#' datasets with quality control and differential expression analysis capabilities.
#'
#' @description
#' GSEexplore provides a user-friendly interface for:
#'   - Fetching GEO datasets by accession ID
#'   - Quality control assessment (PCA, sample distributions)
#'   - Batch effect detection and correction
#'   - Heatmap visualization of expression patterns
#'   - Differential expression analysis using limma
#'   - Interactive data exploration
#'
#' @keywords internal
#'
#' @import shiny
#' @import Biobase
#' @import limma
#' @import sva
#' @import pheatmap
#' @import ggplot2
#' @import plotly
#' @import DT
#' @importFrom GEOquery getGEO
#' @importFrom stats prcomp median sd quantile model.matrix var
#' @importFrom utils write.csv
#'
"_PACKAGE"
