#' Fetch and Load a GEO Dataset by Accession
#'
#' Downloads a GEO dataset from NCBI using GEOquery and returns an ExpressionSet object.
#'
#' @param accession A character string specifying the GSE accession (e.g., "GSE63310").
#'   Must be a valid GEO accession ID.
#'
#' @return An ExpressionSet object containing:
#'   - exprs: Expression matrix (features × samples)
#'   - pData: Sample metadata (phenotypes)
#'   - fData: Feature metadata/annotations
#'
#' @examples
#' \dontrun{
#'   eset <- fetch_geo_dataset("GSE63310")
#'   dim(eset)  # Check dimensions
#' }
#'
#' @importFrom GEOquery getGEO
#' @importFrom Biobase ExpressionSet assayData phenoData featureData
#'
#' @export
fetch_geo_dataset <- function(accession) {
  if (!is.character(accession) || length(accession) != 1) {
    stop("accession must be a single character string")
  }

  # Validate accession format (basic check)
  if (!grepl("^GSE\\d+$", accession)) {
    stop("Invalid GSE accession format. Expected format: GSE followed by digits (e.g., GSE63310)")
  }

  tryCatch(
    {
      message("Fetching GEO dataset: ", accession)
      geo_data <- GEOquery::getGEO(accession, GSEMatrix = TRUE, AnnotGPL = FALSE)

      if (is.list(geo_data)) {
        eset <- geo_data[[1]]
      } else {
        eset <- geo_data
      }

      validate_expressionset(eset)
      return(eset)
    },
    error = function(e) {
      stop("Failed to fetch GEO dataset: ", conditionMessage(e), call. = FALSE)
    }
  )
}

#' Validate an ExpressionSet Object
#'
#' Checks the structure and integrity of an ExpressionSet object.
#'
#' @param eset An ExpressionSet object to validate.
#'
#' @return The same ExpressionSet object (invisibly) if valid.
#' @keywords internal
#'
#' @details
#' Validates:
#'   - Object is an ExpressionSet
#'   - Expression matrix has non-zero dimensions
#'   - Phenotype data exists and matches sample count
#'   - No completely empty features or samples
#'
#' @importFrom Biobase exprs pData fData
#'
validate_expressionset <- function(eset) {
  if (!inherits(eset, "ExpressionSet")) {
    stop("Object is not an ExpressionSet")
  }

  expr_matrix <- Biobase::exprs(eset)

  if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
    stop("Expression matrix has zero dimensions")
  }

  pdata <- Biobase::pData(eset)
  if (nrow(pdata) != ncol(expr_matrix)) {
    stop("Phenotype data row count does not match expression matrix column count")
  }

  if (is.null(colnames(expr_matrix)) || any(colnames(expr_matrix) == "")) {
    stop("Expression matrix has missing or invalid sample names")
  }

  invisible(eset)
}

#' Prepare Data for Analysis
#'
#' Extracts and formats expression matrix and sample metadata from an ExpressionSet.
#'
#' @param eset An ExpressionSet object.
#'
#' @return A list with elements:
#'   - expr_matrix: Numeric matrix of expression values (features × samples)
#'   - phenotypes: Data frame of sample metadata
#'   - feature_names: Character vector of feature/gene names
#'   - sample_names: Character vector of sample names
#'
#' @examples
#' \dontrun{
#'   eset <- fetch_geo_dataset("GSE63310")
#'   data_list <- prepare_analysis_data(eset)
#'   str(data_list)
#' }
#'
#' @importFrom Biobase exprs pData featureNames sampleNames
#'
#' @export
prepare_analysis_data <- function(eset) {
  validate_expressionset(eset)

  expr_matrix <- Biobase::exprs(eset)
  phenotypes <- Biobase::pData(eset)
  feature_names <- Biobase::featureNames(eset)
  sample_names <- Biobase::sampleNames(eset)

  list(
    expr_matrix = expr_matrix,
    phenotypes = phenotypes,
    feature_names = feature_names,
    sample_names = sample_names
  )
}

#' Get Dataset Summary Information
#'
#' Returns basic summary statistics for a dataset.
#'
#' @param eset An ExpressionSet object.
#'
#' @return A list with:
#'   - n_features: Number of genes/probes
#'   - n_samples: Number of samples
#'   - phenotype_cols: Column names in phenotype data
#'   - expression_range: Min and max expression values
#'
#' @examples
#' \dontrun{
#'   eset <- fetch_geo_dataset("GSE63310")
#'   summary_info <- get_dataset_summary(eset)
#'   str(summary_info)
#' }
#'
#' @importFrom Biobase exprs pData featureNames sampleNames
#'
#' @export
get_dataset_summary <- function(eset) {
  validate_expressionset(eset)

  expr_matrix <- Biobase::exprs(eset)
  phenotypes <- Biobase::pData(eset)

  list(
    n_features = nrow(expr_matrix),
    n_samples = ncol(expr_matrix),
    phenotype_cols = colnames(phenotypes),
    expression_range = c(min = min(expr_matrix, na.rm = TRUE), max = max(expr_matrix, na.rm = TRUE))
  )
}
