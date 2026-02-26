#' Compute PCA for Quality Control
#'
#' Performs Principal Component Analysis on the expression matrix.
#'
#' @param expr_matrix A numeric matrix with features (rows) and samples (columns).
#'
#' @return A list with:
#'   - pca: The prcomp object
#'   - scores: Data frame of PC scores with sample names
#'   - var_explained: Vector of variance explained by each PC (in %)
#'
#' @details
#' Uses prcomp() with scaling enabled. Features are standardized before PCA.
#'
#' @examples
#' \dontrun{
#'   data_list <- prepare_analysis_data(eset)
#'   pca_result <- compute_pca(data_list$expr_matrix)
#'   head(pca_result$scores)
#' }
#'
#' @importFrom stats prcomp
#'
#' @export
compute_pca <- function(expr_matrix) {
  if (!is.numeric(expr_matrix) || !is.matrix(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix")
  }

  if (nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) {
    stop("Expression matrix must have at least 2 features and 2 samples for PCA")
  }

  pca_obj <- stats::prcomp(t(expr_matrix), scale. = TRUE)

  # Prepare scores dataframe with sample names
  scores_df <- as.data.frame(pca_obj$x)
  if (!is.null(rownames(expr_matrix))) {
    rownames(scores_df) <- rownames(expr_matrix)
  }

  # Calculate variance explained (%)
  var_explained <- (pca_obj$sdev^2 / sum(pca_obj$sdev^2)) * 100

  list(
    pca = pca_obj,
    scores = scores_df,
    var_explained = var_explained
  )
}

#' Calculate Sample Distribution Statistics
#'
#' Computes summary statistics for each sample to assess data distributions.
#'
#' @param expr_matrix A numeric matrix with features (rows) and samples (columns).
#'
#' @return A data frame with rows as samples and columns:
#'   - sample_name: Sample identifier
#'   - mean: Mean expression
#'   - median: Median expression
#'   - sd: Standard deviation
#'   - q1: First quartile
#'   - q3: Third quartile
#'   - min: Minimum expression
#'   - max: Maximum expression
#'   - n_genes: Number of expressed features (non-NA)
#'
#' @examples
#' \dontrun{
#'   data_list <- prepare_analysis_data(eset)
#'   dist_stats <- calculate_sample_distributions(data_list$expr_matrix)
#'   head(dist_stats)
#' }
#'
#' @importFrom stats median sd quantile
#'
#' @export
calculate_sample_distributions <- function(expr_matrix) {
  if (!is.numeric(expr_matrix) || !is.matrix(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix")
  }

  sample_stats <- data.frame(
    sample_name = colnames(expr_matrix),
    mean = colMeans(expr_matrix, na.rm = TRUE),
    median = apply(expr_matrix, 2, stats::median, na.rm = TRUE),
    sd = apply(expr_matrix, 2, stats::sd, na.rm = TRUE),
    q1 = apply(expr_matrix, 2, stats::quantile, probs = 0.25, na.rm = TRUE),
    q3 = apply(expr_matrix, 2, stats::quantile, probs = 0.75, na.rm = TRUE),
    min = apply(expr_matrix, 2, min, na.rm = TRUE),
    max = apply(expr_matrix, 2, max, na.rm = TRUE),
    n_genes = colSums(!is.na(expr_matrix)),
    row.names = NULL
  )

  sample_stats
}

#' Run Differential Expression Analysis with Limma
#'
#' Performs limma-based differential expression analysis between two groups.
#'
#' @param eset An ExpressionSet object.
#' @param phenotype_col A character string specifying the phenotype column name in pData(eset).
#' @param group1 A character string specifying the first group value.
#' @param group2 A character string specifying the second group value.
#'
#' @return A data frame with rows as genes and columns:
#'   - gene_id: Feature/gene name
#'   - logFC: Log2 fold change (group2 vs group1)
#'   - AveExpr: Average expression across all samples
#'   - t: t-statistic
#'   - P.Value: Raw p-value
#'   - adj.P.Val: Adjusted p-value (BH)
#'   - significance: Character indicator ("***", "**", "*", or "")
#'
#' @details
#' Uses limma's standard workflow:
#'   1. Design matrix from phenotype column
#'   2. Model fit with lmFit
#'   3. Contrast fit comparing group2 vs group1
#'   4. eBayes moderation
#'   5. Results extracted and annotated
#'
#' @examples
#' \dontrun{
#'   eset <- fetch_geo_dataset("GSE63310")
#'   de_results <- run_differential_expression(
#'     eset,
#'     phenotype_col = "treatment",
#'     group1 = "control",
#'     group2 = "treated"
#'   )
#'   head(de_results)
#' }
#'
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @importFrom Biobase exprs pData featureNames
#'
#' @export
run_differential_expression <- function(eset, phenotype_col, group1, group2) {
  # Validation
  if (!inherits(eset, "ExpressionSet")) {
    stop("eset must be an ExpressionSet object")
  }

  pdata <- Biobase::pData(eset)

  if (!(phenotype_col %in% colnames(pdata))) {
    stop("phenotype_col '", phenotype_col, "' not found in pData")
  }

  pheno_values <- pdata[[phenotype_col]]

  if (!(group1 %in% pheno_values) || !(group2 %in% pheno_values)) {
    stop("group1 and/or group2 not found in phenotype column '", phenotype_col, "'")
  }

  # Filter to only samples in selected groups
  sample_indices <- which(pheno_values %in% c(group1, group2))
  eset_subset <- eset[, sample_indices]
  pheno_subset <- pheno_values[sample_indices]

  # Design matrix
  design <- stats::model.matrix(~0 + factor(pheno_subset))
  colnames(design) <- c("group1", "group2")

  # Fit the model
  expr_matrix <- Biobase::exprs(eset_subset)
  fit <- limma::lmFit(expr_matrix, design)

  # Create contrast (group2 - group1)
  contrast_matrix <- limma::makeContrasts(group2 - group1, levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  # Get results
  de_table <- limma::topTable(fit2, n = Inf, sort.by = "p")
  de_table$gene_id <- rownames(de_table)
  de_table <- de_table[c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

  # Add significance column
  de_table$significance <- ""
  de_table$significance[de_table$adj.P.Val < 0.001] <- "***"
  de_table$significance[de_table$adj.P.Val < 0.01 & de_table$adj.P.Val >= 0.001] <- "**"
  de_table$significance[de_table$adj.P.Val < 0.05 & de_table$adj.P.Val >= 0.01] <- "*"

  rownames(de_table) <- NULL
  de_table
}
