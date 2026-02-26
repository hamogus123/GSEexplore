#' Detect Batch Effects Using SVA
#'
#' Detects potential batch effects in the expression data using the sva package.
#'
#' @param expr_matrix A numeric matrix with features (rows) and samples (columns).
#' @param phenotypes A data frame of sample metadata with at least one column representing groups.
#' @param phenotype_col A character string specifying the column name for biological groups.
#'
#' @return A list with:
#'   - n_sv: Number of surrogate variables detected
#'   - sv_matrix: Matrix of surrogate variables (if n_sv > 0)
#'   - has_batch: Logical indicating if batch effects detected
#'
#' @details
#' Uses the ComBat-Seq methodology to estimate the number of surrogate variables
#' representing batch or other unwanted variation.
#'
#' @examples
#' \dontrun{
#'   data_list <- prepare_analysis_data(eset)
#'   batch_info <- detect_batch_effects(
#'     data_list$expr_matrix,
#'     data_list$phenotypes,
#'     "treatment"
#'   )
#'   cat("Detected surrogate variables:", batch_info$n_sv, "\n")
#' }
#'
#' @importFrom sva num.sv
#'
#' @export
detect_batch_effects <- function(expr_matrix, phenotypes, phenotype_col) {
  if (!is.numeric(expr_matrix) || !is.matrix(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix")
  }

  if (nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) {
    stop("Expression matrix must have at least 2 features and 2 samples")
  }

  if (!(phenotype_col %in% colnames(phenotypes))) {
    stop("phenotype_col not found in phenotypes")
  }

  # Create design matrix from phenotype
  design <- stats::model.matrix(~factor(phenotypes[[phenotype_col]]))

  # Estimate number of surrogate variables
  n_sv <- tryCatch(
    {
      sva::num.sv(expr_matrix, design)
    },
    error = function(e) {
      0
    }
  )

  sv_matrix <- NULL
  if (n_sv > 0) {
    sv_obj <- tryCatch(
      {
        sva::svaseq(expr_matrix, design, design)
      },
      error = function(e) {
        NULL
      }
    )
    if (!is.null(sv_obj)) {
      sv_matrix <- sv_obj$sv
    }
  }

  list(
    n_sv = n_sv,
    sv_matrix = sv_matrix,
    has_batch = n_sv > 0
  )
}

#' Correct for Batch Effects with ComBat-Seq
#'
#' Corrects for batch effects in RNA-seq count data using ComBat-Seq.
#'
#' @param expr_matrix A numeric matrix with features (rows) and samples (columns).
#' @param batch A factor or character vector specifying batch for each sample.
#'
#' @return A corrected expression matrix with same dimensions and dimnames as input.
#'
#' @details
#' Uses the sva::ComBat_seq function which is specifically designed for RNA-seq counts.
#' For microarray data, use ComBat (not implemented here).
#'
#' @examples
#' \dontrun{
#'   data_list <- prepare_analysis_data(eset)
#'   batch_variable <- Biobase::pData(eset)$batch_column
#'   corrected <- correct_batch_effects(data_list$expr_matrix, batch_variable)
#'   head(corrected)
#' }
#'
#' @importFrom sva ComBat_seq
#'
#' @export
correct_batch_effects <- function(expr_matrix, batch) {
  if (!is.numeric(expr_matrix) || !is.matrix(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix")
  }

  if (length(batch) != ncol(expr_matrix)) {
    stop("batch length must equal number of samples (columns)")
  }

  # Ensure batch is a factor
  batch <- as.factor(batch)

  # Apply ComBat-Seq correction
  corrected <- tryCatch(
    {
      sva::ComBat_seq(expr_matrix, batch = batch, group = NULL)
    },
    error = function(e) {
      warning("Batch correction failed, returning original matrix: ", conditionMessage(e))
      expr_matrix
    }
  )

  corrected
}

#' Create Expression Heatmap for Top Genes
#'
#' Generates a heatmap of top differentially expressed genes or selected genes.
#'
#' @param expr_matrix A numeric matrix with features (rows) and samples (columns).
#' @param phenotypes A data frame of sample metadata.
#' @param gene_indices Integer vector of row indices to display (optional).
#' @param color_by Character string specifying phenotype column for sample coloring (optional).
#' @param main Character string for plot title.
#' @param scale Character specifying scaling: "row", "column", or "none".
#'
#' @return A pheatmap object (invisibly) with a rendered heatmap.
#'
#' @details
#' Uses pheatmap for hierarchical clustering and heatmap display.
#' If gene_indices not provided, uses all genes.
#'
#' @examples
#' \dontrun{
#'   data_list <- prepare_analysis_data(eset)
#'   de_results <- run_differential_expression(eset, "treatment", "control", "treated")
#'   # Top 50 genes
#'   top_indices <- order(de_results$adj.P.Val)[1:50]
#'   plot_expression_heatmap(
#'     data_list$expr_matrix,
#'     data_list$phenotypes,
#'     gene_indices = top_indices,
#'     color_by = "treatment",
#'     main = "Top DE Genes"
#'   )
#' }
#'
#' @importFrom pheatmap pheatmap
#'
#' @export
plot_expression_heatmap <- function(expr_matrix, phenotypes, gene_indices = NULL,
                                    color_by = NULL, main = "Expression Heatmap",
                                    scale = "row") {
  if (!is.numeric(expr_matrix) || !is.matrix(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix")
  }

  if (is.null(gene_indices)) {
    gene_indices <- seq_len(nrow(expr_matrix))
  }

  if (any(gene_indices < 1) || any(gene_indices > nrow(expr_matrix))) {
    stop("gene_indices out of bounds")
  }

  # Subset matrix
  expr_subset <- expr_matrix[gene_indices, ]

  # Create annotation data frame for columns
  annotation_col <- NULL
  if (!is.null(color_by) && color_by %in% colnames(phenotypes)) {
    annotation_col <- data.frame(
      Group = phenotypes[[color_by]],
      row.names = colnames(expr_matrix)
    )
  }

  # Create heatmap
  p <- pheatmap::pheatmap(
    expr_subset,
    scale = scale,
    annotation_col = annotation_col,
    main = main,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    silent = TRUE
  )

  invisible(p)
}

#' Create Clustered Gene Expression Heatmap
#'
#' Creates an interactive heatmap with sample and gene clustering.
#'
#' @param expr_matrix A numeric matrix with features (rows) and samples (columns).
#' @param feature_names Character vector of feature/gene names (rownames by default).
#' @param top_n Integer specifying number of most variable genes to include.
#'
#' @return A data frame with variance information for top genes.
#'
#' @details
#' Filters to the top N most variable genes, useful for large datasets.
#' Returns variance statistics for selected genes.
#'
#' @examples
#' \dontrun{
#'   data_list <- prepare_analysis_data(eset)
#'   top_var_genes <- get_top_variable_genes(data_list$expr_matrix, top_n = 100)
#'   head(top_var_genes)
#' }
#'
#' @export
get_top_variable_genes <- function(expr_matrix, feature_names = NULL, top_n = 100) {
  if (!is.numeric(expr_matrix) || !is.matrix(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix")
  }

  if (is.null(feature_names)) {
    feature_names <- rownames(expr_matrix)
  }

  if (length(feature_names) != nrow(expr_matrix)) {
    stop("feature_names length must match number of rows")
  }

  # Calculate variance for each gene
  gene_vars <- apply(expr_matrix, 1, stats::var, na.rm = TRUE)

  # Handle NA variances
  gene_vars[is.na(gene_vars)] <- 0

  # Get top N
  top_n <- min(top_n, nrow(expr_matrix))
  top_indices <- order(gene_vars, decreasing = TRUE)[1:top_n]

  # Create results data frame
  result_df <- data.frame(
    gene_id = feature_names[top_indices],
    variance = gene_vars[top_indices],
    row_index = top_indices,
    row.names = NULL
  )

  result_df[order(result_df$variance, decreasing = TRUE), ]
}
