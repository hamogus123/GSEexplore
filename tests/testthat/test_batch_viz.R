test_that("detect_batch_effects identifies batch structure", {
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  phenotypes <- data.frame(
    treatment = c(rep("control", 5), rep("treated", 5)),
    row.names = paste0("Sample", 1:10)
  )

  batch_info <- detect_batch_effects(expr_matrix, phenotypes, "treatment")

  expect_true(is.list(batch_info))
  expect_true(all(c("n_sv", "sv_matrix", "has_batch") %in% names(batch_info)))
  expect_true(is.numeric(batch_info$n_sv))
  expect_true(is.logical(batch_info$has_batch))
})

test_that("detect_batch_effects validates inputs", {
  expect_error(
    detect_batch_effects(data.frame(), NULL, "col"),
    "must be a numeric matrix"
  )

  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  phenotypes <- data.frame(treatment = rep("A", 10))

  expect_error(
    detect_batch_effects(expr_matrix, phenotypes, "nonexistent"),
    "not found in phenotypes"
  )
})

test_that("correct_batch_effects returns corrected matrix", {
  # ComBat_seq requires count data (non-negative integers)
  set.seed(42)
  expr_matrix <- matrix(rpois(500, lambda = 10), nrow = 50, ncol = 10)
  rownames(expr_matrix) <- paste0("Gene", 1:50)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  batch <- c(rep("batch1", 5), rep("batch2", 5))

  corrected <- correct_batch_effects(expr_matrix, batch)

  expect_equal(dim(corrected), dim(expr_matrix))
  expect_equal(colnames(corrected), colnames(expr_matrix))
  expect_equal(rownames(corrected), rownames(expr_matrix))
  expect_true(all(corrected >= 0))  # ComBat_seq preserves non-negative values
})

test_that("correct_batch_effects validates inputs", {
  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  batch_wrong_length <- c("a", "b", "c")

  expect_error(
    correct_batch_effects(expr_matrix, batch_wrong_length),
    "batch length must equal"
  )
})

test_that("plot_expression_heatmap creates heatmap", {
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr_matrix) <- paste0("Gene", 1:50)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  phenotypes <- data.frame(
    treatment = c(rep("control", 5), rep("treated", 5)),
    row.names = paste0("Sample", 1:10)
  )

  # Test without gene selection
  p <- plot_expression_heatmap(expr_matrix, phenotypes)
  expect_true(is.list(p))

  # Test with gene indices
  gene_idx <- 1:20
  p2 <- plot_expression_heatmap(expr_matrix, phenotypes, gene_indices = gene_idx)
  expect_true(is.list(p2))
})

test_that("plot_expression_heatmap validates inputs", {
  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  phenotypes <- data.frame(group = rep("A", 10))

  expect_error(
    plot_expression_heatmap(expr_matrix, phenotypes, gene_indices = c(0, 100)),
    "out of bounds"
  )
})

test_that("get_top_variable_genes selects most variable genes", {
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr_matrix) <- paste0("Gene", 1:50)

  top_genes <- get_top_variable_genes(expr_matrix, top_n = 20)

  expect_equal(nrow(top_genes), 20)
  expect_true(all(c("gene_id", "variance", "row_index") %in% colnames(top_genes)))
  expect_true(all(top_genes$variance > 0))
  # Variance should be in descending order
  expect_true(all(diff(top_genes$variance) <= 0))
})

test_that("get_top_variable_genes handles edge cases", {
  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr_matrix) <- paste0("Gene", 1:10)

  # Request more genes than available
  top_genes <- get_top_variable_genes(expr_matrix, top_n = 1000)
  expect_equal(nrow(top_genes), 10)

  # Custom feature names
  custom_names <- paste0("Feature", 1:10)
  top_genes2 <- get_top_variable_genes(expr_matrix, feature_names = custom_names, top_n = 5)
  expect_equal(nrow(top_genes2), 5)
  expect_true(all(startsWith(top_genes2$gene_id, "Feature")))
})
