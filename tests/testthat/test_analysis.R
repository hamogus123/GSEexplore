test_that("compute_pca returns correct structure", {
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  pca_result <- compute_pca(expr_matrix)

  expect_true(is.list(pca_result))
  expect_true(all(c("pca", "scores", "var_explained") %in% names(pca_result)))
  expect_equal(nrow(pca_result$scores), 10)
  expect_equal(length(pca_result$var_explained), 10)
  expect_true(all(pca_result$var_explained > 0))
  expect_true(sum(pca_result$var_explained) > 99)  # Should sum to ~100
})

test_that("compute_pca validates input", {
  expect_error(
    compute_pca(data.frame()),
    "must be a numeric matrix"
  )

  expect_error(
    compute_pca(matrix(1, nrow = 1, ncol = 1)),
    "at least 2 features and 2 samples"
  )
})

test_that("compute_pca scores have sample names", {
  expr_matrix <- matrix(rnorm(200), nrow = 20, ncol = 10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  pca_result <- compute_pca(expr_matrix)

  expect_equal(nrow(pca_result$scores), 10)
  expect_equal(rownames(pca_result$scores), colnames(expr_matrix))
})

test_that("calculate_sample_distributions returns correct structure", {
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  dist_stats <- calculate_sample_distributions(expr_matrix)

  expected_cols <- c("sample_name", "mean", "median", "sd", "q1", "q3", "min", "max", "n_genes")
  expect_equal(colnames(dist_stats), expected_cols)
  expect_equal(nrow(dist_stats), 10)
  expect_equal(dist_stats$sample_name, colnames(expr_matrix))
  expect_true(all(dist_stats$n_genes == 50))
})

test_that("calculate_sample_distributions handles missing values", {
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)
  expr_matrix[1:10, 1] <- NA

  dist_stats <- calculate_sample_distributions(expr_matrix)

  expect_equal(dist_stats$n_genes[1], 40)  # First sample has 10 NAs
  expect_equal(dist_stats$n_genes[2], 50)  # Other samples have no NAs
})

test_that("run_differential_expression returns correct structure", {
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 20)
  rownames(expr_matrix) <- paste0("Gene", 1:50)
  colnames(expr_matrix) <- paste0("Sample", 1:20)

  pdata <- data.frame(
    sample_id = paste0("Sample", 1:20),
    group = c(rep("control", 10), rep("treated", 10)),
    row.names = paste0("Sample", 1:20)
  )

  eset <- Biobase::ExpressionSet(
    assayData = expr_matrix,
    phenoData = Biobase::AnnotatedDataFrame(pdata)
  )

  de_results <- run_differential_expression(
    eset,
    phenotype_col = "group",
    group1 = "control",
    group2 = "treated"
  )

  expected_cols <- c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "significance")
  expect_equal(colnames(de_results), expected_cols)
  expect_equal(nrow(de_results), 50)
  expect_true(all(!is.na(de_results$logFC)))
})

test_that("run_differential_expression validates inputs", {
  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr_matrix) <- paste0("Gene", 1:10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  pdata <- data.frame(
    group = rep(c("A", "B"), 5),
    row.names = paste0("Sample", 1:10)
  )

  eset <- Biobase::ExpressionSet(
    assayData = expr_matrix,
    phenoData = Biobase::AnnotatedDataFrame(pdata)
  )

  # Invalid phenotype column
  expect_error(
    run_differential_expression(eset, "nonexistent", "A", "B"),
    "not found in pData"
  )

  # Invalid group values
  expect_error(
    run_differential_expression(eset, "group", "X", "Y"),
    "not found in phenotype column"
  )
})

test_that("run_differential_expression produces meaningful results", {
  # Create data with actual difference
  set.seed(42)
  control_matrix <- matrix(rnorm(250, mean = 5, sd = 1), nrow = 50, ncol = 5)
  treated_matrix <- matrix(rnorm(250, mean = 5, sd = 1), nrow = 50, ncol = 5)

  # Make first 10 genes highly expressed in treated group
  treated_matrix[1:10, ] <- treated_matrix[1:10, ] + 3

  expr_matrix <- cbind(control_matrix, treated_matrix)
  rownames(expr_matrix) <- paste0("Gene", 1:50)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  pdata <- data.frame(
    group = c(rep("control", 5), rep("treated", 5)),
    row.names = paste0("Sample", 1:10)
  )

  eset <- Biobase::ExpressionSet(
    assayData = expr_matrix,
    phenoData = Biobase::AnnotatedDataFrame(pdata)
  )

  de_results <- run_differential_expression(eset, "group", "control", "treated")

  # First 10 genes should have low p-values
  expect_true(mean(de_results$P.Value[1:10]) < 0.05)
})
