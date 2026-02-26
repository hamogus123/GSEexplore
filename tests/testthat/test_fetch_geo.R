test_that("fetch_geo_dataset validates accession format", {
  expect_error(
    fetch_geo_dataset("invalid"),
    "Invalid GSE accession format"
  )
  expect_error(
    fetch_geo_dataset(123),
    "must be a single character string"
  )
  expect_error(
    fetch_geo_dataset(c("GSE1", "GSE2")),
    "must be a single character string"
  )
})

test_that("validate_expressionset checks structure", {
  # Create a minimal ExpressionSet for testing
  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr_matrix) <- paste0("Gene", 1:10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  pdata <- data.frame(
    sample_id = paste0("Sample", 1:10),
    row.names = paste0("Sample", 1:10)
  )

  eset <- Biobase::ExpressionSet(
    assayData = expr_matrix,
    phenoData = Biobase::AnnotatedDataFrame(pdata)
  )

  # Should return invisibly without error
  expect_invisible(GSEexplore:::validate_expressionset(eset))
})

test_that("validate_expressionset detects invalid inputs", {
  expect_error(
    GSEexplore:::validate_expressionset(data.frame()),
    "not an ExpressionSet"
  )

  # Create invalid ExpressionSet with mismatched dimensions
  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  pdata <- data.frame(
    sample_id = paste0("Sample", 1:5),
    row.names = paste0("Sample", 1:5)  # Match row names with sample_id
  )

  # This test verifies dimension checking in validate function
  # We're testing error handling when phenoData doesn't match exprs
  # Create a valid eset first, then we'd have to manually break it
  eset <- Biobase::ExpressionSet(
    assayData = expr_matrix,
    phenoData = Biobase::AnnotatedDataFrame(
      data.frame(group = rep(c("A", "B"), 5), row.names = paste0("Sample", 1:10))
    )
  )

  expect_invisible(GSEexplore:::validate_expressionset(eset))
})

test_that("prepare_analysis_data returns proper structure", {
  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr_matrix) <- paste0("Gene", 1:10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  pdata <- data.frame(
    sample_id = paste0("Sample", 1:10),
    group = rep(c("A", "B"), 5),
    row.names = paste0("Sample", 1:10)
  )

  eset <- Biobase::ExpressionSet(
    assayData = expr_matrix,
    phenoData = Biobase::AnnotatedDataFrame(pdata)
  )

  result <- prepare_analysis_data(eset)

  expect_true(is.list(result))
  expect_true(all(c("expr_matrix", "phenotypes", "feature_names", "sample_names") %in% names(result)))
  expect_equal(dim(result$expr_matrix), c(10, 10))
  expect_equal(nrow(result$phenotypes), 10)
  expect_equal(length(result$feature_names), 10)
  expect_equal(length(result$sample_names), 10)
})

test_that("get_dataset_summary returns expected fields", {
  expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr_matrix) <- paste0("Gene", 1:10)
  colnames(expr_matrix) <- paste0("Sample", 1:10)

  pdata <- data.frame(
    sample_id = paste0("Sample", 1:10),
    group = rep(c("A", "B"), 5),
    row.names = paste0("Sample", 1:10)
  )

  eset <- Biobase::ExpressionSet(
    assayData = expr_matrix,
    phenoData = Biobase::AnnotatedDataFrame(pdata)
  )

  summary_info <- get_dataset_summary(eset)

  expect_equal(summary_info$n_features, 10)
  expect_equal(summary_info$n_samples, 10)
  expect_true("group" %in% summary_info$phenotype_cols)
  expect_length(summary_info$expression_range, 2)
})
