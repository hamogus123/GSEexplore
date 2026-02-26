# Batch Correction & Heatmap Feature Guide

## Quick Start

### Batch Effect Detection

```r
library(GSEexplore)

# Load data
eset <- fetch_geo_dataset("GSE63310")
data <- prepare_analysis_data(eset)

# Detect batch effects
batch_info <- detect_batch_effects(
  data$expr_matrix,
  data$phenotypes,
  phenotype_col = "treatment"  # biological grouping
)

if (batch_info$has_batch) {
  cat("Found", batch_info$n_sv, "batch variables\n")
}
```

### Batch Correction

```r
# Correct for batch (requires count data)
batch_var <- data$phenotypes$batch_column
corrected <- correct_batch_effects(data$expr_matrix, batch_var)

# Use corrected matrix for downstream analysis
```

### Heatmap Visualization

```r
# Top DE genes heatmap
if (!is.null(de_results)) {
  top_de <- order(de_results$adj.P.Val)[1:50]
  plot_expression_heatmap(
    data$expr_matrix,
    data$phenotypes,
    gene_indices = top_de,
    color_by = "treatment",
    main = "Top 50 DE Genes"
  )
}

# Most variable genes
top_var <- get_top_variable_genes(data$expr_matrix, top_n = 100)
plot_expression_heatmap(
  data$expr_matrix,
  data$phenotypes,
  gene_indices = top_var$row_index,
  color_by = "treatment",
  main = "100 Most Variable Genes"
)
```

## Interactive Shiny App Usage

### Batch Correction Tab

1. **Load Dataset** (Data Loading tab) with GSE accession
2. **Select Batch Variable** from phenotype columns
3. **Click "Check for Batch Effects"**
4. If batch detected:
   - Review "Surrogate Variables" count
   - Click "Correct Batch Effects" to apply ComBat-Seq
   - Status shows "[Corrected]" when complete
5. Corrected data automatically used in downstream analyses

### Heatmap Tab

1. **Gene Source Selection:**
   - "Top DE genes": Uses DE results (requires running DE analysis first)
   - "Most variable": Selects genes with highest variance
   - "All genes": Uses all genes up to the count limit

2. **Number of Genes:** Slider 10-500

3. **Color Samples By:** Select any phenotype column for annotations

4. **Generate Heatmap:** Creates clustered heatmap with:
   - Hierarchical clustering (correlation-based)
   - Row scaling for standardized colors
   - Sample phenotype annotations
   - Gene/sample names

## Function Reference

### Batch Detection
```r
detect_batch_effects(expr_matrix, phenotypes, phenotype_col)
# Returns: list(n_sv, sv_matrix, has_batch)
```

### Batch Correction
```r
correct_batch_effects(expr_matrix, batch)
# Returns: corrected expression matrix (same dimensions)
```

### Heatmap Plotting
```r
plot_expression_heatmap(
  expr_matrix,
  phenotypes,
  gene_indices = NULL,      # Use all genes if NULL
  color_by = NULL,          # Phenotype column for annotation
  main = "Expression Heatmap",
  scale = "row"             # "row", "column", or "none"
)
# Returns: pheatmap object (invisibly)
```

### Variable Gene Selection
```r
get_top_variable_genes(expr_matrix, feature_names = NULL, top_n = 100)
# Returns: data.frame with columns: gene_id, variance, row_index
```

## Example Workflow

```r
# Complete analysis pipeline
library(GSEexplore)

# 1. Fetch data
eset <- fetch_geo_dataset("GSE63310")
data <- prepare_analysis_data(eset)

# 2. QC: Check distributions
dist_stats <- calculate_sample_distributions(data$expr_matrix)
print(summary(dist_stats$mean))

# 3. QC: PCA
pca_result <- compute_pca(data$expr_matrix)
plot(pca_result$var_explained[1:10])

# 4. Batch assessment
batch_info <- detect_batch_effects(
  data$expr_matrix, data$phenotypes, "batch"
)
if (batch_info$has_batch) {
  data$expr_matrix <- correct_batch_effects(
    data$expr_matrix, data$phenotypes$batch
  )
  cat("Applied batch correction\n")
}

# 5. DE analysis
de_results <- run_differential_expression(
  eset, "treatment", "control", "treated"
)

# 6. Visualize results
plot_expression_heatmap(
  data$expr_matrix,
  data$phenotypes,
  gene_indices = order(de_results$adj.P.Val)[1:100],
  color_by = "treatment",
  main = "Top 100 DE Genes"
)
```

## Troubleshooting

**Issue:** Batch correction fails with "Negative counts not allowed"
- **Solution:** ComBat_seq requires count data. Use only with RNA-seq counts, not transformed/normalized data

**Issue:** Heatmap takes a long time with all genes
- **Solution:** Use "Most variable" option or limit gene count with the slider

**Issue:** No batch effects detected but data looks unusual
- **Solution:** Check for confounding variables; batch column might not capture the main source of variation

**Issue:** Heatmap colors all look similar
- **Solution:** Try different scale options in the code (column scale, no scale) or use a subset of genes

## References

- **Surrogate Variable Analysis:** https://www.bioconductor.org/packages/sva/
- **ComBat-Seq:** Zhang et al. (2020) Nature Communications
- **pheatmap:** https://CRAN.R-project.org/package=pheatmap
