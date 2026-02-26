# GSEexplore: Interactive GEO Expression Data Explorer

An R Shiny application for interactive exploration of Gene Expression Omnibus (GEO) datasets with quality control assessment and differential expression analysis.

## Features

- **Data Fetching**: Download GEO datasets directly by GSE accession ID
- **Quality Control**: 
  - Principal Component Analysis (PCA) for sample relationships
  - Sample distribution analysis with summary statistics
- **Batch Effect Detection & Correction**:
  - Surrogate variable analysis (SVA)
  - ComBat-Seq batch correction
- **Visualization**:
  - Hierarchical clustered heatmaps with sample/gene annotations
  - Top DE genes or most variable genes
- **Differential Expression Analysis**:
  - Limma-based statistical testing
  - Volcano plot visualization
  - Downloadable results
- **Interactive Exploration**:
  - Gene expression viewer
  - Phenotype-based sample coloring
  - Interactive plots with plotly
  - Sortable/searchable results tables

## Installation

```r
# Install dependencies
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("Biobase", "GEOquery", "limma"))
install.packages(c("shiny", "ggplot2", "plotly", "DT", "golem", "config"))

# Install GSEexplore
devtools::install_github("user/GSEexplore")
```

## Quick Start

### Interactive Shiny App

```r
library(GSEexplore)
run_app()
```

Then:
1. Enter a GSE accession (e.g., "GSE63310")
2. Click "Fetch Dataset"
3. Explore QC views, run differential expression, browse genes

### Programmatic Usage

```r
library(GSEexplore)

# Fetch dataset
eset <- fetch_geo_dataset("GSE63310")

# Prepare data
data_list <- prepare_analysis_data(eset)

# Quality control
pca_result <- compute_pca(data_list$expr_matrix)
dist_stats <- calculate_sample_distributions(data_list$expr_matrix)

# Differential expression
de_results <- run_differential_expression(
  eset,
  phenotype_col = "treatment",
  group1 = "control",
  group2 = "treated"
)

# View results
head(de_results)
```

## Core Functions

### Data Handling

- `fetch_geo_dataset(accession)` - Download GEO dataset by accession
- `prepare_analysis_data(eset)` - Extract and format data
- `get_dataset_summary(eset)` - Get dataset overview

### Quality Control

- `compute_pca(expr_matrix)` - Principal Component Analysis
- `calculate_sample_distributions(expr_matrix)` - Sample statistics

### Batch Effects

- `detect_batch_effects(expr_matrix, phenotypes, phenotype_col)` - Surrogate variable analysis
- `correct_batch_effects(expr_matrix, batch)` - ComBat-Seq batch correction

### Visualization

- `plot_expression_heatmap(expr_matrix, phenotypes, gene_indices, color_by)` - Clustered heatmaps
- `get_top_variable_genes(expr_matrix, top_n)` - Select most variable genes

### Analysis

- `run_differential_expression(eset, phenotype_col, group1, group2)` - Limma-based DE testing

## Documentation

- **Vignette**: `vignette("gse63310_walkthrough")` - Complete walkthrough with GSE63310
- **Function Help**: `?function_name` for detailed documentation on each function

## Example Dataset

The package vignette uses **GSE63310**, a microarray study with:
- ~48,000 probes
- Multiple treatment groups
- Well-documented phenotype data

Perfect for learning the full workflow.

## Package Structure

```
GSEexplore/
├── R/
│   ├── app_ui.R              # Shiny UI
│   ├── app_server.R          # Server logic
│   ├── fetch_geo.R           # Data fetching
│   ├── analysis.R            # QC and DE functions
│   └── run_app.R             # App launcher
├── tests/testthat/           # Unit tests
├── vignettes/                # Documentation
├── DESCRIPTION               # Package metadata
└── README.md                 # This file
```

## Dependencies

**Required:**
- shiny (≥ 1.7.0)
- Biobase
- GEOquery
- limma
- sva (batch effect detection/correction)
- pheatmap (heatmap visualization)
- ggplot2
- plotly
- DT
- golem (≥ 0.3.0)
- config

**Testing:**
- testthat (≥ 3.0.0)

## Testing

Run tests with:

```r
devtools::test()
```

Tests cover:
- Data fetching and validation
- Analysis function outputs
- Edge cases and error handling

## Citation

If you use GSEexplore in your research, please cite:

```
GSEexplore: Interactive GEO Expression Data Explorer (2026)
https://github.com/user/GSEexplore
```

## License

MIT License - see LICENSE file

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## Issues & Support

For bugs, feature requests, or questions:
- Open an issue on GitHub
- Check existing documentation and vignettes first

## Authors

GSEexplore Contributors (2026)

## References

- **GEO**: https://www.ncbi.nlm.nih.gov/geo/
- **Biobase**: https://bioconductor.org/packages/Biobase/
- **GEOquery**: https://bioconductor.org/packages/GEOquery/
- **limma**: https://bioconductor.org/packages/limma/
- **Shiny**: https://shiny.posit.co/
