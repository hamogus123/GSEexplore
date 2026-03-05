# Generate manifest.json for Posit Connect Cloud deployment
# This script maps Bioconductor repositories explicitly to avoid resolution errors

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("rsconnect", quietly = TRUE)) install.packages("rsconnect")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 1. Capture the standard Bioconductor repositories
repos <- BiocManager::repositories()

# 2. MANUALLY INJECT THE MAPPING
# This forces the scanner to recognize the string "Bioconductor 3.21"
repos["Bioconductor 3.21"] <- "https://bioconductor.org/packages/3.21/bioc"
repos["BioCann"] <- "https://bioconductor.org/packages/3.21/data/annotation"

# 3. Apply these repositories to the global options
options(repos = repos)

# 4. Write the manifest.json
# This file will now contain the explicit URLs for all Bioconductor dependencies
rsconnect::writeManifest(
  appDir = ".",
  appPrimaryDoc = "app.R"
)

cat("manifest.json generated successfully!\n")
