# Plan: Fixing Bioconductor Repository Resolution in Posit Cloud

## 1. The Problems Encountered

During the deployment of the `GeneSetEnrichmentExplorer` to Posit Connect Cloud, we encountered two critical errors related to Bioconductor dependencies in the `renv.lock` file.

### Issue A: The "URL Resolution" Bug
**Error:** `repository Bioconductor 3.21 cannot be resolved to a URL.`
*   **Cause:** Modern versions of `renv` generate a repository name string called `"Bioconductor 3.21"`. The internal scanner used by Posit Connect Cloud has a hardcoded lookup table of standard repository names (like "CRAN" or "BioCsoft"). It does not recognize the versioned string `"Bioconductor 3.21"`, causing it to fail the "pre-flight" check before the build even begins.

### Issue B: The "Annotation 404" Error
**Error:** `Error downloading ... GenomeInfoDbData_1.2.14.tar.gz: status_code=404`
*   **Cause:** Bioconductor splits its packages into multiple repositories: **Software** (`/bioc/`) and **Annotation** (`/data/annotation/`). By default, deployment scanners look for all Bioconductor packages in the Software repository. Packages like `org.Hs.eg.db` and `GenomeInfoDbData` do not exist there, leading to a 404 error during the download phase.

---

## 2. The Resolution Strategy: "Surgical Lockfile Remapping"

To fix these issues without breaking the R package structure, we performed a surgical update to the `renv.lock` file using a JSON-aware R script.

### Step 1: Standardizing Repository Names
We renamed every instance of `"Bioconductor 3.21"` to **`"BioCsoft"`**. 
*   **Why:** `BioCsoft` is a "legacy" name that is universally recognized by Posit scanners. By renaming the repository, we forced the scanner to correctly map the packages to the Bioconductor software URL.

### Step 2: Explicit Annotation Mapping
We manually injected a `"Repository": "BioCann"` field into the entry for every annotation package in the `renv.lock` file.
*   **Target Packages:** `org.Hs.eg.db`, `org.Mm.eg.db`, and `GenomeInfoDbData`.
*   **Why:** This explicitly tells the Posit Cloud server: *"Do not look for this package in the software repo; go to the Annotation (BioCann) repo instead."*

### Step 3: Ensuring JSON Integrity
Because `renv.lock` is a strict JSON file, manual text replacement (like `sed` or `grep`) can easily break commas or brackets, causing further deployment crashes. We used the R package `jsonlite` to:
1. Parse the lockfile into an R list.
2. Programmatically update the `Repository` fields.
3. Write the list back to JSON format with "pretty" formatting.

---

## 3. Implementation Script
For future use or sharing, here is the R code used to perform this fix:

```r
library(jsonlite)
lock <- fromJSON("renv.lock", simplifyVector = FALSE)

# Define packages that belong in the Annotation repo
ann_pkgs <- c("org.Hs.eg.db", "org.Mm.eg.db", "GenomeInfoDbData")

for (pkg_name in names(lock$Packages)) {
  pkg <- lock$Packages[[pkg_name]]
  if (identical(pkg$Source, "Bioconductor")) {
    # 1. Map to BioCann if it's an annotation package, otherwise BioCsoft
    repo <- if (pkg_name %in% ann_pkgs) "BioCann" else "BioCsoft"
    lock$Packages[[pkg_name]]$Repository <- repo
  }
}

# 2. Update the header mapping
# Ensure the top-level Repositories list uses the names we just assigned
lock$R$Repositories <- list(
  list(Name = "BioCsoft", URL = "https://bioconductor.org/packages/3.21/bioc"),
  list(Name = "BioCann", URL = "https://bioconductor.org/packages/3.21/data/annotation"),
  list(Name = "CRAN", URL = "https://cloud.r-project.org")
)

write(toJSON(lock, auto_unbox = TRUE, pretty = TRUE), "renv.lock")
```

---

## 4. Final Deployment Outcome
By standardizing these names, the Posit Publisher scanner passed the "pre-flight" check, and the server was able to download and compile all Bioconductor dependencies successfully.
