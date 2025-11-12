# CSEGA
Cell-type Specific Enrichment by Gene-based Annotation

# Download CSEGA

```r
library(remotes)
remotes::install_github("Snnmmk/CSEGA")
```

# CSEGA R

A R function to evaluate whether certain cell types are enriched in the High vs Low expression group of a target gene in Single-Cell data.

## 🔧 Usage

```r
source("CSEGA.R")
result <- CSEGA(sce_data, "celltype", "gene_name", threshold_method = "median")

⚙️ Threshold options

threshold_method can be one of the following:
	•	"median" (default) — split by median expression
	•	"mean" — split by mean expression
	•	"q0.xx" — split by a given quantile (e.g. "q0.25", "q0.75")
	•	numeric value — use a fixed cutoff (e.g. 1.2)

result$summary
draw(result$heatmap)
result$bubbleplot
