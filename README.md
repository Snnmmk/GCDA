# CSEGA
Cell-type Specific Enrichment by Gene-based Annotation
# geneEnrichR

A R function to evaluate whether certain cell types are enriched in the High vs Low expression group of a target gene in Single-Cell data.

## 🔧 Usage

```r
source("CSEGA.R")
result <- gene_enrich_analysis(sce_data, "celltype", "gene_name")
result$summary
draw(result$heatmap)
