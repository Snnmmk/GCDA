# CSEGA
Cell-type Specific Enrichment by Gene-based Annotation

`CSEGA` performs **cell-type composition analysis** based on the expression of a single **GENE** in single-cell data.
It divides cells into *High* and *Low* expression groups, compares the relative abundance of each cell type,
and provides visual and statistical outputs including heatmaps, bubble plots, and regression summaries.

The advanced version introduces **machine learning–based threshold learning** and **explainable logistic regression**.

---
## Key Features
- 🔹 **Automatic threshold learning**  
  - Default: *k-means (2 clusters)* for adaptive High/Low partitioning  
  - Optional `"learn"` mode scans thresholds (1–99% quantiles) to find the cutoff that maximizes group separation  

- 🔹 **Explainable logistic regression**
  - For each cell type *ct*, the model fits:  
    **logit[P(celltype = ct)] = α + β × expression**
  - β > 0 indicates enrichment, β < 0 indicates depletion.  
  - Returns β and p-values for interpretability.
 
- 🔹 **Rich visualization outputs**  
  - `ComplexHeatmap` of cell-type proportions  
  - Bubble plot comparing High vs. Low expression groups  

---

# Download CSEGA

```r
library(remotes)
remotes::install_github("Snnmmk/CSEGA")
```

# CSEGA R


## 🔧 Usage

```r
library(CSEGA)
result <- CSEGA(sce_data, "celltype", "gene_name", threshold_method = "kmeans")

⚙️ Threshold options

threshold_method can be one of the following:
	•	numeric value — use a fixed cutoff (e.g. 1.2)
	•	"median"  — split by median expression
	•	"mean" — split by mean expression
	•	"q0.xx" — split by a given quantile (e.g. "q0.25", "q0.75")
	•	"kmeans"(default) - Apply two-cluster k-means on expression values. The cluster with the higher mean is labeled High.
    •	"learn" - Perform an adaptive search across quantile thresholds (from 1% to 99%) and automatically select the cutoff that maximizes cell-type compositional divergence. This mode enables machine-learned threshold discovery.


# View summary
head(result$summary)

# Draw heatmap
draw(result$heatmap)

# Bubble plot
print(result$bubbleplot)

# the glm of a celltype 
result$models[["celltypename"]] (e.g. "Epithelial")

glm_obj <- result$models$celltypename  
summary(glm_obj)       
coef(glm_obj)


```
## Citation
If you use **CSEGA** in your research, please cite it as:

Zidong Feng (2025). *CSEGA: Cell-type Specific Expression-Guided Analysis with Adaptive and Explainable Learning*. https://github.com/Snnmmk/CSEGA
