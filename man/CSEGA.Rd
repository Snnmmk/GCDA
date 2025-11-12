% Please edit documentation in R/CSEGA.R
\name{CSEGA}
\alias{CSEGA}
\title{Cell-type Specific Expression-Guided Analysis}
\usage{
CSEGA(seu, celltype_var, gene, threshold_method = "kmeans")
}
\arguments{
\item{seu}{A Seurat object.}

\item{group_var}{Metadata column name (e.g., "celltype" ).}

\item{gene}{Target gene symbol (e.g., "DDIT3").}

\item{threshold_method}{Method to split expression into High/Low ("kmeans"(default),"learn","","median", "q0.xx" or a numeric cutoff).}
}
\value{
A list with:
  \item{summary}{Data frame of proportions, fold changes, odds ratios, logistic β values, and adjusted p-values.}
  \item{heatmap}{A ComplexHeatmap object visualizing cell-type composition differences.}
  \item{bubbleplot}{A ggplot2 bubble plot showing High/Low proportions and logistic β significance.}
  \item{adaptive_cutoff}{The applied or learned cutoff threshold.}
  \item{models}{List of per-celltype logistic regression model objects.}
}
\description{
Given a Seurat object, a metadata column (cell type or cluster), and a gene name,
this function compares the composition of each cell type between High and Low
expression groups of the specified gene. High/Low groups are defined by the
specified threshold method, optionally learned automatically. Logistic regression
is fitted per cell type to quantify the effect size (β) and significance of
gene expression on the probability of belonging to that cell type.
}
