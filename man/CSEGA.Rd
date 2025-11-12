% Please edit documentation in R/CSEGA.R
\name{CSEGA}
\alias{CSEGA}
\title{Cell-type Specific Expression-Guided Analysis}
\usage{
CSEGA(seu, celltype_var, gene, threshold_method = "median")
}
\arguments{
\item{seu}{A Seurat object.}

\item{group_var}{Metadata column name (e.g., "celltype" ).}

\item{gene}{Target gene symbol (e.g., "DDIT3").}

\item{threshold_method}{Method to split expression into High/Low ("kmeans"(default),"learn","","median", "q0.xx" or a numeric cutoff).}
}
\value{
A list with:
  \item{summary}{Data frame of proportions, fold change, odds ratio, and adjusted p-values.}
  \item{heatmap}{A ComplexHeatmap object visualizing proportions.}
  \item{bubbleplot}{A bubble plot object visualizing proportions.}
  \item{models}{The glm of a celltype.}
}
\description{
Given a Seurat object, a metadata column (cell type or cluster), and a gene name,
this function compares the composition of each cell type between High and Low
expression groups of the specified gene.
}
