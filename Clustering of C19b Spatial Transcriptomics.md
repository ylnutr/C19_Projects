---
title: 'Clustering of C19b Spatial Transcriptomics'
author: "Yujin L"
date: "`r Sys.Date()`"
output: md_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 0. Background and Information

## 0-1. Project Backgound

* Early onset colorectal cancer (EOCRC) is increasing yet the mechanism for the accelerated shift to a younger age is unknown. We hypothesize that sporadic gene mutations that are necessary for colon carcinogenesis but not sufficient to cause cancer, create a preneoplastic cell state for increased risk of transformation.

## 0-2. Method

* To identify transitional cell phenotypes, we used Spatial Transcriptomics (ST) on distal colon samples from Car1-Cre (CAC) mice with Cre expression limited to ~6% of the distal colon/rectum epithelium, where EOCRC is more prevalent.  

  **(A) Mice:** Car1-Cre (CAC) transgenic mice express Cre recombinase under the control of the mouse Car1 gene promoter (Xue, Fleet 2010 Mol Cancer Res.8:1095). This limits Cre expression to the distal colon and rectum, where EOCRC is more prevalent. CAC mice were crossed to those with floxed alleles to generate mice with one (Apc+/-, ACR) or two (Apc+/-; KrasG12D, AKCR) non-transforming, oncogene mutations. 
  **(B) Sample information and tissue processing:** Distal colon and tumors with phenotypes ranging from benign to malignant were collected from 10-wk CR, 15-20 wk ACR, or 4 wk and 8-10 wk AKCR mice. Tumors and Swiss Rolls (with the distal colon in the center) were embedded in Optimal Cutting Temperature compound (OCT) on dry ice and stored at −80 °C. OCT blocks were cut with a pre-cooled cryostat at 10 um thickness, and sections were transferred onto 6.5 mm2 oligo-barcoded capture areas on the Visium 10x Genomics slide. 
  **(C) 10X Genomics Visium Spatial Gene Expression System:** Tissue sections on oligo-barcoded capture slides were briefly fixed and stained with hematoxylin and eosin (H&E) and imaged. The Visium Spatial Gene Expression protocol was used to process the slides. The permeabilization enzymes were added to the slides for RNA capture, followed by reverse transcription for cDNA synthesis. cDNA libraries were constructed and sequenced on an Illumina NextSeq 500. Sequencing yielded data for 15,277 capture areas under tissue, with 84,554 mean reads and an average of 3,583 genes per capture area.
  
## 0-3. Set working directory

In this step, we verify and update the working directory to ensure that all subsequent file operations (e.g., reading or writing data) are performed in the correct project folder. This is particularly important for maintaining reproducibility and avoiding file path errors. The directory is set to the folder containing the results of the functional annotation analysis for spatial transcriptomics clusters.
  
```{r setwd, echo=TRUE, warning=TRUE}
getwd()
setwd("/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/05.Analysis_in_R/3.Clustering/")
getwd()
```
  
# 1. Load packages

* Besides Seurat, we need to import some useful packages.

* dplyr and ggplot2 are part of tidyverse, a readable and flexible R language for solving data science challenges. I personally prefer the coding style with tidyverse, but you may use base R too. patchwork combines separate ggplots into the same graphic with easy access to control layouts. limma is a Bioconductor package to analyze microarray data. It is not essential for our analysis, but may provide a more efficient implementation of the Wilcoxon rank sum test in differential expression analysis.

```{r library, echo=TRUE, warning=TRUE}
library(tidyverse)
library(Seurat)
library(dplyr) # data manipulation
library(ggplot2)
library(patchwork)
library(clustree)
library(scran)
library(viridis)
library(ggforce)
library(gghalves)
library(ggridges)
library(future)
```

## 1-1. Configure Parallel Processing and Memory Allocation

To improve performance and avoid memory-related errors during computationally intensive steps (e.g., SCTransform()), we configure the future package to:

  * Use multisession parallelism, allowing tasks to run on multiple CPU cores.
  * Increase the maximum allowed size for global variables to 4 GB, which is necessary for handling large Seurat objects during parallelized operations.
  
```{r multisession, echo=TRUE, warning=TRUE}
plan("multisession") # Enables parallel processing across multiple cores
options(future.globals.maxSize = 10 * 1024^3) # Increase memory limit to 10 GB
```

# 2. Load Data

## 2-1. Load 12 10X Genomics Visium dataset into Seurat
A spatial gene expression dataset of mouse distal colon and tumor tissue collected by Space Ranger 2.0.0. will
be analyzed throughout the tutorial. Both the gene expression matrix and spatial imaging data are
necessary for the computational analysis.
The data files we will be using today include:

  * a (filtered) feature / cell matrix HDF5 file (.h5)
  * a spatial imaging data folder (.tar.gz)

```{r load, echo=TRUE, warning=TRUE}
E296_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/E296_C-A1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "CR_1", filter.matrix = TRUE, to.upper = FALSE)
C507_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C507_A-A1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "ACR_1", filter.matrix = TRUE, to.upper = FALSE)
C569_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C569_C-B1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "ACR_2", filter.matrix = TRUE, to.upper = FALSE)
C822_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C822_C-C1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "AKCR_1", filter.matrix = TRUE, to.upper = FALSE)
C498_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C498_A-B1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "AKCR_2", filter.matrix = TRUE, to.upper = FALSE)
C531_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C531_A-D1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "AKCR_3", filter.matrix = TRUE, to.upper = FALSE)
C583_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C583_A-C1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "AKCR_4", filter.matrix = TRUE, to.upper = FALSE)
C520_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C520_B-A1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "AKCR_5", filter.matrix = TRUE, to.upper = FALSE)
C522_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C522_B-B1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "AKCR_6", filter.matrix = TRUE, to.upper = FALSE)
C565_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C565_C-D1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "ACR_TU1", filter.matrix = TRUE, to.upper = FALSE)
C515_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C515_B-C1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "ACR_TU2", filter.matrix = TRUE, to.upper = FALSE)
C533_data <- Seurat::Load10X_Spatial(
  data.dir = "/stor/work/Fleet/BCG/2024_04.Fleet_visium_customgenome/02.count/C533_B-D1/outs/",
  filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "AKCR_TU1", filter.matrix = TRUE, to.upper = FALSE)
```


## 2-2. View each dataset

```{r dataview, echo=TRUE, warning=TRUE}
E296_data
C507_data
C569_data
C822_data
C498_data
C583_data
C531_data
C520_data
C522_data
C565_data
C515_data
C533_data
```
# 3. Preprocess data

## 3-1. Quality control

A few common QC metrics include

  * The number of unique genes detected in each sample (nFeature_Spatial).
  * The total number of molecules detected within a sample (nCount_Spatial).
  * The percentage of reads that map to the mitochondrial genome (percent.mt).

```{r qcmetrics, echo=TRUE, warning=TRUE}
E296_data[["percent.mt"]] <- PercentageFeatureSet(E296_data, pattern = "^mt-")
VlnPlot(
  E296_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C507_data[["percent.mt"]] <- PercentageFeatureSet(C507_data, pattern = "^mt-")
VlnPlot(
  C507_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C569_data[["percent.mt"]] <- PercentageFeatureSet(C569_data, pattern = "^mt-")
VlnPlot(
  C569_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C822_data[["percent.mt"]] <- PercentageFeatureSet(C822_data, pattern = "^mt-")
VlnPlot(
  C822_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C498_data[["percent.mt"]] <- PercentageFeatureSet(C498_data, pattern = "^mt-")
VlnPlot(
  C498_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C531_data[["percent.mt"]] <- PercentageFeatureSet(C531_data, pattern = "^mt-")
VlnPlot(
  C531_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C583_data[["percent.mt"]] <- PercentageFeatureSet(C583_data, pattern = "^mt-")
VlnPlot(
  C583_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C520_data[["percent.mt"]] <- PercentageFeatureSet(C520_data, pattern = "^mt-")
VlnPlot(
  C520_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C522_data[["percent.mt"]] <- PercentageFeatureSet(C522_data, pattern = "^mt-")
VlnPlot(
  C522_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C565_data[["percent.mt"]] <- PercentageFeatureSet(C565_data, pattern = "^mt-")
VlnPlot(
  C565_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C515_data[["percent.mt"]] <- PercentageFeatureSet(C515_data, pattern = "^mt-")
VlnPlot(
  C515_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
C533_data[["percent.mt"]] <- PercentageFeatureSet(C533_data, pattern = "^mt-")
VlnPlot(
  C533_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```
## 3-2. Subsetting samples based on observation

After observing these plots from each sample, I decided to filter out molecules < 1500 & a percentage of reads that map to the mitochondrial genome < 20%.


```{r subsetting, echo=TRUE, warning=TRUE}
E296_subset <- subset(E296_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(E296_data) - ncol(E296_subset), 
            "samples because of the outlier QC metrics, with", ncol(E296_subset), "samples left."))
C507_subset <- subset(C507_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C507_data) - ncol(C507_subset), 
            "samples because of the outlier QC metrics, with", ncol(C507_subset), "samples left."))
C569_subset <- subset(C569_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C569_data) - ncol(C569_subset), 
            "samples because of the outlier QC metrics, with", ncol(C569_subset), "samples left."))
C822_subset <- subset(C822_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C822_data) - ncol(C822_subset), 
            "samples because of the outlier QC metrics, with", ncol(C822_subset), "samples left."))
C498_subset <- subset(C498_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C498_data) - ncol(C498_subset), 
            "samples because of the outlier QC metrics, with", ncol(C498_subset), "samples left."))
C531_subset <- subset(C531_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C531_data) - ncol(C531_subset), 
            "samples because of the outlier QC metrics, with", ncol(C531_subset), "samples left."))
C583_subset <- subset(C583_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C583_data) - ncol(C583_subset), 
            "samples because of the outlier QC metrics, with", ncol(C583_subset), "samples left."))
C520_subset <- subset(C520_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C520_data) - ncol(C520_subset), 
            "samples because of the outlier QC metrics, with", ncol(C520_subset), "samples left."))
C522_subset <- subset(C522_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C522_data) - ncol(C522_subset), 
            "samples because of the outlier QC metrics, with", ncol(C522_subset), "samples left."))
C565_subset <- subset(C565_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C565_data) - ncol(C565_subset), 
            "samples because of the outlier QC metrics, with", ncol(C565_subset), "samples left."))
C515_subset <- subset(C515_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C515_data) - ncol(C515_subset), 
            "samples because of the outlier QC metrics, with", ncol(C515_subset), "samples left."))
C533_subset <- subset(C533_data, subset = nCount_Spatial > 1500 & percent.mt < 20)
print(paste("Filter out", ncol(C533_data) - ncol(C533_subset), 
            "samples because of the outlier QC metrics, with", ncol(C533_subset), "samples left."))
```

## 3-3. Gene Expression Visualization

The SpatialFeaturePlot() function in Seurat extends FeaturePlot(), and can overlay molecular data on top of tissue histology. I used this feature to see the subsetted out spots overlayed on tissue image.

```{r postqcview,fig.height=8, echo=TRUE, warning=TRUE}
plot <- SpatialFeaturePlot(
  E296_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("E296_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/E296_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C507_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C507_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C507_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C569_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C569_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C569_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C822_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C822_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C822_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C498_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C498_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C498_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C531_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C531_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C531_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C583_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C583_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C583_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C520_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C520_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C520_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C522_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C522_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C522_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C565_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C565_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C565_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C515_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C515_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C515_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
plot <- SpatialFeaturePlot(
  C533_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "right") & ggtitle("C533_subset")
plot
png(filename = "../3.Clustering/Subsetted plots/C533_subset.png", height = 2000, width = 6000, res = 300)
print(plot)
dev.off()
```
We can see that the filtering successfully discarded spots that appeared as low quality.

# 4. Normalization

## 4-1. Relabel orig.ident

The variance of molecular counts expresses spatial heterogeneity which cannot be solely explained by technical noise. Satija Lab and Collaborators recommends normalization using SCTransform (Hafemeister and Satija, 2019) in order to account for technical bias while preserving true biological differences. We will also change the "orig.ident" of each sample to unique sample number (genotype+number, e.g. AKCR_3). By changing their orig.ident to unique names, we can still call up a single sample when they are pooled.

```{r changeorigident, echo=TRUE, warning=TRUE}
E296_subset$orig.ident <- "CR_1"
C507_subset$orig.ident <- "ACR_1"
C569_subset$orig.ident <- "ACR_2"
C822_subset$orig.ident <- "AKCR_1"
C498_subset$orig.ident <- "AKCR_2"
C531_subset$orig.ident <- "AKCR_3"
C583_subset$orig.ident <- "AKCR_4"
C520_subset$orig.ident <- "AKCR_5"
C522_subset$orig.ident <- "AKCR_6"
C565_subset$orig.ident <- "ACR_TU1"
C515_subset$orig.ident <- "ACR_TU2"
C533_subset$orig.ident <- "AKCR_TU1"
```

## 4-2. Merge Spatial Samples and Visualize Pre-Normalization Metrics

In this step, we merge multiple subsetted spatial transcriptomics datasets into a single Seurat object using the merge() function. This pooled object allows for integrated downstream analysis across all samples. We also define a specific order for the orig.ident factor to maintain consistent sample labeling in visualizations.

To assess the quality and distribution of the merged data before normalization, we generate violin plots of two key metrics:

  * nCount_Spatial: Total number of transcripts per spot
  * nFeature_Spatial: Number of detected genes per spot

These plots help identify variability across samples and provide a baseline for evaluating the impact of normalization.

```{r poolspatialsamples, echo=TRUE, warning=TRUE}
spatial_pooled <- merge(x = E296_subset, y = c(C507_subset, C569_subset, C822_subset, C498_subset, C531_subset, C583_subset, C520_subset, C522_subset, C565_subset, C515_subset, C533_subset))
order <- c("CR_1","ACR_1","ACR_2","AKCR_1","AKCR_2","AKCR_3","AKCR_4","AKCR_5", "AKCR_6", "ACR_TU1", "ACR_TU2", "AKCR_TU1")
spatial_pooled[["orig.ident"]] <- factor(x = spatial_pooled@meta.data$orig.ident, levels = order)
p<- VlnPlot(spatial_pooled, group.by = c("orig.ident"), features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 1)
p
png(filename = "../3.Clustering/presctransform.png", height = 2350, width = 1500, res = 300)
print(p)
dev.off()
```
## 4-3. Join Spatial Layers

This step uses the JoinLayers() function to combine multiple spatial layers within the merged Seurat object (spatial_pooled). This is useful when working with spatial transcriptomics data that includes multiple tissue slices or capture areas, allowing for unified downstream analysis.

```{r join layer, echo=TRUE, warning=TRUE}
spatial_pooled
spatial_joined <- JoinLayers(spatial_pooled)
spatial_joined
```
## 4-4. Normalize and Identify 3,000 Variable Features

Here, we apply SCTransform() to the joined spatial object to normalize the data and identify the top 3,000 variable features. This normalization method models technical noise while preserving biological variation, which is essential for accurate clustering and dimensionality reduction.

```{r 3000genes, echo=TRUE, warning=TRUE}
spatial_3000 <- SCTransform(spatial_joined, assay = "Spatial", ncells = 14799, residual.features = NULL, variable.features.n = 3000,  verbose = TRUE)
```

Now, we can check the impact of normalization by generating violin plots of two key metrics:

```{r checkpost, echo=TRUE, warning=TRUE}
p<- VlnPlot(spatial_3000, group.by = c("orig.ident"), features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 1)
p
png(filename = "../3.Clustering/postsctransform.png", height = 2350, width = 1500, res = 300)
print(p)
dev.off()
```

# 5. Clustering

## 5-1. Clustering and UMAP Visualization

This chunk performs the core steps of clustering and dimensionality reduction on the normalized dataset (spatial_3000) using the top 20 principal components:

  * FindNeighbors(): Constructs a shared nearest neighbor (SNN) graph based on PCA-reduced data.
  * FindClusters(): Detects transcriptionally similar groups (clusters) using a graph-based modularity optimization algorithm.
  * RunUMAP(): Projects the data into a 2D space for intuitive visualization of the clusters.

These steps are essential for identifying and visualizing distinct cell populations in the spatial transcriptomics data.

```{r 3000clustering, echo=TRUE, warning=TRUE}
spatial_3000 <- RunPCA(spatial_3000, verbose = FALSE)
spatial_3000 <- FindNeighbors(spatial_3000, dims = 1:20)
spatial_3000 <- FindClusters(spatial_3000, verbose = FALSE)
spatial_3000 <- RunUMAP(spatial_3000, dims = 1:20)
```
## 5-2. Select Clustering Resolution and Visualize UMAP by Sample

In this step, we finalize the clustering resolution by setting it to 0.7, which was likely chosen based on prior exploration (e.g., clustree analysis). This resolution determines the granularity of the clusters identified in the dataset.

The resulting clusters are then visualized using UMAP, with cells grouped by both their cluster assignment (SCT_snn_res.0.7) and original sample identity (orig.ident). Labels are added to enhance interpretability.

This visualization helps assess how well the clusters are distributed across samples and whether they reflect meaningful biological patterns.

```{r resolution, echo=TRUE, warning=TRUE}
spatial_3000 <- FindClusters(spatial_3000, assay = "SCT", resolution = 0.7)
plot <- DimPlot(spatial_3000, reduction = "umap", group.by = c("SCT_snn_res.0.7", "orig.ident"), label = TRUE) 
plot
```

We then visualize the clusters using UMAP, splitting the plot by sample (orig.ident) to assess how well the clusters are distributed across different tissue sections. This helps evaluate whether clustering is driven by biological variation rather than technical artifacts.

```{r origumap, echo=TRUE, warning=TRUE}
plot <- DimPlot(spatial_3000, reduction = "umap", group.by = "SCT_snn_res.0.7", split.by = "orig.ident", label = TRUE, ncol = 3, label.size = 3) 
plot
```
This chunk generates a spatial UMAP plot using SpatialDimPlot() to visualize the spatial distribution of clusters identified at resolution 0.7. Each spot is colored by its assigned cluster, and labels are added for clarity. The plot is split across samples and arranged in a grid layout, allowing for easy comparison of spatial patterns across tissue sections.

This visualization is particularly useful for assessing how transcriptionally defined clusters map onto the physical tissue architecture.

```{r spatialumapindividualsamplessave, echo=TRUE, warning=TRUE}
plot <- SpatialDimPlot(spatial_3000, group.by = "SCT_snn_res.0.7", label = TRUE, label.size = 3, ncol = 3)
plot
```

I will explain all the metrics and variables I put into consideration when deciding on this specific clustering in the separate Rmd file.