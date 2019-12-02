#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
   stop("need to specify names of filtered_feature_bc_matrix dir, CITE-seq-count_res file, & output ID", 
        call.=FALSE)
} 

# HB 1st July 2019
# Run this script with, for example:
# Rscript --vanilla cellHash.R /home/hamid.bolouri/__Data/cellHashing/cellRanger_miseq/filtered_feature_bc_matrix /home/hamid.bolouri/__Data/cellHashing/CITEseqCount_res_D703.csv testHash &


## Allocate memory  and CPU usage

library(future)
plan(strategy = "multiprocess", workers = 8)
options(future.globals.maxSize = 20 * 1024 ^ 3)

library("future.apply")
library("stats")


#############################################################
## https://satijalab.org/seurat/v3.0/hashing_vignette.html ##
#############################################################

library(Seurat)
library(Matrix)

#### Get command-line arguments ####
filtered_feature_bc_matrix = args[1]
CITEseqCount_res = args[2]
runID = args[3] # e.g. == "test" for a test run

print(filtered_feature_bc_matrix)
print(CITEseqCount_res)
print(runID)

# Load in the UMI matrix
pbmc.umis <- Read10X(filtered_feature_bc_matrix)

# Load in the HTO count matrix #### NOW FROM COMMAND-LINE INPUT ####
setwd("~/__Data/cellHashing")
pbmc.htos <- read.csv(CITEseqCount_res, row.names=1, as.is=TRUE)
pbmc.htos <- pbmc.htos[!grepl("bad_struct|no_match|total_reads", rownames(pbmc.htos)), ]

# Select cell barcodes detected by both RNA and HTO:
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos)) # 8962

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Check that the HTO have the correct names
# rownames(pbmc.htos)

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, using centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

# For very large dataset we suggest using k_function = 'clara'
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99,
                         seed = 123)

# Global classification results
doubletRateTbl <- table(pbmc.hashtag$HTO_classification.global)
# miseq test:
# Doublet Negative  Singlet 
#    1138       16     7808 
setwd("~/__Data/cellHashing")
save(doubletRateTbl, file=paste0(runID, "_doubletRate.RData"))


# Group cells based on the max HTO signal
Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:2], ncol = 2)

# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
setwd("~/__Data/cellHashing")
pdf(paste0(runID, "_featureScatterPlot.pdf"))
FeatureScatter(pbmc.hashtag, feature1 = "hto_POP1", feature2 = "hto_POP2")
dev.off()

# Compare number of UMIs for singlets, doublets and negative cells
Idents(pbmc.hashtag) <- "HTO_classification.global"
setwd("~/__Data/cellHashing")
pdf(paste0(runID, "_nUMIs_doubletsEtcVlnPlot.pdf"))
VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()


# Generate a tSNE embedding for HTOs. Here, grouping cells by singlets and doublets

# First, we will remove negative cells from the object
pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = pbmc.hashtag.subset, assay = "HTO"))))

# Calculate tSNE embeddings using above distance matrix
pbmc.hashtag.subset <- RunTSNE(pbmc.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
setwd("~/__Data/cellHashing")
pdf(paste0(runID, "_hashtag_tSNE.pdf"))
DimPlot(pbmc.hashtag.subset)
dev.off()

# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
setwd("~/__Data/cellHashing")
pdf(paste0(runID, "_heatmapHTO.pdf"))
HTOHeatmap(pbmc.hashtag, assay = "HTO", ncells = 5000)
dev.off()

# Cluster and visualize cells using the usual scRNA-seq workflow

# Extract the singlets
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

# Select the top 1000 most variable features
pbmc.singlet <- FindVariableFeatures(pbmc.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
pbmc.singlet <- ScaleData(pbmc.singlet, features = VariableFeatures(pbmc.singlet))


# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
pbmc.singlet <- FindNeighbors(pbmc.singlet, reduction = "pca", dims = 1:10)
pbmc.singlet <- FindClusters(pbmc.singlet, resolution = 0.6)

pbmc.singlet <- RunTSNE(pbmc.singlet, reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
setwd("~/__Data/cellHashing")
pdf(paste0(runID, "_tSNE_clusters.pdf"))
DimPlot(pbmc.singlet, group.by = "HTO_classification")
dev.off()



