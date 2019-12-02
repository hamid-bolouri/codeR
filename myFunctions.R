FPKM2TPM <- function(fpkm) {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
#TPM = FPKM / (sum of FPKM over all genes/transcripts) * 10^6


#------------------------------- tSNE -----------------------------------------#
library(Rtsne)

tsneTbl <- # matrix of samples X features

# Simple filtering below. More generally, use library(irlba)
# truncated SVD as in scATAC-seq below.
M <- apply(tsneTbl, 2, max)
tsneTbl <- tsneTbl[ , which(M >= 5)]

MAD <- unlist(apply(tsneTbl, 2, mad))
tsneTbl <- tsneTbl[ , which(MAD > 0.5)]

pees = 10*(1:10); allTSNE = list(); i =0
for (P in pees) { # testing different perplexity para values
	i <- i + 1
	tsneRes <- Rtsne(tsneTbl, perplexity=P) 
	allTSNE[[i]] <- tsneRes
	plot(tsneRes$Y[ , 1], tsneRes$Y[ , 2], cex=0.5, 
		 col=rgb(0,0.5,1,0.5), pch=19, xaxt='n', 
		 yaxt='n', main=P, xlab="", ylab="")
	print(P)
}
dev.off()


#------------------------------ scATAC-seq tSNE -------------------------------#
library(Matrix)
library(Rtsne)
library(irlba)
# https://github.com/shendurelab/mouse-atac/blob/master/dim_reduction/dim_reduction.R
atac_dim_reduction = function(atac_matrix, site_frequency_threshold=0.03) {
               num_cells_ncounted = Matrix::rowSums(atac_matrix)
               threshold = ncol(atac_matrix) * site_frequency_threshold
			   
			   print(threshold)

               ncounts = atac_matrix[num_cells_ncounted >= threshold,]

               ## Normalize the data with TF-IDF
			   ## term frequency?inverse document frequency _penalizes_ common entries
			   ## N.B. scATAC-seq specific!
               nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
               tf_idf_counts = nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts))

               ## Do SVD
               set.seed(0)
               SVDtsne = irlba(tf_idf_counts, 50, 50, maxit=1000) # 50 right/left Singular Vectors
               d_diagtsne = matrix(0, nrow=length(SVDtsne$d), ncol=length(SVDtsne$d))
               diag(d_diagtsne) = SVDtsne$d
               SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v)) # SVD reduced matrix
               rownames(SVDtsne_vd) = colnames(atac_matrix)
               colnames(SVDtsne_vd) = paste0('pca_', 1:ncol(SVDtsne_vd))

               ## Run TSNE to 2 dimensions
               tsnetfidf = Rtsne(SVDtsne_vd, pca=F, perplexity=30, max_iter=5000)

               tsne_coords = as.data.frame(tsnetfidf$Y)
               colnames(tsne_coords) = c('tsne_1', 'tsne_2')
               rownames(tsne_coords) = colnames(ncounts)

               return(list('tsne_coords'=tsne_coords, 'pca_coords'=as.data.frame(SVDtsne_vd)))
}


#---------------------------------------- UMAP --------------------------------#

library(umap)
X <- # matrix of samples X features
umapRes <- umap(X, min_dist=2, n_neighbors=30)
pdf("xxx.pdf")
plot(umapRes$layout, xaxt='n', yaxt='n', xlab="UMAP1", 
	 ylab="UMAP2", col="gray", main="Coral=<3, blue=MLLT10")
points(umapRes$layout[LT3, ], col=rgb(0,0,0,0.5), pch=19)
legend("bottomright", legend=c("<3"),
		col=c(rgb(0,0,0,0.5)), pch=c(19), bty='n')
dev.off()


#----------------------------------------- PCA --------------------------------#
X <- # matrix of samples X features
PCA <- prcomp(X)
varPC1 <- signif(100*summary(PCA)$importance['Proportion of Variance', 1], 3)
varPC2 <- signif(100*summary(PCA)$importance['Proportion of Variance', 2], 3)
pdf("PCA_IF_DE.miRs.pdf")
plot(PCA$x, col=rgb(0,0.5,1,0.2), pch=19, xaxt='n', yaxt='n')
points(PCA$x[grp1, 1], PCA$x[grp1, 2], pch=19, col="coral")
dev.off()


#----------------------------------------- NMF --------------------------------#

TBL <- # matrix of features X samples

library(NMF)

## From NMFTutorial:
resSVD <- svd(TBL)
H = hist(resSVD$d, breaks=50, plot=FALSE) 
sum(H$counts[H$mids > 200]) # If 5 SVD diag mat singular vals are outliers => K=5

## From NMF Vignette:
## WARNING: NMF will use ALL AVAILABLE CPUs !!!!!
Sys.time()
estimRank <- nmf(TBL, rank=(K-1):(K+1), nrun=100, seed=123)
save(estimRank, file="NMF_bySample_estimRank.RData")
Sys.time()

pdf("NMF.rankEsimtationLinePlot.pdf")
plot(estimRank)
dev.off()

pdf("NMF.rankEsimtationHeatmap.pdf")
consensusmap(estimRank, labCol=NA, labRow=NA) # annCol=colAnnot, 
dev.off()

Sys.time()
resNN <- nmf(TBL, rank=5, seed='nndsvd')
resICA <- nmf(TBL, rank=5, seed='ica')
Sys.time()

scoresNN <- extractFeatures(resNN)
scoresICA <- extractFeatures(resICA)

save(resNN, file="resNN_NMF.Rdata")
save(resICA, file="resICA_NMF.Rdata")

resNN <- get(load("resNN_NMF.Rdata"))
resICA <- get(load("resICA_NMF.Rdata"))

pdf("basisMaps_NN.ICA_rankX.pdf")
basismap(resNN, color="YlGnBu:50", main="nndsvd")
basismap(resICA, color="YlGnBu:50", main="ICA")
dev.off()

# Assuming selected RANK == 5
pdf("basisMaps_estimRank_rank5s.pdf")
basismap(estimRank$fit$'5', color="YlGnBu:50", main="estimRank=5")
dev.off()


#----------------------------- Convert CPM to TPM -----------------------------#
# http://luisvalesilva.com/datasimple/rna-seq_units.html
# Start with feature summarization data (y == an edgeR DGEList object)
y <- edgeR::calcNormFactors(y)

# Calculate TPM from normalized counts
calc_tpm <- function(x, gene.length) {
  x <- as.matrix(x)
  len.norm.lib.size <- colSums(x / gene.length)
  return((t(t(x) / len.norm.lib.size) * 1e06) / gene.length)
}

tpm_man <- calc_tpm(y, gene.length = y$genes$Length)

