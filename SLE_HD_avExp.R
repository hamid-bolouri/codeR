library(future)
plan(strategy = "multiprocess", workers = 8)    ############## set number of CPUs #################
options(future.globals.maxSize = 20 * 1024 ^ 3)

library("future.apply")
library("stats")


library("Matrix")
library('dplyr')
library('Seurat')


setwd("/home/hamid.bolouri/__Data/Aldan.SLE.HD")
normalizedPBMC <- get(load("SeuratNormalizedCountsALL_GRCh38_filteredMT.RData"))

# setwd("~/__Data/shelfLife_LSE.v.HD")
# annot <- read.csv("library_metadata.csv", as.is=TRUE)
setwd("~/__Data/shelfLifeSLE/")
meta <- get(load("Aggr_Human_Immunology_RTX-850-856_GRCh38-3.0.0_premrna_samp.dat_v3_reads.RData"))

tags <- unlist(lapply(sapply(meta[ ,1], function(x) strsplit(x, "BRISL")), function(y) y[1])) # 276,505

iIDs <- match(tags, colnames(normalizedPBMC))
normalizedPBMC <- normalizedPBMC[ , iIDs[which(!is.na(iIDs))]] # 24,801 x 271,635
colnames(normalizedPBMC) <- names(tags)[which(!is.na(iIDs))] # loss of ~ 500 cells. AOK.

setwd("~/__Data/Aldan.SLE.HD/clustering/kmeans_10_clusters")
CLs <- read.csv("clusters.csv", header=TRUE, as.is=TRUE)

iIDs <- match(tags, CLs$Barcode)
CLs <- CLs[iIDs[which(!is.na(iIDs))], ]
CLs$Barcode <- names(tags)[which(!is.na(iIDs))] # 274,797 x 2

iIDs <- match(CLs$Barcode, colnames(normalizedPBMC))
normalizedPBMC <- normalizedPBMC[ , iIDs[which(!is.na(iIDs))]] # 24,801 x 27,1635 # AOK
CLs <- CLs[iIDs[which(!is.na(iIDs))], ] # 271,635 x 2 # AOK


IDs = c("BRISL2", "BRISL3", "BRISL4", "BRISL5", "BRISL6", "BRISL7")
HRs = c("-2h$|-2hr$", "-4h$|-4hr$", "-6h$|-6hr$", "-8h$|-8hr$", "-18h$|-18hr$")


clusterNos = 1:10

Sys.time()
avExpList <- list(); indx = 0
for (i in IDs) {
   print(i)
   i.ID <- grep(i, colnames(normalizedPBMC))
   for (cl in clusterNos) {
      i.CL <- grep(cl, CLs$Cluster)
      resMtrx = colID = NULL
      for (h in HRs) {
         i.HR <- grep(h, colnames(normalizedPBMC))
         j <- intersect(intersect(i.ID, i.HR), i.CL)
         if (length(j) > 2) { # skip clusters with few cells in this ID+HR
            resMtrx <- cbind(resMtrx, apply(normalizedPBMC[ , j], 1, median))
            colID <- c(colID, paste0(i, ":", "C", cl, ":", h))
         }
      }
      if (!is.null(resMtrx)) {
         colnames(resMtrx) <- colID
         indx <- indx + 1
         avExpList[[indx]] <- resMtrx
         names(avExpList)[indx] <- paste0(i, ":", "C", cl)
      }
   }
}
Sys.time()

# selGenes <- lapply(avExpList, 
#                    function(x) which(unlist(apply(x, 1, max)) > 1))
# for (i in 1:length(selGenes)) {
#    avExpList[[i]] <- avExpList[[i]][selGenes[[i]], ]
# }

setwd("~/__Data")
save(avExpList, file="avExpList_Aldan.GRCh38.SLE.HD_noFilter.RData")


