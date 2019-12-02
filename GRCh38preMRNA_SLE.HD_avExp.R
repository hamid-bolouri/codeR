library(future)
plan(strategy = "multiprocess", workers = 16)    ############## set number of CPUs #################
options(future.globals.maxSize = 20 * 1024 ^ 3)

library("future.apply")
library("stats")


library("Matrix")
library('dplyr')
library('Seurat')


setwd("/home/hamid.bolouri/__Data/")
pbmc <- get(load("Aggr_Human_Immunology_RTX-850-856_GRCh38-3.0.0_premrna_mat_v3_reads.RData"))

setwd("~/__Data/shelfLifeSLE/")
meta <- get(load("Aggr_Human_Immunology_RTX-850-856_GRCh38-3.0.0_premrna_samp.dat_v3_reads.RData"))

tags <- unlist(lapply(sapply(meta[ ,1], function(x) strsplit(x, "BRISL")), function(y) y[1])) # 276,505

iIDs <- match(names(tags), colnames(pbmc))
pbmc <- pbmc[ , iIDs[which(!is.na(iIDs))]] # 33538 x 76505

setwd("~/__Data/")
CLs <- read.csv("Aggr_Human_Immunology_RTX-850-856_GRCh38-3.0.0_premrna_clustersK8.csv", header=TRUE, as.is=TRUE)

iIDs <- match(tags, CLs$Barcode)
CLs <- CLs[iIDs[which(!is.na(iIDs))], ]
CLs$Barcode <- names(tags)[which(!is.na(iIDs))] # 276505 x 2

iIDs <- match(CLs$Barcode, colnames(pbmc))
pbmc <- pbmc[ , iIDs[which(!is.na(iIDs))]] # 33538 x 276505 # AOK
CLs <- CLs[iIDs[which(!is.na(iIDs))], ] # 276505 x 2 # AOK


IDs = c("BRISL2", "BRISL3", "BRISL4", "BRISL5", "BRISL6", "BRISL7")
HRs = c("-2h$|-2hr$", "-4h$|-4hr$", "-6h$|-6hr$", "-8h$|-8hr$", "-18h$|-18hr$")


clusterNos = 1:8

Sys.time()
avExpList <- list(); indx = 0
for (i in IDs) {
   print(i)
   i.ID <- grep(i, colnames(pbmc))
   for (cl in clusterNos) {
      i.CL <- grep(cl, CLs$Cluster)
      resMtrx = colID = NULL
      for (h in HRs) {
         i.HR <- grep(h, colnames(pbmc))
         j <- intersect(intersect(i.ID, i.HR), i.CL)
         if (length(j) > 2) { # skip clusters with few cells in this ID+HR
            resMtrx <- cbind(resMtrx, apply(pbmc[ , j], 1, median))
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
save(avExpList, file="avExpList_JeffGoldy.GRCh38.preMRNA.SLE.HD_noFilter.RData")

# To plot average expression of a gene for just (say) HD samples, do:
avHD <- avExpList[c("BRISL2:C6", "BRISL3:C6", "BRISL4:C6")]
selGene <- colSums(do.call(rbind, lapply(avHD, function(x) x["NFKBIA", ])))
names(selGene) <- c("2h", "4h", "6h", "8h", "18h")
plot(selGene)




