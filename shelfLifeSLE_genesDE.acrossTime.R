
### Allocate memory  and use all CPUs

library(future)
plan(strategy = "multiprocess", workers = 64)
options(future.globals.maxSize = 20 * 1024 ^ 3)

library("future.apply")
library("stats")


library("Matrix")
library("dplyr")
library("Seurat")
library("erer")
library("purrr")


setwd("~/__Data/shelfLifeSLE")
normalizedPBMC <- get(load("SeuratNormalizedCountsSLE_GRCh38_filteredMT.RData"))
PBMC.IDs <- unlist(lapply(strsplit(colnames(normalizedPBMC), "-"), function(x) x[1]))
 
setwd("/home/hamid.bolouri/__Data/shelfLifeSLE/outs/analysis/clustering/kmeans_10_clusters/")
CLs <- read.csv("clusters.csv", as.is=TRUE)
clusterIDs <- unlist(lapply(strsplit(CLs$Barcode, "-"), function(x) x[1]))
CLs <- CLs[match(PBMC.IDs, clusterIDs), ]
CLs$Barcode <- colnames(normalizedPBMC)

colnames(CLs) <- c("ID", "group")
CLs$group<- as.factor(CLs$group)
library(S4Vectors)
CLs <- DataFrame(CLs)

IDs = c("BRISL5", "BRISL6", "BRISL7")
HRs = c("-2hr", "-4hr", "-6hr", "-8hr", "-18hr")


##################################################################################################
####         per Adam's 27Sept2019 request for expression values as well as p-vals            ####
##################################################################################################

runDE2 <- function(time1, time2) {
   
   testDE2 <- function(x) {
      x1 <- x[grep(paste0("-", time1, "hr$"), names(x))]
      x2 <- x[grep(paste0("-", time2, "hr$"), names(x))]
      x1 <- x1[x1 > 0]
      x2 <- x2[x2 > 0]
      # Consider only clutsers with > 100 cells
      if (length(x1) > 100 & length(x2)>100) {
         pW <- wilcox.test(x1, x2)$p.value
         resTuplet <- data.frame(T1=as.integer(time1), T2=as.integer(time2), 
                                 medianT1=median(x1), medianT2=median(x2), 
                                 meanT1=mean(x1), meanT2=mean(x2), 
                                 maxT1=max(x1), maxT2=max(x2), P=pW)
         return(resTuplet)
      } else return(NULL)
   }
   
   adjDF <- function(DF) {
      tmpL <- t(apply(DF, 1, function(x) {
         x["P"] <- x["P"] * nrow(clTbl)
         return(x) }))
      return(tmpL)
   }
   
   clusterNos = 1:10
   res = NULL; N = 0
   resList = resLabels = list()
   Sys.time()
   for (cl in clusterNos) {
      i.CL <- grep(cl, CLs$group)
      clTbl <- normalizedPBMC[ , i.CL]
      
      stepSize = 1000
      M <- dim(clTbl)[1] - stepSize
      i1 = i2 = 0; res = NULL
      while (i2 < M) {
         i1 <- i2 + 1; i2 <- i2 + stepSize
         iTbl <- clTbl[i1:i2, ]
         res1 <- future_apply(iTbl, 1, FUN = testDE2)
         if (!is.null(res1)) {
            res1 <- do.call(rbind, res1)
            iCheck <- which(res1[ ,"P"] < 0.01 & !is.na(res1[ , "P"]))
            if (length(iCheck) > 0) {
               res1 <- res1[iCheck, ]
               res <- rbind(res, res1)	   
            }
         }
         print(i1)
      }
      
      remaining <- (i2+1):(dim(clTbl)[1])
      iTbl <- clTbl[remaining, ]
      # genes <- rownames(iTbl)
      res2 <- future_apply(iTbl, 1, FUN = testDE2)
      if (!is.null(res2)) {
         res2 <- do.call(rbind, res2)
         iCheck <- which(res2[ ,"P"] < 0.01 & !is.na(res2[ , "P"]))
         if (length(iCheck) > 0) {
            res2 <- res2[res2[ ,"P"] < 0.01 & !is.na(res2[ , "P"]), ]
            res <- rbind(res, res2)
         }
      }
      if (!is.null(res)) {
         N <- N +1
         resList[[N]] <- res
         resLabels[[N]] <- paste0("cluster", cl)
         print(paste0("cluster ", cl, "  time:", Sys.time()))
      }
   }
   
   # res.adj <- lapply(resList, function(x) p.adjust(x, method="fdr")) 
   # NB. Previous method above only corrected for number of tests KEPT!
   if (length(resList) > 0) {
      res.adj <- lapply(resList, adjDF)
      names(res.adj) <- unlist(resLabels)
      res.adj <- lapply(res.adj, function(x) {
         x <- x[x[ , "P"] < 0.05, ]
         return(x)
      })
      res.adj <- compact(res.adj)
      res.adj <- lapply(res.adj, function(x) {
         if (is.matrix(x)) {
            x <- x[order(x[ ,"P"]), ]
         }
      })
   }
   return(res.adj)   
}

T1 = 2
for (T2 in c(4, 6, 8, 18)) {
  print(Sys.time())
  res.adj <- runDE2(time1=T1, time2=T2)
  res.adj <- compact(res.adj)
  setwd("~/__Data/shelfLifeSLE")
  saveRDS(res.adj, file=paste0("shelfLifeSLE_DE", T1, "to", T2, "hr_adjFDR_details.RDS"))
  write.list(z=res.adj, 
             file=paste0("shelfLifeSLE_DE", T1, "to", T2, "hr_pVals_adjFDR_details.CSV"), 
             row.names=TRUE)
}
