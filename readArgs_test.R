#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
   stop("need to specify names of filtered_feature_bc_matrix, CITE-seq-count_res, & outputFile", call.=FALSE)
} 

df = read.csv(args[1], header=TRUE)
num_vars = which(sapply(df, class)=="numeric")
df_out = df[ ,num_vars]
write.table(df_out, file=args[2], row.names=FALSE)

# Rscript --vanilla readArgs_test.R /home/hamid.bolouri/__Data/signatures/angelova.2015.markers.csv /home/hamid.bolouri/__Data/readArgsTest_out.txt