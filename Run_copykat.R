#/usr/bin/env Rscript

library(copykat)
library(Seurat)

setwd("/volumes/lab/users/aschalck/TENX/analysis/BCR/BCR_4/copykat")
#patients = c("BCR01","BCR03","BCR04","BCR05","BCR06","BCR09","BCR13","BCR15","BCR17","BCR18","BCR20")
patients = c("BCR05","BCR13")


for(i in patients){
  print(i)
  temp = Read10X(paste0("/volumes/lab/users/aschalck/TENX/analysis/BCR/BCR_4/data/fil_mats/",i,"/"))
  cktemp = copykat(temp, sam.name = i, n.cores = 60)
  saveRDS(cktemp, paste0("./Rda/",i,".rds"),compress=FALSE)
}

print("done")





