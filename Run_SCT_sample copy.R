#/usr/bin/env Rscript

require(Seurat)
library(tidyverse)
library(stringr)
require(future)
require(future.apply)


plan("multiprocess", workers = 30)
plan()
options(future.globals.maxSize = 100000 * 1024^2)


BCR_fil_seurat = readRDS(file = "/volumes/lab/users/aschalck/TENX/analysis/BCR/BCR_4/Rda/fil_seurats_combined.Rda")

BCR_fil_seurat = SCTransform(BCR_fil_seurat, new.assay.name = "SCT_all", variable.features.n = 5000, vars.to.regress = c("sample"))
BCR_fil_seurat = FindVariableFeatures(BCR_fil_seurat, nfeatures = 5000, assay = "SCT_all")
BCR_fil_seurat@assays$SCT_all@var.features = BCR_fil_seurat@assays$SCT_all@var.features[!str_detect(BCR_fil_seurat@assays$SCT_all@var.features, "^TR[A|B|D|G][V|D|J]")]


BCR_fil_seurat = RunPCA(BCR_fil_seurat, npcs = 200,assay = "SCT_all")


BCR_fil_seurat = RunUMAP(BCR_fil_seurat, dims = 1:30, reduction = "pca", seed.use = 77, reduction.name = "UMAP_30", reduction.key = "UMAP_30_")
BCR_fil_seurat = RunUMAP(BCR_fil_seurat, dims = 1:50, reduction = "pca", seed.use = 77, reduction.name = "UMAP_50", reduction.key = "UMAP_50_")

saveRDS(BCR_fil_seurat, "/volumes/lab/users/aschalck/TENX/analysis/BCR/BCR_4/Rda/BCR_fil_seurat.sample.Rda", compress = FALSE)
