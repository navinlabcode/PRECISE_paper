
library(DropletUtils)
library(Seurat)
library(tidyverse)
library(wesanderson)
library(reshape2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(cowplot)
library(devtools)
library(beyonce, lib.loc = "/volumes/lab/users/aschalck/R/")
library(ggalluvial)
require(foreign)
require(nnet)
library(lme4)
library(GGally)
library(fgsea)
library(future)
library(future.apply)
library(tidyverse)
library(fuzzyjoin)
source("~/TENX/analysis/MP/MP_project_4/TCR_functions_script.R")
options(future.globals.maxSize = 60000 * 1024^2)



####################################################
####read in the emptydorp filtered matrix files#####
####################################################

parent_dir = "/volumes/USR1/aschalck/TENX/analysis/BCR/BCR_4/empty_drops/"

sample_dirs = list.files(parent_dir)
sample_dirs = sample_dirs[str_detect(sample_dirs,"^BCR")]
sample_dirs = paste0(parent_dir, sample_dirs)
sample_names = str_sub(basename(sample_dirs), 1, 7)


raw_mats = lapply(sample_dirs, function(x){return(Read10X(data.dir = x))})
names(raw_mats) = sample_names

raw_seurats = lapply(raw_mats, function(x){return(CreateSeuratObject(x,project = "BCR"))})


raw_seurats = mapply(function(x,y) {x@meta.data$sample = y;return(x)}, raw_seurats, sample_names)

raw_seurats = lapply(raw_seurats, function(x){x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-");return(x)})
ribo.genes = grep(x=row.names(raw_seurats[[1]]@assays$RNA$counts), pattern = "^(RP[LS])|(MRP[LS])",value = TRUE)
raw_seurats = lapply(raw_seurats, function(x){x[["percent.rp"]] <- PercentageFeatureSet(x, features = ribo.genes);return(x)})
raw_seurats = lapply(raw_seurats, function(x){x@meta.data$umi.gene.ratio = x@meta.data$nCount_RNA/x@meta.data$nFeature_RNA;return(x)})

raw_seurats_combined = merge(x = raw_seurats[[1]], y = raw_seurats[2:length(raw_seurats)], add.cell.ids = names(raw_seurats))

raw_seurats_combined$patient = str_sub(raw_seurats_combined$sample, 1, 5)

raw_seurats_combined$patient.2 = patients_key[raw_seurats_combined$patient] %>% setNames(NULL)
raw_seurats_combined$patient.2 = factor(raw_seurats_combined$patient.2, levels = unique(raw_seurats_combined$patient.2) %>% sort())

raw_seurats_combined$timepoint = c("pre","post")[as.numeric(str_sub(raw_seurats_combined$sample, 7, 7))]
raw_seurats_combined$timepoint = factor(raw_seurats_combined$timepoint, levels =  c("pre","post"))

VlnPlot(raw_seurats_combined, features = "percent.mt", 
        cols = treatment_colors, split.by = "timepoint", 
        group.by = "patient.2", pt.size = 0, split.plot = T) +
  geom_hline(yintercept = 20) + theme(axis.title.x = element_blank()) -> p
  ggsave(p,filename = "/volumes/USR1/aschalck/projects/Breast_cancer_radiation/revision/percent.mt.raw.png", height = 2.5, width = 5, units = "in", device = "png")

#dir.create(paste0(parent_dir,"vln_plots"))

for(i in 1:length(raw_seurats)){
  temp_plot = VlnPlot(object = raw_seurats[[i]], features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp","umi.gene.ratio"), ncol=5)
  ggsave(temp_plot, filename = paste0(parent_dir,"vln_plots/",sample_names[i],"_prefilter.png"), width = 16, height = 6, units = "in", device = "png") 
}

fil_seurats = lapply(raw_seurats, function(x){return(subset(x, subset = percent.mt < 20 & nCount_RNA < 40000 & nCount_RNA > 200 & umi.gene.ratio < 10 & umi.gene.ratio > 1.3 & percent.rp < 40 & nFeature_RNA <7500))})

for(i in 1:length(fil_seurats)){
  temp_plot = VlnPlot(object = fil_seurats[[i]], features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rp","umi.gene.ratio"), ncol=5)
  ggsave(temp_plot, filename = paste0(parent_dir,"vln_plots/",sample_names[i],"_postfilter.png"), width = 16, height = 6, units = "in", device = "png") 
}


saveRDS(fil_seurats, file = "/volumes/lab/users/aschalck/TENX/analysis/BCR/BCR_4/Rda/fil_seurats.Rda", compress = FALSE)


fil_seurats_combined = merge(x = fil_seurats[[1]], y = fil_seurats[2:length(fil_seurats)], add.cell.ids = names(fil_seurats))

saveRDS(fil_seurats_combined, file = "/volumes/lab/users/aschalck/TENX/analysis/BCR/BCR_4/Rda/fil_seurats_combined.Rda", compress = FALSE)



