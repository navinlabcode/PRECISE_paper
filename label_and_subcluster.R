
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



patient_colors = c("#90CC95","#A882B7","#AA936F","#5167c9","#FC6F60","#ED9D2D","#5E92A4","#a3a3a3","#d8db6e")
sample_colors = c("#90CC95","#2F5932","#A882B7","#3C2544","#AA936F","#685230","#5167c9","#26305e","#FC6F60","#872219","#ED9D2D","#774601","#5E92A4","#112630","#a3a3a3","#4a4a4a","#d8db6e","#5d5e30")
treatment_colors = c("black","#D3BD30")


BCR_all = readRDS("/volumes/lab/users/aschalck/TENX/analysis/BCR/BCR_4/Rda/BCR_fil_seurat.mt.clustered.Rda")


Idents(BCR_all) = "SCT_all_snn_res.0.15"

BCR_all$treatment = c("pre","post")[as.numeric(str_sub(BCR_all$sample, -1))]
BCR_all$treatment = factor(BCR_all$treatment, levels = c("pre","post"))
BCR_all$patient = str_sub(BCR_all$sample, 1, 5)




###############
#####myeloid###
###############
myeloid = SubsetData(BCR_all, ident.use = "4")

VlnPlot(myeloid, assay = "RNA",features = c("IFNA1", "IFNB1","IL6","CXCL10","CD68","PTPRC","CASP1","EPCAM","AIM2","IL1B","IL18","TNF","TMEM173","NFKB1","NFKB2","IRF3"), group.by = "sample", ncol = 4)
VlnPlot(myeloid, assay = "SCT_all",features = c("IFNA1", "IFNB1","IL6","CXCL10","CD68","PTPRC","CASP1","EPCAM","AIM2","IL1B","IL18","TNF","TMEM173","NFKB1","NFKB2","IRF3"), group.by = "sample", ncol = 4)

VlnPlot(myeloid, assay = "SCT_all",features = c("IFNA1", "IFNB1","IL6","CXCL10","CD68","PTPRC","CASP1","EPCAM","AIM2","IL1B","IL18","TNF","TMEM173","NFKB1","NFKB2","IRF3"), group.by = "treatment", ncol = 4)

myeloid = SCTransform(myeloid, new.assay.name = "SCT_myeloid", return.only.var.genes = FALSE, vars.to.regress = "percent.mt")
myeloid_split = SplitObject(myeloid, split.by = "sample")
myeloid_split = lapply(myeloid_split, function(x) return(SCTransform(x, new.assay.name = "SCT_myeloid", return.only.var.genes = FALSE)))
myeloid = merge(myeloid_split[[1]], myeloid_split[2:length(myeloid_split)])

myeloid = FindVariableFeatures(myeloid, assay = "SCT_myeloid")
myeloid = RunPCA(myeloid, npcs = 200, assay = "SCT_myeloid", features = myeloid@assays$SCT_myeloid@var.features)
ElbowPlot(myeloid, ndims = 200)
myeloid = RunUMAP(myeloid, dims = 1:30)
DimPlot(myeloid, group.by = "sample")
myeloid = FindNeighbors(myeloid, dims = 1:30)
myeloid = FindClusters(myeloid, resolution = c(0.6,0.4,0.3,0.2,0.15,0.1,0.05))
DimPlot(myeloid, group.by = "SCT_myeloid_snn_res.0.2")
FeaturePlot(myeloid,slot = "data", features = c("CD68", "PTPRC", "EPCAM","ARG1"))
FeaturePlot(myeloid, features = c("percent.mt","percent.rp","nCount_RNA"))

Idents(myeloid) = "SCT_myeloid_snn_res.0.2"
myeloid_SCT_all_snn_res.0.2_markers = FindAllMarkers(myeloid, only.pos = TRUE)

View(myeloid_SCT_all_snn_res.0.2_markers)

VlnPlot(myeloid, assay = "SCT_myeloid",slot = "data", features = c("IFNA1", "IFNB1","IL6","CXCL10","CD68","PTPRC","CASP1","EPCAM","AIM2","IL1B","IL18","TNF","TMEM173","NFKB1","NFKB2","IRF3","CGAS"), group.by = "sample", ncol = 4)


plot_list = list()
for(gene in c("TMEM173","CGAS","IRF3","IFNA1","IFNB1","NFKB1","NFKB2","TNF","IL6","CASP1","AIM2","IL1B","IL18","STAT6","CXCL10")){
  plot_list[[gene]] = VlnPlot(myeloid, assay = "RNA",slot = "data", features = gene, group.by = "treatment", ncol = 5, pt.size = 0.2, cols = c("grey20", "goldenrod2")) + theme(axis.title = element_blank(), legend.position = "none")
}
plot_grid(plotlist = plot_list, ncol = 5)


plot_list = list()
for(gene in c("TMEM173","CGAS","IRF3","NFKB1","NFKB2","TNF","IL6","CASP1","AIM2","IL1B","IL18","STAT6","CXCL10")){
  plot_list[[gene]] = VlnPlot(myeloid, assay = "SCT_myeloid",slot = "data", features = gene, group.by = "sample", ncol = 5, pt.size = 0.2, cols = sample_colors) + theme(axis.title = element_blank(), legend.position = "none")
}
plot_grid(plotlist = plot_list, ncol = 3)


#VlnPlot(myeloid, assay = "RNA",slot = "data", features = c("IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA11", "IFNA12", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA20", "IFNA21", "IFNA22"), group.by = "treatment", ncol = 5, pt.size = 0.2)
#VlnPlot(BCR_all, assay = "RNA",slot = "data", features = c("IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA11", "IFNA12", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA20", "IFNA21", "IFNA22"), group.by = "SCT_all_snn_res.0.15", ncol = 5, pt.size = 0.2)





Idents(SCTall_fil2_myeloid) = "minor_clusters2"
SCTall_fil2_myeloid_minor_clusters2_RNA_markers = FindAllMarkers(SCTall_fil2_myeloid, assay = "RNA", only.pos = TRUE, logfc.threshold = 2)
SCTall_fil2_myeloid_minor_clusters2_SCT_markers = FindAllMarkers(SCTall_fil2_myeloid, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.1)


SCTall_fil2_myeloid_minor_clusters2_RNA_markers %>% 
  arrange(cluster, -avg_logFC) %>% 
  View()
  

SCTall_fil2_myeloid_minor_clusters2_SCT_markers %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  filter(cluster != "Myeloid - 4") %>% 
  filter(cluster != "Mast Cells") %>% 
  arrange(cluster, -avg_logFC) %>% 
  #filter(avg_logFC >= 0.3) %>% 
  #group_by(cluster) %>% 
  top_n(30, avg_logFC) %>% 
  pull(gene) %>% 
  unique() ->
  test




Idents(SCTall_fil2) = "minor_clusters2"
SCTall_fil2_bcells_cells = WhichCells(SCTall_fil2, idents = c("B - Cells", "Plasma Cells"))
SCTall_fil2_bcells_minor_clusters2_RNA_markers = FindAllMarkers(SCTall_fil2[,SCTall_fil2_bcells_cells], assay = "RNA", only.pos = TRUE, logfc.threshold = 2)
SCTall_fil2_bcells_minor_clusters2_SCT_markers = FindAllMarkers(SCTall_fil2[,SCTall_fil2_bcells_cells], assay = "SCT", only.pos = TRUE, logfc.threshold = 0.1)


SCTall_fil2_bcells_minor_clusters2_SCT_markers %>% 
  arrange(cluster, -avg_logFC) %>% 
  View()



SCTall_fil2_bcells_minor_clusters2_RNA_markers %>% 
  arrange(cluster, -avg_logFC) %>% 
  View()

DimPlot(SCTall_fil2_myeloid)

####

SCTall_fil2_tcells

VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters2", split.by = "timepoint", features = c("CD3G","CD3E","CD8A","CD4"), ncol = 2)

SCTall_fil2_tcells = FindVariableFeatures(SCTall_fil2_tcells, assay = "SCT")
#SCTall_fil2_tcells = SCTransform(SCTall_fil2_tcells, vars.to.regress = "percent.mt")
SCTall_fil2_tcells@assays$SCT@var.features = SCTall_fil2_tcells@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells@assays$SCT@var.features, "^IG[HKL][VDJ]")]
SCTall_fil2_tcells@assays$SCT@var.features = SCTall_fil2_tcells@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells@assays$SCT@var.features, "^TR[ABGD][VDJ]")]
SCTall_fil2_tcells@assays$SCT@var.features = SCTall_fil2_tcells@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells@assays$SCT@var.features, "^HB")]

SCTall_fil2_tcells = RunPCA(SCTall_fil2_tcells, npcs = 200)
ElbowPlot(SCTall_fil2_tcells, ndims = 200)
SCTall_fil2_tcells = RunUMAP(SCTall_fil2_tcells, dims = 1:30, min.dist = 0.5, n.neighbors = 70)

DimPlot(SCTall_fil2_tcells, group.by = "sample", cols = sample_colors)
DimPlot(SCTall_fil2_tcells, group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2_tcells, group.by = "minor_clusters2",cols = NewSet1(11))


SCTall_fil2_tcells$CD4CD8 = "DN"
SCTall_fil2_tcells$CD4CD8[(SCTall_fil2_tcells@assays$RNA@counts["CD8A",] > 0 | SCTall_fil2_tcells@assays$RNA@counts["CD8B",] > 0) & SCTall_fil2_tcells@assays$RNA@counts["CD4",] < 1] = "CD8"
SCTall_fil2_tcells$CD4CD8[(SCTall_fil2_tcells@assays$RNA@counts["CD8A",] < 1 & SCTall_fil2_tcells@assays$RNA@counts["CD8B",] < 1) & SCTall_fil2_tcells@assays$RNA@counts["CD4",] > 0] = "CD4"
SCTall_fil2_tcells$CD4CD8[(SCTall_fil2_tcells@assays$RNA@counts["CD8A",] > 0 | SCTall_fil2_tcells@assays$RNA@counts["CD8B",] > 0) & SCTall_fil2_tcells@assays$RNA@counts["CD4",] > 0] = "DP"


Idents(SCTall_fil2_tcells) = "CD4CD8"
SCTall_fil2_tcells_CD8 = subset(SCTall_fil2_tcells, idents = "CD8")
SCTall_fil2_tcells_CD4 = subset(SCTall_fil2_tcells, idents = "CD4")
SCTall_fil2_tcells_DN = subset(SCTall_fil2_tcells, idents = "DN")
SCTall_fil2_tcells_DP = subset(SCTall_fil2_tcells, idents = "DP")


################
###CD8##########
################

SCTall_fil2_tcells_CD8 = FindVariableFeatures(SCTall_fil2_tcells_CD8, assay = "SCT")
#SCTall_fil2_tcells_CD8 = SCTransform(SCTall_fil2_tcells_CD8, vars.to.regress = "percent.mt")
SCTall_fil2_tcells_CD8@assays$SCT@var.features = SCTall_fil2_tcells_CD8@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_CD8@assays$SCT@var.features, "^IG[HKL][VDJ]")]
SCTall_fil2_tcells_CD8@assays$SCT@var.features = SCTall_fil2_tcells_CD8@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_CD8@assays$SCT@var.features, "^TR[ABGD][VDJ]")]
SCTall_fil2_tcells_CD8@assays$SCT@var.features = SCTall_fil2_tcells_CD8@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_CD8@assays$SCT@var.features, "^HB")]

SCTall_fil2_tcells_CD8 = RunPCA(SCTall_fil2_tcells_CD8, npcs = 200)
ElbowPlot(SCTall_fil2_tcells_CD8, ndims = 200)
SCTall_fil2_tcells_CD8 = RunUMAP(SCTall_fil2_tcells_CD8, dims = 1:30, min.dist = 0.5, n.neighbors = 70)

DimPlot(SCTall_fil2_tcells_CD8, group.by = "sample", cols = sample_colors)
DimPlot(SCTall_fil2_tcells_CD8, group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2_tcells_CD8, group.by = "minor_clusters2",cols = NewSet1(11))

SCTall_fil2_tcells_CD8 = FindNeighbors(SCTall_fil2_tcells_CD8, dims = 1:30, k.param = 30)
SCTall_fil2_tcells_CD8 = FindClusters(SCTall_fil2_tcells_CD8, resolution = c(2,1.5,1,0.8,0.6,0.4))
DimPlot(SCTall_fil2_tcells_CD8, group.by = "SCT_snn_res.2",cols = NewSet1(23))

SCTall_fil2_tcells_CD8 = AddModuleScore(SCTall_fil2_tcells_CD8, name = "naive", features = list(c("CCR7","LEF1","TCF7")), nbin = 15)


SCTall_fil2_tcells_CD8$minor_clusters3 = as.character(SCTall_fil2_tcells_CD8$SCT_snn_res.2)
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters2 %in% c("NK Cells")] = "NK Cells"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters2 %in% c("Tcells - GD")] = "Tcells - GD"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("12","3")] = "Tcells - CD8:Naive"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("3","12","15","6","11","9","13","10") & SCTall_fil2_tcells_CD8$naive1 > 0.25] = "Tcells - CD8:Naive"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("11","13","6")] = "Tcells - CD8:ZNF683"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("15")] = "Tcells - CD8:ISG"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("4", "7", "16")] = "Tcells - CD8:CXCL13"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("14","9")] = "Tcells - CD8:CCL4"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("5","2")] = "Tcells - CD8:GZMK"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("0","1","20","8","18")] = "Tcells - CD8:FCGR3A"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("19")] = "Tcells - Cycling"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("17")] = "Tcells - CD8:MAIT"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("10") & SCTall_fil2_tcells_CD8$minor_clusters2 == "Tcells - CXCL13"] = "Tcells - CD8:CXCL13"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("10") & SCTall_fil2_tcells_CD8$minor_clusters2 == "Tcells - GZMK"] = "Tcells - CD8:GZMK"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("10") & SCTall_fil2_tcells_CD8$minor_clusters2 == "Tcells - Tcm/Naive"] = "Tcells - CD8:Naive"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("10") & SCTall_fil2_tcells_CD8$minor_clusters2 == "Tcells - GZMB"] = "Tcells - CD8:FCGR3A"
SCTall_fil2_tcells_CD8$minor_clusters3[SCTall_fil2_tcells_CD8$minor_clusters3 %in% c("10") & SCTall_fil2_tcells_CD8$minor_clusters2 == "Tcells - MAIT"] = "Tcells - CD8:MAIT"





VlnPlot(SCTall_fil2_tcells_CD8, group.by = "SCT_snn_res.2", features = "CCR7", pt.size = 0)
VlnPlot(SCTall_fil2_tcells_CD8, group.by = "minor_clusters3", features = "TRAV1-2", pt.size = 0)



CD8_clusters = unique(SCTall_fil2_tcells_CD8$minor_clusters3)
CD8_clusters = CD8_clusters[c(6,7,1,4,2,3,11,8,10,9,5)]
SCTall_fil2_tcells_CD8$minor_clusters3 = factor(SCTall_fil2_tcells_CD8$minor_clusters3, levels = CD8_clusters)

DimPlot(SCTall_fil2_tcells_CD8, group.by = "minor_clusters3",cols = NewSet1(11))
DimPlot(SCTall_fil2_tcells_CD8, group.by = "minor_clusters2",cols = NewSet1(10))
DimPlot(SCTall_fil2_tcells_CD8, group.by = "patient",cols = patient_colors)


FeaturePlot(SCTall_fil2_tcells_CD8, features = c("FCGR3A","GNLY"), blend = TRUE, max.cutoff = 1, min.cutoff = 0)

Idents(SCTall_fil2_tcells_CD8) = "minor_clusters3"
SCTall_fil2_tcells_CD8_cells = WhichCells(SCTall_fil2_tcells_CD8, idents = "NK Cells", invert = TRUE)


Idents(SCTall_fil2_tcells_CD8) = "minor_clusters3"
SCTall_fil2_tcells_CD8_cluster3_markers = FindAllMarkers(SCTall_fil2_tcells_CD8, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.1)


SCTall_fil2_tcells_CD8_cluster3_markers %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  filter(avg_logFC > 0.2) %>% 
  filter(pct.2 < 0.85) %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>%
  top_n(5, avg_logFC) %>% 
  ungroup() %>% 
  mutate(cluster = factor(cluster, levels = CD8_clusters)) %>% 
  arrange(cluster) %>% 
  pull(gene) ->
  SCTall_fil2_tcells_CD8_cluster3_markers_top
  
DoHeatmap(SCTall_fil2_tcells_CD8, cells = sample(SCTall_fil2_tcells_CD8$cell.names), features = c("IFNG","TNF","IL2","CCL3","CCL5","ITGB1","PDCD1","LAG3",SCTall_fil2_tcells_CD8_cluster3_markers_top), group.by = "minor_clusters3", group.colors = NewSet1(11)) + scale_fill_gradient2(low="steelblue", mid="white", high="firebrick")  


VlnPlot(SCTall_fil2_tcells_CD8, group.by = "minor_clusters3", features = "GZMK")
VlnPlot(SCTall_fil2_tcells_CD8, group.by = "minor_clusters3", features = "GNLY")
VlnPlot(SCTall_fil2_tcells_CD8, group.by = "minor_clusters3", features = "TCF7")
DimPlot(SCTall_fil2_tcells_CD8, group.by = "SCT_snn_res.1.5",cols = NewSet1(18))
FeaturePlot(SCTall_fil2_tcells_CD8, features = c("CCR7","TCF7"), blend = TRUE, min.cutoff = 0, max.cutoff = 1)



SCTall_fil2_tcells_CD8_cluster3_markers %>% 
  filter(cluster == "Tcells - CD8:ISG") %>% 
  arrange(-avg_logFC) %>% 
  top_n(25, avg_logFC) %>% 
  pull(gene) ->
  ISGs
  



SCTall_fil2_tcells_CD8 = AddModuleScore(SCTall_fil2_tcells_CD8, features = list(ISGs), name = "ISGs", nbin = 15)
hist(SCTall_fil2_tcells_CD8$ISGs1, breaks = 100)


SCTall_fil2_tcells_CD8$stim = "no"
SCTall_fil2_tcells_CD8$stim[SCTall_fil2_tcells_CD8$ISGs1 >= 0.25] = "yes"


DoHeatmap(SCTall_fil2_tcells_CD8, cells=sample(SCTall_fil2_tcells_CD8$cell.names), features = ISGs, group.by = "minor_clusters3", group.colors = NewSet1(13)) + scale_fill_gradient2(low="steelblue", mid="white", high="firebrick")  

FeaturePlot(SCTall_fil2_tcells_CD8, features = "ISGs1", min.cutoff = 0.25, max.cutoff = 0.5)





################
###CD4##########
################

SCTall_fil2_tcells_CD4 = FindVariableFeatures(SCTall_fil2_tcells_CD4, assay = "SCT")
#SCTall_fil2_tcells_CD4 = SCTransform(SCTall_fil2_tcells_CD4, vars.to.regress = "percent.mt")
SCTall_fil2_tcells_CD4@assays$SCT@var.features = SCTall_fil2_tcells_CD4@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_CD4@assays$SCT@var.features, "^IG[HKL][VDJ]")]
SCTall_fil2_tcells_CD4@assays$SCT@var.features = SCTall_fil2_tcells_CD4@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_CD4@assays$SCT@var.features, "^TR[ABGD][VDJ]")]
SCTall_fil2_tcells_CD4@assays$SCT@var.features = SCTall_fil2_tcells_CD4@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_CD4@assays$SCT@var.features, "^HB")]

SCTall_fil2_tcells_CD4 = RunPCA(SCTall_fil2_tcells_CD4, npcs = 200)
ElbowPlot(SCTall_fil2_tcells_CD4, ndims = 200)
SCTall_fil2_tcells_CD4 = RunUMAP(SCTall_fil2_tcells_CD4, dims = 1:30, min.dist = 0.5, n.neighbors = 70)

DimPlot(SCTall_fil2_tcells_CD4, group.by = "sample", cols = sample_colors)
DimPlot(SCTall_fil2_tcells_CD4, group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2_tcells_CD4, group.by = "minor_clusters2",cols = NewSet1(11))

SCTall_fil2_tcells_CD4 = FindNeighbors(SCTall_fil2_tcells_CD4, dims = 1:30, k.param = 30)
SCTall_fil2_tcells_CD4 = FindClusters(SCTall_fil2_tcells_CD4, resolution = c(2,1.5,1,0.8,0.6,0.4))
DimPlot(SCTall_fil2_tcells_CD4, group.by = "SCT_snn_res.0.6",cols = NewSet1(6))
DimPlot(SCTall_fil2_tcells_CD4, group.by = "SCT_snn_res.2",cols = NewSet1(20))


SCTall_fil2_tcells_CD4 = AddModuleScore(SCTall_fil2_tcells_CD4, name = "naive", features = list(c("CCR7","LEF1","TCF7")), nbin = 15)

VlnPlot(SCTall_fil2_tcells_CD4, group.by = "SCT_snn_res.0.6", features = "naive1")

FeaturePlot(SCTall_fil2_tcells_CD4,features = "IL17A")


SCTall_fil2_tcells_CD4$minor_clusters3 = as.character(SCTall_fil2_tcells_CD4$SCT_snn_res.0.6)
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$minor_clusters2 %in% c("Tcells - Cycling")] = "Tcells - Cycling"
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$SCT_snn_res.2 %in% c("13")] = "Tcells - CD4:ISG"
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$SCT_snn_res.2 %in% c("12","2")] = "Tcells - CD4:TREG1"
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$SCT_snn_res.2 %in% c("5")] = "Tcells - CD4:TREG2"
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$minor_clusters3 %in% c("0","1")] = "Tcells - CD4:Naive"
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$minor_clusters3 %in% c("5")] = "Tcells - CD4:GZMK"
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$minor_clusters3 %in% c("2")] = "Tcells - CD4:TREG1"
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$minor_clusters3 %in% c("4")] = "Tcells - CD4:TFH"
SCTall_fil2_tcells_CD4$minor_clusters3[SCTall_fil2_tcells_CD4$minor_clusters3 %in% c("3")] = "Tcells - CD4:TEMRA"

VlnPlot(SCTall_fil2_tcells_CD4, group.by = "SCT_snn_res.2", features = "CCR7", pt.size = 0)
VlnPlot(SCTall_fil2_tcells_CD4, group.by = "minor_clusters3", features = "IL1R1", pt.size = 0)
VlnPlot(SCTall_fil2_tcells_CD4, group.by = "SCT_snn_res.2", features = "CD4", pt.size = 0)


SCTall_fil2_tcells_CD4 = AddModuleScore(SCTall_fil2_tcells_CD4, features = list(ISGs), name="ISGs", nbin = 15)


CD4_clusters = unique(SCTall_fil2_tcells_CD4$minor_clusters3)
CD4_clusters = CD4_clusters[c(1,8,3,7,4,5,6,2)]
SCTall_fil2_tcells_CD4$minor_clusters3 = factor(SCTall_fil2_tcells_CD4$minor_clusters3, levels = CD4_clusters)

DimPlot(SCTall_fil2_tcells_CD4, group.by = "minor_clusters3",cols = NewSet1(8))
DimPlot(SCTall_fil2_tcells_CD4, group.by = "minor_clusters2",cols = NewSet1(10))
DimPlot(SCTall_fil2_tcells_CD4, group.by = "patient",cols = patient_colors)




Idents(SCTall_fil2_tcells_CD4) = "minor_clusters3"
SCTall_fil2_tcells_CD4_cluster3_markers = FindAllMarkers(SCTall_fil2_tcells_CD4, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.1)


SCTall_fil2_tcells_CD4_cluster3_markers %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  #filter(avg_logFC > 0.2) %>% 
  filter(pct.2 < 0.90) %>% 
  filter(pct.1 > 0.5) %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>%
  top_n(10, avg_logFC) %>% 
  ungroup() %>% 
  mutate(cluster = factor(cluster, levels = CD4_clusters)) %>% 
  arrange(cluster) %>% 
  pull(gene) ->
  SCTall_fil2_tcells_CD4_cluster3_markers_top

DoHeatmap(SCTall_fil2_tcells_CD4, cells = sample(SCTall_fil2_tcells_CD4$cell.names), features = c("TCF7","CCR7",SCTall_fil2_tcells_CD4_cluster3_markers_top), group.by = "minor_clusters3", group.colors = NewSet1(11)) + scale_fill_gradient2(low="steelblue", mid="white", high="firebrick")  


SCTall_fil2_tcells_CD4$stim = "no"
SCTall_fil2_tcells_CD4$stim[SCTall_fil2_tcells_CD4$ISGs1 >= 0.25] = "yes"




################
###DN##########
################

SCTall_fil2_tcells_DN = FindVariableFeatures(SCTall_fil2_tcells_DN, assay = "SCT")
#SCTall_fil2_tcells_DN = SCTransform(SCTall_fil2_tcells_DN, vars.to.regress = "percent.mt")
SCTall_fil2_tcells_DN@assays$SCT@var.features = SCTall_fil2_tcells_DN@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_DN@assays$SCT@var.features, "^IG[HKL][VDJ]")]
SCTall_fil2_tcells_DN@assays$SCT@var.features = SCTall_fil2_tcells_DN@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_DN@assays$SCT@var.features, "^TR[ABGD][VDJ]")]
SCTall_fil2_tcells_DN@assays$SCT@var.features = SCTall_fil2_tcells_DN@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_DN@assays$SCT@var.features, "^HB")]

SCTall_fil2_tcells_DN = RunPCA(SCTall_fil2_tcells_DN, npcs = 200)
ElbowPlot(SCTall_fil2_tcells_DN, ndims = 200)
SCTall_fil2_tcells_DN = RunUMAP(SCTall_fil2_tcells_DN, dims = 1:30, min.dist = 0.5, n.neighbors = 70)

DimPlot(SCTall_fil2_tcells_DN, group.by = "sample", cols = sample_colors)
DimPlot(SCTall_fil2_tcells_DN, group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2_tcells_DN, group.by = "minor_clusters2",cols = NewSet1(11))

SCTall_fil2_tcells_DN = FindNeighbors(SCTall_fil2_tcells_DN, dims = 1:30, k.param = 30)
SCTall_fil2_tcells_DN = FindClusters(SCTall_fil2_tcells_DN, resolution = c(2,1.5,1,0.8,0.6,0.4))
DimPlot(SCTall_fil2_tcells_DN, group.by = "SCT_snn_res.0.6",cols = NewSet1(11))
DimPlot(SCTall_fil2_tcells_DN, group.by = "SCT_snn_res.1",cols = NewSet1(14))


SCTall_fil2_tcells_DN = AddModuleScore(SCTall_fil2_tcells_DN, name = "naive", features = list(c("CCR7","LEF1","TCF7")), nbin = 15)

VlnPlot(SCTall_fil2_tcells_DN, group.by = "SCT_snn_res.0.6", features = "naive1", pt.size = 0)

FeaturePlot(SCTall_fil2_tcells_DN,features = "CXCL13")


SCTall_fil2_tcells_DN$minor_clusters3 = paste0("DN: ",as.character(SCTall_fil2_tcells_DN$SCT_snn_res.1))


#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$minor_clusters2 %in% c("Tcells - Cycling")] = "Tcells - Cycling"
#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$SCT_snn_res.2 %in% c("13")] = "Tcells - DN:ISG"
#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$SCT_snn_res.2 %in% c("12","2")] = "Tcells - DN:TREG1"
#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$SCT_snn_res.2 %in% c("5")] = "Tcells - DN:TREG2"
#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$minor_clusters3 %in% c("0","1")] = "Tcells - DN:Naive"
#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$minor_clusters3 %in% c("5")] = "Tcells - DN:GZMK"
#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$minor_clusters3 %in% c("2")] = "Tcells - DN:TREG1"
#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$minor_clusters3 %in% c("4")] = "Tcells - DN:TFH"
#SCTall_fil2_tcells_DN$minor_clusters3[SCTall_fil2_tcells_DN$minor_clusters3 %in% c("3")] = "Tcells - DN:TEMRA"

VlnPlot(SCTall_fil2_tcells_DN, group.by = "SCT_snn_res.1", features = "LEF1", pt.size = 0)
VlnPlot(SCTall_fil2_tcells_DN, group.by = "minor_clusters3", features = "IL1R1", pt.size = 0)
VlnPlot(SCTall_fil2_tcells_DN, group.by = "SCT_snn_res.1", features = "D", pt.size = 0)


SCTall_fil2_tcells_DN = AddModuleScore(SCTall_fil2_tcells_DN, features = list(ISGs), name="ISGs", nbin = 15)


DN_clusters = unique(SCTall_fil2_tcells_DN$minor_clusters3)
DN_clusters = DN_clusters[c(1,8,3,7,4,5,6,2)]
SCTall_fil2_tcells_DN$minor_clusters3 = factor(SCTall_fil2_tcells_DN$minor_clusters3, levels = DN_clusters)

DimPlot(SCTall_fil2_tcells_DN, group.by = "minor_clusters3",cols = NewSet1(8))
DimPlot(SCTall_fil2_tcells_DN, group.by = "minor_clusters2",cols = NewSet1(10))
DimPlot(SCTall_fil2_tcells_DN, group.by = "patient",cols = patient_colors)




Idents(SCTall_fil2_tcells_DN) = "minor_clusters3"
SCTall_fil2_tcells_DN_cluster3_markers = FindAllMarkers(SCTall_fil2_tcells_DN, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.1)


SCTall_fil2_tcells_DN_cluster3_markers %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  #filter(avg_logFC > 0.2) %>% 
  #filter(pct.2 < 0.90) %>% 
  #filter(pct.1 > 0.5) %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>%
  top_n(10, avg_logFC) %>% 
  ungroup() %>% 
  #mutate(cluster = factor(cluster, levels = DN_clusters)) %>% 
  arrange(cluster) %>% 
  pull(gene) ->
  SCTall_fil2_tcells_DN_cluster3_markers_top

DoHeatmap(SCTall_fil2_tcells_DN, cells = sample(SCTall_fil2_tcells_DN$cell.names), features = c("TCF7","CCR7",SCTall_fil2_tcells_DN_cluster3_markers_top), group.by = "minor_clusters3", group.colors = NewSet1(11)) + scale_fill_gradient2(low="steelblue", mid="white", high="firebrick")  


SCTall_fil2_tcells_DN$stim = "no"
SCTall_fil2_tcells_DN$stim[SCTall_fil2_tcells_DN$ISGs1 >= 0.25] = "yes"





################
###DP##########
################

SCTall_fil2_tcells_DP = FindVariableFeatures(SCTall_fil2_tcells_DP, assay = "SCT")
#SCTall_fil2_tcells_DP = SCTransform(SCTall_fil2_tcells_DP, vars.to.regress = "percent.mt")
SCTall_fil2_tcells_DP@assays$SCT@var.features = SCTall_fil2_tcells_DP@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_DP@assays$SCT@var.features, "^IG[HKL][VDJ]")]
SCTall_fil2_tcells_DP@assays$SCT@var.features = SCTall_fil2_tcells_DP@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_DP@assays$SCT@var.features, "^TR[ABGD][VDJ]")]
SCTall_fil2_tcells_DP@assays$SCT@var.features = SCTall_fil2_tcells_DP@assays$SCT@var.features[!str_detect(SCTall_fil2_tcells_DP@assays$SCT@var.features, "^HB")]

SCTall_fil2_tcells_DP = RunPCA(SCTall_fil2_tcells_DP, npcs = 200)
ElbowPlot(SCTall_fil2_tcells_DP, ndims = 200)
SCTall_fil2_tcells_DP = RunUMAP(SCTall_fil2_tcells_DP, dims = 1:30, min.dist = 0.5, n.neighbors = 20)

DimPlot(SCTall_fil2_tcells_DP, group.by = "sample", cols = sample_colors)
DimPlot(SCTall_fil2_tcells_DP, group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2_tcells_DP, group.by = "minor_clusters2",cols = NewSet1(11))

SCTall_fil2_tcells_DP = FinDPeighbors(SCTall_fil2_tcells_DP, dims = 1:30, k.param = 30)
SCTall_fil2_tcells_DP = FindClusters(SCTall_fil2_tcells_DP, resolution = c(2,1.5,1,0.8,0.6,0.4))
DimPlot(SCTall_fil2_tcells_DP, group.by = "SCT_snn_res.0.6",cols = NewSet1(11))
DimPlot(SCTall_fil2_tcells_DP, group.by = "SCT_snn_res.1",cols = NewSet1(14))


SCTall_fil2_tcells_DP = AddModuleScore(SCTall_fil2_tcells_DP, name = "naive", features = list(c("CCR7","LEF1","TCF7")), nbin = 15)

VlnPlot(SCTall_fil2_tcells_DP, group.by = "SCT_snn_res.0.6", features = "naive1", pt.size = 0)

FeaturePlot(SCTall_fil2_tcells_DP,features = "CXCL13")


SCTall_fil2_tcells_DP$minor_clusters3 = paste0("DP: ", as.character(SCTall_fil2_tcells_DP$minor_clusters2))


#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$minor_clusters2 %in% c("Tcells - Cycling")] = "Tcells - Cycling"
#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$SCT_snn_res.2 %in% c("13")] = "Tcells - DP:ISG"
#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$SCT_snn_res.2 %in% c("12","2")] = "Tcells - DP:TREG1"
#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$SCT_snn_res.2 %in% c("5")] = "Tcells - DP:TREG2"
#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$minor_clusters3 %in% c("0","1")] = "Tcells - DP:Naive"
#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$minor_clusters3 %in% c("5")] = "Tcells - DP:GZMK"
#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$minor_clusters3 %in% c("2")] = "Tcells - DP:TREG1"
#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$minor_clusters3 %in% c("4")] = "Tcells - DP:TFH"
#SCTall_fil2_tcells_DP$minor_clusters3[SCTall_fil2_tcells_DP$minor_clusters3 %in% c("3")] = "Tcells - DP:TEMRA"

VlnPlot(SCTall_fil2_tcells_DP, group.by = "SCT_snn_res.1", features = "LEF1", pt.size = 0)
VlnPlot(SCTall_fil2_tcells_DP, group.by = "minor_clusters3", features = "IL1R1", pt.size = 0)
VlnPlot(SCTall_fil2_tcells_DP, group.by = "SCT_snn_res.1", features = "D", pt.size = 0)


SCTall_fil2_tcells_DP = AddModuleScore(SCTall_fil2_tcells_DP, features = list(ISGs), name="ISGs", nbin = 15)


DP_clusters = unique(SCTall_fil2_tcells_DP$minor_clusters3)
DP_clusters = DP_clusters[c(1,8,3,7,4,5,6,2)]
SCTall_fil2_tcells_DP$minor_clusters3 = factor(SCTall_fil2_tcells_DP$minor_clusters3, levels = DP_clusters)

DimPlot(SCTall_fil2_tcells_DP, group.by = "minor_clusters3",cols = NewSet1(8))
DimPlot(SCTall_fil2_tcells_DP, group.by = "minor_clusters2",cols = NewSet1(10))
DimPlot(SCTall_fil2_tcells_DP, group.by = "patient",cols = patient_colors)




Idents(SCTall_fil2_tcells_DP) = "minor_clusters3"
SCTall_fil2_tcells_DP_cluster3_markers = FindAllMarkers(SCTall_fil2_tcells_DP, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.1)


SCTall_fil2_tcells_DP_cluster3_markers %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  #filter(avg_logFC > 0.2) %>% 
  #filter(pct.2 < 0.90) %>% 
  #filter(pct.1 > 0.5) %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>%
  top_n(10, avg_logFC) %>% 
  ungroup() %>% 
  #mutate(cluster = factor(cluster, levels = DP_clusters)) %>% 
  arrange(cluster) %>% 
  pull(gene) ->
  SCTall_fil2_tcells_DP_cluster3_markers_top

DoHeatmap(SCTall_fil2_tcells_DP, cells = sample(SCTall_fil2_tcells_DP$cell.names), features = c("TCF7","CCR7",SCTall_fil2_tcells_DP_cluster3_markers_top), group.by = "minor_clusters3", group.colors = NewSet1(11)) + scale_fill_gradient2(low="steelblue", mid="white", high="firebrick")  


SCTall_fil2_tcells_DP$stim = "no"
SCTall_fil2_tcells_DP$stim[SCTall_fil2_tcells_DP$ISGs1 >= 0.25] = "yes"


############################
####combine#################
############################



SCTall_fil2_tcells$minor_clusters3 = NA
SCTall_fil2_tcells$minor_clusters3[SCTall_fil2_tcells_CD8$cell.names] = as.character(SCTall_fil2_tcells_CD8$minor_clusters3)
SCTall_fil2_tcells$minor_clusters3[SCTall_fil2_tcells_CD4$cell.names] = as.character(SCTall_fil2_tcells_CD4$minor_clusters3)
SCTall_fil2_tcells$minor_clusters3[SCTall_fil2_tcells_DN$cell.names] = SCTall_fil2_tcells_DN$minor_clusters3
SCTall_fil2_tcells$minor_clusters3[SCTall_fil2_tcells_DP$cell.names] = SCTall_fil2_tcells_DP$minor_clusters3


Idents(SCTall_fil2_tcells) = "minor_clusters3"
temp = subset(SCTall_fil2_tcells, idents=c("DN: 9","DP: Tcells - Treg"))
temp = FindVariableFeatures(temp, assay = "SCT")
#temp = SCTransform(temp, vars.to.regress = "percent.mt")
temp@assays$SCT@var.features = temp@assays$SCT@var.features[!str_detect(temp@assays$SCT@var.features, "^IG[HKL][VDJ]")]
temp@assays$SCT@var.features = temp@assays$SCT@var.features[!str_detect(temp@assays$SCT@var.features, "^TR[ABGD][VDJ]")]
temp@assays$SCT@var.features = temp@assays$SCT@var.features[!str_detect(temp@assays$SCT@var.features, "^HB")]
temp = RunPCA(temp, npcs = 200)
ElbowPlot(temp, ndims = 200)
temp = RunUMAP(temp, dims = 1:15, min.dist = 0.5, n.neighbors = 20)
DimPlot(temp, group.by = "sample", cols = sample_colors)
DimPlot(temp, group.by = "patient", cols = patient_colors)
DimPlot(temp, group.by = "minor_clusters2",cols = NewSet1(11))
temp = FindNeighbors(temp, dims = 1:15)
temp = FindClusters(temp, resolution = c(2,1.5,1,0.8,0.6,0.4))
DimPlot(temp, group.by = "SCT_snn_res.0.6",cols = NewSet1(5))
DimPlot(temp, group.by = "SCT_snn_res.1",cols = NewSet1(14))
Idents(temp) = "SCT_snn_res.0.6"
temp_markers = FindAllMarkers(temp, only.pos = TRUE, logfc.threshold = 0.1)
temp$minor_clusters3 = as.character(temp$SCT_snn_res.0.6) 
temp$minor_clusters3[temp$minor_clusters3 %in% c("0","1","3")] = "0"
Idents(temp) = "minor_clusters3"
temp_markers = FindAllMarkers(temp, only.pos = TRUE, logfc.threshold = 0.1)

FeaturePlot(temp, features = "LAG3")





SCTall_fil2_tcells$minor_clusters4 = as.character(SCTall_fil2_tcells$minor_clusters3)
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 == "DP: Tcells - Cycling"] = "Tcells - Cycling"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 1","DN: 12","DN: 3","DP: Tcells - Tcm/Naive","DN: 2","DN: 6")] = "Tcells - CD4:Naive"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 5")] = "Tcells - CD8:Naive"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DP: Tcells - CXCL13")] = "Tcells - CD8:CXCL13"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 13")] = "Tcells - CD4:TFH"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 0")] = "NK Cells"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 7") & SCTall_fil2_tcells$minor_clusters2 %in% c("NK Cells")] = "NK Cells"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 10","DP: Tcells - GD")] = "Tcells - GD"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 8","DP: Tcells - GZMB")] = "Tcells - CD4:TEMRA"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 9","DP: Tcells - Treg")] = "Tcells - CD4:TREG1"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 11","DP: Tcells - MAIT")] = "Tcells - CD8:MAIT"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 4") & SCTall_fil2_tcells@assays$RNA@counts["FCGR3A",] > 0] = "Tcells - CD8:FCGR3A"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 4")] = "Tcells - CD4:TEMRA"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 7") & SCTall_fil2_tcells$minor_clusters2 %in% c("Tcells - GZMB")] = "NK Cells"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DN: 7")] = "Tcells - CD8:CCL4"
SCTall_fil2_tcells$minor_clusters4[SCTall_fil2_tcells$minor_clusters3 %in% c("DP: Tcells - GZMK")] = "Tcells - CD8:GZMK"


DimPlot(SCTall_fil2_tcells, group.by = "minor_clusters4",cols = NewSet1(18), pt.size = 0.5)

Idents(SCTall_fil2_tcells) = "minor_clusters4"
SCTall_fil2_tcells_minor_clusters4_markers = FindAllMarkers(SCTall_fil2_tcells, logfc.threshold = 0.1, assay = "SCT", only.pos = TRUE)
SCTall_fil2_tcells_minor_clusters4_markers %>% 
  filter(!str_detect(gene, "^R[LS]")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  group_by(cluster) %>% 
  top_n(20, avg_logFC) %>% 
  pull(gene) ->
  SCTall_fil2_tcells_minor_clusters4_markers_top

DoHeatmap(SCTall_fil2_tcells, cells = sample(SCTall_fil2_tcells$cell.names), features = c("CD8A","CD4",SCTall_fil2_tcells_minor_clusters4_markers_top), group.by = "minor_clusters4") + 
  scale_fill_gradient2(low="steelblue", mid="white", high="firebrick")  


SCTall_fil2_tcells = AddModuleScore(SCTall_fil2_tcells, features = list(ISGs), name = "ISGs", nbin = 20)
hist(SCTall_fil2_tcells$ISGs1, breaks = 100)
SCTall_fil2_tcells$stim = "no"
SCTall_fil2_tcells$stim[SCTall_fil2_tcells$ISGs1 > 0.25] = "yes"

table(SCTall_fil2_tcells$minor_clusters4,SCTall_fil2_tcells$sample)

SCTall_fil2_tcells$CD4CD8_cluster = "other"
SCTall_fil2_tcells$CD4CD8_cluster[str_detect(SCTall_fil2_tcells$minor_clusters4, "CD4")] = "CD4"
SCTall_fil2_tcells$CD4CD8_cluster[str_detect(SCTall_fil2_tcells$minor_clusters4, "CD8")] = "CD8"





##################
####neutrophils###
##################

Idents(SCTall_fil2) = "minor_clusters2"
SCTall_fil2_neu = subset(SCTall_fil2, idents =c("Myeloid - 3","Myeloid - 4","Myeloid - 5"))

SCTall_fil2_neu = FindVariableFeatures(SCTall_fil2_neu, assay = "SCT")
#SCTall_fil2_neu = SCTransform(SCTall_fil2_neu, vars.to.regress = "percent.mt")
SCTall_fil2_neu@assays$SCT@var.features = SCTall_fil2_neu@assays$SCT@var.features[!str_detect(SCTall_fil2_neu@assays$SCT@var.features, "^IG[HKL][VDJ]")]
SCTall_fil2_neu@assays$SCT@var.features = SCTall_fil2_neu@assays$SCT@var.features[!str_detect(SCTall_fil2_neu@assays$SCT@var.features, "^TR[ABGD][VDJ]")]
SCTall_fil2_neu@assays$SCT@var.features = SCTall_fil2_neu@assays$SCT@var.features[!str_detect(SCTall_fil2_neu@assays$SCT@var.features, "^HB")]

SCTall_fil2_neu = RunPCA(SCTall_fil2_neu, npcs = 200)
ElbowPlot(SCTall_fil2_neu, ndims = 200)
SCTall_fil2_neu = RunUMAP(SCTall_fil2_neu, dims = 1:30, min.dist = 0.5, n.neighbors = 20)

DimPlot(SCTall_fil2_neu, group.by = "sample", cols = sample_colors)
DimPlot(SCTall_fil2_neu, group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2_neu, group.by = "minor_clusters2",cols = NewSet1(3))

SCTall_fil2_neu = FindNeighbors(SCTall_fil2_neu, dims = 1:30, k.param = 30)
SCTall_fil2_neu = FindClusters(SCTall_fil2_neu, resolution = c(2,1.5,1,0.8,0.6,0.5,0.4))
DimPlot(SCTall_fil2_neu, group.by = "SCT_snn_res.0.6",cols = NewSet1(9))
DimPlot(SCTall_fil2_neu, group.by = "SCT_snn_res.1",cols = NewSet1(14))


SCTall_fil2_neu$minor_clusters3 = as.character(SCTall_fil2_neu$SCT_snn_res.0.6)
SCTall_fil2_neu$minor_clusters3[SCTall_fil2_neu$minor_clusters3 %in% c("0","4","5","2","1")] = "Neutrophils - VCAN"
SCTall_fil2_neu$minor_clusters3[SCTall_fil2_neu$minor_clusters3 %in% c("3","6")] = "Neutrophils - FCGR3A"


DimPlot(SCTall_fil2_neu, group.by = "minor_clusters3",cols = NewSet1(6))

VlnPlot(SCTall_fil2_neu, group.by = "minor_clusters3", features = "CSF3R")
FeaturePlot(SCTall_fil2_neu, features = c("VMO1"),  min.cutoff = 0, max.cutoff = 0.1)
FeaturePlot(SCTall_fil2_neu, features = c("ICAM1"),  min.cutoff = 0, max.cutoff = 0.1)


Idents(SCTall_fil2_neu) = "minor_clusters3"
SCTall_fil2_neu_minor_clusters3_markers = FindAllMarkers(SCTall_fil2_neu, only.pos = TRUE, assay = "SCT", logfc.threshold = 0.1)
SCTall_fil2_neu_minor_clusters3_markers %>% 
  group_by(cluster) %>%
  filter(pct.2 < 0.85) %>% 
  arrange(cluster, -avg_logFC) %>% 
  top_n(20, avg_logFC) %>% 
  View()


Idents(SCTall_fil2_myeloid) = "minor_clusters2"
SCTall_fil2_myeloid_minor_clusters2_markers = FindAllMarkers(SCTall_fil2_myeloid, only.pos = TRUE, assay = "SCT", logfc.threshold = 0.1)
SCTall_fil2_myeloid_minor_clusters2_markers %>% 
  arrange(cluster, -avg_logFC) %>% 
  View()

FeaturePlot(SCTall_fil2_myeloid, features = "ARG1")
DimPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", cols = NewSet1(11))


