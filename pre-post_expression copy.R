


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
NewSet1 = function(n){return(colorRampPalette(brewer.pal(9,"Set1")[c(1,5,6,3,2,4,8,7,9)])(n))}


#############################
#######myeloid pre-post######
#############################

Idents(SCTall_fil2) = "minor_clusters2"
myeloid_clusters = unique(SCTall_fil2$minor_clusters2)[str_detect(unique(SCTall_fil2$minor_clusters2), "Myeloid")]
SCTall_fil_myeloid = subset(SCTall_fil2, idents="Myeloid")



Idents(SCTall_fil_myeloid) = "timepoint"
SCTall_fil_myeloid_timepoint_SCT_markers = FindAllMarkers(SCTall_fil_myeloid, logfc.threshold = 0.2, assay = "SCT", only.pos = TRUE)
SCTall_fil_myeloid_timepoint_RNA_markers = FindAllMarkers(SCTall_fil_myeloid, logfc.threshold = 2, assay = "RNA", only.pos = TRUE)


myeloid_prepost_SCT_markers = data.frame()
for(i in patients){
  Idents(SCTall_fil_myeloid) = "patient"
  cells = WhichCells(SCTall_fil_myeloid, idents = i)
  Idents(SCTall_fil_myeloid) = "timepoint"
  marks = FindAllMarkers(SCTall_fil_myeloid[,cells], logfc.threshold = 0.2, assay = "SCT", only.pos = TRUE)
  marks$patient = i
  myeloid_prepost_SCT_markers = rbind(myeloid_prepost_SCT_markers, marks)
}


myeloid_prepost_SCT_markers %>% 
  filter(cluster == "post") %>% 
  arrange(patient, -avg_logFC) %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  filter(n>=8) %>%  
  pull(gene) ->
  myeloid_post_markers_top


myeloid_prepost_SCT_markers %>% 
  filter(cluster == "pre") %>% 
  arrange(patient, -avg_logFC) %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  filter(n>=8) %>%  
  pull(gene) ->
  myeloid_pre_markers_top

myeloid_heatmaps = list()
for(i in patients){
  Idents(SCTall_fil_myeloid) = "patient"
  cells = WhichCells(SCTall_fil_myeloid, idents = i)
  plot = DoHeatmap(SCTall_fil_myeloid, cells = sample(cells),
                   group.by = "timepoint", assay = "SCT", 
                   features = c(myeloid_post_markers_top, myeloid_pre_markers_top)) + 
                   scale_fill_gradient2(low="steelblue",mid="white",high = "firebrick") +
                   theme(legend.position = "none")
  myeloid_heatmaps[[i]] = plot
}

plot_grid(plotlist = myeloid_heatmaps, nrow = 3, labels = names(myeloid_heatmaps))





#############################
#######fibro pre-post######
#############################

Idents(SCTall_fil) = "major_clusters"
SCTall_fil_fibro = subset(SCTall_fil, idents="Fibroblasts")

Idents(SCTall_fil_fibro) = "timepoint"
SCTall_fil_fibro_timepoint_SCT_markers = FindAllMarkers(SCTall_fil_fibro, logfc.threshold = 0.2, assay = "SCT", only.pos = TRUE)
SCTall_fil_fibro_timepoint_RNA_markers = FindAllMarkers(SCTall_fil_fibro, logfc.threshold = 2, assay = "RNA", only.pos = TRUE)


fibro_prepost_SCT_markers = data.frame()
for(i in patients[c(2,7,8,11)]){
  Idents(SCTall_fil_fibro) = "patient"
  cells = WhichCells(SCTall_fil_fibro, idents = i)
  Idents(SCTall_fil_fibro) = "timepoint"
  marks = FindAllMarkers(SCTall_fil_fibro[,cells], logfc.threshold = 0.2, assay = "SCT", only.pos = TRUE)
  marks$patient = i
  fibro_prepost_SCT_markers = rbind(fibro_prepost_SCT_markers, marks)
}


fibro_prepost_SCT_markers %>% 
  filter(cluster == "post") %>% 
  arrange(patient, -avg_logFC) %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  filter(n>=3) %>%  
  pull(gene) ->
  fibro_post_markers_top


fibro_prepost_SCT_markers %>% 
  filter(cluster == "pre") %>% 
  arrange(patient, -avg_logFC) %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  filter(n>=3) %>%  
  pull(gene) ->
  fibro_pre_markers_top

fibro_heatmaps = list()
for(i in patients[c(2,7,8,11)]){
  Idents(SCTall_fil_fibro) = "patient"
  cells = WhichCells(SCTall_fil_fibro, idents = i)
  plot = DoHeatmap(SCTall_fil_fibro, cells = sample(cells),
                   group.by = "timepoint", assay = "SCT", 
                   features = c(fibro_post_markers_top, fibro_pre_markers_top)) + 
    scale_fill_gradient2(low="steelblue",mid="white",high = "firebrick") +
    theme(legend.position = "none")
  fibro_heatmaps[[i]] = plot
}

plot_grid(plotlist = fibro_heatmaps, nrow = 2, labels = names(fibro_heatmaps))





#############################
#######tcells pre-post#######
#############################

#Idents(SCTall_fil) = "major_clusters"
#SCTall_fil_tcells = subset(SCTall_fil, idents="T - Cells")

Idents(SCTall_tcells_fil) = "timepoint"
SCTall_fil_tcells_timepoint_SCT_markers = FindAllMarkers(SCTall_tcells_fil, logfc.threshold = 0.2, assay = "SCT", only.pos = TRUE)
SCTall_fil_tcells_timepoint_RNA_markers = FindAllMarkers(SCTall_tcells_fil, logfc.threshold = 2, assay = "RNA", only.pos = TRUE)


tcells_prepost_SCT_markers = data.frame()
for(i in patients){
  Idents(SCTall_tcells_fil) = "patient"
  cells = WhichCells(SCTall_tcells_fil, idents = i)
  Idents(SCTall_tcells_fil) = "timepoint"
  marks = FindAllMarkers(SCTall_tcells_fil[,cells], logfc.threshold = 0.2, assay = "SCT", only.pos = TRUE)
  marks$patient = i
  tcells_prepost_SCT_markers = rbind(tcells_prepost_SCT_markers, marks)
}


tcells_prepost_SCT_markers %>% 
  filter(cluster == "post") %>% 
  arrange(patient, -avg_logFC) %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  filter(n>=8) %>%  
  pull(gene) ->
  tcells_post_markers_top


tcells_prepost_SCT_markers %>% 
  filter(cluster == "pre") %>% 
  arrange(patient, -avg_logFC) %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  filter(n>=8) %>%  
  pull(gene) ->
  tcells_pre_markers_top

tcells_heatmaps = list()
for(i in patients){
  Idents(SCTall_tcells_fil) = "patient"
  cells = WhichCells(SCTall_tcells_fil, idents = i)
  plot = DoHeatmap(SCTall_tcells_fil, cells = sample(cells),
                   group.by = "timepoint", assay = "SCT", 
                   features = c(tcells_post_markers_top, tcells_pre_markers_top)) + 
    scale_fill_gradient2(low="steelblue",mid="white",high = "firebrick") +
    theme(legend.position = "none")
  tcells_heatmaps[[i]] = plot
}

plot_grid(plotlist = tcells_heatmaps, nrow = 4, labels = names(tcells_heatmaps))





#############################
#######bcells pre-post#######
#############################

Idents(SCTall_fil) = "major_clusters"
SCTall_fil_bcells = subset(SCTall_fil, idents="B - Cells")

Idents(SCTall_fil_bcells) = "timepoint"
SCTall_fil_bcells_timepoint_SCT_markers = FindAllMarkers(SCTall_fil_bcells, logfc.threshold = 0.2, assay = "SCT", only.pos = TRUE)
SCTall_fil_bcells_timepoint_RNA_markers = FindAllMarkers(SCTall_fil_bcells, logfc.threshold = 2, assay = "RNA", only.pos = TRUE)


bcells_prepost_SCT_markers = data.frame()
for(i in patients[c(2,4,5,6,7,9,10,11)]){
  Idents(SCTall_fil_bcells) = "patient"
  cells = WhichCells(SCTall_fil_bcells, idents = i)
  Idents(SCTall_fil_bcells) = "timepoint"
  marks = FindAllMarkers(SCTall_fil_bcells[,cells], logfc.threshold = 0.2, assay = "SCT", only.pos = TRUE)
  marks$patient = i
  bcells_prepost_SCT_markers = rbind(bcells_prepost_SCT_markers, marks)
}


bcells_prepost_SCT_markers %>% 
  filter(cluster == "post") %>% 
  arrange(patient, -avg_logFC) %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  filter(n>=5) %>%  
  pull(gene) ->
  bcells_post_markers_top


bcells_prepost_SCT_markers %>% 
  filter(cluster == "pre") %>% 
  arrange(patient, -avg_logFC) %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  filter(n>=5) %>%  
  pull(gene) ->
  bcells_pre_markers_top

bcells_heatmaps = list()
for(i in patients[c(2,4,5,6,7,9,10,11)]){
  Idents(SCTall_fil_bcells) = "patient"
  cells = WhichCells(SCTall_fil_bcells, idents = i)
  plot = DoHeatmap(SCTall_fil_bcells, cells = sample(cells),
                   group.by = "timepoint", assay = "SCT", 
                   features = c(bcells_post_markers_top, bcells_pre_markers_top)) + 
    scale_fill_gradient2(low="steelblue",mid="white",high = "firebrick") +
    theme(legend.position = "none")
  bcells_heatmaps[[i]] = plot
}

plot_grid(plotlist = bcells_heatmaps, nrow = 3, labels = names(bcells_heatmaps))



##################################
####pre-post aneuploid diplpod####
##################################

Idents(SCTall_fil2) = "minor_clusters2"
EpiBCR15_cells = WhichCells(SCTall_fil2, idents = "Epithelial - BCR15")
Idents(SCTall_fil2) = "timepoint"
pre_cells = WhichCells(SCTall_fil2, idents = "pre")
post_cells = WhichCells(SCTall_fil2, idents = "post")
EpiBCR15_pre_cells = EpiBCR15_cells[EpiBCR15_cells %in% pre_cells]
EpiBCR15_post_cells = EpiBCR15_cells[EpiBCR15_cells %in% post_cells]

Idents(SCTall_fil2) = "ck_new"
SCTall_fil2_EpiBCR15_pre_ck_new_SCT_markers = FindAllMarkers(SCTall_fil2[,EpiBCR15_pre_cells], assay = "SCT", logfc.threshold = 0.2, only.pos = TRUE)
SCTall_fil2_EpiBCR15_post_ck_new_SCT_markers = FindAllMarkers(SCTall_fil2[,EpiBCR15_post_cells], assay = "SCT", logfc.threshold = 0.2, only.pos = TRUE)

SCTall_fil2_EpiBCR15_pre_ck_new_SCT_markers %>% 
  arrange(cluster, -avg_logFC) %>%
  filter(!cluster=="undetermined") %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  group_by(cluster) %>% 
  top_n(50, -p_val_adj) %>% 
  View()

SCTall_fil2_EpiBCR15_post_ck_new_SCT_markers %>% 
  arrange(cluster, -avg_logFC) %>%
  filter(!cluster=="undetermined") %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  group_by(cluster) %>% 
  top_n(50, -p_val_adj) %>% 
  View()

Idents(SCTall_fil2) = "ck_new"
aneuploid_cells = WhichCells(SCTall_fil2, idents = "aneuploid")
diploid_cells = WhichCells(SCTall_fil2, idents = "diploid")
EpiBCR15_anu_cells = EpiBCR15_cells[EpiBCR15_cells %in% aneuploid_cells]
EpiBCR15_dip_cells = EpiBCR15_cells[EpiBCR15_cells %in% diploid_cells]

Idents(SCTall_fil2) = "timepoint"
SCTall_fil2_EpiBCR15_anu_ck_new_SCT_markers = FindAllMarkers(SCTall_fil2[,EpiBCR15_anu_cells], assay = "SCT", logfc.threshold = 0.2, only.pos = TRUE)
SCTall_fil2_EpiBCR15_dip_ck_new_SCT_markers = FindAllMarkers(SCTall_fil2[,EpiBCR15_dip_cells], assay = "SCT", logfc.threshold = 0.2, only.pos = TRUE)

SCTall_fil2_EpiBCR15_anu_ck_new_SCT_markers %>% 
  arrange(cluster, -avg_logFC) %>%
  filter(!cluster=="undetermined") %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!gene %in% SCTall_fil2_EpiBCR15_dip_ck_new_SCT_markers$gene) %>% 
  View()

SCTall_fil2_EpiBCR15_dip_ck_new_SCT_markers %>% 
  arrange(cluster, -avg_logFC) %>%
  filter(!cluster=="undetermined") %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!gene %in% SCTall_fil2_EpiBCR15_anu_ck_new_SCT_markers$gene) %>% 
  View()



##################
####macrophages###
##################

Idents(SCTall_fil2) = "minor_clusters2"
SCTall_fil2_macrophages = subset(SCTall_fil2, idents ="Myeloid - 1")

Idents(SCTall_fil2_macrophages) = "timepoint"

SCTall_fil2_macrophages_timepoint_markers_RNA = FindAllMarkers(SCTall_fil2_macrophages, assay = "RNA", only.pos = TRUE, logfc.threshold = 2)
SCTall_fil2_macrophages_timepoint_markers_SCT = FindAllMarkers(SCTall_fil2_macrophages, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.1)

SCTall_fil2_macrophages_timepoint_markers_SCT %>% 
  arrange(cluster, -avg_logFC) %>% 
  View()



###############
###tcells      ####
###############


Idents(SCTall_fil2) = "minor_clusters2"
tcell_clusters = names(table(SCTall_fil2$minor_clusters2))[str_detect(names(table(SCTall_fil2$minor_clusters2)), "^Tcells")]
SCTall_fil2_tcells = subset(SCTall_fil2, idents =c(tcell_clusters, "NK/NKT/GD"))


Idents(SCTall_fil2_tcells) = "minor_clusters2"
tcell_pre_post_cluster_markers = c()
for(i in levels(Idents(SCTall_fil2_tcells))){
  Idents(SCTall_fil2_tcells) = "minor_clusters2"
  cells = WhichCells(SCTall_fil2_tcells, idents = i)
  Idents(SCTall_fil2_tcells) = "timepoint"
  temp = FindAllMarkers(SCTall_fil2_tcells[,cells], only.pos = TRUE, assay = "SCT", logfc.threshold = 0.1)
  temp$minor_clusters2 = i
  tcell_pre_post_cluster_markers = rbind(tcell_pre_post_cluster_markers,temp)
}

tcell_pre_post_cluster_markers %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  arrange(-avg_logFC) %>% 
  View()
  



Idents(SCTall_fil2_tcells) = "patient"
tcell_pre_post_patient_markers = c()
for(i in levels(Idents(SCTall_fil2_tcells))){
  Idents(SCTall_fil2_tcells) = "patient"
  cells = WhichCells(SCTall_fil2_tcells, idents = i)
  Idents(SCTall_fil2_tcells) = "timepoint"
  temp = FindAllMarkers(SCTall_fil2_tcells[,cells], only.pos = TRUE, assay = "SCT", logfc.threshold = 0.1)
  temp$patient = i
  tcell_pre_post_patient_markers = rbind(tcell_pre_post_patient_markers,temp)
}


tcell_pre_post_patient_markers %>% 
  filter(patient == "BCR15") %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  arrange(cluster, -avg_logFC) %>% 
  View()




##################
####tumor cells###
##################

Idents(SCTall_fil2) = "minor_clusters2"
SCTall_fil2_epi = subset(SCTall_fil2, idents = paste0("Epithelial - ",patients[-9]))


Idents(SCTall_fil2_epi) = "patient"
tumor_prepost_list = list()
for(i in patients[-9]){
  temp = subset(SCTall_fil2_epi, idents=i)
  Idents(temp) = "timepoint"
  tumor_prepost_list[[i]] = FindAllMarkers(temp, only.pos = FALSE, assay = "SCT", logfc.threshold = 0)
}


Idents(SCTall_fil2_epi) = "timepoint"
cells = WhichCells(SCTall_fil2_epi, idents = "post")
Idents(SCTall_fil2_epi) = "patient"
tumor_patient_post_markers = FindAllMarkers(SCTall_fil2_epi[,cells], only.pos = TRUE, assay = "SCT",logfc.threshold = 0.1)


tumor_patient_post_markers %>% 
  filter(cluster == "BCR15") %>% 
  filter(pct.2 < 0.1) %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>% 
  top_n(15, avg_logFC) %>% 
  pull(gene) ->
  BCR15_epi_post_markers_top



VlnPlot(SCTall_fil2_epi, features = c("ESR1","PGR","ERBB2"),
        group.by = "patient", split.plot = TRUE, split.by = "timepoint", 
        pt.size = 0, ncol = 1, cols = rev(treatment_colors)) -> temp
ggsave(temp, filename = "./plots/plots10/receptor_expression.png", width = 5, height = 8)



tumor_prepost_df = bind_rows(tumor_prepost_list, .id = "patient")



tumor_prepost_df %>% 
  filter(cluster=="post") %>% 
  filter(avg_logFC >= 0.3) %>% 
  group_by(gene) %>% 
  summarise(n = n(), pct.2 = mean(pct.2)) %>% 
  arrange(-n) %>% 
  filter(n>=5) %>% 
  arrange(pct.2) %>%
  top_n(20, -pct.2) %>% 
  arrange(-n, pct.2) %>% 
  pull(gene) ->
  tumor_post_avg_top


Idents(SCTall_fil2_epi) = "patient"
tumor_prepost_avg = AverageExpression(SCTall_fil2_epi, features = tumor_post_avg_top, return.seurat = T, add.ident = "timepoint", assays = "SCT", slot = "data")

DoHeatmap(tumor_prepost_avg, features = tumor_post_avg_top, assay = "SCT", slot = "data", disp.max = 2, lines.width = 1, size = 5, 
          group.colors = patient_colors[c(1:5,7:11)]) + 
  scale_fill_gradientn(colors = c("black","firebrick","orange","yellow"), na.value = "white") -> prepost_epi_avg_plot
#ggsave(prepost_epi_avg_plot, filename = "./plots/plots11/prepost_epi_avg.png", width = 6, height = 6)


plot_grid(NULL, NULL, prepost_epi_avg_plot, tumor_prepost_df_gsea_H_plot, 
          labels = c("","","A","B"), nrow = 2, label_size = 20, rel_widths = c(0.6,1), rel_heights = c(0.05,1)) -> tum_expr_fig
ggsave(tum_expr_fig, filename = "./plots/plots11/tum_expr_fig.png", width = 16, height = 8)




#####
#patient specific tumor changes


tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  filter(avg_logFC > 0.1) %>% 
  group_by(gene) %>% 
  mutate(n = n()) %>% 
  filter(n >= 8) %>%
  pull(gene) %>%  
  unique() ->
  tum_shared_genes



tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  filter(avg_logFC > 0.1) %>% 
  group_by(gene) %>% 
  mutate(n = n()) %>% 
  filter(n < 8) %>%
  filter(n >= 4) %>%
  pull(gene) %>%  
  unique() ->
  tum_semi_shared_genes

tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  filter(avg_logFC > 0.1) %>% 
  group_by(gene) %>% 
  mutate(n = n()) %>% 
  filter(n < 4) %>% 
  filter(avg_logFC > 0.4) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  group_by(patient) %>% 
  top_n(20, avg_logFC) %>% 
  arrange(patient, -avg_logFC) %>% 
  View()


#pre_post_uniq_gsea_H = list()
#pre_post_uniq_gsea_C2 = list()
#pre_post_uniq_gsea_C5 = list()
pre_post_uniq_gsea_C6 = list()
pre_post_uniq_gsea_C7 = list()
for(i in patients[c(1,2,3,5,6,8,10,11)]){
  tumor_prepost_df %>% 
    filter(patient ==i) %>% 
    filter(cluster == "post") %>% 
    filter(!gene %in% c(tum_shared_genes)) ->
    temp
  #pre_post_uniq_gsea_H[[i]] = rungsea(temp, fc_cutoff = 0.1, cats = "H", mincat = 3)
  #pre_post_uniq_gsea_C2[[i]] = rungsea(temp, fc_cutoff = 0.1, cats = "C2", mincat = 3)
  #pre_post_uniq_gsea_C5[[i]] = rungsea(temp, fc_cutoff = 0.1, cats = "C5", mincat = 3)
  pre_post_uniq_gsea_C6[[i]] = rungsea(temp, fc_cutoff = 0.1, cats = "C6", mincat = 3)
  pre_post_uniq_gsea_C7[[i]] = rungsea(temp, fc_cutoff = 0.1, cats = "C7", mincat = 3)
  
}
#pre_post_uniq_gsea_H = bind_rows(pre_post_uniq_gsea_H, .id="patient")
#pre_post_uniq_gsea_C2 = bind_rows(pre_post_uniq_gsea_C2, .id="patient")
#pre_post_uniq_gsea_C5 = bind_rows(pre_post_uniq_gsea_C5, .id="patient")
pre_post_uniq_gsea_C6 = bind_rows(pre_post_uniq_gsea_C6, .id="patient")
pre_post_uniq_gsea_C7 = bind_rows(pre_post_uniq_gsea_C7, .id="patient")


pre_post_uniq_gsea_C7 %>% 
  filter(!str_detect(pathway, "TCELL")) %>% 
  filter(!str_detect(pathway, "NK_")) %>% 
  filter(!str_detect(pathway, "MONO")) %>% 
  filter(!str_detect(pathway, "MACRO")) %>% 
  filter(!str_detect(pathway, "BCELL")) %>% 
  filter(!str_detect(pathway, "DC_")) %>%
  filter(!str_detect(pathway, "PBMC")) %>% 
  filter(pval <= 0.05) %>% 
  group_by(patient) %>% 
  arrange(patient, -NES) %>% 
  filter(NES > 0) %>% 
  top_n(5, NES) %>% 
  View()


temp2_H = rungsea(temp, fc_cutoff = 0.2, cats = "H", mincat = 5)
temp2_C2 = rungsea(temp, fc_cutoff = 0.2, cats = "C2", mincat = 5)


temp2_H %>% 
  arrange(-NES) %>% 
  filter(str_detect(pathway, "^WP"))






epi_myeloid_genes = tumor_post_avg_top[tumor_post_avg_top %in% myeloid_post_markers_top]
tumor_post_avg_top[!tumor_post_avg_top %in% myeloid_post_markers_top]


tumor_prepost_df %>% 
  filter(gene %in% ifna_only_genes_top) %>% 
  filter(cluster == "post") %>% 
  View()


epi_myeloid_genes



tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  left_join(selection_df) %>% 
  filter(avg_logFC > 0) %>% 
  group_by(gene, group) %>% 
  summarise(n=n()) %>% 
  spread(key=group, value=n, fill=0) %>%
  filter(low_selection <= 1) %>% 
  filter(high_selection >= 3) %>% 
  arrange(-high_selection, low_selection) %>% 
  View()
  
tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  left_join(selection_df) %>% 
  filter(avg_logFC > 0) %>% 
  group_by(gene, group) %>% 
  summarise(n=n()) %>% 
  spread(key=group, value=n, fill=0) %>%
  filter(low_selection >= 3) %>% 
  filter(high_selection <= 1) %>% 
  arrange(high_selection, -low_selection) %>% 
  pull(gene) ->
  low_sel_tumor_change


tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  complete(gene, patient, fill = list(avg_logFC = 0)) %>% 
  left_join(selection_df) %>%
  filter(group != "unclassified") %>% 
  #filter(avg_logFC > 0) %>% 
  group_by(gene) %>% do(tidy(t.test(avg_logFC ~ group, data = .))) %>% 
  arrange(statistic) %>% 
  dplyr::select(gene, statistic) ->
  temp

temp = pull(temp, statistic)  %>%set_names(pull(temp, gene))


m_df = msigdbr(species = "Homo sapiens", category = "H")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
test = fgsea(pathways = m_list, stats = temp, nperm = 1000)

test %>% 
  arrange(-NES) %>% 
  pull(pathway)

plotEnrichment(pathway = m_list[["HALLMARK_APOPTOSIS"]], stats = temp[!is.na(temp)]) + labs(title="HALLMARK_APOPTOSIS")
plotGseaTable(m_list[test %>% arrange(-NES) %>% pull(pathway)], temp[!is.na(temp)], test, 
              gseaParam=0.5)

tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  left_join(selection_df) %>% 
  mutate(avg_logFC = as.numeric(avg_logFC)) %>% 
  group_by(gene, group) %>% 
  dplyr::summarise(m=mean(avg_logFC, na.rm=T)) %>% 
  spread(key=group, value=m, fill=0) %>%
  arrange(-high_selection, low_selection) %>% 
  View()

tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  left_join(selection_df) %>% 
  mutate(avg_logFC = as.numeric(avg_logFC)) %>% 
  group_by(gene, group) %>% 
  dplyr::summarise(m=mean(avg_logFC, na.rm=T)) %>% 
  spread(key=group, value=m, fill=0) %>%
  arrange(-low_selection,high_selection) %>% 
  View()

tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  left_join(selection_df) %>% 
  mutate(avg_logFC = as.numeric(avg_logFC)) %>% 
  group_by(gene, group) %>% 
  dplyr::summarise(m=mean(avg_logFC, na.rm=T)) %>% 
  spread(key=group, value=m, fill=0) %>%
  arrange(-low_selection,high_selection) %>% 
  mutate(h = high_selection-low_selection, 
         l = low_selection-high_selection) %>% 
  View()


tumor_prepost_df %>% 
  filter(cluster == "pre") %>% 
  left_join(selection_df) %>% 
  mutate(avg_logFC = as.numeric(avg_logFC)) %>% 
  group_by(gene, group) %>% 
  dplyr::summarise(m=mean(avg_logFC, na.rm=T)) %>% 
  spread(key=group, value=m, fill=0) %>%
  arrange(-low_selection,high_selection) %>% 
  View()




tumor_prepost_df %>% 
  filter(cluster == "post") %>% 
  left_join(selection_df) %>% 
  filter(gene == "STAT1")



##########################
###HER2 increase genes####
##########################


SCTall_fil2_epi$HER2

Idents(SCTall_fil2_epi) = "HER2"

her2_tum_post_df = c()
for(i in patients[c(1:6,8,10,11)]){
  marks = FindAllMarkers(SCTall_fil2_epi[,SCTall_fil2_epi$timepoint=="post" & SCTall_fil2_epi$patient == i], only.pos = T, logfc.threshold = 0)
  marks$patient = i
  her2_tum_post_df = rbind(her2_tum_post_df, marks)
}

her2_tum_pre_df = c()
for(i in patients[c(1,2,3,5,6,8,11)]){
  marks = FindAllMarkers(SCTall_fil2_epi[,SCTall_fil2_epi$timepoint=="pre" & SCTall_fil2_epi$patient == i], only.pos = T, logfc.threshold = 0)
  marks$patient = i
  her2_tum_pre_df = rbind(her2_tum_pre_df, marks)
}




her2_tum_post_df %>% 
  filter(cluster == "+") %>% 
  group_by(gene) %>% 
  summarise(n = n(), m = mean(avg_logFC)) %>% 
  filter(n >= 4) %>% 
  arrange(-m) %>% 
  View()



her2_tum_post_df %>% 
  filter(cluster == "+") %>% 
  group_by(gene) %>% 
  filter(n() >= 4) %>% 
  summarise(n = n(), m = mean(avg_logFC), sd = sd(avg_logFC, na.rm = T)) %>% 
  mutate(t = m/sd) %>% 
  mutate(p = 2*pt(t, df = n, lower.tail = F)) %>% 
  arrange(p) %>% 
  pull(gene) -> 
  her2_tum_post_top


her2_tum_post_df %>% 
  filter(cluster == "-") %>% 
  group_by(gene) %>% 
  filter(n() >= 4) %>% 
  summarise(n = n(), m = mean(avg_logFC), sd = sd(avg_logFC, na.rm = T)) %>% 
  mutate(t = m/sd) %>% 
  mutate(p = 2*pt(t, df = n, lower.tail = F)) %>%
  arrange(p) %>% 
  View()


her2_tum_pre_df %>% 
  filter(cluster == "+") %>% 
  group_by(gene) %>% 
  summarise(n = n(), m = mean(avg_logFC)) %>% 
  filter(n >= 3) %>% 
  arrange(-m) %>% 
  View()


her2_tum_pre_df %>% 
  filter(cluster == "+") %>% 
  group_by(gene) %>% 
  filter(n() >= 3) %>% 
  summarise(n = n(), m = mean(avg_logFC), sd = sd(avg_logFC, na.rm = T)) %>% 
  mutate(t = m/sd) %>% 
  mutate(p = 2*pt(t, df = n, lower.tail = F)) %>% 
  arrange(p) %>% 
  pull(gene) -> her2_tum_pre_top


her2_tum_pre_df %>% 
  filter(cluster == "-") %>% 
  group_by(gene) %>% 
  filter(n() >= 3) %>% 
  summarise(n = n(), m = mean(avg_logFC), sd = sd(avg_logFC, na.rm = T)) %>% 
  mutate(t = m/sd) %>% 
  mutate(p = 2*pt(t, df = n, lower.tail = F)) %>% 
  arrange(p) %>% 
  View()

