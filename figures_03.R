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
library(ggrepel)
library(GSVA)
library(GSVAdata)
library(GSEABase)
data(c2BroadSets)
library(org.Hs.eg.db)
data(leukemia)
library(GO.db)
library(msigdbr)
library(fgsea)
library(ggpubr)
source("~/TENX/analysis/MP/MP_project_4/TCR_functions_script.R")
options(future.globals.maxSize = 60000 * 1024^2)
NewSet1 = function(n){return(colorRampPalette(brewer.pal(9,"Set1")[c(1,5,6,3,2,4,8,7,9)])(n))}

gran_clusters2_order = c("Neutrophils - VCAN","Neutrophils - FCGR3A","Mast Cells")
dc_clusters2_order = c("Dendritic Cells - mDCs","Dendritic Cells - cDC1","Dendritic Cells - cDC2")
minor_clusters2_order = c(names(table(SCTall_fil2_epi$minor_clusters2)), "Epithelial - Cycling", "Basal Epithelial", "Luminal Epithelial", "Lymphatic Endothelial", "Vascular Endothelial", "Perivascular Cells", "Fibroblasts", "B - Cells", "Plasma Cells",names(tcell_colors), names(macs_colors), dc_clusters2_order, gran_clusters2_order)
SCTall_fil2$minor_clusters2 = factor(SCTall_fil2$minor_clusters2, levels = minor_clusters2_order)


theme_one_line <- theme_classic() + theme(axis.text = element_text(size=12), 
                                          axis.ticks = element_line(size=0.5),
                                          axis.line = element_line(size=0.5),
                                          axis.title = element_text(size=16),
                                          axis.title.x = element_blank(),
                                          axis.line.x = element_blank(),
                                          axis.ticks.x = element_blank())
theme_two_line <- theme_classic() + theme(axis.text = element_text(size=12), 
                                          axis.ticks = element_line(size=0.5),
                                          axis.line = element_line(size=0.5),
                                          axis.title = element_text(size=16))





sample_colors = c("#94E3D3", "#517D74", "#A882B7", "#3C2544", "#AA936F", "#685230", "#5167c9", "#26305E", "#FC6F60", "#872219", "#ED9D2D", "#774601", "#5E92A4", "#112630", "#a3a3a3", "#4a4a4a", "#d8db6e", "#5d5e30", "#B9FF94", "#5B7D48", "#FF9CB9", "#6B414E")
patient_colors = sample_colors[c(1,3,5,7,9,11,13,15,17,19,21)]

tcell_colors_df = data.frame(cluster = unique(SCTall_fil2_tcells$minor_clusters4), color = c("#2E99C7","#2E64C7","#FF0000","#FF9E9E","#005717","#757067","#00FF44","black","#BF8828","#6E4444","#487554","#C7C2B9","#A13030","#6E0000","#9CFFB6","#36BF5A","#523CAB","#7F3CAB"), stringsAsFactors = F)
#tcell_colors_df = tcell_colors_df[c(2,1,14,13,3,4,10,11,5,16,7,15,17,18,9,6,12,8),]
tcell_colors_df = tcell_colors_df[c(2,14,13,3,4,10,17,9,1,11,5,16,7,15,18,6,12,8),]
tcell_colors = tcell_colors_df$color
names(tcell_colors) = tcell_colors_df$cluster



macs_colors_df = data.frame(cluster = unique(SCTall_fil2_macs$minor_clusters2), color = c("#8ce043","#43e0db","#43a4e0","#4353e0","#a143e0","#e043d0"), stringsAsFactors = F)
macs_colors_df = macs_colors_df[c(1,3,2,6,4,5),]
macs_colors = macs_colors_df$color
names(macs_colors) = macs_colors_df$cluster


gran_colors_df = data.frame(cluster = gran_clusters2_order, color = c("#c9ac8f","#a6c98f","#bb8fc9"), stringsAsFactors = F)
gran_colors = gran_colors_df$color
names(gran_colors) = gran_colors_df$cluster

dc_colors_df = data.frame(cluster = dc_clusters2_order, color = c("#e04343","#e09443","#e0ce43"), stringsAsFactors = F)
dc_colors = dc_colors_df$color
names(dc_colors) = dc_colors_df$cluster



cytokines_df = read_csv("~/code/resources/gene_lists/cytokines.csv")
cytokines$gene = UpdateSymbolList(cytokines_df$gene)
cytokines = cytokines_df$gene[cytokines_df$gene %in% row.names(SCTall_fil2)]
#cytokines = c(cytokines, "GZMA","GZMB","GZMH","GZMK","GZMM","NKG7","GNLY","PRF1","TNF","IFNG")
#cytokines = cytokines[cytokines %in% row.names(SCTall_fil2)]

chemlig_df = read_csv("~/code/resources/gene_lists/chem_ligands.csv")
chemlig_df$gene = UpdateSymbolList(chemlig_df$gene)
chemlig = chemlig_df$gene[chemlig_df$gene %in% row.names(SCTall_fil2)]
chemlig = unique(chemlig)

chemrec_df = read_csv("~/code/resources/gene_lists/chem_receptors.csv")
chemrec_df$gene = UpdateSymbolList(chemrec_df$gene)
chemrec = chemrec_df$gene[chemrec_df$gene %in% row.names(SCTall_fil2)]



receptors_df = read_csv("~/code/resources/gene_lists/greceptors.csv")
greceptors = receptors_df$gene[receptors_df$gene %in% row.names(SCTall_fil2)]

lectins_df = read_csv("~/code/resources/gene_lists/lectins.csv")
lectins = lectins_df$gene[lectins_df$gene %in% row.names(SCTall_fil2)]

CD_df = read_csv("~/code/resources/gene_lists/CDs.csv")
CDs = CD_df$gene[CD_df$gene %in% row.names(SCTall_fil2)]
CDs = CDs[!CDs %in% c("CD8A","CD8B","CD4","CD3E","CD3D","CD3G")]

ISGs = row.names(SCTall_fil2)[str_detect(row.names(SCTall_fil2),"^ISG")]
ISGs = c(ISGs,row.names(SCTall_fil2)[str_detect(row.names(SCTall_fil2),"^IFI")])
ISGs = c(ISGs,"MX1","MX2")


hydrolases_df = read_csv("~/code/resources/gene_lists/hydrolases.csv")
hydrolases = hydrolases_df$gene[hydrolases_df$gene %in% row.names(SCTall_fil2)]


metallopeptidases_df = read_csv("~/code/resources/gene_lists/metallopeptidases.csv")
metallopeptidases = metallopeptidases_df$gene[metallopeptidases_df$gene %in% row.names(SCTall_fil2)]

complement_df = read_csv("~/code/resources/gene_lists/complement.csv")
complement = complement_df$gene[complement_df$gene %in% row.names(SCTall_fil2)]

cystatins_df = read_csv("~/code/resources/gene_lists/cystatins.csv")
cystatins = cystatins_df$gene[cystatins_df$gene %in% row.names(SCTall_fil2)]


metallothioneins_df = read_csv("~/code/resources/gene_lists/Metallothioneins.csv")
metallothioneins = metallothioneins_df$gene[metallothioneins_df$gene %in% row.names(SCTall_fil2)]

RLs_df = read_csv("~/code/resources/gene_lists/receptor_ligands.csv")
RLs = RLs_df$gene[RLs_df$gene %in% row.names(SCTall_fil2)]

#unique(SingleCellSignalR::LRdb$ligand))



revlog_trans <- function(base = exp(1)) {
  trans <- function(x) -log1p(x)
  inv <- function(x) expm1(-x)
  scales::trans_new("revlog1p", trans, inv, domain = c(0, Inf))
}





#################################
#Myeloid minor cluster markers###
##################################

Idents(SCTall_fil2) = "minor_clusters2"
SCTall_fil2_minor_clusters2_macro_markers_SCT = FindAllMarkers(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "^Macro")], assay = "SCT", logfc.threshold = 0.2,only.pos = T)
SCTall_fil2_minor_clusters2_dc_markers_SCT = FindAllMarkers(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "^Dend")], assay = "SCT", logfc.threshold = 0.2,only.pos = T)
SCTall_fil2_minor_clusters2_gran_markers_SCT = FindAllMarkers(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "^Neu|Mast")], assay = "SCT", logfc.threshold = 0.2,only.pos = T)

rbind(SCTall_fil2_minor_clusters2_macro_markers_SCT,
      SCTall_fil2_minor_clusters2_dc_markers_SCT,
      SCTall_fil2_minor_clusters2_gran_markers_SCT) %>% 
  arrange(cluster,-avg_logFC) %>% 
  write_csv("./data/myeloid_minor_cluster_markers.csv")

SCTall_fil2_major_clusters_macro_markers_SCT = FindMarkers(SCTall_fil2, ident.1 = minor_clusters2[str_detect(minor_clusters2, "^Macro")], logfc.threshold = 0.2,only.pos = T)
SCTall_fil2_major_clusters_dc_markers_SCT = FindMarkers(SCTall_fil2, ident.1 = minor_clusters2[str_detect(minor_clusters2, "^Dend")], logfc.threshold = 0.2,only.pos = T)
SCTall_fil2_major_clusters_gran_markers_SCT = FindMarkers(SCTall_fil2, ident.1 = minor_clusters2[str_detect(minor_clusters2, "^Neu|Mast")], logfc.threshold = 0.2,only.pos = T)

SCTall_fil2_major_clusters_macro_markers_SCT$cluster = "Macrophages"
SCTall_fil2_major_clusters_macro_markers_SCT$gene = row.names(SCTall_fil2_major_clusters_macro_markers_SCT)

SCTall_fil2_major_clusters_dc_markers_SCT$cluster = "Dendritic Cells"
SCTall_fil2_major_clusters_dc_markers_SCT$gene = row.names(SCTall_fil2_major_clusters_dc_markers_SCT)

SCTall_fil2_major_clusters_gran_markers_SCT$cluster = "Granulocytes"
SCTall_fil2_major_clusters_gran_markers_SCT$gene = row.names(SCTall_fil2_major_clusters_gran_markers_SCT)



rbind(SCTall_fil2_major_clusters_macro_markers_SCT,
      SCTall_fil2_major_clusters_dc_markers_SCT,
      SCTall_fil2_major_clusters_gran_markers_SCT) %>% 
  arrange(cluster,-avg_logFC) %>% 
  #filter(avg_logFC >= 1) %>% group_by(cluster) %>%  top_n(10, avg_logFC) %>% View()
  write_csv("./data/myeloid_major_cluster_markers.csv")

rbind(SCTall_fil2_major_clusters_macro_markers_SCT,
      SCTall_fil2_major_clusters_dc_markers_SCT,
      SCTall_fil2_major_clusters_gran_markers_SCT) %>% 
  arrange(cluster,-avg_logFC) %>% 
  filter(avg_logFC >= 1) %>% group_by(cluster) %>%  top_n(20, avg_logFC) %>% 
  pull(gene) %>%  unique() ->
  myeloid_major_clusters_markers_top


#######################
####macro lineplot####
#######################


SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  #filter(major_clusters %in% c("Myeloid")) %>% 
  #filter(str_detect(minor_clusters2, "Macro")) %>% 
  group_by(sample, timepoint,patient,cohort, minor_clusters2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  filter(str_detect(minor_clusters2, "Macro")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = minor_clusters2, value = p, fill=0) %>% 
  gather(-timepoint, -patient, key = minor_clusters2, value =p) %>% 
  #mutate(minor_clusters2 = factor(minor_clusters2, levels = )) %>% 
  ggplot() + aes(x = timepoint, y=p, group=patient, color=patient) + 
  geom_point(size=2.5) + geom_line(size=0.8) + 
  scale_color_manual(values = patient_colors) +
  facet_wrap(~minor_clusters2,  nrow = 1) + 
  theme(panel.background = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title=element_blank(), legend.title = element_blank(),
        strip.background = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text = element_text(size=12), legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        panel.spacing.y = unit(2, "lines")) + 
  guides(color = guide_legend(override.aes = list(size = 2, shape=15))) ->
  macro_lineplot_grid

y.grob <- textGrob("Fraction of Immune Cells", 
                   gp=gpar(col="black", fontsize=20), rot=90)

grid.arrange(arrangeGrob(macro_lineplot_grid, left = y.grob, padding = unit(1,"cm"))) -> macro_lineplot
#ggsave(macro_lineplot, filename = "./plots/fig_pieces_02/macro_lineplot.pdf", width=9, height=4.5, device = cairo_pdf)
ggsave(macro_lineplot, filename = "./plots/fig_pieces_02/macro_lineplot_wide.pdf", width=14, height=4, device = cairo_pdf)



SCTall_fil2@meta.data %>%
  #filter(str_detect(minor_clusters2, "Macro")) %>% 
  #filter(major_clusters %in% c("Myeloid")) %>% 
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, minor_clusters2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  filter(str_detect(minor_clusters2, "Macro")) %>% 
  filter(str_detect(minor_clusters2, "Macrophages - CCL18")) %>% 
  ungroup() %>% 
  dplyr::select(-sample,-cohort,-n) %>% 
  spread(key = timepoint, value = p, fill = 0) -> temp
  #temp$pre = 0
  t.test(temp$pre, temp$post, paired = T)


  

  SCTall_fil2_tcells@meta.data %>%
    #filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
    group_by(sample, timepoint,patient,cohort, minor_clusters4) %>%
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(sample) %>% 
    mutate(p = n/sum(n)) %>% 
    filter(str_detect(minor_clusters4, "Tcells")) %>% 
    filter(str_detect(minor_clusters4, "Tcells - CD8:CCL4")) %>% 
    ungroup() %>% 
    select(-sample,-cohort,-n) %>% 
    spread(key = timepoint, value = p, fill = 0) -> temp 
  t.test(temp$pre, temp$post, paired = T)
  
  
  
#Idents(SCTall_fil2_tcells) = "minor_clusters4"
#temp = AverageExpression(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD8")], assays = c("RNA","SCT"), features = tfs[tfs%in%row.names(SCTall_fil2_tcells)], return.seurat = T)







SCTall_fil2_macs_minor_clusters2_markers = FindAllMarkers(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2,"Macrophages")], only.pos = T, assay = "SCT", logfc.threshold = 0.1)

SCTall_fil2_macs_minor_clusters2_markers %>% 
  filter(avg_logFC >= 0.1) %>% 
  pull(gene) %>% 
  unique() ->
  SCTall_fil2_macs_minor_clusters2_markers_high



################
###macro umap###
################

#SCTall_fil2_macs = RunUMAP(SCTall_fil2_macs, n.neighbors = 150,
#                             features = c(SCTall_fil2_macs_minor_clusters2_markers_high), 
#                             reduction.name = "UMAP_top", reduction.key = "umap_")
#DimPlot(SCTall_fil2_macs, reduction = "UMAP_top", cols = NewSet1(6), pt.size = 1)
#DimPlot(SCTall_fil2_macs, reduction = "UMAP_top", cols = patient_colors, pt.size = 1, group.by="patient")


SCTall_fil2_macs = FindVariableFeatures(SCTall_fil2_macs)
SCTall_fil2_macs = RunPCA(SCTall_fil2_macs, npcs = 200)
SCTall_fil2_macs = RunUMAP(SCTall_fil2_macs, dims=1:70, n.neighbors = 50, spread = 0.25)
DimPlot(SCTall_fil2_macs, reduction = "umap", cols = macs_colors, pt.size = 1)
DimPlot(SCTall_fil2_macs, reduction = "umap",group.by = "patient", cols = patient_colors, pt.size = 1)



macs_umap_df = data.frame(umap1 = SCTall_fil2_macs@reductions$umap@cell.embeddings[,1], umap2 = SCTall_fil2_macs@reductions$umap@cell.embeddings[,2], cluster = SCTall_fil2_macs$minor_clusters2, timepoint = SCTall_fil2_macs$timepoint, stringsAsFactors = F)
macs_umap_df = left_join(macs_umap_df, macs_colors_df)
macs_umap_df$color = as.character(macs_umap_df$color)



macs_umap_df %>% 
  ggplot() + aes(x=umap1, y=umap2, color=color)  +
  geom_point(aes(alpha=timepoint), size=0.01) + 
  theme_void() + theme(axis.title = element_text(size=16), axis.title.y = element_text(angle = 90)) +
  scale_color_identity("cluster", labels = names(macs_colors), breaks = macs_colors,
                    guide = "legend") + scale_alpha_manual(values = c(0.05,1)) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 20, 10)) ->
  plot_macs_umap_post


macs_umap_df %>% 
  ggplot() + aes(x=umap1, y=umap2, color=color)  +
  geom_point(aes(alpha=timepoint), size=0.01) + 
  theme_void() + theme(axis.title = element_text(size=16), axis.title.y = element_text(angle = 90)) + scale_alpha_discrete(guide = FALSE) +
  scale_color_identity("cluster", labels = names(macs_colors), breaks = macs_colors, guide = "legend") + 
  guides(color = guide_legend(override.aes = list(size=5))) + theme(legend.title = element_blank()) ->
  plot_macs_umap_legend
plot_macs_umap_legend = get_legend(plot_macs_umap_legend)
plot_grid(plot_macs_umap_legend) -> plot_macs_umap_legend_plot
ggsave(plot = plot_macs_umap_legend_plot, filename = "./plots/fig_pieces_02/macs_umap_legend.pdf", width = 2, height = 4.5, dpi = 300, device = cairo_pdf)



macs_umap_df %>% 
  mutate(timepoint = factor(timepoint, levels = c("post","pre"))) %>% 
  ggplot() + aes(x=umap1, y=umap2, color=color)  +
  geom_point(aes(alpha=timepoint), size=0.005) + scale_alpha_manual(values = c(0.05,1)) +
  theme_void() + theme(axis.title = element_text(size=16), axis.title.y = element_text(angle = 90)) +
  scale_color_identity("cluster", labels = names(macs_colors), breaks = macs_colors,
                       guide = "legend") + theme(legend.position = "none", plot.margin = margin(0, 0, 20, 10)) ->
  plot_macs_umap_pre

plot_grid(plot_macs_umap_pre,plot_macs_umap_post) -> plot_macs_umap
ggsave(plot = plot_macs_umap, filename = "./plots/fig_pieces_02/macs_umap.png", width = 6, height = 2.5, dpi = 300)


DimPlot(SCTall_fil2_macs, group.by = "patient", cols = patient_colors)-> plot_macs_patient_umap
ggsave(plot = plot_macs_patient_umap, filename = "~/projects/Breast_cancer_radiation/fig_pieces_02/plot_macs_patient_umap.png", width = 6, height = 4, dpi = 300)

DimPlot(SCTall_fil2_macs, group.by = "timepoint", cols = treatment_colors)-> plot_macs_timepoint_umap
ggsave(plot = plot_macs_timepoint_umap, filename = "./plots/fig_pieces_02/plot_macs_timepoint_umap.png", width = 6, height = 2.5, dpi = 300)



###########################
###myeloid umap_barchart###
###########################

SCTall_fil2_macs@meta.data %>% 
  ggplot() + aes(x=timepoint, fill=minor_clusters2) + geom_bar(position="fill") + 
  scale_fill_manual(values = macs_colors) +   theme_minimal() +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position = "none") + theme(panel.grid = element_blank()) + 
  theme(axis.title = element_blank())  + theme(axis.text.y=element_blank()) +
  theme(axis.text.x = element_text(size=12)) -> myeloid_umap_bar

ggsave(plot = myeloid_umap_bar, filename = "./plots/fig_pieces_02/myeloid_umap_bar.pdf", width = 1.5, height = 4)







####################
###mac heatmap######
####################


rbind(SCTall_fil2_macs_minor_clusters2_markers) %>% 
  mutate(cluster = factor(cluster, levels = as.character(macs_colors_df$cluster))) %>% 
  arrange(cluster, avg_logFC) %>% 
  group_by(cluster) %>% 
  top_n(15, avg_logFC) %>% 
  pull(gene) %>% 
  unique() -> SCTall_fil2_macs_minor_clusters2_markers_top
  
  #filter(!gene %in% RLs) %>% 
  filter(!gene %in% sort(unique(SingleCellSignalR::LRdb$ligand))) %>% 
  #filter(!gene %in% cytokines) %>% 
  #filter(!gene %in% tfs) %>% 
  filter(!gene %in% ISGs) %>% 
  #filter(!gene %in% sort(unique(SingleCellSignalR::LRdb$receptor))) %>% 
  #View()
  #filter(!gene %in% c(tfs,cytokines,checkpoint_genes)) %>% 
  #filter(!gene %in% greceptors) %>% 
  #filter(!gene %in% lectins) %>%  
  #filter(!gene %in% CDs) %>% 
  #filter(!gene %in% hydrolases) %>% 
  #filter(!gene %in% metallopeptidases) %>%
  #filter(!gene %in% complement) %>% 
  #filter(!gene %in% cystatins) %>% 
  #filter(!gene %in% metallothioneins) %>% 
  View()



#top_gene_filter_macs = function(genes_list,n=5){
  
  SCTall_fil2_macs_minor_clusters2_markers %>% 
    mutate(cluster = factor(cluster, levels = as.character(macs_colors_df$cluster))) %>% 
    arrange(cluster) %>% 
    #filter(cluster %in% c("NK Cells", "Tcells - Cycling","Tcells - GD")) %>% 
    filter(gene %in% genes_list) %>% 
    group_by(cluster) %>% 
    top_n(n,avg_logFC) %>% 
    pull(gene) %>% 
    unique() ->
    return_genes
  
  #SCTall_fil2_tcells_minor_clusters4_markers_CD4 %>% 
  #  mutate(cluster = factor(cluster, levels = as.character(tcell_colors_df$cluster))) %>% 
  #  arrange(cluster) %>% 
  #  filter(gene %in% genes_list) %>% 
  #  group_by(cluster) %>% 
  #  top_n(n,avg_logFC) %>% 
  #  pull(gene) %>% 
  #  unique() %>% 
   # c(return_genes) %>% 
   # unique() ->
   # return_genes
  
  #SCTall_fil2_tcells_minor_clusters4_markers_CD8 %>% 
  #  mutate(cluster = factor(cluster, levels = as.character(tcell_colors_df$cluster))) %>% 
  #  arrange(cluster) %>% 
  #  filter(gene %in% genes_list) %>% 
  #  group_by(cluster) %>% 
  #  top_n(n,avg_logFC) %>% 
   # pull(gene) %>% 
  #  c(return_genes) %>% 
   # unique() ->
   # return_genes
  
#  return(return_genes)
#}

#
#mac_genes_RLs = top_gene_filter_macs(unique(SingleCellSignalR::LRdb$ligand), 10)
#mac_genes_ISGs = top_gene_filter_macs(ISGs, 10)


#mac_genes_allgroups = c(mac_genes_RLs, mac_genes_ISGs)
#mac_genes_allgroups = unique(mac_genes_allgroups)



#tcell_genes_tfs = rev(tcell_genes_tfs)

SCTall_fil2_macs$minor_clusters2 = factor(SCTall_fil2_macs$minor_clusters2, levels = as.character(macs_colors_df$cluster))
Idents(SCTall_fil2_macs) = "minor_clusters2"

SCTall_fil2_macs_avg = AverageExpression(SCTall_fil2_macs, features = SCTall_fil2_macs_minor_clusters2_markers_top, return.seurat = T)

DoHeatmap(SCTall_fil2_macs_avg, features = SCTall_fil2_macs_minor_clusters2_markers_top, assay = "SCT", slot = "scale.data", label = F,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = macs_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm"))-> macs_heatmap


library(gtable)
g <- ggplotGrob(macs_heatmap)
s <- gtable_filter(g, 'axis-l|ylab', trim=F)  # use trim depending on need
grid.newpage()
grid.draw(s)


DoHeatmap(SCTall_fil2_macs_avg, features = SCTall_fil2_macs_minor_clusters2_markers_top, assay = "SCT", slot = "scale.data", label = F,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = macs_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.text.y = element_blank()) -> macs_heatmap


plot_grid(s,macs_heatmap, nrow=1, align = "h", rel_widths = c(5,8)) -> macs_heatmap_bottom
ggsave(macs_heatmap_bottom, filename = "./plots/fig_pieces_02/macs_heatmap_bottom.pdf", width=2, height=9)
#ggsave(tcells_heatmap_bottom, filename = "./plots/fig_pieces_01/tcells_heatmap_bottom.png", width=2, height=9)





myeloid_major_clusters_markers_top
Idents(SCTall_fil2) = "minor_clusters2"
SCTall_fil2_myeloid_avg = AverageExpression(SCTall_fil2[,(SCTall_fil2$major_clusters=="Myeloid")], features = myeloid_major_clusters_markers_top, return.seurat = T)

myeloid_major_clusters_markers_top_clust = myeloid_major_clusters_markers_top[hclust(dist(SCTall_fil2_myeloid_avg@assays$SCT@scale.data))$order]

DoHeatmap(SCTall_fil2_myeloid_avg, features = myeloid_major_clusters_markers_top_clust, assay = "SCT", slot = "scale.data", label = T,
          disp.min = -2, disp.max = 2, draw.lines = F,  raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.text.y = element_blank()) -> myeloid_heatmap




Idents(SCTall_fil2) = "minor_clusters2"
SCTall_fil2_myeloid_avg_cyto = AverageExpression(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2,"Mac")], features = c(cytokines), return.seurat = T)

Idents(SCTall_fil2_myeloid_avg_cyto) = factor(Idents(SCTall_fil2_myeloid_avg_cyto), levels =  names(macs_colors)[c(1,3,2,4,5,6)])
SCTall_fil2_myeloid_avg_cyto = SCTall_fil2_myeloid_avg_cyto[names(which(rowSums(SCTall_fil2_myeloid_avg_cyto@assays$SCT@data) >= 0.01)),]
#h = hclust(dist(SCTall_fil2_myeloid_avg_cyto@assays$SCT@data))
#cyto_fil = names(cutree(h, k = 6))
#[cutree(h, k = 6) != 1]
#myeloid_cytokines_clust = cyto_fil[hclust(dist(SCTall_fil2_myeloid_avg_cyto@assays$SCT@data[cyto_fil,]))$order]

data.frame(SCTall_fil2_myeloid_avg_cyto@assays$SCT@data) %>% 
  rownames_to_column("gene") %>% 
  gather(-gene, key = "cluster", value = "value") %>% 
  mutate(cluster = factor(cluster, levels = str_replace(names(macs_colors), "Macrophages - ","Macrophages...")[c(1,3,2,4,5,6)])) %>% 
  group_by(gene) %>% 
  slice(which.max(value)) %>% 
  arrange(cluster, -value) %>% 
  pull(gene) ->
  myeloid_cytokines_top

DoHeatmap(SCTall_fil2_myeloid_avg_cyto, features = myeloid_cytokines_top, assay = "SCT", slot = "scale.data", label = T,
          disp.min = -2, disp.max = 2, draw.lines = F,  raster = F, group.colors = macs_colors[c(1,3,2,4,5,6)]) + NoLegend() + 
  scale_x_discrete(limits = names(macs_colors)[c(1,3,2,4,5,6)]) + 
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) -> mac_cytokine_heatmap





SCTall_fil2_myeloid_avg_hla2 = AverageExpression(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2,"Mac")], features = c(HLA_2_genes), return.seurat = T)

Idents(SCTall_fil2_myeloid_avg_hla2) = factor(Idents(SCTall_fil2_myeloid_avg_hla2), levels =  names(macs_colors)[c(1,3,2,4,5,6)])
SCTall_fil2_myeloid_avg_hla2 = SCTall_fil2_myeloid_avg_hla2[names(which(rowSums(SCTall_fil2_myeloid_avg_hla2@assays$SCT@data) >= 0.01)),]
#h = hclust(dist(SCTall_fil2_myeloid_avg_hla2@assays$SCT@data))
#hla2_fil = names(cutree(h, k = 6))
#[cutree(h, k = 6) != 1]
#myeloid_hla2kines_clust = hla2_fil[hclust(dist(SCTall_fil2_myeloid_avg_hla2@assays$SCT@data[hla2_fil,]))$order]

data.frame(SCTall_fil2_myeloid_avg_hla2@assays$SCT@data) %>% 
  rownames_to_column("gene") %>% 
  gather(-gene, key = "cluster", value = "value") %>% 
  mutate(cluster = factor(cluster, levels = str_replace(names(macs_colors), "Macrophages - ","Macrophages...")[c(1,3,2,4,5,6)])) %>% 
  group_by(gene) %>% 
  slice(which.max(value)) %>% 
  arrange(cluster, -value) %>% 
  pull(gene) ->
  myeloid_hla2_top

DoHeatmap(SCTall_fil2_myeloid_avg_hla2, features = myeloid_hla2_top, assay = "SCT", slot = "scale.data", label = T,
          disp.min = -2, disp.max = 2, draw.lines = F,  raster = F, group.colors = macs_colors[c(1,3,2,4,5,6)]) + NoLegend() + 
  scale_x_discrete(limits = names(macs_colors)[c(1,3,2,4,5,6)]) + 
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) -> mac_hla2_heatmap


mac_hla2_heatmap + theme(plot.margin = margin(t = 3, b = 0, l = 0, r = 4, "cm"))


#####dotplots
chemrec_top = names(which(rowSums(SCTall_fil2_macs[chemrec,]) >= 250))
DotPlot(SCTall_fil2_macs, assay = "SCT", features = chemrec_top, scale = F) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) -> temp
ggsave(temp, filename = "~/projects/Breast_cancer_radiation/sup_figure_pieces/mac_dotplot_rec.pdf", device = cairo_pdf, width = 6, height = 4)

chemlig_top = names(which(rowSums(SCTall_fil2_macs[chemlig,]) >= 500))
DotPlot(SCTall_fil2_macs, assay = "SCT", features = chemlig_top, scale = F) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))-> temp
ggsave(temp, filename = "~/projects/Breast_cancer_radiation/sup_figure_pieces/mac_dotplot_lig.pdf", device = cairo_pdf, width = 9.5, height = 4)



DotPlot(SCTall_fil2_macs, assay = "SCT",scale = F, features = HLA_2_genes) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) -> temp
ggsave(temp, filename = "~/projects/Breast_cancer_radiation/sup_figure_pieces/mac_dotplot_hla2.pdf", device = cairo_pdf, width = 7.5, height = 4)





##############
####HLA#######
##############

treatment_colors = c("royalblue4","goldenrod")
names(treatment_colors) = c("pre","post")

SCTall_fil2@meta.data %>% 
  filter(major_clusters == "Myeloid") %>% 
  filter(str_detect(minor_clusters2, "Macrophage") | str_detect(minor_clusters2, "Dendritic")) %>% 
  ggplot() + aes(x=patient, y=HLA.II.1, color = timepoint, fill=NULL) + ylab("HLA Class II Expression") +
  geom_boxplot(outlier.color = NA, size=1) + theme_one_line + scale_color_manual(values = c("royalblue4","goldenrod")) + 
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 12)) -> myeloid_HLA_plot

ggsave(myeloid_HLA_plot, filename = "./plots/fig_pieces_02/myeloid_HLA.pdf",  width = 12, height = 4.5, dpi = 300, device = cairo_pdf)


SCTall_fil2@meta.data %>% 
  filter(major_clusters == "Myeloid") %>% 
  filter(str_detect(minor_clusters2, "Macrophage") | str_detect(minor_clusters2, "Dendritic")) %>% 
  dplyr::select(patient.2, timepoint, HLA.II.1) %>%
  group_by(patient.2, timepoint) %>% 
  summarise(HLA.II.1 = list(HLA.II.1)) %>% 
  spread(timepoint, HLA.II.1) %>% 
  group_by(patient.2) %>% 
  mutate(p_value = wilcox.test(unlist(pre), unlist(post))$p.value,
         t_value = wilcox.test(unlist(pre), unlist(post))$statistic) %>% 
  View()
  

  

################################
####IHC/vectra/freq CD68########
################################



IMT_TME_merged %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  filter(marker == "CD68") %>% 
  tidyr::spread(timepoint, density) %>% 
  summarise(ttest = list(wilcox.test(pre, post,paired = TRUE))) %>% 
  pull(ttest) ->
  temp

#p = 0.0002185

IMT_TME_merged %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  filter(marker == "CD68") %>% 
  ggplot() + aes(x=timepoint, group=patient, y = density) + 
  geom_point(size=2) + geom_line(size=0.5) + theme_one_line + ylab("Density of CD68+ (cells/mm*2)")
  
  





IMT_TME_merged %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  filter(marker == "CD68") %>% 
  ggplot() + aes(x=timepoint, group=patient, y = density) + 
  geom_point(size=3) + geom_line(size=0.75) + 
  theme_one_line + ylab("Density of CD68+ (cells/mm^2)") +
  scale_x_discrete(expand = c(0,0.3)) -> IHC_cd68_density_plot

#left_join(id_key[,c("IMT_id", "patient")])

id_key %>% 
  gather(key = "timepoint", value = "Accession", pre_bx_id, post_surgery_id) %>% 
  mutate(timepoint = str_split(timepoint, "_")[[1]][1]) %>% 
  dplyr::select(c("Accession", "patient","timepoint")) %>% 
  right_join(vectra_merged) %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  ggplot() + aes(x=timepoint, group=patient, y = `Total Macrophages (n/mm2)`) + 
  #ggplot() + aes(x=timepoint, group=patient, y = mac_fraction) + 
  geom_point(size=3) + geom_line(size=0.75) + theme_one_line + ylab("Density of CD68+ (cells/mm^2)") +
  scale_x_discrete(expand = c(0,0.3)) -> vectra_cd68_density_plot

id_key %>% 
  gather(key = "timepoint", value = "Accession", pre_bx_id, post_surgery_id) %>% 
  mutate(timepoint = str_split(timepoint, "_")[[1]][1]) %>% 
  dplyr::select(c("Accession", "patient","timepoint")) %>% 
  right_join(vectra_merged) %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  #ggplot() + aes(x=timepoint, group=patient, y = `Total Macrophages (n/mm2)`) + 
  ggplot() + aes(x=timepoint, group=patient, y = mac_fraction) + 
  geom_point(size=3) + geom_line(size=0.75) +
  theme_one_line + ylab("Frequency of CD68+ cells") +
  scale_x_discrete(expand = c(0,0.3)) -> vectra_cd68_freq_plot



SCTall_fil2@meta.data %>%
  group_by(patient, timepoint, minor_clusters2) %>% 
  summarise(n = n()) %>% 
  mutate(f = n/sum(n)) %>% 
  #filter(str_detect(minor_clusters2, "Macro|Dendritic")) %>% 
  filter(str_detect(minor_clusters2, "Macro")) %>% 
  #filter(str_detect(minor_clusters2, "Neu")) %>% 
  summarise(f = sum(f)) %>% 
  ungroup() %>% 
  complete(patient, nesting(timepoint), fill = list(f=0)) %>% 
  ggplot() + aes(x=timepoint, group=patient, y = f) + 
  geom_point(size=3) + geom_line(size=0.75) +
  theme_one_line + ylab("Frequency of CD68+ clusters") +
  scale_x_discrete(expand = c(0,0.3)) -> sc_cd68_clus_freq_plot
  
  

SCTall_fil2@meta.data %>%
  group_by(patient, timepoint, minor_clusters2) %>% 
  summarise(n = n()) %>% 
  mutate(f = n/sum(n)) %>% 
  #filter(str_detect(minor_clusters2, "Macro|Dendritic")) %>% 
  filter(str_detect(minor_clusters2, "Macro")) %>% 
  #filter(str_detect(minor_clusters2, "Neu")) %>% 
  summarise(f = sum(f)) %>% 
  ungroup() %>% 
  complete(patient, nesting(timepoint), fill = list(f=0)) %>% 
  spread(key = timepoint, value = f) %>% 
  summarise(ttest = list(wilcox.test(pre, post,paired = TRUE))) %>% 
  pull(ttest)
#0.042

#SCTall_fil2$cd68 = if_else(SCTall_fil2@assays$RNA@counts["CD68",] >= 1, 1, 0)
SCTall_fil2@meta.data %>%
  group_by(patient, timepoint, cd68) %>% 
  summarise(n = n()) %>% 
  mutate(f = n/sum(n)) %>% 
  filter(cd68 == 1) %>% 
  dplyr::select(patient, timepoint, f) %>% 
  ggplot() + aes(x=timepoint, group=patient, y = f) + 
  geom_point(size=3) + geom_line(size=0.75) + 
  theme_one_line + ylab("Frequency of CD68+ cells") +
  scale_x_discrete(expand = c(0,0.3)) -> 
  sc_cd68_pos_freq_plot


  
plot_grid(sc_cd68_clus_freq_plot,sc_cd68_pos_freq_plot, vectra_cd68_freq_plot, vectra_cd68_density_plot,IHC_cd68_density_plot, nrow=1) %>% 
  ggsave(filename = "./plots/fig_pieces_02/myeloid_freq.pdf",  width = 18, height = 3.5, dpi = 300, device = cairo_pdf)




id_key %>% 
  gather(key = "timepoint", value = "Accession", pre_bx_id, post_surgery_id) %>% 
  mutate(timepoint = str_split(timepoint, "_")[[1]][1]) %>% 
  dplyr::select(c("Accession", "patient","timepoint")) %>% 
  right_join(vectra_merged) %>% 
  dplyr::select(patient, `Total Macrophages (n/mm2)`, timepoint) %>% 
  tidyr::spread(timepoint, `Total Macrophages (n/mm2)`) %>%
  summarise(ttest = list(wilcox.test(pre, post,paired = TRUE))) %>% 
  pull(ttest)

#p = 0.00618 wilcox
#p = 0.012

id_key %>% 
  gather(key = "timepoint", value = "Accession", pre_bx_id, post_surgery_id) %>% 
  mutate(timepoint = str_split(timepoint, "_")[[1]][1]) %>% 
  dplyr::select(c("Accession", "patient","timepoint")) %>% 
  right_join(vectra_merged) %>% 
  dplyr::select(patient, mac_fraction, timepoint) %>% 
  tidyr::spread(timepoint, mac_fraction) %>%
  summarise(ttest = list(t.test(pre, post,paired = TRUE))) %>% 
  pull(ttest)





































######################
###cluster markers####
######################


SCTall_fil2_minor_clusters2_macro_markers_SCT = FindAllMarkers(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "^Macro")], assay = "SCT", logfc.threshold = 0.2,only.pos = T)


SCTall_fil2_minor_clusters2_macro_markers_SCT_gsea_C7 = rungsea(SCTall_fil2_minor_clusters2_macro_markers_SCT, cats = "C7", mincat = 5)
SCTall_fil2_minor_clusters2_macro_markers_SCT_gsea_C7 %>% 
  filter((str_detect(pathway, "TCELL")|str_detect(pathway, "TREG")|str_detect(pathway, "EFF")|str_detect(pathway, "CONV")|str_detect(pathway, "MEM")|str_detect(pathway, "CD[48]"))) %>%  
  filter(!(str_detect(pathway, "DC"))) %>%  
  filter(!(str_detect(pathway, "BCELL") & (!str_detect(pathway,"TCELL")))) %>% 
  arrange(-NES) %>% 
  #filter(size <= 10) %>% 
  #filter(!str_detect(pathway, "CD4")) %>% 
  #filter(!str_detect(pathway, "TREG")) %>% 
  View()






          



##################
###M1 M2 scores###
##################



bagaev_sigs_df = read_csv("~/code/resources/gene_lists/bagaev_mac_sigs.csv")
cheng_sigs = read_csv("~/code/resources/gene_lists/cheng_mac_sigs.csv")
yang_sigs = read_csv("~/code/resources/gene_lists/yang_mac_sigs.csv")

names(bagaev_sigs_df) = paste0("bagaev.",names(bagaev_sigs_df))
names(cheng_sigs) = paste0("cheng.",names(cheng_sigs))
names(yang_sigs) = paste0("yang.",names(yang_sigs))

bagaev_sigs_df %>% 
  gather(key = "signature", value = "gene") %>% 
  rbind(cheng_sigs %>% 
          gather(key = "signature", value = "gene") ) %>% 
  rbind(yang_sigs %>% 
          gather(key = "signature", value = "gene")) %>% 
  filter(!is.na(gene)) %>% 
  filter(gene %in% row.names(SCTall_fil2_myeloid)) %>% 
  base::split(f = as.factor(.$signature), drop = T) ->
  all_mac_sigs_list

all_mac_sigs_list = lapply(all_mac_sigs_list, function(x) return(x %>% pull(gene)))

mac_sig_names = names(all_tcell_sigs_list)[c(10,11,13,16:25,28,29,31)]
all_mac_sigs_list = c(all_mac_sigs_list, all_tcell_sigs_list[mac_sig_names])


mac_sig_scores = AddModuleScore(SCTall_fil2_myeloid[,str_detect(SCTall_fil2_myeloid$minor_clusters2, "Macrophages")], features = all_mac_sigs_list, nbin = 18)
names(mac_sig_scores@meta.data)[str_detect(names(mac_sig_scores@meta.data), "^Cluster")] = names(all_mac_sigs_list)
mac_sig_scores = mac_sig_scores@meta.data





mac_sig_scores %>% 
  filter(!str_detect(minor_clusters2,"DC")) %>% 
  group_by(patient, timepoint) %>% 
  mutate(n=n()) %>% 
  filter(n >= 30) %>% 
  dplyr::select(-n) %>% 
  group_by(patient) %>% 
  filter(length(unique(timepoint)) == 2) %>% 
  ungroup() %>% 
  dplyr::select(cell.names, patient, timepoint, minor_clusters2,starts_with("bagaev"), starts_with("cheng"), starts_with("yang"),starts_with("azizi"),starts_with("meta")) %>% 
  gather(-cell.names, -patient, -timepoint , -minor_clusters2, key="signature", value = "score") %>% 
  #filter(signature %in% mac_sig_names) %>% 
  group_by(patient, timepoint, signature) %>% 
  summarise(m = median(score)) %>% 
  ungroup() %>% 
  spread(key = timepoint, value = m, fill = 0) %>% 
  mutate(d = post - pre) %>% 
  group_by(signature) %>% 
  mutate(p_value = t.test(unlist(d))$p.value,
         t_value = t.test(unlist(d))$statistic) %>%
  filter(p_value <= 0.05) %>%  
  mutate(m = median(d), sd = sd(d)) %>% 
  dplyr::select(-pre, -post, -d, -patient) %>% 
  mutate(dir = if_else(m >=0, "pos","neg")) %>% 
  ungroup() %>% 
  unique() %>% 
  top_n(18, (abs(m)-sd)) %>% 
  top_n(15, abs(m)) %>% 
  #top_n(15, -p_value) %>% 
  mutate(signature = str_split(signature, "\\.", simplify = T)[,2]) %>% 
  ggplot() + aes(x=reorder(signature, m), y = m, fill=dir) + geom_bar(stat = "identity") + 
  geom_linerange(aes(ymin = m-sd, ymax = m+sd)) +
  scale_fill_manual(values = c("steelblue","firebrick")) +
  theme_minimal() + theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(-0.5,0.5), expand = c(0,0)) +
  ylab("Median Change in Score") +
  theme(axis.line.x = element_line(size=0.5)) +
  theme(axis.title.x = element_text(size=16), axis.text=element_text(size=12)) +
  theme(panel.grid = element_blank(), legend.position = "none") + coord_flip() -> macs_sig_plot
ggsave(macs_sig_plot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_02/macs_sig_plot.pdf", width=8, height=6)




mac_sig_scores %>% 
  filter(!str_detect(minor_clusters2,"DC")) %>% 
  group_by(minor_clusters2) %>% 
  mutate(n=n()) %>% 
  filter(n >= 30) %>% 
  dplyr::select(-n) %>% 
  dplyr::select(cell.names, patient, timepoint, minor_clusters2,starts_with("azizi"), starts_with("meta"), starts_with("andreatta"),starts_with("yang"),starts_with("cheng"), starts_with("bagaev")) %>% 
  gather(-cell.names, -patient, -timepoint , -minor_clusters2, key="signature", value = "score") %>% 
  filter(signature %in% names(all_mac_sigs_list)[c(9,6,42,43,2,40,37,11,20)]) %>% 
  group_by(minor_clusters2, signature) %>% 
  mutate(m = median(score), sd = sd(score)) %>%
  ungroup() %>% 
  filter() %>% 
  ggplot() + aes(x=minor_clusters2 , y = score, fill=minor_clusters2) + geom_boxplot(outlier.colour = NA) + 
  facet_wrap(~signature, scales = "free_y") + theme_real_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust =1, vjust = 0.5)) + 
  scale_fill_manual(values = macs_colors) 




yang_clus_sigs_df = read_csv("~/code/resources/gene_lists/yang_cluster_sigs.csv")

yang_clus_sigs_df %>% 
  filter(!is.na(gene)) %>% 
  filter(gene %in% row.names(SCTall_fil2_myeloid)) %>% 
  mutate(signature = paste0("yang.",signature)) %>% 
  base::split(f = as.factor(.$signature), drop = T) ->
  yang_clus_sigs_list

yang_clus_sigs_list = lapply(yang_clus_sigs_list, function(x) return(x %>% pull(gene)))


mac_clus_sig_scores = AddModuleScore(SCTall_fil2_myeloid[,str_detect(SCTall_fil2_myeloid$minor_clusters2, "Macrophages")], features = yang_clus_sigs_list, nbin = 18)
names(mac_clus_sig_scores@meta.data)[str_detect(names(mac_clus_sig_scores@meta.data), "^Cluster")] = names(yang_clus_sigs_list)
mac_clus_sig_scores = mac_clus_sig_scores@meta.data


mac_clus_sig_scores %>% 
  filter(!str_detect(minor_clusters2,"DC")) %>% 
  group_by(minor_clusters2) %>% 
  mutate(n=n()) %>% 
  filter(n >= 30) %>% 
  dplyr::select(-n) %>% 
  dplyr::select(cell.names, patient, timepoint, minor_clusters2,starts_with("yang")) %>% 
  gather(-cell.names, -patient, -timepoint , -minor_clusters2, key="signature", value = "score") %>% 
  #filter(signature %in% names(all_mac_sigs_list)[c(9,42,43,2,40,37,11,20)]) %>% 
  group_by(minor_clusters2, signature) %>% 
  mutate(m = median(score), sd = sd(score)) %>%
  ungroup() %>% 
  filter() %>% 
  ggplot() + aes(x=minor_clusters2 , y = score) + geom_boxplot(outlier.colour = NA) + 
  facet_wrap(~signature, scales = "free_y") + theme(axis.text.x = element_text(angle=90, hjust =1, vjust = 0.5))







M1_M2_df = read_csv("~/code/resources/gene_lists/M1_M2_sigs.csv")


M1_M2_df %>% 
  gather(key = "sig", value = "gene") %>% 
  filter(!is.na(gene)) %>% 
  filter(gene %in% row.names(SCTall_fil2)) %>% 
  base::split(f = as.factor(.$sig)) ->
  M1_M2_sigs_list

M1_M2_sigs_list = lapply(M1_M2_sigs_list, function(x) return(x %>% pull(gene)))

SCTall_fil2 = AddModuleScore(SCTall_fil2, features = M1_M2_sigs_list, name = "M1M2_sigs.")

M1M2_sig_rename = c("M1M2_sigs.1"="M1",
                     "M1M2_sigs.2"="M2")

names(SCTall_fil2@meta.data)[match(names(M1M2_sig_rename), names(SCTall_fil2@meta.data))] = paste0(M1M2_sig_rename,".sig")


VlnPlot(SCTall_fil2, ncol = 2, pt.size = 0, 
        split.by = "timepoint", split.plot = T, 
        group.by = "medium_clusters", 
        features = paste0(M1M2_sig_rename,".sig"))


VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="Macrophages"], ncol = 2, pt.size = 0, 
        group.by = "minor_clusters2", 
        features = paste0(M1M2_sig_rename,".sig"))

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="Macrophages"], ncol = 2, pt.size = 0, 
        group.by = "minor_clusters2", 
        features = "CX3CR1")
#


VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="Macrophages" & SCTall_fil2$timepoint == "pre"], ncol = 2, pt.size = 0, 
        group.by = "patient", same.y.lims = T, 
        features = paste0(M1M2_sig_rename,".sig"))


VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="Macrophages" & SCTall_fil2$timepoint == "post"], ncol = 2, pt.size = 0, 
        group.by = "patient", same.y.lims = T, 
        features = paste0(M1M2_sig_rename,".sig"))



SCTall_fil2@meta.data %>% 
  filter(medium_clusters=="Macrophages") %>% 
  group_by(patient, timepoint) %>% 
  filter(n() >= 50) %>% 
  summarise(m = median(M1.sig)) %>% 
  spread(key=timepoint, value=m) %>% 
  mutate(d = post-pre) %>% 
  arrange(d)


SCTall_fil2@meta.data %>% 
  filter(medium_clusters=="Macrophages") %>% 
  group_by(patient, timepoint) %>% 
  filter(n() >= 50) %>% 
  summarise(m = median(M2.sig)) %>% 
  spread(key=timepoint, value=m) %>% 
  mutate(d = post-pre) %>% 
  arrange(d)







SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, minor_clusters2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  filter(str_detect(minor_clusters2, "Mac|Den")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(p = post - pre) %>% 
  group_by(minor_clusters2) %>% 
  summarise(m = median(p)) %>% 
  arrange(-m) %>% 
  pull(minor_clusters2) %>% 
  as.character() ->
  myeloid_minorclusters_change_order


#########
##myeloid boxplot
###

library(ggallin)

SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  #filter(major_clusters %in% c("T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient.2,cohort, minor_clusters2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  filter(str_detect(minor_clusters2, "Mac")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(p = post - pre) %>% 
  group_by(minor_clusters2) %>% 
  mutate(p_value = t.test(unlist(p))$p.value,
         t_value = t.test(unlist(p))$statistic) %>% 
  #mutate(minor_clusters2 = factor(minor_clusters2, levels = myeloid_minorclusters_change_order)) %>% 
  ggplot() + aes(x = reorder(minor_clusters2,-t_value), y=p) + 
  geom_hline(yintercept = 0, linetype=2, color="grey") +
  geom_boxplot(outlier.color = NA) + 
  geom_point(size=1.25, aes(color=patient.2), position = position_jitter(width = 0.25)) + 
  ylab("Change in Immune Fraction") +
  scale_y_continuous(breaks = c(-0.6,-0.4,-0.2,0,0.2,0.4)) +
  #scale_y_continuous(trans = pseudolog10_trans) +
  scale_color_manual(values = patients.2_colors) +
  theme(panel.background = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        axis.title.x=element_blank(), 
        axis.title.y = element_text(hjust = 0.5, size=16),
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10)) + 
  guides(color = guide_legend(override.aes = list(size = 2, shape=15))) ->
  myeloid_change_boxplot
ggsave(myeloid_change_boxplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_02/myeloid_change_boxplot.pdf", width=8, height=4, device = cairo_pdf)





for(i in myeloid_minorclusters_change_order[-4]){
  print(i)
  SCTall_fil2@meta.data %>%
    filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
    group_by(sample, timepoint,patient,cohort, minor_clusters2) %>%
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(sample) %>% 
    mutate(p = n/sum(n)) %>%
    ungroup() %>% 
    mutate(minor_clusters2 = as.character(minor_clusters2)) %>% 
    filter(minor_clusters2 == i) %>% 
    #filter(str_detect(minor_clusters2, "Tcells - CD8:CCL4")) %>% 
    # ungroup() %>% 
    dplyr::select(-sample,-cohort,-n) %>% 
    spread(key = timepoint, value = p, fill = 0) -> temp
  print(t.test(temp$pre, temp$post, paired = T)$p.value)
  
}




