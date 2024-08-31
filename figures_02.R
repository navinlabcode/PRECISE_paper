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
NewSet1 = function(n){return(RColorBrewer::colorRampPalette(brewer.pal(9,"Set1")[c(1,5,6,3,2,4,8,7,9)])(n))}


patients.2_colors = patients.2_colors[sort(names(patients.2_colors))]

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


SCTall_fil2_tcells$minor_clusters5 = as.character(SCTall_fil2_tcells$minor_clusters4)
SCTall_fil2_tcells$minor_clusters5[SCTall_fil2_tcells$minor_clusters5 == "Tcells - CD4:Naive"] = "Tcells - CD4:CCR7"
SCTall_fil2_tcells$minor_clusters5[SCTall_fil2_tcells$minor_clusters5 == "Tcells - CD8:Naive"] = "Tcells - CD8:CCR7"
SCTall_fil2_tcells$minor_clusters5[SCTall_fil2_tcells$minor_clusters5 == "Tcells - CD4:TREG1"] = "Tcells - CD4:FOXP3"
SCTall_fil2_tcells$minor_clusters5[SCTall_fil2_tcells$minor_clusters5 == "Tcells - CD4:TREG2"] = "Tcells - CD4:LAG3"
SCTall_fil2_tcells$minor_clusters5[SCTall_fil2_tcells$minor_clusters5 == "Tcells - CD4:TEMRA"] = "Tcells - CD4:GNLY"
SCTall_fil2_tcells$minor_clusters5[SCTall_fil2_tcells$minor_clusters5 == "Tcells - CD4:TFH"] = "Tcells - CD4:GNLY"

  
  

sample_colors = c("#94E3D3", "#517D74", "#A882B7", "#3C2544", "#AA936F", "#685230", "#5167c9", "#26305E", "#FC6F60", "#872219", "#ED9D2D", "#774601", "#5E92A4", "#112630", "#a3a3a3", "#4a4a4a", "#d8db6e", "#5d5e30", "#B9FF94", "#5B7D48", "#FF9CB9", "#6B414E")
patient_colors = sample_colors[c(1,3,5,7,9,11,13,15,17,19,21)]

tcell_colors_df = data.frame(cluster = unique(SCTall_fil2_tcells$minor_clusters4), color = c("#2E99C7","#2E64C7","#FF0000","#FF9E9E","#005717","#757067","#00FF44","black","#BF8828","#6E4444","#487554","#C7C2B9","#A13030","#6E0000","#9CFFB6","#36BF5A","#523CAB","#7F3CAB"), stringsAsFactors = F)
#tcell_colors_df = tcell_colors_df[c(2,1,14,13,3,4,10,11,5,16,7,15,17,18,9,6,12,8),]
tcell_colors_df = tcell_colors_df[c(2,14,13,3,4,10,17,9,1,11,5,16,7,15,18,6,12,8),]
tcell_colors_df = left_join(tcell_colors_df, data.frame("cluster" = SCTall_fil2_tcells$minor_clusters4, "clusters2" = SCTall_fil2_tcells$minor_clusters5, row.names = NULL) %>% unique())
tcell_colors = tcell_colors_df$color
names(tcell_colors) = tcell_colors_df$clusters2


cytokines_df = read_csv("~/code/resources/gene_lists/cytokines.csv")
cytokines = cytokines_df$gene[cytokines_df$gene %in% row.names(SCTall_fil2)]
cytokines = c(cytokines, "GZMA","GZMB","GZMH","GZMK","GZMM","NKG7","GNLY","PRF1","TNF","IFNG")
cytokines = cytokines[cytokines %in% row.names(SCTall_fil2)]

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

revlog_trans <- function(base = exp(1)) {
  trans <- function(x) -log1p(x)
  inv <- function(x) expm1(-x)
  scales::trans_new("revlog1p", trans, inv, domain = c(0, Inf))
}






SCTall_fil2_tcells@meta.data %>%
  group_by(sample, timepoint,patient,cohort, minor_clusters4) %>%
  summarise(n = n()) %>% 
  #group_by(sample, timepoint,patient,cohort) %>% 
  ungroup() %>% 
  tidyr::complete(minor_clusters4,nesting(sample, timepoint,patient,cohort), fill=list(n=0)) %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  filter(!str_detect(minor_clusters4, "Cycling")) %>% 
  filter(!str_detect(minor_clusters4, "NK")) %>% 
  ungroup() %>% 
  select(-sample,-n,-cohort) %>% 
  spread(key = minor_clusters4, value = p, fill=0) %>% 
  gather(-timepoint, -patient, key = minor_clusters4, value =p) %>% 
  ggplot() + aes(x = timepoint, y=p, group=patient, color=patient) + 
  geom_point() + geom_line() + scale_color_manual(values = patient_colors) +
  facet_wrap(~minor_clusters4,  nrow = 4) + 
  theme(panel.background = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title=element_blank(), legend.title = element_blank(),
        strip.background = element_blank())





#######################
####tcells lineplot####
#######################

SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, minor_clusters2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  filter(str_detect(minor_clusters2, "Tcells")|str_detect(minor_clusters2, "^NK")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = minor_clusters2, value = p, fill=0) %>% 
  gather(-timepoint, -patient, key = minor_clusters2, value =p) %>% 
  mutate(minor_clusters2 = factor(minor_clusters2, levels = names(tcell_colors))) %>% 
  ggplot() + aes(x = timepoint, y=p, group=patient, color=patient) + 
  geom_point(size=2.5) + geom_line(size=0.8) + scale_color_manual(values = patient_colors) +
  facet_wrap(~minor_clusters2,  nrow = 3) + 
  theme(panel.background = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title=element_blank(), legend.title = element_blank(),
        strip.background = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text = element_text(size=12), legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        panel.spacing.y = unit(2, "lines")) + 
  guides(color = guide_legend(override.aes = list(size = 2, shape=15))) ->
  tcells_lineplot_grid

y.grob <- textGrob("Fraction of Immune Cells", 
                   gp=gpar(col="black", fontsize=20), rot=90)

grid.arrange(arrangeGrob(tcells_lineplot_grid, left = y.grob, padding = unit(1,"cm"))) -> tcells_lineplot
ggsave(tcells_lineplot, filename = "./plots/fig_pieces_01/tcells_lineplot.pdf", width=13, height=7, device = cairo_pdf)





  
temp<-plot_grid(tcells_lineplot)
ggdraw(add_sub(plot = tcells_lineplot, label = "Fraction of Immune Cells", angle = 90, vpadding=grid::unit(0,"lines"),y=0.5, x=0, vjust=0))


for(i in names(tcell_colors)){
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
    dplyr::select(-sample,-cohort,-n) %>% 
    spread(key = timepoint, value = p, fill = 0) -> temp 
  t.test(temp$pre, temp$post, paired = T)
  
  
  
Idents(SCTall_fil2_tcells) = "minor_clusters4"
temp = AverageExpression(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD8")], assays = c("RNA","SCT"), features = tfs[tfs%in%row.names(SCTall_fil2_tcells)], return.seurat = T)







SCTall_fil2_tcells_minor_clusters4_markers_CD8 = FindAllMarkers(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD8")], only.pos = T, assay = "SCT", logfc.threshold = 0.1)
SCTall_fil2_tcells_minor_clusters4_markers_CD4 = FindAllMarkers(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD4")], only.pos = T, assay = "SCT", logfc.threshold = 0.1)

SCTall_fil2_tcells_minor_clusters4_markers_CD8 %>% 
  filter(avg_logFC >= 0.5) %>% 
  pull(gene) %>% 
  unique() ->
  SCTall_fil2_tcells_minor_clusters4_markers_CD8_high

SCTall_fil2_tcells_minor_clusters4_markers_CD4 %>% 
  filter(avg_logFC >= 0.5) %>% 
  pull(gene) %>% 
  unique() ->
  SCTall_fil2_tcells_minor_clusters4_markers_CD4_high


SCTall_fil2_tcells_minor_clusters4_markers %>% 
  filter(cluster %in% c("NK Cells", "Tcells - Cycling","Tcells - GD")) %>% 
  filter(avg_logFC >= 0.5) %>% 
  pull(gene) %>% 
  unique() ->
  SCTall_fil2_tcells_minor_clusters4_markers_high




################
###tcell umap###
################

SCTall_fil2_tcells = RunUMAP(SCTall_fil2_tcells, n.neighbors = 30,
                             features = c(SCTall_fil2_tcells_minor_clusters4_markers_CD8_high,
                                          SCTall_fil2_tcells_minor_clusters4_markers_CD4_high,
                                          SCTall_fil2_tcells_minor_clusters4_markers_high), 
                             reduction.name = "UMAP_top", reduction.key = "umap_")
DimPlot(SCTall_fil2_tcells, reduction = "UMAP_top", cols = tcell_colors, pt.size = 1)
DimPlot(SCTall_fil2_tcells, reduction = "umap", cols = tcell_colors, pt.size = 0.5)

SCTall_fil2_tcells = RunUMAP(SCTall_fil2_tcells, n.neighbors = 50,
                             dims = 1:30, 
                             reduction.name = "UMAP_all", reduction.key = "umap_all_")
DimPlot(SCTall_fil2_tcells, reduction = "UMAP_all", cols = tcell_colors, pt.size = 0.5)


tcell_umap_df = data.frame(umap1 = SCTall_fil2_tcells@reductions$UMAP_top@cell.embeddings[,1], umap2 = SCTall_fil2_tcells@reductions$UMAP_top@cell.embeddings[,2], cluster = SCTall_fil2_tcells$minor_clusters4, timepoint = SCTall_fil2_tcells$timepoint, stringsAsFactors = F)
tcell_umap_df = left_join(tcell_umap_df, tcell_colors_df)
tcell_umap_df$color = as.character(tcell_umap_df$color)



tcell_umap_df %>% 
  ggplot() + aes(x=umap1, y=umap2, color=color)  +
  geom_point(aes(alpha=timepoint), size=0.005) + 
  theme_void() + theme(axis.title = element_text(size=16), axis.title.y = element_text(angle = 90)) +
  scale_color_identity("cluster", labels = names(tcell_colors), breaks = tcell_colors,
                    guide = "legend") + scale_alpha_manual(values = c(0.05,1)) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 20, 10)) ->
  plot_tcell_umap_post


tcell_umap_df %>% 
  ggplot() + aes(x=umap1, y=umap2, color=color)  +
  geom_point(aes(alpha=timepoint), size=0.01) + 
  theme_void() + theme(axis.title = element_text(size=16), axis.title.y = element_text(angle = 90)) + scale_alpha_discrete(guide = FALSE) +
  scale_color_identity("cluster", labels = names(tcell_colors), breaks = tcell_colors, guide = "legend") + 
  guides(color = guide_legend(override.aes = list(size=5))) + theme(legend.title = element_blank()) ->
  plot_tcell_umap_legend
plot_tcell_umap_legend = get_legend(plot_tcell_umap_legend)
plot_grid(plot_tcell_umap_legend) -> plot_tcell_umap_legend_plot
ggsave(plot = plot_tcell_umap_legend_plot, filename = "./plots/fig_pieces_01/tcell_umap_legend.pdf", width = 2, height = 4.5, dpi = 300)



tcell_umap_df %>% 
  mutate(timepoint = factor(timepoint, levels = c("post","pre"))) %>% 
  ggplot() + aes(x=umap1, y=umap2, color=color)  +
  geom_point(aes(alpha=timepoint), size=0.005) + scale_alpha_manual(values = c(0.05,1)) +
  theme_void() + theme(axis.title = element_text(size=16), axis.title.y = element_text(angle = 90)) +
  scale_color_identity("cluster", labels = names(tcell_colors), breaks = tcell_colors,
                       guide = "legend") + theme(legend.position = "none", plot.margin = margin(0, 0, 20, 10)) ->
  plot_tcell_umap_pre

plot_grid(plot_tcell_umap_pre,plot_tcell_umap_post) -> plot_tcell_umap
ggsave(plot = plot_tcell_umap, filename = "./plots/fig_pieces_01/tcell_umap.png", width = 6, height = 2.5, dpi = 300)



#####################
###tcell umap all####
#####################


tcell_umap_df = data.frame(umap1 = SCTall_fil2_tcells@reductions$umap@cell.embeddings[,1], umap2 = SCTall_fil2_tcells@reductions$umap@cell.embeddings[,2], cluster = SCTall_fil2_tcells$minor_clusters4, timepoint = SCTall_fil2_tcells$timepoint, stringsAsFactors = F)
tcell_umap_df = left_join(tcell_umap_df, tcell_colors_df)
tcell_umap_df$color = as.character(tcell_umap_df$color)


tcell_umap_df %>% 
  ggplot() + aes(x=umap1, y=umap2, color=color)  +
  geom_point(size=0.05) + 
  theme_void() + theme(axis.title = element_text(size=16), axis.title.y = element_text(angle = 90)) +
  scale_color_identity("cluster", labels = names(tcell_colors), breaks = tcell_colors,guide = "legend")  +
  theme(legend.position = "none", plot.margin = margin(0, 0, 20, 10)) ->
  plot_tcell_umap_all



#########################
###tcell umap_barchart###
#########################

SCTall_fil2_tcells@meta.data %>% 
  ggplot() + aes(x=timepoint, fill=minor_clusters4) + geom_bar(position="fill") + 
  scale_fill_manual(values = tcell_colors) +   theme_minimal() +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position = "none") + theme(panel.grid = element_blank()) + 
  theme(axis.title = element_blank())  + theme(axis.text.y=element_blank()) +
  theme(axis.text.x = element_text(size=12)) -> tcell_umap_bar

ggsave(plot = tcell_umap_bar, filename = "./plots/fig_pieces_01/tcell_umap_bar.pdf", width = 1.5, height = 4)






####################
###tcell heatmap####
####################


rbind(SCTall_fil2_tcells_minor_clusters4_markers,SCTall_fil2_tcells_minor_clusters4_markers_CD4,SCTall_fil2_tcells_minor_clusters4_markers_CD8) %>% 
  group_by(cluster) %>% 
  top_n(5, avg_logFC) %>% 
  filter(!gene %in% c(tfs,cytokines,checkpoint_genes)) %>% 
  filter(!gene %in% greceptors) %>% 
  filter(!gene %in% lectins) %>%  
  filter(!gene %in% CDs) %>% 
  filter(!gene %in% ISGs) %>% 
  View()



top_gene_filter_tcells = function(genes_list,n=5){
  
  SCTall_fil2_tcells_minor_clusters4_markers %>% 
    mutate(cluster = factor(cluster, levels = as.character(tcell_colors_df$cluster))) %>% 
    filter(cluster %in% c("NK Cells", "Tcells - Cycling","Tcells - GD")) %>% 
    filter(gene %in% genes_list) %>% 
    group_by(cluster) %>% 
    top_n(n,avg_logFC) %>% 
    pull(gene) %>% 
    unique() ->
    return_genes
  
  SCTall_fil2_tcells_minor_clusters4_markers_CD4 %>% 
    mutate(cluster = factor(cluster, levels = as.character(tcell_colors_df$cluster))) %>% 
    arrange(cluster) %>% 
    filter(gene %in% genes_list) %>% 
    group_by(cluster) %>% 
    top_n(n,avg_logFC) %>% 
    pull(gene) %>% 
    unique() %>% 
    c(return_genes) %>% 
    unique() ->
    return_genes
  
  SCTall_fil2_tcells_minor_clusters4_markers_CD8 %>% 
    mutate(cluster = factor(cluster, levels = as.character(tcell_colors_df$cluster))) %>% 
    arrange(cluster) %>% 
    filter(gene %in% genes_list) %>% 
    group_by(cluster) %>% 
    top_n(n,avg_logFC) %>% 
    pull(gene) %>% 
    c(return_genes) %>% 
    unique() ->
    return_genes
  
  return(return_genes)
}


tcell_genes_tfs = top_gene_filter_tcells(tfs, 5)
tcell_genes_cyto = top_gene_filter_tcells(cytokines, 5)
tcell_genes_ckpts = top_gene_filter_tcells(checkpoint_tcells, 5)
recepts = c(greceptors,lectins,CDs)
tcell_genes_recepts = top_gene_filter_tcells(recepts, 5)
tcell_genes_ISGs = top_gene_filter_tcells(ISGs, 5)

tcell_genes_allgroups = c(tcell_genes_tfs,tcell_genes_cyto,tcell_genes_ckpts,figures_01,tcell_genes_ISGs)
tcell_genes_allgroups = unique(tcell_genes_allgroups)



#tcell_genes_tfs = rev(tcell_genes_tfs)

SCTall_fil2_tcells$minor_clusters4 = factor(SCTall_fil2_tcells$minor_clusters4, levels = as.character(tcell_colors_df$cluster))
Idents(SCTall_fil2_tcells) = "minor_clusters4"

SCTall_fil2_tcells_avg_CD8 = AverageExpression(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD8")], features = tcell_genes_allgroups, return.seurat = T)
DoHeatmap(SCTall_fil2_tcells_avg_CD8, features = tcell_genes_allgroups, assay = "SCT", slot = "scale.data", label = FALSE,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=6)) -> tcell_heatmap_CD8


library(gtable)
g <- ggplotGrob(tcell_heatmap_CD8)
s <- gtable_filter(g, 'axis-l|ylab', trim=F)  # use trim depending on need
grid.newpage()
grid.draw(s)


DoHeatmap(SCTall_fil2_tcells_avg_CD8, features = tcell_genes_allgroups, assay = "SCT", slot = "scale.data", label = FALSE,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.text.y = element_blank()) -> tcell_heatmap_CD8




SCTall_fil2_tcells_avg_CD4 = AverageExpression(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD4")], features = tcell_genes_allgroups, return.seurat = T)
DoHeatmap(SCTall_fil2_tcells_avg_CD4, features = tcell_genes_allgroups, assay = "SCT", label = FALSE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  theme(axis.text.y = element_blank()) -> tcell_heatmap_CD4


SCTall_fil2_tcells_avg_all = AverageExpression(SCTall_fil2_tcells, features = tcell_genes_allgroups, return.seurat = T)
SCTall_fil2_tcells_avg_other = SCTall_fil2_tcells_avg_all[,!str_detect(colnames(SCTall_fil2_tcells_avg_all),"CD[48]")]
DoHeatmap(SCTall_fil2_tcells_avg_other, features = tcell_genes_allgroups, assay = "SCT",label = FALSE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  theme(axis.text.y = element_blank()) -> tcell_heatmap_other


plot_grid(s,tcell_heatmap_CD8, tcell_heatmap_CD4, tcell_heatmap_other, nrow=1, align = "h", rel_widths = c(4,8,7,3)) -> tcells_heatmap_bottom
ggsave(tcells_heatmap_bottom, filename = "./plots/fig_pieces_01/tcells_heatmap_bottom.pdf", width=2, height=9)
#ggsave(tcells_heatmap_bottom, filename = "./plots/fig_pieces_01/tcells_heatmap_bottom.png", width=2, height=9)





tcell_genes_lineage = c("CD3E","CD3D","CD3G","CD8A","CD4","TRAC","TRBC2","TRDC","NCAM1")

SCTall_fil2_tcells_avg_all_lin = AverageExpression(SCTall_fil2_tcells, features = tcell_genes_lineage, return.seurat = T)
DoHeatmap(SCTall_fil2_tcells_avg_all_lin, features = tcell_genes_lineage, assay = "SCT", label = FALSE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=6)) -> tcell_heatmap_markers

ggsave(tcell_heatmap_markers, filename = "./plots/fig_pieces_01/tcell_heatmap_top.pdf", width=2, height=1.1)


DoHeatmap(SCTall_fil2_tcells_avg_all_lin, features = tcell_genes_lineage, assay = "SCT", label = FALSE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=6)) -> temp
tcell_heatmap_legend = get_legend(temp)
plot_grid(tcell_heatmap_legend) -> tcell_heatmap_legend 
ggsave(tcell_heatmap_legend, filename = "./plots/fig_pieces_01/tcell_heatmap_legend.pdf", width=2, height=9)

##################################
###tcell heatmap ungrouped########
##################################

SCTall_fil2_tcells_minor_clusters4_markers_CD4 %>% 
  filter(!str_detect(gene, "MT-")) %>% 
  filter(!str_detect(gene, "RP[LS]")) %>% 
  filter(!str_detect(gene, "^HSP")) %>% 
  filter(!gene %in% c("FOS","FOSB","JUN")) %>% 
  filter(!gene %in% tcell_genes_lineage) %>% 
  mutate(pct = pct.1-pct.2) %>% 
  filter(pct > 0.12) %>% 
  filter(pct.1 > 0.25) %>%
  filter(pct.2 < 0.5) %>%
  group_by(gene) %>% 
  top_n(1, pct.1) %>% 
  #group_by(gene) %>% 
  #top_n(2, pct.1) %>% 
  #arrange(cluster, -pct) %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC) %>% 
  pull(gene) %>% 
  unique() ->
  tcell_genes_CD4
  


SCTall_fil2_tcells_minor_clusters4_markers_CD8 %>% 
  filter(!str_detect(gene, "MT-")) %>% 
  filter(!str_detect(gene, "RP[LS]")) %>% 
  filter(!str_detect(gene, "^HSP")) %>% 
  filter(!str_detect(gene, "^HB")) %>% 
  filter(!str_detect(gene, "^TR[ABDG][VDJ]")) %>% 
  filter(!gene %in% c("FOS","FOSB","JUN")) %>% 
  filter(!gene %in% tcell_genes_lineage) %>% 
  mutate(pct = pct.1-pct.2) %>% 
  #filter(pct > 0.12) %>% 
  #filter(pct.1 > 0.25) %>%
  group_by(gene) %>% 
  top_n(1, pct.1) %>% 
  arrange(cluster, -avg_logFC) %>%  
  group_by(cluster) %>% 
  top_n(10, avg_logFC) %>% 
  pull(gene) %>% 
  unique() ->
  tcell_genes_CD8


SCTall_fil2_tcells_minor_clusters4_markers %>% 
  mutate(cluster = factor(cluster, levels = tcell_colors_df$cluster)) %>% 
  filter(cluster %in% c("NK Cells", "Tcells - Cycling","Tcells - GD")) %>% 
  filter(!str_detect(gene, "MT-")) %>% 
  filter(!str_detect(gene, "RP[LS]")) %>% 
  filter(!str_detect(gene, "^HSP")) %>% 
  filter(!gene %in% c("FOS","FOSB","JUN")) %>% 
  filter(!gene %in% tcell_genes_lineage) %>% 
  mutate(pct = pct.1-pct.2) %>% 
  filter(pct > 0.12) %>% 
  filter(pct.1 > 0.25) %>%
  filter(pct.2 < 0.5) %>%
  group_by(gene) %>% 
  top_n(1, pct.1) %>% 
  #group_by(gene) %>% 
  #top_n(2, pct.1) %>% 
  #arrange(cluster, -pct) %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC) %>% 
  pull(gene) %>% 
  unique() ->
  tcell_genes_other

tcell_genes_all= c(tcell_genes_CD8, tcell_genes_CD4,tcell_genes_other)

SCTall_fil2_tcells$minor_clusters4 = factor(SCTall_fil2_tcells$minor_clusters4, levels = as.character(tcell_colors_df$cluster))
Idents(SCTall_fil2_tcells) = "minor_clusters4"

#SCTall_fil2_tcells_avg_CD8 = AverageExpression(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[48]")], features = tcell_genes_CD8, return.seurat = T)
SCTall_fil2_tcells_avg_CD8 = AverageExpression(SCTall_fil2_tcells, features = tcell_genes_CD8, return.seurat = T)

DoHeatmap(SCTall_fil2_tcells_avg_CD8, features = tcell_genes_CD8, assay = "SCT", slot = "scale.data", label = FALSE,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=6)) -> tcell_heatmap_CD8


#SCTall_fil2_tcells_avg_CD4 = AverageExpression(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD4")], features = tcell_genes_all, return.seurat = T)
SCTall_fil2_tcells_avg_CD4 = AverageExpression(SCTall_fil2_tcells, features = tcell_genes_CD4, return.seurat = T)

DoHeatmap(SCTall_fil2_tcells_avg_CD4, features = tcell_genes_CD4, assay = "SCT", label = FALSE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  theme(axis.text.y = element_text(size=6)) -> tcell_heatmap_CD4


SCTall_fil2_tcells_avg_all = AverageExpression(SCTall_fil2_tcells, features = tcell_genes_other, return.seurat = T)

DoHeatmap(SCTall_fil2_tcells_avg_all, features = tcell_genes_other, assay = "SCT", label = FALSE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  theme(axis.text.y = element_text(size=6)) -> tcell_heatmap_other

tcell_genes_lineage = c("CD3E","CD3D","CD3G","CD8A","CD4","TRAC","TRBC2","TRDC","NCAM1")

SCTall_fil2_tcells_avg_all_lin = AverageExpression(SCTall_fil2_tcells, features = tcell_genes_lineage, return.seurat = T)
DoHeatmap(SCTall_fil2_tcells_avg_all_lin, features = tcell_genes_lineage, assay = "SCT", slot = "scale.data", label = FALSE,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.text.y = element_text(size=6)) -> tcell_heatmap_markers



library(gtable)
g <- ggplotGrob(tcell_heatmap_CD8)
s1 <- gtable_filter(g, 'axis-l|ylab', trim=F)  # use trim depending on need
grid.newpage()
grid.draw(s1)

library(gtable)
g <- ggplotGrob(tcell_heatmap_CD4)
s2 <- gtable_filter(g, 'axis-l|ylab', trim=F)  # use trim depending on need
grid.newpage()
grid.draw(s2)

library(gtable)
g <- ggplotGrob(tcell_heatmap_other)
s3 <- gtable_filter(g, 'axis-l|ylab', trim=F)  # use trim depending on need
grid.newpage()
grid.draw(s3)

g <- ggplotGrob(tcell_heatmap_markers)
s4 <- gtable_filter(g, 'axis-l|ylab', trim=F)  # use trim depending on need
grid.newpage()
grid.draw(s4)


DoHeatmap(SCTall_fil2_tcells_avg_CD8, features = tcell_genes_CD8, assay = "SCT", slot = "scale.data", label = FALSE,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.text.y = element_blank()) -> tcell_heatmap_CD8

DoHeatmap(SCTall_fil2_tcells_avg_CD4, features = tcell_genes_CD4, assay = "SCT", slot = "scale.data", label = FALSE,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.text.y = element_blank()) -> tcell_heatmap_CD4

DoHeatmap(SCTall_fil2_tcells_avg_all, features = tcell_genes_other, assay = "SCT", slot = "scale.data", label = FALSE,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.text.y = element_blank()) -> tcell_heatmap_other

DoHeatmap(SCTall_fil2_tcells_avg_all_lin, features = tcell_genes_lineage, assay = "SCT", slot = "scale.data", label = FALSE,
          disp.min = -2, disp.max = 2, draw.lines = F, group.colors = tcell_colors, raster = F) + NoLegend() +
  scale_fill_gradient2(low = "#0054b5",mid = "#dbdbdb", high = "#ebac00", na.value = "00000000") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.text.y = element_blank())  -> tcell_heatmap_markers

plot_grid(s4, tcell_heatmap_markers, s1,tcell_heatmap_CD8, s2, tcell_heatmap_CD4, s3, tcell_heatmap_other, nrow=4, align = "h", rel_widths = c(0.1,1), rel_heights = c(0.15,1,1,0.43)) -> tcells_heatmap_bottom

ggsave(tcells_heatmap_bottom, filename = "./plots/fig_pieces_01/tcells_heatmap_ungroup.pdf", device = cairo_pdf, width=4.5, height=24)


  
















c("CD3E","CD3D","CD3G","CD8A","CD4","TRAC","TRBC1","TRDC","TRGC1","B3GAT1",checkpoint_tcells, "IFNG","TNF","CCR7","TCF7","LEF1", SCTall_fil2_t_CD4_markers_top,SCTall_fil2_t_CD8_markers_top,"MKI67","TUBB","ACTB")


rbind(SCTall_fil2_t_CD4_markers,SCTall_fil2_t_CD8_markers) %>% 
  filter(gene %in% cytokines) %>% View()
  group_by(cluster) %>% 
  top_n(5,avg_logFC) %>% 
  pull(gene) %>% 
  unique() ->
  tcell_genes_cyto


  
pal = colorRampPalette(c("#0079c9", "#cedbcc","#dbb700"))
scales::show_col(pal(25))

##########################
##########################
###expresison channges####
##########################
##########################

Idents(SCTall_fil2_tcells) = "timepoint"
SCTall_fil2_tcells_timepoint_SCT_markers = FindAllMarkers(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[48]")], logfc.threshold = 0, assay = "SCT", min.pct = 0.01)
SCTall_fil2_tcells_timepoint_RNA_markers = FindAllMarkers(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[48]")], logfc.threshold = 0, assay = "RNA", min.pct = 0.01)

SCTall_fil2_tcells_timepoint_SCT_markers %>% 
  filter(!str_detect(gene, "^HB[AB]")) %>% 
  filter(cluster == "post") ->
  SCTall_fil2_tcells_timepoint_SCT_markers_fil

SCTall_fil2_tcells_timepoint_RNA_markers %>% 
  filter(!is.infinite(avg_logFC)) %>% 
  filter(!str_detect(gene, "^HB[AB]")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(cluster == "post") ->
  SCTall_fil2_tcells_timepoint_RNA_markers_fil

#SCTall_fil2_tcells_patient_timepoint_SCT_markers = c()
#SCTall_fil2_tcells_patient_timepoint_RNA_markers = c()
#for(i in patients[c(2,4,5,6,7,9,10,11)]){
#  temp = FindAllMarkers(SCTall_fil2_tcells[,(str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[48]") & SCTall_fil2_tcells$patient==i)], logfc.threshold = 0, assay = "SCT", only.pos = TRUE, min.pct = 0.01)
#  temp$patient = i
#  SCTall_fil2_tcells_patient_timepoint_SCT_markers = rbind(temp, SCTall_fil2_tcells_patient_timepoint_SCT_markers)
#  temp = FindAllMarkers(SCTall_fil2_tcells[,(str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[48]") & SCTall_fil2_tcells$patient==i)], logfc.threshold = 0, assay = "RNA", only.pos = TRUE, min.pct = 0.01)
#  temp$patient = i
#  SCTall_fil2_tcells_patient_timepoint_RNA_markers = rbind(temp, SCTall_fil2_tcells_patient_timepoint_RNA_markers)
#}





####################
###cd4cd8 prepost###
####################


Idents(SCTall_fil2_tcells) = "timepoint"
SCTall_fil2_tcells_timepoint_SCT_markers_CD8 = FindAllMarkers(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[8]")], logfc.threshold = 0, assay = "SCT", min.pct = 0.01)
SCTall_fil2_tcells_timepoint_SCT_markers_CD4 = FindAllMarkers(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[4]")], logfc.threshold = 0, assay = "SCT", min.pct = 0.01)


SCTall_fil2_tcells_timepoint_SCT_markers_CD8 %>% 
  filter(!str_detect(gene, "^HB[AB]")) %>% 
  filter(cluster == "post") %>% 
  arrange(-avg_logFC) ->
  SCTall_fil2_tcells_timepoint_SCT_markers_CD8_fil

SCTall_fil2_tcells_timepoint_SCT_markers_CD4 %>% 
  filter(!str_detect(gene, "^HB[AB]")) %>% 
  filter(cluster == "post") %>% 
  arrange(-avg_logFC) ->
  SCTall_fil2_tcells_timepoint_SCT_markers_CD4_fil


SCTall_fil2_tcells_timepoint_SCT_CD8_gsea_C7 = rungsea(SCTall_fil2_tcells_timepoint_SCT_markers_CD8_fil, cats = "C7")
SCTall_fil2_tcells_timepoint_SCT_CD4_gsea_C7 = rungsea(SCTall_fil2_tcells_timepoint_SCT_markers_CD4_fil, cats = "C7")

SCTall_fil2_tcells_timepoint_SCT_CD8_gsea_C7 %>% 
  filter((str_detect(pathway, "TCELL")|str_detect(pathway, "TREG")|str_detect(pathway, "EFF")|str_detect(pathway, "CONV")|str_detect(pathway, "MEM")|str_detect(pathway, "CD[48]"))) %>%  
  filter(!(str_detect(pathway, "DC"))) %>%  
  filter(!(str_detect(pathway, "BCELL") & (!str_detect(pathway,"TCELL")))) %>% 
  arrange(-NES) ->
  SCTall_fil2_tcells_timepoint_SCT_CD8_gsea_C7_fil


SCTall_fil2_tcells_timepoint_SCT_CD4_gsea_C7 %>% 
  filter((str_detect(pathway, "TCELL")|str_detect(pathway, "TREG")|str_detect(pathway, "EFF")|str_detect(pathway, "CONV")|str_detect(pathway, "MEM")|str_detect(pathway, "CD[48]"))) %>%  
  filter(!(str_detect(pathway, "DC"))) %>%  
  filter(!(str_detect(pathway, "BCELL") & (!str_detect(pathway,"TCELL")))) %>% 
  arrange(-NES) ->
  SCTall_fil2_tcells_timepoint_SCT_CD4_gsea_C7_fil




######################
###cluster markers####
######################


Idents(SCTall_fil2_tcells) = "minor_clusters4"
SCTall_fil2_tcells_SCT_CD4_naive_markers = FindMarkers(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[4]")], ident.1 = "Tcells - CD4:Naive", logfc.threshold = 0, assay = "SCT", min.pct = 0.01)
SCTall_fil2_tcells_SCT_CD4_naive_markers$gene = row.names(SCTall_fil2_tcells_SCT_CD4_naive_markers)
SCTall_fil2_tcells_SCT_CD4_naive_markers$cluster = "Tcells - CD4:Naive"

SCTall_fil2_tcells_SCT_CD8_CCL4_markers = FindMarkers(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters4,"CD[8]")], ident.1 = "Tcells - CD8:CCL4", logfc.threshold = 0, assay = "SCT", min.pct = 0.01)
SCTall_fil2_tcells_SCT_CD8_CCL4_markers$gene = row.names(SCTall_fil2_tcells_SCT_CD8_CCL4_markers)
SCTall_fil2_tcells_SCT_CD8_CCL4_markers$cluster = "Tcells - CD8:CCL4"


SCTall_fil2_tcells_SCT_CD8_CCL4_gsea_C7 = rungsea(SCTall_fil2_tcells_SCT_CD8_CCL4_markers, cats = "C7", mincat = 5)
SCTall_fil2_tcells_SCT_CD8_CCL4_gsea_C7 %>% 
  filter((str_detect(pathway, "TCELL")|str_detect(pathway, "TREG")|str_detect(pathway, "EFF")|str_detect(pathway, "CONV")|str_detect(pathway, "MEM")|str_detect(pathway, "CD[48]"))) %>%  
  filter(!(str_detect(pathway, "DC"))) %>%  
  filter(!(str_detect(pathway, "BCELL") & (!str_detect(pathway,"TCELL")))) %>% 
  arrange(-NES) %>% 
  #filter(size <= 10) %>% 
  #filter(!str_detect(pathway, "CD4")) %>% 
  #filter(!str_detect(pathway, "TREG")) %>% 
  View()



SCTall_fil2_tcells_SCT_CD4_naive_gsea_C7 = rungsea(SCTall_fil2_tcells_SCT_CD4_naive_markers, cats = "C7", mincat = 5)
SCTall_fil2_tcells_SCT_CD4_naive_gsea_C7 %>% 
  filter((str_detect(pathway, "TCELL")|str_detect(pathway, "TREG")|str_detect(pathway, "EFF")|str_detect(pathway, "CONV")|str_detect(pathway, "MEM")|str_detect(pathway, "CD[48]"))) %>%  
  filter(!(str_detect(pathway, "DC"))) %>%  
  filter(!(str_detect(pathway, "BCELL") & (!str_detect(pathway,"TCELL")))) %>% 
  arrange(-NES) %>% 
  #filter(size <= 10) %>% 
  #filter(!str_detect(pathway, "CD4")) %>% 
  #filter(!str_detect(pathway, "TREG")) %>% 
  View()






                                                
#################          
###volcanoplot###
#################

anno_genes = c("CXCR4","NFKBIA","JUNB","CD69","TSC22D3","FOS","KLF2","KLF6","NFKBIA","FGFBP2","CCL4","CCR7","GZMB","GNLY","GZMH")
anno_genes = c("CXCR4","CD69","FOS","KLF2","FGFBP2","CCL4","CCR7","GZMB","GNLY","GZMH")

SCTall_fil2_tcells_timepoint_SCT_markers_fil = readRDS("./Rda/SCTall_fil2_tcells_timepoint_SCT_markers_fil.rds")
SCTall_fil2_tcells_timepoint_SCT_markers_fil %>% 
  mutate(p_val_adj = if_else(p_val_adj==0, 10^-300, p_val_adj)) %>% 
  mutate(anno = if_else(gene %in% anno_genes, "yes", "no")) -> temp
  temp %>% 
    mutate(dir = "no_change") %>% 
    mutate(dir = if_else(avg_logFC >= 0.25 & p_val_adj <= 0.01, "up", dir)) %>% 
    mutate(dir = if_else(avg_logFC <= -0.25 & p_val_adj <= 0.01, "down", dir)) %>% 
  ggplot(aes(x = avg_logFC, y = -log10(p_val_adj))) +
  theme_minimal() +
  scale_x_continuous(limits = c(-1.3,1.3),expand = c(0,.1)) +
  scale_y_continuous(limits = c(0,300),expand = c(0,5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 2, linetype = "dashed") +
  #geom_hline(yintercept = 0) +
  #geom_hline(yintercept = 25, linetype = "dashed") +
  #geom_vline(xintercept = c(-.25, .25), linetype = "dashed") +
  theme_linedraw() + 
  geom_point(aes(color = dir), alpha = 1) +
  scale_color_manual(values = c("royalblue4","grey50","goldenrod")) +
  theme(panel.grid = element_blank()) +
  xlab("Average Log Fold Change") +
  ylab("-log10(Adjusted P-Value)") +
  theme(legend.position = "none",axis.text = element_text(size=16), axis.title = element_text(size = 20)) +
  theme(axis.ticks = element_line(size=0.75)) +
  theme(panel.border = element_blank()) +
  theme(axis.line.x = element_line(size=0.75), axis.line.y = element_line(size=0.75)) +
  geom_text_repel(
    data = subset(temp, anno=="yes"),
    aes(label = gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    min.segment.length = 0.1) -> tcell_volcano_plot
  ggsave(tcell_volcano_plot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_01/tcells_volcano_plot.pdf", width=6, height=5, device = cairo_pdf)
  
  device = cairo_pdf



#################
###gsea##########
#################
  

#SCTall_fil2_tcells_timepoint_SCT_gsea_C2 = rungsea(SCTall_fil2_tcells_timepoint_SCT_markers_fil, cats = "C2")
#SCTall_fil2_tcells_timepoint_SCT_gsea_C5 = rungsea(SCTall_fil2_tcells_timepoint_SCT_markers_fil, cats = "C5")
SCTall_fil2_tcells_timepoint_SCT_gsea_C7 = rungsea(SCTall_fil2_tcells_timepoint_SCT_markers_fil, cats = "C7")
#SCTall_fil2_tcells_timepoint_RNA_gsea_C2 = rungsea(SCTall_fil2_tcells_timepoint_RNA_markers_fil, cats = "C2")
#SCTall_fil2_tcells_timepoint_RNA_gsea_C5 = rungsea(SCTall_fil2_tcells_timepoint_RNA_markers_fil, cats = "C5")
#SCTall_fil2_tcells_timepoint_RNA_gsea_C7 = rungsea(SCTall_fil2_tcells_timepoint_RNA_markers_fil, cats = "C7")
  

SCTall_fil2_tcells_timepoint_SCT_gsea_C7 %>% 
  filter((str_detect(pathway, "TCELL")|str_detect(pathway, "TREG")|str_detect(pathway, "EFF")|str_detect(pathway, "CONV")|str_detect(pathway, "MEM")|str_detect(pathway, "CD[48]"))) %>%  
  filter(!(str_detect(pathway, "DC"))) %>%  
  filter(!(str_detect(pathway, "BCELL") & (!str_detect(pathway,"TCELL")))) ->
  SCTall_fil2_tcells_timepoint_SCT_gsea_C7_fil


#SCTall_fil2_tcells_timepoint_RNA_gsea_C7 %>% 
#  filter((str_detect(pathway, "TCELL")|str_detect(pathway, "TREG")|str_detect(pathway, "EFF")|str_detect(pathway, "CONV")|str_detect(pathway, "MEM")|str_detect(pathway, "CD[48]"))) %>%  
#  filter(!(str_detect(pathway, "DC"))) %>% 
#  filter(!(str_detect(pathway, "BCELL") & (!str_detect(pathway,"TCELL")))) %>%  
#  View()


#SCTall_fil2_tcells_timepoint_SCT_gsea_C7_fil %>% 
#  mutate(avg_logFC = if_else(cluster=="pre", avg_logFC*-1, avg_logFC)) %>% 
#  mutate(cluster = "post") ->
#  SCTall_fil2_tcells_timepoint_SCT_gsea_C7_fil



SCTall_fil2_tcells_timepoint_SCT_gsea_C7_fil %>%
  filter(dense_rank(NES) <= 10 | dense_rank(desc(NES)) <= 10) %>% 
  mutate(dir = if_else(NES >=0, "pos","neg")) %>% 
  ggplot() + aes(x=reorder(pathway, NES), y = NES, fill=dir) + geom_bar(stat = "identity") +
  scale_fill_manual(values = c("steelblue","firebrick")) +
  theme_minimal() + theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(-3,3), expand = c(0,0)) +
  theme(axis.line.x = element_line(size=0.5)) +
  theme(axis.title.x = element_text(size=12), axis.text=element_text(size=8)) +
  theme(panel.grid = element_blank(), legend.position = "none") + coord_flip() ->
  tcells_gsea_plot
ggsave(tcells_gsea_plot, filename = "./plots/fig_pieces_01/tcells_gsea_plot.pdf", width=10, height=5)



SCTall_fil2_tcells_timepoint_SCT_CD4_gsea_C7_fil %>%
  filter(dense_rank(NES) <= 10 | dense_rank(desc(NES)) <= 10) %>% 
  mutate(dir = if_else(NES >=0, "pos","neg")) %>% 
  ggplot() + aes(x=reorder(pathway, NES), y = NES, fill=dir) + geom_bar(stat = "identity") +
  scale_fill_manual(values = c("steelblue","firebrick")) +
  theme_minimal() + theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(-3,3), expand = c(0,0)) +
  theme(axis.line.x = element_line(size=0.5)) +
  theme(axis.title.x = element_text(size=12), axis.text=element_text(size=8)) +
  theme(panel.grid = element_blank(), legend.position = "none") + coord_flip()


#SCTall_fil2_tcells_patient_timepoint_SCT_gsea_C2 = lapply(split(SCTall_fil2_tcells_patient_timepoint_SCT_markers_fix, SCTall_fil2_tcells_patient_timepoint_SCT_markers_fix$patient), function(x){rungsea(x, cats = "C2", fc_cutoff = 0.2)})
#SCTall_fil2_tcells_patient_timepoint_SCT_gsea_C5 = lapply(split(SCTall_fil2_tcells_patient_timepoint_SCT_markers_fix, SCTall_fil2_tcells_patient_timepoint_SCT_markers_fix$patient), function(x){rungsea(x, cats = "C5", fc_cutoff = 0.2)})
#SCTall_fil2_tcells_patient_timepoint_SCT_gsea_C7 = lapply(split(SCTall_fil2_tcells_patient_timepoint_SCT_markers_fix, SCTall_fil2_tcells_patient_timepoint_SCT_markers_fix$patient), function(x){rungsea(x, cats = "C7", fc_cutoff = 0.2)})


#bind_rows(SCTall_fil2_tcells_patient_timepoint_SCT_gsea_C7, .id = "patient") %>% 
#  filter(NES <= -2) %>% 
#  filter((str_detect(pathway, "TCELL")|str_detect(pathway, "TREG")|str_detect(pathway, "EFF")|str_detect(pathway, "CONV")|str_detect(pathway, "MEM")|str_detect(pathway, "CD[48]"))) %>%  
#  filter(!(str_detect(pathway, "DC"))) %>% 
#  filter(!(str_detect(pathway, "BCELL") & (!str_detect(pathway,"TCELL")))) %>%  
#  group_by(pathway, cluster) %>% 
#  summarise(n=n()) %>% 
#  View()

#SCTall_fil2_tcells_patient_timepoint_SCT_markers_fix %>% 
#  group_by(gene) %>% 
#  summarise(m = mean(avg_logFC), s = sd(avg_logFC)) %>% 
#  View()

#SCTall_fil2_tcells_timepoint_SCT_markers_fil %>% 
#  ggplot() + aes(x = avg_logFC, y = -log10(p_val_adj)) + geom_point()


#SCTall_fil2_tcells_timepoint_RNA_markers %>% 
#  ggplot() + aes(x = avg_logFC, y = -log10(p_val_adj)) + geom_point()



##############
####CD4 imt###
##############

#p = 0.00105

#IMT_TME_merged %>% 
#  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
#  filter(marker == "CD4") %>% 
#  ggplot() + aes(x=timepoint, y=density) + 
##  ylab(expression('CD4'^'+'~'/mm'^2)) +
#  geom_boxplot() +
#  geom_point(color="grey40", size=2, aes(group=patient)) + 
#  geom_line(color="grey40", size=0.5, aes(group=patient)) +
#  scale_y_continuous(limits = c(0,2000)) +
#  theme_one_line 
IMT_TME_merged %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  filter(marker == "CD4") %>% 
  ggpaired(x = "timepoint", y = "density",
           color = "timepoint", line.color = "black", point.size = 2,
           line.size = 0.5,palette = c("black","black")) +
  #stat_compare_means(paired = TRUE) + 
  ylab(expression('CD4'^'+'*'/mm'^2)) +
  theme_one_line + theme(legend.position = "none") ->
  CD4_IHC_plot

ggsave(CD4_IHC_plot, filename = ".fig/plots/fig_pieces_01/CD4_IHC_plot.pdf", width=2.5, height=3, device = cairo_pdf)


  

IMT_TME_merged %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  filter(marker == "CD4") %>% 
  tidyr::spread(timepoint, density) %>% 
  summarise(ttest = list(t.test(pre, post,paired = TRUE))) %>% 
  pull(ttest) ->
  temp



############################
####tcell lineplot, vert####
############################




SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, minor_clusters2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  filter(str_detect(minor_clusters2, "Tcells")|str_detect(minor_clusters2, "^NK")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(p = post - pre) %>% 
  group_by(minor_clusters2) %>% 
  summarise(m = median(p)) %>% 
  arrange(-m) %>% 
  pull(minor_clusters2) ->
  tcell_minorclusters_change_order


SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  #filter(major_clusters %in% c("T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, minor_clusters2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  filter(str_detect(minor_clusters2, "Tcells")|str_detect(minor_clusters2, "^NK")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(p = post - pre) %>% 
  mutate(minor_clusters2 = factor(minor_clusters2, levels = tcell_minorclusters_change_order)) %>% 
  ggplot() + aes(x = minor_clusters2, y=p) + 
  geom_hline(yintercept = 0, linetype=2, color="grey") +
  geom_boxplot(outlier.color = NA) + 
  geom_point(size=1.25, aes(color=patient), position = position_jitter(width = 0.25)) + 
  ylab("Change in Immune Fraction") +
  scale_color_manual(values = patient_colors) +
  theme(panel.background = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.x=element_blank(), 
        axis.title.y = element_text(hjust = 0.5, size=16),
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10)) + 
  guides(color = guide_legend(override.aes = list(size = 2, shape=15))) ->
  tcell_change_boxplot
ggsave(tcell_change_boxplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_01/tcell_change_boxplot.pdf", width=8.4, height=4, device = cairo_pdf)









SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  #filter(major_clusters %in% c("T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, minor_clusters2, group) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  filter(str_detect(minor_clusters2, "Tcells")|str_detect(minor_clusters2, "^NK")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(p = post - pre) %>% 
  mutate(minor_clusters2 = factor(minor_clusters2, levels = tcell_minorclusters_change_order)) %>% 
  ggplot() + aes(x =group, y=p,color=group) + 
  geom_hline(yintercept = 0, linetype=2, color="grey") +
  geom_boxplot(outlier.color = NA) + 
  geom_point(size=1.25, position = position_jitter(width = 0.25), aes(color=patient)) + 
  ylab("Change in Immune Fraction") +
  scale_color_manual(values = group_colors) +
  theme(panel.background = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.x=element_blank(), 
        axis.title.y = element_text(hjust = 0.5, size=16),
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10)) + 
  guides(color = guide_legend(override.aes = list(size = 2, shape=15))) +
  facet_wrap(~minor_clusters2)





SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  #filter(major_clusters %in% c("T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, minor_clusters2, group) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  filter(str_detect(minor_clusters2, "Tcells")|str_detect(minor_clusters2, "^NK")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  filter(timepoint == "post") %>% 
  mutate(minor_clusters2 = factor(minor_clusters2, levels = tcell_minorclusters_change_order)) %>% 
  complete(minor_clusters2,nesting(patient, group,timepoint), fill=list(p=0)) %>% 
  ggplot() + aes(x =group, y=p,color=group) + 
  geom_hline(yintercept = 0, linetype=2, color="grey") +
  geom_boxplot(outlier.color = NA) + 
  geom_point(size=1.25, position = position_jitter(width = 0.25), aes(color=patient)) + 
  ylab("Change in Immune Fraction") +
  scale_color_manual(values = group_colors) +
  theme(panel.background = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.x=element_blank(), 
        axis.title.y = element_text(hjust = 0.5, size=16),
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10)) + 
  guides(color = guide_legend(override.aes = list(size = 2, shape=15))) +
  facet_wrap(~minor_clusters2)


############################
####tcell lineplot, vert####
############################





SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("T - Cells")) %>% 
  filter(str_detect(minor_clusters2, "Tcells") | str_detect(minor_clusters2, "^NK")) %>% 
  mutate(CD4CD8 = "other") %>% 
  mutate(CD4CD8 = if_else(str_detect(minor_clusters2, "CD8"),"CD8",CD4CD8)) %>% 
  mutate(CD4CD8 = if_else(str_detect(minor_clusters2, "CD4"),"CD4",CD4CD8)) %>% 
  filter(!str_detect(CD4CD8,"other")) %>% 
  group_by(timepoint,patient,CD4CD8) %>%
  summarise(n = n()) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  dplyr::select(-n) %>% 
  spread(key = CD4CD8, value = p, fill=0) %>% 
  #mutate(dif = CD4/CD8) %>% 
  dplyr::select(-CD8) %>% 
  spread(key = timepoint, value = CD4, fill=0) %>%
  mutate(dif = post-pre) %>% 
  ggplot() + aes(x = patient, y = dif, group=patient, fill = patient) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = patient_colors)




SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("T - Cells")) %>% 
  filter(str_detect(minor_clusters2, "Tcells") | str_detect(minor_clusters2, "^NK")) %>% 
  mutate(CD4CD8 = "other") %>% 
  mutate(CD4CD8 = if_else(str_detect(minor_clusters2, "CD8"),"CD8",CD4CD8)) %>% 
  mutate(CD4CD8 = if_else(str_detect(minor_clusters2, "CD4"),"CD4",CD4CD8)) %>% 
  filter(!str_detect(CD4CD8,"other")) %>% 
  group_by(timepoint,patient,CD4CD8) %>%
  summarise(n = n()) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  dplyr::select(-n) %>% 
  filter(CD4CD8 == "CD4") ->
  temp

IMT_TME_merged %>% 
  filter(marker == "CD4") %>% 
  filter(patient %in% patients) %>% 
  dplyr::select(patient, timepoint, density) %>% 
  left_join(temp) %>% 
  ggplot() + aes(x = p, y = log(density),color = patient, group=patient,shape=timepoint) + 
  geom_point() + geom_line()
  
 
#######
##tcell signatures

tcell_sigs_df = read_csv("~/code/resources/gene_lists/tcell_sigs.csv")
azizi_sigs = read_csv("~/code/resources/gene_lists/azizi_sigs.csv")
andreatta_sigs = read_csv("~/code/resources/gene_lists/Tcell_andreatta.csv")

names(tcell_sigs_df) = paste0("meta.",names(tcell_sigs_df))
names(azizi_sigs) = paste0("azizi.",names(azizi_sigs))
names(andreatta_sigs) = paste0("andreatta.",names(andreatta_sigs))

tcell_sigs_df %>% 
  gather(key = "signature", value = "gene") %>% 
  rbind(azizi_sigs %>% 
          gather(key = "signature", value = "gene") ) %>% 
  rbind(andreatta_sigs %>% 
          gather(key = "signature", value = "gene")) %>% 
  filter(!is.na(gene)) %>% 
  filter(gene %in% row.names(SCTall_fil2_tcells)) %>% 
  base::split(f = as.factor(.$signature), drop = T) ->
  all_tcell_sigs_list

all_tcell_sigs_list = lapply(all_tcell_sigs_list, function(x) return(x %>% pull(gene)))

tcell_sig_scores = AddModuleScore(SCTall_fil2_tcells, features = all_tcell_sigs_list)
names(tcell_sig_scores@meta.data)[str_detect(names(tcell_sig_scores@meta.data), "^Cluster")] = names(all_tcell_sigs_list)
tcell_sig_scores = tcell_sig_scores@meta.data

tcell_sig_scores %>% 
  filter(minor_clusters5 != "NK Cells") %>% 
  group_by(patient, timepoint) %>% 
  mutate(n=n()) %>% 
  filter(n >= 100) %>% 
  dplyr::select(-n) %>% 
  group_by(patient) %>% 
  filter(length(unique(timepoint)) == 2) %>% 
  ungroup() %>% 
  dplyr::select(cell.names, patient, timepoint, minor_clusters5,starts_with("azizi"), starts_with("meta"), starts_with("andreatta")) %>% 
  gather(-cell.names, -patient, -timepoint , -minor_clusters5, key="signature", value = "score") %>% 
  filter(!str_detect(signature, "Mac")) %>% 
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
  top_n(8, abs(m)) %>% 
  mutate(signature = str_split(signature, "\\.", simplify = T)[,2]) %>% 
  ggplot() + aes(x=reorder(signature, m), y = m, fill=dir) + geom_bar(stat = "identity") + 
  geom_linerange(aes(ymin = m-sd, ymax = m+sd)) +
  scale_fill_manual(values = c("steelblue","firebrick")) +
  theme_minimal() + theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(-0.5,0.5), expand = c(0,0)) +
  ylab("Median Change in Score") +
  theme(axis.line.x = element_line(size=0.5)) +
  theme(axis.title.x = element_text(size=16), axis.text=element_text(size=12)) +
  theme(panel.grid = element_blank(), legend.position = "none") + coord_flip() -> tcells_sig_plot
ggsave(tcells_sig_plot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_01/tcells_sig_plot.pdf", width=8, height=5)

 
  #ggplot() + aes(x = reorder(signature,t_value), y = post-pre) + geom_boxplot() + geom_point(aes(color = patient))+
  #scale_color_manual(values = patient_colors) + coord_flip()
  #ggplot() + aes(x = reorder(signature,t_value), y = post-pre, fill = median(post-pre)) + geom_bar(stat="identity") +
  #scale_color_gradient2(low="steelblue", mid="white", high = "firebrick")





tcell_sigs_list = lapply(tcell_sigs_list, function(x) return(x %>% pull(gene)))

SCTall_fil2_tcells = AddModuleScore(SCTall_fil2_tcells, features = tcell_sigs_list, name = "tcell_sigs.")

tcell_sig_rename = c("tcell_sigs.1"="Activation_Dysfunction",
  "tcell_sigs.2"="Cell_stress",
  "tcell_sigs.3"="Exhaustion",
  "tcell_sigs.4"="Inhibitory",
  "tcell_sigs.5"="Progenitor_Exhausted",
  "tcell_sigs.6"="Terminal_differentiation",
  "tcell_sigs.7"="Terminally_Exhausted")

names(SCTall_fil2_tcells@meta.data)[match(names(tcell_sig_rename), names(SCTall_fil2_tcells@meta.data))] = paste0(tcell_sig_rename,".sig")

#score all cells
SCTall_fil2 = AddModuleScore(SCTall_fil2, features = tcell_sigs_list, name = "tcell_sigs.")

tcell_sig_rename = c("tcell_sigs.1"="Activation_Dysfunction",
                     "tcell_sigs.2"="Cell_stress",
                     "tcell_sigs.3"="Exhaustion",
                     "tcell_sigs.4"="Inhibitory",
                     "tcell_sigs.5"="Progenitor_Exhausted",
                     "tcell_sigs.6"="Terminal_differentiation",
                     "tcell_sigs.7"="Terminally_Exhausted")

names(SCTall_fil2@meta.data)[match(names(tcell_sig_rename), names(SCTall_fil2@meta.data))] = paste0(tcell_sig_rename,".sig")


VlnPlot(SCTall_fil2_tcells, ncol = 4, pt.size = 0, group.by = "minor_clusters4", features = paste0(tcell_sig_rename,".sig"), cols = tcell_colors)
VlnPlot(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2,"Tcells")], ncol = 4, pt.size = 0, group.by = "minor_clusters2", features = paste0(tcell_sig_rename,".sig"), cols = tcell_colors)


VlnPlot(SCTall_fil2_tcells, ncol = 4, pt.size = 0,
        group.by = "patient", features = paste0(tcell_sig_rename,".sig"),
        split.by = "timepoint", split.plot = T, cols = treatment_colors)
VlnPlot(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2,"Tcells")], ncol = 4, pt.size = 0,
        group.by = "patient", features = paste0(tcell_sig_rename,".sig"),
        split.by = "timepoint", split.plot = T, cols = treatment_colors)


VlnPlot(SCTall_fil2_tcells, ncol = 4, pt.size = 0,
        group.by = "patient", features = paste0(tcell_sig_rename,".sig"),
        cols = patient_colors)
VlnPlot(SCTall_fil2_tcells[,SCTall_fil2_tcells$timepoint=="pre"], ncol = 4, pt.size = 0,
        group.by = "patient", features = paste0(tcell_sig_rename,".sig"),
        cols = patient_colors)
VlnPlot(SCTall_fil2_tcells[,SCTall_fil2_tcells$timepoint=="pre"], ncol = 4, pt.size = 0,
        group.by = "group", features = paste0(tcell_sig_rename,".sig"),
        cols = group_colors)

#SCTall_fil2_tcells$group = "unclassified"
#SCTall_fil2_tcells$group[SCTall_fil2_tcells$patient %in% c("BCR01","BCR03","BCR04","BCR15")] = "high_selection"
#SCTall_fil2_tcells$group[SCTall_fil2_tcells$patient %in% c("BCR06","BCR09","BCR18","BCR20")] = "low_selection"

SCTall_fil2_tcells$irds.1 = SCTall_fil2@meta.data[colnames(SCTall_fil2_tcells),]$irds.1


SCTall_fil2@meta.data %>% 
  filter(timepoint == "pre") %>% 
  filter(str_detect(minor_clusters2,"Tcells")) %>% 
  group_by(patient) %>% 
  mutate(n = n()) %>% 
  filter(n >=20) %>% 
  summarise(TE = median(Terminally_Exhausted.sig)) %>%
  ungroup() %>% 
  left_join(SCTall_fil2@meta.data %>% 
              filter(timepoint == "pre") %>% 
              group_by(patient) %>% 
              summarise(irds = median(irds.1))) %>% 
  ggplot() + aes(x=irds, y=TE) + geom_point(aes(color=patient)) + 
  scale_color_manual(values = patient_colors) + geom_smooth(method="lm")
  

  irds_order = c("BCR04","BCR03","BCR15","BCR05","BCR13","BCR01","BCR18","BCR20","BCR06","BCR09")

SCTall_fil2_tcells@meta.data %>% 
  filter(str_detect(minor_clusters4,"CD4")) %>% 
  filter(timepoint == "post") %>% 
  #filter(str_detect(minor_clusters2,"Tcells")) %>% 
  group_by(patient) %>% 
  mutate(n = n()) %>% 
  filter(n >=20) %>% 
  ggplot() + aes(x=patient, y=Terminally_Exhausted.sig) + 
  geom_boxplot() + scale_x_discrete(limits = irds_order) + 
  theme_minimal()

  
  
  

    
############################
####tcell dotplot###########
############################


SCTall_fil2_tcells_minor_clusters4_markers_CD8 = readRDS("./Rda/SCTall_fil2_tcells_minor_clusters4_markers_CD8.rds")
SCTall_fil2_tcells_minor_clusters4_markers_CD4 = readRDS("./Rda/SCTall_fil2_tcells_minor_clusters4_markers_CD4.rds")
SCTall_fil2_tcells_minor_clusters4_markers = readRDS("./Rda/SCTall_fil2_tcells_minor_clusters4_markers.rds")


Idents(SCTall_fil2_tcells) = "minor_clusters5"

tcell_dotplot_features = c("IL7R","CCR7","TCF7","LEF1","SELL","GZMK","CMC1","CST7","GZMA","IGFL2","PDCD1","CXCL13","LAG3","ITGAE","GZMB","NKG7","PRF1","GNLY","FGFBP2","CX3CR1","FCGR3A","CCL4","CCL3","CCL4L2","ZNF683","XCL1","IL32","FOXP3","TIGIT","BATF","TNFRSF4","TNFRSF18","CTLA4","IL2RA","LAG3","KLRB1","TRAV1-2","TRBV20-1","STAT1","MX1","ISG15","IFI6")


DotPlot(SCTall_fil2_tcells[,SCTall_fil2_tcells$minor_clusters5 %in% names(tcell_colors[c(1:6,8:14,7,15)])], 
        features = unique(tcell_dotplot_features), assay = "SCT",
        group.by = "minor_clusters5",
        col.min = -1, col.max = 1, dot.min = 0.1) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = names(tcell_colors[c(1,9,6,10,2,14,3,13,4,5,11,12,8,7,15)])) + 
  scale_color_gradient2(low="lightgrey", mid="lightgrey", high="tomato") + 
  theme(axis.title = element_blank()) -> tcell_dotplot_CD4CD8
ggsave(tcell_dotplot_CD4CD8, filename = "~/projects/Breast_cancer_radiation/fig_pieces_01/tcell_dotplot_CD4CD8.pdf", device = cairo_pdf, width = 14, height = 6.5)


tcell_dotplot_features2 = c("CD3E","CD3D","CD3G","CD8A","CD4","TRAC","TRBC2","TRDC","XCL1","XCL2","AREG","GNLY","NKG7","FGFBP2","PRF1","MKI67","TUBA1B","STMN1")

SCTall_fil2_tcells$medium_clusters = SCTall_fil2$medium_clusters[colnames(SCTall_fil2_tcells)]

DefaultAssay(SCTall_fil2_tcells) = "SCT"
Idents(SCTall_fil2_tcells) = "medium_clusters"
DotPlot(SCTall_fil2_tcells[,colnames(SCTall_fil2_tcells) %in% colnames(SCTall_fil2)], 
        features = unique(tcell_dotplot_features2), assay = "SCT",
        col.min = -1, col.max = 1, dot.min = 0.1) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = c("CD4 Tcells", "CD8 Tcells","GD Tcells", "NK Cells","Cycling Immunocytes")) + 
  scale_color_gradient2(low="lightgrey", mid="lightgrey", high="tomato") + 
  theme(axis.title = element_blank()) + theme (legend.position = "none") -> tcell_dotplot_medium
ggsave(tcell_dotplot_medium, filename = "~/projects/Breast_cancer_radiation/fig_pieces_01/tcell_dotplot_medium.pdf", device = cairo_pdf, width = 6, height = 2.5)



#################
####CCL4 vln#####
#################

Idents(SCTall_fil2_tcells) = "patient"
VlnPlot(SCTall_fil2_tcells, assay = "SCT",slot = "data", features = c("sct_CCL4"), y.max = 6, cols = rainbow(18)) + NoLegend()

VlnPlot(SCTall_fil2_tcells[,str_detect(SCTall_fil2_tcells$minor_clusters5,"CD8")], assay = "SCT", features = c("CCL4"), combine = F, group.by = "minor_clusters5", cols = tcell_colors) 


data.frame(cluster = SCTall_fil2_tcells$minor_clusters5, timepoint = SCTall_fil2_tcells$timepoint, 
           IFNG = SCTall_fil2_tcells@assays$SCT["IFNG",][1,],
           TNF = SCTall_fil2_tcells@assays$SCT["TNF",][1,],
           CCL4 = SCTall_fil2_tcells@assays$SCT["CCL4",][1,],
           PDCD1 = SCTall_fil2_tcells@assays$SCT["PDCD1",][1,]) %>%
  gather(-cluster, -timepoint, key = "gene", value = "value") %>% 
  filter(timepoint == "post") %>% 
  filter(str_detect(cluster, "CD8")) %>% 
  filter(!str_detect(cluster, "MAIT|CCR7|ISG")) %>% 
  mutate(cluster = str_replace(cluster, "Tcells \\- ","")) %>% 
  ggplot() + aes(x = cluster, y = value, fill = cluster) + geom_violin(scale = "width") + facet_wrap(~gene) +
  scale_fill_manual(values = tcell_colors[1:8] %>% setNames(str_replace(names(tcell_colors[1:8]), "Tcells \\- ",""))) +
  theme(legend.position = "none") + ylab("Normalized Expression") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), strip.background = element_blank()) +
  theme(axis.title.x = element_blank(), text = element_text(size = 14)) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) -> CCL4_plot
ggsave(CCL4_plot, filename = "~/projects/Breast_cancer_radiation/sup_figure_pieces/CCL4_vln.pdf", device = cairo_pdf, width = 6, height = 4)



####################
##tcell sig table###
####################

tcell_sigs_df %>% 
  gather(key = "signature", value = "gene") %>% 
  rbind(azizi_sigs %>% 
          gather(key = "signature", value = "gene") ) %>% 
  rbind(andreatta_sigs %>% 
          gather(key = "signature", value = "gene")) %>% 
  filter(!is.na(gene)) %>% 
  filter(gene %in% row.names(SCTall_fil2_tcells)) %>% 
  group_by(signature) %>% 
  mutate(r = 1:n()) %>% 
  spread(key = signature, value = gene) %>% 
  write_csv(file = "~/projects/Breast_cancer_radiation/tcell_sigs.csv")


