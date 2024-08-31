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
library(copykat)
library(magick)
library(ComplexHeatmap)
library(circlize)
source("~/TENX/analysis/MP/MP_project_4/TCR_functions_script.R")
options(future.globals.maxSize = 60000 * 1024^2)
NewSet1 = function(n){return(colorRampPalette(brewer.pal(9,"Set1")[c(1,5,6,3,2,4,8,7,9)])(n))}


#patients_key = c("P01","P02","P03","P09","P05","P06","P10","P04","P11","P07","P08") %>% setNames(patients)
#SCTall_fil2$patient.2 = patients_key[SCTall_fil2$patient]
#patients_key_IHC = paste0("P",12:19) %>% setNames(c("BCR00","BCR07","BCR08","BCR10","BCR12","BCR14","BCR16","BCR21"))
#id_key$patient.2 = c(patients_key,patients_key_IHC)[id_key$patient]
#id_key$sc_RNA = if_else(id_key$patient %in% patients, "YES","NO")
#id_key$sc_DNA = if_else(id_key$patient %in% patients[c(1,2,3,5,6,8,10,11)], "YES","NO")
#write_csv(id_key, "./data/id_key.csv")

#clinical_data = read_csv("./data/clinical_data.csv")
#KI67_total 

medium_cluster_colors["Cycling Epithelial"] = "grey40"
medium_cluster_colors["Macrophages"] = "mediumpurple1"
medium_cluster_colors["NK Cells"] = "burlywood3"


major_cluster_colors = c("#743783","#2f9758","#e5a0f3","#f1b063","#42a0f5","#afad2b","#495faf","#f58373")
#"#9c86d7","#9f5f28"
names(major_cluster_colors) = unique(SCTall_fil2$major_clusters)

ck_colors = c("#FC766A","#184A45","#B0B8B4")
names(ck_colors) = c("aneuploid","diploid",NA)

names(patient_colors) = patients

patients.2_colors = patient_colors[clinical_data$patient]
names(patients.2_colors) = clinical_data$patient.2
patients.2_colors = patients.2_colors[!is.na(patients.2_colors)]

#########################
###add ck predictions####
#########################

ck_predictions = list()
for(i in patients){
  path = paste0("./copykat_2/",i,"_copykat_prediction.txt")
  if(file.exists(path)){
    ck_predictions[[i]] = read_tsv(path)
  }
  else{}
}

ck_predictions$BCR17$copykat.pred[ck_predictions$BCR17$copykat.pred == "aneuploid"] = NA
ck_predictions$BCR01$copykat.pred[ck_predictions$BCR01$copykat.pred == "aneuploid"] = "temp"
ck_predictions$BCR01$copykat.pred[ck_predictions$BCR01$copykat.pred == "diploid"] = "aneuploid"
ck_predictions$BCR01$copykat.pred[ck_predictions$BCR01$copykat.pred == "temp"] = "diploid"

ck_predictions_df = bind_rows(ck_predictions)
names(ck_predictions_df$copykat.pred) = ck_predictions_df$cell.names
ck_predictions_vec = c(ck_predictions_df$copykat.pred)


SCTall_fil2$ck_pred = NA
SCTall_fil2@meta.data[names(ck_predictions_vec),"ck_pred"] = ck_predictions_vec

#table(SCTall_fil2$major_clusters, SCTall_fil2$ck_pred, SCTall_fil2$patient)
#table(SCTall_fil2$minor_clusters2, SCTall_fil2$ck_pred, SCTall_fil2$sample)[,,"BCR01_2"]

ck_predictions_vec_epi = ck_predictions_vec[names(ck_predictions_vec) %in% colnames(SCTall_fil2_epi)]
SCTall_fil2_epi$ck_pred = NA
SCTall_fil2_epi@meta.data[names(ck_predictions_vec_epi),"ck_pred"] = ck_predictions_vec_epi



########################
####all cluster umaps###
########################

DimPlot(SCTall_fil2, group.by = "major_clusters", cols = major_cluster_colors) -> major_clusters_umap
ggsave(major_clusters_umap, filename = "./plots/fig_pieces_05/major_clusters_umap.pdf", width=6, height=4, device = cairo_pdf)

DimPlot(SCTall_fil2, group.by = "patient.2", cols = patients.2_colors, shuffle = T) -> all_patient_umap
ggsave(all_patient_umap, filename = "~/projects/Breast_cancer_radiation/fig_pieces_05/all_patient_umap.pdf", width=5.5, height=4, device = cairo_pdf)

set.seed(77)
DimPlot(SCTall_fil2, group.by = "timepoint", cols = treatment_colors, shuffle = T) -> all_timepoint_umap
ggsave(all_timepoint_umap, filename = "~/projects/Breast_cancer_radiation/fig_pieces_05/all_timepoint_umap.pdf", width=6, height=4, device = cairo_pdf)

DimPlot(SCTall_fil2, group.by = "ck_pred", cols = ck_colors) -> all_ck_umap
ggsave(all_ck_umap, filename = "./plots/fig_pieces_05/all_ck_umap.pdf", width=6, height=4, device = cairo_pdf)


plot_grid(major_clusters_umap, all_patient_umap, all_timepoint_umap, all_ck_umap, ncol = 1, align = "v", axis = "r") -> clusters_umap_combined
ggsave(clusters_umap_combined, filename = "./plots/fig_pieces_05/clusters_umap_combined.pdf", width=6, height=16, device = cairo_pdf)


set.seed(77)
DimPlot(SCTall_fil2, group.by = "medium_clusters", cols = medium_cluster_colors, shuffle = T) -> medium_clusters_umap
ggsave(medium_clusters_umap, filename = "~/projects/Breast_cancer_radiation/fig_pieces_05/medium_clusters_umap.pdf", width=6, height=4, device = cairo_pdf)

plot_grid(medium_clusters_umap, all_patient_umap, all_timepoint_umap, all_ck_umap, ncol = 1, align = "v", axis = "r") -> clusters_umap_combined2
ggsave(clusters_umap_combined2, filename = "~/projects/Breast_cancer_radiation/fig_pieces_05/clusters_umap_combined2.pdf", width=6, height=16, device = cairo_pdf)



###########################
####all cluster barplots###
###########################


SCTall_fil2@meta.data %>% 
  ggplot() + aes(x=timepoint, fill=major_clusters) + 
  geom_bar(position="fill") + facet_wrap(~patient, nrow = 1) +
  scale_fill_manual(values = major_cluster_colors) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), strip.background = element_blank()) +
  theme(axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(strip.text = element_text(size=12), axis.text = element_text(size=12)) +
  theme(panel.spacing = unit(0, "lines")) ->
  major_clusters_barplot

ggsave(major_clusters_barplot, filename = "./plots/fig_pieces_05/major_clusters_barplot.pdf", width=7, height=3, device = cairo_pdf)


SCTall_fil2@meta.data %>% 
  #filter(major_clusters == "Epithelial") %>% 
  filter(!is.na(ck_pred)) %>% 
  ggplot() + aes(x=timepoint, fill=ck_pred) + 
  geom_bar(position="fill") + facet_wrap(~patient, nrow = 1) +
  scale_fill_manual(values = ck_colors) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), strip.background = element_blank()) +
  theme(axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(strip.text = element_text(size=12), axis.text = element_text(size=12)) +
  theme(panel.spacing = unit(0, "lines"))





SCTall_fil2@meta.data %>% 
  filter(major_clusters %in% c("B - Cells","NK Cells","T - Cells","Myeloid")) %>% 
  ggplot() + aes(x=timepoint, fill=major_clusters) + 
  geom_bar(position="fill") + facet_wrap(~patient, nrow = 1) +
  scale_fill_manual(values = major_cluster_colors) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), strip.background = element_blank()) +
  theme(axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(strip.text = element_text(size=12), axis.text = element_text(size=12)) +
  theme(panel.spacing = unit(0, "lines")) ->
  immune_clusters_barplot

ggsave(immune_clusters_barplot, filename = "./plots/fig_pieces_05/immune_clusters_barplot.pdf", width=7, height=3, device = cairo_pdf)




SCTall_fil2@meta.data %>% 
  rbind(SCTall_fil2@meta.data %>% mutate(patient.2 = "all")) %>% 
  ggplot() + 
  aes(x=timepoint, 
    fill=factor(medium_clusters, 
    levels = names(medium_cluster_colors))) + 
  geom_bar(position="fill") + facet_wrap(~patient.2, nrow = 1) +
  scale_fill_manual(values = medium_cluster_colors) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust =1)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), strip.background = element_blank()) +
  theme(axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(strip.text = element_text(size=16), axis.text = element_text(size=12)) +
  theme(panel.spacing = unit(0, "lines")) ->
  medium_clusters_barplot

ggsave(medium_clusters_barplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_05/medium_clusters_barplot.pdf", width=7.5, height=2.75, device = cairo_pdf)





SCTall_fil2@meta.data %>% 
  rbind(SCTall_fil2@meta.data %>% mutate(patient.2 = "all")) %>% 
  filter(major_clusters %in% c("B - Cells","NK Cells","T - Cells","Myeloid")) %>% 
  ggplot() + 
  aes(x=timepoint, 
      fill=factor(medium_clusters, 
                  levels = names(medium_cluster_colors))) + 
  geom_bar(position="fill") + facet_wrap(~patient.2, nrow = 1) +
  scale_fill_manual(values = medium_cluster_colors) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust =1)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), strip.background = element_blank()) +
  theme(axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(strip.text = element_text(size=16), axis.text = element_text(size=12)) +
  theme(panel.spacing = unit(0, "lines")) ->
  medium_clusters_immune_barplot

ggsave(medium_clusters_immune_barplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_05/medium_clusters_immune_barplot.pdf", width=7.5, height=2.75, device = cairo_pdf)







DimPlot(SCTall_fil2, group.by = "major_clusters", cols = major_cluster_colors)
DimPlot(SCTall_fil2, group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2, group.by = "timepoint", cols = treatment_colors)
DimPlot(SCTall_fil2, group.by = "nFeature_SCT")


table(SCTall_fil2$minor_clusters2,SCTall_fil2$major_clusters)




DimPlot(SCTall_fil2[,SCTall_fil2$major_clusters == "Epithelial"], group.by = "ck_pred")
DimPlot(SCTall_fil2[,SCTall_fil2$major_clusters == "Epithelial"], group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2_epi, group.by = "ck_pred")
DimPlot(SCTall_fil2_epi, group.by = "patient", cols = patient_colors)
DimPlot(SCTall_fil2_epi, group.by = "timepoint", cols = treatment_colors)




####
###





SCTall_fil2@meta.data %>%
  #filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, medium_clusters) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  #filter(str_detect(minor_clusters2, "Tcells")|str_detect(minor_clusters2, "^NK")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(p = post - pre) %>% 
  group_by(medium_clusters) %>% 
  summarise(m = median(p)) %>% 
  arrange(-m) %>% 
  pull(medium_clusters) ->
  medium_clusters_change_order


SCTall_fil2@meta.data %>%
  #filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  #filter(major_clusters %in% c("T - Cells", "NK Cells")) %>% 
  group_by(sample, timepoint,patient,cohort, medium_clusters) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  #filter(str_detect(minor_clusters2, "Tcells")|str_detect(minor_clusters2, "^NK")) %>% 
  #filter(!str_detect(minor_clusters2, "Cycling")) %>% 
  dplyr::select(-sample,-n,-cohort) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(p = post - pre) %>% 
  mutate(medium_clusters = factor(medium_clusters, levels = medium_clusters_change_order)) %>% 
  ggplot() + aes(x = medium_clusters, y=p) + 
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
  guides(color = guide_legend(override.aes = list(size = 2, shape=15)))



for(i in medium_clusters_change_order){
  SCTall_fil2@meta.data %>%
    #filter(group == "low_selection") %>% 
    mutate(patient = as.character(patient)) %>% 
    #filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
    group_by(sample, timepoint,patient,cohort, medium_clusters) %>%
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(sample) %>% 
    mutate(p = n/sum(n)) %>% 
    ungroup() %>% 
    dplyr::select(-sample,-cohort,-n) %>% 
    complete(timepoint, patient, medium_clusters, fill = list(p=0)) %>% 
    filter(str_detect(medium_clusters, i)) %>% 
    #filter(str_detect(minor_clusters2, "Tcells - CD8:CCL4")) %>% 
    spread(key = timepoint, value = p, fill = 0) -> temp
  print(i)
  print(wilcox.test(temp$pre, temp$post, paired = T)$p.value)
}


####################
##receptor heatmap##
####################

library(circlize)

clinical_data %>% 
  #mutate(ER = if_else(str_detect(receptor_status,"ER\\+"),"+","-")) %>%   
 # mutate(PR = if_else(str_detect(receptor_status,"PR\\+"),"+","-")) %>%  
  mutate(HER2 = if_else(str_detect(receptor_status,"HER2\\+"),"+","-")) %>%
  arrange(patient.2) %>% 
  dplyr::select(patient.2, sc_RNA, sc_DNA, ER_percent, PR_percent, HER2, KI67, clinical_stage) %>% 
  column_to_rownames(var = "patient.2") ->
  clinical_subset 
  

names(clinical_subset)[names(clinical_subset) == "clinical_stage"] = "stage"
names(clinical_subset)[names(clinical_subset) == "ER_percent"] = "ER%"
names(clinical_subset)[names(clinical_subset) == "PR_percent"] = "PR%"


temp_mat = matrix(row.names(clinical_subset), dimnames = list(row.names(clinical_subset)))

ha = rowAnnotation(
  df = clinical_subset,
  col = list(sc_RNA = c("YES" = "coral", "NO" = "grey"),
             sc_DNA = c("YES" = "coral", "NO" = "grey"),
             `ER%` = colorRamp2(c(0,100), c("white", "khaki2")),
             `PR%` = colorRamp2(c(0,100), c("white", "khaki2")),
             HER2 = c("+" = "khaki2", "-" = "grey"),
             KI67 = colorRamp2(c(0,80), c("white", "chartreuse4")),
             stage = c("IA" = "#A0B3F0","IB" = "#4169E1", "IIA" = "#C3B3E6", "IIB" = "#8968CD", "IIIB" = "#E78FC7", "IIIC" = "#D02090")),
  gp = gpar(col = "white")
)

Heatmap(matrix = temp_mat, col = rep("white", nrow(temp_mat)) %>% set_names(temp_mat),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(temp_mat[i, j], x, y, gp = gpar(fontsize = 10))
        }, show_heatmap_legend = F) + ha -> clinical_heatmap

pdf("~/projects/Breast_cancer_radiation/fig_pieces_05/clinical_heatmap_raw.pdf", width=2.6, height = 6)
draw(clinical_heatmap)
dev.off()









########


Idents(SCTall_fil2) = "medium_clusters"
medium_cluster_markers_SCT = FindAllMarkers(SCTall_fil2, assay = "SCT", only.pos = T, logfc.threshold = 1)

medium_cluster_markers_SCT %>% 
  group_by(cluster) %>% 
  top_n(5, avg_logFC) %>%  View()

medium_cluster_markers_select = c("SCGB2A2","EPCAM","KRT81","TUBB","STMN1","COL1A1","DCN","IGFBP7","SPARCL1","TAGLN","MYL9","CD74","HLA-DRA","CD68","C1QC","S100A8","S100A9","GNLY","NKG7","TRDC","TRAC","CD3E","CD8A","CD4","CD79A","MS4A1","IGKC","MZB1","MKI67")

DotPlot(SCTall_fil2, features = medium_cluster_markers_select, 
        group.by = "medium_clusters", 
        cols = c("grey","black"), col.min = 0) + 
  scale_y_discrete(limits = rev(names(medium_cluster_colors))) +
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle=90, hjust =1, vjust = 0.5)) -> medium_cluster_markers_dotplot

ggsave(medium_cluster_markers_dotplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_05/medium_cluster_markers_dotplot.pdf", device = cairo_pdf, width=10, height=4.5)


###########################
###show epcam expression###
###########################


FeaturePlot(SCTall_fil2, features = "EPCAM", cols = c("papayawhip", "orchid"), max.cutoff = "q50", raster = F, order = T) + ylab("UMAP 2") + xlab("UMAP 1") -> temp
ggsave(temp, filename = "~/projects/Breast_cancer_radiation/fig_pieces_05/all_EPCAM_umap.pdf", width=6, height=4, device = cairo_pdf)


#########
##t.test major cluster freq
####

SCTall_fil2@meta.data %>%
  filter(major_clusters %in% c("Endothelial", "Fibroblasts","Perivascular Cells")) %>% 
  #filter(major_clusters %in% c("B - Cells", "Myeloid","T - Cells", "NK Cells")) %>% 
  #filter(major_clusters %in% c("T - Cells", "NK Cells")) %>% 
  group_by(timepoint,patient.2, minor_clusters2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(timepoint,patient.2) %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup() %>% 
  dplyr::select(-n) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(p = post - pre) %>% 
  group_by(minor_clusters2) %>% 
  mutate(p_value = wilcox.test(unlist(p))$p.value,
         t_value = wilcox.test(unlist(p))$statistic) %>% 
  select(minor_clusters2, p_value) %>% 
  unique() %>% 
  ungroup() %>% 
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>% 
  filter(p_value <= 0.05)


