
library(SingleCellSignalR)
library(broom)
library(tidyverse)
library(Seurat)
library(cowplot)


#########
##scsignalr function
RunSCSignalr = function(so, genes, cluster){
  cs = cell_signaling(data = so@assays$RNA@counts[genes,], 
                      int.type = "paracrine", s.score = 0, 
                      genes = genes, 
                      cluster = as.numeric(factor(so@meta.data[,cluster])), 
                      c.names = levels(factor(so@meta.data[,cluster])), 
                      species = "homo sapiens")
  
  cs_df = cs
  for(i in 1:length(cs_df)){
    test = cs_df[[i]]
    test$ligand_cluster = names(test)[1]
    test$receptor_cluster = names(test)[2]
    names(test)[1] = "ligand"
    names(test)[2] = "receptor"
    cs_df[[i]] = test
  }
  
  s_cs_df = do.call(rbind, cs_df)
  row.names(s_cs_df) = NULL
  
  return(s_cs_df)
}



#########
##on medium clusters

genes = row.names(so)[row.names(so) %in% c(LRdb$ligand, LRdb$receptor)]

cs_medium_pre_list = list()
cs_medium_post_list = list()
for(i in patients){
  print(i)
  so = SCTall_fil2[,SCTall_fil2$patient == i & SCTall_fil2$timepoint == "pre"]
  cs_medium_pre_list[[i]] = RunSCSignalr(so, genes, cluster = "medium_clusters")
  so = SCTall_fil2[,SCTall_fil2$patient == i & SCTall_fil2$timepoint == "post"]
  cs_medium_post_list[[i]] = RunSCSignalr(so, genes, cluster = "medium_clusters")
}

cs_medium_pre_df = bind_rows(cs_medium_pre_list, .id="patient")
cs_medium_pre_df$timepoint = "pre"
cs_medium_post_df = bind_rows(cs_medium_post_list, .id="patient")
cs_medium_post_df$timepoint = "post"

cs_medium_df = rbind(cs_medium_pre_df, cs_medium_post_df) 



#cs_medium_df %>% 
#  filter(timepoint == "pre") %>% 
#  dplyr::select(-`interaction type`, -LRscore) %>% 
#  left_join(selection_df) %>% 
#  group_by(ligand_cluster,receptor_cluster, group) %>% 
#  summarise(n = n()) %>% 
#  filter(group != "unclassified") %>% 
#  mutate(group = as.character(group)) %>% 
#  ungroup() %>% 
#  complete(ligand_cluster, receptor_cluster, group, fill = list(n=0)) %>% 
#  spread(key = group, value = n) %>% 
#  mutate(d = high_selection-low_selection) %>% 
#  dplyr::select(-high_selection, -low_selection) %>% 
#  ggplot() + aes(x = ligand_cluster, y = receptor_cluster, fill=d) + geom_tile() +
#  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
#  scale_fill_gradient2(low="steelblue", high="firebrick", mid = "white")



cs_medium_df %>% 
  filter(timepoint == "pre") %>% 
  filter(!str_detect(ligand_cluster, "Cyc")) %>% 
  filter(!str_detect(receptor_cluster, "Cyc")) %>% 
  filter(LRscore >= 0.5) %>% 
  dplyr::select(-`interaction type`, -LRscore) %>% 
  left_join(selection_df) %>% 
  group_by(ligand_cluster,receptor_cluster, group, patient) %>% 
  summarise(n = n()) %>% 
  filter(group != "unclassified") %>% 
  mutate(group = as.character(group)) %>% 
  ungroup() %>% 
  complete(ligand_cluster, receptor_cluster, nesting(group, patient), fill = list(n=0)) %>% 
  group_by(ligand_cluster, receptor_cluster, group) %>% 
  summarise(m = mean(n)) %>% 
  spread(key = group, value = m) %>% 
  mutate(d = low_selection-high_selection) %>% 
  dplyr::select(-high_selection, -low_selection) -> temp

  temp = data.frame(spread(temp, key = receptor_cluster, value = d))
  row.names(temp) = temp$ligand_cluster
  temp = temp[,-1]
  x_temp = row.names(temp)[hclust(dist(temp))$order]
  y_temp = row.names(temp)[hclust(dist(t(temp)))$order]
  
  
  
  cs_medium_df %>% 
    filter(timepoint == "pre") %>% 
    filter(!str_detect(ligand_cluster, "Cyc")) %>% 
    filter(!str_detect(receptor_cluster, "Cyc")) %>% 
    filter(LRscore >= 0.5) %>% 
    dplyr::select(-`interaction type`, -LRscore) %>% 
    left_join(selection_df) %>% 
    group_by(ligand_cluster,receptor_cluster, group, patient) %>% 
    summarise(n = n()) %>% 
    filter(group != "unclassified") %>% 
    mutate(group = as.character(group)) %>% 
    ungroup() %>% 
    complete(ligand_cluster, receptor_cluster, nesting(group, patient), fill = list(n=0)) %>% 
    group_by(ligand_cluster, receptor_cluster, group) %>% 
    summarise(m = mean(n)) %>% 
    spread(key = group, value = m) %>% 
    mutate(d = low_selection-high_selection) %>% 
    dplyr::select(-high_selection, -low_selection) %>% 
  ggplot() + aes(x = ligand_cluster, y = receptor_cluster, fill=d) + geom_tile() +
  theme_real_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_fill_gradient2(low="steelblue", high="firebrick", mid = "white") + 
  theme(legend.title = element_blank()) + scale_x_discrete(limits = x_temp) +
    scale_y_discrete(limits = y_temp)  -> temp
ggsave(temp, filename = "~/projects/Breast_cancer_radiation/fig_pieces_undecided/RL_group_pre.pdf", device = cairo_pdf(), height = 4, width = 4.8)
  
  

cs_medium_df %>% 
  filter(timepoint == "post") %>% 
  filter(!str_detect(ligand_cluster, "Cyc")) %>% 
  filter(!str_detect(receptor_cluster, "Cyc")) %>% 
  filter(LRscore >= 0.5) %>% 
  dplyr::select(-`interaction type`, -LRscore) %>% 
  left_join(selection_df) %>% 
  group_by(ligand_cluster,receptor_cluster, group, patient) %>% 
  summarise(n = n()) %>% 
  filter(group != "unclassified") %>% 
  mutate(group = as.character(group)) %>% 
  ungroup() %>% 
  complete(ligand_cluster, receptor_cluster, nesting(group, patient), fill = list(n=0)) %>% 
  group_by(ligand_cluster, receptor_cluster, group) %>% 
  summarise(m = mean(n)) %>% 
  spread(key = group, value = m) %>% 
  mutate(d = low_selection-high_selection) %>% 
  dplyr::select(-high_selection, -low_selection) %>% 
  ggplot() + aes(x = ligand_cluster, y = receptor_cluster, fill=d) + geom_tile() +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_fill_gradient2(low="steelblue", high="firebrick", mid = "white") + 
  theme(legend.title = element_blank())



cs_medium_df %>% 
  filter(timepoint == "pre") %>% 
  filter(!str_detect(ligand_cluster, "Cyc")) %>% 
  filter(!str_detect(receptor_cluster, "Cyc")) %>% 
  #filter(LRscore >= 0.5) %>% 
  unique() %>% 
  dplyr::select(-`interaction type`, -timepoint) %>% 
  left_join(selection_df) %>% 
  filter(group != "unclassified") %>% 
  group_by(ligand,receptor,ligand_cluster,receptor_cluster,group) %>% 
  summarise(m = mean(LRscore)) %>% 
  ungroup() %>% 
  spread(key = group, value = m, fill = 0) %>% 
  mutate(d = low_selection-high_selection) %>% 
  View()



###########################
####selection heatmap
######################

library(ComplexHeatmap)

selection_order = irds_order[c(1,2,3,6,7:10)]


all_genes_tum_avg = AverageExpression(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "BCR")],
                                      add.ident = "timepoint",return.seurat = T)

all_genes_tum_avg$timepoint = str_sub(colnames(all_genes_tum_avg),7)
all_genes_tum_avg$patient = str_sub(colnames(all_genes_tum_avg),1,5)


#all_genes_tum_avg = readRDS("./Rda/all_genes_tum_avg.rds")


reshape2::melt(all_genes_tum_avg@assays$SCT@data) %>% 
  rename("gene" = Var1, "sample" = Var2) %>% 
  mutate(patient = str_sub(sample, 1, 5), timepoint = str_sub(sample, 7)) %>% 
  left_join(selection_df) %>% 
  filter(timepoint == "pre") %>% 
  mutate(group = as.character(group)) %>% 
  filter(group != "unclassified") %>% 
  #mutate(group = as.factor(group)) %>% 
  dplyr::select(-sample, -timepoint) %>%
  group_by(gene) %>%
  do(tidy(t.test(value ~ group, data = .))) -> sel_tum_pre_avg_ttests


sel_tum_pre_avg_ttests %>%
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  mutate(gene = as.character(gene)) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p.adj <= 0.05) %>% 
  filter(estimate1 >= 0.25 | estimate2 >= 0.25) %>% 
  arrange(p.value)  %>% 
  mutate(dir = if_else(statistic > 0, "high","low")) %>%
  group_by(dir) %>% 
  #top_n(2000, -p.value) %>% 
  arrange(dir, p.value) %>% 
  pull(gene) -> sel_tum_pre_avg_top



reshape2::melt(all_genes_tum_avg@assays$SCT@data) %>% 
  rename("gene" = Var1, "sample" = Var2) %>% 
  mutate(patient = str_sub(sample, 1, 5), timepoint = str_sub(sample, 7)) %>% 
  left_join(selection_df) %>% 
  filter(timepoint == "post") %>% 
  mutate(group = as.character(group)) %>% 
  filter(group != "unclassified") %>% 
  #mutate(group = as.factor(group)) %>% 
  dplyr::select(-sample, -timepoint) %>%
  group_by(gene) %>%
  do(tidy(t.test(value ~ group, data = .))) -> sel_tum_post_avg_ttests

sel_tum_post_avg_ttests %>%
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  mutate(gene = as.character(gene)) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p.adj <= 0.05) %>% 
  filter(estimate1 >= 0.25 | estimate2 >= 0.25) %>% 
  arrange(p.value)  %>% 
  mutate(dir = if_else(statistic > 0, "high","low")) %>%
  group_by(dir) %>% 
  #top_n(2000, -p.value) %>% 
  arrange(dir, p.value) %>% 
  pull(gene) -> sel_tum_post_avg_top




mat1 = all_genes_tum_avg[sel_tum_pre_avg_top[!sel_tum_pre_avg_top %in% sel_tum_post_avg_top],all_genes_tum_avg$timepoint == "pre"]@assays$SCT@scale.data
ha = HeatmapAnnotation(group = rep(names(group_colors[1:2]), each = 4), patient = selection_order, col = list(group = rep(group_colors[1:2], each = 4), patient = patient_colors[selection_order]))
har1 = rowAnnotation(timepoint = rep("pre", length(sel_tum_pre_avg_top[!sel_tum_pre_avg_top %in% sel_tum_post_avg_top])), col = list(timepoint = treatment_colors[1]), show_annotation_name = F, show_legend = F)
hm1 = Heatmap(mat1[,paste0(selection_order,"_pre")], 
              row_order = sel_tum_pre_avg_top[!sel_tum_pre_avg_top %in% sel_tum_post_avg_top], 
              column_order = paste0(selection_order, "_pre"), show_column_names = F, top_annotation = ha,
              left_annotation =  har1,
              heatmap_legend_param = list(at = c(-4,-2, 0, 2, 4),labels = c("-4","-2", "0", "2","4"),title = "Average\nExpression"))







mat2 = all_genes_tum_avg[sel_tum_pre_avg_top[sel_tum_pre_avg_top %in% sel_tum_post_avg_top],all_genes_tum_avg$timepoint == "pre"]@assays$SCT@scale.data
har2 = rowAnnotation(timepoint = rep("pre", length(sel_tum_pre_avg_top[sel_tum_pre_avg_top %in% sel_tum_post_avg_top])), col = list(timepoint = treatment_colors[1]), show_annotation_name = F, show_legend = F)
hm2 = Heatmap(mat2[,paste0(selection_order,"_pre")], row_order = sel_tum_pre_avg_top[sel_tum_pre_avg_top %in% sel_tum_post_avg_top], 
              column_order = paste0(selection_order, "_pre"), show_column_names = F, show_heatmap_legend = F,
              heatmap_legend_param = list(at = c(-4,-2, 0, 2,4),labels = c("-4","-2", "0", "2","4"),
                                          title = "Average\nExpression"), left_annotation = har2)

mat3 = all_genes_tum_avg[sel_tum_pre_avg_top[sel_tum_pre_avg_top %in% sel_tum_post_avg_top],all_genes_tum_avg$timepoint == "post"]@assays$SCT@scale.data
har3 = rowAnnotation(timepoint = rep("post", length(sel_tum_pre_avg_top[sel_tum_pre_avg_top %in% sel_tum_post_avg_top])), col = list(timepoint = treatment_colors[2]), show_annotation_name = F, show_legend = F)
hm3 = Heatmap(mat3[,paste0(selection_order,"_post")], row_order = sel_tum_pre_avg_top[sel_tum_pre_avg_top %in% sel_tum_post_avg_top], 
              column_order = paste0(selection_order, "_post"), show_column_names = F, show_heatmap_legend = F,
              heatmap_legend_param = list(at = c(-4,-2, 0, 2,4),labels = c("-4","-2", "0", "2","4"),
                                          title = "Average\nExpression"), left_annotation = har3)       

mat4 = all_genes_tum_avg[sel_tum_post_avg_top[!sel_tum_post_avg_top %in% sel_tum_pre_avg_top],all_genes_tum_avg$timepoint == "post"]@assays$SCT@scale.data
har4 = rowAnnotation(timepoint = rep("post", length(sel_tum_post_avg_top[!sel_tum_post_avg_top %in% sel_tum_pre_avg_top])), col = list(timepoint = treatment_colors[2]), show_annotation_name = F, show_legend = F)
hm4 = Heatmap(mat4[,paste0(selection_order,"_post")], row_order = sel_tum_post_avg_top[!sel_tum_post_avg_top %in% sel_tum_pre_avg_top], 
              column_order = paste0(selection_order, "_post"), show_column_names = F, show_heatmap_legend = F,
              heatmap_legend_param = list(at = c(-4,-2, 0, 2,4),labels = c("-4","-2", "0", "2","4"),
                                          title = "Average\nExpression"), left_annotation = har4)

hm1 %v% hm2 %v% hm3 %v% hm4 -> tum_sel_heatmap

pdf("~/projects/Breast_cancer_radiation/fig_pieces_06/tum_sel_heatmap.pdf", width = 5, height = 26)
tum_sel_heatmap
dev.off()




#####################
####volcano plot#####
#####################

tumor_genes = read_csv("~/code/resources/gene_lists/tumor_genes.csv")


library(ggrepel)

plot(log2(sel_tum_pre_avg_ttests$estimate1/sel_tum_pre_avg_ttests$estimate2),-log10(sel_tum_pre_avg_ttests$p.value))

add = min(sel_tum_pre_avg_ttests$estimate1[sel_tum_pre_avg_ttests$estimate1 != 0])

sel_tum_genes_top = c("IFITM1","IFIT1","IFI35","SLCO4A1","IFI27","IFIT3","IFIT2","CDK12","ESR1","CDK10","STAT1","MUC2","ASCC1","ZNF7","FOXO1","PTPRE","EGR1","CDK12","SUMO1")

sel_tum_pre_avg_ttests %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  mutate(estimate1 = if_else( estimate1 == 0, add, estimate1)) %>% 
  mutate(estimate2 = if_else( estimate2 == 0, add, estimate2)) %>% 
  mutate(mean_dif = estimate1/estimate2) %>%
  mutate(de = if_else(log2(mean_dif) >= 1 & p.value <=0.05, "HGS","no")) %>% 
  mutate(de = if_else(log2(mean_dif) <= -1 & p.value <=0.05, "LGS",de)) %>%
  mutate(de = na_if(de, "no")) %>% 
  arrange(de,-mean_dif) %>% 
  #mutate(gene = replace(gene, which(abs(log2(mean_dif)) <= 1.42), NA)) %>% 
  #mutate(gene = replace(gene, which(p.value >= 0.011), NA)) %>% 
  mutate(gene = replace(gene, which(!gene %in% sel_tum_genes_top), NA)) %>% 
  ggplot() + aes(x = log2(mean_dif), y=-log10(p.value), color=de, label=gene) + geom_point() +
  xlab("Mean Average Expression Difference") +
  geom_text_repel(min.segment.length = 0.1, max.overlaps = Inf) + theme_real_minimal() + 
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 16)) +
  scale_color_manual(values = group_colors2[1:2], na.value = "#323D44") +
  scale_x_continuous(limits = c(-15,15)) + 
  scale_y_continuous(limits = c(0,4.5)) + 
  ggtitle("PRE") + theme(plot.title = element_text(size=20, hjust = 0.5)) +
  theme(axis.line = element_line(size=0.5)) -> sel_volcano_pre


 ggsave(sel_volcano_pre, filename = "~/projects/Breast_cancer_radiation/fig_pieces_06/sel_volcano_pre2.pdf", device = cairo_pdf, height = 6, width = 7)

 
 
sel_tum_genes_top = c("SRXN1","CLDN4","NFKBIA","CDKN1A","EFHD1","CYC1","ACAA2","IL3RA","GPR85","VSIG8","ACRBP","ADAM28","GYG2","ESR1","SOS1","CLEC2D","C2","MAP3K7","IFNAR1","RB1","MDM4","NLRC3","CDK3")
 
 
sel_tum_post_avg_ttests %>% 
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  mutate(estimate1 = if_else( estimate1 == 0, add, estimate1)) %>% 
  mutate(estimate2 = if_else( estimate2 == 0, add, estimate2)) %>% 
  mutate(mean_dif = estimate1/estimate2) %>% 
  mutate(de = if_else(log2(mean_dif) >= 1 & p.value <=0.05, "HGS","no")) %>% 
  mutate(de = if_else(log2(mean_dif) <= -1 & p.value <=0.05, "LGS",de)) %>% 
  mutate(de = na_if(de, "no")) %>% 
  arrange(de,mean_dif) %>% 
  #mutate(gene = replace(gene, which(abs(log2(mean_dif)) <= 1.88), NA)) %>% 
  #mutate(gene = replace(gene, which(p.value >= 0.0018), NA)) %>% 
  mutate(gene = replace(gene, which(!gene %in% sel_tum_genes_top), NA)) %>%  
  ggplot() + aes(x = log2(mean_dif), y=-log10(p.value), color=de, label=gene) + geom_point() +
  xlab("Mean Average Expression Difference") +
  geom_text_repel(min.segment.length = 0.1, max.overlaps = Inf) + theme_real_minimal() + 
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 16)) +
  scale_color_manual(values = group_colors2[1:2], na.value = "#323D44") +
  scale_x_continuous(limits = c(-15,15)) + 
  scale_y_continuous(limits = c(0,4.5)) + 
  ggtitle("POST") + theme(plot.title = element_text(size=20, hjust = 0.5)) +
  theme(axis.line = element_line(size=0.5)) -> sel_volcano_post
ggsave(sel_volcano_post, filename = "~/projects/Breast_cancer_radiation/fig_pieces_06/sel_volcano_post2.pdf", device = cairo_pdf, height = 6, width = 7)






ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")






###########################
####average diff_groups####
###########################

library(broom)
library(ComplexHeatmap)

selection_order = irds_order[c(1,2,3,6,7:10)]


so = SCTall_fil2[,str_detect(SCTall_fil2$medium_clusters, "Tcells|NK Cells")]
so = SCTall_fil2[,str_detect(SCTall_fil2$medium_clusters, "Dendritic Cells|Macrophages")]


avg_sel_heatmap = function(so){
  
  pre_patients = names(table(so$patient))[table(so$patient, so$timepoint)[,"pre"] >= 25]
  post_patients = names(table(so$patient))[table(so$patient, so$timepoint)[,"post"] >= 25]
  both_patients = pre_patients[pre_patients %in% post_patients]
  
  so = so[,(so$patient %in% pre_patients & so$timepoint == "pre") | (so$patient %in% post_patients & so$timepoint == "post")]
  
  Idents(so) = "patient"
  all_genes_avg = AverageExpression(so, add.ident = "timepoint",return.seurat = T)
  
  all_genes_avg$timepoint = str_sub(colnames(all_genes_avg),7)
  all_genes_avg$patient = str_sub(colnames(all_genes_avg),1,5)
  
  reshape2::melt(all_genes_avg@assays$SCT@data) %>% 
    rename("gene" = Var1, "sample" = Var2) %>% 
    mutate(patient = str_sub(sample, 1, 5), timepoint = str_sub(sample, 7)) %>% 
    left_join(selection_df) %>% 
    filter(timepoint == "pre") %>% 
    mutate(group = as.character(group)) %>% 
    filter(group != "unclassified") %>% 
    #mutate(group = as.factor(group)) %>% 
    dplyr::select(-sample, -timepoint) %>%
    group_by(gene) %>% 
    do(tidy(t.test(value ~ group, data = .))) -> sel_pre_avg_ttests
  
  sel_pre_avg_ttests %>%
    filter(!str_detect(gene, "^RP[LS]")) %>% 
    mutate(gene = as.character(gene)) %>% 
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
    filter(p.adj <= 0.05) %>% 
    filter(estimate1 >= 0.2 | estimate2 >= 0.2) %>% 
    arrange(p.value)  %>% 
    mutate(dir = if_else(statistic > 0, "high","low")) %>%
    group_by(dir) %>% 
    #top_n(2000, -p.value) %>% 
    arrange(dir, p.value) %>% 
    pull(gene) -> sel_pre_avg_top
  
  reshape2::melt(all_genes_avg@assays$SCT@data) %>% 
    rename("gene" = Var1, "sample" = Var2) %>% 
    mutate(patient = str_sub(sample, 1, 5), timepoint = str_sub(sample, 7)) %>% 
    left_join(selection_df) %>% 
    filter(timepoint == "post") %>% 
    mutate(group = as.character(group)) %>% 
    filter(group != "unclassified") %>% 
    #mutate(group = as.factor(group)) %>% 
    dplyr::select(-sample, -timepoint) %>%
    group_by(gene) %>%
    do(tidy(t.test(value ~ group, data = .))) -> sel_post_avg_ttests
  
  sel_post_avg_ttests %>%
    filter(!str_detect(gene, "^RP[LS]")) %>% 
    mutate(gene = as.character(gene)) %>% 
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
    filter(p.adj <= 0.05) %>% 
    filter(estimate1 >= 0.2 | estimate2 >= 0.2) %>% 
    arrange(p.value)  %>% 
    mutate(dir = if_else(statistic > 0, "high","low")) %>%
    group_by(dir) %>% 
    #top_n(2000, -p.value) %>% 
    arrange(dir, p.value) %>% 
    pull(gene) -> sel_post_avg_top
  
  
  
  
  reshape2::melt(all_genes_avg@assays$SCT@data) %>% 
    rename("gene" = Var1, "sample" = Var2) %>% 
    mutate(patient = str_sub(sample, 1, 5), timepoint = str_sub(sample, 7)) %>% 
    filter(patient %in%  both_patients) %>% 
    left_join(selection_df) %>% 
    filter(group != "unclassified") %>% 
    mutate(group = as.character(group)) %>%
    dplyr::select(-sample) %>%
    spread(key =timepoint, value = value, fill = 0) %>% 
    mutate(d = post - pre) %>% 
    dplyr::select(-pre, -post) %>%
    group_by(gene) %>%
    do(tidy(t.test(d ~ group, data = .))) -> sel_dif_avg_ttests
  
  sel_dif_avg_ttests %>%
    filter(!str_detect(gene, "^RP[LS]")) %>% 
    mutate(gene = as.character(gene)) %>% 
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
    filter(p.adj <= 0.05) %>%   
    arrange(statistic)  %>%
    #filter(estimate >= 0.01 | estimate1 >= 0.01) %>% 
    mutate(dir = if_else(statistic > 0, "high","low")) %>% 
    group_by(dir) %>% 
    filter(row_number() > max(row_number()) - 20 | row_number() <= 20)  %>% 
    arrange(dir, p.value) %>%  
    pull(gene) -> sel_dif_avg_top
  
  
  
  reshape2::melt(all_genes_avg@assays$SCT@data) %>% 
    rename("gene" = Var1, "sample" = Var2) %>% 
    mutate(patient = str_sub(sample, 1, 5), timepoint = str_sub(sample, 7)) %>% 
    left_join(selection_df) %>% 
    filter(group != "unclassified") %>% 
    mutate(group = as.character(group)) %>%
    dplyr::select(-sample) %>%
    spread(key =timepoint, value = value, fill = 0) %>% 
    mutate(d = post - pre) %>% 
    dplyr::select(-pre, -post, -group) %>%
    spread(key = patient, value = d, fill = 0) %>% 
    filter(gene %in%  sel_dif_avg_top) ->
    sel_tum_dif_df
  
  
  
}





###########
###gsea###
##########


sel_tum_pre_avg_ttests %>%
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  mutate(gene = as.character(gene)) %>% 
  #mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  #filter(p.adj <= 0.05) %>% 
  filter(estimate1 >= 0.25 | estimate2 >= 0.25) %>% 
  arrange(statistic) -> temp
temp_gsea_order = pull(temp,statistic) %>% set_names(temp$gene)

m_df = msigdbr(species = "Homo sapiens", category = "H")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
sel_tum_pre_gsea_H = fgsea(pathways = m_list,stats = temp_gsea_order, nperm = 1000)
sel_tum_pre_gsea_H %>% arrange(-NES) %>% 
  filter(row_number() > max(row_number()) - 2 | row_number() <= 1) %>% 
  dplyr::select(pathway, NES) ->
  sel_tum_pre_gsea_H_top
  


m_df = msigdbr(species = "Homo sapiens", category = "C2")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
sel_tum_pre_gsea_C2 = fgsea(pathways = m_list,stats = temp_gsea_order, nperm = 1000)
sel_tum_pre_gsea_C2 %>% 
  filter(str_detect(pathway, "^WP|BIOCARTA|REACTOME")) %>% 
  #filter(str_detect(pathway, "^REACTOME")) %>% 
  #filter(str_detect(pathway, "^WP")) %>% 
  #filter(str_detect(pathway, "^KEGG")) %>% 
  arrange(-NES) %>% 
  filter(row_number() > max(row_number()) - 10 | row_number()  <= 10) %>% 
  filter(row_number() %in% c(5,7,8,18,19)) %>% 
  dplyr::select(pathway, NES) ->
  sel_tum_pre_gsea_C2_top

rbind(sel_tum_pre_gsea_C2_top, sel_tum_pre_gsea_H_top) %>% 
  mutate(pathway = str_replace(pathway, "_PATHWAY_PUBERTY_STAGE_2_OF_4","")) %>% 
  mutate(pathway = str_replace(pathway, "INTERFERON","IFN")) %>% 
  arrange(NES) %>% 
  ggplot() + aes(y = reorder(pathway, NES), x = NES, fill = NES) + 
  geom_bar(stat="identity") + theme_real_minimal() + 
  theme(axis.title = element_blank()) +
  #scale_fill_gradient2(low = "#510270",mid = "#db5b6c",  high = "#ffa16e")
scale_fill_gradient(low = "#007477",  high = "#ff8f64") -> sel_tum_pre_gsea_plot

ggsave(sel_tum_pre_gsea_plot,filename =  "~/projects/Breast_cancer_radiation/fig_pieces_06/sel_tum_pre_gsea_plot.pdf", width=6, height = 4)









sel_tum_post_avg_ttests %>%
  filter(!str_detect(gene, "^RP[LS]")) %>% 
  mutate(gene = as.character(gene)) %>% 
  #mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  #filter(p.adj <= 0.05) %>% 
  filter(estimate1 >= 0.25 | estimate2 >= 0.25) %>% 
  arrange(statistic) -> temp
temp_gsea_order = pull(temp,statistic) %>% set_names(temp$gene)

m_df = msigdbr(species = "Homo sapiens", category = "H")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
sel_tum_post_gsea_H = fgsea(pathways = m_list,stats = temp_gsea_order, nperm = 1000)
sel_tum_post_gsea_H %>% arrange(-NES) %>% 
  
  
  filter(row_number() > max(row_number()) - 2 | row_number() <= 1) %>% 
  dplyr::select(pathway, NES) ->
  sel_tum_post_gsea_H_top

m_df = msigdbr(species = "Homo sapiens", category = "C2")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
sel_tum_post_gsea_C2 = fgsea(pathways = m_list,stats = temp_gsea_order, nperm = 1000)
sel_tum_post_gsea_C2 %>% 
  filter(str_detect(pathway, "^WP|BIOCARTA|REACTOME")) %>% 
  #filter(str_detect(pathway, "^REACTOME")) %>% 
  #filter(str_detect(pathway, "^WP")) %>% 
  #filter(str_detect(pathway, "^KEGG")) %>% 
  arrange(-NES) %>% 
  filter(row_number() > max(row_number()) - 10 | row_number()  <= 10) %>% 
  filter(row_number() %in% c(5,7,8,18,19)) %>% 
  dplyr::select(pathway, NES) ->
  sel_tum_post_gsea_C2_top





##########################
######ER/AR/PR/HER2#######
##########################

irds_order2 = patients_key[irds_order] %>% setNames(NULL)
selection_group = data.frame(group = c(rep("HGS",4),rep("LGS",4),rep(NA,3)), patient.2 = paste0("P",str_pad(1:11, side="left", width=2, pad = "0")))
Idents(sct)
AverageExpression(so, add.ident = "timepoint",return.seurat = T)

Idents(SCTall_fil2) = "patient.2"
Idents(SCTall_fil2_epi) = "patient.2"
recetor_tum_average = AverageExpression(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "BCR")],
                                      add.ident = "timepoint",return.seurat = T, features = c("ESR1","AR","PGR","ERBB2"))


data.frame(recetor_tum_average@assays$SCT@data) %>% 
  rownames_to_column(var = "gene") %>% 
  gather(-gene, key = "sample", value = "value") %>% 
  mutate(patient.2 = as.character(str_sub(sample, 1,3)), timepoint = str_sub(sample, 5)) %>% 
  dplyr::select(-sample) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  left_join(selection_group) %>% 
  spread(key = gene, value = value)  ->
  recetor_tum_average_df

recetor_tum_average_df_pre = recetor_tum_average_df %>%  filter(timepoint == "pre") %>% column_to_rownames("patient.2")
recetor_tum_average_df_post = recetor_tum_average_df %>%  filter(timepoint == "post") %>% column_to_rownames("patient.2")
recetor_tum_average_mat_pre = t(scale(as.matrix(recetor_tum_average_df_pre[,c(-1,-2)])))
recetor_tum_average_mat_post = t(scale(as.matrix(recetor_tum_average_df_post[,c(-1,-2)])))


group_colors2 = group_colors %>% setNames(c("HGS","LGS","unclassified"))

#library(circlize)
col_fun = colorRamp2(c(-2, 0,2), c("steelblue","white", "firebrick"))
ca = columnAnnotation(timepoint = recetor_tum_average_df_pre$timepoint,
                       group = recetor_tum_average_df_pre$group,
                      annotation_name_side = "left",
                       col = list(timepoint = treatment_colors, 
                                  group = group_colors2[1:2]))
Heatmap(recetor_tum_average_mat_pre, column_split = recetor_tum_average_df_pre$group,
        show_row_dend = F, show_column_dend = F, 
        row_names_side = "left", column_names_side = "bottom",
        column_order = paste0("P0",1:8), col = col_fun,
        row_order = c("ESR1","AR","PGR","ERBB2"), top_annotation =  ca) -> receptor_group_heatmap_pre
pdf("~/projects/Breast_cancer_radiation/fig_pieces_06/receptor_group_heatmap_pre.pdf", width = 4, height = 2)
draw(receptor_group_heatmap_pre)
dev.off()

col_fun = colorRamp2(c(-2, 0,2), c("steelblue","white", "firebrick"))
ca = columnAnnotation(timepoint = recetor_tum_average_df_post$timepoint,
                      group = recetor_tum_average_df_post$group,
                      annotation_name_side = "left",
                      col = list(timepoint = treatment_colors, 
                                 group = group_colors2[1:2]))
Heatmap(recetor_tum_average_mat_post, column_split = recetor_tum_average_df_post$group,
        show_row_dend = F, show_column_dend = F, 
        row_names_side = "left", column_names_side = "bottom",
        column_order = paste0("P0",1:8), col = col_fun,
        row_order = c("ESR1","AR","PGR","ERBB2"), top_annotation =  ca) -> receptor_group_heatmap_post
pdf("~/projects/Breast_cancer_radiation/fig_pieces_06/receptor_group_heatmap_post.pdf", width = 4, height = 2)
draw(receptor_group_heatmap_post)
dev.off()



recetor_tum_average_df %>% 
  mutate(group = if_else(patient.2 %in% c("P01","P02","P03","P04"), "HGS","LGS")) %>% 
  #filter(timepoint == "post") %>% 
  gather(-patient.2, -timepoint, -group, key = receptor, value = m) %>% 
  group_by(timepoint, receptor) %>%
  summarise(ttest = t.test(m ~ group)$p.value)


data.frame(patient = SCTall_fil2_epi$patient.2, timepoint = SCTall_fil2_epi$timepoint, s = colSums(SCTall_fil2_epi@assays$SCT@data[c("ESR1","ERBB2","PGR","AR"),])) %>% 
  mutate(e = if_else(s >0, "e","ne")) %>% 
  group_by(patient, timepoint, e) %>% 
  summarise(n = n()) %>% 
  mutate(p = n/sum(n)) %>% 
  filter(e == "e") %>% 
  ungroup() %>% 
  mutate(group = if_else(patient %in% c("P01","P02","P03","P04"), "HGS","LGS")) %>% 
  group_by(timepoint) %>% 
  summarise(ttest = t.test(p ~ group, var.equal = TRUE)$p.value)




###################
##irds boxplots####
###################

SCTall_fil2@meta.data %>% 
  filter(medium_clusters == "Tumor Epithelial") %>% 
  filter(timepoint == "pre") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  ggplot() + aes(x = patient.2, y = irds.1, color = group) + geom_boxplot(fill = NA) +
  scale_x_discrete(limits = irds_order2[c(1:3,6:10)]) + theme_real_minimal() +
  ylab("IFN-related DNA damage resistance score") +
  theme(axis.title.x = element_blank()) + 
  theme(legend.title = element_blank()) +
  theme(axis.line.y = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5)) +
  scale_color_manual(values = group_colors[1:2]) -> irds_boxplots 
ggsave(irds_boxplots, filename = "~/projects/Breast_cancer_radiation/fig_pieces_06/irds_boxplots.pdf", width = 6, height = 3.5, device = cairo_pdf)
  
  

#####################
####IRDS vln plots###
#####################

irds_vln_plots = VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial" & SCTall_fil2$timepoint == "pre" & SCTall_fil2$patient.2 %in% paste0("P0",1:8)], 
                         group.by = "patient.2", features = c(irds_genes_top,"ESR1"), same.y.lims = T, combine = F, pt.size = 0, cols = c("#DD7861","#DD7861","#DD7861","#DD7861","#3A8A8A","#3A8A8A","#3A8A8A","#3A8A8A"))
names(irds_vln_plots) = c(irds_genes_top,"ESR1")

irds_vln_plots = lapply(irds_vln_plots, function(x) return(x + scale_y_discrete(expand = c(0,0)) + 
                                                             theme(axis.ticks = element_blank(), 
                                                                   legend.position = "none", 
                                                                   axis.line.y = element_blank(),
                                                                   axis.title = element_blank(), 
                                                                   axis.text.x = element_blank())))

plot_grid(plotlist = irds_vln_plots, ncol=2) -> irds_vln_plots_combined
ggsave(irds_vln_plots_combined, filename = "~/projects/Breast_cancer_radiation/fig_pieces_06/irds_vln_plots_combined_raw.pdf", width = 4, height = 5, device = cairo_pdf)



##################
###checkpoints####
##################

checkpoint_tcells = c("PVRIG","CD200R1","CD244","TMIGD2","CD96","CD226","BTLA","TNFRSF18","CTLA4","TIGIT","PDCD1","HAVCR2","LAG3","CD28","CD40LG","TNFRSF4","TNFRSF9","ICOS","CD27")
checkpoint_tcells_pos = c("CD28","TNFRSF4","TNFRSF9","CD27","ICOS","CD226","CD244","TNFRSF18","CD40LG")
checkpoint_tcells_neg = c("CD96","BTLA","HAVCR2","CD200R1","LAG3","CTLA4","TMIGD2","PDCD1","TIGIT","PVRIG")

checkpoint_genes = c("TNFRSF14","TNFSF18","CD80","CD86","PVR","CD274","HHLA2",
                     "PDCD1LG2","LGALS9","CD40","TNFSF4","TNFSF9","ICOSLG",
                     "CD70","NECTIN2","CD276","VSIR","BTN3A1","CD48","CD200","VTCN1")

checkpoint_neg = c("TNFRSF14","TNFSF18","CD80","CD86","PVR","CD274","HHLA2",
                   "PDCD1LG2","LGALS9","CD40",
                   "NECTIN2","CD276","CD200","VTCN1")

checkpoint_pos = c("TNFSF4","TNFSF9","ICOSLG","CD48","CD70")


tumor_temp = AddModuleScore(SCTall_fil2[,SCTall_fil2$medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")], features = list(checkpoint_neg), search = T, name = "tum.ckpt.neg.")
tumor_temp = AddModuleScore(tumor_temp, features = list(checkpoint_neg), search = T, name = "tum.ckpt.pos.")
tumor_temp = tumor_temp@meta.data



set.seed(77)
tumor_temp %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "pre") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  group_by(patient.2) %>% 
  slice_sample(n=100) %>% 
  mutate(group = if_else(patient.2 %in% paste0("P0",1:4),"HGS","LGS")) %>% 
  ungroup() %>% 
  mutate(pval = wilcox.test(tum.ckpt.neg.1 ~ group, var.equal = TRUE)$p.value) %>% 
  #pull(pval)
  ggplot() + aes(x = group, y = tum.ckpt.neg.1, color = group) + 
  geom_boxplot(outlier.color = NA) + scale_color_manual(values = group_colors2[1:2]) +
  theme_real_minimal() + ggtitle("PRE") + theme(legend.position = "none") +
  theme(axis.title = element_blank()) + theme(text = element_text(size = 16)) +
  theme(axis.line.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(-0.1, 0.2)) -> ckpt_tum_pre

#pre-tumor,neg, wilcox.test p = 2.98217e-21

set.seed(77)
tumor_temp %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "pre") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  mutate(group = if_else(patient.2 %in% paste0("P0",1:4),"HGS","LGS")) %>% 
  group_by(group,patient.2) %>%
  summarise(m = median(tum.ckpt.neg.1)) %>% 
  ungroup() %>% 
  mutate(pval = t.test(m ~ group, var.equal = TRUE)$p.value)

#pre-tumor,neg, t.test p = 0.206


set.seed(77)
tumor_temp %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "post") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  group_by(patient.2) %>% 
  slice_sample(n=100) %>% 
  mutate(group = if_else(patient.2 %in% paste0("P0",1:4),"HGS","LGS")) %>% 
  ungroup() %>% 
  mutate(pval = wilcox.test(tum.ckpt.neg.1 ~ group)$p.value) %>% 
  #pull(pval)
  ggplot() + aes(x = group, y = tum.ckpt.neg.1, color = group) + 
  geom_boxplot(outlier.color = NA) + scale_color_manual(values = group_colors2[1:2]) +
  theme_real_minimal() + ggtitle("POST") + theme(legend.position = "none") +
  theme(axis.title = element_blank()) + theme(text = element_text(size = 16)) +
  theme(axis.line.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(-0.1, 0.2)) -> ckpt_tum_post

#post-tumor, t.test p = 1.679496e-09


set.seed(77)
tumor_temp %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "post") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  mutate(group = if_else(patient.2 %in% paste0("P0",1:4),"HGS","LGS")) %>% 
  group_by(group,patient.2) %>%
  summarise(m = median(tum.ckpt.neg.1)) %>% 
  ungroup() %>% 
  mutate(pval = t.test(m ~ group, var.equal = TRUE)$p.value)

#post-tumor, t.test p = 0.136



tcells_temp = AddModuleScore(SCTall_fil2[,SCTall_fil2$major_clusters %in% c("T - Cells","NK Cells")], features = list(checkpoint_tcells_neg), search = T, name = "tcell.ckpt.neg.", nbin = 20)
tcells_temp = AddModuleScore(tcells_temp, features = list(checkpoint_tcells_pos), search = T, name = "tcell.ckpt.pos.", nbin = 20)
tcells_temp = tcells_temp@meta.data


set.seed(77)
tcells_temp %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "pre") %>% 
  filter(patient.2 %in% paste0("P0",c(2,3,5:8))) %>% 
  group_by(patient.2) %>% 
  slice_sample(n=100, replace = F) %>% 
  mutate(group = if_else(patient.2 %in% paste0("P0",1:4),"HGS","LGS")) %>% 
  ungroup() %>% 
  mutate(pval = wilcox.test(tcell.ckpt.neg.1 ~ group)$p.value) %>%  
  #pull(pval)
  ggplot() + aes(x = group, y = tcell.ckpt.neg.1, color=group) + 
  geom_boxplot(outlier.color = NA) + scale_color_manual(values = group_colors2[1:2]) +
  theme_real_minimal() + ggtitle("PRE") + theme(legend.position = "none") +
  theme(axis.title = element_blank()) + theme(text = element_text(size = 16)) +
  theme(axis.line.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(-0.3, 0.3)) -> ckpt_tcell_pre

  #pre-tcells, wilcox.test p = 0.1150403
  
  
  
set.seed(77)
tcells_temp %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "pre") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  mutate(group = if_else(patient.2 %in% paste0("P0",1:4),"HGS","LGS")) %>% 
  group_by(group,patient.2) %>%
  summarise(m = median(tcell.ckpt.neg.1)) %>% 
  ungroup() %>% 
  mutate(pval = t.test(m ~ group, var.equal = TRUE)$p.value)

#pre-tcells, t.test p = 0.855




set.seed(77)
tcells_temp %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "post") %>% 
  filter(patient.2 %in% paste0("P0",c(1,2,5:8))) %>% 
  group_by(patient.2) %>% 
  slice_sample(n=100, replace = F) %>% 
  mutate(group = if_else(patient.2 %in% paste0("P0",1:4),"HGS","LGS")) %>% 
  ungroup() %>% 
  mutate(pval = t.test(tcell.ckpt.neg.1 ~ group, var.equal = TRUE)$p.value) %>%
  pull(pval) %>% 
  ggplot() + aes(x = group, y = tcell.ckpt.neg.1, color=group) + 
  geom_boxplot(outlier.color = NA) + scale_color_manual(values = group_colors2[1:2]) +
  theme_real_minimal() + ggtitle("POST") + theme(legend.position = "none") +
  theme(axis.title = element_blank()) + theme(text = element_text(size = 16)) +
  theme(axis.line.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(-0.3, 0.3)) -> ckpt_tcell_post


#post-tcells, wilcox.test p = 1.896867e-11



set.seed(77)
tcells_temp %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial","Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "post") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  mutate(group = if_else(patient.2 %in% paste0("P0",1:4),"HGS","LGS")) %>% 
  group_by(group,patient.2) %>%
  summarise(m = median(tcell.ckpt.neg.1)) %>% 
  ungroup() %>% 
  mutate(pval = t.test(m ~ group, var.equal = TRUE)$p.value)

#post-tcells, t.test p = 0.0745


plot_grid(ckpt_tum_pre,ckpt_tum_post, ckpt_tcell_pre, ckpt_tcell_post, nrow=1) -> ckpt_tum_tcells 

ggsave(ckpt_tum_tcells, filename = "~/projects/Breast_cancer_radiation/fig_pieces_06/ckpt_boxplot.pdf", device = cairo_pdf, width = 6, height = 3)



##########
major_percent_group %>% 
  left_join(data.frame(patient = names(patients_key), patient.2 = patients_key)) %>% 
  mutate(group = if_else(patient.2 %in% c("P01","P02","P03","P04"), "HGS", "undetermined")) %>% 
  mutate(group = if_else(patient.2 %in% c("P05","P06","P07","P08"), "LGS", group)) %>% 
  filter(timepoint == "pre") %>%
  ungroup() %>% 
  select(-patient, -timepoint) %>% 
  gather(-group, -patient.2, key = "clus", value = "p") %>% 
  ggplot() + aes(x = clus, y = patient.2, size = p, color = group) + 
  geom_point() + scale_color_manual(values = group_colors2[1:2]) +
  theme_minimal() + scale_size_continuous(range = c(0,6)) +
  theme(axis.title = element_blank(), panel.grid = element_blank()) -> temp
ggsave(temp, filename = "~/projects/Breast_cancer_radiation/sup_figure_pieces/immune_frac_pre.pdf", device = cairo_pdf, width = 4, height = 4)
  
  
  

  
