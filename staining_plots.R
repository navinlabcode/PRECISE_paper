library(tidyverse)
library(corrr)
library(ComplexHeatmap)

setwd("/volumes/USR1/aschalck/TENX/analysis/BCR/BCR_4/")

# SCTall_fil2@meta.data %>% 
#   select(patient, timepoint, medium_clusters) %>% 
#   group_by(patient, timepoint, medium_clusters) %>% 
#   summarise(n = n()) %>% 
#   mutate(percent = n/sum(n)) %>% 
#   select(-n) %>% 
#   mutate(platform = "single cell") %>% 
#   rename("marker" = "medium_clusters") %>% 
#   mutate(tissue = "TME") ->
#   medium_cluster_percent

#write_csv(medium_cluster_percent, "./data/medium_clusters_percent.csv")

# SCTall_fil2@meta.data %>% 
#   select(patient, timepoint,major_clusters, medium_clusters) %>% 
#   mutate(marker = if_else(medium_clusters %in% c("Tumor Epithelial", "Cycling Epithelial"),"Tumor","Stroma")) %>% 
#   mutate(marker = if_else(major_clusters %in% c("B - Cells","Myeloid","NK Cells","T - Cells"),"Immune",marker)) %>% 
#   select(patient, timepoint, marker) %>% 
#   group_by(patient, timepoint, marker) %>% 
#   summarise(n = n()) %>% 
#   mutate(percent = n/sum(n)) %>% 
#   select(-n) %>% 
#   mutate(platform = "single cell") %>% 
#   mutate(tissue = "TME") ->
#   mega_cluster_percent

#write_csv(mega_cluster_percent, "./data/mega_clusters_percent.csv")


# 
# staining_raw = read_csv("./data/staining_data.csv")
# 
# staining_raw %>%
#    mutate(patient = str_replace(patient, "BR","BCR")) %>%
#    mutate(patient = if_else(patient == "N/A","BCR00",patient)) %>%
#    mutate(timepoint = if_else(str_detect(timepoint, "post|POST|Post"),"post",timepoint)) %>%
#    mutate(timepoint = if_else(str_detect(timepoint, "pre|PRE|Pre"),"pre",timepoint)) %>%
#    rename("percent" = "percentage", "density" = "density(n/mm)") %>%
#    mutate(percent = as.numeric(percent)) ->
#    staining_clean
# 
#write_csv(staining_clean, "./data/staining_data_clean.csv")


IHC_1 = read_csv("./data/staining_data_clean.csv")
IHC_2 = read_csv("./data/ki67.csv")
IHC_3 = read_csv("./data/tunel.csv")

sc_1 = read_csv("./data/medium_clusters_percent.csv") %>%  complete(patient, timepoint, marker, fill = list(percent = 0, platform = "single cell", tissue = "TME"))
sc_2 = read_csv("./data/mega_clusters_percent.csv") %>%  complete(patient, timepoint, marker, fill = list(percent = 0, platform = "single cell", tissue = "TME"))


bind_rows(IHC_1, IHC_2,IHC_3,sc_1,sc_2) %>% 
  gather(density, percent, key = "measure", value = "value") %>% 
  filter(timepoint %in% c("pre","post")) %>% 
  group_by(patient,marker,measure,tissue,timepoint,platform) %>% 
  filter(!is.na(value)) %>% 
  summarise(value = as.numeric(mean(value))) %>% 
  ungroup() %>% 
  mutate(platform = str_replace(platform, "IHC Image analysis\\(HALO\\)", "IHC-HALO")) %>% 
  mutate(platform = str_replace(platform, "Manual Microscopy", "IHC-manual")) %>% 
  unite(variable, c("marker","measure","tissue","timepoint","platform")) %>%
  mutate(variable = if_else(str_detect(variable, "single cell"), str_replace(variable, "_percent_TME",""),variable)) %>% 
  mutate(variable = str_replace(variable, " ","-")) %>% 
  mutate(variable = str_replace(variable, "single cell","single-cell")) %>%
  spread(key = variable, value = value) -> 
  IHC_sc_df

write_csv(IHC_sc_df,"./data/IHC_sc_df.csv")

IHC_sc_df %>% 
  column_to_rownames("patient") -> IHC_sc_df

IHC_sc_rho_mat = matrix(NA, nrow = length(names(IHC_sc_df)), ncol = length(names(IHC_sc_df)), dimnames = list(names(IHC_sc_df),names(IHC_sc_df)))
IHC_sc_p_mat = matrix(NA, nrow = length(names(IHC_sc_df)), ncol = length(names(IHC_sc_df)), dimnames = list(names(IHC_sc_df),names(IHC_sc_df)))
for(i in names(IHC_sc_df)){
  print(i)
  for(j in names(IHC_sc_df)){
    #print(j)
    t = cor.test(IHC_sc_df[,i], IHC_sc_df[,j], method = "spearman",use = "pairwise.complete.obs")
    IHC_sc_rho_mat[i,j] = t$estimate[[1]]
    IHC_sc_p_mat[i,j] = t$p.value
  }
}

IHC_sc_padj_mat = matrix(p.adjust(IHC_sc_p_mat, method = "BH"), nrow = length(names(IHC_sc_df)), ncol = length(names(IHC_sc_df)), dimnames = list(names(IHC_sc_df),names(IHC_sc_df)))
#IHC_sc_padj_mat = IHC_sc_p_mat
IHC_sc_star_mat = IHC_sc_padj_mat
IHC_sc_star_mat[IHC_sc_star_mat <= 0.05] = "*"
IHC_sc_star_mat[IHC_sc_star_mat != "*"] = ""

col_fun = circlize::colorRamp2(c(-1, -0.25, 0, 0.25, 1), c("steelblue", "white", "white", "white", "firebrick"))
temp = Heatmap(IHC_sc_rho_mat, col = col_fun, row_names_side = "left", show_row_dend = F, show_column_dend = F,
               cell_fun = function(j,i,x, y, width, height, fill) {
                 grid.text(IHC_sc_star_mat[i,j], x, y, gp = gpar(fontsize = 16, col="white"))
               })

pdf(file = "~/projects/Breast_cancer_radiation/sup_figure_pieces/IHC_cor_heatmap_all.pdf", useDingbats = F, width = 12, height = 12)
draw(temp, padding = unit(c(30, 30, 2, 2), "mm"))
dev.off()


IHC_sc_p_mat2 = IHC_sc_p_mat[!str_detect(row.names(IHC_sc_p_mat), "single-cell"),str_detect(row.names(IHC_sc_p_mat), "single-cell")]
IHC_sc_rho_mat2 = IHC_sc_rho_mat[!str_detect(row.names(IHC_sc_rho_mat), "single-cell"),str_detect(row.names(IHC_sc_rho_mat), "single-cell")]

IHC_sc_padj_mat2 = matrix(p.adjust(IHC_sc_p_mat2, method = "BH"), nrow = nrow(IHC_sc_p_mat2), ncol = ncol(IHC_sc_p_mat2), dimnames = list(rownames(IHC_sc_p_mat2),colnames(IHC_sc_p_mat2)))
#IHC_sc_padj_mat2 = IHC_sc_p_mat2
IHC_sc_star_mat2 = IHC_sc_padj_mat2
IHC_sc_star_mat2[IHC_sc_star_mat2 <= 0.05] = "*"
IHC_sc_star_mat2[IHC_sc_star_mat2 != "*"] = ""
IHC_sc_star_mat2[is.na(IHC_sc_star_mat2)] = ""

col_fun = circlize::colorRamp2(c(-1, -0.25, 0, 0.25, 1), c("steelblue", "white", "white", "white", "firebrick"))
temp = Heatmap(IHC_sc_rho_mat2, col = col_fun, row_names_side = "left", show_row_dend = F, show_column_dend = F,
               cell_fun = function(j,i,x, y, width, height, fill) {
                 grid.text(IHC_sc_star_mat2[i,j], x, y, gp = gpar(fontsize = 16, col="white"))
               })

pdf(file = "~/projects/Breast_cancer_radiation/sup_figure_pieces/IHC_cor_heatmap_sc_IHC.pdf", useDingbats = F, width = 10, height = 7)
draw(temp, padding = unit(c(30, 30, 2, 2), "mm"))
dev.off()



IHC_sc_p_mat3 = IHC_sc_p_mat[!str_detect(row.names(IHC_sc_p_mat), "single-cell"),!str_detect(row.names(IHC_sc_p_mat), "single-cell")]
IHC_sc_rho_mat3 = IHC_sc_rho_mat[!str_detect(row.names(IHC_sc_rho_mat), "single-cell"),!str_detect(row.names(IHC_sc_rho_mat), "single-cell")]

IHC_sc_padj_mat3 = matrix(p.adjust(IHC_sc_p_mat3, method = "BH"), nrow = nrow(IHC_sc_p_mat3), ncol = ncol(IHC_sc_p_mat3), dimnames = list(rownames(IHC_sc_p_mat3),colnames(IHC_sc_p_mat3)))
#IHC_sc_padj_mat3 = IHC_sc_p_mat3
IHC_sc_star_mat3 = IHC_sc_padj_mat3
IHC_sc_star_mat3[IHC_sc_star_mat3 <= 0.05] = "*"
IHC_sc_star_mat3[IHC_sc_star_mat3 != "*"] = ""

col_fun = circlize::colorRamp2(c(-1, -0.25, 0, 0.25, 1), c("steelblue", "white", "white", "white", "firebrick"))
temp = Heatmap(IHC_sc_rho_mat3, col = col_fun, row_names_side = "left", show_row_dend = F, show_column_dend = F,
               cell_fun = function(j,i,x, y, width, height, fill) {
                 grid.text(IHC_sc_star_mat3[i,j], x, y, gp = gpar(fontsize = 16, col="white"))
               })

pdf(file = "~/projects/Breast_cancer_radiation/sup_figure_pieces/IHC_cor_heatmap_IHC.pdf", useDingbats = F, width = 8, height = 7)
draw(temp, padding = unit(c(30, 30, 2, 2), "mm"))
dev.off()




bind_rows(IHC_1,IHC_2,IHC_3) %>% 
  gather(density, percent, key = "measure", value = "value") %>% 
  filter(timepoint %in% c("pre","post")) %>% 
  group_by(patient,marker,measure,tissue,timepoint,platform) %>% 
  filter(!is.na(value)) %>% 
  summarise(value = as.numeric(mean(value))) %>% 
  ungroup() %>% 
  mutate(platform = str_replace(platform, "IHC Image analysis\\(HALO\\)", "IHC-HALO")) %>% 
  mutate(platform = str_replace(platform, "Manual Microscopy", "IHC-manual")) %>% 
  unite(variable, c("marker","measure","tissue","platform")) %>%
  mutate(variable = str_replace(variable, " ","-")) ->
  IHC_prepost_df

  
for(i in unique(IHC_prepost_df$variable)){
  
  IHC_prepost_df %>% 
    filter(variable == i) %>% 
    pull(timepoint) %>% 
    unique() %>% 
    length() -> n
  if(n == 2){
    print(i)
    IHC_prepost_df %>% 
      filter(variable == i) %>% 
      spread(key = timepoint, value = value) ->
      temp_df
    p = wilcox.test(temp_df$pre, temp_df$post, paired = T)$p.value
    #p = round(p, digits = 2)
    print(i)
    print(p)
    IHC_prepost_df %>% 
      filter(variable == i) %>% 
      ggplot() + aes(x = timepoint, y = value, group = patient) + 
      scale_x_discrete(limits = c("pre","post"),expand = c(0.1, 0)) +
      geom_point() + geom_line() + ylab(i) + 
      annotate("text", label = paste0("p = ",p), x = 1.5, y = Inf, vjust = 2, size=6) + 
      theme(panel.background = element_blank(), panel.grid = element_blank()) +
      theme(axis.title.x = element_blank(), axis.line.y = element_line(size = 0.5)) +
      theme(axis.ticks.x = element_blank(), text = element_text(size = 16)) -> temp_plot
    #ggsave(temp_plot, filename = paste0("~/projects/Breast_cancer_radiation/sup_figure_pieces/",i,"_dotplot.pdf"), height = 4, width = 2.5, device = cairo_pdf)
    
  }
  
  
}






IHC_sc_df %>% 
  left_join(clinical_data %>% select(patient, patient.2, cohort)) %>% 
  left_join(selection_group) %>% 
  mutate(dif = CD8_density_TME_post_IMT-CD8_density_TME_pre_IMT) %>% 
  select(dif, cohort) %>% 
  summarise(pval = wilcox.test(dif ~ cohort)$p.value) 

IHC_sc_df %>% 
  left_join(clinical_data %>% select(patient, patient.2, cohort)) %>% 
  left_join(selection_group) %>% 
  mutate(dif = CD3_density_TME_post_IMT-CD3_density_TME_pre_IMT) %>% 
  select(dif, cohort) %>% 
  summarise(pval = wilcox.test(dif ~ cohort)$p.value) 




IHC_sc_df %>% 
  left_join(clinical_data %>% select(patient, patient.2, cohort)) %>% 
  left_join(selection_group) %>% 
  mutate(dif = CD8_density_TME_post_IMT-CD8_density_TME_pre_IMT) %>% 
  ggplot() + aes(x = cohort, y = dif) + 
  geom_boxplot(outlier.colour = NA, fill = NA) +
  geom_point(position = position_jitter(width = 0.2, height = 0), aes(color = group)) + 
  scale_color_manual(values = group_colors2) + theme_minimal() + 
  ggtitle("dif CD8 density")


IHC_sc_df %>% 
  left_join(clinical_data %>% select(patient, patient.2, cohort)) %>% 
  left_join(selection_group) %>% 
  mutate(dif = CD3_density_TME_post_IMT-CD3_density_TME_pre_IMT) %>% 
  ggplot() + aes(x = cohort, y = dif) + 
  geom_boxplot(outlier.colour = NA, fill = NA) +
  geom_point(position = position_jitter(width = 0.2, height = 0), aes(color = group)) + 
  scale_color_manual(values = group_colors2) + theme_minimal() + 
  ggtitle("dif CD3 density")



IHC_sc_df %>% 
  left_join(clinical_data %>% select(patient, patient.2, cohort)) %>% 
  left_join(selection_group) %>% 
  ggplot() + aes(x = cohort, y = CD8_density_TME_post_IMT) + 
  geom_boxplot(outlier.colour = NA, fill = NA) +
  geom_point(position = position_jitter(width = 0.2, height = 0), aes(color = group)) + 
  scale_color_manual(values = group_colors2) + theme_minimal()

IHC_sc_df %>% 
  left_join(clinical_data %>% select(patient, patient.2, cohort)) %>% 
  left_join(selection_group) %>% 
  ggplot() + aes(x = cohort, y = CD3_density_TME_post_IMT) + 
  geom_boxplot(outlier.colour = NA, fill = NA) +
  geom_point(position = position_jitter(width = 0.2, height = 0), aes(color = group)) + 
  scale_color_manual(values = group_colors2) + theme_minimal()


IHC_sc_df %>% 
  left_join(clinical_data %>% select(patient, patient.2, cohort)) %>% 
  left_join(selection_group) %>% 
  summarise(pval = wilcox.test(CD8_density_TME_post_IMT ~ cohort)$p.value) 

IHC_sc_df %>% 
  left_join(clinical_data %>% select(patient, patient.2, cohort)) %>% 
  left_join(selection_group) %>% 
  summarise(pval = wilcox.test(CD3_density_TME_post_IMT ~ cohort)$p.value) 




