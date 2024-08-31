

#SCTall_fil2

radio_genes_df = read_csv("~/code/resources/gene_lists/radio_genes_table.csv", col_select = 1:6)


radio_genes_df %>% 
  filter(!is.na(gene)) %>% 
  filter(gene %in% row.names(SCTall_fil2)) ->
  radio_genes_df
  #group_by(signature_name, coefficient_direction)


Idents(SCTall_fil2) = "medium_clusters"

radio_genes_pre_average = AverageExpression(SCTall_fil2[,SCTall_fil2$timepoint=="pre"],
                                            add.ident = "patient", 
                                        features = unique(radio_genes_df$gene))

temp1 = SCTall_fil2[,SCTall_fil2$timepoint=="pre"]$medium_clusters
temp2 = SCTall_fil2[,SCTall_fil2$timepoint=="pre"]$patient
table(temp1,temp2) -> temp

temp %>% 
  reshape2::melt(variable.name = c("cluster","patient"), value.name = "n") %>% 
  filter(n >= 25) %>% 
  rename("cluster" = temp1, "patient" = temp2) %>% 
  mutate("cluster_patient" = paste0(cluster, "_", patient)) %>% 
  pull(cluster_patient) ->
  cluster_patient_25
  

radio_genes_pre_average$SCT %>% 
  rownames_to_column(var = "gene") %>% 
  reshape2::melt(value.name = "expr", variable.name = "cluster_patient") %>% 
  mutate(cluster = str_split(cluster_patient, pattern = "_", simplify = T)[,1]) %>% 
  mutate(patient = str_split(cluster_patient, pattern = "_", simplify = T)[,2]) %>% 
  filter(cluster_patient %in% cluster_patient_25) %>% 
  group_by(gene, cluster) %>% 
  summarise(m = mean(expr)) %>% 
  filter(m >= 0.1) %>% 
  group_by(cluster) %>% 
  summarise(genes = list(unique(gene))) %>% 
  pull(genes, name = cluster) -> lt

ltm = list_to_matrix(lt)

m1 = make_comb_mat(lt)
m2 = m1[comb_size(m1) >= 3]

#pdf(file = "./plots/plots13/rad_sig_genes_upset.pdf")
#UpSet(m2, comb_order = rev(order(comb_size(m2))))
#dev.off()

tum_genes = row.names(ltm[rowSums(ltm[,c("Tumor Epithelial","Cycling Epithelial")]) > 0,])

radio_genes_df %>% 
  mutate(tum = if_else(gene %in% tum_genes, "yes","no")) %>% 
  group_by(signature_name, coefficient_direction, tum) %>% 
  summarise(n = n()) %>% 
  spread(key = tum, value = n, fill = 0) %>% 
  mutate(p = yes/(yes+no)) %>% 
  arrange(-p) %>% 
  filter(!is.na(coefficient_direction)) %>% 
  mutate(sig_name = paste0(signature_name,"_",coefficient_direction,"_",1)) %>% 
  filter(sig_name %in% large_rad_sigs) %>% 
  pull(p, name = sig_name) -> p_sig_expr
  



#####################
###score the cells###
#####################


radio_genes_df %>% 
  group_by(signature_name, coefficient_direction) %>% 
  summarise(genes = list(unique(gene))) %>% 
  ungroup() %>% 
  filter(length(genes) >= 3) %>% 
  filter(!is.na(coefficient_direction)) %>% 
  mutate(signature = paste0(signature_name, "_", coefficient_direction)) %>% 
  dplyr::select(signature, genes) %>% 
  pull(genes, name = signature) -> rad_sigs_list

SCTall_fil2 = AddModuleScore(SCTall_fil2, features = rad_sigs_list, name = "rad.sig.")
  
for(i in 1:length(names(rad_sigs_list))){
  colnames(SCTall_fil2@meta.data) = gsub(x = colnames(SCTall_fil2@meta.data)
       , pattern = paste0("^rad.sig.",i,"$")
       , replacement = names(rad_sigs_list)[i]
  )
}



radio_genes_df %>% 
  group_by(signature_name, coefficient_direction) %>% 
  summarise(genes = list(unique(gene))) %>% 
  ungroup() %>% 
  mutate(l = lengths(genes)) %>%  
  filter(l >=10) %>% 
  mutate(full_sig = paste0(signature_name, "_", coefficient_direction,"_1")) %>%  
  pull(full_sig) -> 
  large_rad_sigs




SCTall_fil2@meta.data %>% 
  filter(medium_clusters %in% c("Tumor Epithelial", "Cycling Epithelial")) %>% 
  filter(timepoint == "pre") %>% 
  dplyr::select(patient,group,irds.1,80:105) %>% 
  reshape2::melt(value.name = "score", variable.name = "signature") %>% 
  ggplot() + aes(x = patient, y = score, color = group) + geom_boxplot() + 
  facet_wrap(~signature, scales = "free") + 
  scale_color_manual(values = group_colors)
  

SCTall_fil2@meta.data %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial", "Cycling Epithelial")) %>% 
  filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "pre") %>% 
  dplyr::select(patient,group,irds.1,80:105) ->
  rad_sigs_scores

set.seed(77)
samp = sample(nrow(rad_sigs_scores))

rad_mat = rad_sigs_scores[samp,c(-1,-2)]
#rad_mat = t(scale(rad_mat, center = T, scale = F))

rad_mat_colanno = HeatmapAnnotation(df = data.frame(patient = rad_sigs_scores[samp,]$patient, group = rad_sigs_scores[samp,]$group), 
                                    col = list(patient = patient_colors, group = group_colors))

#png(file = "./plots/plots13/rad_sigs_heatmap_all.png", height = 7, width = 35, units = "in", res=72)
#Heatmap(t(rad_mat), top_annotation = rad_mat_colanno, show_column_names = F) 
#dev.off()

########################
###med rad sig scores###
########################


SCTall_fil2@meta.data %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial", "Cycling Epithelial")) %>% 
  filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "pre") %>% 
  dplyr::select(patient.2,group,irds.1,80:105) %>% 
  group_by(patient.2, group) %>% 
  summarise(across(everything(), list(median))) ->
  med_rad_sigs_scores

med_rad_mat = med_rad_sigs_scores[,c(-1,-2)]
med_rad_mat = med_rad_mat[,large_rad_sigs]
#med_rad_mat = med_rad_mat[,names(med_rad_mat)[!str_detect(names(med_rad_mat), "neg")]]

temp = colnames(med_rad_mat)
med_rad_mat = apply(med_rad_mat, 2, scale)
row.names(med_rad_mat) = med_rad_sigs_scores$patient.2

roword = colnames(med_rad_mat)[order(sapply(1:ncol(med_rad_mat), function(x){t.test(med_rad_mat[1:4,x],med_rad_mat[5:8,x])$statistic}))]

prog = colnames(med_rad_mat)
prog[str_detect(prog,"pos") & !str_detect(prog, "Veer")] = "unfavorable"
prog[str_detect(prog,"neg") & !str_detect(prog, "Veer")] = "favorable"
prog[str_detect(prog,"pos")] = "favorable"
prog[str_detect(prog,"neg")] = "unfavorable"

#sapply(1:11, function(x){t.test(med_rad_mat[1:4,x],med_rad_mat[5:8,x])$p.value})
# med_rad_mat_colanno = HeatmapAnnotation(df = data.frame(patient = med_rad_sigs_scores$patient.2, group = med_rad_sigs_scores$group), 
#                                     col = list(patient = patients.2_colors, group = group_colors))
med_rad_mat_colanno = HeatmapAnnotation(df = data.frame(group = med_rad_sigs_scores$group),
                                    col = list(group = group_colors))
library(circlize)
col_fun = colorRamp2(c(0,1), c("white", "darkgoldenrod2"))
med_rad_mat_rowanno = rowAnnotation(df = data.frame(p = p_sig_expr[colnames(med_rad_mat)]),
                                    col = list(p = col_fun))


pdf(file = "~/projects/Breast_cancer_radiation/sup_figure_pieces/rad_sigs_heatmap.pdf", width = 7, height = 4, useDingbats = F)
Heatmap(t(med_rad_mat), top_annotation = med_rad_mat_colanno, show_column_names = T, 
        column_order = med_rad_sigs_scores$patient.2, row_order = roword, 
        column_split = med_rad_sigs_scores$group,right_annotation = med_rad_mat_rowanno,
        row_split = prog) 
dev.off()




med_rad_sigs_scores %>% 
  gather(-patient.2, -group, key = sig, value = score) %>% 
  filter(sig %in% large_rad_sigs) %>% 
  group_by(sig) %>%
  summarise(m = median(score)) %>% 
  arrange(-m)
  







sel_sigs = c("Weichselbaum_IRDS","Sjostrom_ARTIC","Kim_RSS","Cui_RSS","Speers_RSS")


med_rad_mat_vars = apply(med_rad_mat, 2, var)
most_var_rad_sigs = names(med_rad_mat_vars[med_rad_mat_vars >= 0.001])
Heatmap(t(med_rad_mat[,most_var_rad_sigs]), top_annotation = med_rad_mat_colanno, show_column_names = T) 



SCTall_fil2@meta.data %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial", "Cycling Epithelial")) %>% 
  #filter(medium_clusters %in% c("Tumor Epithelial")) %>% 
  filter(timepoint == "pre") %>% 
  dplyr::select(patient,group,irds.1,80:105) %>% 
  group_by(patient, group) %>% 
  summarise(across(everything(), list(median))) %>% 
  filter(group != "unclassified") ->
  med_rad_sigs_scores_all



med_rad_mat_all = med_rad_sigs_scores_all[,c(-1,-2)]
med_rad_mat_all = med_rad_mat_all[,large_rad_sigs]
med_rad_mat_all = med_rad_mat_all[,!str_detect(names(med_rad_mat_all), "neg")]



#med_rad_mat = scale(med_rad_mat)

med_rad_mat_colanno = HeatmapAnnotation(df = data.frame(patient = med_rad_sigs_scores_all$patient, group = med_rad_sigs_scores_all$group), 
                                        col = list(patient = patient_colors, group = group_colors))

Heatmap(t(med_rad_mat_all), top_annotation = med_rad_mat_colanno, show_column_names = T) 















#####################
###subtract scores###
#####################


SCTall_fil2@meta.data %>% 
  dplyr::select(cell.names, timepoint, medium_clusters, patient,group,80:105) %>% 
  reshape2::melt(value.name = "score", variable.name = "signature") %>% 
  filter(signature %in% str_sub(large_rad_sigs, 1, -3)) %>%  
  mutate(signature_name = str_sub(signature, 1, -5)) %>% 
  mutate(coefficient_direction = str_sub(signature, -3)) %>%
  dplyr::select(-signature) %>% 
  group_by(cell.names, timepoint, medium_clusters, patient, group, signature_name) %>% 
  spread(key = coefficient_direction, value = score, fill = 0) %>% 
  ungroup() %>% 
  mutate(score = pos-neg) %>% 
  dplyr::select(cell.names,timepoint, patient, medium_clusters, group, signature_name, score) ->
  rad_sig_all_difs
  
rad_sig_all_difs %>% 
  filter(medium_clusters == "Tumor Epithelial") %>% 
  filter(timepoint=="pre") %>%
  ggplot() + aes(x = patient, y = score, color = group) +
  geom_boxplot() + 
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust =1)) +
  facet_wrap(~signature_name) + 
  scale_color_manual(values = group_colors) + 
  scale_x_discrete(limits = irds_order) + 
  theme_real_minimal() -> temp
ggsave(temp, filename = "./plots/plots13/rad_sigs_diff_tum_pre.pdf")


rad_sig_all_difs %>% 
  filter(medium_clusters == "Tumor Epithelial") %>% 
  filter(timepoint=="pre") %>% 
  group_by(patient, group, signature_name) %>% 
  summarise(score = median(score)) %>% 
  spread(key = signature_name, value = score) ->
  rad_sig_tum_difs_mat

rad_sig_tum_difs_mat_colanno = HeatmapAnnotation(df = data.frame(patient = rad_sig_tum_difs_mat$patient, group = rad_sig_tum_difs_mat$group), 
                                        col = list(patient = patient_colors, group = group_colors))


#pdf(file = "./plots/plots13/rad_sig_dif_tum_pre.pdf")
#Heatmap(t(rad_sig_tum_difs_mat[,c(-1,-2)]), top_annotation = rad_sig_tum_difs_mat_colanno, show_column_names = T) 
#dev.off()


Heatmap(t(scale(rad_sig_tum_difs_mat[,c(-1,-2)], center = T, scale = T)), top_annotation = rad_sig_tum_difs_mat_colanno, show_column_names = T) 

#sel_sigs = c("Weichselbaum_IRDS","Sjostrom_MET141","Kim_RSS")
#pdf(file = "./plots/plots13/rad_sig_dif_sel_tum_pre.pdf")
#Heatmap(t(scale(rad_sig_tum_difs_mat[,sel_sigs],scale = F)), top_annotation = rad_sig_tum_difs_mat_colanno, show_column_names = T)
#dev.off()

sel_sigs = c("Weichselbaum_IRDS","Sjostrom_ARTIC","Kim_RSS","Cui_RSS","Speers_RSS")
Heatmap(t(scale(rad_sig_tum_difs_mat[,sel_sigs],scale = F)), top_annotation = rad_sig_tum_difs_mat_colanno, show_column_names = T)





rad_sig_all_difs %>% 
  filter(timepoint=="pre") %>% 
  group_by(patient, group, signature_name) %>% 
  summarise(score = median(score)) %>% 
  spread(key = signature_name, value = score) ->
  rad_sig_all_difs_mat

rad_sig_all_difs_mat_colanno = HeatmapAnnotation(df = data.frame(patient = rad_sig_all_difs_mat$patient, group = rad_sig_all_difs_mat$group), 
                                                 col = list(patient = patient_colors, group = group_colors))


#pdf(file = "./plots/plots13/rad_sig_dif_all_pre.pdf")
#Heatmap(t(rad_sig_all_difs_mat[,c(-1,-2)]), top_annotation = rad_sig_all_difs_mat_colanno, show_column_names = T) 
#dev.off()



Heatmap(t(scale(rad_sig_all_difs_mat[,c(-1,-2)], center = T, scale = F)), top_annotation = rad_sig_all_difs_mat_colanno, show_column_names = T) 

sel_sigs = c("Weichselbaum_IRDS","Sjostrom_MET141","Kim_RSS")
#pdf(file = "./plots/plots13/rad_sig_dif_sel_all_pre.pdf")
#Heatmap(t(scale(rad_sig_all_difs_mat[,sel_sigs],scale = T)), top_annotation = rad_sig_all_difs_mat_colanno, show_column_names = T)
#dev.off()

sel_sigs = c("Weichselbaum_IRDS","Sjostrom_ARTIC","Kim_RSS","Cui_RSS","Speers_RSS")
Heatmap(t(scale(rad_sig_all_difs_mat[,sel_sigs],scale = F)), top_annotation = rad_sig_all_difs_mat_colanno, show_column_names = T)










#temp = sort( sapply(ls(),function(x){object.size(get(x))}), decreasing = T) 
#head(temp, 50)


#purge = names(head(temp, 50))[c(14:19,21,23,24:26)]

#purge = ls()[str_detect(ls(), "all_genes")]

#purge = "dna_damage_sig_scores_tum_meta"

#for(i in purge){
#  saveRDS(object = get(i), file = paste0("./Rda/",i,".rds"), compress = F)
#}

#rm(list = purge)

