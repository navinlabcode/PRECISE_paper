library(Seurat)
library(circlize)


med_sig_scores = function(so, sig_names){
  sig_scores = setNames(data.frame(matrix(ncol = length(sig_names), nrow = ncol(so))), sig_names)
  for(i in sig_names){
    print(i)
    genes = pull(m_df[m_df$gs_name == i,"gene_symbol"])
    genes_fil = genes[genes %in% row.names(so)]
    so = AddModuleScore(so, features = list(genes_fil), name = i)
    sig_scores[i] = so@meta.data[,paste0(i,"1")]
  }
  row.names(sig_scores) = colnames(so)
  
  cbind(sig_scores, so@meta.data[row.names(sig_scores),]) -> sig_scores_meta
  return(sig_scores_meta)
}

sig_scores_meta = med_sig_scores(so, sig_names)

top_sig_scores = function(sig_scores_meta,sig_names){
  
  
  sig_scores_meta %>% 
    melt(value.name = "sig_score", variable.name = "signature") %>% 
    group_by(patient.2, timepoint, signature) %>% 
    summarise(m = median(sig_score)) %>% 
    filter(signature %in% sig_names) %>% 
    filter(patient.2 %in% paste0("P0",1:8)) %>% 
    #filter(str_detect(signature, "WP|REACTOME_")) %>% 
    group_by(signature) %>% 
    mutate(pval = t.test(m ~ timepoint, paired=T)$p.value, stat = t.test(m ~ timepoint, paired=T)$statistic) %>% 
    filter(pval<= 0.05) %>% 
    pull(signature) %>% 
    unique() -> pathways_top
  
  
  sig_scores_meta %>% 
    melt(value.name = "sig_score", variable.name = "signature") %>% 
    group_by(patient.2, timepoint, signature) %>% 
    summarise(m = median(sig_score)) %>% 
    filter(signature %in% sig_names) %>% 
    filter(patient.2 %in% paste0("P0",1:8)) %>% 
    spread(key = timepoint, value = m) %>% 
    mutate(d = post-pre) %>% 
    filter(signature %in% pathways_top) %>% 
    group_by(signature) %>% 
    mutate(m = median(d)) %>% 
    filter(abs(m) >= 0.01) %>% 
    select(-pre, -post, -m) %>% 
    spread(key = signature, value = d) %>% 
    column_to_rownames("patient.2") %>% 
    as.data.frame() -> sig_difs_top
  
  return(sig_difs_top)
  
}

sig_difs_top = top_sig_scores(sig_scores_meta, sig_names)

all_sig_scores = function(sig_scores_meta,sig_names){

  sig_scores_meta %>% 
    melt(value.name = "sig_score", variable.name = "signature") %>% 
    group_by(patient.2, timepoint, signature) %>% 
    summarise(m = median(sig_score)) %>% 
    filter(signature %in% sig_names) %>% 
    filter(patient.2 %in% paste0("P0",1:8)) %>% 
    spread(key = timepoint, value = m) %>% 
    mutate(d = post-pre) %>% 
    select(-pre, -post) %>% 
    spread(key = signature, value = d) %>% 
    column_to_rownames("patient.2") %>% 
    as.data.frame() -> sig_difs_all
  
  return(sig_difs_all)
  
}

sig_difs_all = all_sig_scores(sig_scores_meta, sig_names)

sig_difs_plot = function(sig_difs_top){
  
  m = max(max(sig_difs_top), abs(min(sig_difs_top)))
  
  col_fun = colorRamp2(c(-round(m, 2), 0, round(m, 2)), c(treatment_colors["pre"], "white", treatment_colors["post"]))
  
  hm = Heatmap(t(sig_difs_top), 
               column_order = paste0("P0",1:8), 
               col = col_fun, row_names_side = "left", 
               show_row_dend = F)
  return(hm)
  
}

draw(sig_difs_plot(sig_difs_top), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))


#######
##run##
#######

m_df = msigdbr(species = "Homo sapiens", category = "C2")

m_df %>% 
  filter(str_detect(gs_name, "DNA_DAMAGE") | 
           str_detect(gs_name, "REPAIR") | 
           str_detect(gs_name, "HDR|HOMOLOGOUS") |
           str_detect(gs_name, "BREAKS|BREAK_") |
           str_detect(gs_name, "STRAND")) %>% 
  pull(gs_name) %>% 
  unique() -> dna_damage_sig_names

dna_damage_sig_names = dna_damage_sig_names[str_detect(dna_damage_sig_names, "^REACTOME_")]

dna_damage_sig_names = dna_damage_sig_names[!str_detect(dna_damage_sig_names, "PROTEIN_REPAIR")]
#dna_damage_sig_names = dna_damage_sig_names[1:3]

dna_damage_sig_scores_meta = med_sig_scores(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"], dna_damage_sig_names)

dna_damage_sig_difs_top = top_sig_scores(dna_damage_sig_scores_meta, dna_damage_sig_names)
draw(sig_difs_plot(dna_damage_sig_difs_top), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))

dna_damage_sig_difs_all = all_sig_scores(dna_damage_sig_scores_meta, dna_damage_sig_names)
draw(sig_difs_plot(dna_damage_sig_difs_all), padding = unit(c(0.1, 8, 0.1, 0.1), "in"))

pdf(file = "~/projects/Breast_cancer_radiation/sup_figure_pieces/reactome_dna_damage_sig_difs_all.pdf", width = 14, height = 7)
draw(sig_difs_plot(dna_damage_sig_difs_all), padding = unit(c(0.1, 8, 0.1, 0.1), "in"))
dev.off()


#write_csv()


m_df %>% 
  filter(str_detect(gs_name, "STRESS") | str_detect(gs_name, "HYPOXIA")) %>% 
  pull(gs_name) %>% 
  unique() -> cell_stress_sig_names

cell_stress_sig_names = cell_stress_sig_names[!str_detect(cell_stress_sig_names, "PROTEIN_REPAIR")]
#cell_stress_sig_names = cell_stress_sig_names[1:3]

cell_stress_sig_scores_meta = med_sig_scores(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"], cell_stress_sig_names)
cell_stress_sig_difs_top = top_sig_scores(cell_stress_sig_scores_meta, cell_stress_sig_names)
draw(sig_difs_plot(cell_stress_sig_difs_top), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))

cell_stress_sig_difs_all = all_sig_scores(cell_stress_sig_scores_meta, cell_stress_sig_names)
draw(sig_difs_plot(cell_stress_sig_difs_all), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))

cell_stress_sig_difs_reactome = cell_stress_sig_difs_all[,str_detect(colnames(cell_stress_sig_difs_all), "REACTOME")]
draw(sig_difs_plot(cell_stress_sig_difs_reactome), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))


pdf(file = "~/projects/Breast_cancer_radiation/sup_figure_pieces/reactome_cell_stress_sig_difs_all.pdf", width = 14, height = 3)
draw(sig_difs_plot(cell_stress_sig_difs_reactome), padding = unit(c(0.1, 8, 0.1, 0.1), "in"))
dev.off()






m_df %>% 
  filter(str_detect(gs_name, "DEATH") | str_detect(gs_name, "TOSIS")) %>% 
  filter(!str_detect(gs_name, "MITOSIS")) %>% 
  filter(!str_detect(gs_name, "DYSOSTOSIS")) %>% 
  filter(!str_detect(gs_name, "SIDS")) %>% 
  filter(!str_detect(gs_name,"LEIOMYOMATOSIS")) %>% 
  filter(!str_detect(gs_name,"_DN$")) %>% 
  pull(gs_name) %>% 
  unique() -> cell_death_sig_names


cell_death_sig_scores_meta = med_sig_scores(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"], cell_death_sig_names)

cell_death_sig_difs_top = top_sig_scores(cell_death_sig_scores_meta, cell_death_sig_names)
draw(sig_difs_plot(cell_death_sig_difs_top), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))

cell_death_sig_difs_all = all_sig_scores(cell_death_sig_scores_meta, cell_death_sig_names)
cell_death_sig_difs_reactome = cell_death_sig_difs_all[,str_detect(colnames(cell_death_sig_difs_all), "REACTOME")]


draw(sig_difs_plot(cell_death_sig_difs_all), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))
draw(sig_difs_plot(cell_death_sig_difs_reactome), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))

pdf(file = "~/projects/Breast_cancer_radiation/sup_figure_pieces/reactome_cell_death_sig_difs_all.pdf", width = 14, height = 5)
draw(sig_difs_plot(cell_death_sig_difs_reactome), padding = unit(c(0.1, 8, 0.1, 0.1), "in"))
dev.off()


m_df %>% 
  filter(gs_name %in% c(cell_death_sig_names, cell_stress_sig_names, dna_damage_sig_names)) %>%
  #filter(str_detect(gs_name, "^WP_") | str_detect(gs_name, "^KEGG_") | str_detect(gs_name, "^REACTOME_")) %>% 
  filter(str_detect(gs_name, "^REACTOME_")) %>% 
  pull(gs_name) %>% 
  unique() -> canon_sig_names 

canon_sig_scores_meta = med_sig_scores(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"], canon_sig_names)

canon_sig_difs_top = top_sig_scores(canon_sig_scores_meta, canon_sig_names)
draw(sig_difs_plot(canon_sig_difs_top), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))


canon_sig_difs_all = all_sig_scores(canon_sig_scores_meta, canon_sig_names)
draw(sig_difs_plot(canon_sig_difs_all), padding = unit(c(0.1, 5, 0.1, 0.1), "in"))





















#########
###old###
#########

dna_damage_sig_scores_tum_meta = read_csv("./data/dna_damage_sig_scores_tum_meta.csv")

dna_damage_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% dna_damage_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  #filter(str_detect(signature, "WP|REACTOME_")) %>% 
  group_by(signature) %>% 
  mutate(pval = t.test(m ~ timepoint, paired=T)$p.value, stat = t.test(m ~ timepoint, paired=T)$statistic) %>% 
  filter(pval<= 0.05) %>% 
  pull(signature) %>% 
  unique() -> DNA_damage_pathways_top

dna_damage_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% dna_damage_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  spread(key = timepoint, value = m) %>% 
  mutate(d = post-pre) %>% 
  filter(signature %in% DNA_damage_pathways_top) %>% 
  group_by(signature) %>% 
  mutate(m = median(d)) %>% 
  filter(abs(m) >= 0.01) %>% 
  #ungroup() %>% 
  #top_n(48,m) %>% 
  select(-pre, -post, -m) %>% 
  spread(key = signature, value = d) %>% 
  column_to_rownames("patient.2") %>% 
  as.data.frame() -> dna_damage_sig_difs_top


#library(circlize)
col_fun = colorRamp2(c(-0.2, 0, 0.2), c(treatment_colors["pre"], "white", treatment_colors["post"]))

hm = Heatmap(t(dna_damage_sig_difs_top), 
        column_order = paste0("P0",1:8), 
        col = col_fun, row_names_side = "left", 
        show_row_dend = F)
draw(hm, padding = unit(c(0.1, 5, 0.1, 0.1), "in"))


m_df %>% 
  filter(gs_name == "AMUNDSON_DNA_DAMAGE_RESPONSE_TP53") %>% 
  pull(gene_symbol) %>% 
  unique() -> temp

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="Tumor Epithelial" & SCTall_fil2$patient.2 %in% paste0("P0",1:8)], 
        group.by = "patient.2",same.y.lims = T,
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors, ncol= 4,
        features = temp) + 
  NoLegend() -> AMUNDSON_vln

ggsave(AMUNDSON_vln, filename = "~/projects/Breast_cancer_radiation/sup_figure_pieces/AMUNDSON_vln.pdf", device = cairo_pdf, height = 8, width = 12)


temp = SCTall_fil2[,SCTall_fil2$medium_clusters=="Tumor Epithelial" & SCTall_fil2$patient.2 %in% paste0("P0",1:8)]


#################
###cell stress###
#################





cell_stress_sig_scores_tum_meta = read_csv("./data/cell_stress_sig_scores_tum_meta.csv")

cell_stress_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% cell_stress_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  spread(key = timepoint, value = m) %>% 
  mutate(d = post-pre) %>% 
  select(-pre, -post) %>% 
  spread(key = signature, value = d) %>% 
  column_to_rownames("patient.2") %>% 
  as.data.frame() -> cell_stress_sig_difs

cell_stress_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% cell_stress_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  #filter(str_detect(signature, "WP|REACTOME_")) %>% 
  group_by(signature) %>% 
  mutate(pval = t.test(m ~ timepoint, paired=T)$p.value, stat = t.test(m ~ timepoint, paired=T)$statistic) %>%
  #filter(patient.2 == "P01") %>% 
  #filter(timepoint == "pre") %>% 
  #pull(pval) %>% 
  #sort() %>% 
  #p.adjust(method = "BH")
  filter(pval<= 0.05) %>% 
  pull(signature) %>% 
  unique() -> cell_stress_pathways_top


cell_stress_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% cell_stress_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  spread(key = timepoint, value = m) %>% 
  mutate(d = post-pre) %>% 
  filter(signature %in% cell_stress_pathways_top) %>% 
  group_by(signature) %>% 
  mutate(m = median(d)) %>% 
  filter(abs(m) >= 0.01) %>% 
  #ungroup() %>% 
  #top_n(48,m) %>% 
  select(-pre, -post, -m) %>% 
  spread(key = signature, value = d) %>% 
  column_to_rownames("patient.2") %>% 
  as.data.frame() -> cell_stress_sig_difs_top



#library(circlize)
col_fun = colorRamp2(c(-0.2, 0, 0.2), c(treatment_colors["pre"], "white", treatment_colors["post"]))

hm = Heatmap(t(cell_stress_sig_difs_top), 
             column_order = paste0("P0",1:8), 
             col = col_fun, row_names_side = "left", 
             show_row_dend = F)
draw(hm, padding = unit(c(0.1, 5, 0.1, 0.1), "in"))

###########
###gh2ax###
###########

GammaH2AX = read_csv("./data/IHC/GammaH2AX.csv")

GammaH2AX %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  ggpaired(x = "timepoint", y = "pH2ax",
           color = "timepoint", line.color = "black", point.size = 2,
           line.size = 0.5,palette = c("black","black")) +
  #stat_compare_means(paired = TRUE) + 
  ylab("GammaH2Ax (H-score)") +
  theme_one_line + theme(legend.position = "none") ->
  ghzax_IHC_plot

ggsave(ghzax_IHC_plot, filename = "~/projects/Breast_cancer_radiation/sup_figure_pieces/ghzax_IHC_plot.pdf", width=2.5, height=3, device = cairo_pdf)

#p-value = 0.02731, paired t.test

GammaH2AX %>% 
  spread(timepoint,pH2ax) %>% 
  summarise(ttest = list(t.test(pre, post,paired = TRUE))) %>% 
  pull(ttest)

   t.test(GammaH2AX$pH2ax[1:13], GammaH2AX$pH2ax[14:26], paired = T)


