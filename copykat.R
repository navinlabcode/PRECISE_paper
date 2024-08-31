
library(copykat)
library(ComplexHeatmap)



####read new copykat calls
ck_files = list.files("./copykat_2/")
ck_pred_files = paste0("./copykat_2/", ck_files[str_detect(ck_files, pattern = "_copykat_prediction.txt")])
ck_pred_list = lapply(ck_pred_files, function(x) read_tsv(x))
ck_pred = do.call(rbind, ck_pred_list)
names(ck_pred)[2] = "ck_new2"


SCTall_fil@meta.data = left_join(SCTall_fil@meta.data, ck_pred)
row.names(SCTall_fil@meta.data) = SCTall_fil$cell.names
SCTall_fil$ck_new[is.na(SCTall_fil$ck_new)] = "undetermined"





BCR15_results = read_tsv("./copykat/BCR15_copykat_CNA_results.txt")
BCR15_results = BCR15_results[,c(-1,-2,-3)]
BCR15_results_mat = t(as.matrix(data.frame(BCR15_results)))
colnames(BCR15_results_mat) = 1:ncol(BCR15_results_mat)
#BCR15_results_m = melt(BCR15_results, variable.name = "cell", value.name = "value", id.vars = "abspos")
#BCR15_results_w = pivot_wider(BCR15_results_m, names_from = abspos, values_from = value)

Heatmap(BCR15_results_mat[1:10,1:1000])


BCR15_copykat_clustering_results = readRDS("./copykat_2/BCR15_copykat_clustering_results.rds")


(BCR15_copykat_clustering_results, col=cutree(BCR15_copykat_clustering_results, k = 6))

labcols = cutree(BCR15_copykat_clustering_results, k = 6)[BCR15_copykat_clustering_results$order]

dend = as.dendrogram(BCR15_copykat_clustering_results)
dend %>% 
  set("labels_col", labcols) %>% 
  plot()

SCTall_fil2_epi@meta.data %>% 
  select(cell.names, patient, timepoint, ck_pred,minor_clusters2) %>% 
  filter(patient == "BCR15") %>% 
  #filter(timepoint == "post") %>% 
  left_join(data.frame(cell.names = SCTall_fil2@meta.data$cell.names, irds = SCTall_fil2@meta.data$irds.1)) %>% 
  left_join(data.frame(cell.names = names(labcols), ck_clus = unname(labcols))) %>% 
  filter(ck_clus %in% c(1,2,4,6)) %>% 
  ggplot() + aes(x = as.factor(ck_clus), y = irds) + 
  geom_boxplot()
    


#ComplexHeatmap::Heatmap(BCR15_results_mat, column_dend_reorder = F, 
#                        show_column_dend = F, row_dend_reorder = F, show_row_names = F)

#BCR15_results_dist = dist(BCR15_results_mat)




ck_hclust_files = paste0("./copykat_2/", ck_files[str_detect(ck_files, pattern = "_copykat_clustering_results.rds")])


copykat_clustering_results = lapply(ck_hclust_files, function(x){readRDS(x)})




