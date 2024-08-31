



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

bind_rows(all_mac_sigs_list) %>% 
  write_csv("/volumes/USR1/aschalck/projects/Breast_cancer_radiation/myeloid_sigs.csv", col_names = T)




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

bind_rows(all_tcell_sigs_list) %>% 
  write_csv("/volumes/USR1/aschalck/projects/Breast_cancer_radiation/tcell_sigs.csv", col_names = T)


