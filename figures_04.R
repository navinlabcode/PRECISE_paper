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
library(org.Hs.eg.db)
library(GO.db)
library(msigdbr)
library(fgsea)
library(ggpubr)
source("~/TENX/analysis/MP/MP_project_4/TCR_functions_script.R")
options(future.globals.maxSize = 60000 * 1024^2)
NewSet1 = function(n){return(colorRampPalette(brewer.pal(9,"Set1")[c(1,5,6,3,2,4,8,7,9)])(n))}
library(ComplexHeatmap)




tumor_prepost_list_gsea_H = lapply(tumor_prepost_list[c(1,2,3,5,6,8,9,10)], function(x){rungsea(x, cats = "H", fc_cutoff = 0.1, mincat = 5)})
tumor_prepost_list_gsea_C2 = lapply(tumor_prepost_list[c(1,2,3,5,6,8,9,10)], function(x){rungsea(x, cats = "C2", fc_cutoff = 0.1, mincat = 5)})
tumor_prepost_list_gsea_C5 = lapply(tumor_prepost_list[c(1,2,3,5,6,8,9,10)], function(x){rungsea(x, cats = "C5", subcat = "MF", fc_cutoff = 0.1, mincat = 5)})
tumor_prepost_list_gsea_C3 = lapply(tumor_prepost_list[c(1,2,3,5,6,8,9,10)], function(x){rungsea(x, cats = "C3", subcat = "TFT:GTRD", fc_cutoff = 0.1, mincat = 5)})


tumor_prepost_df_gsea_H = bind_rows(tumor_prepost_list_gsea_H, .id = "patient")
tumor_prepost_df_gsea_C2 = bind_rows(tumor_prepost_list_gsea_C2, .id = "patient")
tumor_prepost_df_gsea_C5 = bind_rows(tumor_prepost_list_gsea_C5, .id = "patient")
tumor_prepost_df_gsea_C3 = bind_rows(tumor_prepost_list_gsea_C3, .id = "patient")


tumor_prepost_df_gsea_H %>% 
  filter(cluster == "post") %>% 
  group_by(pathway) %>% 
  mutate(order = sum(NES)) %>% 
  arrange(-order) %>%
  ungroup() %>%
  filter(dense_rank(order) <= 1*5 | dense_rank(desc(order)) <= 1*5) %>%
  #top_n(8*5, order) %>% 
  ggplot() + aes(x = reorder(pathway, order), y = NES, fill=NES) + 
  geom_bar(stat = "identity") + 
  facet_grid(cols = vars(patient)) + ylab("Normalized Enrichment Score") +
  coord_flip() + 
  scale_fill_gradient2(low="steelblue", mid="white", high = "firebrick", 
                       midpoint = 0, limits = c(-4,4), breaks=c(-4,-2,0,2,4)) +
  scale_y_continuous(limits = c(-4,4)) +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank()) -> gsea_H_epi_post_plot
#ggsave(gsea_H_epi_post_plot, filename = "./plots/fig_pieces_03/gsea_H_epi_post.pdf", height = 3, width = 20)

tum_gsea_legend = get_legend(gsea_H_epi_post_plot)

gsea_H_epi_post_plot = gsea_H_epi_post_plot + theme(legend.position = "none", axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank(), strip.text = element_text(size=16)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


tumor_prepost_df_gsea_C2 %>% 
  filter(cluster == "post") %>% 
  filter(str_detect(pathway, "^KEGG")) %>% 
  group_by(pathway) %>% 
  mutate(order = sum(NES)) %>% 
  arrange(-order) %>%
  ungroup() %>%
  filter(dense_rank(order) <= 1*5 | dense_rank(desc(order)) <= 1*5) %>% 
  ggplot() + aes(x = reorder(pathway, order), y = NES, fill=NES) + 
  geom_bar(stat = "identity") + 
  facet_grid(cols = vars(patient)) + xlab("Normalized Enrichment Score") +
  coord_flip() +  
  scale_fill_gradient2(low="steelblue", mid="white", high = "firebrick", 
                       midpoint = 0, limits = c(-4,4), breaks=c(-4,-2,0,2,4)) +
  scale_y_continuous(limits = c(-4,4)) +
  theme(panel.background = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank()) -> gsea_C2_KEGG_epi_post_plot
#ggsave(gsea_C2_KEGG_epi_post_plot, filename = "./plots/fig_pieces_03/gsea_C2_KEGG_epi_post.pdf", height = 3, width = 20)


gsea_C2_KEGG_epi_post_plot = gsea_C2_KEGG_epi_post_plot + theme(legend.position = "none", axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank(), strip.text = element_blank()) +
 theme(plot.margin = unit(c(0.1,0,0,0), "cm"))

tumor_prepost_df_gsea_C2 %>% 
  filter(cluster == "post") %>% 
  filter(str_detect(pathway, "^REAC")) %>% 
  group_by(pathway) %>% 
  mutate(order = sum(NES)) %>% 
  arrange(-order) %>%
  ungroup() %>%
  filter(dense_rank(order) <= 1*5 | dense_rank(desc(order)) <= 1*5) %>% 
  ggplot() + aes(x = reorder(pathway, order), y = NES, fill=NES) + 
  geom_bar(stat = "identity") + 
  facet_grid(cols = vars(patient)) + xlab("Normalized Enrichment Score") +
  coord_flip() + 
  scale_fill_gradient2(low="steelblue", mid="white", high = "firebrick", 
                       midpoint = 0, limits = c(-4,4), breaks=c(-4,-2,0,2,4)) +
  scale_y_continuous(limits = c(-4,4)) +
  theme(panel.background = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank()) -> gsea_C2_reactome_epi_post_plot
#ggsave(gsea_C2_reactome_epi_post_plot, filename = "./plots/fig_pieces_03/gsea_C2_reactome_epi_post.pdf", height = 3, width = 20)


gsea_C2_reactome_epi_post_plot = gsea_C2_reactome_epi_post_plot + theme(legend.position = "none", axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank(), strip.text = element_blank()) +
  theme(plot.margin = unit(c(0.1,0,0,0), "cm"))

tumor_prepost_df_gsea_C2 %>% 
  filter(cluster == "post") %>% 
  filter(str_detect(pathway, "^WP_")) %>% 
  group_by(pathway) %>% 
  mutate(order = sum(NES)) %>% 
  arrange(-order) %>%
  ungroup() %>%
  filter(dense_rank(order) <= 1*5 | dense_rank(desc(order)) <= 1*5) %>% 
  ggplot() + aes(x = reorder(pathway, order), y = NES, fill=NES) + 
  geom_bar(stat = "identity") + 
  facet_grid(cols = vars(patient)) + xlab("Normalized Enrichment Score") +
  coord_flip() + 
  scale_fill_gradient2(low="steelblue", mid="white", high = "firebrick", 
                       midpoint = 0, limits = c(-4,4), breaks=c(-4,-2,0,2,4)) +
  scale_y_continuous(limits = c(-4,4)) +theme(panel.background = element_blank(), 
                                              axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                                              strip.background = element_blank()) -> gsea_C2_wikipaths_epi_post_plot
#ggsave(gsea_C2_wikipaths_epi_post_plot, filename = "./plots/fig_pieces_03/gsea_C2_wikipaths_epi_post.pdf", height = 3, width = 20)

gsea_C2_wikipaths_epi_post_plot = gsea_C2_wikipaths_epi_post_plot + theme(legend.position = "none", axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank(), strip.text = element_blank()) +
  theme(plot.margin = unit(c(0.1,0,0,0), "cm"))




tumor_prepost_df_gsea_C3 %>% 
  filter(cluster == "post") %>% 
  #filter(str_detect(pathway, "^WP_")) %>% 
  group_by(pathway) %>% 
  mutate(order = sum(NES)) %>% 
  arrange(-order) %>%
  ungroup() %>%
  filter(dense_rank(order) <= 1*5 | dense_rank(desc(order)) <= 1*5) %>% 
  ggplot() + aes(x = reorder(pathway, order), y = NES, fill=NES) + 
  geom_bar(stat = "identity") + 
  facet_grid(cols = vars(patient)) + xlab("Normalized Enrichment Score") +
  coord_flip() + scale_fill_gradient2(low="steelblue", mid="white", high = "firebrick", 
                                      midpoint = 0, limits = c(-4,4), breaks=c(-4,-2,0,2,4)) +
  scale_y_continuous(limits = c(-4,4)) +theme(panel.background = element_blank(), 
                                              axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                                              strip.background = element_blank()) -> gsea_C3_tfts_epi_post_plot
#ggsave(gsea_C3_tfts_epi_post_plot, filename = "./plots/fig_pieces_03/gsea_C3_tfts_epi_post.pdf", height = 3, width = 20)
#


gsea_C3_tfts_epi_post_plot = gsea_C3_tfts_epi_post_plot + theme(legend.position = "none", strip.text = element_blank(), axis.line.x = element_line(size=0.5), axis.title.x =element_blank()) +
  theme(plot.margin = unit(c(0.1,0,0,0), "cm"))


tum_gsea_left = plot_grid(gsea_H_epi_post_plot, gsea_C2_KEGG_epi_post_plot,gsea_C2_reactome_epi_post_plot,gsea_C2_wikipaths_epi_post_plot,gsea_C3_tfts_epi_post_plot, ncol = 1, align = "v", axis = "lt")




y.grob <- textGrob("Normalized Enrichment Score", 
                   gp=gpar(col="black", fontsize=20),hjust = -0.3, rot=0)

tum_gsea_left2 = grid.arrange(arrangeGrob(tum_gsea_left, bottom = y.grob, padding = unit(1,"cm"))) 


tum_gsea_plot = plot_grid(tum_gsea_left2, tum_gsea_legend, nrow = 1, rel_widths = c(1,0.1))
ggsave(tum_gsea_plot,  filename = "./plots/fig_pieces_03/tum_gsea_plot.pdf", width=18, height = 10)



####################
###pre-pos gsea 2###
####################



tumor_prepost_df_gsea_C2 %>% 
  filter(cluster == "post") %>% 
  filter(str_detect(pathway, "^KEGG|^WP_|^REAC")) %>% 
  group_by(pathway) %>% 
  mutate(order = sum(NES)) %>% 
  arrange(-order) %>%
  ungroup() %>%
  filter(dense_rank(order) <= 1*5 | dense_rank(desc(order)) <= 1*5) %>% 
  ggplot() + aes(x = reorder(pathway, order), y = NES, fill=NES) + 
  geom_bar(stat = "identity") + 
  facet_grid(cols = vars(patient)) + xlab("Normalized Enrichment Score") +
  coord_flip() +  
  scale_fill_gradient2(low="steelblue", mid="white", high = "firebrick", 
                       midpoint = 0, limits = c(-4,4), breaks=c(-4,-2,0,2,4)) +
  scale_y_continuous(limits = c(-4,4)) +
  theme(panel.background = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank()) -> gsea_C2_epi_post_plot

gsea_C2_epi_post_plot = gsea_C2_epi_post_plot + theme(legend.position = "none", axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank(), strip.text = element_blank()) +
  theme(plot.margin = unit(c(0.1,0,0,0), "cm"))



tum_gsea_left = plot_grid(gsea_H_epi_post_plot, gsea_C2_epi_post_plot,gsea_C3_tfts_epi_post_plot, ncol = 1, align = "v", axis = "lt")

y.grob <- textGrob("Normalized Enrichment Score", 
                   gp=gpar(col="black", fontsize=20),hjust = -0.3, rot=0)

tum_gsea_left2 =grid.arrange(arrangeGrob(tum_gsea_left, bottom = y.grob, padding = unit(1,"cm"))) 


tum_gsea_plot = plot_grid(tum_gsea_left2, tum_gsea_legend, nrow = 1, rel_widths = c(1,0.1))
ggsave(tum_gsea_plot,  filename = "./plots/fig_pieces_03/tum_gsea_plot.pdf", width=18, height = 6)






#######################
####pre-post gsea 3####
#######################


tumor_prepost_df_gsea_H = readRDS("./Rda/tumor_prepost_df_gsea_H.rds")

tumor_prepost_df_gsea_H %>% 
  left_join(patients_key_df) %>% 
  dplyr::select(-patient) %>% 
  filter(cluster == "post") %>% 
  filter(padj <= 0.05) %>%
  mutate(pathway = str_replace(pathway, "HALLMARK_","")) %>% 
  group_by(pathway) %>% 
  mutate(order = sum(NES)) %>% 
  filter(abs(order) >= 5) %>% 
  arrange(-order) %>% 
  ungroup() %>% 
  ggplot() + aes(x = patient.2, y = reorder(pathway, order), size=NES, color=NES) +
  geom_point() + 
  scale_color_gradient2(low="royalblue4", mid="white", high = "goldenrod", 
                                      midpoint = 0, limits = c(-4,4), breaks=c(-4,-2,0,2,4)) +
  theme_real_minimal() +
  theme(axis.title = element_blank()) + 
  theme(text = element_text(size=16)) +
  scale_x_discrete(position = "top") -> temp
ggsave(temp, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/tum_gsea_H2.pdf", width = 8, height = 5, device = cairo_pdf)




####################
###tum expression###
####################

Idents(SCTall_fil2_epi) = "patient"
tumor_prepost_avg = AverageExpression(SCTall_fil2_epi, features = tumor_post_avg_top, return.seurat = T, add.ident = "timepoint", assays = "SCT", slot = "data")

DoHeatmap(tumor_prepost_avg, features = tumor_post_avg_top, assay = "SCT", slot = "data", disp.max = 2, lines.width = 1, size = 5, 
          group.colors = patient_colors[c(1:5,7:11)]) +
  scale_fill_gradientn(colors = c("black","firebrick","orange","yellow"), na.value = "white") -> prepost_epi_avg_plot
ggsave(prepost_epi_avg_plot, filename = "./plots/fig_pieces_03/prepost_epi_avg.pdf", width = 6, height = 4)



#####################
####tum cycling######
#####################



cellcycle_colors = c("#FFA177","#E0BC6E","#949C73")

SCTall_fil2@meta.data %>% 
  filter(str_detect(minor_clusters2, "^Epithelial")) %>% 
  filter(patient != "BCR17") %>% 
  dplyr::select(Phase, minor_clusters2, patient, timepoint) %>% 
  mutate(Phase = factor(Phase, levels = c("G2M", "S", "G1"))) %>% 
  ggplot() + aes(x=timepoint, fill=Phase) + 
  geom_bar(position="fill") + scale_fill_manual(values = cellcycle_colors) +
  facet_wrap(~patient, nrow = 2) + ylab("Percent of Tumor Cells") +
  theme(panel.background = element_blank(), 
        axis.title.x= element_blank(), strip.background = element_blank()) -> temp
ggsave(temp, filename = "./plots/fig_pieces_03/phase_epi_bar.pdf", width = 4.5, height = 2.5)



######################
####KI67##############
######################


KI67_total = read_csv("./data/IHC/ki67.csv")
KI67_total %>% 
  left_join(id_key[,c("IMT_id", "patient")], by = c("pa_id" = "IMT_id")) ->
  KI67_total_merged



KI67_total_merged %>% 
  filter(!is.na(patient)) %>% 
  mutate(KI67_pre = if_else(is.na(KI67_pre), KI67_dx,KI67_pre)) %>% 
  dplyr::select(patient,KI67_pre,KI67_post) %>% 
  melt(variable.name = c("timepoint"), value.name = "cycling") %>% 
  mutate(timepoint = str_sub(timepoint, 6)) %>% 
  mutate(assay = "IHC") %>% 
  mutate(cycling = cycling*.01) ->
  cycling_IHC


# cycling_IHC %>% 
#   select(-assay) %>% 
#   spread(key = timepoint, value = cycling, fill=0) %>% 
#   mutate(d = post-pre) %>% 
#   pull(d) %>% wilcox.test(exact=F)
#p-value = 0.3202


SCTall_fil2_epi$mki67 = if_else(SCTall_fil2_epi@assays$RNA@counts["MKI67",] > 0, "pos", "neg")
SCTall_fil2$mki67 = if_else(SCTall_fil2@assays$RNA@counts["MKI67",] > 0, "pos", "neg")



data.frame(SCTall_fil2@meta.data) %>% 
  filter(str_detect(minor_clusters2, "^Epithelial")) %>% 
  group_by(patient, timepoint, mki67) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  complete(patient, timepoint, mki67, fill = list(n=0)) %>% 
  group_by(patient, timepoint) %>% 
  mutate(N = sum(n)) %>% 
  mutate(cycling = n/N, assay="sc") %>% 
  filter(mki67 == "pos") %>% 
  dplyr::select(patient,timepoint,cycling,assay) ->
  mki67_sc


# mki67_sc %>% 
#   select(-assay) %>% 
#   spread(key = timepoint, value = cycling, fill=0) %>% 
#   mutate(d = post-pre) %>% 
#   pull(d) %>% wilcox.test(exact=F)
##p-value = 0.3066

rbind(data.frame(mki67_sc),  data.frame(cycling_IHC)) %>% 
  ggplot() + aes(x = timepoint, y = cycling, group = assay, color=assay) + 
  geom_point() +
  geom_line() +
  facet_wrap(~patient, ncol = 6) +
  scale_color_manual(values = assay_colors) +
  scale_y_continuous(limits = c(0,1)) +
  theme(panel.background = element_blank()) + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=20)) +
  theme(strip.background = element_blank(), strip.text = element_text(size=16)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.line.y = element_line(size = 0.5)) +
  theme(legend.title = element_text(size=16)) +
  theme(legend.text  = element_text(size=14)) +
  ylab("Percent KI67+") -> temp 
ggsave(temp, filename = "ki67_fraction.pdf",path = "./plots/fig_pieces_03/", width = 8, height = 4, device = cairo_pdf)

SCTall_fil2@meta.data %>% 
  filter(str_detect(minor_clusters2, "^Epithelial")) %>% 
  group_by(patient.2, timepoint, Phase) %>% 
  summarise(n = n()) %>% 
  mutate(N = sum(n)) %>% 
  filter(N >= 50) %>% 
  mutate(p = n/N) %>% 
  select(-n,-N) %>% 
  ungroup()  ->
  phase_sc

# phase_sc %>% 
#   filter(Phase %in% c("G2M","S")) %>% 
#   group_by(patient.2, timepoint) %>% 
#   summarise(p = sum(p)) %>% 
#   spread(key = timepoint, value = p, fill=0) %>% 
#   mutate(d = post-pre) %>% 
#   pull(d) %>% wilcox.test(exact=F)
#p = 0.359

phase_sc %>% 
  filter(Phase != "G1") %>% 
  ggplot() + aes(x = timepoint, y = p, group=patient.2) +
  facet_wrap(~Phase) +
  geom_point() +
  geom_line() + scale_y_continuous(limits = c(0,1)) +
  theme_real_minimal() + ylab("Fraction in Phase") +
  theme(axis.title.x = element_blank()) + 
  theme(text = element_text(size=12)) -> 
  phase_lineplots

rbind(data.frame(mki67_sc),  data.frame(cycling_IHC)) %>% 
  left_join(clinical_data[,c("patient","patient.2")]) %>% 
  ggplot() + aes(x = timepoint, y = cycling, group=patient.2) +
  facet_wrap(~assay) +
  geom_point() +
  geom_line() +
  theme_real_minimal() + ylab("Fraction KI67+") +
  theme(axis.title.x = element_blank()) + 
  theme(text = element_text(size=12)) ->
  ki67_lineplots


ggsave(ki67_lineplots, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/ki67_lineplots.pdf", height = 2, width = 4, device = cairo_pdf)


plot_grid(ki67_lineplots,phase_lineplots,  ncol = 1) -> cycle_lineplots

ggsave(cycle_lineplots, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/cycle_lineplots.pdf", height = 3, width = 4, device = cairo_pdf)



###########
####HRs####
###########

VlnPlot(SCTall_fil2_epi, features = c("ESR1","PGR","ERBB2","AR"),
        group.by = "patient", split.plot = TRUE, split.by = "timepoint", 
        pt.size = 0, ncol = 2, cols = treatment_colors,y.max = 3) -> temp
ggsave(temp, filename = "./plots/fig_pieces_03/hormone_receptor_expression.pdf", width = 9, height = 7)



table(SCTall_fil2$sample,SCTall_fil2$major_clusters)[,8]/sum(table(SCTall_fil2$sample,SCTall_fil2$major_clusters)[,8])

table(SCTall_fil2$sample,SCTall_fil2$major_clusters)[,8]



VlnPlot(SCTall_fil2_epi, features = c("ESR1","PGR","ERBB2","AR"),
        group.by = "timepoint", 
        pt.size = 0, ncol = 2, cols = treatment_colors,y.max = 3) -> temp
ggsave(temp, filename = "./plots/fig_pieces_03/hormone_receptor_expression.pdf", width = 9, height = 7)



SCTall_fil2_epi$ESR1 = SCTall_fil2_epi@assays$RNA@counts["ESR1",]
SCTall_fil2_epi$PGR = SCTall_fil2_epi@assays$RNA@counts["PGR",]
SCTall_fil2_epi$ERBB2 = SCTall_fil2_epi@assays$RNA@counts["ERBB2",]
SCTall_fil2_epi$AR = SCTall_fil2_epi@assays$RNA@counts["AR",]


SCTall_fil2_epi@meta.data %>% 
  dplyr::select(patient.2,timepoint,cell.names,ESR1,PGR,ERBB2,AR) %>% 
  gather(-patient.2,-timepoint,-cell.names,key = "receptor", value = "expression") %>% 
  group_by(patient.2,timepoint,receptor) %>% 
  mutate(n = n()) %>% 
  filter(n>=100) %>% 
  group_by(patient.2, timepoint, receptor) %>% 
  summarise(m_expr = mean(expression)) %>% 
  group_by(timepoint, receptor) %>% 
  summarise(m = mean(m_expr), pos = m+sd(m_expr), neg = m-sd(m_expr), n = n()) %>% 
  ggplot() + aes(x=timepoint, y=m, fill=timepoint) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = treatment_colors) +
  geom_linerange(aes(ymin = m, ymax = pos)) +
  facet_wrap(~receptor, scale="free_y") +
  theme_real_minimal() + theme(legend.position = "none") +
  ylab("Mean Expression") + theme(axis.title.x = element_blank()) +
  theme(text = element_text(size=16)) -> recetor_mean_barplot
ggsave(recetor_mean_barplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/recetor_mean_barplot.pdf", height = 4, width = 4, device = cairo_pdf)


SCTall_fil2_epi@meta.data %>% 
  dplyr::select(patient.2,timepoint,cell.names,ESR1,PGR,ERBB2,AR) %>% 
  gather(-patient.2,-timepoint,-cell.names,key = "receptor", value = "expression") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  group_by(patient.2,timepoint,receptor) %>% 
  mutate(n = n()) %>% 
  filter(n>=100) %>% 
  group_by(patient.2, timepoint, receptor) %>% 
  summarise(m_expr = mean(expression)) %>% 
  ggplot() + aes(x = timepoint, y = m_expr, fill=timepoint) + 
  geom_boxplot(outlier.colour = NA) +  
  geom_point() +
  geom_line(aes(group=patient.2)) +
  facet_wrap(~receptor, scale="free_y") + 
  #facet_wrap(~receptor) + 
  theme_real_minimal() + 
  theme(legend.position = "none") +
  scale_fill_manual(values = treatment_colors) +
  ylab("Mean Expression") + theme(axis.title.x = element_blank()) +
  theme(text = element_text(size=16)) -> recetor_mean_boxplot
ggsave(recetor_mean_boxplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/recetor_mean_boxplot.pdf", height = 4, width = 4, device = cairo_pdf)



SCTall_fil2_epi@meta.data %>% 
  dplyr::select(patient.2,timepoint,cell.names,ESR1,PGR,ERBB2,AR) %>% 
  gather(-patient.2,-timepoint,-cell.names,key = "receptor", value = "expression") %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  group_by(patient.2,timepoint,receptor) %>% 
  mutate(n = n()) %>% 
  filter(n>=100) %>% 
  group_by(patient.2, timepoint, receptor) %>% 
  summarise(m_expr = mean(expression)) %>% 
  spread(key = timepoint, value = m_expr, fill=0) %>% 
  mutate(d = post-pre) %>% 
  group_by(receptor) %>% 
  mutate(p_value = t.test(unlist(d))$p.value,
         t_value = t.test(unlist(d))$statistic) %>% 
  arrange(p_value)

##################
###Epi UMAP#######
##################

SCTall_fil2_epi = FindVariableFeatures(SCTall_fil2_epi, assay = "SCT")
SCTall_fil2_epi@assays$SCT@var.features = SCTall_fil2_epi@assays$SCT@var.features[!str_detect(SCTall_fil2_epi@assays$SCT@var.features, "^HB[BA]")]
SCTall_fil2_epi@assays$SCT@var.features = SCTall_fil2_epi@assays$SCT@var.features[!str_detect(SCTall_fil2_epi@assays$SCT@var.features, "^TR[ABGD][VDJ]")]
SCTall_fil2_epi = RunPCA(SCTall_fil2_epi, npcs = 100)
SCTall_fil2_epi = RunUMAP(SCTall_fil2_epi, dims = 1:30, reduction.name = "umap_epi")


###############
###epi umaps###
###############
DimPlot(SCTall_fil2_epi, reduction = "umap_epi", group.by = "patient", cols = patient_colors) -> epi_patient_umap
ggsave(epi_patient_umap, filename = "./plots/fig_pieces_03/epi_patient_umap.pdf", width=6, height=4, device = cairo_pdf)

DimPlot(SCTall_fil2_epi, reduction = "umap_epi", group.by = "timepoint", cols = treatment_colors) -> epi_timepoint_umap
ggsave(epi_timepoint_umap, filename = "./plots/fig_pieces_03/epi_timepoint_umap.pdf", width=6, height=4, device = cairo_pdf)

DimPlot(SCTall_fil2_epi, reduction = "umap_epi", group.by = "ck_pred", cols = ck_colors) -> epi_ck_umap
ggsave(epi_ck_umap, filename = "./plots/fig_pieces_03/epi_ck_umap.pdf", width=6, height=4, device = cairo_pdf)

plot_grid(epi_patient_umap, epi_timepoint_umap, epi_ck_umap, ncol = 1, align = "v", axis = "r") -> epi_umap_combined
ggsave(epi_umap_combined, filename = "./plots/fig_pieces_03/epi_umap_combined.pdf", width=6, height=12, device = cairo_pdf)


#########
#########
#########

SCTall_fil2_epi@meta.data %>% 
  filter(timepoint == "pre") %>% 
  filter(PAM50 != "undetermined") %>% 
  group_by(patient) %>% 
  mutate(n=n()) %>% 
  filter(n>=50) %>% 
  ggplot() + aes(x=patient.2, fill=PAM50) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = PAM50.new_colors) + 
  theme_real_minimal() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12)) -> PAM50_pre_barplot
ggsave(PAM50_pre_barplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/PAM50_pre_barplot.pdf", height = 4, width = 6)


SCTall_fil2_epi@meta.data %>% 
  #filter(timepoint == "pre") %>% 
  filter(PAM50 != "undetermined") %>% 
  group_by(patient.2, timepoint) %>% 
  mutate(n=n()) %>% 
  filter(n>=10) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  ggplot() + aes(x=timepoint, fill=PAM50) + 
  geom_bar(position="fill") + facet_wrap(~patient.2, nrow = 1) +
  scale_fill_manual(values = PAM50.new_colors) + 
  theme_real_minimal() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12)) -> PAM50_prepost_barplot
ggsave(PAM50_prepost_barplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/PAM50_prepost_barplot.pdf", height = 4, width = 6)

####################
###pam50 stat test##
####################

SCTall_fil2_epi@meta.data %>% 
  #filter(timepoint == "pre") %>% 
  filter(PAM50 != "undetermined") %>% 
  group_by(patient.2, timepoint) %>% 
  mutate(n=n()) %>% 
  filter(n>=10) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  group_by(patient.2, timepoint, PAM50) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(patient.2, timepoint) %>% 
  mutate(f = n/sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  spread(key = PAM50, value = f, fill = 0) ->
  PAM50_mn

#PAM50_mn[,c(1,2,5)] %>%  spread(key = timepoint, value = LumA) %>% 
#  mutate(d = post-pre) %>% pull(d) %>% t.test()
#p= 0.0474
#
#PAM50_mn[,c(1,2,5)] %>%  
#spread(key = timepoint, value = LumA) %>% 
#mutate(d = post-pre) %>% 
#pull(d) %>% wilcox.test(exact=F)
#p-value = 0.03461
#
#
#PAM50_mn[,c(1,2,5)] %>% group_by(timepoint) %>%  summarise(m = median(LumA))
#pre  0.600
#post 0.334
#
#

PAM50_mn %>% 
ggplot() + 
  aes(x=timepoint, y = LumA) + 
  geom_point(aes(group=patient.2)) + 
  geom_line(aes(group=patient.2)) + 
  theme_real_minimal() + 
  theme(axis.line.y=element_line(size=0.5)) +
  theme(axis.title.x = element_blank()) + 
  ylab("Luminal A Fraction") + 
  theme(text = element_text(size=16)) +
  scale_x_discrete(expand=c(0.2, 0.2)) +
  scale_y_continuous(limits = c(0,1)) -> 
  lumA_frac
ggsave(lumA_frac, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/lumA_frac.pdf", height = 3, width = 2, device = cairo_pdf)
  


PAM50_mn %>% 
  gather(-patient.2, -timepoint, key=subtype, value = p) %>% 
  spread(key = timepoint, value = p, fill=0) %>% 
  mutate(d = post - pre) %>% 
  ggplot() + aes(x = subtype, y = patient.2, fill = d) + 
  geom_tile() + 
  scale_fill_gradient2(low = "steelblue", mid="white", high ="firebrick") +
  scale_x_discrete(limits = c("LumA","Normal","Her2","Basal","LumB")) +
  scale_y_discrete(limits = rev(c("P08","P02","P05","P06","P03","P07","P01","P04"))) +
  theme_real_minimal() +
  theme(axis.title = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(size=16))



#multiple logistic regression on frequencies, the first one didn't account for the different sample count data  


#mn_test = multinom(data = PAM50_mn, formula = timepoint ~ `Basal` + `Her2` + `LumA` + `LumB`+ `Normal`, maxit=500)
#mn_test = multinom(data = PAM50_mn, formula = timepoint ~ `LumA`, maxit=500)
#summary(mn_test)
#mn_test_z <- summary(mn_test)$coefficients/summary(mn_test)$standard.errors
#mn_test_p <- (1 - pnorm(abs(mn_test_z), 0, 1)) * 2




#####################
###pre_post_heatmap##
#####################
library(ComplexHeatmap)
library(circlize)

tumor_prepost_df = bind_rows(tumor_prepost_list, .id = "patient")

tumor_prepost_df %>% 
  filter(cluster=="post") %>% 
  filter(patient %in% patients[c(1,2,3,5,6,8,10,11)]) %>% 
  filter(avg_logFC >= 0.3) %>% 
  group_by(gene) %>% 
  summarise(n = n(), pct.2 = mean(pct.2)) %>% 
  arrange(-n) %>% 
  filter(n>=5) %>% 
  arrange(pct.2) %>%
  #top_n(20, -pct.2) %>% 
  arrange(pct.2,-n) %>% 
  pull(gene) ->
  tumor_post_avg_top


tumor_prepost_df %>% 
  filter(cluster=="pre") %>% 
  filter(patient %in% patients[c(1,2,3,5,6,8,10,11)]) %>% 
  filter(avg_logFC >= 0.3) %>% 
  group_by(gene) %>% 
  summarise(n = n(), pct.2 = mean(pct.2)) %>% 
  filter(n>=4) %>% 
  arrange(pct.2) %>% 
  #top_n(20, -pct.2) %>% 
  arrange(pct.2,-n) %>% 
  pull(gene) ->
  tumor_pre_avg_top


Idents(SCTall_fil2_epi) = "patient.2"
tumor_prepost_avg = AverageExpression(SCTall_fil2_epi[,SCTall_fil2_epi$patient.2 %in% paste0("P0",1:8)], features = c(tumor_post_avg_top,tumor_pre_avg_top), return.seurat = T, add.ident = "timepoint", assays = "SCT", slot = "data")

tumor_prepost_avg@assays$SCT@data %>% 
  reshape2::melt(c("gene", "sample")) %>% 
  mutate(patient = str_sub(sample, 1, 3),timepoint = str_sub(sample, 5)) %>% 
  dplyr::select(-sample) %>% 
  spread(key=timepoint, value = value) %>% 
  mutate(d = post-pre) %>% 
  dplyr::select(-pre,-post) %>% 
  spread(key=patient, value = d) %>% 
  column_to_rownames(var = "gene") %>% 
  as.matrix() ->
  temp




#col_fun = colorRamp2(c(-1.5,0,0.5,1,1.5), c("black","black","firebrick","orange","yellow"))
col_fun = colorRamp2(c(-1,0,1), c("royalblue4","white","goldenrod"))
Heatmap(temp, row_names_side = "left", column_names_side = "top",
        row_order = c(tumor_post_avg_top,rev(tumor_pre_avg_top)),
        col = col_fun, column_order = colnames(temp)) -> pre_post_tum_heatmap

pdf(file = "~/projects/Breast_cancer_radiation/fig_pieces_03/pre_post_tum_heatmap.pdf", width = 6, height = 6)
draw(pre_post_tum_heatmap)
dev.off()




Idents(SCTall_fil2_tcells) = "patient"
temp = AverageExpression(SCTall_fil2_tcells, features = tumor_post_avg_top, return.seurat = T, add.ident = "timepoint", assays = "SCT", slot = "data")

DoHeatmap(temp, features = tumor_post_avg_top, assay = "SCT", slot = "data", disp.max = 2, lines.width = 1, size = 5, 
          group.colors = patient_colors) + 
  scale_fill_gradientn(colors = c("black","firebrick","orange","yellow"), na.value = "white")


Idents(SCTall_fil2) = "medium_clusters"
TN_dif_genes_list = list()
norm_epi_cells = names(SCTall_fil2$medium_clusters[SCTall_fil2$medium_clusters == "Normal Epithelial"])
for(i in paste0("P", str_pad(1:10, width = 2, side = "left", pad = "0"))){
  temp = SCTall_fil2[,((SCTall_fil2$medium_clusters %in% "Tumor Epithelial") & SCTall_fil2$patient.2 == i) | (SCTall_fil2$cell.names %in% norm_epi_cells)]
  TN_dif_genes = FindAllMarkers(temp, assay = "SCT",logfc.threshold = 0.5, only.pos = T)
  TN_dif_genes$patient.2 = i
  TN_dif_genes_list[[i]] = TN_dif_genes
}

TN_dif_genes_all = FindAllMarkers(SCTall_fil2[,(SCTall_fil2$medium_clusters %in% c("Tumor Epithelial","Normal Epithelial"))], assay = "SCT",logfc.threshold = 0.5, only.pos = T)

TN_dif_genes_all %>% 
  filter(pct.1 >= 0.8) %>% 
  filter(pct.2 <= 0.5) %>% 
  filter(cluster == "Tumor Epithelial") %>% 
  View()


DotPlot(SCTall_fil2_epi, col.min = -1, col.max = 1, cols = c("lightgrey","black"),
        features = c("TFF3","COX6C","XBP1","CST3","H2AFJ",
                     "GATA3","S100A13","KRT18","KRT19","FTL"), 
        group.by = "patient.2") + theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


SCTall_fil2$TN = as.character(SCTall_fil2$medium_clusters)
SCTall_fil2$TN[!SCTall_fil2$TN %in% c("Tumor Epithelial","Normal Epithelial")] = "non-epithelial"
SCTall_fil2$TN = str_replace(SCTall_fil2$TN, "Tumor Epithelial", "Tumor")
SCTall_fil2$TN[str_detect(SCTall_fil2$TN, "Tumor")] = paste0(SCTall_fil2@meta.data[str_detect(SCTall_fil2$TN, "Tumor"),"patient.2"], " Tumor")

DotPlot(SCTall_fil2[,SCTall_fil2$TN != "non-epithelial"],scale = F,
        col.min = 0.1, col.max = 2, cols = c("grey90","black"),
        features = c("TFF3","COX6C","XBP1","CST3","H2AFJ",
                     "GATA3","S100A13","KRT18","KRT19"), 
        group.by = "TN") + theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) ->
  TN_epi_dotplot

ggsave(TN_epi_dotplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/TN_epi_dotplot.pdf", height = 4, width = 6)

# 
# VlnPlot(SCTall_fil2_epi[,SCTall_fil2_epi$patient.2 == "P04"], 
#         group.by = "ck_pred", 
#         features = c("nFeature_RNA","nCount_RNA","percent.mt"), 
#         pt.size = 0, cols = ck_colors) -> dip_aneu_P04_vln
# ggsave(dip_aneu_P04_vln, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/dip_aneu_P04_vln.pdf", height = 3, width = 5)
# 
# 
# 
# temp = SCTall_fil2_epi[,SCTall_fil2_epi$patient.2 == "P04"]
# Idents(temp) = "ck_pred"
# temp = temp[,!is.na(temp$ck_pred)]
# temp = temp[,temp$timepoint == "pre"]
# 
# dip_anu_P04_genes = FindAllMarkers(temp, only.pos = T, logfc.threshold = 0.2, assay = "SCT")
# 
# 
# DotPlot(temp, cols = c("lightgrey","black"),scale = F,
#         features = c("TFF3","COX6C","XBP1","CST3","H2AFJ",
#                      "GATA3","S100A13","KRT18","KRT19","FTL"), 
#         group.by = "ck_pred") + theme(axis.title = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) -> 
#   dip_aneu_P04_dotplot
# ggsave(dip_aneu_P04_dotplot, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/dip_aneu_P04_dotplot.pdf", height = 4, width = 6)
# 
# 

read_csv("~/TENX/analysis/BCR/BCR_4/data/HER2.csv") %>% 
  ggplot() +aes(x = timepoint, y=patient.2,fill=HER2) + geom_tile() + 
  scale_fill_manual(values = receptor_colors[c(3,4,2,1)]) + 
  scale_x_discrete(limits = c("pre","post")) +
  scale_y_discrete(limits = paste0("P", str_pad(as.character(1:19), width = 2, side = "left", pad = "0"))) +
  theme_real_minimal() + theme(axis.title = element_blank()) +
  theme(text = element_text(size=16)) -> her2_heatmap
ggsave(her2_heatmap, filename = "~/projects/Breast_cancer_radiation/fig_pieces_03/her2_heatmap.pdf", height = 3.5, width = 2.5)


