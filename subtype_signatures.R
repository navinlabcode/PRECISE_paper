

library(DropletUtils)
library(Seurat)
library(tidyverse)
library(wesanderson)
library(reshape2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(cowplot) #moo
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
source("~/TENX/analysis/MP/MP_project_4/TCR_functions_script.R")
options(future.globals.maxSize = 60000 * 1024^2)
NewSet1 = function(n){return(colorRampPalette(brewer.pal(9,"Set1")[c(1,5,6,3,2,4,8,7,9)])(n))}
NewSet3 = function(n){return(colorRampPalette(brewer.pal(10,"Set3")[c(4,6,7,1,5,3,8,9)])(n))}



PAM50 = read_csv("./data/PAM50.csv")

PAM50 %>% 
  dplyr::select(group, group.name) %>% 
  unique() %>% 
  data.frame() ->
  PAM50_groups 

for(i in PAM50_groups$group){
  set.seed(77)
  PAM50 %>% 
    filter(group == i) %>% 
    pull(gene) ->
    genes
  SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, 
                   features = list(c(genes)), name=i)
}

set.seed(77)
PAM50 %>% 
  filter(group %in% c('b','c')) %>% 
  pull(gene) ->
  genes
SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, 
                                features = list(c(genes)), name="bc")



SCTall_fil2_epi$a = FALSE
SCTall_fil2_epi$a[SCTall_fil2_epi$a1>0] = TRUE
SCTall_fil2_epi$b = FALSE
SCTall_fil2_epi$b[SCTall_fil2_epi$b1>0] = TRUE
SCTall_fil2_epi$c = FALSE
SCTall_fil2_epi$c[SCTall_fil2_epi$c1>0] = TRUE
SCTall_fil2_epi$d = FALSE
SCTall_fil2_epi$d[SCTall_fil2_epi$d1>0] = TRUE
SCTall_fil2_epi$e = FALSE
SCTall_fil2_epi$e[SCTall_fil2_epi$e1>0] = TRUE
SCTall_fil2_epi$bc = FALSE
SCTall_fil2_epi$bc[SCTall_fil2_epi$bc1>0] = TRUE

SCTall_fil2_epi$PAM50.new = "undetermined"
SCTall_fil2_epi$PAM50.new[SCTall_fil2_epi$a & SCTall_fil2_epi$b & SCTall_fil2_epi$c & !SCTall_fil2_epi$d & !SCTall_fil2_epi$e] = "Basal/Myo"
SCTall_fil2_epi$PAM50.new[!SCTall_fil2_epi$a & SCTall_fil2_epi$c & SCTall_fil2_epi$e] = "Myo/Luminal/Her2"
SCTall_fil2_epi$PAM50.new[!SCTall_fil2_epi$a & SCTall_fil2_epi$b & SCTall_fil2_epi$c & !SCTall_fil2_epi$e] = "Myo/Luminal B"
SCTall_fil2_epi$PAM50.new[!SCTall_fil2_epi$a & !SCTall_fil2_epi$b & SCTall_fil2_epi$c & !SCTall_fil2_epi$e] = "Myo/Luminal A"
SCTall_fil2_epi$PAM50.new[SCTall_fil2_epi$a & SCTall_fil2_epi$e & !SCTall_fil2_epi$d] = "Her2/Basal"
SCTall_fil2_epi$PAM50.new[SCTall_fil2_epi$a & SCTall_fil2_epi$d & !SCTall_fil2_epi$b & !SCTall_fil2_epi$c & !SCTall_fil2_epi$e] = "Basal/Luminal"
SCTall_fil2_epi$PAM50.new[!SCTall_fil2_epi$a & SCTall_fil2_epi$d & !SCTall_fil2_epi$b & !SCTall_fil2_epi$c & !SCTall_fil2_epi$e] = "Luminal"

SCTall_fil2_epi$PAM50 = "undetermined"
SCTall_fil2_epi$PAM50[SCTall_fil2_epi$a & SCTall_fil2_epi$bc & !SCTall_fil2_epi$d & !SCTall_fil2_epi$e] = "Basal"
SCTall_fil2_epi$PAM50[SCTall_fil2_epi$a & !SCTall_fil2_epi$bc & !SCTall_fil2_epi$d & SCTall_fil2_epi$e] = "Her2"
SCTall_fil2_epi$PAM50[!SCTall_fil2_epi$a & !SCTall_fil2_epi$bc & SCTall_fil2_epi$d & !SCTall_fil2_epi$e] = "LumA"
SCTall_fil2_epi$PAM50[SCTall_fil2_epi$a & !SCTall_fil2_epi$bc & SCTall_fil2_epi$d & !SCTall_fil2_epi$e] = "LumB"
SCTall_fil2_epi$PAM50[SCTall_fil2_epi$a & !SCTall_fil2_epi$bc & SCTall_fil2_epi$d & SCTall_fil2_epi$e] = "LumB"
SCTall_fil2_epi$PAM50[!SCTall_fil2_epi$a & SCTall_fil2_epi$bc & SCTall_fil2_epi$d & !SCTall_fil2_epi$e] = "Normal"

PAM50.new_colors = aislynpalette::aislyn_palette("food2", 15)[c(1,4,11,3,6,12,15)]

SCTall_fil2_epi = readRDS("./Rda/SCTall_fil2_epi.rds")

Idents(SCTall_fil2_epi) = "PAM50.new"
DimPlot(SCTall_fil2_epi[,WhichCells(SCTall_fil2_epi, idents = "undetermined", invert = TRUE)], 
        pt.size = 0.5, group.by = "PAM50.new", split.by = "PAM50.new",
        cols = PAM50.new_colors) %>% 
  ggsave(filename = "PAM50.new_UMAP.png", path = "./plots/plots06/", width=16, height=3)


Idents(SCTall_fil2_epi) = "PAM50"
DimPlot(SCTall_fil2_epi[,WhichCells(SCTall_fil2_epi, idents = "undetermined", invert = TRUE)], 
        pt.size = 0.5, group.by = "PAM50", split.by = "PAM50",
        cols = PAM50.new_colors) %>% 
  ggsave(filename = "PAM50_UMAP.png", path = "./plots/plots06/", width=12, height=3)


Idents(SCTall_fil2_epi) = "PAM50.new"
DimPlot(SCTall_fil2_epi[,WhichCells(SCTall_fil2_epi, idents = "undetermined", invert = TRUE)], 
        pt.size = 0.5, group.by = "PAM50.new", split.by = "PAM50.new", ncol = 4,
        cols = PAM50.new_colors)

Idents(SCTall_fil2_epi) = "PAM50.new"
DimPlot(SCTall_fil2_epi[,WhichCells(SCTall_fil2_epi, idents = "undetermined", invert = TRUE)], 
        pt.size = 0.5, group.by = "PAM50.new", cols = PAM50.new_colors) %>% 
  ggsave(filename = "PAM50.new_UMAP_all.png", path = "./plots/plots06/", width=5, height=3)


Idents(SCTall_fil2_epi) = "PAM50"
DimPlot(SCTall_fil2_epi[,WhichCells(SCTall_fil2_epi, idents = "undetermined", invert = TRUE)], 
        pt.size = 0.5, group.by = "PAM50", cols = PAM50.new_colors) %>% 
  ggsave(filename = "PAM50_UMAP_all.png", path = "./plots/plots06/", width=5, height=3)



SCTall_fil2_epi$receptor_status = "ER+PR+HER2-"
SCTall_fil2_epi$receptor_status[SCTall_fil2_epi$patient %in% c("BCR06","BCR17")] = "ER+PR-HER2-"
#SCTall_fil2$receptor_status

VlnPlot(SCTall_fil2_epi, group.by = "PAM50", features = c("ESR1","PGR","ERBB2","MKI67"), cols = c(PAM50.new_colors,"grey40"), ncol = 2) %>% 
ggsave(filename = "PAM50_vln.png", path = "./plots/plots06/", width=8, height=8)

VlnPlot(SCTall_fil2_epi, group.by = "PAM50.new", features = c("ESR1","PGR","ERBB2","MKI67"), cols = c(PAM50.new_colors,"grey40"), ncol = 2) %>% 
ggsave(filename = "PAM50.new_vln.png", path = "./plots/plots06/", width=10, height=8)


SCTall_fil2_epi@meta.data %>% 
  ggplot() + aes(x = timepoint, fill=PAM50) + geom_bar(position = "fill") + 
  facet_wrap(~patient, nrow = 2) + scale_fill_manual(values = c(PAM50.new_colors[1:5],"grey40")) -> temp
  ggsave(temp, filename = "PAM50_bar.png", path = "./plots/plots06/", width=8, height=4)


SCTall_fil2_epi@meta.data %>% 
  filter(!PAM50 == "undetermined") %>% 
  ggplot() + aes(x = timepoint, fill=PAM50) + geom_bar(position = "fill") + 
  facet_wrap(~patient, nrow = 2) + scale_fill_manual(values = c(PAM50.new_colors[1:5])) -> temp
ggsave(temp, filename = "PAM50_bar_det.png", path = "./plots/plots06/", width=8, height=4)

SCTall_fil2_epi@meta.data %>% 
  ggplot() + aes(x = timepoint, fill=PAM50.new) + geom_bar(position = "fill") + 
  facet_wrap(~patient, nrow = 2) + scale_fill_manual(values = c(PAM50.new_colors,"grey40")) ->temp
  ggsave(temp, filename = "PAM50.new_bar.png", path = "./plots/plots06/", width=8, height=4)


SCTall_fil2_epi@meta.data %>% 
  filter(!PAM50.new == "undetermined") %>% 
  ggplot() + aes(x = timepoint, fill=PAM50.new) + geom_bar(position = "fill") + 
  facet_wrap(~patient, nrow = 2) + scale_fill_manual(values = PAM50.new_colors) -> temp
  ggsave(temp, filename = "PAM50.new_bar_det.png", path = "./plots/plots06/", width=8, height=4)

  
SCTall_fil2_epi@meta.data %>% 
  select(a1,b1,c1,d1,e1) %>% 
  melt(variable.name = "gene.group", value.name = "score") %>% 
  ggplot() + aes(x = score, fill = gene.group) + 
  geom_density() + 
  scale_fill_manual(values = NewSet1(11)[c(1,3,5,7,9)]) +
  geom_vline(xintercept = 0) +
  facet_wrap(~gene.group, ncol = 1) + 
  theme(panel.background = element_blank(), legend.position = "none") -> temp
  ggsave(temp, filename = "gene.groups_hist.png", path = "./plots/plots06/", width=4, height=8)

  
    
  
  
  
##################
#####HR status####
##################
  

SCTall_fil2_epi@assays$RNA@counts["ESR1",]
SCTall_fil2_epi@assays$RNA@counts["PGR",]
SCTall_fil2_epi@assays$RNA@counts["ERBB2",]
  

SCTall_fil2_epi$ER = "-"
SCTall_fil2_epi$ER[SCTall_fil2_epi@assays$RNA@counts["ESR1",] > 0] = "+"

SCTall_fil2_epi$PR = "-"
SCTall_fil2_epi$PR[SCTall_fil2_epi@assays$RNA@counts["PGR",] > 0] = "+"

SCTall_fil2_epi$HER2 = "-"
SCTall_fil2_epi$HER2[SCTall_fil2_epi@assays$RNA@counts["ERBB2",] > 0] = "+"


SCTall_fil2_epi$ERPRHER2 = "ER-PR-HER2-"
SCTall_fil2_epi$ERPRHER2[SCTall_fil2_epi$ER == "+" & SCTall_fil2_epi$PR == "+" & SCTall_fil2_epi$HER2 == "+"] = "ER+PR+HER2+"
SCTall_fil2_epi$ERPRHER2[SCTall_fil2_epi$ER == "+" & SCTall_fil2_epi$PR == "+" & SCTall_fil2_epi$HER2 == "-"] = "ER+PR+HER2-"
SCTall_fil2_epi$ERPRHER2[SCTall_fil2_epi$ER == "-" & SCTall_fil2_epi$PR == "+" & SCTall_fil2_epi$HER2 == "+"] = "ER-PR+HER2+"
SCTall_fil2_epi$ERPRHER2[SCTall_fil2_epi$ER == "-" & SCTall_fil2_epi$PR == "+" & SCTall_fil2_epi$HER2 == "-"] = "ER-PR+HER2-"
SCTall_fil2_epi$ERPRHER2[SCTall_fil2_epi$ER == "+" & SCTall_fil2_epi$PR == "-" & SCTall_fil2_epi$HER2 == "+"] = "ER+PR-HER2+"
SCTall_fil2_epi$ERPRHER2[SCTall_fil2_epi$ER == "+" & SCTall_fil2_epi$PR == "-" & SCTall_fil2_epi$HER2 == "-"] = "ER+PR-HER2-"
SCTall_fil2_epi$ERPRHER2[SCTall_fil2_epi$ER == "-" & SCTall_fil2_epi$PR == "-" & SCTall_fil2_epi$HER2 == "+"] = "ER-PR-HER2+"


receptor_colors = c("grey60", "#DB5B5E","#DBD05E", "#db965e", 
                    "#6198DB", "#9e7a9b", "#9eb49c","#b2967d")
                    

SCTall_fil2_epi@meta.data %>% 
  filter(str_detect(minor_clusters2, "^Epithelial")) %>% 
  filter(patient != "BCR17") %>% 
  dplyr::select(minor_clusters2, patient, timepoint, ERPRHER2) %>% 
  ggplot() + aes(x=timepoint, fill=ERPRHER2) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = receptor_colors) +
  facet_wrap(~patient, nrow = 2) + ylab("Percent of Tumor Cells") +
  theme(panel.background = element_blank(),
        axis.title.x= element_blank(), strip.background = element_blank()) -> temp
ggsave(temp, filename = "./plots/plots10/receptor_percent.png", width = 6, height = 4)


SCTall_fil2_epi@meta.data %>% 
  filter(str_detect(minor_clusters2, "^Epithelial")) %>% 
  filter(patient != "BCR17") %>% 
  dplyr::select(minor_clusters2, patient, timepoint, ERPRHER2) %>% 
  mutate(ERPRHER2 = if_else(ERPRHER2 %in% c("ER-PR+HER2+","ER-PR-HER2+","ER-PR+HER2-","ER-PR-HER2-"), "ER-",ERPRHER2)) %>% 
  ggplot() + aes(x=timepoint, fill=ERPRHER2) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = receptor_colors[c(1,5,6,7,8)]) +
  facet_wrap(~patient, nrow = 2) + ylab("Percent of Tumor Cells") +
  theme(panel.background = element_blank(), 
        axis.title.x= element_blank(), strip.background = element_blank()) -> temp
ggsave(temp, filename = "./plots/plots10/ER_percent.png", width = 6, height = 3)




SCTall_fil2_epi@meta.data %>% 
  filter(str_detect(minor_clusters2, "^Epithelial")) %>% 
  filter(patient != "BCR17") %>% 
  dplyr::select(minor_clusters2, patient, timepoint, ERPRHER2) %>% 
  mutate(ERPRHER2 = if_else(ERPRHER2 %in% c("ER+PR-HER2+","ER-PR-HER2+","ER+PR-HER2-","ER-PR-HER2-"), "PR-",ERPRHER2)) %>% 
  mutate(ERPRHER2 = factor(ERPRHER2, levels = c("PR-","ER-PR+HER2-","ER-PR+HER2+","ER+PR+HER2-","ER+PR+HER2+"))) %>% 
  ggplot() + aes(x=timepoint, fill=ERPRHER2) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = receptor_colors[c(1,3,4,7,8)]) +
  facet_wrap(~patient, nrow = 2) + ylab("Percent of Tumor Cells") +
  theme(panel.background = element_blank(), 
        axis.title.x= element_blank(), strip.background = element_blank()) -> temp
ggsave(temp, filename = "./plots/plots10/PR_percent.png", width = 6, height = 3)




SCTall_fil2_epi@meta.data %>% 
  filter(str_detect(minor_clusters2, "^Epithelial")) %>% 
  filter(patient != "BCR17") %>% 
  dplyr::select(minor_clusters2, patient, timepoint, ERPRHER2) %>% 
  mutate(ERPRHER2 = if_else(ERPRHER2 %in% c("ER+PR-HER2-","ER+PR+HER2-","ER-PR+HER2-","ER-PR-HER2-"), "HER2-",ERPRHER2)) %>% 
  mutate(ERPRHER2 = factor(ERPRHER2, levels = c("HER2-","ER-PR-HER2+","ER-PR+HER2+","ER+PR-HER2+","ER+PR+HER2+"))) %>% 
  ggplot() + aes(x=timepoint, fill=ERPRHER2) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = receptor_colors[c(1,2,4,6,8)]) +
  facet_wrap(~patient, nrow = 2) + ylab("Percent of Tumor Cells") +
  theme(panel.background = element_blank(), 
        axis.title.x= element_blank(), strip.background = element_blank()) -> temp
ggsave(temp, filename = "./plots/plots10/HER2_percent.png", width = 6, height = 3)




SCTall_fil2_epi@meta.data %>% 
  filter(timepoint == "pre") %>% 
  ggplot() + aes(x=d1, y=a1, color=patient) + geom_density_2d( bins=10) +
  scale_color_manual(values = patient_colors)

  
SCTall_fil2_epi@meta.data %>% 
  filter(timepoint == "pre") %>% 
  filter(PAM50 != "undetermined") %>% 
  ggplot() + aes(x=patient, fill=PAM50) + geom_bar(position="fill") 
  
#SCTall_fil2_epi$patient.2 = SCTall_fil2$patient.2[colnames(SCTall_fil2_epi)]
  
  
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
  


  
  
