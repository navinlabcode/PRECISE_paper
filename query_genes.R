

library(DropletUtils)
library(Seurat, lib.loc = "~/R/x86_64-redhat-linux-gnu-library/3.6/")
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
library(quantreg)
source("~/TENX/analysis/MP/MP_project_4/TCR_functions_script.R")
options(future.globals.maxSize = 60000 * 1024^2)
NewSet1 = function(n){return(colorRampPalette(brewer.pal(9,"Set1")[c(1,5,6,3,2,4,8,7,9)])(n))}
library(ComplexHeatmap)

VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", split.plot = T, col=treatment_colors,
        pt.size = 0, features = c("rad_resp.1"))


VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("IL18"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("IL1B"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("IL33"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("NLRP3"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TNFRSF1A"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("IL1R1"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TNFRSF1B"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("MYD88"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("IRAK2"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("GSDMD"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("KCNK6"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("CASP4"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("CASP5"))




VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TLR3"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TICAM1"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("JUN"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("NFKB1"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("NFKB2"))

VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TLR4"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("MYD88"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("CXCL8"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("MMP9"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("HIF1A"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("NFKBIA"))

VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TLR2"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TLR3"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TLR1"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("TLR4"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, features = c("IFNG"))
VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, 
        features = c(paste0("TLR",1:10)), ncol = 5)


VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, 
        features = c("CXCL8"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, 
        features = c("IL13"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("sct_IFNG"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, 
        features = c("CD80"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL1R1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL1RAP"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("NFKBIA"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL1RL1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL4"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL5"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("EOMES"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("TBX21"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("RORC"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("GATA3"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("MX1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("OAS1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IRF1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("SIN3A"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("ISG15"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("STAT2"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("STAT1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("STAT3"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IFNGR1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("TNF"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IFNGR2"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IFNAR1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IFNAR2"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IRF1"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL4","IL5","IL9","IL13"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("sct_IFNG"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("CD40LG"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("FASLG"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
         pt.size = 0, assay = "SCT",slot = "data",
        features = c("TNF","sct_IFNG","IL2"))

VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("ITGAE","EOMES","TIGIT","ENTPD1","CTLA4","PDCD1","TOX"))

VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("LAYN","TOX","ID2","HIF1A","ZEB2","TBX21"))

VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("GZMA","GZMK","GNLY","GZMB"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("GATA3","IL4","IL10","RORGT"))



VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("CD69"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("CD69"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL2RA"))
VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", 
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("LAG3"))



VlnPlot(SCTall_fil2_epi, group.by = "patient", split.plot = T,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("S.Score","G2M.Score") ,cols = treatment_colors)













VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL5"))
VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL4"))
VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL13"))

VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL1B"))
VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL18"))
VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("IL33"))

VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("TGFB1"))
VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("RNASEH2A"))


VlnPlot(SCTall_fil2, group.by = "patient", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("CD68"))



VlnPlot(SCTall_fil2, group.by = "major_clusters", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("HLA-A"))

SCTall_fil2 = AddModuleScore(SCTall_fil2, features = ISGs, name = "ISGs")

VlnPlot(SCTall_fil2, group.by = "major_clusters", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("ISGs1"))

VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("ISGs1"))

VlnPlot(SCTall_fil2, group.by = "minor_clusters2", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("LIG1"))


#saveRDS(SCTall_fil2_epi, file = "./Rda/tumor_RNA.rds", compress = FALSE)



VlnPlot(SCTall_fil2_tcells, group.by = "patient", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("CCL4L2","CCL4"), 
        split.plot = TRUE, cols = rev(treatment_colors))


VlnPlot(SCTall_fil2_epi, group.by = "patient", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("NR4A1","ZFP36"), 
        split.plot = TRUE, cols = rev(treatment_colors))



VlnPlot(SCTall_fil2_epi, group.by = "patient", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = BCR15_epi_post_markers_top, ncol = 5,
        split.plot = TRUE, cols = rev(treatment_colors))




VlnPlot(SCTall_fil2_epi, group.by = "patient", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("PRKCZ"), 
        split.plot = TRUE, cols = treatment_colors)

VlnPlot(SCTall_fil2, group.by = "major_clusters", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("PRKCZ"), 
        split.plot = TRUE, cols = treatment_colors)


VlnPlot(SCTall_fil2_epi, group.by = "patient", ncol = 2,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("HLA-A","HLA-B","HLA-C","B2M"), 
        split.plot = TRUE, cols = treatment_colors)

VlnPlot(SCTall_fil2_myeloid, group.by = "patient", ncol = 2,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("HLA-A","HLA-B","HLA-C","B2M"), 
        split.plot = TRUE, cols = treatment_colors)




VlnPlot(SCTall_fil2_epi, group.by = "patient", ncol = 2,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("HLA-A","HLA-B","HLA-C","B2M"), 
        split.plot = TRUE, cols = treatment_colors)


VlnPlot(SCTall_fil2, group.by = "major_clusters", ncol = 2,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("HLA-A","HLA-B","HLA-C","B2M"), 
        split.plot = TRUE, cols = treatment_colors)



SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(c("HLA-A","HLA-B","HLA-C","B2M")), name = "HLA.")


VlnPlot(SCTall_fil2_epi, group.by = "patient", ncol = 1,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = "HLA.1", 
        split.plot = TRUE, cols = treatment_colors)


VlnPlot(SCTall_fil2_epi, group.by = "patient", ncol = 2,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("CD274","CDK4","CDK6"), 
        split.plot = TRUE, cols = treatment_colors)

VlnPlot(SCTall_fil2, group.by = "major_clusters", ncol = 2,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = c("CD274","CDK4","CDK6"), 
        split.plot = TRUE, cols = treatment_colors)




checkpoint_genes = c("TNFRSF14","TNFSF18","CD80","CD86","PVR","CD274","PDCD1LG2","LGALS9","CD40","TNFSF4","TNFSF9","ICOSLG","CD70")
checkpoint_tcells = c("BTLA","TNFRSF18","CTLA4","TIGIT","PDCD1","HAVCR2","LAG3","CD28","CD40LG","TNFRSF4","TNFRSF9","ICOS","CD27")



VlnPlot(SCTall_fil2_epi, group.by = "patient", ncol = 4,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = checkpoint_genes, 
        split.plot = TRUE, cols = treatment_colors)



VlnPlot(SCTall_fil2_epi, group.by = "patient", ncol = 4,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = "PDCD1", 
        split.plot = TRUE, cols = treatment_colors)






VlnPlot(SCTall_fil2_myeloid, group.by = "patient", ncol = 4,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = checkpoint_genes, 
        split.plot = TRUE, cols = treatment_colors)



VlnPlot(SCTall_fil2_myeloid, group.by = "minor_clusters2", ncol = 2,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = checkpoint_genes[c(1,4,8,9)], 
        split.plot = TRUE, cols = treatment_colors)




VlnPlot(SCTall_fil2_tcells, group.by = "patient", ncol = 4,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = checkpoint_tcells, 
        split.plot = TRUE, cols = treatment_colors)





VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4", ncol = 4,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = checkpoint_tcells, 
        split.plot = TRUE, cols = treatment_colors)


VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4",
         pt.size = 0, assay = "SCT",slot = "data",
        features = "BCL6")




tcga_bc_genes = pull(read_tsv("/volumes/seq/database/CancerGenes/TCGA_breastCancer_genelist_45.txt", col_names = "gene"))


VlnPlot(SCTall_fil2[,SCTall_fil2$minor_clusters2 %in% tumor_clusters], group.by = "patient", ncol = 10,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = tcga_bc_genes, 
        split.plot = TRUE, cols = treatment_colors) -> temp
ggsave(temp, filename = "~/TENX/analysis/BCR/BCR_4/plots/plots11/BC_tcga_genes_vln.png", width = 30, height = 15)


BC_bed_genes = pull(read_tsv("~/code/resources/gene_lists/BC_genes.sorted.bed", col_names = F), "X4")



VlnPlot(SCTall_fil2[,SCTall_fil2$minor_clusters2 %in% tumor_clusters], group.by = "patient", ncol = 10,
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = BC_bed_genes, 
        split.plot = TRUE, cols = treatment_colors) -> temp
ggsave(temp, filename = "~/TENX/analysis/BCR/BCR_4/plots/plots11/BC_bed_genes_vln.png", width = 30, height = 25)


VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4",
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("ENTPD1","PDCD1","ITGAE"))

VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4",
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("CX3CR1","FGFBP2","FCGR3A"))

VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4",
        pt.size = 0, assay = "SCT",slot = "data",
        features = c("sct_IFNG","TNF"))


SCTall_fil2_tcells$PD1 = c("-","+")[as.numeric(SCTall_fil2_tcells@assays$RNA@counts["PDCD1",] > 0) + 1]

SCTall_fil2_tcells@meta.data %>% 
        ggplot() + aes(x=patient, fill=PD1) + 
        geom_bar(position="fill") + theme(axis.text.x = element_text(angle = 90)) + 
        facet_wrap(~timepoint)



VlnPlot(SCTall_fil2_tcells, group.by = "patient",split.by = "timepoint",
        pt.size = 0, assay = "SCT",slot = "data",split.plot = T,
        features = c("ENTPD1","PDCD1","ITGAE"))

VlnPlot(SCTall_fil2_tcells, group.by = "patient",split.by = "timepoint",
        pt.size = 0, assay = "SCT",slot = "data",split.plot = T,
        features = c("CX3CR1","FGFBP2","FCGR3A"))


VlnPlot(SCTall_fil2_tcells, group.by = "minor_clusters4",
        pt.size = 0, assay = "SCT",slot = "data",
        features = "BCL6")

pro_apop_genes = c("BCL2L11","BID","BAD","BMF","BBC3","PMAIP1","BAK1","BAX","BOK")
anti_apop_genes = c("MCL1","BCL2","BCL2L1","BCL2L2","BCL2A1")


VlnPlot(SCTall_fil2[,SCTall_fil2$minor_clusters2 %in% tumor_clusters], group.by = "patient", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = pro_apop_genes, 
        split.plot = TRUE, cols = treatment_colors)



VlnPlot(SCTall_fil2[,SCTall_fil2$minor_clusters2 %in% tumor_clusters], group.by = "patient", 
        split.by = "timepoint", pt.size = 0, assay = "SCT",slot = "data",
        features = anti_apop_genes, 
        split.plot = TRUE, cols = treatment_colors)

VlnPlot(SCTall_fil2, group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data",
        features = pro_apop_genes)

VlnPlot(SCTall_fil2, group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data",split.by = "timepoint",
        features = "GLB1")



SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(pro_apop_genes), name = "pro.apop.")
SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(anti_apop_genes), name = "anti.apop.")


VlnPlot(SCTall_fil2[,SCTall_fil2$minor_clusters2 %in% tumor_clusters], 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "pro.apop.1") + NoLegend()

VlnPlot(SCTall_fil2[,SCTall_fil2$minor_clusters2 %in% tumor_clusters], 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "anti.apop.1") + NoLegend()



SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(pro_apop_genes), name = "pro.apop.")
SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(anti_apop_genes), name = "anti.apop.")

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "pro.apop.1") + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "anti.apop.1") + NoLegend()


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = anti_apop_genes) + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = pro_apop_genes) + NoLegend()

VlnPlot(SCTall_fil2, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "MCL1") + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "NFKB1") + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "PITPNM3") + NoLegend()


atm, atr, dnapcs 

SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(pro_apop_genes), name = "pro.apop.")



DAS_genes = c("CDKN2A","HMGA1","HMGA2","HIRA","ASF1A","H2AFY","H2AFY2","PML","MMP3","TP53BP1","CDKN1A")

SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(DAS_genes), name = "pro.DAS.")


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = DAS_genes) + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "pro.DAS.1") + NoLegend()


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "SLCO2B1") + NoLegend()


HLA_genes = sort(row.names(SCTall_fil2_macs)[str_detect(row.names(SCTall_fil2_macs), "^HLA-")])
HLA_2_genes = c(HLA_genes[str_detect(HLA_genes, "^HLA-[D]")],"CD74")

VlnPlot(SCTall_fil2_macs, 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = HLA_2_genes) + NoLegend()


VlnPlot(SCTall_fil2[,SCTall_fil2$major_clusters == "Myeloid"], 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = HLA_2_genes) + NoLegend()


VlnPlot(SCTall_fil2, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "CD68") + NoLegend()


VlnPlot(SCTall_fil2, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "CD68") + NoLegend()



SCTall_fil2_macs = AddModuleScore(SCTall_fil2_macs, features = list(HLA_2_genes[c(1,2,3,5,6,7,9,12,13,14,15)]), name = "HLA.II.", nbin = 15)

VlnPlot(SCTall_fil2_macs, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "HLA.II.1") + NoLegend()

for(i in patients){
        pre = SCTall_fil2_macs$HLA.II.1[SCTall_fil2_macs$patient == i & SCTall_fil2_macs$timepoint == "pre"]
        post = SCTall_fil2_macs$HLA.II.1[SCTall_fil2_macs$patient == i & SCTall_fil2_macs$timepoint == "post"]
        data.frame(score = pre, timepoint="pre") %>% 
                rbind( data.frame(score = post, timepoint="post")) %>% 
                ggplot(aes(x=score, color = timepoint, fill = NULL)) + 
                geom_boxplot() + coord_flip() + theme_minimal() + ggtitle(i) -> g
        print(g)
        p = wilcox.test(pre, post)$p.value
        print(i)
        print(p)
}


SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(HLA_2_genes[c(1,2,3,5,6,7,9,12,13,14,15)]), name = "HLA.II.", nbin = 15)

VlnPlot(SCTall_fil2, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "HLA.II.1") + NoLegend()

for(i in patients){
        pre = SCTall_fil2$HLA.II.1[SCTall_fil2$patient == i & SCTall_fil2$timepoint == "pre" & str_detect(SCTall_fil2$minor_clusters2 , "Macro|Dend")]
        post = SCTall_fil2$HLA.II.1[SCTall_fil2$patient == i & SCTall_fil2$timepoint == "post" & str_detect(SCTall_fil2$minor_clusters2 , "Macro|Dend")]
        print(g)
        p = wilcox.test(pre, post)$p.value
        print(i)
        print(p)
}




VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = c( "CEP55", "MCM7", "CDC45", "MCM5", "KIF18B"," CIT", "EZH2")) + NoLegend()


 

SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(sgs), name = "sen.")


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "sen.1") + NoLegend()


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = sgs[1:16]) + NoLegend()




C2_df = msigdbr(species = "Homo sapiens", category = "C2")

C2_df %>% 
        filter(gs_id %in% c("M9143")) %>% 
        filter(gene_symbol %in% row.names(SCTall_fil2_epi)) %>% 
        pull(gene_symbol) ->
        temp

temp = temp[rowSums(SCTall_fil2_epi[temp,]) > 5000]


SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(temp), name = "M9143.")


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "M9143.1") + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = temp[1:16]) + NoLegend()

so = SCTall_fil2_epi
gsea_df = C2_df
path.id = "M27188"
n=5000
module_plot = function(so, gsea_df, path.id, n, clus = "minor_clusters2", split = "timepoint"){
        gsea_df %>% 
                filter(gs_id %in% path.id) %>% 
                filter(gene_symbol %in% row.names(so)) %>% 
                pull(gene_symbol) ->
                temp
        temp = temp[rowSums(so[temp,]) >= n]
        so = AddModuleScore(so, features = list(temp), name = paste0(path.id,"."))
        plotlist = list()
        p = VlnPlot(so, 
                group.by = clus,
                pt.size = 0, assay = "SCT",slot = "data", 
                split.by = split,split.plot = T,
                features = paste0(path.id,".1")) + NoLegend()
        plotlist[["score"]] = p
        for(i in temp){
                plotlist[[i]] = VlnPlot(so, 
                        group.by = clus,
                        pt.size = 0, assay = "SCT",slot = "data", 
                        split.by = split,split.plot = T,
                        features = i) + NoLegend()    
                
        }
        
        return(plotlist)
}

genes = c("PVRL4","GPR172B","DAO","CCDC74B")
custom_module_plot = function(so, genes, n, clus = "minor_clusters2", split = "timepoint"){
        print(genes[!genes %in% row.names(so)])
        temp = genes[genes %in% row.names(so)] 
        temp = temp[rowSums(so[temp,]) >= n]
        so = AddModuleScore(so, features = list(temp), name = "score.")
        plotlist = list()
        p = VlnPlot(so, 
                    group.by = clus,
                    pt.size = 0, assay = "SCT",slot = "data", 
                    split.by = split,split.plot = T,
                    features = "score.1") + NoLegend()
        plotlist[["score"]] = p
        for(i in temp){
                plotlist[[i]] = VlnPlot(so, 
                                        group.by = clus,
                                        pt.size = 0, assay = "SCT",slot = "data", 
                                        split.by = split,split.plot = T,
                                        features = i) + NoLegend()    
                
        }
        
        return(plotlist)
}



vln_sen_M27188 = module_plot(SCTall_fil2_epi, gsea_df, "M27188",5000)
vln_sen_M9143 = module_plot(SCTall_fil2_epi, gsea_df, "M9143",5000)
vln_sen_M27189 = module_plot(SCTall_fil2_epi, gsea_df, "M27189",5000)
vln_sen_M27191 = module_plot(SCTall_fil2_epi, gsea_df, "M27191",5000)
vln_sen_M27190 = module_plot(SCTall_fil2_epi, gsea_df, "M27190",5000)
vln_sen_M27187 = module_plot(SCTall_fil2_epi, gsea_df, "M27187",1000)


nagano_1 = c("NECTIN4","SLC52A1","DAO","CCDC74B","LOXL4","EVL","PRODH","E2F7","LY6D","IGFBP2","CRABP2","EPN3","APOBEC3B","IER5","ANGPTL2","SLC48A1","METTL27","E2F2","NXPH4","PPM1D","CDKN1A","BTG2","SULF2")
vln_nagano_1 = custom_module_plot(SCTall_fil2_epi, nagano_1,500)

nagano_3 = c("NECTIN4","PRODH","LY6D","DAO","EPN3","SLC52A1")
nagano_2 = c("NECTIN4","CDKN1A","IGFBP2","E2F7","SLC48A1","APOBEC3B","IER5","E2F2","LOXL4")

vln_nagano_2 = custom_module_plot(SCTall_fil2_epi, nagano_2,500)
vln_nagano_3 = custom_module_plot(SCTall_fil2_epi, nagano_3,100)


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "CDKN2B") + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        features = "FOXM1") + NoLegend()

DotPlot(SCTall_fil2_epi, assay = "SCT", features = c("FOXM1","FOS","JUN"), split.by = "timepoint", cols = c("red","blue"))


tumor_prepost_df %>% 
        filter(gene %in% tfs) %>% 
        arrange(patient, -avg_logFC) %>% 
        filter(cluster=="post") %>% 
        filter(avg_logFC > 0) %>% 
        mutate(pct.dif = pct.1-pct.2) %>% 
        group_by(gene) %>% 
        summarise(m = mean(avg_logFC), pct = mean(pct.dif)) %>% View()


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",adjust = 1,
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = c("PIK3CA")) + NoLegend()




VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",adjust = 1,
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = c("PTK2","ERBB2")) + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",adjust = 1,
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = c("ESR1")) + NoLegend()


mtor_genes = c("HIF1A", "NFE2L2", "HSF1", "NFAT5", "TP53", "FOXO3")

mtor_genes %in% row.names(SCTall_fil2_epi)

SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(mtor_genes), name = "mtor.")

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = mtor_genes) + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = "mtor.1") + NoLegend()


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = mtor_genes) + NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = tumor_post_avg_top[1:6]) + NoLegend()


FOXO3A_targets = read_tsv("./data/gene_lists/FOXO3A_motifs.txt", skip = 2, col_names = "gene") %>%  pull(gene)

FOXO3A_targets = FOXO3A_targets[FOXO3A_targets %in% row.names(SCTall_fil2_epi)]

temp = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[FOXO3A_targets,]), decreasing = T)
#FOXO3A_targets_high = names(temp[temp >= 5000])
FOXO3A_targets_high = FOXO3A_targets[FOXO3A_targets %in% SCTall_fil2_epi@assays$SCT@var.features]


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = FOXO3A_targets_high[1:6]) + NoLegend()

SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(FOXO3A_targets_high), name = "foxo3.top.")

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = "foxo3.top.1") + NoLegend()




p53_targets = read_tsv("./data/gene_lists/p53_motifs.txt", skip = 2, col_names = "gene") %>%  pull(gene)

p53_targets = p53_targets[p53_targets %in% row.names(SCTall_fil2_epi)]

temp = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[p53_targets,]), decreasing = T)
#p53_targets_high = names(temp[temp >= 5000])
p53_targets_high = p53_targets[p53_targets %in% SCTall_fil2_epi@assays$SCT@var.features]



SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(p53_targets_high), name = "p53.top.")

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = "p53.top.1") + NoLegend()




mdm2_genes = read_tsv("./data/gene_lists/mdm2_genes.txt", skip = 2, col_names = "gene") %>%  pull(gene)

mdm2_genes = mdm2_genes[mdm2_genes %in% row.names(SCTall_fil2_epi)]

temp = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[mdm2_genes,]), decreasing = T)
#mdm2_genes_high = names(temp[temp >= 5000])
mdm2_genes_high = mdm2_genes[mdm2_genes %in% SCTall_fil2_epi@assays$SCT@var.features]


SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(mdm2_genes_high), name = "mdm2.top.")

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = "mdm2.top.1") + NoLegend()





p300_targets = read_tsv("./data/gene_lists/p300_motifs.txt", skip = 2, col_names = "gene") %>%  pull(gene)

p300_targets = p300_targets[p300_targets %in% row.names(SCTall_fil2_epi)]

temp = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[p300_targets,]), decreasing = T)
#p300_targets_high = names(temp[temp >= 5000])
p300_targets_high = p300_targets[p300_targets %in% SCTall_fil2_epi@assays$SCT@var.features]


SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(p300_targets_high), name = "p300.top.")

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = c("p300.top.1","mdm2.top.1","p53.top.1","foxo3.top.1","mtor.1")) + 
        NoLegend()


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = c("PARP1","NFKBIA","TNF")) + 
        NoLegend()


NFKB_targets = read_tsv("./data/gene_lists/NFKB_targets.txt", skip = 2, col_names = "gene") %>%  pull(gene)

NFKB_targets = NFKB_targets[NFKB_targets %in% row.names(SCTall_fil2_epi)]

temp = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[NFKB_targets,]), decreasing = T)
#NFKB_targets_high = names(temp[temp >= 5000])
NFKB_targets_high = NFKB_targets[NFKB_targets %in% SCTall_fil2_epi@assays$SCT@var.features]

SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(NFKB_targets_high), name = "nfkb.top.")

VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,cols = treatment_colors,
        features = c("p300.top.1","mdm2.top.1","p53.top.1","foxo3.top.1","mtor.1","nfkb.top.1")) + 
        NoLegend()



m_df %>% 
        filter(gs_name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_P53_PATHWAY","HALLMARK_HYPOXIA",
                              "HALLMARK_UV_RESPONSE_UP","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")) %>% 
        group_by(gene_symbol) %>% 
        summarise(n=n()) %>% View()
        filter(n <= 2) %>% 
        pull(gene_symbol) ->
        H_spec


m_df %>% 
filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>% 
        pull(gene_symbol) -> temp

temp = temp[temp %in% row.names(SCTall_fil2_epi)]

temp2 = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[temp,]), decreasing = T)
temp3 = names(temp2[temp2 >= 5000])
temp4 = temp3[temp3 %in% tumor_post_avg_top]


#VlnPlot(SCTall_fil2_epi, 
#        group.by = "minor_clusters2",
#        pt.size = 0, assay = "SCT",slot = "data", 
#        split.by = "timepoint",split.plot = T,
#        cols = treatment_colors,
#        features = temp3[11:20]) + 
#        NoLegend()


#HALLMARK_TNFA_SIGNALING_VIA_NFKB
VlnPlot(SCTall_fil2_epi, 
        group.by = "orig.ident",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,
        features = c("IER2","EGR1","NFKBIA")) + 
        NoLegend()




m_df %>% 
        filter(gs_name == "HALLMARK_P53_PATHWAY") %>% 
        pull(gene_symbol) -> temp

temp = temp[temp %in% row.names(SCTall_fil2_epi)]

temp2 = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[temp,]), decreasing = T)
temp3 = names(temp2[temp2 >= 3000])
temp4 = temp3[temp3 %in% tumor_post_avg_top]
temp4[temp4 %in% H_spec]



m_df %>% 
        filter(gs_name == "HALLMARK_HYPOXIA") %>% 
        pull(gene_symbol) -> temp

temp = temp[temp %in% row.names(SCTall_fil2_epi)]

temp2 = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[temp,]), decreasing = T)
temp3 = names(temp2[temp2 >= 3000])
temp4 = temp3[temp3 %in% tumor_post_avg_top]
temp4[temp4 %in% H_spec]


#HALLMARK_HYPOXIA
VlnPlot(SCTall_fil2_epi, 
        group.by = "orig.ident",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,
        features = c("DUSP1","ZFP36","KLF6")) + 
        NoLegend()


m_df %>% 
        filter(gs_name == "HALLMARK_UV_RESPONSE_UP") %>% 
        pull(gene_symbol) -> temp

temp = temp[temp %in% row.names(SCTall_fil2_epi)]

temp2 = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[temp,]), decreasing = T)
temp3 = names(temp2[temp2 >= 3000])
temp4 = temp3[temp3 %in% tumor_post_avg_top]
temp4[temp4 %in% H_spec]


"HALLMARK_UV_RESPONSE_UP"
VlnPlot(SCTall_fil2_epi, 
        group.by = "orig.ident",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,
        features = c("JUNB","FOSB","NR4A1")) + 
        NoLegend()



#"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"



m_df %>% 
        filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% 
        pull(gene_symbol) -> temp

temp = temp[temp %in% row.names(SCTall_fil2_epi)]

temp2 = sort(rowSums(SCTall_fil2_epi@assays$RNA@counts[temp,]), decreasing = T)
temp3 = names(temp2[temp2 >= 3000])
temp4 = temp3[temp3 %in% tumor_post_avg_top2]
temp4[temp4 %in% H_spec]
temp3[temp3 %in% H_spec]

emt_df = read_csv("./data/gene_lists/emt.csv")
emt_df %>% filter(dir=="up") %>% pull(gene) -> emt_up
emt_df %>% filter(dir=="down") %>% pull(gene) -> emt_down

 

tumor_prepost_df %>% 
        filter(cluster=="post") %>% 
        filter(avg_logFC >= 0.1) %>% 
        group_by(gene) %>% 
        summarise(n = n()) %>% 
        arrange(-n) %>% 
        filter(n>=5) %>% 
        #arrange(pct.2) %>%
        #top_n(20, -pct.2) %>% 
        #arrange(-n, pct.2) %>% 
        pull(gene) ->
        tumor_post_avg_top2


#"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
VlnPlot(SCTall_fil2_epi, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,
        features = c("VIM","EPCAM","GADD45B")) + 
        NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "orig.ident",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,
        features = c("CDH1","CTNNB1","OCLN","CLDN1","CDH2","GADD45B","VIM","FN1","S100A4")) + 
        NoLegend()




VlnPlot(SCTall_fil2_epi[,cells], 
        group.by = "orig.ident",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,adjust = 1,
        features = c("SAT1","TPM1","GADD45B")) + 
        NoLegend()

VlnPlot(SCTall_fil2_epi[,cells], 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 4,
        features = temp3[temp3 %in% emt_up]) + 
        NoLegend()

cells = c()
n = 300
set.seed(77)
for(i in patients[-9]){
        pre = row.names(SCTall_fil2_epi@meta.data)[SCTall_fil2_epi$patient == i & SCTall_fil2_epi$timepoint == "pre"]
        post = row.names(SCTall_fil2_epi@meta.data)[SCTall_fil2_epi$patient == i & SCTall_fil2_epi$timepoint == "post"]

        if(length(pre)>=n){
                cells = c(cells, sample(pre, size = n))
        }
        else{
                cells = c(cells, pre)
        }
        if(length(post)>=n){
                cells = c(cells, sample(post, size = n))
        }
        else{
                cells = c(cells, post)
        }
        
}

Idents(SCTall_fil2_epi) = "timepoint"
test = FindAllMarkers(SCTall_fil2_epi[,cells], assay = "SCT", logfc.threshold = 0)

test %>% 
filter(cluster == "post") %>% 
arrange(-avg_logFC) %>% 
top_n(25, avg_logFC) %>% 
        pull(gene) -> test2

VlnPlot(SCTall_fil2_epi[,cells], 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 5,
        features = test2) + 
        NoLegend()

SCTall_fil2_epi = AddModuleScore(SCTall_fil2_epi, features = list(c("HLA-A","HLA-B","HLA-C")), name = "HLA1.")

VlnPlot(SCTall_fil2_epi, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 4,
        features = c("HLA-A","HLA-B","HLA-C")) + 
        NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 1,
        features = c("HLA1.1")) + 
        NoLegend()


rad_genes = read_csv(file = "~/TENX/analysis/BCR/BCR_4/data/gene_lists/rad_genes.csv")


c("CRABP2","EIF3L","H2AFZ","LEF1","RUNX1","JUN","STAT1","SUMO1","IRF1")

VlnPlot(SCTall_fil2, 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 3,
        features = c("CRABP2","EIF3L","H2AFZ","LEF1","RUNX1","JUN","STAT1","SUMO1","IRF1")) + 
        NoLegend()


VlnPlot(SCTall_fil2, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 3,
        features = c("CRABP2","EIF3L","H2AFZ","LEF1","RUNX1","JUN","STAT1","SUMO1","IRF1")) + 
        NoLegend()


rad_genes %>% 
        filter(signature == "ARCTIC") %>% 
        filter(dir == "pos") %>% 
        pull(gene) -> arctic_pos


rad_genes %>% 
        filter(signature == "ARCTIC") %>% 
        filter(dir == "neg") %>% 
        pull(gene) -> arctic_neg

SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(arctic_pos), name = "arctic.pos.")
SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(arctic_neg), name = "arctic.neg.")
SCTall_fil2$arctic = SCTall_fil2$arctic.pos.1 - SCTall_fil2$arctic.neg.1



VlnPlot(SCTall_fil2, 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 1,
        features = c("arctic")) + 
        NoLegend()


VlnPlot(SCTall_fil2, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 1,
        features = c("arctic")) + 
        NoLegend()


VlnPlot(SCTall_fil2, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors,ncol = 1,
        features = c("arctic")) + 
        NoLegend()


VlnPlot(SCTall_fil2[,SCTall_fil2$timepoint=="pre"], 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        cols = patient_colors,ncol = 1,
        features = c("arctic")) + 
        NoLegend()



m_df %>% 
        filter(gs_id == "M9221") %>% 
        pull(gene_symbol) ->
        brwn_ifn_genes

SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(brwn_ifn_genes), name = "b.ifn.")


SCTall_fil2@meta.data %>% 
        ggplot() + aes(x=b.ifn.1, color=patient) + 
        geom_density() + scale_color_manual(values = patient_colors) +
        theme_minimal()
        
SCTall_fil2@meta.data %>% 
        #filter(timepoint=="post") %>% 
        ggplot() + aes(x=b.ifn.1, color=major_clusters) + 
        geom_density() + scale_color_manual(values = major_cluster_colors) +
        theme_minimal() + facet_wrap(~patient, nrow = 3)



rbind(SCTall_fil2_tcells_minor_clusters4_markers_CD4, SCTall_fil2_tcells_minor_clusters4_markers_CD8) %>% 
        filter(str_detect(cluster, "ISG")) %>% 
        filter(avg_logFC > 0.5) %>%
        arrange(-avg_logFC) %>%
        group_by(gene) %>% 
        mutate(n = n()) %>% 
        filter(n >1) %>% 
        pull(gene) %>% 
        unique() ->
        ifn_genes2


rbind(SCTall_fil2_tcells_minor_clusters4_markers_CD4, SCTall_fil2_tcells_minor_clusters4_markers_CD8) %>% 
        filter(str_detect(cluster, "ISG")) %>% 
        #filter(avg_logFC > 0.5) %>%
        arrange(-avg_logFC) %>%
        group_by(gene) %>% 
        mutate(n = n()) %>% 
        filter(n < 2) %>% View()
        pull(gene) %>% 
        unique() ->
        ifn_genes2


SCTall_fil2_tcells_minor_clusters4_markers %>% 
        filter(str_detect(cluster, "ISG")) %>% 
        filter(avg_logFC > 0.5) %>% 
        group_by(gene) %>% 
        mutate(n = n()) %>% 
        filter(n > 1) %>%
        arrange(-avg_logFC) %>% 
        pull(gene) %>% 
        unique() ->
        ifn_genes


SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(ifn_genes), name = "bcr.ifn.")

SCTall_fil2@meta.data %>% 
        ggplot() + aes(x=bcr.ifn.1, color=patient) + 
        geom_density() + scale_color_manual(values = patient_colors) +
        theme_minimal() + facet_wrap(~timepoint)

SCTall_fil2@meta.data %>% 
        #filter(timepoint=="post") %>% 
        ggplot() + aes(x=bcr.ifn.1, color=major_clusters) + 
        geom_density() + scale_color_manual(values = major_cluster_colors) +
        theme_minimal() + facet_wrap(~timepoint)

SCTall_fil2@meta.data %>% 
        #filter(timepoint=="post") %>% 
        ggplot() + aes(x=bcr.ifn.1, color=timepoint) + 
        geom_density() + scale_color_manual(values = treatment_colors) +
        theme_minimal() + facet_wrap(~major_clusters)


SCTall_fil2@meta.data %>% 
        #filter(major_clusters=="T - Cells") %>% 
        #filter(major_clusters=="NK Cells") %>% 
        filter(major_clusters=="Myeloid") %>%
        ggplot() + aes(x=bcr.ifn.1, color=timepoint) + 
        geom_density() + scale_color_manual(values = treatment_colors) +
        theme_minimal() + facet_wrap(~patient)


SCTall_fil2@meta.data %>% 
        filter(bcr.ifn.1 > 0.25) %>% 
        ggplot() + aes(x=patient, fill = major_clusters) + 
        geom_bar() + facet_wrap(~timepoint)


VlnPlot(SCTall_fil2, 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 5,
        features = checkpoint_genes) + 
        NoLegend()

VlnPlot(SCTall_fil2, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 5,
        features = checkpoint_genes) + 
        NoLegend()


#tnfrsf14
#lgals9

VlnPlot(SCTall_fil2[,SCTall_fil2$patient == "BCR06"], 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 5,
        features = c("LGALS9","TNFRSF14")) + 
        NoLegend()

VlnPlot(SCTall_fil2, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        features = c("HAVCR2")) + 
        NoLegend()
VlnPlot(SCTall_fil2_myeloid, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,
        features = c("HAVCR2")) + 
        NoLegend()
VlnPlot(SCTall_fil2_myeloid, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,
        features = c("HAVCR2")) + 
        NoLegend()


VlnPlot(SCTall_fil2, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,
        features = c("bcr.ifn.1")) + 
        NoLegend()



Idents(SCTall_fil2_macs) = "minor_clusters2"

M1_mac_markers_1 = FindMarkers(SCTall_fil2_macs, ident.1 = "Macrophages - CCL18", ident.2 = c("Macrophages - CD52","Macrophages - TAMs"),assay = "SCT", logfc.threshold = 0.1)
M1_mac_markers_2 = FindMarkers(SCTall_fil2_macs, ident.1 = "Macrophages - CPB1", ident.2 = c("Macrophages - CD52","Macrophages - TAMs"),assay = "SCT", logfc.threshold = 0.1)
M1_mac_markers_3 = FindMarkers(SCTall_fil2_macs, ident.1 = "Macrophages - CXCL3", ident.2 = c("Macrophages - CD52","Macrophages - TAMs"),assay = "SCT", logfc.threshold = 0.1)
M1_mac_markers_4 = FindMarkers(SCTall_fil2_macs, ident.1 = "Macrophages - CXCL10", ident.2 = c("Macrophages - CD52","Macrophages - TAMs"),assay = "SCT", logfc.threshold = 0.1)

M1_mac_markers_1$gene = row.names(M1_mac_markers_1)
M1_mac_markers_2$gene = row.names(M1_mac_markers_2)
M1_mac_markers_3$gene = row.names(M1_mac_markers_3)
M1_mac_markers_4$gene = row.names(M1_mac_markers_4)

rbind(M1_mac_markers_1,M1_mac_markers_2,M1_mac_markers_3,M1_mac_markers_4) -> temp

temp %>% 
filter(avg_logFC > 0.25) %>% 
        group_by(gene) %>% 
        summarise(n = n(), m = mean(avg_logFC)) %>% 
        filter(n == 4) %>% 
        arrange(-m) %>% 
        pull(gene) %>% 
        cat()

ifna_ifng_df = read_csv("~/code/resources/gene_lists/ifng_ifna_genes.csv")
ifna_ifng_df$gene = toupper(ifna_ifng_df$gene)
ifna_ifng_df %>% 
        filter(gene %in% row.names(SCTall_fil2)) %>% 
        filter(stim == "IFNG") %>% 
        pull(gene) ->
        ifng_only_genes
ifna_ifng_df %>% 
        filter(gene %in% row.names(SCTall_fil2)) %>% 
        filter(stim == "IFNA") %>% 
        pull(gene) ->
        ifna_only_genes


VlnPlot(SCTall_fil2, 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 6,
        features =   ifng_only_genes) + 
        NoLegend()


VlnPlot(SCTall_fil2_epi, 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 5,
        features =   ifng_only_genes_top) + 
        NoLegend()


ifna_only_genes_top = ifna_only_genes[rowSums(SCTall_fil2_epi@assays$RNA@counts[ifna_only_genes,]) > 1000]
ifng_only_genes_top = ifng_only_genes[rowSums(SCTall_fil2_epi@assays$RNA@counts[ifng_only_genes,]) > 1000]

VlnPlot(SCTall_fil2, 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 4,
        features =   epi_myeloid_genes) + 
        NoLegend()


SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(epi_myeloid_genes), name = "rad.")

SCTall_fil2@meta.data %>% 
        #filter(timepoint=="post") %>% 
        ggplot() + aes(x=rad.1, color=timepoint) + 
        geom_density() + scale_color_manual(values = treatment_colors) +
        theme_minimal() + facet_wrap(~major_clusters, nrow = 3)


SCTall_fil2@meta.data %>% 
        #filter(patient == "BCR03") %>% 
        filter(major_clusters=="Myeloid") %>% 
        ggplot() + aes(x=rad.1, color=timepoint) + 
        geom_density() + scale_color_manual(values = treatment_colors) +
        theme_minimal() + facet_wrap(~minor_clusters2, nrow = 3)



IRDS_df = read_csv("~/code/resources/gene_lists/IRDS_genes.csv")

IRDS_df %>% 
   filter(gene %in% row.names(SCTall_fil2)) %>% 
   pull(gene) ->
   IRDS_genes
        

SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(IRDS_genes), name = "irds.")

VlnPlot(SCTall_fil2, 
        group.by = "major_clusters",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 5,
        features =   IRDS_genes[31:40]) + 
        NoLegend()


SCTall_fil2@meta.data %>% 
        #filter(timepoint=="post") %>% 
        ggplot() + aes(x=irds.1, color=timepoint) + 
        geom_density() + scale_color_manual(values = treatment_colors) +
        theme_minimal() + facet_wrap(~major_clusters, nrow = 3)

SCTall_fil2@meta.data %>% 
        filter(major_clusters=="Epithelial") %>% 
        ggplot() + aes(x=irds.1, color=timepoint) + 
        geom_density() + scale_color_manual(values = treatment_colors) +
        theme_minimal() + facet_wrap(~minor_clusters2, nrow = 3)


SCTall_fil2@meta.data %>% 
        filter(major_clusters=="Epithelial") %>% 
        ggplot() + aes(x=irds.1, color=timepoint) + 
        geom_density() + scale_color_manual(values = treatment_colors) +
        theme_minimal() + facet_wrap(~minor_clusters2, nrow = 3)

SCTall_fil2$irds.group = "mid"
SCTall_fil2$irds.group[SCTall_fil2$irds.1 >= 0.1] = "high"
SCTall_fil2$irds.group[SCTall_fil2$irds.1 <= 0] = "low"


BCR15_irds_high_cells = names(SCTall_fil2$cell.names[SCTall_fil2$patient == "BCR15" & SCTall_fil2$irds.group == "high"])
BCR15_irds_low_mid_cells = names(SCTall_fil2$cell.names[SCTall_fil2$patient == "BCR15" & SCTall_fil2$irds.group != "high"])

Idents(SCTall_fil2) = "irds.group"
BCR15_irds_markers = FindAllMarkers(SCTall_fil2[,SCTall_fil2$patient == "BCR15" & SCTall_fil2$minor_clusters2 == "Epithelial - BCR15"], only.pos = T)
BCR15_irds_pre_markers = FindAllMarkers(SCTall_fil2[,SCTall_fil2$patient == "BCR15" & SCTall_fil2$minor_clusters2 == "Epithelial - BCR15" & SCTall_fil2$timepoint == "pre"], only.pos = T, logfc.threshold = 0.1)

BCR01_irds_markers = FindAllMarkers(SCTall_fil2[,SCTall_fil2$patient == "BCR01" & SCTall_fil2$minor_clusters2 == "Epithelial - BCR01"], only.pos = T)
BCR01_irds_pre_markers = FindAllMarkers(SCTall_fil2[,SCTall_fil2$patient == "BCR01" & SCTall_fil2$minor_clusters2 == "Epithelial - BCR01" & SCTall_fil2$timepoint == "pre"], only.pos = T)


VlnPlot(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "BCR")], 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 5,
        features =   IRDS_genes[21:30]) + 
        NoLegend()


VlnPlot(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "BCR")], 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,
        features = "CD274") + 
        NoLegend()


VlnPlot(SCTall_fil2[,str_detect(SCTall_fil2$minor_clusters2, "BCR")], 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,ncol = 1,
        features = "irds.1") + 
        NoLegend()




SCTall_fil2@meta.data %>% 
        filter(major_clusters=="Epithelial") %>% 
        #filter(timepoint=="post") %>% 
        filter(str_detect(minor_clusters2, "BCR")) %>% 
        #mutate(patient = as.character(patient)) %>% 
        group_by(patient) %>% 
        mutate(m = median(irds.1)) %>% 
        ungroup() %>% 
        arrange(m) %>% 
        mutate(patient = factor(patient, levels=unique(patient))) %>%
        ggplot() + aes(x=irds.1, color=group, group=timepoint) + 
        geom_boxplot(outlier.size = 0.1) + 
        #geom_point(position = position_jitter(),size=0.1) + 
        scale_color_manual(values = group_colors) +
        theme_minimal() + 
        facet_wrap(~patient, ncol = 1, strip.position = "left") +
        theme(panel.background = element_blank(), 
              panel.grid = element_blank(),
              axis.text = element_text(size=12),
              axis.text.y = element_blank(),
              strip.text = element_text(size=12),
              axis.title = element_text(size=16)) +
        #ylab("Density") + 
        xlab("IFN-related DNA damage resistance score")

VlnPlot(SCTall_fil2_tcells, group.by="minor_clusters4", features="LAMP1")


FeaturePlot(SCTall_fil2_tcells, features=c("BCL6","CXCR3"))




VlnPlot(SCTall_fil2_macs, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors,
        features =   "CXCL9") + 
        NoLegend()


Idents(SCTall_fil2) = "minor_clusters2"
irds_avg_all = AverageExpression(SCTall_fil2, features = IRDS_genes, return.seurat = T)

Idents(SCTall_fil2) = "major_clusters"
irds_avg_maj = AverageExpression(SCTall_fil2, features = IRDS_genes, return.seurat = T)


DoHeatmap(irds_avg_all, features = irds_gene_epi_order)  + 
        scale_fill_gradient(low="white", high="firebrick") + NoLegend()
DoHeatmap(irds_avg_maj, features = IRDS_genes, slot = "data") + 
        scale_fill_gradient(low="white", high="firebrick") + NoLegend()


irds_avg_maj@assays$SCT@data %>% 
        reshape2::melt() %>% 
        rename("gene" = Var1, "cluster" = Var2, "avg" = value) %>% 
        arrange(cluster, -avg) %>% 
        filter(cluster == "Epithelial") %>% 
        pull(gene) -> irds_gene_epi_order


VlnPlot(SCTall_fil2_epi, 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=3,
        features =   HLA_genes[!HLA_genes %in% HLA_2_genes]) + 
        NoLegend()


SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","B2M"),c("HLA-A","HLA-B","HLA-C","B2M")), name = "HLA_1.")

HLA_1_genes = c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","B2M")

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="Tumor Epithelial"], 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=1,
        features = "HLA_1.1") + 
        NoLegend()


VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="Tumor Epithelial"], 
        group.by = "patient",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=1,
        features = "HLA_1.2") + 
        NoLegend()

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="CD8 Tcells"], 
        group.by = "minor_clusters2",
        pt.size = 0.1, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=2,
        features = c("PDCD1")) + 
        NoLegend()


SCTall_fil2_tcells@meta.data %>% 
        filter(str_detect(minor_clusters4, "CD4")) %>% 
        group_by(patient, PD1, timepoint) %>% 
        summarise(n = n()) %>% 
        group_by(patient) %>% 
        filter(sum(n) > 100) %>% 
        filter(patient != "BCR17") %>% 
        #filter(timepoint == "post") %>% 
        ggplot() + aes(x = timepoint, fill=PD1, y = n) + 
        geom_bar(stat="identity", position="fill") + 
        facet_wrap(~factor(patient, levels=irds_order))
        
        
        
VlnPlot(SCTall_fil2, 
        group.by = "medium_clusters",
        pt.size = 0.1, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=1,
        features = c("IL8")) + 
        NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "patient",
        pt.size = 0.1, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=2,
        features = c("HLA_1.2")) + 
        NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "patient",
        pt.size = 0.1, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=3,
        features = HLA_1_genes) + 
        NoLegend()

VlnPlot(SCTall_fil2_epi, 
        group.by = "patient",
        pt.size = 0.1, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=1,
        features = "VEGFA") + 
        NoLegend()


VlnPlot(SCTall_fil2_macs, 
        group.by = "patient",
        pt.size = 0.1, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=3,
        features = row.names(tumor_prepost_avg)[1:9]) + 
        NoLegend()



VlnPlot(SCTall_fil2_macs, 
        group.by = "minor_clusters2",
        pt.size = 0.1, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=3,
        features = "CFP") + 
        NoLegend()


SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list(row.names(tumor_prepost_avg)), name = "rad_resp.")

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "CD4 Tcells"], 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=1,
        features = "rad_resp.1") + 
        NoLegend()




HLA.II.1

SCTall_fil2@meta.data %>% 
        filter(medium_clusters == "Tumor Epithelial") %>% 
        ggplot() + aes(x= irds.1, y = HLA_1.2) +
        geom_density_2d(aes(color=patient)) +
        facet_wrap(~timepoint)

SCTall_fil2@meta.data %>% 
        filter(medium_clusters == "Tumor Epithelial") %>% 
        ggplot() + aes(x= irds.1, y = HLA_1.1,color=patient) +
        geom_point(size=0.1) +
        scale_color_manual(values = patient_colors) +
        facet_wrap(~timepoint)


plot(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"]$irds.1, 
     SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"]$HLA_1.2,
     col=as.numeric(SCTall_fil2$patient))





#SCTall_fil2_epi$irds.1 = SCTall_fil2$irds.1[colnames(SCTall_fil2_epi)]


SCTall_fil2_epi@meta.data %>% 
        filter(timepoint == "pre") %>% 
        ggplot() + aes(x = ERPRHER2, y = irds.1) + 
        geom_violin() + facet_wrap(~patient) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

SCTall_fil2_epi@meta.data %>% 
        filter(timepoint == "pre") %>% 
        left_join(selection_df) %>% 
        filter(group == "low_selection") %>% 
        mutate(ERPR = paste0("ER",ER,"PR",PR)) %>% 
        ggplot() + aes(x = ERPRHER2, y = irds.1, color=patient) + 
        geom_boxplot(outlier.color = NA) 



TIS_genes = read_csv("~/code/resources/gene_lists/TIS_genes.csv") %>% pull(gene)

Idents(SCTall_fil2) = "patient"

TIS_average = AverageExpression(SCTall_fil2, return.seurat = T, features = TIS_genes, add.ident = "timepoint")

DoHeatmap(object = TIS_average, features = TIS_genes)

plot(colSums(TIS_average@assays$SCT@data[c(1:4,8:13),]),colSums(TIS_average@assays$SCT@data[c(5:7,16:18),]))




colSums(TIS_average@assays$SCT@data)


VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "CD8 Tcells"], 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=1,
        features = "rad_resp.1") + 
        NoLegend()


VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "CD4 Tcells" & SCTall_fil2$timepoint == "pre"], 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        cols = tcell_colors, ncol=1,
        features = "rad_resp.1") + 
        NoLegend()

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "CD4 Tcells" & SCTall_fil2$timepoint == "post"], 
        group.by = "minor_clusters2",
        pt.size = 0, assay = "SCT",slot = "data", 
        cols = tcell_colors, ncol=1,
        features = "rad_resp.1") + 
        NoLegend()


VlnPlot(SCTall_fil2[,SCTall_fil2$major_clusters == "T - Cells"], 
        group.by = "minor_clusters2",
         assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=2,
        features = c("IL17A","IL17RA","IL17F","TRDC")) + 
        NoLegend()



irds_genes_top = c("STAT1", "MX1", "ISG15", "OAS1", "IFIT1", "IFIT3", "IFI44")


VlnPlot(SCTall_fil2, 
        group.by = "major_clusters",
        assay = "SCT",slot = "data", 
        split.by = "timepoint", split.plot = T,
        cols = treatment_colors, ncol=4,
        features = irds_genes_top) + 
        NoLegend()


VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial" & 
                            SCTall_fil2$timepoint == "pre"], 
        group.by = "patient",pt.size = 0,
        assay = "SCT",slot = "data", combine = F, 
        cols = patient_colors, ncol=2,
        features = irds_genes_top) + 
        NoLegend()


irds_order = c("BCR04","BCR03","BCR15","BCR05","BCR13","BCR01","BCR18","BCR20","BCR06","BCR09")

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial" & 
                            SCTall_fil2$timepoint == "pre"], 
        group.by = "patient",pt.size = 0,
        assay = "SCT",slot = "data", combine = F,
        cols = patient_colors, ncol=1,
        features = irds_genes_top) -> irds_genes_top_vln_list 

names(irds_genes_top_vln_list) = irds_genes_top

for(i in names(irds_genes_top_vln_list)){
        irds_genes_top_vln_list[[i]] = irds_genes_top_vln_list[[i]] + scale_x_discrete(limits = irds_order) + 
                theme(legend.position = "none", axis.title = element_blank())  
}

plot_grid(plotlist = irds_genes_top_vln_list, ncol = 2) -> temp
ggsave(temp, filename = "./plots/plots12/irds7_vln.pdf", height = 7, width = 6)





irds_genes_top_avg = AverageExpression(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial" & 
                                                           SCTall_fil2$timepoint == "pre"], 
                                       features = irds_genes_top, 
                                       return.seurat = T)

irds_genes_top_avg@meta.data$group = selection_df[match(colnames(irds_genes_top_avg),selection_df$patient),"group"] 
        

DoHeatmap(irds_genes_top_avg, features = irds_genes_top, group.by = "group", slot="data") + 
        scale_fill_gradient(low = "white", high = "firebrick")

Idents(irds_genes_top_avg) = factor(Idents(irds_genes_top_avg), levels = irds_order)



DoHeatmap(irds_genes_top_avg, features = irds_genes_top, slot="scale.data", group.by = "group", disp.max = 1, disp.min = -1, group.colors = group_colors)  +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", na.value = "white") -> temp
ggsave(temp, filename = "./plots/plots12/irds7_heatmap.png", height = 7, width = 6)




################
##DNAdamage_signatures
################

m_df = msigdbr(species = "Homo sapiens", category = "C2")

dna_damage_sig_names = c()
m_df %>% 
filter(str_detect(gs_name, "DNA_DAMAGE")) %>% 
        pull(gs_name) %>% 
        unique() -> temp
dna_damage_sig_names = c(dna_damage_sig_names, temp)

m_df %>% 
        filter(str_detect(gs_name, "REPAIR")) %>% 
        pull(gs_name) %>% 
        unique() -> temp
dna_damage_sig_names = c(dna_damage_sig_names, temp)

m_df %>% 
        filter(str_detect(gs_name, "HDR|HOMOLOGOUS")) %>% 
        pull(gs_name) %>% 
        unique() -> temp
dna_damage_sig_names = c(dna_damage_sig_names, temp)

m_df %>% 
        filter(str_detect(gs_name, "BREAKS|BREAK_")) %>% 
        pull(gs_name) %>% 
        unique() -> temp
dna_damage_sig_names = c(dna_damage_sig_names, temp)

m_df %>% 
        filter(str_detect(gs_name, "STRAND")) %>% 
        pull(gs_name) %>% 
        unique() -> temp
dna_damage_sig_names = c(dna_damage_sig_names, temp)

dna_damage_sig_names = unique(dna_damage_sig_names)
dna_damage_sig_names = dna_damage_sig_names[!str_detect(dna_damage_sig_names, "PROTEIN_REPAIR")]

"REACTOME_PROTEIN_REPAIR"
dna_damage_sig_scores = setNames(data.frame(matrix(ncol = length(dna_damage_sig_names), nrow = ncol(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"]))), dna_damage_sig_names)
temp1 = SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"]
for(i in dna_damage_sig_names){
  print(i)
        genes = pull(m_df[m_df$gs_name == i,"gene_symbol"])
        genes_fil = genes[genes %in% row.names(SCTall_fil2)]
        temp = AddModuleScore(temp1, features = list(genes_fil), name = i)
        dna_damage_sig_scores[i] = temp@meta.data[,paste0(i,"1")]
}
row.names(dna_damage_sig_scores) = colnames(temp1)



dna_damage_sig_scores_tum = dna_damage_sig_scores[colnames(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"]),]

set.seed(77)
cells = sample(row.names(dna_damage_sig_scores_tum))[1:1500]




ha = HeatmapAnnotation(
        patient = SCTall_fil2@meta.data[cells,"patient"], 
        timepoint = SCTall_fil2@meta.data[cells,"timepoint"],
        group = SCTall_fil2@meta.data[cells,"group"],
        col = list(patient = patient_colors,
                   timepoint = treatment_colors,
                   group = group_colors
        )
)



Heatmap(t(dna_damage_sig_scores_tum[cells,]), top_annotation = ha,show_column_names = F)





set.seed(77)
cells = sample(row.names(dna_damage_sig_scores))[1:5000]




ha = HeatmapAnnotation(
        patient = SCTall_fil2@meta.data[cells,"patient"], 
        timepoint = SCTall_fil2@meta.data[cells,"timepoint"],
        group = SCTall_fil2@meta.data[cells,"group"],
        medium_clusters = SCTall_fil2@meta.data[cells,"medium_clusters"],
        col = list(patient = patient_colors,
                   timepoint = treatment_colors,
                   group = group_colors,
                   medium_clusters = medium_cluster_colors
        )
)

Heatmap(t(dna_damage_sig_scores[cells,]), top_annotation = ha,show_column_names = F)











cbind(dna_damage_sig_scores_tum, SCTall_fil2@meta.data[row.names(dna_damage_sig_scores_tum),]) -> dna_damage_sig_scores_tum_meta
        
dna_damage_sig_boxplots = list()
for(i in dna_damage_sig_names){
        dna_damage_sig_scores_tum_meta %>% 
                ggplot() + aes_string(x = "patient", y = i, fill = "timepoint") + 
                geom_hline(yintercept = 0, color = "grey") +
                geom_boxplot()  +
                scale_fill_manual(values = treatment_colors) + 
                theme_real_minimal() -> dna_damage_sig_boxplots[[i]]
        ggsave(dna_damage_sig_boxplots[[i]], filename = paste0("./plots/plots12/dna_damage_boxplots/",i,".pdf"), height = 6, width =8)
}




dna_damage_sig_scores_tum_meta %>% 
        #dplyr::select(AMUNDSON_DNA_DAMAGE_RESPONSE_TP53:REACTOME_TELOMERE_C_STRAND_SYNTHESIS_INITIATION, patient, timepoint,cell.names) %>% 
        melt(value.name = "sig_score", variable.name = "signature") %>% 
        group_by(patient, timepoint, signature) %>% 
        summarise(m = median(sig_score)) %>% 
        spread(key = signature, value = m, fill = 0) %>% 
        ungroup() %>% 
        as.data.frame() ->
        dna_damage_sig_scores_tum_median
        
 row.names(dna_damage_sig_scores_tum_median) = paste0(dna_damage_sig_scores_tum_median$patient, "_",dna_damage_sig_scores_tum_median$timepoint)       

 
 dna_damage_sig_scores_tum_median_scale = cbind(dna_damage_sig_scores_tum_median[,c(1,2)], apply(dna_damage_sig_scores_tum_median[,c(-1,-2)], 2, scale))    
 
pdf(file = "./plots/plots12/dna_damage_score_heatmap.pdf", width = 12, height = 10) 
Heatmap(dna_damage_sig_scores_tum_median[,c(-1,-2)], row_order = c(paste0(irds_order,"_pre"), paste0(irds_order, "_post")))
dev.off()

pdf(file = "./plots/plots12/dna_damage_score_heatmap_scale.pdf", width = 12, height = 10) 
Heatmap(dna_damage_sig_scores_tum_median_scale[,c(-1,-2)], row_order = c(paste0(irds_order,"_pre"), paste0(irds_order, "_post")))
dev.off()


Heatmap(dna_damage_sig_scores_tum_median_scale[dna_damage_sig_scores_tum_median$timepoint == "pre",c(-1,-2)], row_order = c(paste0(irds_order,"_pre")))
Heatmap(dna_damage_sig_scores_tum_median_scale[dna_damage_sig_scores_tum_median$timepoint == "post",c(-1,-2)], row_order = c(paste0(irds_order,"_post")))



VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "CD8 Tcells" & 
                      SCTall_fil2$timepoint == "pre"], 
        group.by = "patient",pt.size = 0,
        assay = "SCT",slot = "data", 
        cols = patient_colors, ncol=1,
        features = "IFNG")

m_df = msigdbr(species = "Homo sapiens", category = "H")

m_df %>% 
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% 
  pull(gene_symbol) %>%  unique() -> IGR_genes

IGR_genes[!IGR_genes %in% c(IFNG.GS_genes,ISG.RS_genes)]


IFNG.GS_genes = read_csv("~/code/resources/gene_lists/Benci_IFNG.GS.csv") %>% pull(IFNG.GS)
ISG.RS_genes = read_csv("~/code/resources/gene_lists/Benci_ISG.RS.csv") %>% pull(ISG.RS)

IFNG.GS_genes[!IFNG.GS_genes %in% row.names(SCTall_fil2@assays$RNA@counts)]
ISG.RS_genes[!ISG.RS_genes %in% row.names(SCTall_fil2@assays$RNA@counts)]


SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list("IFNG.GS" = IFNG.GS_genes), name = "IFNG.GS.")
SCTall_fil2 = AddModuleScore(SCTall_fil2, features = list("ISG.RS" = ISG.RS_genes), name = "ISG.RS.")


VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial" & 
                      SCTall_fil2$timepoint == "pre"], 
        group.by = "patient",pt.size = 0,
        assay = "SCT",slot = "data", 
        cols = patient_colors, ncol=1,
        features = "ISG.RS.1")

VlnPlot(SCTall_fil2[,SCTall_fil2$major_clusters %in% 
                      c("B - Cells","T - Cells","Myeloid","NK Cells") & 
                      SCTall_fil2$timepoint == "pre"], 
        group.by = "patient",pt.size = 0,
        assay = "SCT",slot = "data", 
        cols = patient_colors, ncol=1,
        features = "IFNG.GS.1")


SCTall_fil2@meta.data %>% 
  filter(timepoint == "pre") %>% 
  mutate(tum_im = if_else(major_clusters %in% c("B - Cells","T - Cells","Myeloid","NK Cells"), "immune", "other")) %>% 
  mutate(tum_im = if_else(medium_clusters %in% c("Tumor Epithelial"), "tumor", tum_im)) %>% 
  filter(tum_im != "other") %>% 
  dplyr::select(patient, tum_im, IFNG.GS.1, ISG.RS.1) %>% 
  gather(key = "sig", value = "score", IFNG.GS.1, ISG.RS.1) %>% 
  filter((tum_im == "tumor" & sig == "ISG.RS.1") | (tum_im == "immune" & sig == "IFNG.GS.1")) %>% 
  group_by(patient, tum_im,sig) %>% 
  summarise(sig_med = median(score)) %>% 
  ungroup() %>% 
  dplyr::select(-tum_im) %>% 
  spread(key=sig, value = sig_med) %>% 
  left_join(selection_df) %>% 
  ggplot() + aes(x=IFNG.GS.1, y = ISG.RS.1, color = patient) + 
  geom_point() + 
  scale_color_manual(values = patient_colors) + 
  theme_minimal() + geom_abline(slope = 1, intercept = 0)
  
  
  
  

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"], 
        group.by = "patient",pt.size = 0,split.by = "timepoint",
        assay = "SCT",slot = "data", split.plot = T,
        cols = treatment_colors, ncol=2,
        features = c("ACSL4","SLC7A11","GPX4","PTGS2"))
  



dna_damage_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% dna_damage_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  spread(key = timepoint, value = m) %>% 
  mutate(d = post-pre) %>% 
  select(-pre, -post) %>% 
  spread(key = signature, value = d) %>% 
  column_to_rownames("patient.2") %>% 
  as.data.frame() -> dna_damage_sig_difs

col_fun = colorRamp2(c(-0.2, 0, 0.2), c(treatment_colors["pre"], "white", treatment_colors["post"]))


Heatmap(t(dna_damage_sig_difs), column_order = paste0("P0",1:8), col = col_fun)


dna_damage_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% dna_damage_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  spread(key = timepoint, value = m) %>% 
  mutate(d = post-pre) %>% 
  group_by(signature) %>% 
  mutate(m = median(d)) %>% 
  ungroup() %>% 
  top_n(48,m) %>% 
  select(-m,-d) %>% 
  gather(pre, post, key = "timepoint", value = "value") %>% 
  mutate(timepoint = factor(timepoint, levels = c("pre","post"))) %>% 
  ggplot() + aes(x = timepoint, y = value, group=patient.2) +
  facet_wrap(~signature) + geom_point() + geom_line() + 
  theme_minimal()
  


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


Heatmap(t(dna_damage_sig_difs_top), 
        column_order = paste0("P0",1:8), 
        col = col_fun, row_names_side = "left", 
        show_row_dend = F)

Heatmap(t(dna_damage_sig_difs), 
        column_order = paste0("P0",1:8), 
        col = col_fun, row_names_side = "left", 
        show_row_dend = F)

dna_damage_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% dna_damage_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  spread(key = timepoint, value = m) %>% 
  mutate(d = post-pre) %>% 
  filter(!signature %in% DNA_damage_pathways_top) %>% 
  #group_by(signature) %>% 
  #mutate(m = median(d)) %>% 
  #ungroup() %>% 
  #top_n(48,m) %>% 
  select(-pre, -post) %>% 
  spread(key = signature, value = d) %>% 
  column_to_rownames("patient.2") %>% 
  as.data.frame() -> dna_damage_sig_difs_bottom


col_fun = colorRamp2(c(-0.25, 0, 0.25), c(treatment_colors["pre"], "white", treatment_colors["post"]))


Heatmap(t(dna_damage_sig_difs_bottom), col = col_fun)




m_df %>% 
  filter(gs_name == "AMUNDSON_DNA_DAMAGE_RESPONSE_TP53") %>% 
  pull(gene_symbol) %>% 
  unique() -> temp

VlnPlot(SCTall_fil2[,SCTall_fil2$medium_clusters=="Tumor Epithelial"], 
        group.by = "patient.2",same.y.lims = T,
        pt.size = 0, assay = "SCT",slot = "data", 
        split.by = "timepoint",split.plot = T,
        cols = treatment_colors, ncol= 4,
        features = temp) + 
  NoLegend() -> AMUNDSON_vln





dna_damage_sig_scores_tum_meta %>% 
  melt(value.name = "sig_score", variable.name = "signature") %>% 
  group_by(patient.2, timepoint, signature) %>% 
  summarise(m = median(sig_score)) %>% 
  filter(signature %in% dna_damage_sig_names) %>% 
  filter(patient.2 %in% paste0("P0",1:8)) %>% 
  filter(str_detect(signature, "REACTOME_")) %>% 
  group_by(signature) %>% 
  mutate(pval = t.test(m ~ timepoint, paired=T)$p.value) %>% 
  filter(pval<= 0.05) %>% 
  ggplot() + aes(x = timepoint, y = m, group=patient.2) +
  facet_wrap(~signature) + geom_point() + geom_line() + 
  theme_minimal()

  


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
  ggplot() + aes(x = timepoint, y = m, group=patient.2) +
  facet_wrap(~reorder(signature,stat), ncol=4) + geom_point() + geom_line() + 
  theme_minimal()



#################
###cell stress###
#################



m_df = msigdbr(species = "Homo sapiens", category = "C2")

cell_stress_sig_names = c()
m_df %>% 
  filter(str_detect(gs_name, "STRESS")) %>% 
  pull(gs_name) %>% 
  unique() -> temp
cell_stress_sig_names = c(cell_stress_sig_names, temp)

m_df %>% 
  filter(str_detect(gs_name, "HYPOXIA")) %>% 
  pull(gs_name) %>% 
  unique() -> temp
cell_stress_sig_names = c(cell_stress_sig_names, temp)


cell_stress_sig_names = unique(cell_stress_sig_names)


cell_stress_sig_scores = setNames(data.frame(matrix(ncol = length(cell_stress_sig_names), nrow = ncol(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"]))), cell_stress_sig_names)
temp1 = SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"]
for(i in cell_stress_sig_names){
  print(i)
  genes = pull(m_df[m_df$gs_name == i,"gene_symbol"])
  genes_fil = genes[genes %in% row.names(SCTall_fil2)]
  temp = AddModuleScore(temp1, features = list(genes_fil), name = i)
  cell_stress_sig_scores[i] = temp@meta.data[,paste0(i,"1")]
}
row.names(cell_stress_sig_scores) = colnames(temp1)


cell_stress_sig_scores_tum = cell_stress_sig_scores[colnames(SCTall_fil2[,SCTall_fil2$medium_clusters == "Tumor Epithelial"]),]


cbind(cell_stress_sig_scores_tum, SCTall_fil2@meta.data[row.names(cell_stress_sig_scores_tum),]) -> cell_stress_sig_scores_tum_meta


