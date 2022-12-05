#####################################################################
#Application using gseGO using microglial genes
#####################################################################
library(edgeR)
library(ggplot2)
library(ggpubr)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("mixOmics")
library(RColorBrewer)
library(mixOmics)
library(dplyr)
library(clusterProfiler)
library(enrichplot) 
require(org.Mm.eg.db)  
require(org.Hs.eg.db)
library(AnnotationDbi)

#FEMALES
micro <- read.table("all_microglial_genes_DTHFvsCTF.txt",sep="\t", header=TRUE,  fill=TRUE, quote = "")
micro <- micro[order(-micro$logFC_DHTFvsCTF), ]
head(micro)
# 0   Lars2       0.4766317
# 56 Mir6236       0.4649637
# 58   Psat1       0.4133512
# 1      Ank       0.3307614
# 98   Timp2       0.2600000
# 96    Rxra      -0.2000000
dim(micro) #105   2
gene_list <- micro$logFC_DHTFvsCTF
names(gene_list) <- micro$ID
head(gene_list)
length(gene_list) #105
#      Lars2    Mir6236      Psat1        Ank      Timp2       Rxra 
#  0.4766317  0.4649637  0.4133512  0.3307614  0.2600000 -0.2000000 

gseMicro <- gseGO(gene= gene_list,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  keyType = "SYMBOL",
  eps = 1e-300) 
head(summary(as.data.frame(gseMicro)))
dim(as.data.frame(gseMicro))# 219  11

write.table(x= gseMicro, file='GSEA_microglial_genes_DTHFvsCTF.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 

dotplot(gseMicro, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("DotPlot_GSEA_microglial_genes_DTHFvsCTF.pdf", width = 11, height = 8)

gseaplot(gseMicro,title = gseMicro$Description[4], geneSetID = 4)
ggsave("gseaplot_4_microDHTFvsCTF.pdf", width = 9, height = 8)
gseaplot(gseMicro,title = gseMicro$Description[6], geneSetID = 6)
ggsave("gseaplot_6_microDHTFvsCTF.pdf", width = 9, height = 8)
gseaplot(gseMicro,title = gseMicro$Description[7], geneSetID = 7)
ggsave("gseaplot_7_microDHTFvsCTF.pdf", width = 9, height = 8)
gseaplot(gseMicro,title = gseMicro$Description[15], geneSetID = 15)
ggsave("gseaplot_15_microDHTFvsCTF.pdf", width = 9, height = 8)
gseaplot(gseMicro,title = gseMicro$Description[20], geneSetID = 20)
ggsave("gseaplot_20_microDHTFvsCTF.pdf", width = 9, height = 8)
gseaplot(gseMicro,title = gseMicro$Description[35], geneSetID = 35)
ggsave("gseaplot_35_microDHTFvsCTF.pdf", width = 9, height = 8)


#MALES
micro <- read.table("GSEA_microglial_genes_DTHMvsCTM.txt",sep="\t", header=TRUE,  fill=TRUE, quote = "")
micro <- micro[order(-micro$logFC_DHTMvsCTM), ]
head(micro)
#       ID logFC_DHTMvsCTM
# 14   Nos2       1.3400000
# 60   Il1b       1.2045669
# 25   Csf2       0.8100000
# 6    Ccl2       0.6500000
# 12   Igf1       0.2700000
# 105 Psat1       0.1546996
dim(micro) #105   2
gene_list <- micro$logFC_DHTMvsCTM
names(gene_list) <- micro$ID
head(gene_list)
length(gene_list) #105
#      Nos2      Il1b      Csf2      Ccl2      Igf1     Psat1 
# 1.3400000 1.2045669 0.8100000 0.6500000 0.2700000 0.1546996 

gseMicroM <- gseGO(gene= gene_list,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  keyType = "SYMBOL",
  eps = 1e-300) 
head(summary(as.data.frame(gseMicroM)))
dim(as.data.frame(gseMicroM))# 2 11

write.table(x= gseMicroM, file='GSEA_microglial_genes_DTHMvsCTM.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
dotplot(gseMicroM, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("DotPlot_GSEA_microglial_genes_DTHMvsCTM.pdf", width = 6, height = 8)






