install.packages("dplyr")
install.packages("Seurat")
install.packages("sctransform")
install.packages("ggplot2")


library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(sctransform)

setwd("C:/.../")
##################################################################################################################################
#loading Rdata object obtained sent by David Van Bruggen from Gonçalo Castelo-Branco's lab
###################################################################################################################################
load("C:.../EAEraw.RData")

objects()
#[1] "annotationEAESS2" "erccs"            "expressions"     
dim(erccs) #92 2304
dim(expressions) #24490  2304


##################################################################################################
# Initialize the Seurat object with the raw (non-normalized data).  Keep all genes expressed in
#>= 5 cells (~0.1% of the data). Keep all cells with at least 100 detected genes
##################################################################################################

EAE_OLglia <- CreateSeuratObject(counts = expressions, project = "Falcao_2018", 
  min.cells = 10, min.features = 200)
dim(expressions) #24490  2304
dim(EAE_OLglia)      #16900  2298

# store mitochondrial percentage in object meta data    
EAE_OLglia <- PercentageFeatureSet(EAE_OLglia, pattern = "^MT-", col.name = "percent.mt")

#Apply sctransform normalization
EAE_OLglia <- SCTransform(EAE_OLglia, vars.to.regress = "percent.mt", verbose = FALSE)

#Run PCA
EAE_OLglia <- RunPCA(EAE_OLglia, verbose = FALSE)

#Visulize loadings 1 to 4
VizDimLoadings(EAE_OLglia, dims = 1:2, reduction = "pca")

#Visulize components PC1-PC2
DimPlot(EAE_OLglia, dims = c(1, 2), reduction = "pca")
ggsave("EAE_OLglia_PCs.pdf", width = 8.27, height = 6.00)

#DETERMINE THE DIMENSIONALITY OF THE DATASET

ElbowPlot(EAE_OLglia) #here, one can observe an 'elbow' around PC16-17, suggesting that the majority of true signal is captured in the first 10 PCs.
ggsave("EAE_OLglia_ElbowPlot_PCs.pdf", width = 8.27, height = 6.00)


#Find Neighbours
EAE_OLglia <- FindNeighbors(EAE_OLglia, dims = 1:30, verbose = FALSE)

# Find Clusters
EAE_OLglia <- FindClusters(EAE_OLglia, resolution = 0.8, verbose = FALSE)
#it produces 15 clusters



##########################################
#FIND NEW MARKERS FOR CLUSTERS
#########################################
EAE_OLglia.markers <- FindAllMarkers(object = EAE_OLglia, only.pos = TRUE, min.pct = 0.25, 
  thresh.use = 0.25)
dim(EAE_OLglia.markers)# 23983     7

#Top 10 markers per cluster
markers10_EAE_OLglia <- EAE_OLglia.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
head(markers10_EAE_OLglia)
print(markers10_EAE_OLglia)[c(30:40),]

#heatmap of the 10 markers of each cluster by expression levels      
library(viridis)
DoHeatmap(EAE_OLglia, features = markers10_EAE_OLglia$gene, slot = "scale.data", size = 3.5, draw.lines = TRUE) +
  scale_fill_viridis(option = "D")
ggsave("heatmap_EAE_OLglia_clusters.pdf", width = 14, height = 20)



ElbowPlot(EAE_OLglia, ndims = 25) #entre 12 et 23

#define best dimension for newly generated UMAP

dim = 19 # Dimension chosen

l = dim-5 # how many dimension before chosen dimension

h = dim+6# how many dimension after chosen dimension

p <- lapply(l:h, FUN = function(d) { # loop through dimension l to h)
  
  print(d) # print current dimension
  
  EAE_OLglia <- RunUMAP(EAE_OLglia, dims = 1:d, verbose = FALSE) # run umap using current dimension
  
  p1 <- DimPlot(EAE_OLglia, label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle(paste0("dimension ", d)) # Generate DimPlot using umap with current dim
  
  p1
  
})

CombinePlots(p)

ggsave(filename = "umap_dims_5_18_seurat.png", plot = gridExtra::arrangeGrob(grobs = p, ncol = 6),
  device = "png", width = 6000, height = 1800, units = "px", dpi = 200)


################################################# 
#Run  dimensionality reduction UMAP        
#################################################       
#Run UMAP
EAE_OLglia <- RunUMAP(EAE_OLglia, dims = 1:21, verbose = FALSE)

DimPlot(EAE_OLglia, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5) + NoLegend()
ggsave("EAE_OLglia_UMAP.pdf", width = 8.27, height = 6.00)


table(Idents(EAE_OLglia))
#  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 
# 439 219 215 210 193 174 172 128 110 101  93  73  69  60  42 

#   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14 
# 446 269 224 219 216 159 145 134  99  93  76  63  61  52  42 

##########################################
#Top 50 markers per cluster
#########################################

markers50_EAE_OLglia <- EAE_OLglia.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
head(markers50_EAE_OLglia)
dim(markers50_EAE_OLglia) # 750   7
write.table(markers50_EAE_OLglia, file='markers50_EAE_OLglia.txt', sep="\t", quote=F, col.names=T, row.names=F)

##################################################################
#Markers of different stages and cell types per cluster
#################################################################

OPCmarkers <- c( "Pdgfra", "Cspg5",  "Cspg4",  "Ascl1","Emid1","Tmem255b","Ncan", "Rlbp1",  "Kcnip3", "Rgcc", "Sdc3","Ednrb", "Pcdh15", "Sox11", "Tmem100", 
  "Spon1", "Gria3", "Traf4", "Matn4", "3110035E14Rik" )
p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  OPCmarkers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_OPC_markers_EAE_OLglia.pdf", width = 12, height = 7)

COPmarkers <-c("Tmem108", "Ptpre","Bmp4", "Pdcd4", "AI414108", "Tubb2b", "Ppfibp1", "Neu4","Trio", "Tnr",  "Arsb","Ust", "Lrrc42", "Epb4.1l2", "Sgk1", "Cyfip2","Gpr17", "Fyn")

p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  COPmarkers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_EAE_OLglia_COP_markers.pdf", width = 10, height = 7)


COPiOLmarkers <- c("Slc1a1", "Tcf7l2", "Frmd4a", "Bcan","Mpzl1", "Cadm1", "Cadm2", "Tns3", "Enpp6", "Klhl5", "Grlf1",
  "Dynll2", "Pik3r3", "Map1b", "Gpr17", "Marcks", "Fyn", "Epb4.1l2", "Cyp2j6", "Bcas1")

p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  COPiOLmarkers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_EAE_OLglia_COPiOL_markers.pdf", width = 12, height = 7)

iOLmarkers <- c("Rras2", "Rims2", "Fam107b",  "Prom1",
  "Fkbp9",  "Rap2a", "Ubl3", "Sema6a", "Tmem163", "Rell1",
  "Arpc1b",  "Add1", "Eml1", "Pdia6", "Kndc1", "Hmgcr", "Reep5", "Mt1")

p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  iOLmarkers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_EAE_OLglia_iOL_markers.pdf", width = 11, height = 7)


astromarkers <- c("Rarres2", "Ccdc153", "Gja1","Vim", "Aqp4","Mlc1", "Mia", "Enkur", "Sox9", 
  "Pltp", "Tm4sf1", "Sdc4", "Fgfr3", "Slc4a4", "Agt", "S1pr1", "Ezr", "Aldoc")
p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  astromarkers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_EAE_OLglia_astro_markers.pdf", width = 11, height = 7)

microglialmarkers <- c("C1qa", "Csf1r", "Cx3cr1", "Lyz2", "Pf4", "C1qb", "Cd14", "Tyrobp", 
  "Fcgr3",  "Laptm5", "Aif1",  "Ctss", "Ly86", "F13a1", "Ccl12", "C1qc", 
  "Ccl7", "Hpgds", "Mrc1", "Rnase4", "Ccl3", "Rgs1", "Gpr34", "Ctsc", "Apoe")

p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  microglialmarkers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_EAE_OLglia_microglial_markers.pdf", width = 14, height = 7)


MFOLmarkers= c("Tmod1","Slc9a3r2","Igsf8","Hs3st1", "Mal", "Pdlim2","Opalin", "Grb14", "Lgi3",  "Cldn11" )

p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  MFOLmarkers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_EAE_OLglia_MFOL_markers.pdf", width = 9, height = 7)


mOLmarkers= c("Apod", "Trf",  "Gpm6b", "Ugt8a",  "Fth1", "Gng11",   "Dbndd2",  "Hapln2", "Acot7", "Gstm5", "Elovl5", 
  "Spock3",  "Gjb1", "B3galt5", "Cpd", "Slco3a1", "Itgb4", "Anln", "Klk6", "Cd59a","Mgst3","Pmp22","Anxa5")

p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  mOLmarkers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_EAE_OLglia_mOL_markers.pdf", width = 12, height = 7)

mOL2markers= c("Car2", "Qdpr", "Gatm", "Tubb4a", "Stmn4", "Ndrg1", "Aplp1", "Neat1", "Il33","Dpy19l1", 
  "Serpinb1a", "Ank2", "Pls1", "Zdhhc20", "Dnajb2", "Scd1", "Trim59", "Pex5l", "Tmem229a", 
  "Rassf2", "Ppp1cc")

p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features =  mOL2markers) 
q<- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) +
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") +     labs(x=NULL, y="Cluster")
ggsave("DotPlot_EAE_OLglia_mOL2_markers.pdf", width = 12, height = 7)




VLMC_markers <- c( "Dcn", "Serping1", "Cp",  "Egr1", "Jun", "Lamc1")


FeaturePlot(EAE_OLglia, features = VLMC_markers, label = TRUE, label.size = 3,
  pt.size = 0.2, ncol = 3, cols = c("lightgrey", "red"))
ggsave("FeaturePlot_EAE_OLglia_VLMC_markers.pdf", width = 16, height = 9)

#################################################################
#Renaming the  clusters
#################################################################
new.cluster.ids <- c("MFOL" ,  "micro"  ,"mOL2" ,  "mOL6" ,  "mOL3" ,  "OPC1" , "mOL4" ,  "mOL5" ,  "OPC2"  ,  "OPC3"  , "iOL" , "mOL1" ,  "COP" ,"mOL7" , "mOL8" )
names(new.cluster.ids) <- levels(EAE_OLglia)
#rename the identifiers (cluster names)
EAE_OLglia <- RenameIdents(EAE_OLglia, new.cluster.ids)
levels(EAE_OLglia)
# "MFOL"  "micro" "mOL2"  "mOL6"  "mOL3"  "OPC1"  "mOL4"  "mOL5"  "OPC2"  "OPC3"  "iOL"   "mOL1"  "mOL7"  "COP"   "mOL8"   

DimPlot(EAE_OLglia, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5) + NoLegend()
ggsave("EAE_OLglia_UMAP.pdf", width = 8.27, height = 6.00)

table(Idents(EAE_OLglia))
#  "MFOL"  "micro" "mOL2"  "mOL6"  "mOL3"  "OPC1"  "mOL4"  "mOL5"  "OPC2"  "OPC3"  "iOL"   "mOL1"  "mOL7"  "COP"   "mOL8"   
#   439     219     215     210      193     174    172      128    110      101     93      73      69      60      42 


#################################################
#Visualizing the EAE and Control cells 
#################################################
library(stringr)
control <- colnames(EAE_OLglia)[grepl("counts_290|counts_291|counts_295" , colnames(EAE_OLglia))]  
length(control) #1146
tail(control)

EAE <- colnames(EAE_OLglia)[grepl("counts_292|counts_293|counts_294" , colnames(EAE_OLglia))]
length(EAE) #1152
head(EAE)
length(colnames(EAE_OLglia))# 2298

p <- DimPlot(EAE_OLglia, label=T, label.size = 7, pt.size = 0.6, cells.highlight= EAE, cols.highlight = "palevioletred1",  cols= "skyblue1")
p +    guides(colour = guide_legend("EAE pink, ctrl blue")) $  ggrepel
ggsave("DimPlot_EAE_ctrl_cells.pdf", width = 12, height = 8)

########################
#Reordering the clusters
########################
EAE_OLglia@active.ident <- factor(EAE_OLglia@active.ident, levels= 
    c("OPC1",  "OPC2",  "OPC3", "COP", "iOL", "MFOL", "mOL1", "mOL2",  "mOL3" , "mOL4" , "mOL5", "mOL6", "micro", "VLMC1",  "VLMC2"))
levels(EAE_OLglia)
#"OPC1"  "OPC2"  "OPC3"  "COP"   "iOL"   "MFOL"  "mOL1"  "mOL2"  "mOL3"  "mOL4"  "mOL5"  "mOL6"  "micro" "VLMC1" "VLMC2"


#################################################
#Renaming the clusters for  EAE OLglial 
#################################################
new.cluster.ids <- c("OPC1C",  "OPC2E",  "OPC3E", "COPE","iOLad", "MFOLC", "mOL1C", "mOL2E",  "mOL3E" , "mOL4C" , "mOL5E", "mOL6", "microE", "VLMC1C",  "VLMC2E")
names(new.cluster.ids) <- levels(EAE_OLglia)
#rename the identifiers (cluster names)
EAE_OLglia <- RenameIdents(EAE_OLglia, new.cluster.ids)
levels(EAE_OLglia)
#"OPC1C"  "OPC2E"  "OPC3E"  "COPE"   "iOLad"  "MFOLC"  "mOL1C"  "mOL2E"  "mOL3E"  "mOL4C"  "mOL5E"  "mOL6"   "microE" "VLMC1C" "VLMC2E" 

DimPlot(EAE_OLglia, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5) + NoLegend()
ggsave("EAE_OLglia_UMAP.pdf", width = 8.27, height = 6.00)




############################################################################################
#Dot plots for downregulated genes in female microglia
############################################################################################

downregFmicroglia <- "Cd33, Fcgr2b, Tlr4, Tnf, B2m, Ccl6, Cd9, Clec7a, Csf1, Csf2, Itgax, Anxa2, Capg, Cd52, Crip1, Ctss, Cybb, H2-K1, Ifitm3, 
Lgals1, Lgals3, Spp1, Tspo, Vim, Il1b, Rpl14, Rpl21, Rpl35, Rpsa, Tmsb4x, Ank, Fn1, Psat1"
downregFmicroglia <- strsplit(downregFmicroglia, split=", ")
downregFmicroglia <- unlist(downregFmicroglia)
downregFmicroglia <- sort(downregFmicroglia)
 downregFmicroglia
 # [1] "Ank"    "Anxa2"  "B2m"    "Capg"   "Ccl6"   "Cd33"   "Cd52"   "Cd9"    "Clec7a" "Crip1"  "Csf1"   "Csf2"   "Ctss"  
 # [14] "Cybb"   "Fcgr2b" "Fn1"    "H2-K1"  "Ifitm3" "Il1b"   "Itgax"  "Lgals1" "Lgals3" "Psat1"  "Rpl14"  "Rpl21"  "Rpl35" 
 # [27] "Rpsa"   "Spp1"   "Tlr4"   "Tmsb4x" "Tnf"    "Tspo"   "Vim"

p <- DotPlot(object = EAE_OLglia,  dot.scale = 12, features = downregFmicroglia) 
q <- p + theme_grey() 
q + theme(axis.text.x = element_text(size=12, face = "bold", angle = 30,color = "black",  hjust = 1),
  axis.text.y = element_text(size=12,face = "bold")) + 
  scale_color_viridis(alpha = 1, begin = 0.2, end = 1, direction = -1,
    discrete = FALSE, option = "D") + 
  labs(x=NULL, y=NULL)
ggsave("FeaturePlot_downregFmicroglia_genes.pdf", width = 18, height = 7)



