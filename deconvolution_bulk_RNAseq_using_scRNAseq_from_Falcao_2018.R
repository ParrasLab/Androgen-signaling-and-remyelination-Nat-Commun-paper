
library(data.table)
library(Seurat)
library(readr)
#download datasets from Falcao et al., Nat Med 2018 deposited at https://cells.ucsc.edu/?ds=oligo-lineage-ms

mat <- read_tsv('exprMatrix_2.tsv')
meta <- read_tsv('meta_2.tsv')

mat <- as.data.frame(mat)
rownames(mat) <- mat$...1
mat$...1 <- NULL

meta <- as.data.frame(meta)
rownames(meta) <- meta$...1
meta$...1 <- NULL

obj.seurat = CreateSeuratObject(counts = mat)
meta <- meta[row.names(obj.seurat@meta.data),]
obj.seurat@meta.data <- meta

# Get CPM avec NormalizeData
# RC = Relative counts. For CPM, sclae.factor = 1e6. No log-transformation is applied
# Inutile car counts contient les donnees en cpm
#so = NormalizeData(so,normalization.method="RC",scale.factor=1e6) 
#cpm_data = GetAssayData(so,slot = "data")

# Print final matrix 
Idents(obj.seurat) = obj.seurat@meta.data$Renamed_clusternames
list_clusters = levels(as.factor(obj.seurat@meta.data$Renamed_clusternames))
final_data = list()

for(i in 1:length(list_clusters)){
    group.cells <- WhichCells(object = obj.seurat, idents = list_clusters[i])
    data_to_write_out <- as.data.frame(GetAssayData(obj.seurat,slot = 'data')[,group.cells]) #matrice cpm gene cellule pour le cluster i
    new_row = (rep(list_clusters[i],length(group.cells))) # cell ident
    data_to_write_out = rbind(new_row,data_to_write_out)
    final_data[[i]] = data_to_write_out
}

final = do.call(cbind,final_data)
colnames(final) = final[1,]

write.table(final,file = " ExprMatrix_2.txt",col.names = F,sep = "\t",quote = F,row.names = T)

avg = AverageExpression(obj.seurat,slot = "counts") # counts = donnees brutes
write.table(avg$RNA,file = "AverageExpression_2.txt",col.names = T,sep = "\t",quote = F,row.names = T)

#####

DHTFvsCTF <- read.table(file = "DHTFvsCTF_stats_EdgeR.txt", header= TRUE, sep = "\t")
DHTFvsCTF <- DHTFvsCTF[,c("CT1F", "CT2F", "CT3F", "CT4F", "DHT1F", "DHT2F", "DHT3F")]
 
DHTFvsCTF <- DHTFvsCTF %>%  mutate(across(c("CT1F", "CT2F", "CT3F", "CT4F", "DHT1F", "DHT2F", "DHT3F"), function(x) exp(log(2)*x)))

DHTFvsCTF$Blank <- DHTFvsCTF$CT1F
DHTFvsCTF <- DHTFvsCTF[,c(8,1:7)]

write.table(DHTFvsCTF ,file = "DHTFvsCTF.txt",col.names = T, sep = "\t", quote = F, row.names = T)


DHTMvsCTM <- read.table(file = "DHTMvsCTM_stats_EdgeR.txt", header= TRUE, sep = "\t")
DHTMvsCTM <- DHTMvsCTM[,c("CT1M", "CT2M", "CT3M", "CT4M", "DHT2M", "DHT3M", "DHT4M")]

DHTMvsCTM <- DHTMvsCTM %>%  mutate(across(c("CT1M", "CT2M", "CT3M", "CT4M", "DHT2M", "DHT3M", "DHT4M"), function(x) exp(log(2)*x)))

DHTMvsCTM$Blank <- DHTMvsCTM$CT1M
DHTMvsCTM <- DHTMvsCTM[,c(8,1:7)]

write.table(DHTMvsCTM ,file = "DHTMvsCTM.txt",col.names = T, sep = "\t", quote = F, row.names = T)

#############################################################################################
#run deconvolution in CIBERSORTx  at https://cibersortx.stanford.edu/
#############################################################################################

#run the following code in CIBERSORTx (sudo docker run -v <dossier input>:/src/data -v /<dossier output>:/src/outdir cibersortx/fractions --username <login> --token <token> --single_cell TRUE --refsample ./ExprMatrix.txt --mixture ./DHTMvsCTM.txt --perm 100)


### Result

CIBERSORTx_cell_type_sourceGEP <- read.table(file = "CIBERSORTx_cell_type_sourceGEP.txt", header= TRUE, sep = "\t")

CIBERSORTx_inferred_phenoclasses_inferred_refsample <- read.table(file = "CIBERSORTx_ExprMatrix_inferred_phenoclasses.CIBERSORTx_ExprMatrix_inferred_refsample.bm.K999.txt", header= TRUE, sep = "\t")

CIBERSORTx_ExprMatrix_inferred_phenoclasses <- read.table(file = "CIBERSORTx_ExprMatrix_inferred_phenoclasses.txt", header= TRUE, sep = "\t")

CIBERSORTx_ExprMatrix_inferred_refsample <- read.table(file = "CIBERSORTx_ExprMatrix_inferred_refsample.txt", header= TRUE, sep = "\t")

CIBERSORTx_Results <- read.table(file = "CIBERSORTx_Results.txt", header= TRUE, sep = "\t")

