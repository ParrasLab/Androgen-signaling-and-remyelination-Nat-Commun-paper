###########################################################
## RNAseq analysis with EdgeR from Chen & Smyth 2020
###########################################################
setwd("...")
library(edgeR)
library(ggplot2)
library(ggpubr)
#if (!requireNamespace("BiocManager", quietly = TRUE))
 #  install.packages("BiocManager")
#BiocManager::install("mixOmics")
library(RColorBrewer)
library(mixOmics)
library(viridis)
library(dplyr)
library( clusterProfiler )
library(enrichplot) 
require(org.Mm.eg.db)  
require(org.Hs.eg.db)
library(AnnotationDbi)



###########################################################
#importing the RNAseq Feature Counts 
###########################################################

data1 <- read.table("featureCounts_79_Vehicule.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(data1)
dim(data1) #55401     1

data2 <- read.table("featureCounts_80_Vehicule.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(data2)
dim(data2) #55401     1
class(data2)

data3 <- read.table("featureCounts_81_Vehicule.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(data3)
dim(data3) #55401     1
data4 <- read.table("featureCounts_82_Vehicule.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(data4)
dim(data4) #55401     1

data5 <- read.table("featureCounts_83_DHT.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(data5)
dim(data5) #55401     1

data6 <- read.table("featureCounts_84_DHT.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(data6)
dim(data6) #55401     1


data7 <- read.table("featureCounts_85_DHT.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(data7)
dim(data7) #55401     1



data <- cbind(data1, data2, data3, data4, data5, data6, data7)
head(data) 
dim(data)#55401     7
colnames(data)<- c("CT1F",   "CT2F",   "CT3F",   "CT4F", "DHT1F", "DHT2F", "DHT3F")

write.table(x=data, file='dataF.txt',
            row.names=T, col.names=T, quote=F, sep="\t", dec='.')


dat1 <- read.table("Ctrl_1.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(dat1)
dim(dat1) #55401     1

dat2 <- read.table("Ctrl_2.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(dat2)
dim(dat2) #55401     1
class(dat2)

dat3 <- read.table("Ctrl_3.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(dat3)
dim(dat3) #55401     1
dat4 <- read.table("Ctrl_8.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(dat4)
dim(dat4) #55401     1


dat5 <- read.table("DHT_15.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(dat5)
dim(dat5) #55401     1


dat6 <- read.table("DHT_16.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(dat6)
dim(dat6) #55401     1

dat7 <- read.table("DHT_18.tabular",sep="\t", header=TRUE, row.names=1, fill=TRUE, quote = "")
head(dat7)
dim(dat7) #55401     1

dat <- cbind(dat1, dat2, dat3, dat4, dat5,  dat6, dat7)
head(dat) 
dim(dat)#55401     7
colnames(dat)<- c("CT1M",   "CT2M",   "CT3M",   "CT4M", "DHT1M", "DHT2M", "DHT3M")


write.table(x=dat, file='dataM.txt',
            row.names=T, col.names=T, quote=F, sep="\t", dec='.')
write.table(x=data, file='dataF.txt',
            row.names=T, col.names=T, quote=F, sep="\t", dec='.')
head(data)
data_merged <- cbind(data, dat)
head(data_merged)
dim(data_merged) # 55401    14

data <- data_merged

#####################################################################       
#Transforming ENSEMBL IDs with .XX into SYMBOL IDs
#####################################################################   
ID <- as.character(rownames(data))
head(ID)
ID <- sub('\\.[0-9]*$', '', ID)
head(ID)
length(ID) # 55401
data <- data.frame(ID,data)
dim(data) #55401       15
head(data)

eg = bitr(ID, fromType <- "ENSEMBL", toType <- "SYMBOL", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
#Warning message:  50.81% of of input gene IDs are fail to map...
head(eg)
dim(eg) #32386     2
colnames(eg) <- c("ID", "SYMBOL")
head(eg)

dat1 = left_join(data, eg,   by = c( "ID" = "ID"))
dim(dat1) # 60533    16
head(dat1)
dat1 <- na.omit(dat1)
dim(dat1) #32386    16
head(dat1)

write.table(x=dat1, file= "rowcounts_symbols_IDs_genotypes.txt", row.names=T, col.names=T, quote=F, sep="\t", dec='.')

library(dplyr)
dim(dat1)#32386    16
dat1.nona = na.omit(dat1) #on vire le NA
dim(dat1.nona) #32386    16
dat1.nona$mean = rowMeans(dat1.nona[, 2:15]) 
head(dat1.nona)
dim(dat1.nona) #32386    17

library(data.table)
dat2 = setDT(dat1.nona)[order(-mean), head(.SD, 1) , SYMBOL] 
# Convert the 'data.frame' to 'data.table' (setDT(dat3.nona)), grouped by "SYMBOL", order the 'mean' in descending and subset the first observation with head.
head(dat2)
class(dat2)
dim(dat2) # 27369    17
dat5 <- data.frame(row.names =dat2$SYMBOL, dat2[ , c(3:16)] )
head(dat5)
dim(dat5) #27369    14

rawCountTable <- dat5 #counts file
head(rawCountTable)


######################################################################################################################  
# Normalization using house keeping genes
######################################################################################################################  
# Output of seqRUVg is a list containg both objects: 
# 1/ "W": vector of normalization factors (one factor per sample) 
# 2/ "normalizedCounts": matrix of normalized expression data 

lhg4781 <- as.character( as.matrix( read.table("List_4781HG_PMID28646208.txt", header = FALSE, sep = "\t") ) )
lhg <- intersect(rownames(rawCountTable), lhg4781)
length(lhg) # 4468

#BiocManager::install("RUVSeq")
library("RUVSeq")

rc <- as.matrix(rawCountTable)
is(rawCountTable)#dataframe
seqRUVg <- RUVg(x = rc, cIdx = lhg, k = 1)

######################################################################
# Matrix of normalized counts
######################################################################

normCount <- seqRUVg$normalizedCounts
dim(normCount) #27369    14
head(rc)
#        CT1F   CT2F   CT3F   CT4F  DHT1F  DHT2F  DHT3F   CT1M   CT2M   CT3M   CT4M  DHT1M  DHT2M  DHT3M
# COX1 712839 607003 638211 519103 668700 749155 615475 590169 498880 515203 541976 710470 594692 649654
# CYTB 353046 304389 307163 241303 375290 420031 326561 287545 243108 241255 263526 396688 354200 344112
# Mbp  264421 251241 254332 183928 407566 408552 296899 302951 285504 256075 262947 233131 292778 336583
# ND1  260348 231606 204515 161044 274689 326297 245562 171292 144248 148754 151215 308386 276980 257614
# ND5  237069 213499 210347 176211 258358 310639 222910 182931 151905 150627 161811 300482 254427 215360
# Fth1 205962 193468 235055 218210 126451 157958 174234 233203 220363 235222 217067 216844 125775 172490

#-- Plot distributions of housekeeping gene expression levels before and after normalization by RUVg --#

par(mar=c(7.1,3.1,1.1,1.1), las = 2, mfrow = c(2,1))
# HG expression before RUVg normalization
boxplot(log2(rc[lhg,]+1), col = c(rep("palegreen3",7), rep("salmon",7)), outline = FALSE)
# HG expression after RUVg normalization
boxplot(log2(normCount[lhg,]+1), col = c(rep("palegreen3",7), rep("salmon",7)), outline = FALSE)

par(las = 1, mfrow = c(1,1))
outpca <- prcomp(t(log2(normCount+1)))
plot(outpca$x, col = c(rep("palegreen3",7), rep("salmon",7)), pch = 16)
text(outpca$x, labels = rownames(outpca$x), cex = .8, pos = 4)



# Relative log expression (RLE) plots are boxplots that can be used to determine the overall quality of a data set and, 
# in particular, identify bad chips ??????
par(las = 1)
plotRLE(normCount, outline=FALSE, ylim=c(-.15, .15))


###########################################################
## PCA with ExPosition  PRINCIPAL COMPONENT ANALYSIS
###########################################################
library('ExPosition')

MyIQR=apply(normCount, 1, IQR)
##  barplot to visualize the quantiles from 0 to 1 by 10%
barplot(quantile(MyIQR, probs=seq(0,1,.01)))

##  Filter the data for > or = to .5 (50%) quantiles
MyFiltData50=normCount[MyIQR >= quantile(MyIQR, .5), ]
head(MyFiltData50)
#       CT1F   CT2F   CT3F   CT4F  DHT1F  DHT2F  DHT3F   CT1M   CT2M   CT3M   CT4M  DHT1M  DHT2M  DHT3M
# COX1 707361 600013 611199 513394 668212 736283 612878 585475 499455 516750 550305 712530 638099 656926
# CYTB 354810 306676 315869 243034 375467 424762 327455 289033 242927 240788 260941 395946 338435 341645
# Mbp  266041 253556 263192 185546 407802 414197 297895 304873 285244 255468 259791 232597 276893 333629
# ND1  263114 235311 217000 163503 274964 334139 246989 173177 144021 148144 148088 307165 251493 253715
# ND5  238030 214800 215173 177236 258457 313476 223405 183699 151813 150390 160522 300026 245197 214106
# Fth1 198806 183463 192781 207423 126028 145893 170889 224819 221529 238478 232788 219742 173725 181521


MyEpPCA = epPCA(t(MyFiltData50), scale = TRUE, center = TRUE, DESIGN = NULL,
  make_design_nominal = TRUE, graphs = F, k = 0)

mycol=rownames(MyEpPCA$ExPosition.Data$fi) #mycol is an object that is the rownames  of MyEPCA..
mycol=gsub("^CT.*$", 'red', mycol) #gsub function: substitutes everthying starting with O (^O) follows by whatever (.*$)
mycol=gsub("^DHT.*$", 'blue', mycol)
mycol <- c("blue", "blue","blue","blue", "red", "red", "red", "cadetblue3","cadetblue3","cadetblue3","cadetblue3", "coral", "coral", "coral") 

epGraphs = epGraphs(MyEpPCA, x_axis = 1, y_axis = 2, epPlotInfo = NULL, DESIGN=NULL,
  fi.col = as.matrix(mycol), fi.pch = NULL, fj.col = NULL, fj.pch = NULL, col.offset = 10,
  constraints = NULL, xlab = NULL, ylab = NULL, main = NULL,
  contributionPlots = TRUE, correlationPlotter = F, 
  graphs = TRUE)

 ######################################################################################################## 
 #Importing the "SampleInfo" and "rawCounts" files & Create a DGEList data object 
 ########################################################################################################

v1 <- colnames(normCount)
v2 <- c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "DHT", "DHT", "DHT", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "DHT", "DHT", "DHT")
v3 <- c( "F", "F", "F", "F", "F", "F", "F", "M", "M", "M", "M", "M", "M", "M" )
v4 <- c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "DHTF", "DHTF", "DHTF", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "DHTM", "DHTM", "DHTM")
v5 <- c("CTF","CTF","CTF","CTF", "DHTF", "DHTF", "DHTF", "CTM", "CTM","CTM","CTM","DHTM", "DHTM", "DHTM")
sampleInfo <- data.frame(v1,v2, v3, v4, v5)
colnames(sampleInfo) = c("sample", "condition", "gender", "merge", "allgroups")
sampleInfo

#    sample condition gender merge allgroups
# 1    CT1F      Ctrl      F  Ctrl     CTF
# 2    CT2F      Ctrl      F  Ctrl     CTF
# 3    CT3F      Ctrl      F  Ctrl     CTF
# 4    CT4F      Ctrl      F  Ctrl     CTF
# 5   DHT1F       DHT      F  DHTF    DHTF
# 6   DHT2F       DHT      F  DHTF    DHTF
# 7   DHT3F       DHT      F  DHTF    DHTF
# 8    CT1M      Ctrl      M  Ctrl     CTM
# 9    CT2M      Ctrl      M  Ctrl     CTM
# 10   CT3M      Ctrl      M  Ctrl     CTM
# 11   CT4M      Ctrl      M  Ctrl     CTM
# 12  DHT1M       DHT      M  DHTM    DHTM
# 13  DHT2M       DHT      M  DHTM    DHTM
# 14  DHT3M       DHT      M  DHTM    DHTM
 
 dgeFull <- DGEList(normCount, group=sampleInfo$merge)
 head(dgeFull)
 
 
 
 # y1 <- dgeFull
 # y1$samples$group <- 1
 # y0 <- estimateDisp(y1[lhg, ], trend="none", tagwise=F)
 # dgeFull$common.dispersion <- y0$common.dispersion
 
 
 ####################################################################################################################################################################################
 # Design matrix for DHTM vs CTM 
 ####################################################################################################################################################################################
 
 sampleInfo$allgroups <- as.factor(sampleInfo$allgroups)
 levels(sampleInfo$allgroups) #"CTF"  "CTM"  "DHTF" "DHTM"
 
 group <- as.factor(sampleInfo$allgroups)
 
 levels(group) # #"CTF"  "CTM"  "DHTF" "DHTM"
 design <- model.matrix(~0+group) 
 colnames(design) <- levels(group) 
 design
 #    CTF CTM DHTF DHTM
 # 1    1   0    0    0
 # 2    1   0    0    0
 # 3    1   0    0    0
 # 4    1   0    0    0
 # 5    0   0    1    0
 # 6    0   0    1    0
 # 7    0   0    1    0
 # 8    0   1    0    0
 # 9    0   1    0    0
 # 10   0   1    0    0
 # 11   0   1    0    0
 # 12   0   0    0    1
 # 13   0   0    0    1
 # 14   0   0    0    1
 # attr(,"assign")
 # [1] 1 1 1 1
 # attr(,"contrasts")
 # attr(,"contrasts")$group
 # [1] "contr.treatment"
 # 
 ############################################################
 # Filtering to remove low counts 
 ############################################################
 keep <- filterByExpr(dgeFull, group = group) 
 table(keep)
 #keep
 #FALSE  TRUE 
 #10666 16703 
 head(keep)
 
 
 #filtering by keep
 dgeFull <- dgeFull[keep, , keep.lib.sizes=FALSE]
 
 dim(dgeFull$counts)  # 16703    14
 
 ############################################################
 # Dispersion estimation
 ############################################################ 
 dgeFull <- estimateDisp(dgeFull, design, robust=TRUE) 
 
 fit <- glmQLFit(dgeFull, design, robust=TRUE) 
 head(fit$coefficients)
 head(fit)
 ############################################################ 
 # Differential expression analysis DHTM vs CTM
 ############################################################
 
 DHTMvsCTM  <- makeContrasts(DHTM-CTM, levels=design)
 
 res <- glmQLFTest(fit, contrast=DHTMvsCTM)
 topTags(res, n = 10)
 # Coefficient:  -1*CTM 1*DHTM 
 #            logFC   logCPM        F       PValue         FDR
 # Adgrg1 -0.4711078 6.911596 93.43650 2.856944e-07 0.002041894
 # Gpbp1   0.5393574 5.977704 84.62565 5.016173e-07 0.002041894
 # Ugcg    0.5127577 6.344091 84.13458 5.183794e-07 0.002041894
 # Cacybp  0.6724970 6.157181 83.47500 5.419223e-07 0.002041894
 # Srebf1 -0.9152669 6.310685 77.84996 8.018302e-07 0.002041894
 # Fmr1    0.7479595 5.896388 77.43255 8.262864e-07 0.002041894
 # Zfp948  0.5992285 3.651637 75.12954 9.777752e-07 0.002041894
 # Ptges3  0.6241881 6.822780 75.12676 9.779770e-07 0.002041894
 # Ube2d3  0.5465935 7.995843 69.06682 1.557306e-06 0.002858795
 # Zfp467 -0.5903106 4.767036 67.88927 1.711546e-06 0.002858795
 # 
 
 
 #The total number of DE genes identified at an FDR of 5% can be shown with decideTestsDGE 
 is.de <- decideTestsDGE(res) 
 summary(is.de)
 #        -1*CTM 1*DHTM
 # Down            1720
 # NotSig         12922
 # Up              2061
 
 
 
 #################################################################################################  
 # Creating the fusion of stat table and the logCPM for all expressed genes
 ################################################################################################# 
 
 DEGs_FDR005  <- topTags(res, n= 16703)
 head(DEGs_FDR005)
 # Coefficient:  -1*CTM 1*DHTM 
 #             logFC   logCPM        F       PValue         FDR
 # Adgrg1 -0.4711078 6.911596 93.43650 2.856944e-07 0.002041894
 # Gpbp1   0.5393574 5.977704 84.62565 5.016173e-07 0.002041894
 # Ugcg    0.5127577 6.344091 84.13458 5.183794e-07 0.002041894
 # Cacybp  0.6724970 6.157181 83.47500 5.419223e-07 0.002041894
 # Srebf1 -0.9152669 6.310685 77.84996 8.018302e-07 0.002041894
 # Fmr1    0.7479595 5.896388 77.43255 8.262864e-07 0.002041894
 
 dim(DEGs_FDR005) #16703     5
 
 DEGstat <-  data.frame(ID = rownames(DEGs_FDR005), DEGs_FDR005)
 head(DEGstat)
 #           ID      logFC   logCPM        F       PValue         FDR
 # Adgrg1 Adgrg1 -0.4711078 6.911596 93.43650 2.856944e-07 0.002041894
 # Gpbp1   Gpbp1  0.5393574 5.977704 84.62565 5.016173e-07 0.002041894
 # Ugcg     Ugcg  0.5127577 6.344091 84.13458 5.183794e-07 0.002041894
 # Cacybp Cacybp  0.6724970 6.157181 83.47500 5.419223e-07 0.002041894
 # Srebf1 Srebf1 -0.9152669 6.310685 77.84996 8.018302e-07 0.002041894
 # Fmr1     Fmr1  0.7479595 5.896388 77.43255 8.262864e-07 0.002041894
 
 DEGstat <- DEGstat[order(DEGstat$ID, decreasing = F), ]
 head(DEGstat)
 #                         ID       logFC    logCPM          F       PValue        FDR
 # 0610009B22Rik 0610009B22Rik  0.08397274 3.5245657  0.4843735 0.4988035894 0.64776212
 # 0610009E02Rik 0610009E02Rik -0.43048842 0.3013286  3.0067583 0.1067485877 0.22516050
 # 0610009L18Rik 0610009L18Rik -0.71095334 0.6483960 13.7953093 0.0026337597 0.02148032
 # 0610010F05Rik 0610010F05Rik  0.29478516 5.4776492 18.0253161 0.0009723056 0.01313536
 # 0610010K14Rik 0610010K14Rik -0.56435093 0.7906806  7.1847535 0.0189984823 0.06906021
 # 0610012G03Rik 0610012G03Rik -0.35600099 3.9017953 20.0144454 0.0006394552 0.01071368
 dim(DEGstat)
 
 
 dgeFullCPM <- cpm(dgeFull, prior.count=2, log=T) 
 head(dgeFullCPM)
 #          CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
 # COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
 # CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
 # Mbp  13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
 # ND1  13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
 # ND5  12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
 # Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
 dim(dgeFullCPM) # 16703    14
 
 DEGsCPM <-  data.frame(ID = rownames(dgeFullCPM), dgeFullCPM)
 head(DEGsCPM)
 #        ID     CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
 # COX1 COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
 # CYTB CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
 # Mbp   Mbp 13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
 # ND1   ND1 13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
 # ND5   ND5 12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
 # Fth1 Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
 dim(DEGsCPM)#  16703   15
 
 DEGsCPM <- DEGsCPM[order(DEGsCPM$ID, decreasing = F), ]
 head(DEGsCPM)
 #                         ID      CT1F      CT2F      CT3F       CT4F     DHT1F     DHT2F     DHT3F      CT1M      CT2M       CT3M      CT4M       DHT1M     DHT2M      DHT3M
 # 0610009B22Rik 0610009B22Rik 3.3267319 3.2889966 3.7131574  3.6338225 3.6424848 3.6638320 3.4858948 3.3368192 3.3731138 3.68190779 3.4511837  3.45351361 3.6395959  3.5531386
 # 0610009E02Rik 0610009E02Rik 0.1494915 0.4662663 0.1753633 -0.5622992 0.5107245 0.4909081 0.5145337 0.3761698 0.8093237 0.09139392 0.5083199  0.09431313 0.2408394 -0.1858660
 # 0610009L18Rik 0610009L18Rik 0.4200818 0.7406566 0.2941583  0.4788918 0.9736110 0.6790344 0.5145337 0.9618237 0.9853484 0.85521187 0.9800794 -0.13508851 0.4238192  0.4328421
 # 0610010F05Rik 0610010F05Rik 5.5497149 5.5320587 5.4783899  5.4256583 5.6863199 5.7140562 5.6436749 5.1764783 5.2653428 5.35964449 5.1259288  5.54069860 5.6359610  5.4006756
 # 0610010K14Rik 0610010K14Rik 0.6478472 0.8731971 0.6011559  0.9557145 0.6212767 0.5193323 0.8236557 1.4303555 0.9853484 0.69156560 1.0493971  0.13612002 0.7321061  0.6212411
 # 0610012G03Rik 0610012G03Rik 3.8912906 3.8484805 3.8093833  4.0480171 3.5788348 3.4916204 3.6157369 4.2235889 4.0473192 4.15224861 4.2402180  3.87678357 3.7187710  3.8391264
 tail(DEGstat)
 
 datDHTMvsCTM <- cbind(DEGstat, DEGsCPM[, c(2:15)] )
 head(datDHTMvsCTM)
 dim(datDHTMvsCTM) #16703    20
 datDHTMvsCTM[datDHTMvsCTM$ID =="Nipal2", ]
 # 
 #                     ID    logFC   logCPM        F      PValue        FDR     CT1F     CT2F    CT3F     CT4F DHT1F  DHT2F  DHT3F     CT1M     CT2M     CT3M  CT4M DHT1M    DHT2M    DHT3M
 #Nipal2 Nipal2 0.3199073 3.933548 12.17715 0.004039552 0.02691234 3.744199 3.785954 3.98154 3.573086 4.432566 4.34795 4.3563 3.641062 3.610329 3.724875 3.577799 3.845623 4.071208 3.948748
 #  
 
 
 
 write.table(x=datDHTMvsCTM, file='AEGs_DHTMvsCTM_logCPM_EdgeR.txt',
   row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 
 
 #################################################################################################  
 # Creating the fusion of stat table and the logCPM only for only the FDR005 genes
 ################################################################################################# 
 
 DEGsDHTMvsCTM <- subset(datDHTMvsCTM, FDR <0.05)
 dim(DEGsDHTMvsCTM) #3781   20
 
 write.table(x=DEGsDHTMvsCTM, file='DEGs_DHTMvsCTM_logCPM_EdgeR.txt',
   row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 
 ############################################################ 
 #Cluster profiler
 ############################################################
 
 # BiocManager::install("clusterProfiler")
 # BiocManager::install("org.Mm.eg.db")
 # BiocManager::install("org.Hs.eg.db")
 # BiocManager::install("pathview")
 # BiocManager::install("RDAVIDWebService")
 # browseVignettes("clusterProfiler")
 
 library( clusterProfiler )
 library(ggplot2) #Librairie pour le graph
 library(enrichplot) #Librairie pour les graphs de ClusterProfiler
 require(org.Mm.eg.db)  
 require(org.Hs.eg.db)  
 
 
 # if (!requireNamespace("BiocManager", quietly = TRUE))
 #  install.packages("BiocManager")
 # BiocManager::install(version = "3.11")
 
 # if (!requireNamespace("BiocManager", quietly = TRUE))
 #    install.packages("BiocManager")
 # BiocManager::install("AnnotationDbi")
 
 
 dat <- DEGsDHTMvsCTM
 dim(dat) #3781   20
 head(dat)
 #           CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M      logFC   logCPM        F       PValue          FDR
 #Zfp365  6.137217 6.137217 6.137217 6.137217 6.611770 6.611770 6.611770 6.137217 6.137217 6.137217 6.137217 6.611770 6.611770 6.611770  0.4745529 6.360923 87.55694 1.300372e-07 0.0008067208
 #B3galt5 7.161934 7.161934 7.161934 7.161934 7.661093 7.661093 7.661093 7.161934 7.161934 7.161934 7.161934 7.661093 7.661093 7.661093  0.4991589 7.397738 78.54227 2.595115e-07 0.0008067208
 #Inava   2.720966 2.720966 2.720966 2.720966 1.589203 1.589203 1.589203 2.720966 2.720966 2.720966 2.720966 1.589203 1.589203 1.589203 -1.1317634 2.355709 73.41449 3.964932e-07 0.0008067208
 #Sigirr  2.422466 2.422466 2.422466 2.422466 1.252664 1.252664 1.252664 2.422466 2.422466 2.422466 2.422466 1.252664 1.252664 1.252664 -1.1698014 2.051387 69.86394 5.397130e-07 0.0008067208
 #Lamp5   6.945556 6.945556 6.945556 6.945556 7.533163 7.533163 7.533163 6.945556 6.945556 6.945556 6.945556 7.533163 7.533163 7.533163  0.5876066 7.227608 68.66582 6.006816e-07 0.0008067208
 #Map9    5.808191 5.808191 5.808191 5.808191 6.465249 6.465249 6.465249 5.808191 5.808191 5.808191 5.808191 6.465249 6.465249 6.465249  0.6570581 6.128063 66.50237 7.316731e-07 0.0008067208
 

 
 dat <- dat[order(dat$PValue),]
 
 updat<- dat[dat$logFC >0, ]
 dim(updat) #2061   20
 head(updat)
 
 downdat<- dat[dat$logFC <0, ]
 dim(downdat) #1720   20
 head(downdat)
 
 # write.table(x= dat, file='DHT_EAE_merge_DEGs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 # write.table(x= updat, file='DHT_EAE_merge_upregs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 # write.table(x= downdat, file='DHT_EAE_merge_downregs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 
 ID <- as.character(dat$ID)
 eg = bitr(ID, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
 head(eg)
 dim(eg)  #3733    2
 geneList = eg$ENTREZID
 head(geneList)
 
 upID <- as.character(updat$ID)
 upeg <- bitr(upID, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
 # Warning message:     In bitr(ID, fromType <- "ENSEMBL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") :   1.7%  of input gene IDs are fail to map...
 head(upeg)
 dim(upeg) # 2033    2
 upgeneList = upeg$ENTREZID
 head(upgeneList)
 
 downID <- as.character(downdat$ID)
 downeg <- bitr(downID, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
 # Warning message:     In bitr(ID, fromType <- "ENSEMBL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") : 1.36% of input gene IDs are fail to map...
 head(downeg)
 dim(downeg) # 1700    2
 downgeneList = downeg$ENTREZID
 head(downgeneList)
 
 
 
 #############################
 #Application using enrichGO
 #############################
 
 ego <- enrichGO(gene= geneList,
   universe      = names(geneList),
   OrgDb         = org.Mm.eg.db,
   ont           = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.01,
   qvalueCutoff  = 0.05,
   readable      = TRUE) 
 head(summary(as.data.frame(ego)))
 dim(as.data.frame(ego))# 474   9
 
 
 upego <- enrichGO(gene= upgeneList,
   universe      = names(upgeneList),
   OrgDb         = org.Mm.eg.db,
   ont           = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.01,
   qvalueCutoff  = 0.05,
   readable      = TRUE) 
 head(summary(as.data.frame(upego)))
 dim(as.data.frame(upego))#  432   9
 
 downego <- enrichGO(gene= downgeneList,
   universe      = names(downgeneList),
   OrgDb         = org.Mm.eg.db,
   ont           = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.01,
   qvalueCutoff  = 0.05,
   readable      = TRUE) 
 head(summary(as.data.frame(downego)))
 dim(as.data.frame(downego))#68  9
 
 #write.table
 
 write.table(x= ego, file='GO_DEGsDHTMvsCTM_DEGs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 write.table(x= upego, file='GO_DEGsDHTMvsCTM_upregs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 write.table(x= downego, file='GO_DEGsDHTMvsCTM_downregs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 
 
 
 ############## VIZUALISATION FROM CLUSTER PROFILER #####################
 # Barplot
 # barplot(ego, showCategory=20)
 # ggsave("BarPlot_DHTFvsCT_EAE_DEGs_EdgeR_FDR005.pdf", width = 10, height = 6)
 # 
 # 
 # barplot(upego, showCategory=20)
 # ggsave("BarPlot_DHTFvsCT_EAE_upregs_EdgeR_FDR005.pdf", width = 10, height = 8)
 # 
 # barplot(downego, showCategory=20)
 # ggsave("BarPlot_DHTFvsCT_EAE_downregs_EdgeR_FDR005.pdf", width = 10, height = 8)
 
 # DotPlot count
 dotplot(ego, x = "Count", orderBy = "Count",
   color = "pvalue", font.size = 10, title = "DEGs GOby counts", showCategory=20)
 ggsave("DotPlot_DEGsDHTMvsCTM_DEGs_EdgeR_FDR005.pdf", width = 8, height = 6)
 
 dotplot(upego, x = "Count", orderBy = "Count",
   color = "pvalue", font.size = 10, title = "DEGs GOby counts", showCategory=20)
 ggsave("DotPlot_DEGsDHTMvsCTM_upreg_EdgeR_FDR005.pdf", width = 6, height = 7)
 
 dotplot(downego, x = "Count", orderBy = "Count",
   color = "pvalue", font.size = 10, title = "DEGs GOby counts", showCategory=20)
 ggsave("DotPlot_DEGsDHTMvsCTM_downreg_EdgeR_FDR005.pdf", width = 6, height = 7)
 
 
 
 # Dotplot pvalue
 # dotplot(ego, x = "pvalue", orderBy = "pvalue", 
 #   color = "pvalue", font.size = 12, title = "DEGs GO by pvalue", showCategory=40)
 # ggsave("DotPlot_GOpval_DHTFvsCT_EAE_DEGs_EdgeR_FDR005.pdf", width = 10, height = 12)
 # 
 # dotplot(upego, x = "pvalue", orderBy = "pvalue", 
 #   color = "pvalue", font.size = 12, title = "DEGs GO by pvalue", showCategory=40)
 # ggsave("DotPlot_GOpval_DHTFvsCT_EAE_upregs_EdgeR_FDR005.pdf", width = 10, height = 12)
 # 
 # dotplot(downego, x = "pvalue", orderBy = "pvalue", 
 #   color = "pvalue", font.size = 12, title = "DEGs GO by pvalue", showCategory=40)
 # ggsave("DotPlot_GOpval_DHTFvsCT_EAE_downregs_EdgeR_FDR005.pdf", width = 10, height = 12)
 
 
 # CnetPlot (depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network)
 # cnetplot(ego, categorySize="pvalue", font.size = 16, foldChange=geneList ,  circular = FALSE, colorEdge = T, showCategory=10)
 # ggsave("cnetPlot_DHTFvsCT_EAE_DEGs_EdgeR_FDR005.pdf", width = 30, height = 20)
 
 cnetplot(upego, categorySize="pvalue", font.size = 16, foldChange=geneList ,  circular = FALSE, colorEdge = T, showCategory=10)
 ggsave("cnetPlot_DEGsDHTMvsCTM_upregs_EdgeR_FDR005.pdf", width = 30, height = 20)
 
 cnetplot(downego, categorySize="pvalue", font.size = 16, foldChange=geneList ,  circular = FALSE, colorEdge = T, showCategory=10)
 ggsave("cnetPlot_DEGsDHTMvsCTM_downregs_EdgeR_FDR005.pdf", width = 30, height = 20)
 
 # 
 # # Enrichment Map (GO networks)
 # emapplot(ego, pie_scale=1.5) 
 # ggsave("emapplot2_DHTFvsCT_EAE_DEGs_EdgeR_FDR005.pdf", width = 15, height = 10)
 # emapplot(ego, showCategory = 20, color = "p.adjust")
 # ggsave("emapplot_DHTFvsCT_EAE_DEGs_EdgeR_FDR005.pdf", width = 15, height = 10)
 # 
 # 
 # emapplot(downego, pie_scale=1.5) 
 # ggsave("emapplot2_DHTFvsCT_EAE_downregs_EdgeR_FDR005.pdf", width = 15, height = 10)
 # emapplot(downego, showCategory = 20, color = "p.adjust")
 # ggsave("emapplot_DHTFvsCT_EAE_downregs_EdgeR_FDR005.pdf", width = 15, height = 10)
 # 
 
 
 
 ####################################################################################################################################################################################
 #  Differential expression analysis for DHTF vs CTF 
 ####################################################################################################################################################################################
 
 DHTFvsCTF  <- makeContrasts(DHTF-CTF, levels=design)
 
 res <- glmQLFTest(fit, contrast=DHTFvsCTF)
 topTags(res, n = 10)
 # Coefficient:  -1*CTF 1*DHTF 
 #                    logFC     logCPM        F       PValue         FDR
 # Chrm5          1.3197612  0.3061247 72.02873 1.235364e-06 0.005735475
 # 1810037I17Rik -0.4021326  5.5654145 60.68229 3.150993e-06 0.005735475
 # Nab2          -0.9061339  4.5717830 60.02690 3.340776e-06 0.005735475
 # Riox2         -0.6570302  3.7956422 59.99098 3.351553e-06 0.005735475
 # Zfp36l1       -1.1380056  6.7538767 50.84225 8.295015e-06 0.005735475
 # Swap70        -0.7439644  4.5989216 48.40727 1.040847e-05 0.005735475
 # Slc35d3        0.8296476  2.8284327 48.31812 1.050829e-05 0.005735475
 # Pnrc1         -0.6418558  5.7220823 47.57965 1.137909e-05 0.005735475
 # Creb3l4       -1.6853344 -0.1263452 45.68923 1.401592e-05 0.005735475
 # Apol8          0.7363703  2.8775567 45.46403 1.437482e-05 0.005735475
 # 
 
 
 #The total number of DE genes identified at an FDR of 5% can be shown with decideTestsDGE 
 is.de <- decideTestsDGE(res) 
 summary(is.de)
 #         -1*CTF 1*DHTF
 # Down            4185
 # NotSig          9233
 # Up              3285
 
 
 
 #################################################################################################  
 # Creating the fusion of stat table and the logCPM for all expressed genes
 ################################################################################################# 
 
 DEGs_FDR005  <- topTags(res, n= 16703)
 head(DEGs_FDR005)
 # Coefficient:  -1*CTM 1*DHTM 
 #                    logFC    logCPM        F       PValue         FDR
 # Chrm5          1.3197612 0.3061247 72.02873 1.235364e-06 0.005735475
 # 1810037I17Rik -0.4021326 5.5654145 60.68229 3.150993e-06 0.005735475
 # Nab2          -0.9061339 4.5717830 60.02690 3.340776e-06 0.005735475
 # Riox2         -0.6570302 3.7956422 59.99098 3.351553e-06 0.005735475
 # Zfp36l1       -1.1380056 6.7538767 50.84225 8.295015e-06 0.005735475
 # Swap70        -0.7439644 4.5989216 48.40727 1.040847e-05 0.005735475
 
 dim(DEGs_FDR005) #16703     5
 
 DEGstat <-  data.frame(ID = rownames(DEGs_FDR005), DEGs_FDR005)
 head(DEGstat)
 #                         ID      logFC    logCPM        F       PValue         FDR
 # Chrm5                 Chrm5  1.3197612 0.3061247 72.02873 1.235364e-06 0.005735475
 # 1810037I17Rik 1810037I17Rik -0.4021326 5.5654145 60.68229 3.150993e-06 0.005735475
 # Nab2                   Nab2 -0.9061339 4.5717830 60.02690 3.340776e-06 0.005735475
 # Riox2                 Riox2 -0.6570302 3.7956422 59.99098 3.351553e-06 0.005735475
 # Zfp36l1             Zfp36l1 -1.1380056 6.7538767 50.84225 8.295015e-06 0.005735475
 # Swap70               Swap70 -0.7439644 4.5989216 48.40727 1.040847e-05 0.005735475
 
 DEGstat <- DEGstat[order(DEGstat$ID, decreasing = F), ]
 head(DEGstat)
 #                          ID       logFC    logCPM          F       PValue        FDR
 # 0610009B22Rik 0610009B22Rik  0.09778662 3.5245657  0.6600299 0.4312965846 0.52077037
 # 0610009E02Rik 0610009E02Rik  0.41981462 0.3013286  2.9796472 0.1081866265 0.17120239
 # 0610009L18Rik 0610009L18Rik  0.25302986 0.6483960  1.8293554 0.1994518963 0.28122953
 # 0610010F05Rik 0610010F05Rik  0.18460126 5.4776492  7.0818002 0.0197065168 0.04557712
 # 0610010K14Rik 0610010K14Rik -0.12198064 0.7906806  0.3410115 0.5693233748 0.64954975
 # 0610012G03Rik 0610012G03Rik -0.34071230 3.9017953 18.2045910 0.0009351809 0.00736382
 dim(DEGstat) #16703     6
 
 
 dgeFullCPM <- cpm(dgeFull, prior.count=2, log=T) 
 head(dgeFullCPM)
 #          CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT2M    DHT2M    DHT3M
 # COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
 # CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
 # Mbp  13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
 # ND1  13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
 # ND5  12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
 # Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
 dim(dgeFullCPM) # 16703    14
 
 DEGsCPM <-  data.frame(ID = rownames(dgeFullCPM), dgeFullCPM)
 head(DEGsCPM)
 #        ID     CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT2M    DHT3M    DHT3M
 # COX1 COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
 # CYTB CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
 # Mbp   Mbp 13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
 # ND1   ND1 13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
 # ND5   ND5 12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
 # Fth1 Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
 dim(DEGsCPM)#  16703   15
 DEGsCPM <- DEGsCPM[order(DEGsCPM$ID, decreasing = F), ]
 head(DEGsCPM)
 #                         ID      CT1F      CT2F      CT3F       CT4F     DHT1F     DHT2F     DHT3F      CT1M      CT2M       CT3M      CT4M       DHT2M     DHT3M      DHT3M
 # 0610009B22Rik 0610009B22Rik 3.3267319 3.2889966 3.7131574  3.6338225 3.6424848 3.6638320 3.4858948 3.3368192 3.3731138 3.68190779 3.4511837  3.45351361 3.6395959  3.5531386
 # 0610009E02Rik 0610009E02Rik 0.1494915 0.4662663 0.1753633 -0.5622992 0.5107245 0.4909081 0.5145337 0.3761698 0.8093237 0.09139392 0.5083199  0.09431313 0.2408394 -0.1858660
 # 0610009L18Rik 0610009L18Rik 0.4200818 0.7406566 0.2941583  0.4788918 0.9736110 0.6790344 0.5145337 0.9618237 0.9853484 0.85521187 0.9800794 -0.13508851 0.4238192  0.4328421
 # 0610010F05Rik 0610010F05Rik 5.5497149 5.5320587 5.4783899  5.4256583 5.6863199 5.7140562 5.6436749 5.1764783 5.2653428 5.35964449 5.1259288  5.54069860 5.6359610  5.4006756
 # 0610010K14Rik 0610010K14Rik 0.6478472 0.8731971 0.6011559  0.9557145 0.6212767 0.5193323 0.8236557 1.4303555 0.9853484 0.69156560 1.0493971  0.13612002 0.7321061  0.6212411
 # 0610012G03Rik 0610012G03Rik 3.8912906 3.8484805 3.8093833  4.0480171 3.5788348 3.4916204 3.6157369 4.2235889 4.0473192 4.15224861 4.2402180  3.87678357 3.7187710  3.8391264
 tail(DEGstat)
 
 datDHTFvsCTF <- cbind(DEGstat, DEGsCPM[, c(2:15)] )
 head(datDHTFvsCTF)
 dim(datDHTFvsCTF) #16703    20
 datDHTFvsCTF[datDHTFvsCTF$ID =="Nipal2", ]
 
 #           ID     logFC   logCPM        F       PValue         FDR     CT1F     CT2F    CT3F     CT4F    DHT1F   DHT2F  DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
 #Nipal2 Nipal2 0.6028827 3.933548 43.96044 1.706281e-05 0.005735475 3.744199 3.785954 3.98154 3.573086 4.432566 4.34795 4.3563 3.641062 3.610329 3.724875 3.577799 3.845623 4.071208 3.948748  
 
 
 
 write.table(x=datDHTFvsCTF, file='AEGs_DHTFvsCTF_logCPM_EdgeR.txt',
   row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 
 
 #################################################################################################  
 # Creating the fusion of stat table and the logCPM only for only the FDR005 genes
 ################################################################################################# 
 
 DEGsDHTFvsCTF <- subset(datDHTFvsCTF, FDR <0.05)
 dim(DEGsDHTFvsCTF) #7470   20
 
 write.table(x=DEGsDHTFvsCTF, file='DEGsDHTFvsCTF_logCPM_EdgeR.txt',
   row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 
 
 
 
 ############################################################ 
 #Cluster profiler
 ############################################################
 
 # BiocManager::install("clusterProfiler")
 # BiocManager::install("org.Mm.eg.db")
 # BiocManager::install("org.Hs.eg.db")
 # BiocManager::install("pathview")
 # BiocManager::install("RDAVIDWebService")
 # browseVignettes("clusterProfiler")
 
 library( clusterProfiler )
 library(ggplot2) #Librairie pour le graph
 library(enrichplot) #Librairie pour les graphs de ClusterProfiler
 require(org.Mm.eg.db)  
 require(org.Hs.eg.db)  
 
 
 # if (!requireNamespace("BiocManager", quietly = TRUE))
 #  install.packages("BiocManager")
 # BiocManager::install(version = "3.11")
 
 # if (!requireNamespace("BiocManager", quietly = TRUE))
 #    install.packages("BiocManager")
 # BiocManager::install("AnnotationDbi")
 
 
 dat <- DEGsDHTFvsCTF
 dim(dat) #7470   20
 head(dat)
 #CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M      logFC   logCPM        F       PValue          FDR
 #Zfp365  6.137217 6.137217 6.137217 6.137217 6.611770 6.611770 6.611770 6.137217 6.137217 6.137217 6.137217 6.611770 6.611770 6.611770  0.4745529 6.360923 87.55694 1.300372e-07 0.0008067208
 #B3galt5 7.161934 7.161934 7.161934 7.161934 7.661093 7.661093 7.661093 7.161934 7.161934 7.161934 7.161934 7.661093 7.661093 7.661093  0.4991589 7.397738 78.54227 2.595115e-07 0.0008067208
 #Inava   2.720966 2.720966 2.720966 2.720966 1.589203 1.589203 1.589203 2.720966 2.720966 2.720966 2.720966 1.589203 1.589203 1.589203 -1.1317634 2.355709 73.41449 3.964932e-07 0.0008067208
 #Sigirr  2.422466 2.422466 2.422466 2.422466 1.252664 1.252664 1.252664 2.422466 2.422466 2.422466 2.422466 1.252664 1.252664 1.252664 -1.1698014 2.051387 69.86394 5.397130e-07 0.0008067208
 #Lamp5   6.945556 6.945556 6.945556 6.945556 7.533163 7.533163 7.533163 6.945556 6.945556 6.945556 6.945556 7.533163 7.533163 7.533163  0.5876066 7.227608 68.66582 6.006816e-07 0.0008067208
 #Map9    5.808191 5.808191 5.808191 5.808191 6.465249 6.465249 6.465249 5.808191 5.808191 5.808191 5.808191 6.465249 6.465249 6.465249  0.6570581 6.128063 66.50237 7.316731e-07 0.0008067208
 

 
 dat <- dat[order(dat$PValue),]
 
 updat<- dat[dat$logFC >0, ]
 dim(updat) #3285   20
 head(updat)
 
 downdat<- dat[dat$logFC <0, ]
 dim(downdat) #4185   20
 head(downdat)
 
 # write.table(x= dat, file='DHT_EAE_merge_DEGs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 # write.table(x= updat, file='DHT_EAE_merge_upregs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 # write.table(x= downdat, file='DHT_EAE_merge_downregs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 
 ID <- as.character(dat$ID)
 eg = bitr(ID, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
 head(eg)
 dim(eg)  #7377    2
 geneList = eg$ENTREZID
 head(geneList)
 
 upID <- as.character(updat$ID)
 upeg <- bitr(upID, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
 # Warning message:     In bitr(ID, fromType <- "ENSEMBL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") :   1.7%  of input gene IDs are fail to map...
 head(upeg)
 dim(upeg) # 3243    2
 upgeneList = upeg$ENTREZID
 head(upgeneList)
 
 downID <- as.character(downdat$ID)
 downeg <- bitr(downID, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
 # Warning message:     In bitr(ID, fromType <- "ENSEMBL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") : 1.36% of input gene IDs are fail to map...
 head(downeg)
 dim(downeg) # 4134    2
 downgeneList = downeg$ENTREZID
 head(downgeneList)
 
 
 #############################
 #Application using enrichGO
 #############################
 
 ego <- enrichGO(gene= geneList,
   universe      = names(geneList),
   OrgDb         = org.Mm.eg.db,
   ont           = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.01,
   qvalueCutoff  = 0.05,
   readable      = TRUE) 
 head(summary(as.data.frame(ego)))
 dim(as.data.frame(ego))# 2139    9
 
 
 upego <- enrichGO(gene= upgeneList,
   universe      = names(upgeneList),
   OrgDb         = org.Mm.eg.db,
   ont           = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.01,
   qvalueCutoff  = 0.05,
   readable      = TRUE) 
 head(summary(as.data.frame(upego)))
 dim(as.data.frame(upego))#  742   9
 
 downego <- enrichGO(gene= downgeneList,
   universe      = names(downgeneList),
   OrgDb         = org.Mm.eg.db,
   ont           = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.01,
   qvalueCutoff  = 0.05,
   readable      = TRUE) 
 head(summary(as.data.frame(downego)))
 dim(as.data.frame(downego))#1982    9
 
 #write.table
 
 write.table(x= ego, file='GO_DEGsDHTFvsCTF_DEGs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 write.table(x= upego, file='GO_DEGsDHTFvsCTF_upregs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 write.table(x= downego, file='GO_DEGsDHTFvsCTF_downregs_EdgeR_FDR005.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
 
 
 
 ############## VIZUALISATION FROM CLUSTER PROFILER #####################
 # Barplot
 # barplot(ego, showCategory=20)
 # ggsave("BarPlot_DHTFvsCT_EAE_DEGs_EdgeR_FDR005.pdf", width = 10, height = 6)
 # 
 # 
 # barplot(upego, showCategory=20)
 # ggsave("BarPlot_DHTFvsCT_EAE_upregs_EdgeR_FDR005.pdf", width = 10, height = 8)
 # 
 # barplot(downego, showCategory=20)
 # ggsave("BarPlot_DHTFvsCT_EAE_downregs_EdgeR_FDR005.pdf", width = 10, height = 8)
 
 # DotPlot count
 # dotplot(ego, x = "Count", orderBy = "Count",
 #   color = "pvalue", font.size = 10, title = "DEGs GOby counts", showCategory=20)
 # ggsave("DotPlot_DEGsDHTFvsCTF_DEGs_EdgeR_FDR005.pdf", width = 8, height = 6)
 
 dotplot(upego, x = "Count", orderBy = "Count",
   color = "pvalue", font.size = 10, title = "DEGs GOby counts", showCategory=20)
 ggsave("DotPlot_DEGsDHTFvsCTF_upreg_EdgeR_FDR005.pdf", width = 7, height = 7)
 
 dotplot(downego, x = "Count", orderBy = "Count",
   color = "pvalue", font.size = 10, title = "DEGs GOby counts", showCategory=20)
 ggsave("DotPlot_DEGsDHTFvsCTF_downreg_EdgeR_FDR005.pdf", width = 6, height = 7)
 
 
 
 # Dotplot pvalue
 # dotplot(ego, x = "pvalue", orderBy = "pvalue", 
 #   color = "pvalue", font.size = 12, title = "DEGs GO by pvalue", showCategory=40)
 # ggsave("DotPlot_GOpval_DHTFvsCT_EAE_DEGs_EdgeR_FDR005.pdf", width = 10, height = 12)
 # 
 # dotplot(upego, x = "pvalue", orderBy = "pvalue", 
 #   color = "pvalue", font.size = 12, title = "DEGs GO by pvalue", showCategory=40)
 # ggsave("DotPlot_GOpval_DHTFvsCT_EAE_upregs_EdgeR_FDR005.pdf", width = 10, height = 12)
 # 
 # dotplot(downego, x = "pvalue", orderBy = "pvalue", 
 #   color = "pvalue", font.size = 12, title = "DEGs GO by pvalue", showCategory=40)
 # ggsave("DotPlot_GOpval_DHTFvsCT_EAE_downregs_EdgeR_FDR005.pdf", width = 10, height = 12)
 
 
 # CnetPlot (depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network)
 # cnetplot(ego, categorySize="pvalue", font.size = 16, foldChange=geneList ,  circular = FALSE, colorEdge = T, showCategory=10)
 # ggsave("cnetPlot_DHTFvsCT_EAE_DEGs_EdgeR_FDR005.pdf", width = 30, height = 20)
 
 cnetplot(upego, categorySize="pvalue", font.size = 16, foldChange=geneList ,  circular = FALSE, colorEdge = T, showCategory=10)
 ggsave("cnetPlot_DEGsDHTFvsCTF_upregs_EdgeR_FDR005.pdf", width = 30, height = 20)
 
 cnetplot(downego, categorySize="pvalue", font.size = 16, foldChange=geneList ,  circular = FALSE, colorEdge = T, showCategory=10)
 ggsave("cnetPlot_DEGsDHTFvsCTF_downregs_EdgeR_FDR005.pdf", width = 30, height = 20)
 

####################################################################################################################################################################################
# Design matrix for DHTF vs CTs 
####################################################################################################################################################################################

   sampleInfo$merge <- as.factor(sampleInfo$merge)
   levels(sampleInfo$merge) #"Ctrl" "DHTF" "DHTM"
 
   group <- as.factor(sampleInfo$merge)
   #mouse <- as.factor(sampleInfo$mouse)
   levels(group) # "Ctrl" "DHTF" "DHTM"
   design <- model.matrix(~0+group) 
   colnames(design) <- levels(group) 
   design
   #     Ctrl DHTF DHTM
   # 1     1    0    0
   # 2     1    0    0
   # 3     1    0    0
   # 4     1    0    0
   # 5     0    1    0
   # 6     0    1    0
   # 7     0    1    0
   # 8     1    0    0
   # 9     1    0    0
   # 10    1    0    0
   # 11    1    0    0
   # 12    0    0    1
   # 13    0    0    1
   # 14    0    0    1
   # attr(,"assign")
   # [1] 1 1 1
   # attr(,"contrasts")
   # attr(,"contrasts")$group
   # [1] "contr.treatment"
   # 
   ############################################################
   # Filtering to remove low counts 
   ############################################################
   keep <- filterByExpr(dgeFull, group = group) 
   table(keep)
   #keep
   #FALSE  TRUE 
   #10666 16703 
   head(keep)
   
   
   dim(dgeFull) 
   #27369    14
   
   #filtering by keep
   dgeFull <- dgeFull[keep, , keep.lib.sizes=FALSE]
   
   dim(dgeFull$counts)  # 16703    14
   
   ############################################################
   # Dispersion estimation
   ############################################################ 
   dgeFull <- estimateDisp(dgeFull, design, robust=TRUE) 
   
   fit <- glmQLFit(dgeFull, design, robust=TRUE) 
   head(fit$coefficients)
   head(fit)
   ############################################################ 
   # Differential expression analysis DHTF vs Ctrl
   ############################################################
   
   DHTFvsCtrl  <- makeContrasts(DHTF-Ctrl, levels=design)
   
   res <- glmQLFTest(fit, contrast=DHTFvsCtrl)
   topTags(res, n = 10)
   # Coefficient:  -1*Ctrl 1*DHTF 
   #                   logFC     logCPM        F       PValue         FDR
   # Chrm5          1.2597317  0.3052936 89.37033 1.940994e-07 0.001174503
   # Nab2          -1.0120484  4.5720092 71.19785 7.580945e-07 0.001174503
   # Tdrd6          1.6076184  1.9763454 69.24889 8.925863e-07 0.001174503
   # Tenm4          0.8629905  6.8701874 61.86246 1.720604e-06 0.001174503
   # Zfp36l1       -1.1142984  6.7540460 60.56091 1.944549e-06 0.001174503
   # Nipal2         0.6710960  3.9336007 60.32448 1.988715e-06 0.001174503
   # Aqp6           1.8361465  0.8964240 60.28354 1.996478e-06 0.001174503
   # Pnrc1         -0.6533985  5.7222562 59.94721 2.061588e-06 0.001174503
   # Creb3l4       -1.7297530 -0.1253640 59.17785 2.219875e-06 0.001174503
   # 1810037I17Rik -0.3860097  5.5655784 58.30992 2.415384e-06 0.001174503
  


   #The total number of DE genes identified at an FDR of 5% can be shown with decideTestsDGE 
   is.de <- decideTestsDGE(res) 
   summary(is.de)
   #        -1*Ctrl 1*DHTF
   # Down             4618
   # NotSig           8244
   # Up               3841
   
   #1st way is to keep only the FDR005 genes

   #################################################################################################  
   # Creating the fusion of stat table and the logCPM for all expressed genes
   ################################################################################################# 
   
   DEGs_FDR005  <- topTags(res, n= 16703)
   head(DEGs_FDR005)
   #              logFC    logCPM        F       PValue        FDR
   # Chrm5    1.2597247 0.3050825 87.86976 2.087586e-07 0.00118427
   # Nab2    -1.0120388 4.5718659 71.25205 7.351703e-07 0.00118427
   # Tdrd6    1.6076335 1.9761480 69.26567 8.686926e-07 0.00118427
   # Tenm4    0.8630078 6.8700419 61.98367 1.661347e-06 0.00118427
   # Zfp36l1 -1.1142864 6.7538992 60.69636 1.875348e-06 0.00118427
   # Aqp6     1.8361654 0.8961909 60.21322 1.963674e-06 0.00118427
   
   dim(DEGs_FDR005) #16703     5
  
   
   dgeFullCPM <- cpm(dgeFull, prior.count=2, log=T) 
   head(dgeFullCPM)
   #          CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
   # COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
   # CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
   # Mbp  13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
   # ND1  13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
   # ND5  12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
   # Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
   dim(dgeFullCPM) # 16703    14
   
   DEGsCPM <-  data.frame(ID = rownames(dgeFullCPM), dgeFullCPM)
   head(DEGsCPM)
   #        ID     CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
   # COX1 COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
   # CYTB CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
   # Mbp   Mbp 13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
   # ND1   ND1 13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
   # ND5   ND5 12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
   # Fth1 Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
   dim(DEGsCPM)#  16703   15
   DEGstat <-  data.frame(ID = rownames(DEGs_FDR005), DEGs_FDR005)
   head(DEGstat)
  
   #              ID      logFC    logCPM        F       PValue        FDR
   # Chrm5     Chrm5  1.2597247 0.3050825 87.86976 2.087586e-07 0.00118427
   # Nab2       Nab2 -1.0120388 4.5718659 71.25205 7.351703e-07 0.00118427
   # Tdrd6     Tdrd6  1.6076335 1.9761480 69.26567 8.686926e-07 0.00118427
   # Tenm4     Tenm4  0.8630078 6.8700419 61.98367 1.661347e-06 0.00118427
   # Zfp36l1 Zfp36l1 -1.1142864 6.7538992 60.69636 1.875348e-06 0.00118427
   # Aqp6       Aqp6  1.8361654 0.8961909 60.21322 1.963674e-06 0.00118427
   dim(DEGstat) #16703     6
   DEGstat <- DEGstat[order(DEGstat$ID, decreasing = F), ]
   head(DEGstat)
   #                          ID       logFC    logCPM            F       PValue         FDR
   # 0610009B22Rik 0610009B22Rik  0.11539049 3.5245444 1.191452e+00 0.2934865909 0.376969127
   # 0610009E02Rik 0610009E02Rik  0.21605058 0.3011016 8.927590e-01 0.3607710083 0.446698158
   # 0610009L18Rik 0610009L18Rik -0.00174543 0.6483236 7.473577e-05 0.9932245047 0.995250399
   # 0610010F05Rik 0610010F05Rik  0.31015849 5.4776340 1.165302e+01 0.0042034192 0.012502679
   # 0610010K14Rik 0610010K14Rik -0.27756858 0.7909004 1.961683e+00 0.1831305768 0.255648142
   # 0610012G03Rik 0610012G03Rik -0.48026239 3.9018570 2.221954e+01 0.0003334297 0.002384108
   # 

   
      DEGsCPM <- DEGsCPM[order(DEGsCPM$ID, decreasing = F), ]
   head(DEGsCPM)
   dim(DEGsCPM)
   #                          ID      CT1F      CT2F      CT3F       CT4F     DHT1F     DHT2F     DHT3F      CT1M      CT2M       CT3M      CT4M       DHT1M     DHT2M      DHT3M
   # 0610009B22Rik 0610009B22Rik 3.3267319 3.2889966 3.7131574  3.6338225 3.6424848 3.6638320 3.4858948 3.3368192 3.3731138 3.68190779 3.4511837  3.45351361 3.6395959  3.5531386
   # 0610009E02Rik 0610009E02Rik 0.1494915 0.4662663 0.1753633 -0.5622992 0.5107245 0.4909081 0.5145337 0.3761698 0.8093237 0.09139392 0.5083199  0.09431313 0.2408394 -0.1858660
   # 0610009L18Rik 0610009L18Rik 0.4200818 0.7406566 0.2941583  0.4788918 0.9736110 0.6790344 0.5145337 0.9618237 0.9853484 0.85521187 0.9800794 -0.13508851 0.4238192  0.4328421
   # 0610010F05Rik 0610010F05Rik 5.5497149 5.5320587 5.4783899  5.4256583 5.6863199 5.7140562 5.6436749 5.1764783 5.2653428 5.35964449 5.1259288  5.54069860 5.6359610  5.4006756
   # 0610010K14Rik 0610010K14Rik 0.6478472 0.8731971 0.6011559  0.9557145 0.6212767 0.5193323 0.8236557 1.4303555 0.9853484 0.69156560 1.0493971  0.13612002 0.7321061  0.6212411
   # 0610012G03Rik 0610012G03Rik 3.8912906 3.8484805 3.8093833  4.0480171 3.5788348 3.4916204 3.6157369 4.2235889 4.0473192 4.15224861 4.2402180  3.87678357 3.7187710  3.8391264

   
   datDHTFvsCT <- cbind(DEGstat, DEGsCPM[, c(2:15)] )
   head(datDHTFvsCT)
   dim(datDHTFvsCT) #16703    20
  
    datDHTFvsCT[datDHTFvsCT$ID =="Nipal2", ]
   #          ID     logFC  logCPM        F       PValue        FDR     CT1F     CT2F    CT3F     CT4F    DHT1F   DHT2F  DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
   # Nipal2 Nipal2 0.6711111 3.93345 60.15782 1.974102e-06 0.00118427 3.744199 3.785954 3.98154 3.573086 4.432566 4.34795 4.3563 3.641062 3.610329 3.724875 3.577799 3.845623 4.071208 3.948748
   
  
   write.table(x=datDHTFvsCT, file='AEGs_DHTFvsCT_logCPM_EdgeR.txt',
     row.names=F, col.names=T, quote=F, sep="\t", dec='.') 

   
   #################################################################################################  
   # Creating the fusion of stat table and the logCPM only for only the FDR005 genes
   ################################################################################################# 
   
   DEGsDHTFvsCT <- subset(datDHTFvsCT, FDR <0.05)
   dim(DEGsDHTFvsCT) #8459   20
   
   write.table(x=DEGsDHTFvsCT, file='DEGs_DHTFvsCT_logCPM_EdgeR.txt',
      row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
   
    
   ############################################################
   # DHTargets
   ############################################################
 DHTargets <- toupper(c("Ar", "Bag1", "Carm1", "Crebbp", "CSF2", "CYP11A1", "Cyp11b1", "CYP11B2", "Cyp17a1", "CYP19A1", "CYP21A2", "Dnaja2", "DNAJB1", "Ep300", "ESRP1", "ESRP2", "FKBP5",
                "Hip", "Hopx", "HSD", "HSD17B2", "HSD3B1", "HSD3B2", "Hsp90aa1", "Hsp90ab1", "HSPA4", "Hspa8", "Ncoa1", "Ncoa2", "Ncoa3", "Ncor1", "Ncor2", "Pias2", "POR", "Prmt1", 
                "Shbg", "SRD53A", "SRD5A1", "SRD5A2", "ST13", "STAR", "Stip1"))


 allCPM <- cpm(res, prior.count=2, log=F) 
 dim(allCPM) # 16452     7
 head(allCPM)
 colnames(allCPM) <- colnames(res$fitted.values)
 head(allCPM) 
 allCPM <-  data.frame(ID = rownames(allCPM), allCPM)
 head(allCPM)
 
 all <- topTags(res, n =  16452)
  dim(all)
  head(all)
 
 
 allstat <-  data.frame(ID = rownames(all), all)
 head(allstat)
 dim(allstat)
 IDS=intersect(allstat$ID, allCPM$ID)
 length(IDS) #6051
 head(IDS)
 allCPM <- allCPM[match(IDS, allCPM$ID), ]
 head(allstat)
 head(allCPM)
 DEGsCPM<- as.data.frame(DEGsCPM)[, -1]
 head(DEGsCPM)
 DEGstat<- as.data.frame(DEGstat)[, -1]
 head(DEGstat)
 
 all_cpm_EdgeR <- cbind(allCPM, allstat)
 head(all_cpm_EdgeR)
 dim(all_cpm_EdgeR) #16452    14
 
DHTargetsall <- all_cpm_EdgeR[ toupper(rownames(all_cpm_EdgeR)) %in% DHTargets, ]
dim(DHTargetsall) #31  7
head(DHTargetsall)

  

write.table(x=DHTargetsall, file='DHTargets_regulation.txt',
            row.names=F, col.names=T, quote=F, sep="\t", dec='.')




 
  
 


####################################################################################################################################################################################
# Design matrix for DHTM vs CTs 
####################################################################################################################################################################################
  
  sampleInfo$merge <- as.factor(sampleInfo$merge)
  levels(sampleInfo$merge) #"Ctrl" "DHTM" "DHTM"
  
  group <- as.factor(sampleInfo$merge)
  #mouse <- as.factor(sampleInfo$mouse)
  levels(group) # "Ctrl" "DHTF" "DHTM"
  design <- model.matrix(~0+group) 
  colnames(design) <- levels(group) 
  design
  #     Ctrl DHTF DHTM
  # 1     1    0    0
  # 2     1    0    0
  # 3     1    0    0
  # 4     1    0    0
  # 5     0    1    0
  # 6     0    1    0
  # 7     0    1    0
  # 8     1    0    0
  # 9     1    0    0
  # 10    1    0    0
  # 11    1    0    0
  # 12    0    0    1
  # 13    0    0    1
  # 14    0    0    1
  # attr(,"assign")
  # [1] 1 1 1
  # attr(,"contrasts")
  # attr(,"contrasts")$group
  # [1] "contr.treatment"
  # 
  ############################################################
  # Filtering to remove low counts 
  ############################################################
  keep <- filterByExpr(dgeFull, group = group) 
  table(keep)
  #keep
  #FALSE  TRUE 
  #10666 16703 
  head(keep)
  
  
  dim(dgeFull) 
  #27369    14
  
  #filtering by keep
  dgeFull <- dgeFull[keep, , keep.lib.sizes=FALSE]
  
  dim(dgeFull$counts)  # 16703    14
  
  ############################################################
  # Dispersion estimation
  ############################################################ 
  dgeFull <- estimateDisp(dgeFull, design, robust=TRUE) 
  
  fit <- glmQLFit(dgeFull, design, robust=TRUE) 
  head(fit$coefficients)
  head(fit)
  ############################################################ 
  # Differential expression analysis DHTM vs Ctrl
  ############################################################
  
  DHTMvsCtrl  <- makeContrasts(DHTM-Ctrl, levels=design)
  
  res <- glmQLFTest(fit, contrast=DHTMvsCtrl)
  topTags(res, n = 10)
  # Coefficient:  -1*Ctrl 1*DHTM 
  #             logFC     logCPM        F       PValue       FDR
  # Zfp365   0.3911191  6.3618217 23.42200 0.0002631443 0.2893611
  # Ttn      2.1021740  0.7225588 19.41292 0.0006345130 0.2893611
  # Ckm      2.9983089 -0.7867570 18.86178 0.0007220057 0.2893611
  # Tnnt3    2.3948821  0.3689932 18.48074 0.0008010942 0.2893611
  # Myh4     4.4438378  0.3437379 18.17169 0.0009861639 0.2893611
  # Sigirr  -0.7888094  2.0467742 16.57624 0.0011468558 0.2893611
  # Gstm7    0.5233526  4.7121049 16.02306 0.0013113945 0.2893611
  # B3galt5  0.3313440  7.3982941 15.83843 0.0013722265 0.2893611
  # Atp2a1   3.0131618  0.6519693 16.21172 0.0014834622 0.2893611
  # Fgd3    -0.3034843  5.4537744 15.39091 0.0015335945 0.2893611
  
  
  
  #The total number of DE genes identified at an FDR of 5% can be shown with decideTestsDGE 
  is.de <- decideTestsDGE(res) 
  summary(is.de)
  #        -1*Ctrl 1*DHTM
  # Down                0
  # NotSig          16703
  # Up                  0
  

 
  #################################################################################################  
  # Creating the fusion of stat table and the logCPM for all expressed genes
  ################################################################################################# 
  
  DEGs_FDR005  <- topTags(res, n= 16703)
  head(DEGs_FDR005)
  #             logFC     logCPM        F       PValue       FDR
  # Zfp365  0.3911191  6.3618217 23.42200 0.0002631443 0.2893611
  # Ttn     2.1021740  0.7225588 19.41292 0.0006345130 0.2893611
  # Ckm     2.9983089 -0.7867570 18.86178 0.0007220057 0.2893611
  # Tnnt3   2.3948821  0.3689932 18.48074 0.0008010942 0.2893611
  # Myh4    4.4438378  0.3437379 18.17169 0.0009861639 0.2893611
  # Sigirr -0.7888094  2.0467742 16.57624 0.0011468558 0.2893611
  
  dim(DEGs_FDR005) #16703     5
  
  DEGstat <-  data.frame(ID = rownames(DEGs_FDR005), DEGs_FDR005)
  head(DEGstat)
  #            ID      logFC     logCPM        F       PValue       FDR
  # Zfp365 Zfp365  0.3911191  6.3618217 23.42200 0.0002631443 0.2893611
  # Ttn       Ttn  2.1021740  0.7225588 19.41292 0.0006345130 0.2893611
  # Ckm       Ckm  2.9983089 -0.7867570 18.86178 0.0007220057 0.2893611
  # Tnnt3   Tnnt3  2.3948821  0.3689932 18.48074 0.0008010942 0.2893611
  # Myh4     Myh4  4.4438378  0.3437379 18.17169 0.0009861639 0.2893611
  # Sigirr Sigirr -0.7888094  2.0467742 16.57624 0.0011468558 0.2893611
  
  DEGstat <- DEGstat[order(DEGstat$ID, decreasing = F), ]
  head(DEGstat)
  #                          ID       logFC    logCPM         F     PValue       FDR
  # 0610009B22Rik 0610009B22Rik  0.06612572 3.5245444 0.3877315 0.54351552 0.7502140
  # 0610009E02Rik 0610009E02Rik -0.25204772 0.3011016 1.1065969 0.31066178 0.5705315
  # 0610009L18Rik 0610009L18Rik -0.49500519 0.6483236 5.4790064 0.03458738 0.2895048
  # 0610010F05Rik 0610010F05Rik  0.15729201 5.4776340 2.9459498 0.10815132 0.3639896
  # 0610010K14Rik 0610010K14Rik -0.42411728 0.7909004 4.4219596 0.05406955 0.3019614
  # 0610012G03Rik 0610012G03Rik -0.22884624 3.9018570 5.1938639 0.03888489 0.2895048
  dim(DEGstat)
  
  
  dgeFullCPM <- cpm(dgeFull, prior.count=2, log=T) 
  head(dgeFullCPM)
  #          CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  # COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
  # CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
  # Mbp  13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
  # ND1  13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
  # ND5  12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
  # Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
  dim(dgeFullCPM) # 16703    14
  
  DEGsCPM <-  data.frame(ID = rownames(dgeFullCPM), dgeFullCPM)
  head(DEGsCPM)
  #        ID     CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  # COX1 COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
  # CYTB CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
  # Mbp   Mbp 13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
  # ND1   ND1 13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
  # ND5   ND5 12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
  # Fth1 Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
  dim(DEGsCPM)#  16703   15
  DEGsCPM <- DEGsCPM[order(DEGsCPM$ID, decreasing = F), ]
  head(DEGsCPM)
  #                         ID      CT1F      CT2F      CT3F       CT4F     DHT1F     DHT2F     DHT3F      CT1M      CT2M       CT3M      CT4M       DHT1M     DHT2M      DHT3M
  # 0610009B22Rik 0610009B22Rik 3.3267319 3.2889966 3.7131574  3.6338225 3.6424848 3.6638320 3.4858948 3.3368192 3.3731138 3.68190779 3.4511837  3.45351361 3.6395959  3.5531386
  # 0610009E02Rik 0610009E02Rik 0.1494915 0.4662663 0.1753633 -0.5622992 0.5107245 0.4909081 0.5145337 0.3761698 0.8093237 0.09139392 0.5083199  0.09431313 0.2408394 -0.1858660
  # 0610009L18Rik 0610009L18Rik 0.4200818 0.7406566 0.2941583  0.4788918 0.9736110 0.6790344 0.5145337 0.9618237 0.9853484 0.85521187 0.9800794 -0.13508851 0.4238192  0.4328421
  # 0610010F05Rik 0610010F05Rik 5.5497149 5.5320587 5.4783899  5.4256583 5.6863199 5.7140562 5.6436749 5.1764783 5.2653428 5.35964449 5.1259288  5.54069860 5.6359610  5.4006756
  # 0610010K14Rik 0610010K14Rik 0.6478472 0.8731971 0.6011559  0.9557145 0.6212767 0.5193323 0.8236557 1.4303555 0.9853484 0.69156560 1.0493971  0.13612002 0.7321061  0.6212411
  # 0610012G03Rik 0610012G03Rik 3.8912906 3.8484805 3.8093833  4.0480171 3.5788348 3.4916204 3.6157369 4.2235889 4.0473192 4.15224861 4.2402180  3.87678357 3.7187710  3.8391264
  tail(DEGstat)
  
  datDHTMvsCT <- cbind(DEGstat, DEGsCPM[, c(2:15)] )
  head(datDHTMvsCT)
  dim(datDHTMvsCT) #16703    20
  dat[dat$ID =="Nipal2", ]
  # 
  #            ID    logFC  logCPM          F    PValue       FDR         CT1F     CT2F    CT3F     CT4F    DHT1F   DHT2F  DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  # Nipal2 Nipal2 0.248365 3.93345  7.862668  0.0140802    0.2893611 3.744199  3.785954 3.98154 3.573086 4.432566 4.34795 4.3563 3.641062 3.610329 3.724875 3.577799 3.845623 4.071208 3.948748
  #  

  
  
  write.table(x=datDHTMvsCT, file='AEGs_DHTMvsCT_logCPM_EdgeR.txt',
     row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
  
  
  #################################################################################################  
  # Creating the fusion of stat table and the logCPM only for only the FDR005 genes
  ################################################################################################# 
  
  DEGsDHTMvsCT <- subset(datDHTMvsCT, FDR <0.05)
  dim(DEGsDHTMvsCT) #0   20 > no DEGs in this comparison
  
  # write.table(x=DEGsDHTFvsCT, file='DEGs_DHTFvsCT_logCPM_EdgeR.txt',
  #    row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
  
  
  
  
  
 

  
  
####################################################################################################################################################################################
#  Differential expression analysis for DHTF vs DHTM
####################################################################################################################################################################################

  
  DHTFvsDHTM  <- makeContrasts(DHTF-DHTM, levels=design)
  
  res <- glmQLFTest(fit, contrast=DHTFvsDHTM)
  topTags(res, n = 10)
  # Coefficient:  1*DHTF -1*DHTM 
  #                    logFC     logCPM        F       PValue       FDR
  # Chrm5          1.0478837  0.3061247 39.73389 2.835611e-05 0.1114594
  # Sf3b5         -0.5229268  3.7717099 31.02361 9.358475e-05 0.1114594
  # Tmem258       -0.5435836  4.3747399 29.64570 1.155812e-04 0.1114594
  # Sf3b4         -0.3744896  5.0493043 29.59917 1.164224e-04 0.1114594
  # Gm41611       -2.6498175 -0.5660913 29.27201 1.225413e-04 0.1114594
  # Riox2         -0.4755971  3.7956422 28.08725 1.480370e-04 0.1114594
  # Ccdc130       -0.4550647  3.3732230 27.42906 1.648341e-04 0.1114594
  # Nfkbib        -0.4552927  4.4129297 25.62956 2.232697e-04 0.1114594
  # Trir          -0.3816206  5.7820426 24.21451 2.864544e-04 0.1114594
  # 1810037I17Rik -0.2689972  5.5654145 24.08963 2.929638e-04 0.1114594
  
  
  #The total number of DE genes identified at an FDR of 5% can be shown with decideTestsDGE 
  is.de <- decideTestsDGE(res) 
  summary(is.de)
  #        1*DHTF -1*DHTM
  # Down                0
  # NotSig          16703
  # Up                  0
  
  
  
  #################################################################################################  
  # Creating the fusion of stat table and the logCPM for all expressed genes
  ################################################################################################# 
  
  DEGs_FDR005  <- topTags(res, n= 16703)
  head(DEGs_FDR005)
  #              logFC     logCPM        F       PValue       FDR
  # Chrm5    1.0478837  0.3061247 39.73389 2.835611e-05 0.1114594
  # Sf3b5   -0.5229268  3.7717099 31.02361 9.358475e-05 0.1114594
  # Tmem258 -0.5435836  4.3747399 29.64570 1.155812e-04 0.1114594
  # Sf3b4   -0.3744896  5.0493043 29.59917 1.164224e-04 0.1114594
  # Gm41611 -2.6498175 -0.5660913 29.27201 1.225413e-04 0.1114594
  # Riox2   -0.4755971  3.7956422 28.08725 1.480370e-04 0.1114594
  
  dim(DEGs_FDR005) #16703     5
  
  DEGstat <-  data.frame(ID = rownames(DEGs_FDR005), DEGs_FDR005)
  head(DEGstat)
  #              ID      logFC     logCPM        F       PValue       FDR
  # Chrm5     Chrm5  1.0478837  0.3061247 39.73389 2.835611e-05 0.1114594
  # Sf3b5     Sf3b5 -0.5229268  3.7717099 31.02361 9.358475e-05 0.1114594
  # Tmem258 Tmem258 -0.5435836  4.3747399 29.64570 1.155812e-04 0.1114594
  # Sf3b4     Sf3b4 -0.3744896  5.0493043 29.59917 1.164224e-04 0.1114594
  # Gm41611 Gm41611 -2.6498175 -0.5660913 29.27201 1.225413e-04 0.1114594
  # Riox2     Riox2 -0.4755971  3.7956422 28.08725 1.480370e-04 0.1114594
  
  DEGstat <- DEGstat[order(DEGstat$ID, decreasing = F), ]
  head(DEGstat)
  #                          ID      logFC    logCPM         F     PValue       FDR
  # 0610009B22Rik 0610009B22Rik  0.0493005 3.5245657 0.1466301 0.70801610 0.8044893
  # 0610009E02Rik 0610009E02Rik  0.4681297 0.3013286 3.1958126 0.09734428 0.2217596
  # 0610009L18Rik 0610009L18Rik  0.4934448 0.6483960 5.9216946 0.03027511 0.1325181
  # 0610010F05Rik 0610010F05Rik  0.1528711 5.4776492 4.2269419 0.06062225 0.1725730
  # 0610010K14Rik 0610010K14Rik  0.1464326 0.7906806 0.4235603 0.52659605 0.6634530
  # 0610012G03Rik 0610012G03Rik -0.2514394 3.9017953 8.7670804 0.01112335 0.1114594
  dim(DEGstat) #16703     6
  
  
  dgeFullCPM <- cpm(dgeFull, prior.count=2, log=T) 
  head(dgeFullCPM)
  #          CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  # COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
  # CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
  # Mbp  13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
  # ND1  13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
  # ND5  12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
  # Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
  dim(dgeFullCPM) # 16703    14
  
  DEGsCPM <-  data.frame(ID = rownames(dgeFullCPM), dgeFullCPM)
  head(DEGsCPM)
  #        ID     CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  # COX1 COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
  # CYTB CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
  # Mbp   Mbp 13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
  # ND1   ND1 13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
  # ND5   ND5 12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
  # Fth1 Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
  dim(DEGsCPM)#  16703   15
  DEGsCPM <- DEGsCPM[order(DEGsCPM$ID, decreasing = F), ]
  head(DEGsCPM)
  #                         ID      CT1F      CT2F      CT3F       CT4F     DHT1F     DHT2F     DHT3F      CT1M      CT2M       CT3M      CT4M       DHT1M     DHT2M      DHT3M
  # 0610009B22Rik 0610009B22Rik 3.3267319 3.2889966 3.7131574  3.6338225 3.6424848 3.6638320 3.4858948 3.3368192 3.3731138 3.68190779 3.4511837  3.45351361 3.6395959  3.5531386
  # 0610009E02Rik 0610009E02Rik 0.1494915 0.4662663 0.1753633 -0.5622992 0.5107245 0.4909081 0.5145337 0.3761698 0.8093237 0.09139392 0.5083199  0.09431313 0.2408394 -0.1858660
  # 0610009L18Rik 0610009L18Rik 0.4200818 0.7406566 0.2941583  0.4788918 0.9736110 0.6790344 0.5145337 0.9618237 0.9853484 0.85521187 0.9800794 -0.13508851 0.4238192  0.4328421
  # 0610010F05Rik 0610010F05Rik 5.5497149 5.5320587 5.4783899  5.4256583 5.6863199 5.7140562 5.6436749 5.1764783 5.2653428 5.35964449 5.1259288  5.54069860 5.6359610  5.4006756
  # 0610010K14Rik 0610010K14Rik 0.6478472 0.8731971 0.6011559  0.9557145 0.6212767 0.5193323 0.8236557 1.4303555 0.9853484 0.69156560 1.0493971  0.13612002 0.7321061  0.6212411
  # 0610012G03Rik 0610012G03Rik 3.8912906 3.8484805 3.8093833  4.0480171 3.5788348 3.4916204 3.6157369 4.2235889 4.0473192 4.15224861 4.2402180  3.87678357 3.7187710  3.8391264
  tail(DEGstat)
  
  datDHTFvsDHTM <- cbind(DEGstat, DEGsCPM[, c(2:15)] )
  head(datDHTFvsDHTM)
  dim(datDHTFvsDHTM) #16703    20
  datDHTFvsDHTM[datDHTFvsDHTM$ID =="Nipal2", ]
  
  #           ID     logFC   logCPM        F       PValue       FDR     CT1F     CT2F    CT3F     CT4F    DHT1F   DHT2F  DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  #Nipal2 Nipal2 0.4227549 3.933548 18.68717 0.0008431334 0.1114594 3.744199 3.785954 3.98154 3.573086 4.432566 4.34795 4.3563 3.641062 3.610329 3.724875 3.577799 3.845623 4.071208 3.948748
  
  
  write.table(x=datDHTFvsDHTM, file='AEGs_DHTFvsDHTM_logCPM_EdgeR.txt',
     row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
  
  
  #################################################################################################  
  # Creating the fusion of stat table and the logCPM only for only the FDR005 genes
  ################################################################################################# 
  
  DEGsDHTFvsDHTM <- subset(datDHTFvsDHTM, FDR <0.05)
  dim(DEGsDHTFvsDHTM) #0   20
  
  # write.table(x=DEGsDHTFvsCTF, file='DEGsDHTFvsCTF_logCPM_EdgeR.txt',
  #    row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
  
  
  
  
  
  ####################################################################################################################################################################################
  #  Differential expression analysis for CTF vs CTM
  ####################################################################################################################################################################################
  
  
  CTFvsCTM  <- makeContrasts(CTF-CTM, levels=design)
  
  res <- glmQLFTest(fit, contrast=CTFvsCTM)
  topTags(res, n = 10)
  # Coefficient:  1*CTF -1*CTM 
  #              logFC   logCPM        F       PValue         FDR
  # Adgrg1  -0.4410666 6.911596 97.05735 2.297839e-07 0.002192422
  # Gpbp1    0.5049153 5.977704 85.01982 4.886149e-07 0.002192422
  # Ugcg     0.4676207 6.344091 80.30188 6.739704e-07 0.002192422
  # Metrn   -0.7118352 5.860048 78.55347 7.624780e-07 0.002192422
  # Zfp948   0.5688350 3.651637 77.39104 8.287655e-07 0.002192422
  # Fmr1     0.6944720 5.896388 76.06021 9.129869e-07 0.002192422
  # Ube2d3   0.5354158 7.995843 75.94336 9.208428e-07 0.002192422
  # Ptges3   0.5779628 6.822780 73.67297 1.090146e-06 0.002192422
  # Slc24a4 -0.8491843 4.444590 72.61309 1.181332e-06 0.002192422
  # Rabac1  -0.5813836 5.822430 67.95529 1.702445e-06 0.002843594
  
  
  #The total number of DE genes identified at an FDR of 5% can be shown with decideTestsDGE 
  is.de <- decideTestsDGE(res) 
  summary(is.de)
  #        1*CTF -1*CTM
  # Down            869
  # NotSig        14647
  # Up             1187
  
  
  
  #################################################################################################  
  # Creating the fusion of stat table and the logCPM for all expressed genes
  ################################################################################################# 
  
  DEGs_FDR005  <- topTags(res, n= 16703)
  head(DEGs_FDR005)
  #           logFC   logCPM        F       PValue         FDR
  # Adgrg1 -0.4410666 6.911596 97.05735 2.297839e-07 0.002192422
  # Gpbp1   0.5049153 5.977704 85.01982 4.886149e-07 0.002192422
  # Ugcg    0.4676207 6.344091 80.30188 6.739704e-07 0.002192422
  # Metrn  -0.7118352 5.860048 78.55347 7.624780e-07 0.002192422
  # Zfp948  0.5688350 3.651637 77.39104 8.287655e-07 0.002192422
  # Fmr1    0.6944720 5.896388 76.06021 9.129869e-07 0.002192422
  
  dim(DEGs_FDR005) #16703     5
  
  DEGstat <-  data.frame(ID = rownames(DEGs_FDR005), DEGs_FDR005)
  head(DEGstat)
  #           ID      logFC   logCPM        F       PValue         FDR
  # Adgrg1 Adgrg1 -0.4410666 6.911596 97.05735 2.297839e-07 0.002192422
  # Gpbp1   Gpbp1  0.5049153 5.977704 85.01982 4.886149e-07 0.002192422
  # Ugcg     Ugcg  0.4676207 6.344091 80.30188 6.739704e-07 0.002192422
  # Metrn   Metrn -0.7118352 5.860048 78.55347 7.624780e-07 0.002192422
  # Zfp948 Zfp948  0.5688350 3.651637 77.39104 8.287655e-07 0.002192422
  # Fmr1     Fmr1  0.6944720 5.896388 76.06021 9.129869e-07 0.002192422
  
  DEGstat <- DEGstat[order(DEGstat$ID, decreasing = F), ]
  head(DEGstat)
  #                          ID       logFC    logCPM          F      PValue        FDR
  # 0610009B22Rik 0610009B22Rik  0.03548662 3.5245657  0.1004021 0.756420389 0.86502052
  # 0610009E02Rik 0610009E02Rik -0.38217332 0.3013286  2.7970948 0.118514251 0.28527793
  # 0610009L18Rik 0610009L18Rik -0.47053843 0.6483960  7.3286299 0.018058369 0.09124405
  # 0610010F05Rik 0610010F05Rik  0.26305502 5.4776492 16.5780792 0.001343571 0.02254203
  # 0610010K14Rik 0610010K14Rik -0.29593765 0.7906806  2.3864680 0.146576111 0.32341622
  # 0610012G03Rik 0610012G03Rik -0.26672813 3.9017953 13.2762050 0.003011483 0.03333386
  dim(DEGstat) #16703     6
  
  
  dgeFullCPM <- cpm(dgeFull, prior.count=2, log=T) 
  head(dgeFullCPM)
  #          CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  # COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
  # CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
  # Mbp  13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
  # ND1  13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
  # ND5  12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
  # Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
  dim(dgeFullCPM) # 16703    14
  
  DEGsCPM <-  data.frame(ID = rownames(dgeFullCPM), dgeFullCPM)
  head(DEGsCPM)
  #        ID     CT1F     CT2F     CT3F     CT4F    DHT1F    DHT2F    DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  # COX1 COX1 14.49719 14.23604 14.26942 14.09648 14.21047 14.32967 14.21530 14.21573 13.98610 14.02834 14.12056 14.44912 14.31591 14.33175
  # CYTB CYTB 13.50180 13.26776 13.31711 13.01758 13.37885 13.53607 13.31100 13.19736 12.94628 12.92665 13.04406 13.60147 13.40101 13.38852
  # Mbp   Mbp 13.08640 12.99335 13.05390 12.62820 13.49804 13.49973 13.17451 13.27433 13.17795 13.01203 13.03769 12.83401 13.11146 13.35427
  # ND1   ND1 13.07044 12.88562 12.77548 12.44575 12.92941 13.18987 12.90415 12.45838 12.19204 12.22589 12.22680 13.23518 12.97265 12.95923
  # ND5   ND5 12.92590 12.75405 12.76329 12.56210 12.84010 13.09777 12.75937 12.54347 12.26806 12.24760 12.34311 13.20126 12.93608 12.71435
  # Fth1 Fth1 12.66612 12.52654 12.60475 12.78900 11.80393 11.99434 12.37277 12.83489 12.81325 12.91274 12.87935 12.75199 12.43895 12.47616
  dim(DEGsCPM)#  16703   15
  DEGsCPM <- DEGsCPM[order(DEGsCPM$ID, decreasing = F), ]
  head(DEGsCPM)
  #                         ID      CT1F      CT2F      CT3F       CT4F     DHT1F     DHT2F     DHT3F      CT1M      CT2M       CT3M      CT4M       DHT1M     DHT2M      DHT3M
  # 0610009B22Rik 0610009B22Rik 3.3267319 3.2889966 3.7131574  3.6338225 3.6424848 3.6638320 3.4858948 3.3368192 3.3731138 3.68190779 3.4511837  3.45351361 3.6395959  3.5531386
  # 0610009E02Rik 0610009E02Rik 0.1494915 0.4662663 0.1753633 -0.5622992 0.5107245 0.4909081 0.5145337 0.3761698 0.8093237 0.09139392 0.5083199  0.09431313 0.2408394 -0.1858660
  # 0610009L18Rik 0610009L18Rik 0.4200818 0.7406566 0.2941583  0.4788918 0.9736110 0.6790344 0.5145337 0.9618237 0.9853484 0.85521187 0.9800794 -0.13508851 0.4238192  0.4328421
  # 0610010F05Rik 0610010F05Rik 5.5497149 5.5320587 5.4783899  5.4256583 5.6863199 5.7140562 5.6436749 5.1764783 5.2653428 5.35964449 5.1259288  5.54069860 5.6359610  5.4006756
  # 0610010K14Rik 0610010K14Rik 0.6478472 0.8731971 0.6011559  0.9557145 0.6212767 0.5193323 0.8236557 1.4303555 0.9853484 0.69156560 1.0493971  0.13612002 0.7321061  0.6212411
  # 0610012G03Rik 0610012G03Rik 3.8912906 3.8484805 3.8093833  4.0480171 3.5788348 3.4916204 3.6157369 4.2235889 4.0473192 4.15224861 4.2402180  3.87678357 3.7187710  3.8391264
  tail(DEGstat)
  
  datCTFvsCTM <- cbind(DEGstat, DEGsCPM[, c(2:15)] )
  head(datCTFvsCTM)
  dim(datCTFvsCTM) #16703    20
  datCTFvsCTM[datCTFvsCTM$ID =="Nipal2", ]
  
  #           ID     logFC   logCPM        F       PValue       FDR     CT1F     CT2F    CT3F     CT4F    DHT1F   DHT2F  DHT3F     CT1M     CT2M     CT3M     CT4M    DHT1M    DHT2M    DHT3M
  #Nipal2 Nipal2 0.1397795 3.933548  2.67527    0.1260858 0.2961449 3.744199 3.785954 3.98154 3.573086 4.432566 4.34795 4.3563 3.641062 3.610329 3.724875 3.577799 3.845623 4.071208 3.948748
  
  
  write.table(x=datCTFvsCTM, file='AEGs_CTFvsCTM_logCPM_EdgeR.txt',
     row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
  
  
  #################################################################################################  
  # Creating the fusion of stat table and the logCPM only for only the FDR005 genes
  ################################################################################################# 
  
  DEGsCTFvsCTM <- subset(datCTFvsCTM, FDR <0.05)
  dim(DEGsCTFvsCTM) #2056   20
  
   write.table(x=DEGsCTFvsCTM, file='DEGsCTFvsCTM_logCPM_EdgeR.txt',
      row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
  
  
  
##################################################################################################################################################################  
# Creating the fusion of all stats  and the logCPM 
################################################################################################################################################################## 
  
   alldat <- cbind(datDHTFvsCTF[, c(1:6)], datDHTFvsCT[, c(1:6)], datDHTMvsCTM[, c(1:6)], datDHTMvsCT[, c(1:6)], datDHTFvsDHTM[, c(1:6)], datCTFvsCTM )
  dim(alldat) #16703    50
  
  colnames(alldat)
  # [1] "ID"     "logFC"  "logCPM" "F"      "PValue" "FDR"    "ID"     "logFC"  "logCPM" "F"      "PValue" "FDR"    "ID"     "logFC"  "logCPM" "F"      "PValue" "FDR"    "ID"     "logFC"  "logCPM"
  # [22] "F"      "PValue" "FDR"    "ID"     "logFC"  "logCPM" "F"      "PValue" "FDR"    "ID"     "logFC"  "logCPM" "F"      "PValue" "FDR"    "CT1F"   "CT2F"   "CT3F"   "CT4F"   "DHT1F"  "DHT2F" 
  # [43] "DHT3F"  "CT1M"   "CT2M"   "CT3M"   "CT4M"   "DHT2M"  "DHT3M"  "DHT4M" 
  
  colnames(alldat) <- c(
     "ID", "logFC_DHTFvsCTF",  "logCPM_DHTFvsCTF", "F_DHTFvsCTF", "PValue_DHTFvsCTF", "FDR_DHTFvsCTF",
     "ID", "logFC_DHTFvsCT",  "logCPM_DHTFvsCT", "F_DHTFvsCT", "PValue_DHTFvsCT", "FDR_DHTFvsCT",
     "ID", "logFC_DHTMvsCTM",  "logCPM_DHTMvsCTM", "F_DHTMvsCTM", "PValue_DHTMvsCTM", "FDR_DHTMvsCTM",
     "ID", "logFC_DHTMvsCT",  "logCPM_DHTMvsCT", "F_DHTMvsCT", "PValue_DHTMvsCT", "FDR_DHTMvsCT",
     "ID", "logFC_DHTFvsDHTM",  "logCPM_DHTFvsDHTM", "F_DHTFvsDHTM", "PValue_DHTFvsDHTM", "FDR_DHTFvsDHTM",
     "ID", "logFC_CTFvsCTM",  "logCPM_CTFvsCTM", "F_CTFvsCTM", "PValue_CTFvsCTM", "FDR_CTFvsCTM",
     "CT1F",   "CT2F",   "CT3F",   "CT4F",   "DHT1F",  "DHT2F", "DHT3F",  "CT1M",   "CT2M",   "CT3M",   "CT4M",   "DHT1M", "DHT2M",  "DHT3M") 
  colnames(alldat) 
  
  # [1] "ID"                "logFC_DHTFvsCTF"   "logCPM_DHTFvsCTF"  "F_DHTFvsCTF"       "PValue_DHTFvsCTF"  "FDR_DHTFvsCTF"     "ID"                "logFC_DHTFvsCT"    "logCPM_DHTFvsCT"  
  # [10] "F_DHTFvsCT"        "PValue_DHTFvsCT"   "FDR_DHTFvsCT"      "ID"                "logFC_DHTMvsCTM"   "logCPM_DHTMvsCTM"  "F_DHTMvsCTM"       "PValue_DHTMvsCTM"  "FDR_DHTMvsCTM"    
  # [19] "ID"                "logFC_DHTMvsCT"    "logCPM_DHTMvsCT"   "F_DHTMvsCT"        "PValue_DHTMvsCT"   "FDR_DHTMvsCT"      "ID"                "logFC_DHTFvsDHTM"  "logCPM_DHTFvsDHTM"
  # [28] "F_DHTFvsDHTM"      "PValue_DHTFvsDHTM" "FDR_DHTFvsDHTM"    "ID"                "logFC_CTFvsCTM"    "logCPM_CTFvsCTM"   "F_CTFvsCTM"        "PValue_CTFvsCTM"   "FDR_CTFvsCTM"     
  # [37] "CT1F"              "CT2F"              "CT3F"              "CT4F"              "DHT1F"             "DHT2F"             "DHT3F"             "CT1M"              "CT2M"      
  # 
  write.table(x=alldat, file='AEGs_alldat_logCPM_EdgeR.txt',
     row.names=F, col.names=T, quote=F, sep="\t", dec='.') 
  
  dim(alldat)
  head(alldat)
#################################################################################################  
# Creating Curated lists for DHTFvsCTF
#################################################################################################  
  alldat <- read.table("AEGs_alldat_logCPM_EdgeR.txt",sep="\t", header=TRUE,  fill=TRUE, quote = "")
  head(alldat)
  dim(alldat) #16703    50
  
  
  
  #filter for FDR_DHTFvsCTF <0.05
  DHTFvsCTF_DEGs <- subset(alldat, FDR_DHTFvsCTF <0.05)
  dim(DHTFvsCTF_DEGs) #7470   50
  head(DHTFvsCTF_DEGs)
  
  DHTMvsCTM_DEGs <- subset(alldat, FDR_DHTMvsCTM <0.05)
  dim(DHTMvsCTM_DEGs) #3781   50
  #
 
  

  OL_curated2 = read.table('OLcurated_445.txt', sep = "\t", header=TRUE, row.names=1)
  head(OL_curated2)
  #       Specif Prolif Migrat Survival Differ Myelin Remyel
  # Aatk       0      0      0        0      1      0      0
  # Abcd1      0      0      0        1      0      0      0
  # Abhd6      0      0      0        0      0      0     -1
  # Abl1       0      1      0        0      0      0      0
  # Acaca      0      0      0        0      0      1      0
  # Acox1      0      0      0        0      0      1      0
  dim(OL_curated2) #444   7
  
  
  
  DHTFvsCTF_DEGs_OLcurated = OL_curated2[rownames(OL_curated2) %in% DHTFvsCTF_DEGs$ID , ]
  head(DHTFvsCTF_DEGs_OLcurated)
  #        Specif Prolif Migrat Survival Differ Myelin Remyel
  # Aatk        0      0      0        0      1      0      0
  # Abcd1       0      0      0        1      0      0      0
  # Abhd6       0      0      0        0      0      0     -1
  # Acss2       0      0      0        1      0      0      1
  # Adam15      0      0      0        2      0      0      0
  # Adam17      0      2      0        2      1      1      0
  dim(DHTFvsCTF_DEGs_OLcurated) # 210   7
  
  
  DHTFvsCTF_curated_logFC <- DHTFvsCTF_DEGs[DHTFvsCTF_DEGs$ID %in% rownames(DHTFvsCTF_DEGs_OLcurated),  ]
  head(DHTFvsCTF_curated_logFC)
  dim(DHTFvsCTF_curated_logFC) #210  50

  DHTFvsCTF_curated_logFC2 <- subset(DHTFvsCTF_curated_logFC, select = c(1:6))
  DHTFvsCTF_curated_logFC2 <- DHTFvsCTF_curated_logFC2[order(DHTFvsCTF_curated_logFC2$ID, decreasing = F), ]
  head(DHTFvsCTF_curated_logFC2)
   #         ID logFC_DHTFvsCTF logCPM_DHTFvsCTF F_DHTFvsCTF PValue_DHTFvsCTF FDR_DHTFvsCTF
  # 576   Aatk       0.5015684         8.368987   14.089697     0.0024440935   0.011660543
  # 609  Abcd1      -0.8766540         5.262127   23.007199     0.0003725117   0.005735475
  # 639  Abhd6       0.4463744         5.075470   16.994018     0.0012222148   0.008132975
  # 736  Acss2       0.4555666         6.234625    7.624692     0.0164444377   0.039906990
  # 773 Adam15      -0.4260926         7.189044   17.470342     0.0010985355   0.007778227
  # 774 Adam17      -1.2036608         6.453859   23.471667     0.0004054137   0.005873048
  dim(DHTFvsCTF_curated_logFC2) # 210  6
  
  DHTFvsCTF_DEGs_OLcurated_logFC <- DHTFvsCTF_DEGs_OLcurated[ , c(1:6)]*DHTFvsCTF_curated_logFC2$logFC_DHTFvsCTF
  head(DHTFvsCTF_DEGs_OLcurated_logFC)
  #        Specif    Prolif Migrat   Survival     Differ    Myelin
  # Aatk        0  0.000000      0  0.0000000  0.5015684  0.000000
  # Abcd1       0  0.000000      0 -0.8766540  0.0000000  0.000000
  # Abhd6       0  0.000000      0  0.0000000  0.0000000  0.000000
  # Acss2       0  0.000000      0  0.4555666  0.0000000  0.000000
  # Adam15      0  0.000000      0 -0.8521852  0.0000000  0.000000
  # Adam17      0 -2.407322      0 -2.4073215 -1.2036608 -1.203661
  

  DHTFvsCTF_stats_OLcurated_logFC <- cbind(DHTFvsCTF_DEGs_OLcurated_logFC, DHTFvsCTF_curated_logFC2)
  dim(DHTFvsCTF_stats_OLcurated_logFC) #210  12
  head(DHTFvsCTF_stats_OLcurated_logFC)
  #        Specif    Prolif Migrat   Survival     Differ    Myelin     ID logFC_DHTFvsCTF logCPM_DHTFvsCTF F_DHTFvsCTF PValue_DHTFvsCTF FDR_DHTFvsCTF
  # Aatk        0  0.000000      0  0.0000000  0.5015684  0.000000   Aatk       0.5015684         8.368987   14.089697     0.0024440935   0.011660543
  # Abcd1       0  0.000000      0 -0.8766540  0.0000000  0.000000  Abcd1      -0.8766540         5.262127   23.007199     0.0003725117   0.005735475
  # Abhd6       0  0.000000      0  0.0000000  0.0000000  0.000000  Abhd6       0.4463744         5.075470   16.994018     0.0012222148   0.008132975
  # Acss2       0  0.000000      0  0.4555666  0.0000000  0.000000  Acss2       0.4555666         6.234625    7.624692     0.0164444377   0.039906990
  # Adam15      0  0.000000      0 -0.8521852  0.0000000  0.000000 Adam15      -0.4260926         7.189044   17.470342     0.0010985355   0.007778227
  # Adam17      0 -2.407322      0 -2.4073215 -1.2036608 -1.203661 Adam17      -1.2036608         6.453859   23.471667     0.0004054137   0.005873048
  
  
  
  write.table(x=DHTFvsCTF_stats_OLcurated_logFC, file='DHTFvsCTF_stats_OLcurated_logFC.txt', col.names=T, row.names=T, 
    sep="\t", dec='.')
  
  
  
  
  ############################################################ 
  #Creating Curated lists for DHTMvsCTM
  ############################################################ 

  DHTMvsCTM_DEGs_OLcurated = OL_curated2[rownames(OL_curated2) %in% DHTMvsCTM_DEGs$ID , ]
  head(DHTMvsCTM_DEGs_OLcurated)
  #        Specif Prolif Migrat Survival Differ Myelin Remyel
  # Abcd1       0      0      0        1      0      0      0
  # Abl1        0      1      0        0      0      0      0
  # Acox1       0      0      0        0      0      1      0
  # Adam15      0      0      0        2      0      0      0
  # Adam22      0      0      0        0      1      0      0
  # Ago2        0      0      0        0      0      1      0
  dim(DHTMvsCTM_DEGs_OLcurated) # 94   7
  
  
  DHTMvsCTM_curated_logFC <- DHTMvsCTM_DEGs[DHTMvsCTM_DEGs$ID %in% rownames(DHTMvsCTM_DEGs_OLcurated),  ]
  dim(DHTMvsCTM_curated_logFC) #94  50
  
  DHTMvsCTM_curated_logFC2 <- subset(DHTMvsCTM_curated_logFC, select = c(13:18))
  DHTMvsCTM_curated_logFC2 <- DHTMvsCTM_curated_logFC2[order(DHTMvsCTM_curated_logFC2$ID, decreasing = F), ]
  head(DHTMvsCTM_curated_logFC2)
  #       ID.2 logFC_DHTMvsCTM logCPM_DHTMvsCTM F_DHTMvsCTM PValue_DHTMvsCTM FDR_DHTMvsCTM
  # 609  Abcd1      -0.6068864         5.262127   11.301826      0.005258009   0.031476684
  # 646   Abl1      -0.3967174         6.551083   29.116638      0.000125577   0.005805802
  # 714  Acox1      -0.1584439         7.514612    8.838281      0.010870613   0.048810679
  # 773 Adam15      -0.3647606         7.189044   12.856839      0.003363087   0.024508573
  # 779 Adam22       0.4140273         7.472845   11.202796      0.005305608   0.031615971
  # 926   Ago2       0.1972504         7.270222    8.935250      0.010537115   0.047930673
  dim(DHTMvsCTM_curated_logFC2) # 94  6


  
  DHTMvsCTM_DEGs_OLcurated_logFC <- DHTMvsCTM_DEGs_OLcurated[ , c(1:6)]*DHTMvsCTM_curated_logFC2$logFC_DHTMvsCTM
  head(DHTMvsCTM_DEGs_OLcurated_logFC)
  #       Specif     Prolif Migrat   Survival    Differ     Myelin
  # Abcd1       0  0.0000000      0 -0.6068864 0.0000000  0.0000000
  # Abl1        0 -0.3967174      0  0.0000000 0.0000000  0.0000000
  # Acox1       0  0.0000000      0  0.0000000 0.0000000 -0.1584439
  # Adam15      0  0.0000000      0 -0.7295212 0.0000000  0.0000000
  # Adam22      0  0.0000000      0  0.0000000 0.4140273  0.0000000
  # Ago2        0  0.0000000      0  0.0000000 0.0000000  0.1972504
  
  
  DHTMvsCTM_stats_OLcurated_logFC <- cbind(DHTMvsCTM_DEGs_OLcurated_logFC, DHTMvsCTM_curated_logFC2)
  dim(DHTMvsCTM_stats_OLcurated_logFC) #210  12
  head(DHTMvsCTM_stats_OLcurated_logFC)
  #        Specif     Prolif Migrat   Survival    Differ     Myelin   ID.2 logFC_DHTMvsCTM logCPM_DHTMvsCTM F_DHTMvsCTM PValue_DHTMvsCTM FDR_DHTMvsCTM
  # Abcd1       0  0.0000000      0 -0.6068864 0.0000000  0.0000000  Abcd1      -0.6068864         5.262127   11.301826      0.005258009   0.031476684
  # Abl1        0 -0.3967174      0  0.0000000 0.0000000  0.0000000   Abl1      -0.3967174         6.551083   29.116638      0.000125577   0.005805802
  # Acox1       0  0.0000000      0  0.0000000 0.0000000 -0.1584439  Acox1      -0.1584439         7.514612    8.838281      0.010870613   0.048810679
  # Adam15      0  0.0000000      0 -0.7295212 0.0000000  0.0000000 Adam15      -0.3647606         7.189044   12.856839      0.003363087   0.024508573
  # Adam22      0  0.0000000      0  0.0000000 0.4140273  0.0000000 Adam22       0.4140273         7.472845   11.202796      0.005305608   0.031615971
  # Ago2        0  0.0000000      0  0.0000000 0.0000000  0.1972504   Ago2       0.1972504         7.270222    8.935250      0.010537115   0.047930673
  
  
  write.table(x=DHTMvsCTM_stats_OLcurated_logFC, file='DHTMvsCTM_stats_OLcurated_logFC.txt', col.names=T, row.names=T, 
    sep="\t", dec='.')
  
  
  ############################################################ 
  # microglia and astrocyte subsets
  ############################################################ 
  
  genes <- c("Nos2", "Tnf", "Ccl2", "Arg1", "Cd206", "Igf1", 
    "C1qb", "C1qc", "Cx3cr1", "Mertk", "Cd33", "Apoe", "Fcgr2b", "Trem2", "Axl", "Tlr4", "Rxra", "Cd68")
  
  microglial_genes <- alldat[alldat$ID %in% genes, ]
  dim(microglial_genes) #17 50
  head(microglial_genes)
  
  write.table(microglial_genes, file='microglial_genes_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T)
  
  microglia_homeostatic <- c("Hexb", "Cst3", "Cx3cr1", "Ctsd", "Csf1r", "Ctss", "Sparc", "Tmsb4x", "P2ry12", "C1qa", "C1qb", "TMEM119", "P2RY13", 
    "SLC2A5", "1700017B05Rik", "Ccr5", "Cd164", "Cmtm6", "Ecscr", "Fcrls", "Gpr34", "Ifngr1", "Laptm5", "Ldhb", "Lgmn", 
    "Lpcat2", "Lrrc3", "Mgat4a", "Olfml3", "Ptgs1", "Selplg", "Siglech", "Sirpa", "Srgap2", "Ubqln1", "Upk1b", "Vdac1", 
    "Vsir", "Bhlhe41", "Sall1", "Serpine2")
  
  
  microglia_homeostatic <- alldat[ alldat$ID %in% microglia_homeostatic, ]
  dim(microglia_homeostatic) #37 50
  head(microglia_homeostatic)
  
  write.table(microglia_homeostatic, file='microglia_homeostatics_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T)
  
  microglia_DAM <- c("Tyropb", "Ctsb", "Ctsd", "Apoe", "B2m", "Fth1", "Lyz2", "Trem2", "Axl", "Csl7", "Ctsl", "Lpl", "Cd9", "Csf1", "Ccl6", "Itgax", 
    "Clec7a", "Lilrb4", "Timp2")
  microglia_DAM <- alldat[ alldat$ID %in% microglia_DAM, ]
  dim(microglia_DAM) #16 50
  head(microglia_DAM)
  
  write.table(microglia_DAM, file='microglia_DAM_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T)
  
  
  microglia_WAM <- c("Anxa2", "Anxa5", "Apoe", "Atp6v0c", "B2m", "C1qb", "Capg", "Cd52", "Cd63", "Cd63-ps", "Cd74", "Clec7a", "Crip1", "Cst7", "Ctsb", "Ctss", 
    "Ctsz", "Cybb", "Fam20c", "Fth1", "Ftl1", "Ftl1-ps1", "Ftl1-ps2", "Gm12164", "Gm7541", "H2-D1", "H2-K1", "Ifitm3", "Lgals1", "Lgals3", "Lyz2",
    "Mir692-1", "Spp1", "Tspo", "Vim")
  
  microglia_WAM<- alldat[ alldat$ID %in% microglia_WAM, ]
  dim(microglia_WAM) #29 50
  head(microglia_WAM)
  
  write.table(microglia_WAM, file='microglia_WAM_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T)
  
  
  microglia_activated <- c("CT010467.1", "Eef1a1", "Eef1b2", "Eef2", "Fau", "Gm10076", "Gm10269", "Gm10443", "Gm10736", "Gm14303", "Gm2000", 
    "Gm23935", "Gm3511", "Gm5735", "Gm7618", "Lars2", "Mir6236", "mt-Co1", "mt-Cytb", "mt-Rnr2", "Rpl10", "Rpl13", "Rpl14", 
    "Rpl17", "Rpl18a", "Rpl21", "Rpl27a", "Rpl28", "Rpl31", "Rpl32", "Rpl35", "Rpl36", "Rpl37", "Rpl37a", "Rpl38", "Rpl41", 
    "Rpl5", "Rpl9", "Rplp0", "Rplp1", "Rplp2", "Rps10-ps2", "Rps11", "Rps14", "Rps16", "Rps17", "Rps18", "Rps19", "Rps23", "Rps25", 
    "Rps27", "Rps27a", "Rps28", "Rps29", "Rps3", "Rps6", "Rps9", "Rpsa", "Tmsb4x")
  microglia_activated <- alldat[ alldat$ID %in% microglia_activated, ]
  dim(microglia_activated) #44 50
  head(microglia_activated)
  
  write.table(microglia_activated, file='microglia_activated_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T)
  
  
  microglia_demyelination <- c("Apoe", "Axl", "Igf1", "Lyz2", "Itgax", "Gpnmb", "Apoc1")
  microglia_demyelination <- alldat[ alldat$ID %in% microglia_demyelination, ]
  dim(microglia_demyelination) #7 50  
  head(microglia_demyelination)
  
  write.table(microglia_demyelination, file='microglia_demyelination_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T)
  
  
  microglia_remyelination <- c("Fam20c", "Cst7", "Ccl6", "Fn1", "Ank", "Psat1", "Spp1")
  microglia_remyelination <- alldat[ alldat$ID %in% microglia_remyelination, ]
  dim(microglia_remyelination) #7 50
  head(microglia_remyelination)
  
  write.table(microglia_remyelination, file='microglia_remyelination_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T)
  
  
  astroglia_homeostatic1 <- c("Hsd11b1", "Chst1", "Kcnip3", "Phactr3", "Tspan7", "Chrdl1", "Fam212b", "St6galnac5", "Cth", "Gabbr2", "Efhd2", "AW011738", "Grm3",
    "0610040J01Rik", "Abcb9", "Ptprz1", "Akr1b10", "Timp4", "Eps8", "Slco1c1", "Aplp1", "Slc7a10", "Grm5", "Dhcr7", "Ccnd1", "S100b", "Gpc5",
    "Cadm1", "Arpp21", "Rnf215", "Slc13a5", "Etv4", "Dbx2", "Olig2", "Cbr3", "Hbegf", "Cdo1", "Vldlr")
  
  astroglia_homeostatic1 <- alldat[ alldat$ID %in% astroglia_homeostatic1, ]
  dim(astroglia_homeostatic1) #37 50
  head(astroglia_homeostatic1)
  
  write.table(astroglia_homeostatic1 , file='astroglia_homeostatic1_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T)
  
  
  astroglia_homeostatic2 <- c("Prex2", "Col9a1", "Arhgef4", "Bmpr2", "Nrp2", "Mgat5", "Adora1", "Rgs2", "Glul", "Tnr", "Myoc", "Ccdc3", "Nsmf", "Ggta1", "Nr4a2", 
    "Pkp4", "Ssfa2", "Slc1a2", "Cd44", "Serf2", "Slc20a1", "Btbd3", "Ptprt", "Matn4", "Elmo2", "Sulf2", "Usp9x", "Slc6a8", "Plp1", "Gpm6b",
    "Trim2", "Slc16a1", "Unc5c", "Bmpr1b", "Grin3a", "Rspo1", "Man1c1", "Ptchd2", "Ski", "Agrn", "Prkag2", "Hopx", "Miat", "Hspb8", "Clip2",
    "Tsc22d4", "D630045J12Rik", "Gadd45a", "Atoh8", "Tmem150a", "Tgoln1", "Ctnna2", "Slc6a6", "Adamts9", "A2m", "Kcna6", "Zfp36", "Clip3", 
    "Nav2", "Arhgef17", "Tab2", "Marcks", "Sgpl1", "Slc7a2", "Mfap3l", "Ncan", "Adgrl1", "Zfp423", "Mt3", "Mt2", "Mt1", "Cadps", "Opcml",
    "Sorl1", "Cryab", "Islr", "6030419C18Rik", "Smad6", "Igdcc4", "Pth1r", "Tmie", "Csrnp1", "Tmem158", "Nsg2", "Wwc1", "Col23a1", "Slc36a2",
    "Serpinf1", "Taok1", "Nog", "Stat3", "Dusp3", "Prkca", "Fasn", "Aldh5a1", "Mylip", "Msx2", "Erbb2ip", "Map3k1", "Id2", "Rgs6", "Lifr", 
    "Ank", "Enpp2", "Plec", "Maff", "Shisa8", "Sept3", "Fam19a5", "Slc38a1", "Snai2", "Nfkbiz", "Robo2", "Adamts1", "Tulp4", "Paqr4", "Fgd2", 
    "H2-DMa", "H2-DMb1", "H2-Ab1", "H2-Aa", "H2-Eb1", "Ddr1", "Npc1", "Aqp4", "Nrep", "Cd74", "Cidea", "Smad7", "Fth1", "Aldh1a1", "Kank1", "Scd2")
  
  
  astroglia_homeostatic2 <- alldat[ alldat$ID %in% astroglia_homeostatic2, ]
  dim(astroglia_homeostatic2) #128  50
  head(astroglia_homeostatic2)
  
  write.table(astroglia_homeostatic2 , file='astroglia_homeostatic2_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T) 
  
  astroglia_activation <- c("Bgn", "S100a4", "S100a11", "S100a10", "Igfbp7", "Ifitm3", "C1ql1", "Timp2", "S1pr3", "Ntrk2", "Actn1", "Serpina3n", "Ccdc74a", "C4b")
  
  astroglia_activation  <- alldat[ alldat$ID %in% astroglia_activation, ]
  dim(astroglia_activation) #14 50
  head(astroglia_activation)
  
  write.table(astroglia_activation , file='astroglia_activation_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T) 
  
  astroglia_A1 <-  c("H2.T23", "Serping1", "H2.D1", "Ggta1", "Iigp1", "Gpb2", "Fbln5", "Ugt1a", "Fkbp5", "Psmb8", "Srgn", "Amigo2")
  astroglia_A1  <- alldat[ alldat$ID %in% astroglia_A1, ]
  dim(astroglia_A1) #8 50
  head(astroglia_A1)
  write.table(astroglia_A1 , file='astroglia_A1_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T) 
  
  
  astroglia_A2 <-  c("Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd109", "Ptgs2", "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14")
  astroglia_A2  <- alldat[ alldat$ID %in% astroglia_A2, ]
  dim(astroglia_A2) #12 50
  head(astroglia_A2)
  write.table(astroglia_A2 , file='astroglia_A2_alldat.txt', sep="\t",
    quote=F, col.names=T, row.names=T) 
  
  
  
  
