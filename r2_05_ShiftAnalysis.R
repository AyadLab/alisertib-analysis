library(Seurat)
library(GSEABase)
library(ggsignif)
library(ggpubr)
library(Rmisc)
library(cowplot)
library(singscore)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(cowplot)

# For raincloud plots...
library(scrabble)
library(readr)
library(tidyr)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(PupillometryR)
library(singscore)
library(viridis)
library(future)

# change the current plan to access parallelization
plan(strategy = "sequential")
# plan()
options(future.globals.maxSize = 40000 * 1024^2)
obj.integrated <- readRDS("ali_obj.int1.RDS")

# Define raincloud theme: from https://neuroconscience.wordpress.com/2018/03/15/introducing-raincloud-plots/
raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

##############################################################################

# obj.integrated <- readRDS(file = "AlisertibAnalysisPipeline/04_output/04_CellTagHu_CytoTRACE.RDS")
# 
# Scale RNA assay
DefaultAssay(obj.integrated) <- "RNA"
obj.integrated <- NormalizeData(obj.integrated, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
obj.integrated <- CellCycleScoring(obj.integrated, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)


##----
### setting up vector of genes to keep for scaling
# variableFeatures
variableFeatures <- obj.integrated@assays$RNA@var.features
# suvaFeatures
suvasigs <- read.csv(file = "SuvaMetaSignatures.csv", header = T, skip = 4)
SuvaNPC1selected <- suvasigs$NPC1
SuvaNPC2selected <- suvasigs$NPC2
SuvaOPCselected <- suvasigs$OPC
SuvaACselected <- suvasigs$AC
SuvaMES1selected <- suvasigs$MES1
SuvaMES2selected <- suvasigs$MES2
suva_genes <- c(SuvaNPC1selected, SuvaNPC2selected, SuvaOPCselected, SuvaACselected, SuvaMES1selected, SuvaMES2selected)
#L1000Features
l1000signatures <- read.delim(file = "matPH3_2_1_0.2_0.3_L1000_Batch2017_Regina_removed.txt", row.names = 1)
L1000Features <- rownames(l1000signatures)

toScale <- unique(c(L1000Features, suva_genes, variableFeatures))
##----
obj.integrated <- ScaleData(obj.integrated, 
                            assay = "RNA", 
                            features = toScale, # only scale features of interest
                            vars.to.regress = c("G2M.Score", "S.Score"))
saveRDS(obj.integrated,"interestgenes.rds")
dim(obj.integrated@assays$RNA@scale.data)
dim(obj.integrated@assays$RNA@data)
# obj.integrated@assays$RNA@scale.data

saveRDS(obj.integrated, file = "alisertib_r2_cellCycleRegressedInt.RDS")
obj.integrated <- readRDS(file = "alisertib_r2_cellCycleRegressedInt.RDS")
# Split into seperate objects

head(obj.integrated$arm)

obj.integrated <- RenameCells(object = obj.integrated, new.names = paste0(obj.integrated$arm, "_", colnames(x = obj.integrated)))
Idents(obj.integrated) <- obj.integrated$arm
Alisertib <- subset(obj.integrated, idents = "alisertib")
metadata <- Alisertib@meta.data
subtypes_alisertib <- metadata_alisertib$arm


head(metadata)
table(subtypes_alisertib)


#saveRDS(Alisertib,"Alisertib.rds")
DMSO <- subset(obj.integrated, idents = "vehicle")
#saveRDS(DMSO,"dmso.rds)

##############################################################################
##############################################################################
##############################################################################

# singscore

############################

# # Trimmed signatures from Neftel et al...
# 
# SuvaNPC1selected <- c("SOX4", "MARCKSL1", "STMN1", "TCF12", "BEX1", "MLLT11", "MEST", "PFN2", "CHD7", "TUBA1A", "PCBP4", "AMOTL2", "DBN1", "SOX11")
# SuvaNPC2selected <- c("RND3", "MLLT11", "SOX11", "NREP", "FNBP1L", "SOX4", "MAP1B", "RBFOX2", "STMN1", "DPYSL3", "ATP1B1", "UCHL1", "DDAH2", "TUBB2A", "LBH", "TCF4", "NFIB",
#                       "DBN1", "NFIX", "CEP170", "BLCAP")
# SuvaOPCselected <- c("PSAT1", "SIRT2", "THY1", "VCAN", "DBI", "LIMA1", "SCD5", "FGF12", "RNF13", "ALCAM", "PGRMC1", "FABP5", "RAB31", "EPB41L2")
# SuvaACselected <- c("CST3", "PON2", "SPARC", "RAMP1", "DBI", "CLU", "PMP22", "S100A16", "PRCP", "PLTP", "F3", "RAB31", "ANXA5")
# SuvaMES1selected <- c("ANXA2", "ANXA1", "CD44", "VIM", "MT2A", "C1S", "NAMPT", "C1R", "SOD2", "IFITM3", "TIMP1", "S100A11", "MT1X", "S100A10", "FN1", "LGALS1", "S100A16", "CLIC1", "MGST1", "RCAN1", "TAGLN2", "NPC2", "SERPING1", "EMP1", "CTSB", "LGALS3", "MT1E", "EMP3", "ACTN1", "PRDX6", "IGFBP7", "SERPINE1", "PLP2", "CLIC4", "GSN", "NNMT", "TUBA1C", "GJA1", "TNFRSF1A", "WWTR1")
# SuvaMES2selected <- c("HILPDA", "ADM", "DDIT3", "NDRG1", "HERPUD1", "DNAJB9", "AKAP12", "SQSTM1", "MT1X", "ATF3", "NAMPT", "NRN1", "SLC2A1", "BNIP3", "LGALS3", "INSIG2", "PPP1R15A", "VIM", "PLOD2", "GBE1", "SLC2A3", "FTL", "WARS", "XPOT", "HSPA5", "ANXA2", "EPAS1", "LDHA", "P4HA1", "SERTAD1", "PFKP", "PGK1", "EGLN3", "SLC6A6", "BNIP3L", "RPL21", "TRAM1", "UFM1", "ASNS", "GOLT1B", "SLC39A14", "HSPA9")
# SuvaG2_Mselected <- c("CCNB1", "CDC20", "CCNB2", "PLK1", "CCNA2", "CKAP2", "KNSTRN", "RACGAP1", "TROAP", "KIF2C", "AURKA", "CENPF", "KPNA2", "ECT2", "BUB1", "CDCA8", "TACC3", "TTK", "TUBA1C", "NCAPD2", "ARL6IP1", "MZT1", "ANP32E", "KIF11", "TUBB4B", "SMC4", "CDC25B", "REEP4", "FOXM1", "TMPO", "GPSM2", "HMGB3", "RANGAP1", "H2AFZ")

# What if we use full suva sigs...

suvasigs <- read.csv(file = "SuvaMetaSignatures.csv", header = T, skip = 4)
suvasigs

SuvaNPC1selected <- suvasigs$NPC1
SuvaNPC2selected <- suvasigs$NPC2
SuvaOPCselected <- suvasigs$OPC
SuvaACselected <- suvasigs$AC
SuvaMES1selected <- suvasigs$MES1
SuvaMES2selected <- suvasigs$MES2

# Extract data from seurat obj

DefaultAssay(obj.integrated) <- "RNA"
intData <- GetAssayData(obj.integrated, assay = "RNA", slot = "scale.data")
datfram <- as.data.frame(intData)
rankData <- singscore::rankGenes(datfram)

# NPC1
NPC1GS <- GSEABase::GeneSet()
NPC1GS@geneIds <- as.character(SuvaNPC1selected)
ScoredNPC1int <- singscore::simpleScore(rankData = rankData, upSet = NPC1GS)
# for (i in 1:length(ScoredNPC1int$TotalScore)){
#   ScoredNPC1int$Arm[i] <- unlist(strsplit(rownames(ScoredNPC1int)[i], split = "_", fixed = TRUE))[1]
# }
# Alisertibscores <- subset(ScoredNPC1int, ScoredNPC1int$Arm == "alisertib")
# DMSOscores <- subset(ScoredNPC1obj.integrated, ScoredNPC1int$Arm == "vehicle")
# singNPC1box <- ggpubr::ggboxplot(ScoredNPC1obj.integrated, x = "Arm", y = "TotalScore",
#                          color = "Arm", palette = c("aquamarine", "maroon"),
#                          order = c("vehicle", "alisertib"),
#                          ylab = "NPC1 enrichment score", xlab = "Treatment")
# plot(singNPC1box + ggpubr::stat_compare_means(method = "t.test"))
# NPC1ttestres <- t.test(TotalScore ~ Arm, data = ScoredNPC1int, paired = FALSE)

# NPC2
NPC2GS <- GeneSet()
NPC2GS@geneIds <- as.character(SuvaNPC2selected)
ScoredNPC2int <- singscore::simpleScore(rankData = rankData, upSet = NPC2GS)
# for (i in 1:length(ScoredNPC2int$TotalScore)){
#   ScoredNPC2int$Arm[i] <- unlist(strsplit(rownames(ScoredNPC2int)[i], split = "_", fixed = TRUE))[1]
# }
# Alisertibscores <- subset(ScoredNPC2obj.integrated, ScoredNPC2int$Arm == "Alisertib")
# DMSOscores <- subset(ScoredNPC2obj.integrated, ScoredNPC2int$Arm == "vehicle")
# singNPC2box <- ggpubr::ggboxplot(ScoredNPC2obj.integrated, x = "Arm", y = "TotalScore",
#                          color = "Arm", palette = c("aquamarine", "maroon"),
#                          order = c("vehicle", "alisertib"),
#                          ylab = "NPC2 enrichment score", xlab = "Treatment")
# plot(singNPC2box + ggpubr::stat_compare_means(method = "t.test"))
# NPC2ttestres <- t.test(TotalScore ~ Arm, data = ScoredNPC2obj.integrated, paired = FALSE)

# OPC
OPCGS <- GeneSet()
OPCGS@geneIds <- as.character(SuvaOPCselected)
ScoredOPCint <- singscore::simpleScore(rankData = rankData, upSet = OPCGS)
# for (i in 1:length(ScoredOPCint$TotalScore)){
#   ScoredOPCint$Arm[i] <- unlist(strsplit(rownames(ScoredOPCint)[i], split = "_", fixed = TRUE))[1]
# }
# Alisertibscores <- subset(ScoredOPCobj.integrated, ScoredOPCint$Arm == "Alisertib")
# DMSOscores <- subset(ScoredOPCobj.integrated, ScoredOPCint$Arm == "vehicle")
# singOPCbox <- ggpubr::ggboxplot(ScoredOPCobj.integrated, x = "Arm", y = "TotalScore",
#                         color = "Arm", palette = c("aquamarine", "maroon"),
#                         order = c("vehicle", "alisertib"),
#                         ylab = "OPC enrichment score", xlab = "Treatment")
# plot(singOPCbox + ggpubr::stat_compare_means(method = "t.test"))
# OPCttestres <- t.test(TotalScore ~ Arm, data = ScoredOPCobj.integrated, paired = FALSE)

# AC
ACGS <- GeneSet()
ACGS@geneIds <- as.character(SuvaACselected)
ScoredACint <- singscore::simpleScore(rankData = rankData, upSet = ACGS)
# for (i in 1:length(ScoredACint$TotalScore)){
#   ScoredACint$Arm[i] <- unlist(strsplit(rownames(ScoredACint)[i], split = "_", fixed = TRUE))[1]
# }
# Alisertibscores <- subset(ScoredACobj.integrated, ScoredACint$Arm == "Alisertib")
# DMSOscores <- subset(ScoredACobj.integrated, ScoredACint$Arm == "vehicle")
# singACbox <- ggpubr::ggboxplot(ScoredACobj.integrated, x = "Arm", y = "TotalScore",
#                        color = "Arm", palette = c("aquamarine", "maroon"),
#                        order = c("vehicle", "alisertib"),
#                        ylab = "AC enrichment score", xlab = "Treatment")
# plot(singACbox + ggpubr::stat_compare_means(method = "t.test"))
# ACttestres <- t.test(TotalScore ~ Arm, data = ScoredACobj.integrated, paired = FALSE)

# MES1
MES1GS <- GeneSet()
MES1GS@geneIds <- as.character(SuvaMES1selected)
ScoredMES1int <- singscore::simpleScore(rankData = rankData, upSet = MES1GS)
# for (i in 1:length(ScoredMES1int$TotalScore)){
#   ScoredMES1int$Arm[i] <- unlist(strsplit(rownames(ScoredMES1int)[i], split = "_", fixed = TRUE))[1]
# }
# Alisertibscores <- subset(ScoredMES1obj.integrated, ScoredMES1int$Arm == "Alisertib")
# DMSOscores <- subset(ScoredMES1obj.integrated, ScoredMES1int$Arm == "vehicle")
# singMES1box <- ggpubr::ggboxplot(ScoredMES1obj.integrated, x = "Arm", y = "TotalScore",
#                          color = "Arm", palette = c("aquamarine", "maroon"),
#                          order = c("vehicle", "alisertib"),
#                          ylab = "MES1 enrichment score", xlab = "Treatment")
# plot(singMES1box + ggpubr::stat_compare_means(method = "t.test"))
# MES1ttestres <- t.test(TotalScore ~ Arm, data = ScoredMES1obj.integrated, paired = FALSE)

# MES2
MES2GS <- GeneSet()
MES2GS@geneIds <- as.character(SuvaMES2selected)
ScoredMES2int <- singscore::simpleScore(rankData = rankData, upSet = MES2GS)
# for (i in 1:length(ScoredMES2int$TotalScore)){
#   ScoredMES2int$Arm[i] <- unlist(strsplit(rownames(ScoredMES2int)[i], split = "_", fixed = TRUE))[1]
# }
# Alisertibscores <- subset(ScoredMES2obj.integrated, ScoredMES2int$Arm == "Alisertib")
# DMSOscores <- subset(ScoredMES2obj.integrated, ScoredMES2int$Arm == "vehicle")
# singMES2box <- ggpubr::ggboxplot(ScoredMES2obj.integrated, x = "Arm", y = "TotalScore",
#                          color = "Arm", palette = c("aquamarine", "maroon"),
#                          order = c("vehicle", "alisertib"),
#                          ylab = "MES2 enrichment score", xlab = "Treatment")
# plot(singMES2box + ggpubr::stat_compare_means(method = "t.test"))
# MES2ttestres <- t.test(TotalScore ~ Arm, data = ScoredMES2obj.integrated, paired = FALSE)

# # G2_M
# G2_MGS <- GeneSet()
# G2_MGS@geneIds <- as.character(SuvaG2_Mselected)
# ScoredG2_Mobj.integrated <- singscore::simpleScore(rankData = rankData, upSet = G2_MGS)
# for (i in 1:length(ScoredG2_Mint$TotalScore)){
#   ScoredG2_Mint$Arm[i] <- unlist(strsplit(rownames(ScoredG2_Mint)[i], split = "_", fixed = TRUE))[1]
# }
# Alisertibscores <- subset(ScoredG2_Mobj.integrated, ScoredG2_Mint$Arm == "Alisertib")
# DMSOscores <- subset(ScoredG2_Mobj.integrated, ScoredG2_Mint$Arm == "vehicle")
# singG2_Mbox <- ggpubr::ggboxplot(ScoredG2_Mobj.integrated, x = "Arm", y = "TotalScore",
#                          color = "Arm", palette = c("aquamarine", "maroon"),
#                          order = c("vehicle", "alisertib"),
#                          ylab = "G2_M enrichment score", xlab = "Treatment")
# plot(singG2_Mbox + ggpubr::stat_compare_means(method = "t.test"))

##################################################################################################

########## Aggregate data into one dataframe...

ScoredNPC1int$Sig <- "NPC1"
ScoredNPC2int$Sig <- "NPC2"
ScoredOPCint$Sig <- "OPC"
ScoredACint$Sig <- "AC"
ScoredMES1int$Sig <- "MES1"
ScoredMES2int$Sig <- "MES2"

singscoreDat <- rbind(ScoredNPC1int,
                      ScoredNPC2int,
                      ScoredOPCint,
                      ScoredACint,
                      ScoredMES1int,
                      ScoredMES2int)

singscoreDat

armDF <- as.data.frame(obj.integrated$arm)
armDF

singscoreDat <- merge(x = singscoreDat, y = armDF, by.x = "row.names", by.y = "row.names")
rownames(singscoreDat) <- singscoreDat$Row.names
singscoreDat$Row.names <- NULL
colnames(singscoreDat) <- c("TotalScore", "TotalDispersion", "Sig", "Arm")
singscoresum <- summarySE(singscoreDat, measurevar = "TotalScore", groupvars = c("Arm", "Sig"))


# ggplot(singscoresum, aes(x = Sig, y=TotalScore, fill=Arm)) +
#   geom_bar(position = position_dodge(), stat = "identity") +
#   geom_errorbar(aes(ymin=TotalScore-se, ymax=TotalScore+se),
#                 width = 0.2,
#                 position = position_dodge(.9))

#Side by side of DMSO vs Alisertib mean single-cell enrichment for Neftel et al signatures
#pdf(file = "AlisertibAnalysisPipeline/05_output/Singscore_NeftelEtAlTranscriptionalStateShiftFig_CI.pdf")
pdf("Neftelsigfor DMSO vs Alisertlib.pdf")
ggplot(singscoresum, aes(x=Sig, y=TotalScore, fill=Arm)) +
  geom_bar(aes(fill = Arm), position = position_dodge(), stat = "identity",
           colour="black",
           size = 0.3) +
  geom_errorbar(aes(ymin=TotalScore-ci, ymax = TotalScore+ci),
                width = 0.2,
                position = position_dodge(.9)) +
  scale_fill_manual(values = c("gray80", "chartreuse"),
                    name = "Treatment",
                    breaks = c("vehicle", "alisertib"),
                    labels = c("vehicle", "Alisertib")) +
  xlab("Transcriptional State Signature") +
  ylab("Mean Signature Enrichment Score") +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.6), angle = 90),
        axis.title.x = element_text(size = rel(1.6), angle = 0))
dev.off()

# With standard error instead of confidence interval. 
#pdf(file = "AlisertibAnalysisPipeline/05_output/Singscore_NeftelEtAlTranscriptionalStateShiftFig_SE.pdf")
pdf("Neftelsigfor DMSO vs Alisertlib.pdf_withstandarderror")

ggplot(singscoresum, aes(x=Sig, y=TotalScore, fill=Arm)) +
  geom_bar(aes(fill = Arm), position = position_dodge(), stat = "identity",
           colour="black",
           size = 0.3) +
  geom_errorbar(aes(ymin=TotalScore-se, ymax = TotalScore+se),
                width = 0.2,
                position = position_dodge(.9)) +
  scale_fill_manual(values = c("gray80", "chartreuse"),
                    name = "Treatment",
                    breaks = c("vehicle", "alisertib"),
                    labels = c("vehicle", "UM002")) +
  xlab("Transcriptional State Signature") +
  ylab("Mean Signature Enrichment Score") +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.6), angle = 90),
        axis.title.x = element_text(size = rel(1.6), angle = 0))
dev.off()

# UMAP plots with neftel state enrichment projected - sorry I need to make this one page... 
#pdf(file = "AlisertibAnalysisPipeline/05_output/CellTagHu_NeftelSingScoreFeaturePlots.pdf")
pdf("CellTagHu_NeftelSingScoreFeaturePlots_dim.pdf")
obj.integrated <- AddMetaData(obj.integrated, metadata = ScoredNPC1int[,"TotalScore"], col.name = "NPC1singscore")
FeaturePlot(obj.integrated, features = c("NPC1singscore"), min.cutoff = 0, cols = viridis::viridis(256, option = "B"))
obj.integrated <- AddMetaData(obj.integrated, metadata = ScoredNPC2int[,"TotalScore"], col.name = "NPC2singscore")
FeaturePlot(obj.integrated, features = c("NPC2singscore"), min.cutoff = 0, cols = viridis::viridis(256, option = "B"))
obj.integrated <- AddMetaData(obj.integrated, metadata = ScoredOPCint[,"TotalScore"], col.name = "OPCsingscore")
FeaturePlot(obj.integrated, features = c("OPCsingscore"), min.cutoff = 0, cols = viridis::viridis(256, option = "B"))
obj.integrated <- AddMetaData(obj.integrated, metadata = ScoredACint[,"TotalScore"], col.name = "ACsingscore")
FeaturePlot(obj.integrated, features = c("ACsingscore"), min.cutoff = 0, cols = viridis::viridis(256, option = "B"))
obj.integrated <- AddMetaData(obj.integrated, metadata = ScoredMES1int[,"TotalScore"], col.name = "MES1singscore")
FeaturePlot(obj.integrated, features = c("MES1singscore"), min.cutoff = 0, cols = viridis::viridis(256, option = "B"))
obj.integrated <- AddMetaData(obj.integrated, metadata = ScoredMES2int[,"TotalScore"], col.name = "MES2singscore")
FeaturePlot(obj.integrated, features = c("MES2singscore"), min.cutoff = 0, cols = viridis::viridis(256, option = "B"))
DimPlot(obj.integrated, split.by = "arm")
dev.off()

#### Can we bin cells based on singscores? 

singscoredf <- data.frame(row.names = rownames(obj.integrated@meta.data),
                          NPC1 = obj.integrated@meta.data$NPC1singscore,
                          NPC2 = obj.integrated@meta.data$NPC2singscore,
                          OPC = obj.integrated@meta.data$OPCsingscore,
                          AC = obj.integrated@meta.data$ACsingscore,
                          MES1 = obj.integrated@meta.data$MES1singscore,
                          MES2 = obj.integrated@meta.data$MES2singscore
)


###

singscoredf$assignment <- "Ambiguous"
for (i in 1:length(rownames(singscoredf))){
  if (singscoredf$NPC1[i] > singscoredf$NPC2[i]){
    if (singscoredf$NPC1[i] > singscoredf$OPC[i]){
      if (singscoredf$NPC1[i] > singscoredf$AC[i]){
        if (singscoredf$NPC1[i] > singscoredf$MES1[i]){
          if (singscoredf$NPC1[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "NPC1"
          }
        }
      }
    }
  }
  if (singscoredf$NPC2[i] > singscoredf$NPC1[i]){
    if (singscoredf$NPC2[i] > singscoredf$OPC[i]){
      if (singscoredf$NPC2[i] > singscoredf$AC[i]){
        if (singscoredf$NPC2[i] > singscoredf$MES1[i]){
          if (singscoredf$NPC2[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "NPC2"
          }
        }
      }
    }
  }
  if (singscoredf$OPC[i] > singscoredf$NPC1[i]){
    if (singscoredf$OPC[i] > singscoredf$NPC2[i]){
      if (singscoredf$OPC[i] > singscoredf$AC[i]){
        if (singscoredf$OPC[i] > singscoredf$MES1[i]){
          if (singscoredf$OPC[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "OPC"
          }
        }
      }
    }
  }
  if (singscoredf$AC[i] > singscoredf$NPC1[i]){
    if (singscoredf$AC[i] > singscoredf$NPC2[i]){
      if (singscoredf$AC[i] > singscoredf$OPC[i]){
        if (singscoredf$AC[i] > singscoredf$MES1[i]){
          if (singscoredf$AC[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "AC"
          }
        }
      }
    }
  }
  if (singscoredf$MES1[i] > singscoredf$NPC1[i]){
    if (singscoredf$MES1[i] > singscoredf$NPC2[i]){
      if (singscoredf$MES1[i] > singscoredf$OPC[i]){
        if (singscoredf$MES1[i] > singscoredf$AC[i]){
          if (singscoredf$MES1[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "MES1"
          }
        }
      }
    }
  }
  if (singscoredf$MES2[i] > singscoredf$NPC1[i]){
    if (singscoredf$MES2[i] > singscoredf$NPC2[i]){
      if (singscoredf$MES2[i] > singscoredf$OPC[i]){
        if (singscoredf$MES2[i] > singscoredf$AC[i]){
          if (singscoredf$MES2[i] > singscoredf$MES1[i]){
            singscoredf$assignment[i] <- "MES2"
          }
        }
      }
    }
  }
}

singscoredf$cell <- rownames(singscoredf)
singscoredf

################################################################################
# Conver 6 states to 4 here... will also have to convert in another area?

for (i in 1:length(rownames(singscoredf))){
  print(singscoredf$assignment[i])
  if (singscoredf$assignment[i] %in% c("NPC1", "NPC2")){
    singscoredf$Neftel_State[i] <- "NPC"
  } else if (singscoredf$assignment[i] %in% c("MES1", "MES2")){
    singscoredf$Neftel_State[i] <- "MES"
  } else {
    singscoredf$Neftel_State[i] <- as.character(singscoredf$assignment[i])
  }
}

singscoredf
################################################################################
singscorebin <- singscoredf[,c("cell", "assignment")]
singscorebin$cell <- NULL

obj.integrated <- AddMetaData(obj.integrated, metadata = singscorebin, col.name = "singscorebin")

singscorebin2 <- singscoredf[,c("cell", "Neftel_State")]
singscorebin2$cell <- NULL

# neftelState <- singscorebin[c("Neftel_State")]
obj.integrated <- AddMetaData(obj.integrated, metadata = singscorebin2)

################################################################################

Idents(obj.integrated) <- obj.integrated$singscorebin 
obj.integrated <- subset(obj.integrated, idents = c("NPC1", "NPC2", "OPC", "AC", "MES1", "MES2"))
stateUMAP <- DimPlot(obj.integrated, group.by = "singscorebin")#)cells = WhichCells(obj.integrated, idents = c("NPC1", "NPC2", "OPC", "AC", "MES1", "MES2")
neftelstateUMAP <- DimPlot(obj.integrated, group.by = "Neftel_State", pt.size = 0.5, split.by = "arm", cols = c("forestgreen", "darkred", "navy", "violet")) + theme_void()
pdf("neftelstateumap")
neftelstateUMAP
dev.off()
################################################################################
# test dim redux
################################################################################
# Is there a way to better show the shift via dim redux? run umap on neftel features?
suvafeatures <- c(suvasigs$MES2, suvasigs$MES1, suvasigs$AC, suvasigs$OPC, suvasigs$NPC1, suvasigs$NPC2)
class(suvafeatures)
suvafeatures <- unique(suvafeatures[suvafeatures != ""])
suvafeatures

test <- obj.integrated
test <- RunPCA(test, features = suvafeatures)
ElbowPlot(test)
test <- RunUMAP(test, dims = 1:20, n.neighbors = 5)
DimPlot(test, group.by = "Neftel_State", split.by = "arm")
################################################################################
# Scrabble hierarchy plot...

# need a dataframe, rows are cells, columns are states and scores...

# need to average scores for combined states...
inputdf <- singscoredf[c("NPC1", "NPC2", "OPC", "AC", "MES1", "MES2")]
inputdf$MES <- rowMeans(inputdf[ , c("MES1", "MES2")], na.rm=TRUE)
inputdf$NPC <- rowMeans(inputdf[ , c("NPC1", "NPC2")], na.rm=TRUE)
inputdf <- inputdf[c("AC", "MES", "NPC", "OPC")]
colnames(inputdf) <- c("AC.like", "MES.like", "NPC.like", "OPC.like")
obj.integrated <- AddMetaData(obj.integrated, metadata = inputdf)
hierarchy <- scrabble::hierarchy(m = inputdf, quadrants = NULL, log.scale = T)
scrabble::plot_hierarchy(hierarchy)

# What if we just add the hierarchy output back to the seurat object as a reduction?
hier <- as.matrix(hierarchy)
colnames(hier) <- c("hierarchy_1", "hierarchy_2")
obj.integrated[["hierarchy"]] <- CreateDimReducObject(embeddings = hier, key = "hierarchy_", assay = "RNA")

pdf("hierarchyModuleDimPlot.pdf")
DimPlot(obj.integrated, reduction = "hierarchy", group.by = "Neftel_State") + 
  xlim(c(-0.4,0.4)) + ylim(c(-0.4,0.4)) + 
  xlab(label = "Relative Meta-Module Score [log2(| SC1 - SC2 | + 1)]") + 
  ylab(label = "Relative Meta-Module Score [log2(| SC1 - SC2 | + 1)]") + 
  theme_minimal()
dev.off()

plist <- FeaturePlot(obj.integrated, reduction = "hierarchy", 
            features = c("ACsingscore", 
                         "MES1singscore",
                         "MES2singscore", 
                         "NPC1singscore",
                         "NPC2singscore",
                         "OPCsingscore"), 
            coord.fixed = T, 
            combine = F, cols = viridis(256, option = "turbo"))
plist <- lapply(X = plist, FUN = function(p) p + xlim(c(-0.4,0.4)) + ylim(c(-0.4,0.4)) + 
                                             xlab(label = "Relative Meta-Module Score [log2(| SC1 - SC2 | + 1)]") + 
                                             ylab(label = "Relative Meta-Module Score [log2(| SC1 - SC2 | + 1)]") + 
                                             theme_minimal())
plist <- CombinePlots(plist)

pdf("hierarchyModuleScorePlots.pdf", width = 21, height = 14)
plist
dev.off()
################################################################################

pdf("CellTagHu_SubtypeBinnedUMAP.pdf")
stateUMAP
dev.off()

singscoreNeftelList <- c(
  "ACsingscore",
  "MES1singscore",
  "MES2singscore",
  "NPC1singscore",
  "NPC2singscore",
  "OPCsingscore"
)

pdf("singscoreCellTagSubtypeDotPlots.pdf",height= 3.5,width = 12)
DotPlot(
  obj.integrated,
  assay = NULL,
  group.by = "singscorebin",
  features = singscoreNeftelList,
  cols = c("black", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 9,
  scale.by = "radius"
)
dev.off()

NeftelList <- c(
  "AC.like", "MES.like", "NPC.like", "OPC.like"
)
neftelStateDotPlot <- DotPlot(obj.integrated,
                              assay = NULL,
                              group.by = "Neftel_State",
                              features = NeftelList,
                              cols = c("white", "red"),
                              col.min = -2.5,
                              col.max = 2.5,
                              dot.min = 0,
                              dot.scale = 14,
                              scale.by = "radius",
                              scale = F) + 
                              coord_flip() + 
  theme_minimal() + theme(text = element_text(size = 20))

pdf("neftelStateDotPlot.pdf", height = 3.5)
neftelStateDotPlot
dev.off()

# Idents(int) <- int$singscorebin
# StateMarkers <- FindAllMarkers(obj.integrated, test.use = "MAST", only.pos = TRUE)
# library(dplyr)
# 
# StateMarkers
# topStateMarkers <- StateMarkers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
# topStateMarkers
# stateHeatmap <- DoHeatmap(obj.integrated, features = topStateMarkers$gene, label = FALSE, assay = "RNA")
# pdf(file = "AlisertibAnalysisPipeline/05_output/alisertib_DMSO_stateHeatmap.pdf")
# stateHeatmap
# dev.off()

Idents(obj.integrated) <- obj.integrated$arm
bincounts <- table(obj.integrated$singscorebin, obj.integrated$arm)
DMSOsum <- sum(bincounts[,"vehicle"])
Alisertibsum <- sum(bincounts[,"alisertib"])
bincounts2 <- as.data.frame(bincounts)
colnames(bincounts2) <- c("State", "Treatment", "Freq")

for (i in 1:length(bincounts2$Freq)){
  if (bincounts2$Treatment[i] == "vehicle"){
    bincounts2$Percent[i] <- bincounts2$Freq[i]/DMSOsum*100
  }
  if (bincounts2$Treatment[i] == "alisertib"){
    bincounts2$Percent[i] <- bincounts2$Freq[i]/Alisertibsum*100
  }
}

library(reshape2)
bincounts3 <- recast(bincounts2, formula = State ~ Treatment, measure.var = "Percent")
rownames(bincounts3) <- bincounts3$State
bincounts3$State <- NULL
bincounts4 <- as.matrix(bincounts3)

#' Plots a series of barplots and connects them
#' Modified from https://stackoverflow.com/questions/22560850/barplot-with-connected-series
#'
#' @param dat NxM matrix with N rows as features and M columns as samples
#' @param color Vector of N colors
#' @param space Space between barplots
#' @param alpha Alpha for area connecting barplots
#'
#' @examples
#' dat <- matrix(rnorm(100),10,10)
#' dat <- abs(matrix(rnorm(100),10,10))
#' connectedBarplot(dat, color=rainbow(nrow(dat)))
#'
connectedBarplot <- function(dat, color=rainbow(nrow(dat)), space=1, alpha=0.5, ...) {
  b <- barplot(dat, col=color, space = space, ...)

  for (i in seq_len(ncol(dat) - 1)) {
    lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## bottom line

    for (j in seq_len(nrow(dat))) {
      if (j == 1) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(0, dat[j,i], dat[j,i+1], 0),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
      if (j == 2) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
      if (j > 2) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], colSums(dat[1:(j-1),])[i+1]),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
    }
  }
}

pdf("NeftelEtAlStateShiftBarplot_singscore_total.pdf")
connectedBarplot(dat = bincounts4)
dev.off()

connectedBarplot(dat = bincounts4)

#######################################################################################
#######################################################################################
#######################################################################################

### Cell Cyucle Shift Analysis

## Part I - Seurat method

## Part II - Singscore enrichment

#######################################################################################
#######################################################################################
#######################################################################################

# Part I - Seurat module scoring

DefaultAssay(obj.integrated) <- "RNA"
#obj.integrated <- CellCycleScoring(obj.integrated, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
obj.integrated$CellCycleIdents <- Idents(obj.integrated)

# RidgePlot(obj.integrated, features = c("G2M.Score", "S.Score"))

DimPlot(obj.integrated, group.by = "CellCycleIdents")

VlnPlot(obj.integrated, features = c("G2M.Score", "S.Score"), group.by = "arm", pt.size = 0)

CellCycleCounts <- table(obj.integrated$CellCycleIdents, obj.integrated$arm)

DMSOsum <- sum(CellCycleCounts[,"vehicle"])
Alisertibsum <- sum(CellCycleCounts[,"alisertib"])
CellCycleCounts2 <- as.data.frame(CellCycleCounts)
CellCycleCounts2
colnames(CellCycleCounts2) <- c("Phase", "Treatment", "Freq")

for (i in 1:length(CellCycleCounts2$Freq)){
  if (CellCycleCounts2$Treatment[i] == "vehicle"){
    CellCycleCounts2$Percent[i] <- CellCycleCounts2$Freq[i]/DMSOsum*100
  }
  if (CellCycleCounts2$Treatment[i] == "alisertib"){
    CellCycleCounts2$Percent[i] <- CellCycleCounts2$Freq[i]/Alisertibsum*100
  }
}

library(reshape2)
CellCycleCounts3 <- recast(CellCycleCounts2, formula = Phase ~ Treatment, measure.var = "Percent")
rownames(CellCycleCounts3) <- CellCycleCounts3$Phase
CellCycleCounts3$Phase <- NULL
CellCycleCounts4 <- as.matrix(CellCycleCounts3)

#' Plots a series of barplots and connects them
#' Modified from https://stackoverflow.com/questions/22560850/barplot-with-connected-series
#'
#' @param dat NxM matrix with N rows as features and M columns as samples
#' @param color Vector of N colors
#' @param space Space between barplots
#' @param alpha Alpha for area connecting barplots
#'
#' @examples
#' dat <- matrix(rnorm(100),10,10)
#' dat <- abs(matrix(rnorm(100),10,10))
#' connectedBarplot(dat, color=rainbow(nrow(dat)))
#'
connectedBarplot <- function(dat, color=rainbow(nrow(dat)), space=1, alpha=0.5, ...) {
  b <- barplot(dat, col=color, space = space, ...)

  for (i in seq_len(ncol(dat) - 1)) {
    lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## bottom line

    for (j in seq_len(nrow(dat))) {
      if (j == 1) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(0, dat[j,i], dat[j,i+1], 0),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
      if (j == 2) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
      if (j > 2) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], colSums(dat[1:(j-1),])[i+1]),
                col=adjustcolor(color[j], alpha.f=alpha))
      }
    }
  }
}

pdf("Alisertib_DMSO_CellCycleShift_connectedBarplot.pdf")
connectedBarplot(dat = CellCycleCounts4)
dev.off()


pdf("CellTag_Hu_CellCycleIdents.pdf")
DimPlot(obj.integrated, group.by = "CellCycleIdents")
dev.off()

# Now, is this significant? growth of G1, loss of G2M and S?

CellCycleScoreDF <- data.frame(G2M.Score = obj.integrated$G2M.Score, S.Score = obj.integrated$S.Score, Treatment = obj.integrated$arm)
mlt <- melt(CellCycleScoreDF)
colnames(mlt) <- c("Treatment", "Phase", "Score")

library(Rmisc)
CellCyclesum <- summarySE(mlt, measurevar = "Score", groupvars = c("Treatment", "Phase"))
CellCyclesum
pdf("CellTTagHu_UM002_CellCyclePhaseShiftPairedBarplot.pdf")
# 95% CI instead of SE...
ggplot(CellCyclesum, aes(x=Phase, y=Score, fill=Treatment)) +
  geom_bar(aes(fill = Treatment), position = position_dodge(), stat = "identity",
           colour="black",
           size = 0.3) +
  geom_errorbar(aes(ymin=Score-ci, ymax = Score+ci),
                width = 0.2,
                position = position_dodge(.9)) +
  # geom_signif(y_position = c(0.33), xmin=c(1.23, 1.77, 2.77, 3.77, 4.77, 5.77),
  #             xmax = c(.77, 2.23, 3.23, 4.23, 5.23, 6.23),
  #             annotations = c("NS", "NS", "NS", "NS", "*", "*"),
  #             tip_length = c(0.55, 0.55, 0.24, 0.3, 1.1, 1.15, 1.05, 1.05, 0.85, 0.73, 0.85, 1.1),
  #             textsize = 5) +
  scale_fill_manual(values = c("gray80", "chartreuse"),
                    name = "Treatment",
                    breaks = c("vehicle", "alisertib"),
                    labels = c("vehicle", "alisertib")) +
  xlab("Cell Cycle Phase Signature") +
  ylab("Module Score") +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.6), angle = 90),
        axis.title.x = element_text(size = rel(1.6), angle = 0))
dev.off()


### How to make normalized plot showing growth of cells in S phase...?

g2mmlt <- subset(mlt, mlt$Phase == "G2M.Score")
smlt <- subset(mlt, mlt$Phase == "S.Score")

g2mres <- t.test(Score ~ Treatment, data = g2mmlt, paired = FALSE)
g2mres

sres <- t.test(Score ~ Treatment, data = smlt, paired = FALSE)
sres

############################################################################
#obj.integrated <- readRDS("targetfeatures.rds")
# Differential expression between DMSO and Alisertib treated cells. 

Idents(obj.integrated) <- obj.integrated$arm
AlisertibMarkers <- FindAllMarkers(obj.integrated, test.use = "t", only.pos = TRUE)
AlisertibMarkers2 <- subset(AlisertibMarkers, AlisertibMarkers$avg_log2FC != "Inf")
AlisertibMarkers2 <- subset(AlisertibMarkers2, AlisertibMarkers$avg_log2FC != "-Inf")
topnmarkers <- AlisertibMarkers2 %>% group_by(cluster) %>% top_n(25, wt = avg_log2FC)
# topnmarkers <- su

topnmarkers$gene
top25heatmap <- DoHeatmap(obj.integrated, features = c(topnmarkers$gene), slot = "scale.data")
#top25heatmap

pdf("top25heatmap.pdf")
top25heatmap
dev.off()

### Make normalized state shift barplot...

singscoreDat <- rbind(ScoredNPC1obj.integrated,
                      ScoredNPC2obj.integrated,
                      ScoredOPCobj.integrated,
                      ScoredACobj.integrated,
                      ScoredMES1obj.integrated,
                      ScoredMES2obj.integrated)

singscoreDat


# Make column of arm from barcodes...

# for (i in 1:length(rownames(singscoreDat))){
#   singscoreDat$Arm[i] <- unlist(strsplit(as.character(rownames(singscoreDat)[i]), split = "_", fixed = TRUE))[1]
# }


singscoreAli <- subset(singscoreDat, singscoreDat$Arm == "alisertib")
singscoreAli
singscoresumDMSO <- subset(singscoresum, singscoresum$Arm == "vehicle")
rownames(singscoresumDMSO) <- singscoresumDMSO$Sig
singscoresumDMSO["AC","TotalScore"]

for (i in 1:length(rownames(singscoreAli))){
  j <- singscoreAli$Sig[i]
  print(j)
  k <- singscoreAli$TotalScore[i]
  print(k)
  singscoreAli$deltaMean[i] <- k - singscoresumDMSO[j, "TotalScore"]
}
singscoreAli

singscoreAlisum <- summarySE(singscoreAli, measurevar = "deltaMean", groupvars = c("Arm", "Sig"))
singscoreAlisum
################################################################################

# Normalized Module Score Plot...
normAliShiftPlot1 <- ggplot(singscoreAlisum, aes(x=Sig, y=deltaMean, fill=Sig)) +
                    geom_bar(aes(fill = Sig), position = position_dodge(), stat = "identity",
                             colour="black",
                             size = 0.3) +
                    geom_errorbar(aes(ymin=deltaMean-ci, ymax = deltaMean+ci),
                                  width = 0.2,
                                  position = position_dodge(.9)) +
                    xlab("Transcriptional State Signature") +
                    ylab("Delta Mean Signature Enrichment Score") +
                    theme_bw() +
                    theme(axis.title.y = element_text(size = rel(1.6), angle = 90),
                          axis.title.x = element_text(size = rel(1.6), angle = 0))
normAliShiftPlot1

################################################################################
# 4 States

nefteldf <- data.frame(cell_barcode = rownames(obj.integrated@meta.data), 
           AC.like = obj.integrated$AC.like,
           MES.like = obj.integrated$MES.like,
           NPC.like = obj.integrated$NPC.like,
           OPC.like = obj.integrated$OPC.like,
           Arm = obj.integrated$arm)

nefteldf <- melt(nefteldf, id.vars = c("cell_barcode", "Arm"))
aliNeftelDF <- subset(nefteldf, nefteldf$Arm == "alisertib")

colnames(aliNeftelDF)
neftelstatesum <- summarySE(nefteldf, measurevar = "value", 
                            groupvars = c("Arm", "variable"))
neftelstatesumDMSO <- subset(neftelstatesum, neftelstatesum$Arm == "alisertib")
rownames(neftelstatesumDMSO) <- neftelstatesumDMSO$variable
neftelstatesumDMSO["AC","value"]

for (i in 1:length(rownames(aliNeftelDF))){
  j <- aliNeftelDF$variable[i]
  print(j)
  k <- aliNeftelDF$value[i]
  print(k)
  aliNeftelDF$deltaMean[i] <- k - neftelstatesumDMSO[j, "value"]
}

aliNeftelDF

neftelAlisum <- summarySE(aliNeftelDF, measurevar = "deltaMean", groupvars = c("Arm", "variable"))
neftelAlisum

write.csv(neftelAlisum, file = "inVivoNeftelShiftSummary2.csv")

normvehShiftPlot_fourStates <- ggplot(neftelAlisum, aes(x=variable, y=deltaMean, fill=variable)) +
  geom_bar(aes(fill = variable), position = position_dodge(), stat = "identity",
           colour="black",
           size = 0.3) +
  geom_errorbar(aes(ymin=deltaMean-ci, ymax = deltaMean+ci),
                width = 0.2,
                position = position_dodge(.9)) +
  xlab("Transcriptional State Signature") +
  ylab("Normalized Enrichment Score") +
  theme_minimal() + NoLegend() +
  theme(axis.title.y = element_text(size = rel(1.0), angle = 90),
        axis.title.x = element_text(size = rel(1.0), angle = 0),
        text = element_text(size = 30))

pdf("normvehShiftPlot_4States1.pdf",width=10)
normvehShiftPlot_fourStates
dev.off()

################################################################################

# Turn this into a raincloud plot...

# plotdat <- int@meta.data
# plotdat1 <- melt(plotdat, id.vars = c("Arm"),
#                  measure.vars = c("oxphog_Signature"),
#                  variable.name = "oxphog_Signature",
#                  value.name = "enrichmentScore")

# Calculate Summary Statistics

lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

singscoreAli$logDeltaMean <- log1p(singscoreAli$deltaMean)
colnames(singscoreAli)
singscoreAli
sumId <- ddply(singscoreAli, ~logDeltaMean, summarise,
               mean = mean(logDeltaMean), median = median(logDeltaMean),
               lower = lb(logDeltaMean), upper = ub(logDeltaMean))

g <- ggplot(data = singscoreAli, aes(y = logDeltaMean, x = Sig, fill = Sig)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = logDeltaMean, color = Sig), position = position_jitter(width = .15), size = .5, alpha = 0.2) +
  geom_boxplot(width = .3, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_discrete() +
  scale_fill_discrete() +
  # coord_flip() +
  theme_bw() +
  raincloud_theme

pdf("StateBin_Raincloud.pdf")
g
dev.off()

################################################################################
# Test for normality of normalized values...
################################################################################

ggdensity(singscoreAli$deltaMean,
          main = "Alisertib vs. DMSO",
          xlab = "Normalized Enrichment")

sample <- sample_n(singscoreAli, 30)
shapiro.test(sample$deltaMean)

ggqqplot(singscoreAli, x = "deltaMean")
ggqqplot(sample, x = "deltaMean")


pdf("stateShift.pdf",width = 14)
plot_grid(plotlist = list(stateUMAP, normAliShiftPlot), ncol = 2, labels = "AUTO")
dev.off()

## 20210606 - Characterize the mesenchymal, alisertib resistant population. 

# Split into DMSO and Alisertib. 

# FindMarkers between Mes1, Mes2 and other cells within single treatment arm

# run through DAVID 

Idents(obj.integrated)
DMSO <- subset(obj.integrated, idents = c("vehicle"))
Idents(DMSO) <- DMSO$singscorebin

Idents(Alisertib) <- Alisertib$singscorebin
# aliMesMarkers <- FindMarkers(Alisertib, ident.1 = c("MES1", "MES2"), test.use = "MAST", logfc.threshold = 0)
# 
# # Enhanced volcano of results!...
# library(EnhancedVolcano)
# 
# head(colnames(aliMesMarkers))
# pdf(file = "AlisertibAnalysisPipeline/05_output/aliTreatedMesVolcano.pdf")
# EnhancedVolcano(aliMesMarkers, lab = rownames(aliMesMarkers), x = "avg_log2FC", y = "p_val", pCutoff = 10e-10, FCcutoff = 0.66)
# dev.off()

# saveRDS(obj.integrated, file = "AlisertibAnalysisPipeline/05_output/05out_int.RDS")

## ggalluvial to show shift in assigned identity?

groupingCounts <- table(obj.integrated@meta.data[,"singscorebin"],
                        obj.integrated@meta.data[,"arm"]
                        
                        )
groupingCounts

dmsoSum <- sum(groupingCounts[,"vehicle"])
aliSum <- sum(groupingCounts[,"alisertib"])
groupingCounts2 <- as.data.frame(groupingCountsSplit)

groupingCounts2
gcSplit<- as.data.frame(table(obj@meta.data[,"singscorebin"],
                                             obj.integrated@meta.data[,"CellCycleIdents"],
                                             obj.integrated@meta.data[,"Replicate"]))
colnames(gcSplit) <- c("State", "Idents", "patient_ID", "Freq") 

colnames(groupingCounts2) <- c("State", "arm", "Freq")
for (i in 1:length(groupingCounts2$Freq)){
  if (groupingCounts2$arm[i] == "vehicle"){
    groupingCounts2$Percent[i] <- groupingCounts2$Freq[i]/dmsoSum*100
  }
  if (groupingCounts2$arm[i] == "alisertib"){
    groupingCounts2$Percent[i] <- groupingCounts2$Freq[i]/aliSum*100
  }
}
groupingCounts2$State <- as.factor(groupingCounts2$State)
groupingCounts2$arm <- factor(groupingCounts2$arm, levels = c("vehicle", "alisertib"))

library(ggalluvial)
# is_alluvia_form(groupingCounts2)

for (i in 1:length(rownames(groupingCounts2))){
  print(groupingCounts2$State[i])
  if (groupingCounts2$State[i] %in% c("NPC1", "NPC2")){
    groupingCounts2$State2[i] <- "NPC"
  } else if (groupingCounts2$State[i] %in% c("MES1", "MES2")){
    groupingCounts2$State2[i] <- "MES"
  } else {
    groupingCounts2$State2[i] <- as.character(groupingCounts2$State[i])
    }
  }

groupingCounts2$State2 <- as.factor(groupingCounts2$State2)
class(groupingCounts2$State)
class(groupingCounts2$State2)

groupingCounts2

groupingCounts3 <- groupingCounts2
groupingCounts3$Arm_State2 <- paste0(groupingCounts3$Arm, "_", groupingCounts3$State2)
groupingCounts3$Arm <- NULL
groupingCounts3$State <- NULL
groupingCounts3$State2 <- NULL
#edited this code due to sum function as it's not accepting factor for calculation
'''groupingCounts3
groupingCounts3 <- as.data.frame(
  groupingCounts3 %>% 
    group_by(Arm_State2) %>% 
    mutate(across(.fns = sum)) %>% 
    distinct()
)'''

groupingCounts3 <- as.data.frame(
  groupingCounts3 %>% 
    group_by(Arm_State2) %>% 
    mutate(across(.fns = ~if(is.numeric(.)) sum(.) else .)) %>% 
    distinct()
)



for (i in 1:length(rownames(groupingCounts3))){
  groupingCounts3$Arm[i] <- unlist(strsplit(groupingCounts3$Arm_State2[i], split = "_", fixed = T))[1]
  groupingCounts3$State2[i] <- unlist(strsplit(groupingCounts3$Arm_State2[i], split = "_", fixed = T))[2]
}

groupingCounts3

alluv <- ggplot(groupingCounts2,
       aes(x = Arm, stratum = State, alluvium = State,
           y = Percent,
           fill = State, label = State)) +
  scale_fill_brewer(type = "qual", palette = "Spectral") +
  # scale_fill_discrete(guide = guide_legend(reverse = TRUE)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") + theme_minimal() +
  ggtitle("Transcriptional State Shift")

alluv2 <- ggplot(groupingCounts3,
       aes(x = Arm, stratum = State2, alluvium = State2,
           y = Percent,
           fill = State2, label = State2)) +
  scale_fill_brewer(type = "qual", palette = "Spectral") +
  # scale_fill_discrete(guide = guide_legend(reverse = TRUE)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum(width = as.numeric(1/3)) +
  theme(legend.position = "bottom") + theme_minimal() +
  ggtitle("Transcriptional State Shift")

pdf("FourState_alisertibInVivo_alluvial.pdf")
alluv2
dev.off()

groupingCounts3

df2 <- data.frame()
for (i in 1:length(unique(groupingCounts3$State2))){
  print(i)
  groupBy_search <- unique(groupingCounts3$State2)[i]
  plotdf_state <- subset(groupingCounts3, groupingCounts3$State2 == groupBy_search)
  print(plotdf_state)
  rownames(plotdf_state) <- plotdf_state$arm
  alipct <- plotdf_state["alisertib", "Percent"]
  print(alipct)
  dmsopct <- plotdf_state["vehicle", "Percent"]
  print(dmsopct)
  df1 <- data.frame(State = plotdf_state$State2[1],
                    deltapct <- alipct - dmsopct)
  df1
  df2 <- rbind(df2, df1)
}

colnames(df2) <- c("State", "deltapct")
library(scales)
percShiftBP <- ggbarplot(df2, x = "State", y = "deltapct", title = "State shift", 
                         ylab = "Change in % Proportion", xlab = "Transcriptional State", color = "State", fill = "State") +
                         # add = c("mean_se", "jitter"), position = position_dodge(.8)) + 
  scale_color_manual(values = c(OPC = "black", NPC = muted("black"), AC = "black", MES = "black")) + theme_minimal() + NoLegend()

pdf(file = "alisertib_vs_dmso_NeftelState_percent_bp.pdf")
percShiftBP
dev.off()

write.csv(groupingCounts2, "inVivoAli_freqPerc_shiftAlluvialData_2.csv")
write.csv(groupingCounts2, file = "~/Documents/GU_Bioinformatics/AyadLab/ISOSCELES/inVivoAli_freqPerc_shiftAlluvialData.csv")



#loading these takes less time in computation (future reference purpose) 

#obj.integrated <- readRDS("out_int.RDS")
#singscoreDat <- readRDS("singscoreDat.rds")
#singscorebin <- readRDS("singscorebin.rds")
#singscoreAli <- readRDS("singscoreAli.rds")

numeric_suff <- as.integer(gsub(".*_(\\d+)$", "\\1", colnames(obj.integrated)))
obj.integrated$Replicate <- numeric_suff
obj.integrated@meta.data$Replicate
unique(obj.integrated$Replicate)

saveRDS(obj.integrated,"integrated_new1.rds")

#############

nefteldf <- data.frame(cell_barcode = rownames(obj.integrated@meta.data), 
                       AC.like = obj.integrated$AC.like,
                       MES.like = obj.integrated$MES.like,
                       NPC.like = obj.integrated$NPC.like,
                       OPC.like = obj.integrated$OPC.like,
                       Arm = obj.integrated$arm)

nefteldf <- melt(nefteldf, id.vars = c("cell_barcode", "Arm"))
aliNeftelDF <- subset(nefteldf, nefteldf$Arm == "alisertib")
colnames(aliNeftelDF)
neftelstatesum <- summarySE(nefteldf, measurevar = "value", 
                            groupvars = c("Arm", "variable"))
neftelstatesumDMSO <- subset(neftelstatesum, neftelstatesum$Arm == "vehicle")
rownames(neftelstatesumDMSO) <- neftelstatesumDMSO$variable
neftelstatesumDMSO["AC","value"]

for (i in 1:length(rownames(aliNeftelDF))){
  j <- aliNeftelDF$variable[i]
  print(j)
  k <- aliNeftelDF$value[i]
  print(k)
  aliNeftelDF$deltaMean[i] <- k - neftelstatesumDMSO[j, "value"]
}
aliNeftelDF

neftelAlisum <- summarySE(aliNeftelDF, measurevar = "deltaMean", groupvars = c("Arm", "variable"))
neftelAlisum
write.csv(neftelAlisum, file = "~/Documents/GU_Bioinformatics/AyadLab/ISOSCELES/inVivoNeftelShiftSummary.csv")

normAliShiftPlot_fourStates <- ggplot(neftelAlisum, aes(x=variable, y=deltaMean, fill=variable)) +
  geom_bar(aes(fill = variable), position = position_dodge(), stat = "identity",
           colour="black",
           size = 0.3) +
  geom_errorbar(aes(ymin=deltaMean-ci, ymax = deltaMean+ci),
                width = 0.2,
                position = position_dodge(.9)) +
  xlab("Transcriptional State Signature") +
  ylab("Normalized Enrichment Score") +
  theme_minimal() + NoLegend() +
  theme(axis.title.y = element_text(size = rel(1.0), angle = 90),
        axis.title.x = element_text(size = rel(1.0), angle = 0),
        text = element_text(size = 30))

pdf("AlisertibAnalysisPipeline_4States.pdf")
normAliShiftPlot_fourStates
dev.off()


#####################
library(dittoSeq)
frequencies <- table(obj.integrated@meta.data$Replicate, obj.integrated$arm)
frequencies_df <- as.data.frame(frequencies)
colnames(frequencies_df) <- c("Replicate", "Ident", "Count")
frequencies_df <- transform(frequencies_df,
                                     Percentage = Count / ave(Count, Replicate, FUN = sum) * 100)

dittobarplot(frequencies_df, x = "Ident", y = "Percentage", group = "Replicate",
             main = "%Comparison", xlab = "Identity")

dittoPlot(obj.integrated,"Replicate",group.by = "arm")
dittoPlot(obj.integrated,"AC.like",group.by = "Replicate")

dpb <- dittoBarPlot(obj.integrated, "AC.like", group.by = "Replicate")

obj.integrated$AC.like <- as.factor(obj.integrated$AC.like)
dittoDimPlot(obj.integrated, "AC.like")
dpg1 <- dittoPlot(obj.integrated, "AC.like", group.by = "Replicate")
dpg



pdf("dpg.pdf")
dpg1
dev.off()









# Subset the data for cells with the identity "alisertlib"
Alisertib<-  obj.integrated[obj.integrated$idents == "alisertlib", ]

# Create a table with the counts of each category within each replicate
table_data <- table(Alisertib$Replicate, Alisertib$AC.like,Alisertib$OPC.like,Alisertib$NPC.like,Alisertib$MES.like)

subset_data <- obj.integrated[obj.integrated$arm == "alisertlib", ]




dittodata <- as.data.frame(obj.integrated)

# Create a dittobarplot to visualize the percentages for each state across replicates
dittoBarPlot(
  data = data,
  x = c("AC.like", "OPC.like", "NPC.like", "MES.like"),
  group.by = "Replicate",
  main = "Percentage of Cell Cycle States by Replicate (alisertlib)"
)

data <- obj.integrated@meta.data

# Create the bar plot
barplot_state <- data %>%
  group_by(Replicate, Neftel_State) %>%
  summarise(count = n()) %>%
  group_by(Replicate) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ggplot(aes(x = Replicate, y = percent, fill = Neftel_State)) +
  geom_bar(stat = "identity") +
  labs(title = "Percentage of Each State Across Replicates",
       x = "Replicate",
       y = "Percentage") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_brewer(palette = "Set3")
pdf("barplot_state.pdf")
barplot_state
dev.off()



barplot_state_percent <- data %>%
  mutate(Replicate = factor(Replicate, levels = unique(Replicate))) %>%
  group_by(Replicate, Neftel_State) %>%
  summarise(count = n()) %>%
  group_by(Replicate) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ggplot(aes(x = Replicate, y = percent, fill = Neftel_State)) +
  geom_bar(stat = "identity") +
  labs(title = "Percentage of Each State Across Replicates",
       x = "Replicate",
       y = "Percentage") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_brewer(palette = "Set3")
pdf("barplot_state_percent.pdf")
barplot_state_percent
dev.off()

tail(obj.integrated@meta.data)
cellcycle_replicate_data1 <- obj.integrated@meta.data[, c("arm", "Replicate")]
