# devtools::install_github(repo = "satijalab/seurat", ref = "develop")

library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggplotify)
library(viridis)
library(Hmisc)
#library(vtable)

options(future.globals.maxSize = 4000 * 1024^2)

################################################################################

# Load 10X data & Create Seurat Object...

# samples <- dir(path = "cellranger_outs/")
# 
# obj.list <- list()
# for (s in samples) {
#   print(s)
#   obj.dat <- Read10X(data.dir = paste0("cellranger_outs/", s, "/outs/grch38_out/filtered_feature_bc_matrix/"))
#   obj <- CreateSeuratObject(counts = obj.dat, project = s)
#   obj.list <- append(obj.list, obj)
#   rm(list = c("obj", "obj.dat"))
# }
# 
# obj.list[[1]]@meta.data
# 
# ###
# 
# # And uniformly filter all datasets...
# for (i in 1:length(obj.list)) {
#   # obj.list[[i]]$percent.mt <- PercentageFeatureSet(obj.list[[i]], pattern = "^MT-")
#   obj.list[[i]][["percent.hg19"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^GRCh38-")
#   obj.list[[i]][["percent.mm10"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^mm10---")
# }
# 
# # for (i in 1:length(obj.list)) {
# #   obj.list[[i]] <- subset(x = obj.list[[i]],
# #                                  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
# # }
# 
# obj.list <- lapply(X = obj.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = FALSE)
#   x <- FindVariableFeatures(x, verbose = FALSE)
# })
# 
# features <- SelectIntegrationFeatures(object.list = obj.list)
# obj.list <- lapply(X = obj.list, FUN = function(x) {
#   x <- ScaleData(x, features = features, verbose = FALSE)
#   x <- RunPCA(x, features = features, verbose = FALSE)
# })
# 
# anchors <- FindIntegrationAnchors(object.list = obj.list, reference = c(1, 2), reduction = "rpca",
#                                   dims = 1:50)
# obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
# rm(anchors)
# 
# obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
# obj.integrated <- RunPCA(obj.integrated, verbose = FALSE)
# obj.integrated <- FindNeighbors(obj.integrated, dims = 1:50)
# obj.integrated <- FindClusters(obj.integrated)
# 
# obj.integrated <- RunUMAP(obj.integrated, dims = 1:50)
# 
# # saveRDS(obj.integrated, file="ali_r2_barnyard_obj.int.RDS")
obj.integrated <- readRDS("ali_r2_barnyard_obj.int.RDS")

ident_fp <- DimPlot(obj.integrated, group.by = "orig.ident") + coord_fixed()
clust_fp <- DimPlot(obj.integrated, group.by = "seurat_clusters", label = T) + coord_fixed()
hg19_fp <- FeaturePlot(obj.integrated, features = c("percent.hg19"), order = F) +
  scale_colour_viridis(option = "D") + coord_fixed()
mm10_fp <- FeaturePlot(obj.integrated, features = c("percent.mm10"), order = T) +
  scale_colour_viridis(option = "D") + coord_fixed()
alignment_fp <- plot_grid(plotlist = list(ident_fp, clust_fp, hg19_fp, mm10_fp),
                          ncol = 2, labels = "AUTO")

pdf("ali_r2_barnyard_alignment_featureplots.pdf", height = 14, width = 14)
alignment_fp
dev.off()

###

Idents(obj.integrated) <- obj.integrated$seurat_clusters
human_clusters <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "10")
mouse_clusters <- as.character(unique(obj.integrated$seurat_clusters)[unique(obj.integrated$seurat_clusters) %nin% human_clusters])

human_sub <- subset(obj.integrated, idents = human_clusters)
mouse_sub <- subset(obj.integrated, idents = mouse_clusters)

hist(human_sub$percent.hg19)
hist(human_sub$percent.mm10)

hist(mouse_sub$percent.hg19)
hist(mouse_sub$percent.mm10)

human_barcodes <- as.data.frame(colnames(human_sub))
colnames(human_barcodes) <- c("barcodes")
saveRDS(human_barcodes, file = "alisertib_r2_human_barcode_df.RDS")
human_barcodes <- readRDS(file = "alisertib_r2_human_barcode_df.RDS")

################################################################################

# Load 10X data & Create Seurat Object... from GRCh38 alignment

samples <- dir(path = "cellranger_outs/")

obj.list <- list()
for (s in samples) {
  print(s)
  obj.dat <- Read10X(data.dir = paste0("cellranger_outs/", s, "/outs/barnyard_out/filtered_feature_bc_matrix/"))
  obj <- CreateSeuratObject(counts = obj.dat, project = s)
  obj.list <- append(obj.list, obj)
  rm(list = c("obj", "obj.dat"))
}

# remove mouse cells via barcodes prior to GRCh38 QC...

for (i in 1:length(obj.list)) {
  print("Pre-barcode filter")
  print(dim(obj.list[[i]]))
  obj.list[[i]] <- RenameCells(object = obj.list[[i]], new.names = paste0(colnames(obj.list[[i]]), "_", i))
  obj.list[[i]] <- subset(obj.list[[i]], cells = c(human_barcodes$barcodes))
  print("Post-barcode filter")
  print(dim(obj.list[[i]]))
}

# And uniformly filter all datasets...

qc_plot_list <- list()
for (i in 1:length(obj.list)) {
  obj.list[[i]]$percent.mt <- PercentageFeatureSet(obj.list[[i]], pattern = "^MT-")
  # obj.list[[i]][["percent.hg19"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^GRCh38-")
  # obj.list[[i]][["percent.mm10"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^mm10---")
  count_feature_fp <- FeatureScatter(obj.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    ggtitle(label = unique(obj.list[[i]]$orig.ident[1]))
  mt_feature_fp <- FeatureScatter(obj.list[[i]], feature1 = "nFeature_RNA", feature2 = "percent.mt") + 
    ggtitle(label = unique(obj.list[[i]]$orig.ident[1]))  
  
  qc_fp <- plot_grid(plotlist = list(count_feature_fp, mt_feature_fp), ncol = 2, labels = "AUTO")
  qc_plot_list[[i]] <- qc_fp
}

pdf(file = "grch38_qc_prefilter.pdf", height = 8*7, width = 14)
plot_grid(plotlist = qc_plot_list, ncol = 1)
dev.off()
##need to change here
for (i in 1:length(obj.list)) {
  cutoff <- quantile(obj.list[[i]]$nCount_RNA)[4]
  obj.list[[i]] <- subset(obj.list[[i]], subset = nCount_RNA > 200 & nCount_RNA < cutoff & percent.mt < 5)
}


##############Feature based scaling(independant calculation)##################################################################################################################################################################################
#nFeature_RNA_min <- 200 + feature_sd_cutoff * sd(nFeature_RNA)
#feature_sd <- apply(obj.list[[i]]$RNA, 2, sd)
#feature_sd_cutoff <- 0.5  


for (i in 1:length(obj.list)){
  #cutoff <- quantile(obj.list[[i]]$nFeature_RNA)[4]
  obj.list[[i]] <- subset(obj.list[[i]], 
                          subset = nFeature_RNA > 200 & percent.mt < 5)
}


for (i in 1:length(obj.list)) {
  avg_nCount_RNA <- mean(obj.list[[i]]$)
  print(avg_nCount_RNA)
  obj.list[[i]]$nFeature_RNA <- obj.list[[i]]$nFeature_RNA[obj.list[[i]]$nFeature_RNA < avg_nCount_RNA & obj.list[[i]]$percent.mt < 5]
}





######################################################################################################################################################################################################################################################################################################################################################################################
qc_plot_list <- list()
for (i in 1:length(obj.list)) {
  # obj.list[[i]]$percent.mt <- PercentageFeatureSet(obj.list[[i]], pattern = "^MT-")
  # obj.list[[i]][["percent.hg19"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^GRCh38-")
  # obj.list[[i]][["percent.mm10"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^mm10---")
  count_feature_fp <- FeatureScatter(obj.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    ggtitle(label = unique(obj.list[[i]]$orig.ident[1]))
  mt_feature_fp <- FeatureScatter(obj.list[[i]], feature1 = "nFeature_RNA", feature2 = "percent.mt") + 
    ggtitle(label = unique(obj.list[[i]]$orig.ident[1]))  
  
  qc_fp <- plot_grid(plotlist = list(count_feature_fp, mt_feature_fp), ncol = 2, labels = "AUTO")
  qc_plot_list[[i]] <- qc_fp
}

pdf( "grch38_qc_postfilter.pdf", height = 8*7, width = 14)
plot_grid(plotlist = qc_plot_list, ncol = 1)
dev.off()

################################################################################

# integrate post-qc filtering

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reference = c(1, 2), reduction = "rpca",
                                  dims = 1:50)

# length(unlist(lapply(obj.list, FUN = "rownames")))
to_integrate <- unique(unlist(lapply(obj.list, FUN = "rownames")))

obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50, features.to.integrate = to_integrate)
rm(obj.list)
rm(anchors)

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE, features = rownames(obj.integrated))
obj.integrated <- RunPCA(obj.integrated, verbose = FALSE)
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:50)
obj.integrated <- FindClusters(obj.integrated)

obj.integrated <- RunUMAP(obj.integrated, dims = 1:50)

FeaturePlot(obj.integrated, features = c("AURKA", "CD44"), min.cutoff = 0, slot = "scale.data", cols = c("grey", "red"), order = T)

metadf <- as.data.frame(obj.integrated$orig.ident)
for (i in 1:length(rownames(metadf))){
  metadf$arm[i] <- unlist(strsplit(metadf$`obj.integrated$orig.ident`[i], split = "_", fixed = T))[2]
}
unique(metadf$arm)
metadf$`obj.integrated$orig.ident` <- NULL
metadf
obj.integrated <- AddMetaData(obj.integrated, metadata = metadf)

DimPlot(obj.integrated, group.by = "arm")
table(obj.integrated$seurat_clusters, obj.integrated$arm)
table(obj.integrated$arm)

VlnPlot(obj.integrated, features = c("CDK4", "APOE", "TOP2A", "CD44"), group.by = "arm", pt.size = 0)

VlnPlot(obj.integrated, features = c("percent.mt", "nCount_RNA", "nFeature_RNA"), group.by = "arm", pt.size = 0)

saveRDS(obj.integrated, file="ali_obj.int1.RDS")
obj.integrated <- readRDS(file="ali_r2_GRCh38_obj.int.RDS")


################################################################################
# mem too big - using rPCA instead - not really important here yet...

# # # Run scTransform on all objects...
# # 
# # for (i in 1:length(obj.list)) {
# #   obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = TRUE)
# # }
# # 
# # # Next, select features for downstream integration, and run PrepSCTIntegration,
# # # which ensures that all necessary Pearson residuals have been calculated.
# # 
# # obj.features <- SelectIntegrationFeatures(
# #   object.list = obj.list,
# #   nfeatures = 3000
# # )
# # 
# # obj.list <- PrepSCTIntegration(
# #   object.list = obj.list,
# #   anchor.features = obj.features,
# #   verbose = TRUE)
# # 
# # # Next, identify anchors and integrate the datasets. Commands are identical to
# # # the standard workflow, but make sure to set normalization.method = "SCT":
# # 
# # obj.anchors <- FindIntegrationAnchors(
# #   object.list = obj.list,
# #   normalization.method = "SCT",
# #   anchor.features = obj.features,
# #   verbose = TRUE
# # )
# # 
# # to_integrate <- Reduce(intersect, lapply(obj.anchors@object.list, rownames))
# # 
# # print("Preview of genes to integrate... ")
# # head(to_integrate)
# # 
# # obj <- IntegrateData(
# #   anchorset = obj.anchors,
# #   #  features.to.integrate = to_integrate,
# #   normalization.method = "SCT",
# #   verbose = TRUE
# # )
# # 
# # rm(obj.anchors)
# # rm(obj.features)
# # rm(obj.list)
# # rm(to_integrate)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # setwd(dir = "/Volumes/Robs SSD2/AlisertibExp_CellRangerOuts/")
# # DMSO_tween.data <- Read10X(data.dir = "DMSO_Tween_barnyard/outs/filtered_feature_bc_matrix/")
# # Ali_RdI.data <- Read10X(data.dir = "Alisertib_1_barnyard/outs/filtered_feature_bc_matrix/")
# # Ali_RdII.data <- Read10X(data.dir = "Alisertib_2_barnyard/outs/filtered_feature_bc_matrix/")
# # setwd(dir = "/Volumes/Robs SSD2/AlisertibExpt/")
# ################################################################################
# 
# DMSO_tween$percent.mt <- PercentageFeatureSet(DMSO_tween, pattern = "^MT-")
# DMSO_tween[["percent.hg19"]] <- PercentageFeatureSet(DMSO_tween, pattern = "^hg19-")
# DMSO_tween[["percent.mm10"]] <- PercentageFeatureSet(DMSO_tween, pattern = "^mm10-")
# 
# pdf(
#   file = "01_output/DMSO_tween_prefilterQC.pdf"
# )
# VlnPlot(DMSO_tween, features = c("percent.hg19", "percent.mm10"))
# FeatureScatter(object = DMSO_tween, feature1 = "nCount_RNA", feature2 = "percent.mt")
# FeatureScatter(object = DMSO_tween, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# dev.off()
# 
# # RdII
# 
# # RdII_DMSO$percent.mt <- PercentageFeatureSet(RdII_DMSO, pattern = "^MT-")
# # RdII_DMSO[["percent.hg19"]] <- PercentageFeatureSet(RdII_DMSO, pattern = "^hg19-")
# # RdII_DMSO[["percent.mm10"]] <- PercentageFeatureSet(RdII_DMSO, pattern = "^mm10-")
# 
# # pdf(
# #   file = "01_output/RdII_DMSO_prefilterQC.pdf"
# # )
# # VlnPlot(RdII_DMSO, features = c("percent.hg19", "percent.mm10"))
# # FeatureScatter(object = RdII_DMSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # FeatureScatter(object = RdII_DMSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # dev.off()
# 
# ################################################################################
# 
# # NM002
# #######
# 
# Ali_I <- CreateSeuratObject(
#   counts = Ali_RdI.data,
#   min.cells = 3,
#   min.features = 200,
#   project = "Ali_I"
# )
# 
# Ali_II <- CreateSeuratObject(
#   counts = Ali_RdII.data,
#   min.cells = 3,
#   min.features = 200,
#   project = "Ali_II"
# )
# 
# Ali_I$Arm <- "Alisertib"
# Ali_I <- RenameCells(
#   object = Ali_I,
#   new.names = paste0("Ali_I_", colnames(x = Ali_I))
# )
# 
# Ali_II$Arm <- "Alisertib"
# Ali_II <- RenameCells(
#   object = Ali_II,
#   new.names = paste0("Ali_II_", colnames(x = Ali_II))
# )
# 
# # Cell Filtering QC
# # RdI
# 
# Ali_I$percent.mt <- PercentageFeatureSet(Ali_I, pattern = "^MT-")
# Ali_I[["percent.hg19"]] <- PercentageFeatureSet(Ali_I, pattern = "^hg19-")
# Ali_I[["percent.mm10"]] <- PercentageFeatureSet(Ali_I, pattern = "^mm10-")
# 
# pdf(
#   file = "01_output/Ali_I_prefilterQC.pdf"
# )
# VlnPlot(Ali_I, features = c("percent.hg19", "percent.mm10"))
# FeatureScatter(object = Ali_I, feature1 = "nCount_RNA", feature2 = "percent.mt")
# FeatureScatter(object = Ali_I, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# dev.off()
# 
# # RdII
# 
# Ali_II$percent.mt <- PercentageFeatureSet(Ali_II, pattern = "^MT-")
# Ali_II[["percent.hg19"]] <- PercentageFeatureSet(Ali_II, pattern = "^hg19-")
# Ali_II[["percent.mm10"]] <- PercentageFeatureSet(Ali_II, pattern = "^mm10-")
# 
# pdf(
#   file = "01_output/Ali_II_prefilterQC.pdf"
# )
# VlnPlot(Ali_II, features = c("percent.hg19", "percent.mm10"))
# FeatureScatter(object = Ali_II, feature1 = "nCount_RNA", feature2 = "percent.mt")
# FeatureScatter(object = Ali_II, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# dev.off()
# 
# ################################################################################
# ################################################################################
# ################################################################################
# 
# ##                         INTEGRATION
# 
# ################################################################################
# ################################################################################
# ################################################################################
# 
# obj.list <- c(DMSO_tween, Ali_I, Ali_II)
# rm(DMSO_tween.data, DMSO_tween, Ali_I, Ali_RdI.data, Ali_II, Ali_RdII.data)
# # And uniformly filter all datasets...
# 
# for (i in 1:length(obj.list)) {
#   obj.list[[i]] <- subset(x = obj.list[[i]],
#                                  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
# }
# 
# 
# # Run scTransform on all objects...
# 
# for (i in 1:length(obj.list)) {
#   obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = TRUE)
# }
# 
# # Next, select features for downstream integration, and run PrepSCTIntegration,
# # which ensures that all necessary Pearson residuals have been calculated.
# 
# obj.features <- SelectIntegrationFeatures(
#   object.list = obj.list,
#   nfeatures = 3000
# )
# 
# obj.list <- PrepSCTIntegration(
#   object.list = obj.list,
#   anchor.features = obj.features,
#   verbose = TRUE)
# 
# # Next, identify anchors and integrate the datasets. Commands are identical to
# # the standard workflow, but make sure to set normalization.method = "SCT":
# 
# obj.anchors <- FindIntegrationAnchors(
#   object.list = obj.list,
#   normalization.method = "SCT",
#   anchor.features = obj.features,
#   verbose = TRUE
# )
# 
# to_integrate <- Reduce(intersect, lapply(obj.anchors@object.list, rownames))
# 
# print("Preview of genes to integrate... ")
# head(to_integrate)
# 
# obj <- IntegrateData(
#   anchorset = obj.anchors,
#   #  features.to.integrate = to_integrate,
#   normalization.method = "SCT",
#   verbose = TRUE
# )
# 
# rm(obj.anchors)
# rm(obj.features)
# rm(obj.list)
# rm(to_integrate)
# 
# # saveRDS(obj, file = "AlisertibExpt/01_output/CellTagInt.RDS")
# 
# # Now proceed with downstream analysis (i.e. visualization, clustering) on the
# # integrated dataset. Commands are identical to the standard workflow, but do
# # not run the ScaleData function after integration. This should have helped to
# # remove any variation caused by technical batch effect.
# 
# ################################################################################
# ################################################################################
# ################################################################################
# 
# obj <- RunPCA(obj, verbose = FALSE)
# 
# pdf(
#   file = "01_output/CellTagInt_ElbowPlot.pdf"
# )
# ElbowPlot(obj)
# dev.off()
# 
# obj <- FindNeighbors(obj, dims = 1:15)
# obj <- FindClusters(obj, dims = 1:15)
# 
# obj <- RunUMAP(
#   obj,
#   dims = 1:15,
#   min.dist = 0.5,
#   n_neighbors = 30,
#   umap.method = "uwot",
#   metric = "correlation",
#   assay = "SCT"
# )
# 
# pdf(
#   file = "01_output/CellTagIntUMAP.pdf"
# )
# DimPlot(obj, group.by = c("Arm"), cols = c("azure4", "chartreuse3"))#, combine = FALSE)
# DimPlot(obj, group.by = "ident", label = TRUE)#, combine = FALSE)
# FeaturePlot(obj, features = "percent.hg19", cols = viridis(100, option = "D"))
# FeaturePlot(obj, features = "percent.mm10", cols = viridis(100, option = "D"))
# 
# dev.off()
# 
# # saveRDS(obj, file = "AlisertibExpt/01_output/CellTagInt_Barnyard_UMAP.Rds")
# 
# ################################################################################
# 
# # Generate supplemental figure showing seperation of pdx and mouse cells. 
# 
# # source umap
# 
# # sourceUMAP <- as.ggplot(DimPlot(obj, group.by = "orig.ident") + ggtitle(label = "Original Identity"))
# armUMAP <- as.ggplot(DimPlot(obj, group.by = c("Arm"), cols = c("violetred3", "azure4"), order = "DMSO"))
# snnUMAP <- as.ggplot(DimPlot(obj, group.by = "ident", label = TRUE))
# hg19UMAP <- as.ggplot(FeaturePlot(obj, features = "percent.hg19", cols = viridis(100, option = "D")))# + ggtitle(label = "Percent alignment to hg19"))
# mm10UMAP <- as.ggplot(FeaturePlot(obj, features = "percent.mm10", cols = viridis(100, option = "D")))# + ggtitle(label = "Percent alignment to mm10"))
# 
# dev.off()
# pdf(file = "01_output/mm10_hg19_seperation.pdf", width = 14)
# plot_grid(plotlist = list(armUMAP, 
#                           hg19UMAP, 
#                           mm10UMAP
# ), ncol = 3, labels = "AUTO")
# dev.off()
# 
# ################################################################################
# ################################################################################
# ################################################################################
# 
# print("Pipeline successfully completed...")
# 
# ################################################################################
# ################################################################################
# ################################################################################