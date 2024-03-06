
library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(parallel)
library(openxlsx)
library(plotly)
library(sctransform)
library(glassoFast)
library(igraph)
library(ggraph)
library(graphlayouts)
library(Signac)
library(patchwork)
library(hdf5r)
library(GenomeInfoDbData)
library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)

source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
## cell cycle perturbed data 

## ID of the KO genes
Crk2 <- 'TGGT1-218220'
Ark3 <- 'TGGT1-203010'

## IDs
prod.desc  <- read.xlsx('../Input_copy/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_copy/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


getExpr <- function(in.file, TGGT1_ME49){
  file.counts <- read.csv(in.file)
  genes <- file.counts$X
  ind <- which(genes %in% TGGT1_ME49$TGME49)
  file.counts <- file.counts[ind, ]
  genes <- genes[ind]
  genes <- TGGT1_ME49$TGGT1[match(genes, TGGT1_ME49$TGME49)]
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  return(expr)
}

# mnk 1 hr

intra.file.csv <- "../Input_copy/toxo_cdc_new_221220/scRNA/1_hr_KD_MNK_expr.csv"
mnk.1hr.counts <- getExpr(intra.file.csv, TGGT1_ME49)
dim(mnk.1hr.counts)

feats <- c("nFeature_RNA","nCount_RNA")
S.O.mnk.1hr <- CreateSeuratObject(counts = mnk.1hr.counts)
S.O.mnk.1hr$orig.ident <- 'scRNA.mnk.1hr'
S.O.mnk.1hr$spp <- "scRNA.mnk.1hr"
S.O.mnk.1hr$spp2 <- "mnk"

VlnPlot(S.O.mnk.1hr, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.mnk.1hr, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.mnk.1hr, expression = nFeature_RNA > 100 & nFeature_RNA < 1500)
selected_f <- rownames(S.O.mnk.1hr)[ Matrix::rowSums(S.O.mnk.1hr) > 5]
selected_f <- unique(c(selected_f, Crk2, Ark3)) ## Make sure Crk2 and Ark 3 are selected
S.O.mnk.1hr  <- subset(S.O.mnk.1hr, features=selected_f, cells=selected_c)
#S.O.mnk.1hr <- subset(S.O.mnk.1hr, downsample = 3000)
S.O.mnk.1hr <- prep_S.O(S.O.mnk.1hr, res = 0.4)
DimPlot(S.O.mnk.1hr, reduction = "pca", dims = c(1,2))
DimPlot(S.O.mnk.1hr, reduction = "umap")

## mnk 3 hr

intra.file.csv <- "../Input_copy/toxo_cdc_new_221220/scRNA/3_hr_KD_MNK_expr.csv"
mnk.3hr.counts <- getExpr(intra.file.csv, TGGT1_ME49)
dim(mnk.3hr.counts)

feats <- c("nFeature_RNA","nCount_RNA")
S.O.mnk.3hr <- CreateSeuratObject(counts = mnk.3hr.counts)
S.O.mnk.3hr$orig.ident <- 'scRNA.3hr.mnk'
S.O.mnk.3hr$spp <- "scRNA.3hr.mnk"
S.O.mnk.3hr$spp2 <- "mnk"

VlnPlot(S.O.mnk.3hr, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.mnk.3hr, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.mnk.3hr, expression = nFeature_RNA > 100 & nFeature_RNA < 1500)
selected_f <- rownames(S.O.mnk.3hr)[ Matrix::rowSums(S.O.mnk.3hr) > 5]
selected_f <- unique(c(selected_f, Crk2, Ark3)) ## Make sure Crk2 and Ark 3 are selected
S.O.mnk.3hr  <- subset(S.O.mnk.3hr, features=selected_f, cells=selected_c)
#S.O.intra <- subset(S.O.intra, downsample = 3000)

S.O.mnk.3hr <- prep_S.O(S.O.mnk.3hr, res = 0.4)
DimPlot(S.O.mnk.3hr, reduction = "pca")
DimPlot(S.O.mnk.3hr, reduction = "umap")


##
S.O.tg.boothroyd <- readRDS('../Input_copy/boothroyd_sc_all_data/rds/S.O.tg_RH_boothroyd.rds')

S.Os <- list(S.O.mnk.1hr = S.O.mnk.1hr, S.O.mnk.3hr = S.O.mnk.3hr)

S.Os <- mclapply(S.Os, function(S.O){
  S.O <- prep_S.O(S.O, res = 0.4)
  anchors <- FindTransferAnchors(reference = S.O.tg.boothroyd, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.tg.boothroyd@meta.data$phase,dims = 1:30)
  predictions$phase <- predictions$predicted.id
  #predictions$phase[which(predictions$prediction.score.max < 0.7)] <- 'NA'
  S.O <- AddMetaData(object = S.O, metadata = predictions)
  return(S.O)
}, mc.cores = num.cores)

spps <- names(S.Os)
names(S.Os) <- spps


S.Os$S.O.mnk.1hr$phase <- factor(S.Os$S.O.mnk.1hr$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))
S.Os$S.O.mnk.3hr$phase <- factor(S.Os$S.O.mnk.3hr$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))

Idents(S.Os$S.O.mnk.1hr) <- "phase"
Idents(S.Os$S.O.mnk.3hr) <- "phase"

DimPlot(S.Os$S.O.mnk.1hr)
DimPlot(S.Os$S.O.mnk.3hr)

## 3D plot 1hr
umap.data <- FetchData(S.Os$S.O.mnk.1hr, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "phase","seurat_clusters"))
umap.data$phase <- factor(umap.data$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))
plot_ly(umap.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,color =  umap.data$phase, size = 0.2)

pca.data <- FetchData(S.Os$S.O.mnk.1hr, vars = c("PC_1", "PC_2", "PC_3", "phase"))
pca.data$phase <- factor(pca.data$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))
plot_ly(pca.data, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3,color =  pca.data$phase, size = 0.2)


## 3D plot 3hr
umap.data <- FetchData(S.Os$S.O.mnk.3hr, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "phase","seurat_clusters"))
umap.data$phase <- factor(umap.data$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))
plot_ly(umap.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,color =  umap.data$phase, size = 0.2)

pca.data <- FetchData(S.Os$S.O.mnk.3hr, vars = c("PC_1", "PC_2", "PC_3", "phase"))
pca.data$phase <- factor(pca.data$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))
plot_ly(pca.data, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3,color =  pca.data$phase, size = 0.2)

## DEGs - mnk1hr
mnk1hr.markers <- FindAllMarkers(S.Os$S.O.mnk.1hr)
mnk1hr.markers.sig <- mnk1hr.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
write.xlsx(mnk1hr.markers.sig, "../Output_copy/toxo_cdc/ME49_59/tables/mnk1hr_markers_sig.xlsx")

mnk3hr.markers <- FindAllMarkers(S.Os$S.O.mnk.3hr)
mnk3hr.markers.sig <- mnk3hr.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
write.xlsx(mnk3hr.markers.sig, "../Output/toxo_cdc/ME49_59/tables/mnk3hr_markers_sig.xlsx")

