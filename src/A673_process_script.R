library(Seurat)
library(scCustomize)
library(ggplot2)
library(dplyr)
library(viridis)
# library(SeuratData)
library(magrittr)
library(patchwork)
library(Signac)
# library(qs)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(plyranges)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(SeuratDisk)
library(SeuratObject)
library(reshape2)
library(eulerr)
library(org.Hs.eg.db)
library(BiocParallel)

# A673 <- readRDS(file = "../working_data/A673_tgfb_seurat.rds")
# 
# # re-sort data, make sure it is processed for the subset
# ## rna data
# DefaultAssay(A673) <- "RNA"
# A673 <- FindVariableFeatures(object = A673)
# A673 <- ScaleData(object = A673)
# A673 <- RunPCA(object = A673)
# A673 <- FindNeighbors(object = A673, dims = 1:30)
# A673 <- FindClusters(object = A673)
# A673 <- RunUMAP(object = A673, dims = 1:30)
# 
# ## atac data
# DefaultAssay(A673) <- "ATAC"
# A673 <- RunTFIDF(A673)
# A673 <- FindTopFeatures(A673, min.cutoff = 'q0')
# A673 <- RunSVD(A673)
# A673 <- RunUMAP(A673, reduction = 'lsi', dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# 
# ## cluster factoring in both modules
# A673 <- FindMultiModalNeighbors(A673, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:20))
# A673 <- RunUMAP(A673, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# A673 <- FindClusters(A673, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


# # load data back in
# A673 <- readRDS("../working_data/A673_reclustered.rds")
# Idents(A673) <- "condition"
# 
# # add features/metadata for cells
# DefaultAssay(A673) <- "RNA"
# 
# ## gene sets
# riggi_up <- read.delim(file="../gene_lists/RIGGI_UP.txt", header=FALSE)
# riggi_dn <- read.delim(file="../gene_lists/RIGGI_DN.txt", header=FALSE)
# cc.genes <- Seurat::cc.genes.updated.2019
# 
# 
# ## cell cycle phases
# A673 <- CellCycleScoring(
#   object = A673,
#   s.features = cc.genes$s.genes,
#   g2m.features = cc.genes$g2m.genes,
#   set.ident = TRUE)
# 
# ## Riggi gene list metadata/subsetting
# A673 <- AddModuleScore(A673, features=riggi_up, name="riggi_up")
# A673 <- AddModuleScore(A673, features=riggi_dn, name="riggi_dn")
# rup_cutoff <- quantile(A673$riggi_up1, probs=0.5)
# rup_cutoff <- unname(rup_cutoff)
# rdn_cutoff <- quantile(A673$riggi_dn1, probs=0.5)
# rdn_cutoff <- unname(rdn_cutoff)
# 
# riggi_up_hi <- subset(x = A673, subset = riggi_up1 > rup_cutoff, slot='data')
# riggi_up_lo <- subset(x = A673, subset = riggi_up1 > rup_cutoff, slot='data', invert=TRUE)
# A673 <- SetIdent(A673, cells=Cells(riggi_up_hi), value='riggi_up_hi')
# A673 <- SetIdent(A673, cells=Cells(riggi_up_lo), value='riggi_up_lo')
# A673$riggi_up_identity <- Idents(A673)
# 
# riggi_dn_hi <- subset(x = A673, subset = riggi_dn1 > rdn_cutoff, slot='data')
# riggi_dn_lo <- subset(x = A673, subset = riggi_dn1 > rdn_cutoff, slot='data', invert=TRUE)
# A673 <- SetIdent(A673, cells=Cells(riggi_dn_hi), value='riggi_dn_hi')
# A673 <- SetIdent(A673, cells=Cells(riggi_dn_lo), value='riggi_dn_lo')
# A673$riggi_dn_identity <- Idents(A673)
# 
# hybrid<- subset(x = A673, subset = riggi_up_identity  == "riggi_up_hi" & riggi_dn_identity == "riggi_dn_hi", slot='data')
# true_low<- subset(x = A673, subset = riggi_up_identity  == "riggi_up_lo" & riggi_dn_identity == "riggi_dn_hi", slot='data')
# true_hi<- subset(x = A673, subset = riggi_up_identity  == "riggi_up_hi" & riggi_dn_identity == "riggi_dn_lo", slot='data')
# low_hybrid<- subset(x = A673, subset = riggi_up_identity  == "riggi_up_lo" & riggi_dn_identity == "riggi_dn_lo", slot='data')
# 
# A673 <- SetIdent(A673, cells=Cells(hybrid), value='hybrid')
# A673 <- SetIdent(A673, cells=Cells(true_low), value='true_low')
# A673 <- SetIdent(A673, cells=Cells(true_hi), value='true_high')
# A673 <- SetIdent(A673, cells=Cells(low_hybrid), value='low_hybrid')
# Idents(A673) <- factor(Idents(A673), levels = c("low_hybrid", "true_low", "true_high", "hybrid"))
# A673$ef_state_identity <- Idents(A673)
# 
# ## ZEB2 subsetting
# zeb2_expr <- GetAssayData(A673, slot = "data")["ZEB2", ]
# ZEB2_positive <- subset(A673, cells = names(zeb2_expr[zeb2_expr > 0]))
# ZEB2_negative <- subset(A673, cells = names(zeb2_expr[zeb2_expr <= 0]))
# 
# A673 <- SetIdent(A673, cells=Cells(ZEB2_positive), value='ZEB2+')
# A673 <- SetIdent(A673, cells=Cells(ZEB2_negative), value='ZEB2-')
# A673$zeb2_expr <- Idents(A673)
# 
# ## ZEB1 subsetting
# zeb1_expr <- GetAssayData(A673, slot = "data")["ZEB1", ]
# ZEB1_positive <- subset(A673, cells = names(zeb1_expr[zeb1_expr > 0]))
# ZEB1_negative <- subset(A673, cells = names(zeb1_expr[zeb1_expr <= 0]))
# 
# A673 <- SetIdent(A673, cells=Cells(ZEB1_positive), value='ZEB1+')
# A673 <- SetIdent(A673, cells=Cells(ZEB1_negative), value='ZEB1-')
# A673$zeb1_expr <- Idents(A673)
# 
# 
# # Fix object annotations
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- "UCSC"
# genome(annotations) <- "hg38"
# Annotation(A673[["ATAC"]]) <- annotations
# 
# 
# # Identify differentially expressed genes
# DefaultAssay(A673) <- 'RNA'
# Idents(A673) <- "condition"
# tgfb1_de_genes <- FindMarkers(A673, ident.1 = 'TGFB1', ident.2 = 'vehicle')
# tgfb1_upreg <- tgfb1_de_genes[tgfb1_de_genes$avg_log2FC > 0 & tgfb1_de_genes$p_val_adj < 0.05, "X"]
# 
# 
# # Identify differentially accessible peaks
# DefaultAssay(A673) <- 'ATAC'
# Idents(A673) <- 'condition'
# tgfb1_da_peaks <- FindMarkers(
#   object = A673, ident.1 = "TGFB1", ident.2 = "vehicle", test.use = 'wilcox', min.pct = 0.1)
# 
# ## Record significant peaks, label them with proximal genes
# open_tgfb1 <- rownames(tgfb1_da_peaks[tgfb1_da_peaks$avg_log2FC > 0 & tgfb1_da_peaks$p_val_adj < 0.05, ])
# closest_genes_tgfb1 <- ClosestFeature(A673, regions = open_tgfb1)
# tgfb1_accessible_genes <- closest_genes_tgfb1$gene_name
# 
# 
# # Correlate peak accessibility with gene expression DIRECTLY -- intensive
# A673 <- LinkPeaks(object = A673, peak.assay = "ATAC", expression.assay = "RNA")
# write.csv(Links(A673), "../results/A673_link_ranges.csv")
# 
# 
# # Identify motifs enriched in differentially accessible ranges/peaks
# upreg_pks <- tgfb1_da_peaks[tgfb1_da_peaks$avg_log2FC > 0 & tgfb1_da_peaks$p_val_adj < 0.05, ]
# valid_peaks <- intersect(as.character(rownames(upreg_pks)), rownames(A673[["ATAC"]]))
# DefaultAssay(A673) <- "ATAC"
# enriched_motifs <- FindMotifs(A673, features = valid_peaks)
# write.csv(enriched_motifs, "../results/A673_enriched_motifs.csv")
# 
# ## Save progress
# saveRDS(A673, "../working_data/A673_processed.rds")


# load data back in
A673 <- readRDS("../working_data/A673_processed.rds")
DefaultAssay(A673) <- "ATAC"


# Calculate motif activity per cell -- very intensive
register(SerialParam())
A673 <- RunChromVAR(A673, genome = BSgenome.Hsapiens.UCSC.hg38)


# If ChromVAR works, save again
saveRDS(A673, "../working_data/A673_processed.rds")

