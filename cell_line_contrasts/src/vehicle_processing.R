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
library(JASPAR2022)
library(TFBSTools)


# data input comes from treatment_splitting.R
# vehicle_treated <- readRDS(file = "../working_data/vehicle_treated.rds")
# 
# # re-sort data, make sure it is processed for the subset
# ## rna data
# DefaultAssay(vehicle_treated) <- "RNA"
# vehicle_treated <- FindVariableFeatures(object = vehicle_treated)
# vehicle_treated <- ScaleData(object = vehicle_treated)
# vehicle_treated <- RunPCA(object = vehicle_treated)
# vehicle_treated <- FindNeighbors(object = vehicle_treated, dims = 1:30)
# vehicle_treated <- FindClusters(object = vehicle_treated)
# vehicle_treated <- RunUMAP(object = vehicle_treated, dims = 1:30)
# 
# ## atac data
# DefaultAssay(vehicle_treated) <- "ATAC"
# vehicle_treated <- RunTFIDF(vehicle_treated)
# vehicle_treated <- FindTopFeatures(vehicle_treated, min.cutoff = 'q0')
# vehicle_treated <- RunSVD(vehicle_treated)
# vehicle_treated <- RunUMAP(vehicle_treated, reduction = 'lsi', dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# 
# ## cluster factoring in both modules
# vehicle_treated <- FindMultiModalNeighbors(vehicle_treated, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:20))
# vehicle_treated <- RunUMAP(vehicle_treated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# vehicle_treated <- FindClusters(vehicle_treated, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# 
# # add features/metadata for cells // load data back in
# vehicle_treated <- readRDS("../working_data/vehicle_reclustered.rds")
# DefaultAssay(vehicle_treated) <- "RNA"
# 
# ## gene sets
# riggi_up <- read.delim(file="../../gene_lists/RIGGI_UP.txt", header=FALSE)
# riggi_dn <- read.delim(file="../../gene_lists/RIGGI_DN.txt", header=FALSE)
# vehicle_response_genes <- read.csv("../../gene_lists/Ewing_TGFB1_response_ACP.txt")
# cc.genes <- Seurat::cc.genes.updated.2019
# 
# 
# ## cell cycle phases
# vehicle_treated <- CellCycleScoring(
#   object = vehicle_treated,
#   s.features = cc.genes$s.genes,
#   g2m.features = cc.genes$g2m.genes,
#   set.ident = TRUE)
# 
# ## Riggi gene list metadata/subsetting
# vehicle_treated <- AddModuleScore(vehicle_treated, features=riggi_up, name="riggi_up")
# vehicle_treated <- AddModuleScore(vehicle_treated, features=riggi_dn, name="riggi_dn")
# rup_cutoff <- quantile(vehicle_treated$riggi_up1, probs=0.5)
# rup_cutoff <- unname(rup_cutoff)
# rdn_cutoff <- quantile(vehicle_treated$riggi_dn1, probs=0.5)
# rdn_cutoff <- unname(rdn_cutoff)
# 
# riggi_up_hi <- subset(x = vehicle_treated, subset = riggi_up1 > rup_cutoff, slot='data')
# riggi_up_lo <- subset(x = vehicle_treated, subset = riggi_up1 > rup_cutoff, slot='data', invert=TRUE)
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(riggi_up_hi), value='riggi_up_hi')
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(riggi_up_lo), value='riggi_up_lo')
# vehicle_treated$riggi_up_identity <- Idents(vehicle_treated)
# 
# riggi_dn_hi <- subset(x = vehicle_treated, subset = riggi_dn1 > rdn_cutoff, slot='data')
# riggi_dn_lo <- subset(x = vehicle_treated, subset = riggi_dn1 > rdn_cutoff, slot='data', invert=TRUE)
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(riggi_dn_hi), value='riggi_dn_hi')
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(riggi_dn_lo), value='riggi_dn_lo')
# vehicle_treated$riggi_dn_identity <- Idents(vehicle_treated)
# 
# hybrid<- subset(x = vehicle_treated, subset = riggi_up_identity  == "riggi_up_hi" & riggi_dn_identity == "riggi_dn_hi", slot='data')
# true_low<- subset(x = vehicle_treated, subset = riggi_up_identity  == "riggi_up_lo" & riggi_dn_identity == "riggi_dn_hi", slot='data')
# true_hi<- subset(x = vehicle_treated, subset = riggi_up_identity  == "riggi_up_hi" & riggi_dn_identity == "riggi_dn_lo", slot='data')
# low_hybrid<- subset(x = vehicle_treated, subset = riggi_up_identity  == "riggi_up_lo" & riggi_dn_identity == "riggi_dn_lo", slot='data')
# 
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(hybrid), value='hybrid')
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(true_low), value='true_low')
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(true_hi), value='true_high')
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(low_hybrid), value='low_hybrid')
# Idents(vehicle_treated) <- factor(Idents(vehicle_treated), levels = c("low_hybrid", "true_low", "true_high", "hybrid"))
# vehicle_treated$ef_state_identity <- Idents(vehicle_treated)
# 
# ## ZEB2 subsetting
# zeb2_expr <- GetAssayData(vehicle_treated, slot = "data")["ZEB2", ]
# zeb2_bot_30 <- quantile(zeb2_expr, 0.30)
# zeb2_top_30 <- quantile(zeb2_expr, 0.70)
# zeb2_bot_30 <- unname(zeb2_bot_30)
# zeb2_top_30 <- unname(zeb2_top_30)
# 
# hi_ZEB2 <- subset(vehicle_treated, cells = names(zeb2_expr[zeb2_expr > zeb2_top_30]))
# low_ZEB2 <- subset(vehicle_treated, cells = names(zeb2_expr[zeb2_expr < zeb2_bot_30]))
# mid_zeb2_cells <- names(zeb2_expr[zeb2_expr >= zeb2_bot_30 & zeb2_expr <= zeb2_top_30])
# mid_ZEB2 <- subset(vehicle_treated, cells = mid_zeb2_cells)
# 
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(hi_ZEB2), value='Hi_ZEB2')
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(mid_ZEB2), value='Mid_ZEB2')
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(low_ZEB2), value='Low_ZEB2')
# 
# Idents(vehicle_treated) <- factor(Idents(vehicle_treated), levels = c("Low_ZEB2", "Mid_ZEB2", "Hi_ZEB2"))
# vehicle_treated$zeb2_lvl <- Idents(vehicle_treated)
# 
# ## ZEB1 subsetting
# zeb1_expr <- GetAssayData(vehicle_treated, slot = "data")["ZEB1", ]
# ZEB1_positive <- subset(vehicle_treated, cells = names(zeb1_expr[zeb1_expr > 0]))
# ZEB1_negative <- subset(vehicle_treated, cells = names(zeb1_expr[zeb1_expr <= 0]))
# 
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(ZEB1_positive), value='ZEB1+')
# vehicle_treated <- SetIdent(vehicle_treated, cells=Cells(ZEB1_negative), value='ZEB1-')
# vehicle_treated$zeb1_expr <- Idents(vehicle_treated)
# 
# ## vehicle_treated response module score
# vehicle_treated <- AddModuleScore(vehicle_treated, features=vehicle_response_genes, name="tgfb_response")
# 
# 
# # Fix object annotations
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- "UCSC"
# genome(annotations) <- "hg38"
# Annotation(vehicle_treated[["ATAC"]]) <- annotations
# 
# 
# # Identify differentially expressed genes -- change more closely -- save DE genelists
# DefaultAssay(vehicle_treated) <- 'RNA'
# Idents(vehicle_treated) <- "cell_line"
# cellLine_de_genes <- FindMarkers(vehicle_treated, ident.1 = 'CHLA10', ident.2 = 'TC71')
# # cellLine_CHLA10up <- cellLine_de_genes[cellLine_de_genes$avg_log2FC > 0 & cellLine_de_genes$p_val_adj < 0.05, "X"]
# cellLine_CHLA10up <- rownames(cellLine_de_genes[cellLine_de_genes$avg_log2FC > 0 & cellLine_de_genes$p_val_adj < 0.05, ])
# # cellLine_TC71up <- cellLine_de_genes[cellLine_de_genes$avg_log2FC < 0 & cellLine_de_genes$p_val_adj < 0.05, "X"]
# cellLine_TC71up <- rownames(cellLine_de_genes[cellLine_de_genes$avg_log2FC < 0 & cellLine_de_genes$p_val_adj < 0.05, ])
# write.csv(cellLine_de_genes, "../results/vehicle_treated_cellLine_de_genes.csv")
# write.csv(cellLine_CHLA10up, "../results/vehicle_treated_CHLA10_upreg_genes.csv")
# write.csv(cellLine_TC71up, "../results/vehicle_treated_TC71_upreg_genes.csv")
# 
# 
# # Identify differentially accessible peaks
# DefaultAssay(vehicle_treated) <- 'ATAC'
# Idents(vehicle_treated) <- 'cell_line'
# cellLine_da_peaks <- FindMarkers(
#   object = vehicle_treated, ident.1 = "CHLA10", ident.2 = "TC71", test.use = 'wilcox', min.pct = 0.1)
# 
# ## Record significant peaks, label them with proximal genes
# open_chla10 <- rownames(cellLine_da_peaks[cellLine_da_peaks$avg_log2FC > 0 & cellLine_da_peaks$p_val_adj < 0.05, ])
# closest_genes_chla10 <- ClosestFeature(vehicle_treated, regions = open_chla10)
# chla10_accessible_genes <- unique(closest_genes_chla10$gene_name)
# 
# open_tc71 <- rownames(cellLine_da_peaks[cellLine_da_peaks$avg_log2FC < 0 & cellLine_da_peaks$p_val_adj < 0.05, ])
# closest_genes_tc71 <- ClosestFeature(vehicle_treated, regions = open_tc71)
# tc71_accessible_genes <- unique(closest_genes_tc71$gene_name) 
# 
# write.csv(chla10_accessible_genes, "../results/vehicle_treated_CHLA10_accessible_genes.csv")
# write.csv(tc71_accessible_genes, "../results/vehicle_treated_TC71_accessible_genes.csv")
# 
# # Correlate peak accessibility with gene expression DIRECTLY -- intensive
# vehicle_treated <- LinkPeaks(object = vehicle_treated, peak.assay = "ATAC", expression.assay = "RNA")
# vehicle_treated_links <-  Links(vehicle_treated)
# write.csv(vehicle_treated_links, "../results/vehicle_treated_link_ranges.csv")

## Save progress
# saveRDS(vehicle_treated, "../working_data/vehicle_treated_processed.rds")


# # load data back in
vehicle_treated <- readRDS("../working_data/vehicle_treated_processed.rds")
DefaultAssay(vehicle_treated) <- "ATAC"


# Calculate motif activity per cell -- very intensive
# pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
# vehicle_treated <- AddMotifs(vehicle_treated, BSgenome.Hsapiens.UCSC.hg38, pfm) # add motif set, takes 20+ min
vehicle_treated@assays[["ATAC"]]@motifs <- readRDS("../working_data/vehicle_jaspar_motifs.rds")
# saveRDS(vehicle_treated@assays[["ATAC"]]@motifs, "../working_data/vehicle_jaspar_motifs.rds") # save motif position matrix
### Calculate motif enrichment
# jaspar_CHLA10 <- FindMotifs(vehicle_treated, features = open_chla10)
# write.csv(jaspar_CHLA10, "../results/vehicle_CHLA10_jaspar_enrichment.csv")
# jaspar_TC71 <- FindMotifs(vehicle_treated, features = open_tc71)
# write.csv(jaspar_TC71, "../results/vehicle_TC71_jaspar_enrichment.csv")
### ChromVAR
register(SerialParam()) # calibrate parallelization
vehicle_treated <- RunChromVAR(vehicle_treated, genome = BSgenome.Hsapiens.UCSC.hg38, new.assay.name = "JASPAR")

pfms <- readRDS("../../gene_lists/cisBP_human_pfms_2021.rds") # new motif set, same process
# vehicle_treated <- AddMotifs(vehicle_treated, BSgenome.Hsapiens.UCSC.hg38, pfms)
# saveRDS(vehicle_treated@assays[["ATAC"]]@motifs, "../working_data/vehicle_cisbp_motifs.rds")
vehicle_treated@assays[["ATAC"]]@motifs <- readRDS("../working_data/vehicle_cisbp_motifs.rds")
# cisbp_CHLA10 <- FindMotifs(vehicle_treated, features = open_chla10)
# write.csv(cisbp_CHLA10, "../results/vehicle_CHLA10_cisbp_enrichment.csv")
# cisbp_TC71 <- FindMotifs(vehicle_treated, features = open_tc71)
# write.csv(cisbp_TC71, "../results/vehicle_TC71_cisbp_enrichment.csv")

vehicle_treated <- RunChromVAR(vehicle_treated, genome = BSgenome.Hsapiens.UCSC.hg38, new.assay.name = "cisBP")


# If ChromVAR works, save again
saveRDS(vehicle_treated, "../working_data/vehicle_treated_processed.rds")

