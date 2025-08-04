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

# CHLA10 <- readRDS(file = "../working_data/chla10_tgfb_seurat.rds")
# 
# # re-sort data, make sure it is processed for the subset
# ## rna data
# DefaultAssay(CHLA10) <- "RNA"
# CHLA10 <- FindVariableFeatures(object = CHLA10)
# CHLA10 <- ScaleData(object = CHLA10)
# CHLA10 <- RunPCA(object = CHLA10)
# CHLA10 <- FindNeighbors(object = CHLA10, dims = 1:30)
# CHLA10 <- FindClusters(object = CHLA10)
# CHLA10 <- RunUMAP(object = CHLA10, dims = 1:30)
# 
# ## atac data
# DefaultAssay(CHLA10) <- "ATAC"
# CHLA10 <- RunTFIDF(CHLA10)
# CHLA10 <- FindTopFeatures(CHLA10, min.cutoff = 'q0')
# CHLA10 <- RunSVD(CHLA10)
# CHLA10 <- RunUMAP(CHLA10, reduction = 'lsi', dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# 
# ## cluster factoring in both modules
# CHLA10 <- FindMultiModalNeighbors(CHLA10, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:20))
# CHLA10 <- RunUMAP(CHLA10, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# CHLA10 <- FindClusters(CHLA10, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


# load data back in
CHLA10 <- readRDS("../working_data/CHLA10_reclustered.rds")
Idents(CHLA10) <- "condition"

# add features/metadata for cells
DefaultAssay(CHLA10) <- "RNA"

## gene sets
riggi_up <- read.delim(file="../gene_lists/RIGGI_UP.txt", header=FALSE)
riggi_dn <- read.delim(file="../gene_lists/RIGGI_DN.txt", header=FALSE)
cc.genes <- Seurat::cc.genes.updated.2019


## cell cycle phases
CHLA10 <- CellCycleScoring(
  object = CHLA10,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE)

## Riggi gene list metadata/subsetting
CHLA10 <- AddModuleScore(CHLA10, features=riggi_up, name="riggi_up")
CHLA10 <- AddModuleScore(CHLA10, features=riggi_dn, name="riggi_dn")
rup_cutoff <- quantile(CHLA10$riggi_up1, probs=0.5)
rup_cutoff <- unname(rup_cutoff)
rdn_cutoff <- quantile(CHLA10$riggi_dn1, probs=0.5)
rdn_cutoff <- unname(rdn_cutoff)

riggi_up_hi <- subset(x = CHLA10, subset = riggi_up1 > rup_cutoff, slot='data')
riggi_up_lo <- subset(x = CHLA10, subset = riggi_up1 > rup_cutoff, slot='data', invert=TRUE)
CHLA10 <- SetIdent(CHLA10, cells=Cells(riggi_up_hi), value='riggi_up_hi')
CHLA10 <- SetIdent(CHLA10, cells=Cells(riggi_up_lo), value='riggi_up_lo')
CHLA10$riggi_up_identity <- Idents(CHLA10)

riggi_dn_hi <- subset(x = CHLA10, subset = riggi_dn1 > rdn_cutoff, slot='data')
riggi_dn_lo <- subset(x = CHLA10, subset = riggi_dn1 > rdn_cutoff, slot='data', invert=TRUE)
CHLA10 <- SetIdent(CHLA10, cells=Cells(riggi_dn_hi), value='riggi_dn_hi')
CHLA10 <- SetIdent(CHLA10, cells=Cells(riggi_dn_lo), value='riggi_dn_lo')
CHLA10$riggi_dn_identity <- Idents(CHLA10)

hybrid<- subset(x = CHLA10, subset = riggi_up_identity  == "riggi_up_hi" & riggi_dn_identity == "riggi_dn_hi", slot='data')
true_low<- subset(x = CHLA10, subset = riggi_up_identity  == "riggi_up_lo" & riggi_dn_identity == "riggi_dn_hi", slot='data')
true_hi<- subset(x = CHLA10, subset = riggi_up_identity  == "riggi_up_hi" & riggi_dn_identity == "riggi_dn_lo", slot='data')
low_hybrid<- subset(x = CHLA10, subset = riggi_up_identity  == "riggi_up_lo" & riggi_dn_identity == "riggi_dn_lo", slot='data')

CHLA10 <- SetIdent(CHLA10, cells=Cells(hybrid), value='hybrid')
CHLA10 <- SetIdent(CHLA10, cells=Cells(true_low), value='true_low')
CHLA10 <- SetIdent(CHLA10, cells=Cells(true_hi), value='true_high')
CHLA10 <- SetIdent(CHLA10, cells=Cells(low_hybrid), value='low_hybrid')
Idents(CHLA10) <- factor(Idents(CHLA10), levels = c("low_hybrid", "true_low", "true_high", "hybrid"))
CHLA10$ef_state_identity <- Idents(CHLA10)

## ZEB2 subsetting
zeb2_expr <- GetAssayData(CHLA10, slot = "data")["ZEB2", ]
zeb2_bot_30 <- quantile(zeb2_expr, 0.30)
zeb2_top_30 <- quantile(zeb2_expr, 0.70)
zeb2_bot_30 <- unname(zeb2_bot_30)
zeb2_top_30 <- unname(zeb2_top_30)

hi_ZEB2 <- subset(CHLA10, cells = names(zeb2_expr[zeb2_expr > zeb2_top_30]))
low_ZEB2 <- subset(CHLA10, cells = names(zeb2_expr[zeb2_expr < zeb2_bot_30]))
mid_zeb2_cells <- names(zeb2_expr[zeb2_expr >= zeb2_bot_30 & zeb2_expr <= zeb2_top_30])
mid_ZEB2 <- subset(CHLA10, cells = mid_zeb2_cells)

CHLA10 <- SetIdent(CHLA10, cells=Cells(hi_ZEB2), value='Hi_ZEB2')
CHLA10 <- SetIdent(CHLA10, cells=Cells(mid_ZEB2), value='Mid_ZEB2')
CHLA10 <- SetIdent(CHLA10, cells=Cells(low_ZEB2), value='Low_ZEB2')

Idents(CHLA10) <- factor(Idents(CHLA10), levels = c("Low_ZEB2", "Mid_ZEB2", "Hi_ZEB2"))
CHLA10$zeb2_lvl <- Idents(CHLA10)

## ZEB1 subsetting
zeb1_expr <- GetAssayData(CHLA10, slot = "data")["ZEB1", ]
ZEB1_positive <- subset(CHLA10, cells = names(zeb1_expr[zeb1_expr > 0]))
ZEB1_negative <- subset(CHLA10, cells = names(zeb1_expr[zeb1_expr <= 0]))

CHLA10 <- SetIdent(CHLA10, cells=Cells(ZEB1_positive), value='ZEB1+')
CHLA10 <- SetIdent(CHLA10, cells=Cells(ZEB1_negative), value='ZEB1-')
CHLA10$zeb1_expr <- Idents(CHLA10)


# Fix object annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(CHLA10[["ATAC"]]) <- annotations


# Identify differentially expressed genes
DefaultAssay(CHLA10) <- 'RNA'
Idents(CHLA10) <- "condition"
tgfb1_de_genes <- FindMarkers(CHLA10, ident.1 = 'TGFB1', ident.2 = 'vehicle')
tgfb1_upreg <- tgfb1_de_genes[tgfb1_de_genes$avg_log2FC > 0 & tgfb1_de_genes$p_val_adj < 0.05, "X"]


# Identify differentially accessible peaks
DefaultAssay(CHLA10) <- 'ATAC'
Idents(CHLA10) <- 'condition'
tgfb1_da_peaks <- FindMarkers(
  object = CHLA10, ident.1 = "TGFB1", ident.2 = "vehicle", test.use = 'wilcox', min.pct = 0.1)

## Record significant peaks, label them with proximal genes
open_tgfb1 <- rownames(tgfb1_da_peaks[tgfb1_da_peaks$avg_log2FC > 0 & tgfb1_da_peaks$p_val_adj < 0.05, ])
closest_genes_tgfb1 <- ClosestFeature(CHLA10, regions = open_tgfb1)
tgfb1_accessible_genes <- closest_genes_tgfb1$gene_name


# Correlate peak accessibility with gene expression DIRECTLY -- intensive
CHLA10 <- LinkPeaks(object = CHLA10, peak.assay = "ATAC", expression.assay = "RNA")
CHLA10_links <-  Links(CHLA10)
# write.csv(CHLA10_links, "../results/CHLA10_link_ranges.csv")


# Identify motifs enriched in differentially accessible ranges/peaks
upreg_pks <- tgfb1_da_peaks[tgfb1_da_peaks$avg_log2FC > 0 & tgfb1_da_peaks$p_val_adj < 0.05, ]
valid_peaks <- intersect(as.character(rownames(upreg_pks)), rownames(CHLA10[["ATAC"]]))
DefaultAssay(CHLA10) <- "ATAC"
enriched_motifs <- FindMotifs(CHLA10, features = valid_peaks)
# write.csv(enriched_motifs, "../results/CHLA10_enriched_motifs.csv")

## Save progress
# saveRDS(CHLA10, "../working_data/CHLA10_processed.rds")

# 
# # load data back in
# CHLA10 <- readRDS("../working_data/CHLA10_processed.rds")
# DefaultAssay(CHLA10) <- "ATAC"


# Calculate motif activity per cell -- very intensive
register(SerialParam())
CHLA10 <- RunChromVAR(CHLA10, genome = BSgenome.Hsapiens.UCSC.hg38)


# If ChromVAR works, save again
# saveRDS(CHLA10, "../working_data/CHLA10_processed.rds")

