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

TC71 <- readRDS(file = "../working_data/tc71_tgfb_seurat.rds")

# re-sort data, make sure it is processed for the subset
## rna data
DefaultAssay(TC71) <- "RNA"
TC71 <- FindVariableFeatures(object = TC71)
TC71 <- ScaleData(object = TC71)
TC71 <- RunPCA(object = TC71)
TC71 <- FindNeighbors(object = TC71, dims = 1:30)
TC71 <- FindClusters(object = TC71)
TC71 <- RunUMAP(object = TC71, dims = 1:30)

## atac data
DefaultAssay(TC71) <- "ATAC"
TC71 <- RunTFIDF(TC71)
TC71 <- FindTopFeatures(TC71, min.cutoff = 'q0')
TC71 <- RunSVD(TC71)
TC71 <- RunUMAP(TC71, reduction = 'lsi', dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

## cluster factoring in both modules
TC71 <- FindMultiModalNeighbors(TC71, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:20))
TC71 <- RunUMAP(TC71, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
TC71 <- FindClusters(TC71, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


# load data back in
# TC71 <- readRDS("../working_data/TC71_reclustered.rds")
# Idents(TC71) <- "condition"

# add features/metadata for cells
DefaultAssay(TC71) <- "RNA"

## gene sets
riggi_up <- read.delim(file="../gene_lists/RIGGI_UP.txt", header=FALSE)
riggi_dn <- read.delim(file="../gene_lists/RIGGI_DN.txt", header=FALSE)
cc.genes <- Seurat::cc.genes.updated.2019


## cell cycle phases
TC71 <- CellCycleScoring(
  object = TC71,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE)

## Riggi gene list metadata/subsetting
TC71 <- AddModuleScore(TC71, features=riggi_up, name="riggi_up")
TC71 <- AddModuleScore(TC71, features=riggi_dn, name="riggi_dn")
rup_cutoff <- quantile(TC71$riggi_up1, probs=0.5)
rup_cutoff <- unname(rup_cutoff)
rdn_cutoff <- quantile(TC71$riggi_dn1, probs=0.5)
rdn_cutoff <- unname(rdn_cutoff)

riggi_up_hi <- subset(x = TC71, subset = riggi_up1 > rup_cutoff, slot='data')
riggi_up_lo <- subset(x = TC71, subset = riggi_up1 > rup_cutoff, slot='data', invert=TRUE)
TC71 <- SetIdent(TC71, cells=Cells(riggi_up_hi), value='riggi_up_hi')
TC71 <- SetIdent(TC71, cells=Cells(riggi_up_lo), value='riggi_up_lo')
TC71$riggi_up_identity <- Idents(TC71)

riggi_dn_hi <- subset(x = TC71, subset = riggi_dn1 > rdn_cutoff, slot='data')
riggi_dn_lo <- subset(x = TC71, subset = riggi_dn1 > rdn_cutoff, slot='data', invert=TRUE)
TC71 <- SetIdent(TC71, cells=Cells(riggi_dn_hi), value='riggi_dn_hi')
TC71 <- SetIdent(TC71, cells=Cells(riggi_dn_lo), value='riggi_dn_lo')
TC71$riggi_dn_identity <- Idents(TC71)

hybrid<- subset(x = TC71, subset = riggi_up_identity  == "riggi_up_hi" & riggi_dn_identity == "riggi_dn_hi", slot='data')
true_low<- subset(x = TC71, subset = riggi_up_identity  == "riggi_up_lo" & riggi_dn_identity == "riggi_dn_hi", slot='data')
true_hi<- subset(x = TC71, subset = riggi_up_identity  == "riggi_up_hi" & riggi_dn_identity == "riggi_dn_lo", slot='data')
low_hybrid<- subset(x = TC71, subset = riggi_up_identity  == "riggi_up_lo" & riggi_dn_identity == "riggi_dn_lo", slot='data')

TC71 <- SetIdent(TC71, cells=Cells(hybrid), value='hybrid')
TC71 <- SetIdent(TC71, cells=Cells(true_low), value='true_low')
TC71 <- SetIdent(TC71, cells=Cells(true_hi), value='true_high')
TC71 <- SetIdent(TC71, cells=Cells(low_hybrid), value='low_hybrid')
Idents(TC71) <- factor(Idents(TC71), levels = c("low_hybrid", "true_low", "true_high", "hybrid"))
TC71$ef_state_identity <- Idents(TC71)

## ZEB2 subsetting
zeb2_expr <- GetAssayData(TC71, slot = "data")["ZEB2", ]
zeb2_bot_30 <- quantile(zeb2_expr, 0.30)
zeb2_top_30 <- quantile(zeb2_expr, 0.70)
zeb2_bot_30 <- unname(zeb2_bot_30)
zeb2_top_30 <- unname(zeb2_top_30)

hi_ZEB2 <- subset(TC71, cells = names(zeb2_expr[zeb2_expr > zeb2_top_30]))
low_ZEB2 <- subset(TC71, cells = names(zeb2_expr[zeb2_expr < zeb2_bot_30]))
mid_zeb2_cells <- names(zeb2_expr[zeb2_expr >= zeb2_bot_30 & zeb2_expr <= zeb2_top_30])
mid_ZEB2 <- subset(TC71, cells = mid_zeb2_cells)

TC71 <- SetIdent(TC71, cells=Cells(hi_ZEB2), value='Hi_ZEB2')
TC71 <- SetIdent(TC71, cells=Cells(mid_ZEB2), value='Mid_ZEB2')
TC71 <- SetIdent(TC71, cells=Cells(low_ZEB2), value='Low_ZEB2')

Idents(TC71) <- factor(Idents(TC71), levels = c("Low_ZEB2", "Mid_ZEB2", "Hi_ZEB2"))
TC71$zeb2_lvl <- Idents(TC71)

## ZEB1 subsetting
zeb1_expr <- GetAssayData(TC71, slot = "data")["ZEB1", ]
ZEB1_positive <- subset(TC71, cells = names(zeb1_expr[zeb1_expr > 0]))
ZEB1_negative <- subset(TC71, cells = names(zeb1_expr[zeb1_expr <= 0]))

TC71 <- SetIdent(TC71, cells=Cells(ZEB1_positive), value='ZEB1+')
TC71 <- SetIdent(TC71, cells=Cells(ZEB1_negative), value='ZEB1-')
TC71$zeb1_expr <- Idents(TC71)


# Fix object annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(TC71[["ATAC"]]) <- annotations


# Identify differentially expressed genes
DefaultAssay(TC71) <- 'RNA'
Idents(TC71) <- "condition"
tgfb1_de_genes <- FindMarkers(TC71, ident.1 = 'TGFB1', ident.2 = 'vehicle')
tgfb1_upreg <- tgfb1_de_genes[tgfb1_de_genes$avg_log2FC > 0 & tgfb1_de_genes$p_val_adj < 0.05, "X"]


# Identify differentially accessible peaks
DefaultAssay(TC71) <- 'ATAC'
Idents(TC71) <- 'condition'
tgfb1_da_peaks <- FindMarkers(
  object = TC71, ident.1 = "TGFB1", ident.2 = "vehicle", test.use = 'wilcox', min.pct = 0.1)

## Record significant peaks, label them with proximal genes
open_tgfb1 <- rownames(tgfb1_da_peaks[tgfb1_da_peaks$avg_log2FC > 0 & tgfb1_da_peaks$p_val_adj < 0.05, ])
closest_genes_tgfb1 <- ClosestFeature(TC71, regions = open_tgfb1)
tgfb1_accessible_genes <- closest_genes_tgfb1$gene_name


# Correlate peak accessibility with gene expression DIRECTLY -- intensive
TC71 <- LinkPeaks(object = TC71, peak.assay = "ATAC", expression.assay = "RNA")
TC71_links <-  Links(TC71)
write.csv(TC71_links, "../results/TC71_link_ranges.csv")


# Identify motifs enriched in differentially accessible ranges/peaks
upreg_pks <- tgfb1_da_peaks[tgfb1_da_peaks$avg_log2FC > 0 & tgfb1_da_peaks$p_val_adj < 0.05, ]
valid_peaks <- intersect(as.character(rownames(upreg_pks)), rownames(TC71[["ATAC"]]))
DefaultAssay(TC71) <- "ATAC"
enriched_motifs <- FindMotifs(TC71, features = valid_peaks)
write.csv(enriched_motifs, "../results/TC71_enriched_motifs.csv")

## Save progress
saveRDS(TC71, "../working_data/TC71_processed.rds")


# load data back in
# TC71 <- readRDS("../working_data/TC71_processed.rds")
# DefaultAssay(TC71) <- "ATAC"


# Calculate motif activity per cell -- very intensive
register(SerialParam())
TC71 <- RunChromVAR(TC71, genome = BSgenome.Hsapiens.UCSC.hg38)


# If ChromVAR works, save again
saveRDS(TC71, "../working_data/TC71_processed.rds")

