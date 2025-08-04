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

vehicle_treated <- readRDS("../working_data/vehicle_treated_processed.rds")
DefaultAssay(vehicle_treated) <- "ATAC"

pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
vehicle_treated <- AddMotifs(vehicle_treated, BSgenome.Hsapiens.UCSC.hg38, pfm) # add motif set, takes 20+ min

register(SerialParam()) # calibrate parallelization
vehicle_treated <- RunChromVAR(vehicle_treated, genome = BSgenome.Hsapiens.UCSC.hg38, new.assay.name = "JASPAR")

vehicle_treated@assays[["ATAC"]]@motifs <- readRDS("../working_data/vehicle_cisbp_motifs.rds")
vehicle_treated <- RunChromVAR(vehicle_treated, genome = BSgenome.Hsapiens.UCSC.hg38, new.assay.name = "cisBP")

saveRDS(vehicle_treated, "../working_data/vehicle_treated_processed_2.rds")
