library(Signac)
library(Seurat)


CHLA10 <- readRDS("/data/hps/assoc/private/ewing_beti/user/jvale8/zeb_singleCell/working_data/CHLA10_subset.rds")

CHLA10 <- LinkPeaks(CHLA10, peak.assay = "ATAC", expression.assay = "RNA")

saveRDS(CHLA10, "/data/hps/assoc/private/ewing_beti/user/jvale8/zeb_singleCell/working_data/CHLA10_linked.rds")