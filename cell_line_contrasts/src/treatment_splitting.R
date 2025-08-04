library(Seurat)
library(scCustomize)
library(dplyr)
library(Signac)

# source dataset -- in separate directory for safekeeping
dataset <- readRDS("../../VTP_SC_seq/ews_tgfb_seurat.rds")

Idents(dataset) <- "condition"

t1_treated <- subset(dataset, idents = c("TGFB1"))
vehicle <- subset(dataset, idents = c("vehicle"))

saveRDS(t1_treated, "../working_data/tgfb1_treated.rds")
saveRDS(vehicle, "../working_data/vehicle_treated.rds")
