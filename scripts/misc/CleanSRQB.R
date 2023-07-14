### LAST VERSION UPDATE 27 APRIL 2023 (v2.1) - UPDATED TO BE COMPARABLE TO THE `EIRA-SRQB.txt` FILE PRODUCED FOR THE PRIMARY COHORT
### THIS SCRIPT SUBSETS THE SRQB DATA TO ONLY THOSE REQUIRED FOR MY MTX PROJECT

.libPaths("/home2/genetics/antobe/software/RLibrary/")

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(haven)))

TARGET_POP <- read_sas("data/raw/srqb_mtx_allra_b2_b3.sas7bdat")
TARGET_POP_CLEAN <- TARGET_POP %>% select(pid, Barcode_deCODE_b2, sex, age, persistence_d365, persistence_d1096) %>% mutate(FID = "FAM001") %>% select(pid, FID, IID = Barcode_deCODE_b2, everything()) %>% arrange(FID, IID)

write.table(TARGET_POP_CLEAN, "TMP/CACHE/SRQB.txt", row.names = F, col.names = T, quote = F, sep = "\t")
