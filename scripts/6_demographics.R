### LAST VERSION UPDATE 2 MAY 2023 (v1.2) - BIG REVISION OF THE STRUCTURE, NOW WITH LESS UNNECESSARY COMPUTATIONS
### THIS SCRIPT LAUNCHES MULTIPLE SUBSCRIPTS FOR WHICH TO EXTRACT COHORT DEMOGRAPHIC TABLES FOR THE STUDY MANUSCRIPT

setwd("H:/Projects/MTX_GWAS/")

### 6.1. EXTRACT DEMOGRAPHICS FOR TABLES 1, 2, S6 AND S7

source("scripts/misc/GetDemographics.R")

for(COHORT in c("RA", "SPOS", "SNEG", "SENS", "VALID")){source("scripts/misc/CleanDemographics.R")}

SPOS <- read.xlsx("TMP/tmp-6/SPOS.xlsx"); SNEG <- read.xlsx("TMP/tmp-6/SNEG.xlsx")
SPOS %>% left_join(SNEG, by = "VAR", suffix = c("_SPOS", "_SNEG")) %>% write.xlsx("data/res/clean/Table 2.xlsx")
rm(SPOS); rm(SNEG)

### 6.2. EXTRACT REASON FOR DISCONTINUATION FOR TABLE S2

for(COHORT in c("RA", "SPOS", "SNEG")){source("scripts/misc_private/Discontinuation.R")}

read_tsv("TMP/tmp-6/DISCONTINUE_RA.txt", show_col_types = F) %>% mutate(TYPE = "RA") %>% 
  bind_rows(read_tsv("TMP/tmp-6/DISCONTINUE_SPOS.txt", show_col_types = F) %>% mutate(TYPE = "SPOS")) %>%
  bind_rows(read_tsv("TMP/tmp-6/DISCONTINUE_SNEG.txt", show_col_types = F) %>% mutate(TYPE = "SNEG")) %>%
  pivot_longer(c("non_persistent_365", "non_persistent_1096"), names_to = "OUTCOME") %>% 
  pivot_wider(names_from = c("TYPE", "OUTCOME"), values_from = value) %>%
  select(REASON = reason_CLEAN, RA_365 = RA_non_persistent_365, SPOS_365 = SPOS_non_persistent_365, SNEG_365 = SNEG_non_persistent_365, RA_1096 = RA_non_persistent_1096, SPOS_1096 = SPOS_non_persistent_1096, SNEG_1096 = SNEG_non_persistent_1096) %>%
  slice(6, 4, 1, 3, 2, 5) %>%
  write.xlsx("data/res/clean/Table S2.xlsx")
