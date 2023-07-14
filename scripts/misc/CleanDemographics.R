### LAST VERSUON UPDATE 2 MAY 2023 (v1.0)
### THIS SCRIPT CLEANS THE OUTPUT FROM `GetDemographics.R` INTO THE RELEVANT TABLES

print(paste0("GETTING DEMOGRAPHICS TABLE FOR THE ", COHORT, " COHORT"))

setwd("H:/Projects/MTX_GWAS/")

suppressMessages(suppressWarnings(library(dplyr))); options(dplyr.summarise.inform = FALSE)    #QUIETS `dplyr` SOMEWHAT
suppressMessages(suppressWarnings(library(haven)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(openxlsx)))

demo_df <- read_tsv("TMP/tmp-6/demo_df.txt", show_col_types = F)

### 6.1.6. LOAD KEYS AND SUBSET INTO TARGET COHORT

if(COHORT != "VALID") fam <- read.table(str_c("TMP/CACHE/", COHORT, "_QC.fam"), header = F) %>% distinct(V2) %>% select(IID = V2) %>% as_tibble()
if(COHORT == "VALID") fam <- read.table("TMP/CACHE/VALID_QC.fam", header = F) %>% distinct(V2) %>% mutate(IID = str_c("FAM001_", V2)) %>% as_tibble()
eira_key <- read_tsv("data/raw/Genotyped_EIRA_from_Leonid_20190823_IDconv.txt", show_col_types = F) %>% select(EIRA_id = EIRA, GWAS = SMP.number)
srqb_key <- read_sas("data/raw/key_20230126.sas7bdat") %>% mutate(SRQ_id = as.character(SRQ_id), BARCODE = as.character(barcode_deCODE_b1)) %>% select(SRQ_id, GWAS = BARCODE)
valid_key <- read_sas("data/raw/srqb_mtx_allra_b2_b3.sas7bdat") %>% select(SRQ_id, Barcode_deCODE_b2)
KEY <- read_tsv("data/raw/KEY_pid-GWAS.txt", show_col_types = F)

KEY_FULL <- KEY %>% left_join(eira_key %>% select(COHORT_ID = EIRA_id, IID = GWAS) %>%
                                bind_rows(srqb_key %>% select(COHORT_ID = SRQ_id, IID = GWAS)) %>%
                                bind_rows(valid_key %>% mutate(COHORT_ID = SRQ_id %>% as.character(), IID = Barcode_deCODE_b2 %>% as.character()) %>% select(COHORT_ID, IID)), by = "COHORT_ID") %>% distinct(pid, IID) %>%
  mutate(IID = str_c("FAM001_", IID))
fam_pid <- fam %>% left_join(KEY_FULL, by = "IID") %>% distinct()
demo_df_sub <- demo_df %>% filter(pid %in% fam_pid$pid)

rm(demo_df); rm(fam); rm(eira_key); rm(srqb_key); rm(valid_key); rm(KEY); rm(KEY_FULL); rm(fam_pid)

### 6.1.7. EXTRACT THE SUMMARY STATISTICS

DEMO_LIST <- list()

for(OUTCOME in c("persistence_d365", "persistence_d1096")){
  
  demo_df_pheno <- demo_df_sub %>% filter(!is.na(get(OUTCOME)))
  
  counts <- c("SEX", "SEROPOS", "EDU_CHAR", "DOSECAT", "ORAL", "ANY_FOLIC_ACID", "ANY_PRED", "RDCI")
  means <- c("AGE", "svullna_leder", "omma_leder", "sr", "crp", "patientens_globala", "das28", "das28CRP")
  
  counts_list <- list()
  means_list <- list()
  na_means_list <- list()
  
  for(i in 1:length(counts)){counts_list[[i]] <- demo_df_pheno %>% group_by(get(OUTCOME), get(counts[i])) %>% summarise(N = n()) %>% ungroup() %>% mutate(TYPE = counts[i])}
  for(i in 1:length(means)){means_list[[i]] <- demo_df_pheno %>% filter(!is.na(get(means[i]))) %>% group_by(get(OUTCOME)) %>% summarise(MU = mean(get(means[i])), SIGMA = sd(get(means[i]))) %>% ungroup() %>% mutate(TYPE = means[i])}
  for(i in 1:length(means)){na_means_list[[i]] <- demo_df_pheno %>% filter(is.na(get(means[i]))) %>% group_by(get(OUTCOME)) %>% summarise(N_NA = n()) %>% ungroup() %>% mutate(TYPE = means[i])}
  
  counts_df <- counts_list %>% bind_rows() %>% group_by(`get(OUTCOME)`, TYPE) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>%
    pivot_wider(names_from = `get(OUTCOME)`, values_from = c(N, N_TOT)) %>%
    mutate(V1 = paste0(N_0, " (", round(100 * N_0 / N_TOT_0, 0), "%)"), V2 = paste0(N_1, " (", round(100 * N_1 / N_TOT_1, 0), "%)")) %>%
    mutate(TYPE = paste0(TYPE, "_", `get(counts[i])`)) %>%
    select(TYPE, V1, V2)
  
  means_df <- means_list %>% bind_rows() %>% group_by(`get(OUTCOME)`, TYPE) %>% ungroup() %>% 
    mutate(MU = ifelse(TYPE %in% c("AGE", "svullna_leder", "omma_leder"), round(MU, 0), MU), SIGMA = ifelse(TYPE %in% c("AGE", "svullna_leder", "omma_leder"), round(SIGMA, 0), SIGMA)) %>%
    pivot_wider(names_from = `get(OUTCOME)`, values_from = c(MU, SIGMA)) %>%
    mutate(V1 = paste0(round(MU_0, 2), " (", round(SIGMA_0, 2), ")"), V2 = paste0(round(MU_1, 2), " (", round(SIGMA_1, 2), ")")) %>%
    select(TYPE, V1, V2)
  
  na_means_df <- na_means_list %>% bind_rows() %>% group_by(`get(OUTCOME)`, TYPE) %>% 
    pivot_wider(names_from = `get(OUTCOME)`, values_from = N_NA) %>% 
    mutate(V1 = as.character(`0`), V2 = as.character(`1`), TYPE = paste0(TYPE, "_NA")) %>%
    select(TYPE, V1, V2)
  
  medians_df <- demo_df_pheno %>% group_by(get(OUTCOME)) %>% summarise(MEDIAN_YEAR = median(YEAR), LQ_YEAR = quantile(YEAR, 0.25) %>% str_extract("\\d+") %>% str_extract("\\d\\d$"), UQ_YEAR = quantile(YEAR, 0.75) %>% str_extract("\\d+") %>% str_extract("\\d\\d$")) %>% mutate(YEAR = str_c(MEDIAN_YEAR, " (", LQ_YEAR, "-", UQ_YEAR, ")"))
  
  res <- counts_df %>% bind_rows(means_df) %>% bind_rows(na_means_df) %>% filter(TYPE != "ANY_FOLIC_ACID_0", TYPE != "ANY_PRED_0", TYPE != "SEX_0", TYPE != "SEROPOS_0", TYPE != "ORAL_0") %>% bind_rows(data.frame(TYPE = "YEAR", V1 = medians_df$YEAR[2], V2 = medians_df$YEAR[1]))
  colnames(res) <- c("VAR", paste0("NON_", OUTCOME), OUTCOME)
  DEMO_LIST[[str_which(c("persistence_d365", "persistence_d1096"), OUTCOME)]] <- res
  
}

rm(OUTCOME); rm(demo_df_pheno); rm(counts); rm(means); rm(counts_list); rm(means_list); rm(na_means_list); rm(counts_df); rm(means_df); rm(na_means_df); rm(medians_df); rm(res)

### 6.1.8. CLEAN UP THE OUTPUT TABLE

FINAL <- DEMO_LIST[[1]] %>% left_join(DEMO_LIST[[2]], by = "VAR") %>% slice(1:2, 22, 37, 4:6, 18:20, 23:24, 26, 25, 27, 28:29, 8:10, 12, 14, 16, 3, 7, 11, 13, 15, 17, 21, 30:36)
FINAL <- c("N", nrow(demo_df_sub %>% filter(persistence_d365 == 0)), nrow(demo_df_sub %>% filter(persistence_d365 == 1)), nrow(demo_df_sub %>% filter(persistence_d1096 == 0)), nrow(demo_df_sub %>% filter(persistence_d1096 == 1))) %>% rbind(FINAL) %>% select(VAR, persistence_d365, NON_persistence_d365, persistence_d1096, NON_persistence_d1096)

if(COHORT == "RA"){
  
  counts <- c("SEX", "SEROPOS", "EDU_CHAR", "DOSECAT", "ORAL", "ANY_FOLIC_ACID", "ANY_PRED", "RDCI")
  means <- c("AGE", "svullna_leder", "omma_leder", "sr", "crp", "patientens_globala", "das28", "das28CRP")
  
  counts_list <- list()
  means_list <- list()
  na_means_list <- list()
  
  for(i in 1:length(counts)){counts_list[[i]] <- demo_df_sub %>% group_by(get(counts[i])) %>% summarise(N = n()) %>% ungroup() %>% mutate(TYPE = counts[i])}
  for(i in 1:length(means)){means_list[[i]] <- demo_df_sub %>% filter(!is.na(get(means[i]))) %>% summarise(MU = mean(get(means[i])), SIGMA = sd(get(means[i]))) %>% mutate(TYPE = means[i])}
  for(i in 1:length(means)){na_means_list[[i]] <- demo_df_sub %>% filter(is.na(get(means[i]))) %>% summarise(N_NA = n()) %>% mutate(TYPE = means[i])}
  
  counts_df <- counts_list %>% bind_rows() %>% group_by(TYPE) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>% mutate(V1 = str_c(N, " (", round(100 * N / N_TOT), "%)")) %>% mutate(TYPE = paste0(TYPE, "_", `get(counts[i])`)) %>% select(TYPE, V1)
  means_df <- means_list %>% bind_rows() %>% mutate(V1 = str_c(round(MU, 2), " (", round(SIGMA, 2), ")")) %>% select(TYPE, V1)
  na_means_df <- na_means_list %>% bind_rows() %>% mutate(TYPE = str_c(TYPE, "_NA"), V1 = as.character(N_NA)) %>% select(TYPE, V1)
  medians_df <- demo_df_sub %>% summarise(MEDIAN_YEAR = median(YEAR), LQ_YEAR = quantile(YEAR, 0.25) %>% str_extract("\\d\\d$"), UQ_YEAR = quantile(YEAR, 0.75) %>% str_extract("\\d\\d$")) %>% mutate(YEAR = str_c(MEDIAN_YEAR, " (", LQ_YEAR, "-", UQ_YEAR, ")"))
  
  res <- counts_df %>% bind_rows(means_df) %>% bind_rows(na_means_df) %>% filter(TYPE != "ANY_FOLIC_ACID_0", TYPE != "ANY_PRED_0", TYPE != "SEX_0", TYPE != "SEROPOS_0", TYPE != "ORAL_0") %>% bind_rows(data.frame(TYPE = "YEAR", V1 = medians_df$YEAR[1]))
  colnames(res) <- c("VAR", "OVERALL")
  res <- c("N", nrow(demo_df_sub)) %>% rbind(res)
  FINAL <- FINAL %>% left_join(res, by = "VAR") %>% select(VAR, OVERALL, everything())
  
  FILEPATH <- "data/res/clean/Table 1.xlsx"
  
}else if(COHORT == "SPOS"){FILEPATH <- "TMP/tmp-6/SPOS.xlsx"
}else if(COHORT == "SNEG"){FILEPATH <- "TMP/tmp-6/SNEG.xlsx"
}else if(COHORT == "SENS"){FILEPATH <- "data/res/clean/Table S6.xlsx"
}else if(COHORT == "VALID"){FILEPATH <- "data/res/clean/Table S7.xlsx"}

write.xlsx(FINAL, FILEPATH)
rm(list = ls())
