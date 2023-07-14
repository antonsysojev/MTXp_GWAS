### LAST VERSION UPDATE APRIL 27 2023 (v2.0) - ADDED CODE FOR HANDLING THE VALIDATION DATA
### THIS SCRIPT EXTRACTS THE GWAS COVARIATES FOR AN INPUT COHORT FILENAME

.libPaths("/home2/genetics/antobe/software/RLibrary")

suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(stringr)))

ARGS <- arg_parser('EMPTY FOR NOW') %>% add_argument('COHORT', help = 'ID OF THE TARGET COHORT FOR WHICH TO GET THE COVARIATES', type = 'character') %>% parse_args()
FILEPATH_PREFIX <- str_c("TMP/CACHE/", ARGS$COHORT, "_QC")
FAM <- read_tsv(str_c(FILEPATH_PREFIX, ".fam"), col_names = c("FID", "IID", "V3", "V4", "SEX", "PHENO"), show_col_types = F)

if(ARGS$COHORT != "VALID") COV <- read_tsv("TMP/CACHE/EIRA-SRQB.txt", show_col_types = F)
if(ARGS$COHORT == "VALID") COV <- read_tsv("TMP/CACHE/SRQB.txt", show_col_types = F)
PCA <- read_tsv(str_c(FILEPATH_PREFIX, ".eigenvec"), show_col_types = F)
FAM_COV_PCA <- FAM %>% inner_join(COV %>% select(-FID), by = "IID") %>% inner_join(PCA, by = "IID")

if(ARGS$COHORT != "VALID") COV.txt <- FAM_COV_PCA %>% select(-V3, -V4, -pid, -persistence_d365, -persistence_d1096, -RA, -SPOS, -SNEG, -SENS, -SEX, -PHENO, -`#FID`) %>% distinct()
if(ARGS$COHORT == "VALID") COV.txt <- FAM_COV_PCA %>% select(FID, IID, sex, age, starts_with("PC")) %>% distinct()
write.table(COV.txt, str_c("TMP/tmp-2/", ARGS$COHORT, "_COV.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

PHENO.txt <- FAM_COV_PCA %>% select(FID, IID, persistence_d365, persistence_d1096) %>%
	mutate(persistence_d365 = persistence_d365 + 1, persistence_d1096 = persistence_d1096 + 1) %>%
	mutate(persistence_d365 = if_else(is.na(persistence_d365), -9, persistence_d365), persistence_d1096 = if_else(is.na(persistence_d1096), -9, persistence_d1096)) %>% distinct()
write.table(PHENO.txt, str_c("TMP/CACHE/", ARGS$COHORT, "_PHENO.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

### TO DO:
# 1.1. Actually, I reckon that the covariates handling of the VALID and the other inputs can be handled the same way, I just didn't have the data available to check it when I wrote the code and thus left it untouched.
### NOTES:
# 2.1. There was a discrete bug within this script that led to duplicates within data. 
#	When creating the `TMP/CACHE/EIRA-SRQB.txt` file, I ignored the fact that duplicates were abundant. 
#	This was fine since any and all duplicates were dealth with during the QC. 
#	However, duplicates remained in the data, and the file used for `COV` here contained several. 
#	Since I didn't get rid of this PRIOR to joining them onto the .fam-file, they multiply lines within the `FAM_COV_PCA` table. 
#	I added `distinct()` to the final files which SHOULD work, but one should clean the files when loading them into R, by selecting only the variables of interest and doing a quick `distinct()` call at this point already. 
#	The current approach works for now and remains valid.
