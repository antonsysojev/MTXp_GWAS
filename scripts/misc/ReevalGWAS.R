### LAST VERSION UPDATE 27 APRIL 2023 (v1.0)
### THIS SCRIPT PERFORMS THE RE EVALUATION OF SUGGESTIVE SNP ASSOCIATIONS IN THE PRIMARY COHORT

.libPaths("/home2/genetics/antobe/software/RLibrary/")

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(argparser)))

ARGS <- arg_parser('EMPTY FN') %>% add_argument('OUTCOME', help = 'OUTCOME OF THE TARGET GWAS TO RE EVALUATE - EITHER P1YR OR P3YR', type = 'character') %>% parse_args()
OUTCOME <- ARGS$OUTCOME
GWAS <- read_tsv(str_c("data/res/RA_", OUTCOME, ".PHENO1.glm.logistic.hybrid"), show_col_types = F)
SENS <- read_tsv(str_c("data/res/SENS_", OUTCOME, ".PHENO1.glm.logistic.hybrid"), show_col_types = F)

SUGGESTIVE <- GWAS %>% filter(P < 5e-5) %>% mutate(CHRPOS = str_extract(ID, "^\\d+\\:\\d+"))
BONFERRONI <- 0.05 / nrow(SUGGESTIVE)

TARGET_SNPS <- SENS %>% mutate(CHRPOS = str_extract(ID, "^\\d+\\:\\d+")) %>% filter(CHRPOS %in% SUGGESTIVE$CHRPOS)
SIGNIFICANT <- TARGET_SNPS %>% filter(P < BONFERRONI)

print(str_c("FOR ", OUTCOME, ", ", nrow(SUGGESTIVE), " SNPS WERE SUGGESTIVELY ASSOCIATED"))
print(str_c("OF THESE, ", nrow(TARGET_SNPS), " WERE AVAILABLE WITHIN THE SENSITIVITY ANALYSIS DATA"))
print(str_c("OF THESE, ", nrow(SIGNIFICANT), " WERE SIGNIFICANTLY ASSOCIATED WITH ", OUTCOME, " AFTER BONFERRONI CORRECTION"))

### TO DO:
### NOTES:

