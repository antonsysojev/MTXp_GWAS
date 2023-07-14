### LAST VERSION UPDATE 4 MAY 2023 (v1.0)
### THIS SCRIPT UPDATES THE NAMES OF GWAS SNP WHILE AVOIDING DUPLICATED IDS

.libPaths("/home2/genetics/antobe/software/RLibrary")
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(argparser)))

args <- arg_parser('EMPTY FN') %>% add_argument('GWAS', help = 'NAME OF THE GWAS FOR WHICH TO CLEAN THE VARIANT NAMES', type = 'character') %>% parse_args()
GWAS <- args$GWAS

GWAS_df <- read_tsv(paste0("data/res/", GWAS, ".PHENO1.glm.logistic.hybrid"), show_col_types = F)
GWAS_df_clean <- GWAS_df %>% mutate(CHRPOS = paste0(`#CHROM`, ":", POS)) %>% 
	group_by(CHRPOS) %>% mutate(N = n()) %>% ungroup() %>%
	mutate(ID = ifelse(N == 1, CHRPOS, paste0(CHRPOS, ":", REF, ":", ALT))) %>%
	select(-CHRPOS)

write.table(GWAS_df_clean, paste0("TMP/tmp-2/", GWAS, "_IDFIX.PHENO1.glm.logistic.hybrid"), col.names = T, row.names = F, quote = F, sep = "\t")
