#!/user/bin/env Rscript
### LAST VERSION UPDATE AUG 31 2023 (v2.1.2) - FIXED AN ISSUE WITH THE GWAS CATALOG VALIDATION, NOW CUTS ALL SNPS WITH A NON-ALPHABETIC ALLELE
### THIS SCRIPT CLEANS THE GWAS DATA FOR ONLINE PUBLICATION

.libPaths("/home2/genetics/antobe/software/RLibrary/")

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(argparser)))

args <- arg_parser('EMPTY FN') %>% add_argument('GWAS', help = 'NAME OF THE GWAS FOR WHICH TO CLEAN THE CONTENTS', type = 'character') %>% parse_args()
GWAS <- args$GWAS

if(str_detect(GWAS, "META")){
	
	GWAS_df <- read_tsv(str_c("TMP/tmp-7/", GWAS, ".tsv"), show_col_types = F)

	GWAS_df_clean <- GWAS_df %>% mutate(chromosome = str_extract(MarkerName, "^\\d+"), base_pair_location = str_extract(MarkerName, "\\d+$")) %>%
	    mutate(effect_allele = str_to_upper(Allele2), other_allele = str_to_upper(Allele1)) %>%
	    mutate(Z = Zscore, P = `P-value`) %>%
	    mutate(variant_id = MarkerName, n = Weight) %>%
	    select(chromosome, base_pair_location, effect_allele, other_allele, Z, P, variant_id, n) %>% arrange(chromosome, base_pair_location)

}

if(!str_detect(GWAS, "META")){

	GWAS_df <- read_tsv(str_c("data/res/", GWAS, ".PHENO1.glm.logistic.hybrid"), show_col_types = F)
	FREQ <- read_tsv(str_c("TMP/tmp-7/", GWAS, "_CTRL.afreq"), show_col_types = F) %>% mutate(variant_id = str_replace_all(ID, "\\:", "_")) %>% select(variant_id, ALT_FREQS)
	SNP <- read_tsv("/home2/genetics/antobe/data/EIRA-SRQB/SNP-ID_linkage.txt", col_names = c("ID", "RSNO"), show_col_types = F) %>% mutate(variant_id = str_replace_all(ID, "\\:", "_")) %>% select(variant_id, RSNO)


	GWAS_df_clean <- GWAS_df %>% mutate(chromosome = `#CHROM`, base_pair_location = POS) %>%
	  mutate(effect_allele = ALT, other_allele = REF) %>%
	  mutate(beta = log(OR), standard_error = `LOG(OR)_SE`, p_value = P) %>%
	  mutate(variant_id = str_replace_all(ID, "\\:", "_"), n = OBS_CT) %>%
	  inner_join(FREQ, by = "variant_id") %>% 
	  inner_join(SNP, by = "variant_id") %>%
	  select(chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency = ALT_FREQS, p_value, variant_id, n, rs_id = RSNO)

}

GWAS_df_clean <- GWAS_df_clean %>% filter(!str_detect(effect_allele, "[^A-Z]")) %>% filter(!str_detect(other_allele, "[^A-Z]")) %>% mutate(base_pair_location = as.integer(base_pair_location))    #AD-HOC SOLUTION TO MAKE THE FILE PASS THROUGH THE `gwas-ssf validate` ISSUE
write.table(GWAS_df_clean, str_c("TMP/tmp-7/", GWAS, ".tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

### TO DO:
# 1.1. Currently does not get the frequency or the rs-numbers due to laziness. Can be updated to add this later.
# 1.2. Order of the meta analysis data is off, since it is sorted by character and not be numericals. Should fix this later.
### NOTES:
