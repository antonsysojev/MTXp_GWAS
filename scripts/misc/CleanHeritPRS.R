### LAST VERSION UPDATE 27 APRIL 2023 (v1.1) - ADDED VALIDATION INTO SCRIPT AND ADJUSTED OUTPUT FILES
### THIS SCRIPT CLEANS AND OUTPUTS PRESENTABLE RESULTS FROM THE HERITABILITY ESTIMATION AND PRS TESTING

.libPaths("/home2/genetics/antobe/software/RLibrary/")

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(argparser)))
suppressMessages(suppressWarnings(library(readr)))

ARGS <- arg_parser('EMPTY FN') %>% add_argument('COHORT', help = 'COHORT OF THE TARGET GWAS FOR WHICH TO GET THE COVARIATES', type = 'character') %>% parse_args()
COHORT <- ARGS$COHORT

### 5.2.1. CLEAN HERITABILITY OUTPUT

if(COHORT != "VALID"){    #NO HERITABILITY TO ACQUIRE FOR THE VALIDATION DATA
	
	herit_P1YR <- readLines(str_c("data/res/", COHORT, "_P1YR.hsq")) %>% str_extract_all("V\\(G\\)/Vp_L.+") %>% unlist() %>% str_extract("\\d.+")
	herit_P3YR <- readLines(str_c("data/res/", COHORT, "_P3YR.hsq")) %>% str_extract_all("V\\(G\\)/Vp_L.+") %>% unlist() %>% str_extract("\\d.+")
	herit_P1YR_ADJ <- readLines(str_c("data/res/", COHORT, "_P1YR_adj.hsq")) %>% str_extract_all("V\\(G\\)/Vp_L.+") %>% unlist() %>% str_extract("\\d.+")
	herit_P3YR_ADJ <- readLines(str_c("data/res/", COHORT, "_P3YR_adj.hsq")) %>% str_extract_all("V\\(G\\)/Vp_L.+") %>% unlist() %>% str_extract("\\d.+")

	herit <- data.frame(RAW = c(herit_P1YR, herit_P1YR_ADJ, herit_P3YR, herit_P3YR_ADJ), ID = c("P1YR", "P1YR_ADJ", "P3YR", "P3YR_ADJ")) %>% as_tibble() %>%
	  mutate(h2 = str_extract(RAW, "^.+?(?=\t)") %>% as.numeric(), ERR = str_extract(RAW, "\t.+$") %>% str_remove("\t") %>% as.numeric()) %>%
	  mutate(CI_L = h2 - ERR * qnorm(1 - 0.05 / 2), CI_U = h2 + ERR * qnorm(1 - 0.05 / 2), h2 = h2 %>% round(digits = 2) %>% format(nsmall = 2), standard_error = ERR %>% round(digits = 2) %>% format(nsmall = 2)) %>%
	  mutate(CI = str_c(CI_L %>% round(digits = 2) %>% format(nsmall = 2), CI_U %>% round(digits = 2) %>% format(nsmall = 2), sep = "-")) %>%
	  mutate(CI = str_trim(CI)) %>%
	  select(ID, h2, h2_standard_error = standard_error, h2_CI = CI) %>%
	  mutate(h2 = str_c(h2, " (", h2_CI, ")")) %>% select(ID, h2, h2_standard_error)

}

### 5.2.2. CLEAN PRS OUTPUT

PRS_df <- read_tsv(str_c("data/res/", COHORT, "_PRS.txt"), show_col_types = F)
PRS <- PRS_df %>% mutate(RR = exp(BETA), CI_L = exp(BETA - SE * qnorm(1 - 0.05 / 2)), CI_U = exp(BETA + SE * qnorm(1 - 0.05 / 2))) %>%
  mutate(RR = RR %>% round(digits = 2) %>% format(nsmall = 2), CI = str_c(CI_L %>% round(digits = 2) %>% format(nsmall = 2), CI_U %>% round(digits = 2) %>% format(nsmall = 2), sep = "-")) %>%
  mutate(RR = str_c(RR, " (", CI, ")")) %>%
  mutate(TrendP = TrendP %>% round(digits = 4) %>% format(nsmall = 2)) %>% mutate(TrendP = ifelse(str_detect(TrendP, "NA"), "", TrendP)) %>% 
  select(ID, RR, RR_N = N, RR_P = P, PForTrend = TrendP)

### 5.2.3. OUTPUT FINAL FILES

if(COHORT != "VALID"){
	herit_PRS <- herit %>% inner_join(PRS, by = "ID") %>% select(-PForTrend) %>% mutate(COHORT = COHORT)
	PRS_quantiles <- PRS %>% filter(str_detect(ID, "_0\\d$"))
	write.table(PRS_quantiles, str_c("TMP/tmp-5/", COHORT, "_TableS5.txt"), col.names = T, row.names = F, quote = T, sep = "\t")
}else{ 
	herit_PRS <- PRS %>% slice(1:4) %>% select(-PForTrend)
}

write.table(herit_PRS, str_c("TMP/tmp-5/", COHORT, "_heritPRS.txt"), col.names = T, row.names = F, quote = T, sep = "\t")

### TO DO:
### NOTES:
