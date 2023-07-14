### LAST VERSION UPDATE 24 APRIL 2023 (v1.1) - ADDED CODE FOR PRS ANALYSIS ACROSS QUANTILES
### THIS SCRIPT PERFORMS THE REGRESSION TESTING OF THE EFFECT OF THE RA PRS ON OUR PERSISTENCE PHENOTYPES

.libPaths("/home2/genetics/antobe/software/RLibrary")

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(bigsnpr)))
suppressMessages(suppressWarnings(library(argparser)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(DescTools)))

ARGS <- arg_parser('EMPTY FN') %>% add_argument('COHORT', help = 'ID OF THE TARGET COHORT FOR WHICH TO GET THE COVARIATES', type = 'character') %>% parse_args()
COHORT <- ARGS$COHORT
FILEPATH <- str_c("TMP/CACHE/", COHORT)

PHENO_df <- read_tsv(str_c(FILEPATH, "_PHENO.txt"), show_col_types = F) %>% mutate(ID = str_c(FID, "_", IID)) %>% select(ID, everything()) %>% select(-FID, -IID)
PRS_df <- read_tsv(str_c(FILEPATH, "_PRS.txt"), show_col_types = F) %>% mutate(ID = str_c(FID, "_", IID)) %>% select(ID, everything()) %>% select(-FID, -IID)
COV_df <- read_tsv(str_c(FILEPATH, "_COVSUB.txt"), show_col_types = F) %>% mutate(ID = str_c(FID, "_", IID)) %>% select(ID, everything()) %>% select(-FID, -IID)

### 4.2.1. CLEAN INPUT DATA

PRS_quantiles <- qnorm(c(0, 0.2, 0.4, 0.6, 0.8, 1))
df <- PHENO_df %>% inner_join(PRS_df, by = "ID") %>% inner_join(COV_df, by = "ID") %>%
	mutate(MU = mean(RA_PRS), SIGMA = sd(RA_PRS)) %>% mutate(RA_PRS_NORM = (RA_PRS - MU) / SIGMA) %>% select(-MU, -SIGMA) %>%
	mutate(RA_PRS_NORM_CAT = cut(RA_PRS_NORM, PRS_quantiles, 0:4)) %>%
	select(ID, RA_PRS, RA_PRS_NORM, RA_PRS_NORM_CAT, everything())

df_P1YR <- df %>% select(-ID, -persistence_d1096) %>% filter(persistence_d365 != -9) %>% mutate(persistence_d365 = persistence_d365 - 1)
df_P3YR <- df %>% select(-ID, -persistence_d365) %>% filter(persistence_d1096 != -9) %>% mutate(persistence_d1096 = persistence_d1096 - 1)

### 4.2.2. PERFORM REGRESSIONS
### ### 4.2.2.1. PERFORM PRIMARY REGRESSIONS

glm_P1YR_crude <- glm(persistence_d365 ~ RA_PRS_NORM, data = df_P1YR, family = binomial(link = "log"))
glm_P3YR_crude <- glm(persistence_d1096 ~ RA_PRS_NORM, data = df_P3YR, family = binomial(link = "log"))

glm_P1YR_adj <- glm(persistence_d365 ~ RA_PRS_NORM + sex + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = df_P1YR, family = binomial(link = "log"), start = c(0, 0, 0, 0 - 1e-6, 0, 0, 0, 0, 0, 0))
glm_P3YR_adj <- glm(persistence_d1096 ~ RA_PRS_NORM + sex + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = df_P3YR, family = binomial(link = "log"), start = c(0, 0, 0, 0 - 1e-6, 0, 0, 0, 0, 0, 0))

### ### 4.2.2.2. PERFORM SECONDARY REGRESSIONS (CONTRASTING ACROSS QUANTILES)

df_P1YR_01 <- df_P1YR %>% filter(RA_PRS_NORM_CAT == 0 | RA_PRS_NORM_CAT == 1) %>% mutate(RA_PRS_NORM_CAT_01 = ifelse(RA_PRS_NORM_CAT == 0, 0, 1))
df_P1YR_02 <- df_P1YR %>% filter(RA_PRS_NORM_CAT == 0 | RA_PRS_NORM_CAT == 2) %>% mutate(RA_PRS_NORM_CAT_02 = ifelse(RA_PRS_NORM_CAT == 0, 0, 1))
df_P1YR_03 <- df_P1YR %>% filter(RA_PRS_NORM_CAT == 0 | RA_PRS_NORM_CAT == 3) %>% mutate(RA_PRS_NORM_CAT_03 = ifelse(RA_PRS_NORM_CAT == 0, 0, 1))
df_P1YR_04 <- df_P1YR %>% filter(RA_PRS_NORM_CAT == 0 | RA_PRS_NORM_CAT == 4) %>% mutate(RA_PRS_NORM_CAT_04 = ifelse(RA_PRS_NORM_CAT == 0, 0, 1))

df_P3YR_01 <- df_P3YR %>% filter(RA_PRS_NORM_CAT == 0 | RA_PRS_NORM_CAT == 1) %>% mutate(RA_PRS_NORM_CAT_01 = ifelse(RA_PRS_NORM_CAT == 0, 0, 1))
df_P3YR_02 <- df_P3YR %>% filter(RA_PRS_NORM_CAT == 0 | RA_PRS_NORM_CAT == 2) %>% mutate(RA_PRS_NORM_CAT_02 = ifelse(RA_PRS_NORM_CAT == 0, 0, 1))
df_P3YR_03 <- df_P3YR %>% filter(RA_PRS_NORM_CAT == 0 | RA_PRS_NORM_CAT == 3) %>% mutate(RA_PRS_NORM_CAT_03 = ifelse(RA_PRS_NORM_CAT == 0, 0, 1))
df_P3YR_04 <- df_P3YR %>% filter(RA_PRS_NORM_CAT == 0 | RA_PRS_NORM_CAT == 4) %>% mutate(RA_PRS_NORM_CAT_04 = ifelse(RA_PRS_NORM_CAT == 0, 0, 1))

glm_P1YR_01_crude <- glm(persistence_d365 ~ RA_PRS_NORM_CAT_01, data = df_P1YR_01, family = binomial(link = "log"))
glm_P1YR_02_crude <- glm(persistence_d365 ~ RA_PRS_NORM_CAT_02, data = df_P1YR_02, family = binomial(link = "log"))
glm_P1YR_03_crude <- glm(persistence_d365 ~ RA_PRS_NORM_CAT_03, data = df_P1YR_03, family = binomial(link = "log"))
glm_P1YR_04_crude <- glm(persistence_d365 ~ RA_PRS_NORM_CAT_04, data = df_P1YR_04, family = binomial(link = "log"))

glm_P3YR_01_crude <- glm(persistence_d1096 ~ RA_PRS_NORM_CAT_01, data = df_P3YR_01, family = binomial(link = "log"))
glm_P3YR_02_crude <- glm(persistence_d1096 ~ RA_PRS_NORM_CAT_02, data = df_P3YR_02, family = binomial(link = "log"))
glm_P3YR_03_crude <- glm(persistence_d1096 ~ RA_PRS_NORM_CAT_03, data = df_P3YR_03, family = binomial(link = "log"))
glm_P3YR_04_crude <- glm(persistence_d1096 ~ RA_PRS_NORM_CAT_04, data = df_P3YR_04, family = binomial(link = "log"))

### ### ### 4.2.2.2.1. TEST THE OUTPUT FROM REGRESSION ACROSS QUANTILES FOR A TREND

TrendTestMat_P1YR <- matrix((df_P1YR %>% group_by(RA_PRS_NORM_CAT, persistence_d365) %>% summarise(N = n()) %>% arrange(persistence_d365))$N, byrow = T, nrow = 2) %>% suppressMessages()
TrendTestMat_P3YR <- matrix((df_P3YR %>% group_by(RA_PRS_NORM_CAT, persistence_d1096) %>% summarise(N = n()) %>% arrange(persistence_d1096))$N, byrow = T, nrow = 2) %>% suppressMessages()

TrendP_P1YR <- CochranArmitageTest(TrendTestMat_P1YR)$p.value
TrendP_P3YR <- CochranArmitageTest(TrendTestMat_P3YR)$p.value

### 4.2.3. OUTPUT RESULTS

RES_df <- data.frame(ID = c("P1YR", "P3YR", "P1YR_ADJ", "P3YR_ADJ"), N = rep(c(nrow(df_P1YR), nrow(df_P3YR)), 2), BETA = NA, SE = NA, P = NA)
RES_df[1, 3:5] <- coef(glm_P1YR_crude %>% summary())[2, c(1, 2, 4)]
RES_df[2, 3:5] <- coef(glm_P3YR_crude %>% summary())[2, c(1, 2, 4)]
RES_df[3, 3:5] <- coef(glm_P1YR_adj %>% summary())[2, c(1, 2, 4)]
RES_df[4, 3:5] <- coef(glm_P3YR_adj %>% summary())[2, c(1, 2, 4)]
RES_df$TrendP <- NA

RES_df_quantiles <- data.frame(ID = str_c(rep(c("P1YR", "P3YR"), each = 4), rep(c("01", "02", "03", "04"), 2), sep = "_"), N = c(nrow(df_P1YR_01), nrow(df_P1YR_02), nrow(df_P1YR_03), nrow(df_P1YR_04), nrow(df_P3YR_01), nrow(df_P3YR_02), nrow(df_P3YR_03), nrow(df_P3YR_04)), BETA = NA, SE = NA, P = NA)
#c(rep("P1YR", 4), rep("P3YR", 4)), rep(c("01", "02", "03", "04"), 2), sep = "_"), 
RES_df_quantiles[1, 3:5] <- coef(glm_P1YR_01_crude %>% summary())[2, c(1, 2, 4)]
RES_df_quantiles[2, 3:5] <- coef(glm_P1YR_02_crude %>% summary())[2, c(1, 2, 4)]
RES_df_quantiles[3, 3:5] <- coef(glm_P1YR_03_crude %>% summary())[2, c(1, 2, 4)]
RES_df_quantiles[4, 3:5] <- coef(glm_P1YR_04_crude %>% summary())[2, c(1, 2, 4)]
RES_df_quantiles[5, 3:5] <- coef(glm_P3YR_01_crude %>% summary())[2, c(1, 2, 4)]
RES_df_quantiles[6, 3:5] <- coef(glm_P3YR_02_crude %>% summary())[2, c(1, 2, 4)]
RES_df_quantiles[7, 3:5] <- coef(glm_P3YR_03_crude %>% summary())[2, c(1, 2, 4)]
RES_df_quantiles[8, 3:5] <- coef(glm_P3YR_04_crude %>% summary())[2, c(1, 2, 4)]
RES_df_quantiles$TrendP <- rep(c(TrendP_P1YR, TrendP_P3YR), each = 4)

RES_df_full <- RES_df %>% bind_rows(RES_df_quantiles)

write.table(RES_df_full, str_c("data/res/", COHORT, "_PRS.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

### TO DO:
### NOTES:
