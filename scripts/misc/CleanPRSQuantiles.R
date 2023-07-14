### LAST VERSION UPDATE 28 APRIL 2023 (v1.0)
### THIS SCRIPT CLEANS THE PRS QUANTILE OUTPUT AND HARMONIZES IT INTO A SINGLE FILE

.libPaths("/home2/genetics/antobe/software/RLibrary")

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(openxlsx)))

df1 <- read_tsv("TMP/tmp-5/RA_TableS5.txt", show_col_types = F)
df2 <- read_tsv("TMP/tmp-5/SPOS_TableS5.txt", show_col_types = F)
df3 <- read_tsv("TMP/tmp-5/SNEG_TableS5.txt", show_col_types = F)

df_max <- df1 %>% select(ID, RA = RR) %>% inner_join(df2 %>% select(ID, SPOS = RR), by = "ID") %>% inner_join(df3 %>% select(ID, SNEG = RR), by = "ID")
df_p1yr <- df_max %>% slice(1:4); df_p3yr <- df_max %>% slice(5:8)
df_wide <- df_p1yr %>% select(RA_P1YR = RA, SPOS_P1YR = SPOS, SNEG_P1YR = SNEG) %>% bind_cols(df_p3yr %>% select(RA_P3YR = RA, SPOS_P3YR = SPOS, SNEG_P3YR = SNEG))
df_complete <- df_wide %>% rbind(c(df1$PForTrend %>% unique(), df2$PForTrend %>% unique(), df3$PForTrend %>% unique()))

write.xlsx(df_complete, "data/res/clean/Table S5.xlsx")
