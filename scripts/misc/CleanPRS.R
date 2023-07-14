### LAST VERSION UPDATE 19 APRIL 2023 (v0.1.1) - CHANGED THE SCRIPT TO OUTPUT THE PRS, NOT THE WEIGHTS
### THIS SCRIPT CLEANS THE PRS WEIGHTS AND COMPUTES THE PRS FOR A REUSABLE SAMPLE OF INDIVIDUALS

.libPaths("/home2/genetics/antobe/software/RLibrary/")

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(bigsnpr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(argparser)))
suppressMessages(suppressWarnings(library(readr)))

ARGS <- arg_parser('EMPTY FOR NOW') %>% add_argument('COHORT', help = 'ID OF THE TARGET COHORT FOR WHICH TO GET THE COVARIATES', type = 'character') %>% add_argument('FAST', help = 'EITHER AN INTEGER BETWEEN 1 AND 15 OR THE CHARACTER "FALSE"', type = 'character') %>% parse_args()
COHORT <- ARGS$COHORT
FAST <- ARGS$FAST
FILEPATH <- str_c("TMP/CACHE/", COHORT, "_QC")

PRS_BETA <- read_tsv("data/raw/PRS_RA_beta.txt", show_col_types = F) %>% mutate(CHRPOS = str_c(chr, ":", pos))
FREQ <- read_tsv("TMP/tmp-4/FREQ.gcount", show_col_types = F) %>% mutate(CHRPOS = str_extract(ID, "^\\d+\\:\\d+"))

if(!file.exists(str_c(FILEPATH, ".rds"))){snp_readBed(str_c(FILEPATH, ".bed"))}
DATA <- snp_attach(str_c(FILEPATH, ".rds"))
GENLINK <- data.frame(GEN_IDX = 1:ncol(DATA$genotypes), CHRPOS = str_c(DATA$map$chromosome, ":", DATA$map$physical.pos), A1 = DATA$map$allele1, A2 = DATA$map$allele2) %>% as_tibble()

MATRIX_FLIP <- function(GENOTYPE_MATRIX){
  GENOTYPE_MATRIX_CP <- GENOTYPE_MATRIX    #MAKE AN INTERNAL COPY TO NOT OVERWRITE THE GLOBAL FILE
  flip_to_0 <- which(GENOTYPE_MATRIX_CP == 2)
  flip_to_2 <- which(GENOTYPE_MATRIX_CP == 0)
  GENOTYPE_MATRIX_CP[flip_to_0] <- 0
  GENOTYPE_MATRIX_CP[flip_to_2] <- 2
  GENOTYPE_MATRIX_CP
}

### 4.1. CLEANING OF DISCOVERY SAMPLE

print(paste0("CLEANING DISCOVERY SAMPLE FOR THE ", COHORT, " COHORT"))

GENLINK_OVERLAPREF <- GENLINK %>% inner_join(FREQ %>% select(CHRPOS, A1_REF = ALT, A2_REF = REF), by = "CHRPOS")
REMOVEABLE_NOOVERLAP <- GENLINK %>% anti_join(FREQ %>% select(CHRPOS, A1_REF = ALT, A2_REF = REF), by = "CHRPOS")

GENLINK_NONALIGNED <- GENLINK_OVERLAPREF %>% filter(!(A1 == A1_REF & A2 == A2_REF))    #ID SNPS NOT CURRENTLY ALIGNED
FLIPABLE <- GENLINK_NONALIGNED %>% filter(A1 == A2_REF & A2 == A1_REF)    #ID SNPS THAT CAN BE FLIPPED
REMOVEABLE_NOTFLIPABLE <- GENLINK_NONALIGNED %>% filter(!(CHRPOS %in% (FLIPABLE %>% distinct(CHRPOS))$CHRPOS))    #ID SNPS THAT CAN NOT BE FLIPPED WHICH WILL THUSLY BE CUT

REMOVEABLE <- (data.frame(GEN_IDX = c(REMOVEABLE_NOOVERLAP$GEN_IDX, REMOVEABLE_NOTFLIPABLE$GEN_IDX)) %>% distinct(GEN_IDX))$GEN_IDX
GENLINK_CLEAN <- GENLINK %>% filter(!(GEN_IDX %in% REMOVEABLE))

rm(FREQ); rm(GENLINK_OVERLAPREF); rm(REMOVEABLE_NOOVERLAP); rm(GENLINK_NONALIGNED); rm(REMOVEABLE_NOTFLIPABLE); rm(REMOVEABLE)

### 4.2. CLEANING OF PRS WEIGHTS
### ### 4.2.1. CLEANING OF GENOTYPE DATA

print(paste0("CLEANING GENOTYPE DATA FOR THE ", COHORT, " COHORT"))

DATA_GENOTYPES <- DATA$genotypes
DATA_FLIPPED_GENOTYPES <- MATRIX_FLIP(DATA_GENOTYPES[, FLIPABLE$GEN_IDX])    #FLIP THOSE IDENTIFIED AS FLIPABLE
DATA_GENOTYPES[, FLIPABLE$GEN_IDX] <- DATA_FLIPPED_GENOTYPES    #UPDATE WITH THE FLIPPED EFFECTS

PRS_BETA_MERGED <- PRS_BETA %>% inner_join(GENLINK_CLEAN %>% select(GEN_IDX, CHRPOS), by = "CHRPOS")
DATA_GENOTYPES <- DATA_GENOTYPES[, PRS_BETA_MERGED$GEN_IDX]    #EXTRACT ONLY THOSE THAT WEIGHTS ARE AVAILABLE FOR
DATA_GENOTYPES[is.na(DATA_GENOTYPES)] <- 0    #SET THE NAS TO ZERO

rm(PRS_BETA); rm(MATRIX_FLIP); rm(FLIPABLE); rm(DATA_FLIPPED_GENOTYPES)

### ### 4.2.2. SELECTING APPROPRIATE PRS WEIGHTS

print(paste0("SELECTING PRS WEIGHTS FOR THE ", COHORT, " COHORT"))

if(FAST == "FALSE"){    #IF FAST == FALSE THEN WE RUN THE FULL, SLOWER PROCESS
  
  IND_SCORES <- vector("list", 15)
  for(i in 1:15){IND_SCORES[[i]] <- DATA_GENOTYPES %*% unlist(PRS_BETA_MERGED[, 2 + i])}
  IND_SCORES_BOUND <- bind_cols(IND_SCORES) %>% invisible()
  colnames(IND_SCORES_BOUND) <- paste0("SCORE_", 1:ncol(IND_SCORES_BOUND))
  
  sc <- apply(IND_SCORES_BOUND, 2, sd)
  PRS_BETA_MERGED_SUB <- PRS_BETA_MERGED[, 3:17]
  print(paste0("OF 15 POSSIBLE VALUES, ", sum(abs(sc - median(sc)) < 3 * mad(sc)), " WERE CHOSEN FOR THE ", COHORT, " COHORT..."))
  PRS_BETA_MERGED_SUB <- PRS_BETA_MERGED_SUB[, abs(sc - median(sc)) < 3 * mad(sc)]    #NOTE THAT THIS GETS THE FIRST 8, SEEMINGLY ALWAYS
  PRS_BETA_MERGED_SUB_ROWMEANS <- rowMeans(PRS_BETA_MERGED_SUB)    #COMPLETED WEIGHTS
  
}else{    #IF FAST != FALSE THEN WE RUN THE FAST PROCESS - NOTE THAT WE NEED IT TO CONTAIN AN INTEGER HERE OR WE CRASH...

  NCOL = as.numeric(FAST)  
  PRS_BETA_MERGED_SUB <- PRS_BETA_MERGED[, 1:NCOL + 2]
  PRS_BETA_MERGED_SUB_ROWMEANS <- rowMeans(PRS_BETA_MERGED_SUB)
  
}

PRS_vec <- DATA_GENOTYPES %*% PRS_BETA_MERGED_SUB_ROWMEANS
PRS_df <- data.frame(FID = DATA$fam$family.ID, IID = DATA$fam$sample.ID, RA_PRS = PRS_vec)
write.table(PRS_df, paste0("TMP/CACHE/", COHORT, "_PRS.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
print(paste0("PRS SUCCESFULLY CONSTRUCTED FOR THE ", COHORT, " COHORT... RESULTS WRITTEN TO 'TMP/CACHE/'", COHORT, "_PRS.txt"))

### TO DO:
# 1.1. Silencing doesn't seem to work as expected with `invisible()`. Not sure how to silence the new names when binding columns...
#	EDIT: This can be done by giving an additional argument when loading `dplyr` I believe. Something similar was done in a different project of yours.
### NOTES:
