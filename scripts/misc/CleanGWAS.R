### LAST VERSION UPDATE 20 APRIL 2023 (v1.0)
### THIS SCRIPT CLEANS AND OUTPUTS PRESENTABLE DATA RELATED TO THE GWAS RESULTS

.libPaths("/home2/genetics/antobe/software/RLibrary/")

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(argparser)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(openxlsx)))
suppressMessages(suppressWarnings(library(ggplot2)))

ARGS <- arg_parser('EMPTY FN') %>% add_argument('PREFIX', help = 'PREFIX OF THE TARGET GWAS FOR WHICH TO GET THE COVARIATES', type = 'character') %>% parse_args()
PREFIX <- ARGS$PREFIX
GWAS <- read_tsv(str_c("data/res/", PREFIX, ".PHENO1.glm.logistic.hybrid"), show_col_types = F)

### 5.1.1. EXTRACTING RELEVANT GWAS HITS
### ### 5.1.1.1. NUMBER OF SIGNIFICANT AND SUGGESTIVE ASSOCIATIONS

SIGNIF <- GWAS %>% filter(P < 5e-8)
SUGGEST <- GWAS %>% filter(P < 5e-5) %>% mutate(CHRPOS = str_extract(ID, "^\\d+\\:\\d+"))

print(str_c("NUMBER OF SNPS PASSING GENOME-WIDE SIGNIFICANCE THRESHOLD, N = ", nrow(SIGNIF)))
print(str_c("NUMBER OF SNPS PASSING SUGGESTIVE THRESHOLD, N = ", nrow(SUGGEST)))

rm(SIGNIF)

### ### 5.1.1.2. NUMBER OF REGIONAL CLUSTERS OF SUGGESTIVE SNPS

tophits <- SUGGEST %>% arrange(P, `#CHROM`, POS)
cluster_list <- list(); list_pos <- 1
while(nrow(tophits) > 0){
  
  target_SNP <- tophits %>% slice(1)    #TARGET THE FIRST SNP IN THE LIST
  cluster_list[[list_pos]] <- tophits %>% filter(`#CHROM` == target_SNP$`#CHROM`) %>% filter(between(POS, target_SNP$POS - 200 * 1000, target_SNP$POS + 200 * 1000)) %>% mutate(CLUSTER_ID = list_pos)    #IDENTIFY SURROUNDING SNPS
 
  tophits <- tophits %>% filter(!(ID %in% cluster_list[[list_pos]]$ID)) %>% arrange(P, `#CHROM`, POS)    #POP OUT THE TARGET SNP AND THOSE IN THE SURROUNDING AREA
  list_pos <- list_pos + 1
}
N_CLUSTER <- length(cluster_list)
N_CHROM <- unlist(sapply(cluster_list, function(.) .$`#CHROM`)) %>% sort() %>% unique() %>% length()

print(str_c("NUMBER OF SUGGESTIVE REGIONAL CLUSTERS, N = ", N_CLUSTER))
print(str_c("NUMBER OF UNIQUE CHROMOSOMES ON WHICH WE HAVE SUGGESTIVE REGIONAL CLUSTERS, N = ", N_CHROM))

rm(tophits); rm(N_CLUSTER); rm(N_CHROM)

### ### 5.1.1.3. LISTS OF SUGGESTIVE ASSOCIATIONS FOR PRIMARY COHORT

if(PREFIX %in% c("RA_P1YR", "RA_P3YR")){
  
  #rsLinkage <- read_tsv("/home2/genetics/antobe/data/EIRA-SRQB/SNP-ID_linkage.txt", show_col_types = F, col_names = c("CHRPOS", "RSID")) %>% mutate(CHRPOS_CLEAN = str_extract(CHRPOS, "^\\d+\\:\\d+"))
  rsLinkage <- read_tsv("TMP/SNP-ID_linkage.txt", show_col_types = F, col_names = c("CHRPOS", "RSID")) %>% mutate(CHRPOS_CLEAN = str_extract(CHRPOS, "^\\d+\\:\\d+"))
  
  SUGGEST_CLEAN <- bind_rows(cluster_list) %>% 
    left_join(rsLinkage, by = c("CHRPOS" = "CHRPOS_CLEAN")) %>%
    mutate(CI_L = exp(log(OR) + -1 * `LOG(OR)_SE` * qnorm(1 - 0.05 / 2)) %>% round(digits = 2) %>% format(nsmall = 2), CI_U = exp(log(OR) + 1 * `LOG(OR)_SE` * qnorm(1 - 0.05 / 2)) %>% round(digits = 2) %>% format(nsmall = 2)) %>%
      mutate(CI = str_c(CI_L, CI_U, sep = "-")) %>%
    mutate(P_ALT = format(P, scientific = T)) %>%
      mutate(P_NUMB = str_extract(P_ALT, "^\\d\\.\\d+"), P_DEC = str_extract(P_ALT, "e-0\\d$")) %>%
      mutate(P_NUMB_2DIG = as.numeric(P_NUMB) %>% round(digits = 2) %>% format(nsmall = 2)) %>%
      mutate(P_CLEAN = str_c(P_NUMB_2DIG, P_DEC, sep = "")) %>%
      mutate(OR = round(OR, 2) %>% as.character()) %>%
    select(Region = CLUSTER_ID, SNP = RSID, Chromosome = `#CHROM`, Position = POS, `Effect Allele` = ALT, N = OBS_CT, OR, CI, P = P_CLEAN) %>%
    arrange(Chromosome, Position)
  ID_LINKAGE <- SUGGEST_CLEAN %>% distinct(Region) %>% mutate(Region_NEW = 1:nrow(.))
  SUGGEST_CLEAN <- SUGGEST_CLEAN %>% left_join(ID_LINKAGE, by = "Region") %>% select(-Region) %>% select(Region = Region_NEW, everything())
  
  if(PREFIX == "RA_P1YR"){write.xlsx(SUGGEST_CLEAN, "data/res/clean/Table S3.xlsx")}    #SEE NOTES
  if(PREFIX == "RA_P3YR"){write.xlsx(SUGGEST_CLEAN, "data/res/clean/Table S4.xlsx")}
  
}

### 5.1.2. BUILDING AND PRODUCING MANHATTAN-PLOT

if(!(PREFIX %>% str_detect("SENS"))){    #IF WE'RE DOING SENSITIVITY ANALYSIS THEN DON'T GET MANHATTAN

  if(PREFIX %>% str_detect("RA")){fig_width = 1400; fig_height = 600      #IF USING THE RA COHORT THEN MAKE A BIGGER FIGURE
  }else{fig_width = 1400 / 2; fig_height = 600 / 2}    #OTHERWISE USE A SMALLER ONE
  
  POS_HELPER <- GWAS %>% group_by(`#CHROM`) %>% summarise(MAX_POS = as.numeric(max(POS))) %>% mutate(POS_ADD = lag(cumsum(MAX_POS), default = 0)) %>% select(CHR = `#CHROM`, POS_ADD)
  GWAS_CLEAN <- GWAS %>% select(CHR = `#CHROM`, everything()) %>% inner_join(POS_HELPER, by = "CHR") %>% mutate(CUMUL_POS = POS + POS_ADD, COL_HELPER = if_else(CHR %% 2 == 0, 1, 0)) %>%
    mutate(CHRPOS = str_c(CHR, ":", POS)) %>% left_join(cluster_list %>% bind_rows() %>% select(CHRPOS, CLUSTER_ID), by = "CHRPOS") %>%
    group_by(CLUSTER_ID) %>% mutate(TOPHIT = if_else(min(P) == P, 1, 0)) %>% ungroup() %>% mutate(COL_HELPER = if_else(TOPHIT == 1, 2, COL_HELPER)) %>% mutate(COL_HELPER = as.factor(COL_HELPER))
  CHR_MEDPOS <- GWAS_CLEAN %>% group_by(CHR) %>% summarise(MEDIAN_POS = round(median(CUMUL_POS), 0))    #HELPER DATA FOR X-AXIS
 
  tiff(filename = str_c("data/res/clean/", PREFIX, ".tiff"), width = fig_width, height = fig_height, units = "px")
  ggplot(GWAS_CLEAN, aes(x = CUMUL_POS, y = -log10(P), col = factor(COL_HELPER), size = -log10(P))) + geom_point(alpha = 0.75) +
    geom_hline(yintercept = -log10(5e-5), alpha = 0.5, linetype = "dashed") +
    geom_hline(yintercept = -log10(5e-8), alpha = 0.5, linetype = "dashed") +
    scale_x_continuous(breaks = CHR_MEDPOS$MEDIAN_POS, labels = CHR_MEDPOS$CHR) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, -log10(5e-9))) +
    scale_size_continuous(range = c(0.75, 3)) + 
    scale_color_manual(values = c("Grey 5", "Grey 40", "Fire Brick 2")) +
    xlab("") + ylab("-log10(p)") +
    theme_minimal() +
    theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  dev.off() %>% invisible()
  
}

### TO DO:
# 1.1. There is some bug in the plotting that breaks down, uncler what's actually going on here but it isn't outputting figures as I had initially expected.
#	EDIT: Been trying to bugfix this for almost an hour now and am officially stumped. Makes no sense for it to break down? 
#	Individual components of the if-statement works fine but when running it from top-to-bottom it fails? 
#	Need to carefully check this one, probably outside of the command-line... Doesn't matter for now and we can move on without it...
#	EDIT: Try using `ggsave()` instead of `tiff()` or similar default R commands; the former seems much more stable.
### NOTES:
# 2.1. I write to .xslx here because I can't figure out how to fix the data so that, instead of the actual number given for `Region` here, I have them in chronological order.
#       That is, instead of the first number being 35, I have a number 1 and so forth. 
#	Should be a simple enough task that I'm ignoring for now since it isn't crucial.
