#!/usr/bin/env Rscript
### LAST UPDATED 14 APRIL 2023 (v1.2.4) - FIXED A BUG THAT GAVE THE WRONG FILEPATH FOR THE OUTPUT
### THIS SCRIPT PRODUCES A SCREE-PLOT FOR CHOOSING THE NUMBER OF PCS FOR GWAS

.libPaths("/home2/genetics/antobe/software/RLibrary/")
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(dplyr)))

args <- arg_parser('EMPTY FOR NOW') %>% add_argument('filepath', help = 'FILEPATH TO THE INPUT .EIGENVAL FILE', type = 'character') %>% parse_args()
filepath <- args$filepath

eigenvals <- read_tsv(filepath, col_names = "EIGENVALS", show_col_type = F)
var_expl <- eigenvals / sum(eigenvals$EIGENVALS)

jpeg(filename = paste0(filepath, "_SCREE.jpeg"))
plot(y = var_expl$V1, x = 1:nrow(var_expl),
	xlab = "NUMBER OF COMPONENTS", ylab = "PROPORTION VARIANCE EXPLAINED",
	main = "SCREE PLOT OF VARIANCE EXPLAINED BY PCS")
lines(y = var_expl$V1, x = 1:nrow(var_expl))
dev.off() %>% invisible()

cat(paste0("SCREE PLOT PRODUCED AT '", filepath, "_SCREE.jpeg', PLEASE VISUALLY INSPECT IT TO DECIDE HOW MANY PRINCIPAL COMPONENTS TO INCLUDE AS COVARIATES... \n"))
cat("NOTE THAT YOU CAN NOT INSPECT IT WITHIN THE LINUX SERVER AND THUS NEED TO MOVE IT TO SOMEWHERE WITH IMAGE READING SOFTWARE VIA WINSCP... \n")
cat("HOW MANY COMPONENTS WOULD YOU LIKE TO INCLUDE? OPTIONS INCLUDE ALL INTEGERS FROM 0 TO 10... \n")

### TO DO:
### NOTES:
# 2.1. Note that slightly different versions may coexist across multiple projects. A single standard should be provided to avoid confusion.
