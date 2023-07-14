#!/usr/bin/env bash
### LAST VERSION UPDATE 27 APRIL 2023 (v1.1) - FIXED SOME MINOR BUGS
### THIS SCRIPT PRODUCES ALL THE COUNTS AND MANUSCRIPT FILES THAT DO NOT REQUIRE ACCESS TO K FOLDERS

if ! [[ pwd == /home2/genetics/antobe/projects/MTX_GWAS ]]; then cd /home2/genetics/antobe/projects/MTX_GWAS; fi

if [[ -f ${LOG}/PLog4.log ]]; then rm ${LOG}/PLog4.log; touch ${LOG}/PLog4.log; fi

### ### ### 5.1. PROCESSING OF GWAS OUTPUT

for COHORT in RA SPOS SNEG; do 

	Rscript scripts/misc/CleanGWAS.R ${COHORT}_P1YR
	Rscript scripts/misc/CleanGWAS.R ${COHORT}_P3YR
	echo "SUCCESFULLY PROCESSED THE GWAS OUTPUT FOR THE ${COHORT} COHORT!"

done

Rscript scripts/misc/ReevalGWAS.R P1YR    #RE EVALUATION OF SENSITIVITY ANALYSIS GWAS RESULTS
Rscript scripts/misc/ReevalGWAS.R P3YR

### ### ### 5.2. PROCESSING HERITABILITY AND PRS OUTPUT

touch data/res/clean/HeritPRS.txt
for COHORT in RA SPOS SNEG SENS VALID; do 

	Rscript scripts/misc/CleanHeritPRS.R ${COHORT}
	cat TMP/tmp-5/${COHORT}_heritPRS.txt >> data/res/clean/HeritPRS.txt 
	echo "SUCCESFULLY PROCESSED THE HERITABILITY AND PRS RESULTS FOR THE ${COHORT} COHORT!"

done

Rscript scripts/misc/CleanPRSQuantiles.R

echo "SUCCESFULLY PROCESSED ALL OUTPUT FROM ASSOCIATION ANALYSIS, HERITABILITY ESTIMATION AND POLYGENIC RISK SCORE TESTING..."
rm TMP/tmp-5/*

### TO DO:
### NOTES:
