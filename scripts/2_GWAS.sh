#!/usr/bin/env bash
### LAST VERSION UPDATE 13 APRIL 2023 (v2.0) - MADE SCRIPT MUCH SHORTER AND MORE READABLE, ADDED R SCRIPTS THAT DO THE BULK OF THE JOINING WHICH WILL BE MORE ROBUST
### THIS SCRIPT TAKES THE OUTPUT FROM `1*.sh` AND PERFORMS THE GWAS ON THE GENOTYPED AND IMPUTED SNP DATA

if ! [[ pwd == /home2/genetics/antobe/projects/MTX_GWAS ]]; then cd /home2/genetics/antobe/projects/MTX_GWAS; fi

LOG=data/info
TMP=TMP/tmp-2
CACHE=TMP/CACHE
RES=data/res
SOFTWARE=../../software

if [[ -f ${LOG}/PLog2.log ]]; then rm ${LOG}/PLog2.log; fi; touch ${LOG}/PLog2.log

### ### ### 2.1. PRE-PROCESSING OF DATA

echo "BEGINNING PRE-PROCESSING OF COVARIATES AND OUTCOMES..."

for COHORT in RA SPOS SNEG SENS VALID; do

	Rscript scripts/misc/CleanVar.R ${COHORT}
	Rscript scripts/misc/PCScree.R ${CACHE}/${COHORT}_QC.eigenval
	read NCOMP; COVFIELDS=$((NCOMP + 4)); cut -f 1-${COVFIELDS} ${TMP}/${COHORT}_COV.txt > ${CACHE}/${COHORT}_COVSUB.txt
	rm ${CACHE}/${COHORT}_QC.eigenval_SCREE.jpeg

done

echo "PRE-PROCESSING OF COVARIATES AND OUTCOMES COMPLETED!"

### ### ### 2.2. RUN ASSOCIATION ANALYSIS

echo "BEGINNING ASSOCIATION ANALYSIS... PLEASE HOLD..."
${SOFTWARE}/plink2 --bfile ${CACHE}/RA_QC --geno-counts --out ${TMP}/FREQ &>> ${LOG}/PLog2.log     #SEE NOTES FOR DETAILS ABOUT THIS

for COHORT in RA SPOS SNEG SENS; do

	tail -n +2 ${CACHE}/${COHORT}_PHENO.txt | cut -f 1,2,3 > ${TMP}/${COHORT}_P1YR.txt
	${SOFTWARE}/plink --bfile ${CACHE}/${COHORT}_QC --pheno ${TMP}/${COHORT}_P1YR.txt --make-bed --out ${TMP}/${COHORT}_P1YR.txt &>> ${LOG}/PLog2.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_P1YR.txt --glm hide-covar --covar ${CACHE}/${COHORT}_COVSUB.txt --covar-variance-standardize --read-freq ${TMP}/FREQ.gcount --out ${RES}/${COHORT}_P1YR --threads 8 &>> ${LOG}/PLog2.log
	echo "FIRST PART OF THE ASSOCIATION ANALYSIS COMPLETED FOR THE ${COHORT} COHORT..."

	tail -n +2 ${CACHE}/${COHORT}_PHENO.txt | cut -f 1,2,4 > ${TMP}/${COHORT}_P3YR.txt
	${SOFTWARE}/plink --bfile ${CACHE}/${COHORT}_QC --pheno ${TMP}/${COHORT}_P3YR.txt --make-bed --out ${TMP}/${COHORT}_P3YR.txt &>> ${LOG}/PLog2.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_P3YR.txt --glm hide-covar --covar ${CACHE}/${COHORT}_COVSUB.txt --covar-variance-standardize --read-freq ${TMP}/FREQ.gcount --out ${RES}/${COHORT}_P3YR --threads 8 &>> ${LOG}/PLog2.log
	echo "SECOND PART OF THE ASSOCIATION ANALYSIS COMPLETED FOR THE ${COHORT} COHORT..."

done

echo "ASSOCIATION ANALYSIS COMPLETED FOR ALL TARGET COHORTS! RESULTS ARE AVAILABLE AT 'data/res/'... CLEANING UP AND EXITING BASH..."
rm ${TMP}/*

### TO DO:
# 1.1. The second line of each chunk within the for loop at the end outputs bfiles with a '.txt' within their name. This is not necessarily wrong but it isn't appropriate and should be changed to solely `${COHORT}_P1YR`.
# 1.2. Selection of the number of PCs via scree plot could be improved by implementing a tool that allows graphical illustration via the command-line. 
#	J gave me a link to such a tool built in R, but I have not tried implementing it here, yet.
#	It would remove the awkward procedure that we currently have implemented.
### NOTES:
