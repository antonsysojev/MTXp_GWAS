#!/user/bin/env bash
### LAST UPDATED 17 APRIL 2023 (v2.0)
### THIS SCRIPT EXTRACTS THE HERITABILITY ESTIMATES VIA GCTA FOR OUR MTX GWAS PHENOTYPES

if ! [[ pwd == /home2/genetics/antobe/projects/MTX_GWAS ]]; then cd /home2/genetics/antobe/projects/MTX_GWAS; fi

LOG=data/info
TMP=TMP/tmp-3
CACHE=TMP/CACHE
RES=data/res
SOFTWARE=../../software

if [[ -f ${LOG}/PLog3.log ]]; then rm ${LOG}/PLog3.log; fi; touch ${LOG}/PLog3.log

### 3.1. PRE-PROCESSING OF GCTA INPUT DATA

echo "BEGINNING PRE-PROCESSING OF DATA..."

for COHORT in RA SPOS SNEG SENS; do if ! [[ -f ${TMP}/${COHORT}_QC.grm ]]; then ${SOFTWARE}/gcta-1.94.1 --bfile ${CACHE}/${COHORT}_QC --make-grm --out ${TMP}/${COHORT}_QC --thread-num 8 &>> ${LOG}/PLog3.log; fi; done

### 3.2. ESTIMATION OF HERITABILITY

echo "BEGINNING HERITABILITY ESTIMATION..."

for COHORT in RA SPOS SNEG SENS; do

	tail -n +2 ${CACHE}/${COHORT}_PHENO.txt | cut -f 1,2,3 > ${TMP}/${COHORT}_P1YR.txt
	tail -n +2 ${CACHE}/${COHORT}_PHENO.txt | cut -f 1,2,4 > ${TMP}/${COHORT}_P3YR.txt
	tail -n +2 ${CACHE}/${COHORT}_COVSUB.txt | cut -f 1,2,4 > ${TMP}/${COHORT}_CCOV.txt
	tail -n +2 ${CACHE}/${COHORT}_COVSUB.txt | cut -f 4 --complement > ${TMP}/${COHORT}_QCOV.txt

	${SOFTWARE}/gcta-1.94.1 --reml --grm ${TMP}/${COHORT}_QC --pheno ${TMP}/${COHORT}_P1YR.txt --prevalence 0.66 --out ${RES}/${COHORT}_P1YR --thread-num 8 &>> ${LOG}/PLog3.log
	${SOFTWARE}/gcta-1.94.1 --reml --grm ${TMP}/${COHORT}_QC --pheno ${TMP}/${COHORT}_P3YR.txt --prevalence 0.50 --out ${RES}/${COHORT}_P3YR --thread-num 8 &>> ${LOG}/PLog3.log

	${SOFTWARE}/gcta-1.94.1 --reml --grm ${TMP}/${COHORT}_QC --pheno ${TMP}/${COHORT}_P1YR.txt --covar ${TMP}/${COHORT}_CCOV.txt --qcovar ${TMP}/${COHORT}_QCOV.txt --prevalence 0.66 --out ${RES}/${COHORT}_P1YR_adj --thread-num 8 &>> ${LOG}/PLog3.log
	${SOFTWARE}/gcta-1.94.1 --reml --grm ${TMP}/${COHORT}_QC --pheno ${TMP}/${COHORT}_P3YR.txt --covar ${TMP}/${COHORT}_CCOV.txt --qcovar ${TMP}/${COHORT}_QCOV.txt --prevalence 0.50 --out ${RES}/${COHORT}_P3YR_adj --thread-num 8 &>> ${LOG}/PLog3.log

	echo "ESTIMATION OF HERITABILITY SUCCESFULLY COMPLETED FOR THE ${COHORT} COHORT..."

done

echo "ALL HERITABILITY ESTIMATION COMPLETED! CLEARING FOLDERS AND EXITING BASH..."
rm ${TMP}/*

### TO DO:
### NOTES:
