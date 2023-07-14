#!/user/bin/env bash
### LAST VERSUPN UPDATE JULY 13 2023 (v1.0)
### THIS SCRIPT READIES DATA FOR PUBLICATION ONLINE

if ! [[ pwd == /home2/genetics/antobe/projects/MTX_GWAS ]]; then cd /home2/genetics/antobe/projects/MTX_GWAS; fi

TMP=TMP/tmp-7
SOFTWARE=/home2/genetics/antobe/software

### 7.1. META-ANALYSIS WITH VALIDATION

for COHORT in P1YR P3YR; do

	#Rscript scripts/misc/CleanVariantNames.R RA_${COHORT}; mv TMP/tmp-2/RA_${COHORT}_IDFIX.PHENO1.glm.logistic.hybrid ${TMP}/RA_${COHORT}_IDFIX.PHENO1.glm.logistic.hybrid
	#Rscript scripts/misc/CleanVariantNames.R VALID_${COHORT}; mv TMP/tmp-2/VALID_${COHORT}_IDFIX.PHENO1.glm.logistic.hybrid ${TMP}/VALID_${COHORT}_IDFIX.PHENO1.glm.logistic.hybrid

	#${SOFTWARE}/generic-metal/metal
	#MARKER ID
	#ALLELE REF ALT
	#EFFECT log(OR)
	#PVALUE P
	#WEIGHT OBS_CT
	#SEPARATOR TAB
	#PROCESS ${TMP}/RA_${COHORT}_IDFIX.PHENO1.glm.logistic.hybrid
	
	#MARKER ID
	#ALLELE REF ALT
	#EFFECT log(OR)
	#PVALUE P
	#WEIGHT OBS_CT
	#SEPARATOR TAB
	#PROCESS ${TMP}/VALID_${COHORT}_IDFIX.PHENO1.glm.logstic.hybrid

	#ANALYZE
	#QUIT

	#mv METAANALYSIS1.TBL ${TMP}/RA_META_P1YR.tsv; rm METAANALYSIS.TBL.info

	echo "WARNING - ALPHA VERSION OF THIS SCRIPT DOES NOT ALLOW CONSTRUCTION OF META-ANALYSIS GWAS WITHIN THE PIPELINE - SEE TO DO FOR SUPPORT"

done

### 7.2. CLEAN VARIABLES AND MAKE TAR BALLS

touch ${TMP}/output.log

for COHORT in RA SPOS SNEG; do

	tail -n +2 TMP/CACHE/${COHORT}_PHENO.txt | cut -f 1,2,3 > ${TMP}/${COHORT}_P1YR.txt
	${SOFTWARE}/plink --bfile TMP/CACHE/${COHORT}_QC --pheno ${TMP}/${COHORT}_P1YR.txt --make-bed --out ${TMP}/${COHORT}_P1YR &>> ${TMP}/output.log		#UPDATE PHENOTYPE
	${SOFTWARE}/plink --bfile ${TMP}/${COHORT}_P1YR --filter-controls --make-bed --out ${TMP}/${COHORT}_P1YR_CTRL &>> ${TMP}/output.log	#EXTRACT CONTROLS
	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_P1YR_CTRL --freq --out ${TMP}/${COHORT}_P1YR_CTRL &>> ${TMP}/output.log	#COMPUTE FREQUNCIES
	Rscript scripts/misc/CleanGWASData.R ${COHORT}_P1YR

	tail -n +2 TMP/CACHE/${COHORT}_PHENO.txt | cut -f 1,2,4 > ${TMP}/${COHORT}_P3YR.txt
	${SOFTWARE}/plink --bfile TMP/CACHE/${COHORT}_QC --pheno ${TMP}/${COHORT}_P3YR.txt --make-bed --out ${TMP}/${COHORT}_P3YR &>> ${TMP}/output.log
	${SOFTWARE}/plink --bfile ${TMP}/${COHORT}_P3YR --filter-controls --make-bed --out ${TMP}/${COHORT}_P3YR_CTRL &>> ${TMP}/output.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_P3YR_CTRL --freq --out ${TMP}/${COHORT}_P3YR_CTRL &>> ${TMP}/output.log
	Rscript scripts/misc/CleanGWASData.R ${COHORT}_P3YR

	echo "FINISHED PROCESSING THE ${COHORT} GWAS DATA..."

done

Rscript scripts/misc/CleanGWASData.R RA_META_P1YR
Rscript scripts/misc/CleanGWASData.R RA_META_P3YR
echo "FINISHED PROCESSING THE META-ANALYSIS GWAS DATA..."

for GWAS in RA_P1YR RA_P3YR RA_META_P1YR RA_META_P3YR SPOS_P1YR SPOS_P3YR SNEG_P1YR SNEG_P3YR; do

	tar -czvf data/public/${GWAS}.tsv.tar.gz ${TMP}/${GWAS}.tsv
	md5sum data/public/${GWAS}.tsv.tar.gz data/public/${GWAS}.md5sum

done

echo "FINISHED PROCESSING PUBLIC DATA FILES! EXITING BASH..."
rm ${TMP}/*

### TO DO:
# 1.1. I couldn't figure out how to get METAL running within the pipeline. 
#	Probably we can create a script that runs METAL but then we need to hassle with taking inputs and other details.
#	I don't have time to fix this right now, so will leave it as a TO-DO for later.
#	Valid code is still available within the loop, and we can input that into the command-line manually to acquire the files.
### NOTES:
