#!/user/bin/env bash
### LAST VERSION UPDATE APRIL 27 2023 (v2.1) - ADDED CODE FOR PROCESSING OF VALIDATION DATA INTO EXISTING PIPELINE
### THIS SCRIPT PERFORMS THE PRE-PROCESSING OF THE RAW GWAS DATA FOR ANALYSIS OF THE MTX GWAS

if ! [[ pwd == /home2/genetics/antobe/projects/MTX_GWAS ]]; then cd /home2/genetics/antobe/projects/MTX_GWAS; fi

LOG=data/info
TMP=TMP/tmp-1
CACHE=TMP/CACHE
SOFTWARE=/home2/genetics/antobe/software
GENOTYPED=/home2/genetics/antobe/data/EIRA-SRQB
VALID=/home2/genetics/antobe/data/EIRA-SRQB/SRQB/IMPUTED/RAW
RSQ=0.70

if [[ -f ${LOG}/PLog1.log ]]; then rm ${LOG}/PLog1.log; fi; touch ${LOG}/PLog1.log
if [[ -f ${LOG}/info.log ]]; then rm ${LOG}/info.log; fi; touch ${LOG}/info.log
if [[ -f ${LOG}/info_VALID.log ]]; then rm ${LOG}/info_VALID.log; fi; touch ${LOG}/info_VALID.log

### ### ### 1.1. CLEANING THE RAW INPUT DATA

if ! [[ -f data/raw/KEY_pid-GWAS.txt ]]; then echo "FAILED TO FIND THE KEY FOR LINKING 'pid' TO 'GWAS' AT 'data/raw/KEY_pid-GWAS.txt', EXITING BASH..."; exit 1; fi     #CLEAN PRIMARY COHORT DATA
Rscript scripts/misc/CleanEIRA-SRQB.R
cut -f 2,3 ${CACHE}/EIRA-SRQB.txt >> ${TMP}/EIRA-SRQB_INDS.txt		#SEE NOTE ABOUT THIS FILTERING
${SOFTWARE}/plink2 --bfile ${GENOTYPED}/eira-plus-others-imputed --keep ${TMP}/EIRA-SRQB_INDS.txt --make-bed --out ${TMP}/EIRA-SRQB &>> ${LOG}/PLog1.log

echo "GENOTY $(wc -l < ${GENOTYPED}/eira-plus-others-imputed.fam) $(wc -l < ${GENOTYPED}/eira-plus-others-imputed.bim)" >> ${LOG}/info.log
echo "COHORT $(tail -n +2 ${CACHE}/EIRA-SRQB.txt | cut -f 1 | sort | uniq | wc -l)" >> ${LOG}/info.log
echo "OVERLP $(cut -f 2 ${TMP}/EIRA-SRQB.fam | sort | uniq | wc -l)" >> ${LOG}/info.log

for CHR in $(seq 1 22); do ${SOFTWARE}/plink2 --bfile ${VALID}/SRQB_IMPRAW_${CHR} --maf 0.00001 --make-pgen --out ${TMP}/SRQB_i${CHR} --threads 8 &>> ${LOG}/PLog1.log; echo "${TMP}/SRQB_i${CHR}" >> ${TMP}/PMergeListFile.txt; done     #MERGE VALIDATION DATA
${SOFTWARE}/plink2 --pmerge-list ${TMP}/PMergeListFile.txt --make-bed --out ${TMP}/SRQB_i &>> ${LOG}/PLog1.log
Rscript scripts/misc/CleanSRQB.R
tail -n +2 ${CACHE}/SRQB.txt | cut -f 2,3 | grep -v NA > ${TMP}/SRQB_INDS.txt
${SOFTWARE}/plink2 --bfile ${TMP}/SRQB_i --keep ${TMP}/SRQB_INDS.txt --make-bed --out ${TMP}/SRQB &>> ${LOG}/PLog1.log

echo "GENOTY $(wc -l < ${TMP}/SRQB_i.fam) $(wc -l < ${TMP}/SRQB_i.bim)" >> ${LOG}/info_VALID.log
echo "COHORT $(tail -n +2 ${CACHE}/SRQB.txt | cut -f 1 | sort | uniq | wc -l)" >> ${LOG}/info_VALID.log
echo "OVERLP $(cut -f 2 ${TMP}/SRQB.fam | sort | uniq | wc -l)" >> ${LOG}/info_VALID.log

### ### ### 1.2. QUALITY CONTROL OF THE CLEANED INPUT DATA

echo "BEGINNING QUALITY CONTROL OF THE CLEANED INPUT DATA"

tail -n +2 ${GENOTYPED}/All_chromosomes_info.txt | cut -f 1,7 | awk '{ if ($2 >= 0.70 ) { print } }' > ${TMP}/RSQ070.txt
${SOFTWARE}/plink2 --bfile ${TMP}/EIRA-SRQB --extract ${TMP}/RSQ070.txt --make-bed --out ${TMP}/EIRA-SRQB_RSQ070 --threads 8 &>> ${LOG}/PLog1.log

tail -n +2 ${VALID}/impqual.txt | awk '{ if ($2 >= 0.70 ) { print } }' > ${TMP}/RSQ070_VALID.txt
${SOFTWARE}/plink2 --bfile ${TMP}/SRQB --extract ${TMP}/RSQ070_VALID.txt --make-bed --out ${TMP}/SRQB_RSQ070 --threads 8 &>> ${LOG}/PLog1.log

for COHORT in RA SPOS SNEG SENS VALID; do if ! [[ -f ${CACHE}/${COHORT}_QC.bed ]]; then

	echo "PERFORMING QUALITY CONTROL FOR THE ${COHORT} COHORT... PLEASE HOLD..."
	if [[ ${COHORT} =~ ^RA|SPOS|SNEG|SENS$ ]]; then
		awk -v col=${COHORT} 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $2,$3,$c} NR>1{print $2,$3,$c}' ${CACHE}/EIRA-SRQB.txt | awk '$3 == 1' > ${TMP}/${COHORT}.txt
		${SOFTWARE}/plink2 --bfile ${TMP}/EIRA-SRQB_RSQ070 --keep ${TMP}/${COHORT}.txt --make-bed --out ${TMP}/${COHORT} &>> ${LOG}/PLog1.log

		${SOFTWARE}/plink --bfile ${TMP}/${COHORT} --check-sex 0.2 0.8 --out ${TMP}/${COHORT}_CHKSX &>> ${LOG}/PLog1.log
		grep PROBLEM ${TMP}/${COHORT}_CHKSX.sexcheck > ${TMP}/${COHORT}_SXCHK_PROBLEM.txt
		${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT} --remove ${TMP}/${COHORT}_SXCHK_PROBLEM.txt --chr 1-22 --make-bed --out ${TMP}/${COHORT}_SXCHK &>> ${LOG}/PLog1.log; fi
	if [[ ${COHORT} =~ ^VALID$ ]]; then ${SOFTWARE}/plink2 --bfile ${TMP}/SRQB_RSQ070 --chr 1-22 --make-bed --out ${TMP}/${COHORT}_SXCHK &>> ${LOG}/PLog1.log; fi    #SEE NOTES ABOUT THE SEPARATION

	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_SXCHK --mind 0.05 --make-bed --out ${TMP}/${COHORT}_SXCHK_MIND &>> ${LOG}/PLog1.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_SXCHK_MIND --geno 0.05 --make-bed --out ${TMP}/${COHORT}_SXCHK_MIND_GENO &>> ${LOG}/PLog1.log
    	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_SXCHK_MIND_GENO --maf 0.01 --make-bed --out ${TMP}/${COHORT}_SXCHK_MIND_GENO_MAF &>> ${LOG}/PLog1.log
	${SOFTWARE}/plink --bfile ${TMP}/${COHORT}_SXCHK_MIND_GENO_MAF --hwe 1e-6 include-nonctrl --make-bed --out ${TMP}/${COHORT}_pQC &>> ${LOG}/PLog1.log

	${SOFTWARE}/plink --bfile ${TMP}/${COHORT}_pQC --indep-pairwise 100 5 0.2 --out ${TMP}/${COHORT}_pQC_LDPRUNE &>> ${LOG}/PLog1.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_pQC --extract ${TMP}/${COHORT}_pQC_LDPRUNE.prune.in --make-bed --out ${TMP}/${COHORT}_pQC_LDPRUNE &>> ${LOG}/PLog1.log
	${SOFTWARE}/plink --bfile ${TMP}/${COHORT}_pQC_LDPRUNE --genome --min 0.125 --out ${TMP}/${COHORT}_pQC_LDPRUNE_REL &>> ${LOG}/PLog1.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_pQC_LDPRUNE --remove ${TMP}/${COHORT}_pQC_LDPRUNE_REL.genome --make-bed --out ${TMP}/${COHORT}_pQC_LDPRUNE_REL &>> ${LOG}/PLog1.log

	PCA_INPUT=${COHORT}_pQC_LDPRUNE_REL
	for ITER in $(seq 1 5); do
		${SOFTWARE}/plink2 --bfile ${TMP}/${PCA_INPUT} --pca 10 --out pca --threads 8 &>> ${LOG}/PLog1.log
		Rscript scripts/misc/pcafiltration.R
		${SOFTWARE}/plink2 --bfile ${TMP}/${PCA_INPUT} --remove pca_outliers.txt --out ${TMP}/${PCA_INPUT}_PCA${ITER} --make-bed &>> ${LOG}/PLog1.log
		rm pca*; PCA_INPUT=${PCA_INPUT}_PCA${ITER}
	done
	${SOFTWARE}/plink2 --bfile ${TMP}/${COHORT}_pQC --keep ${TMP}/${PCA_INPUT}.fam --make-bed --out ${CACHE}/${COHORT}_QC &>> ${LOG}/PLog1.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${PCA_INPUT} --pca 10 --out ${CACHE}/${COHORT}_QC --threads 8 &>> ${LOG}/PLog1.log    #THIS COULD BE SKIPPED FOR THE VALIDATION

	if [[ ${COHORT} =~ ^RA|SPOS|SNEG|SENS$ ]]; then cat ${LOG}/info.log > ${LOG}/${COHORT}_info.log; echo "IMPFLT $(wc -l < ${TMP}/EIRA-SRQB_RSQ070.fam) $(wc -l < ${TMP}/EIRA-SRQB_RSQ070.bim)" >> ${LOG}/${COHORT}_info.log
		echo "COHSUB $(wc -l < ${TMP}/${COHORT}.fam) $(wc -l < ${TMP}/${COHORT}.bim)" >> ${LOG}/${COHORT}_info.log; fi
	if [[ ${COHORT} =~ ^VALID$ ]]; then cat ${LOG}/info_VALID.log > ${LOG}/${COHORT}_info.log; echo "IMPFLT $(wc -l < ${TMP}/SRQB_RSQ070.fam) $(wc -l < ${TMP}/SRQB_RSQ070.bim)" >> ${LOG}/${COHORT}_info.log; fi
	echo "SEXCHK $(wc -l < ${TMP}/${COHORT}_SXCHK.fam) $(wc -l < ${TMP}/${COHORT}_SXCHK.bim)" >> ${LOG}/${COHORT}_info.log
	echo "MINDFL $(wc -l < ${TMP}/${COHORT}_SXCHK_MIND.fam) $(wc -l < ${TMP}/${COHORT}_SXCHK_MIND.bim)" >> ${LOG}/${COHORT}_info.log
	echo "GENOFL $(wc -l < ${TMP}/${COHORT}_SXCHK_MIND_GENO.fam) $(wc -l < ${TMP}/${COHORT}_SXCHK_MIND_GENO.bim)" >> ${LOG}/${COHORT}_info.log
	echo "MAFFLT $(wc -l < ${TMP}/${COHORT}_SXCHK_MIND_GENO_MAF.fam) $(wc -l < ${TMP}/${COHORT}_SXCHK_MIND_GENO_MAF.bim)" >> ${LOG}/${COHORT}_info.log
	echo "HWEFLT $(wc -l < ${TMP}/${COHORT}_pQC.fam) $(wc -l < ${TMP}/${COHORT}_pQC.bim)" >> ${LOG}/${COHORT}_info.log
	echo "LDFILT $(wc -l < ${TMP}/${COHORT}_pQC_LDPRUNE.fam) $(wc -l < ${TMP}/${COHORT}_pQC_LDPRUNE.bim)" >> ${LOG}/${COHORT}_info.log
	echo "RELFLT $(wc -l < ${TMP}/${COHORT}_pQC_LDPRUNE_REL.fam) $(wc -l < ${TMP}/${COHORT}_pQC_LDPRUNE_REL.bim)" >> ${LOG}/${COHORT}_info.log
	echo "PCAFLT $(wc -l < ${CACHE}/${COHORT}_QC.fam) $(wc -l < ${CACHE}/${COHORT}_QC.bim)" >> ${LOG}/${COHORT}_info.log

	rm ${TMP}/${COHORT}*
	echo "COMPLETED QUALITY CONTROL FOR THE ${COHORT} COHORT!"

fi; done

echo "PRE-PROCESSING COMPLETED FOR ALL OF THE TARGET COHORTS, CLEANING UP FOLDERS AND EXITING BASH..."
rm ${TMP}/*; rm ${LOG}/info.log; rm ${LOG}/VALID_info.log

### TO DO:
# 1.1.
### NOTES:
# 2.1. When filtering out individuals, we do not care if they are missing a GWAS ID. In other words, a large set of the filtered cohort will lack a GWAS ID and thus have an NA as FID and IID. 
#	When extracting individuals from the genotyped EIRA-SRQB data, we do not take this into account, meaning that if an individual in the genotyped EIRA-SRQB data had FID and IID set to NA, they would be kept too. 
#	However, no such individual exists and we thus run into no issues at this point.
# 2.2. We have to split the introductory parts into two to accomodate differences between the validation data and the primary cohort data. 
#	Firstly, primary cohort data contains four subsets of data which need to be handled only when the input is of this type. 
#	However, I've also included the filtering on sex mismatch here, since I lack imputed data on the X-chromosome. 
#	This way, I skip the check for a sex mismatch within the validation data and only update the name to match the coming lines of code.
