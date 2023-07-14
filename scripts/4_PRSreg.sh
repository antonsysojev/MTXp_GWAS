#!/usr/bin/env bash
### LAST VERSION UPDATE 17 APRIL 2023 (v1.0)
### THIS SCRIPT TAKES THE PREVIOUSLY CONSTRUCTED PRS WEIGHTS AND PERFORMS THE PRS REGRESSION FOR ALL COHORTS

if ! [[ pwd == /home2/genetics/antobe/projects/MTX_GWAS ]]; then cd /home2/genetics/antobe/projects/MTX_GWAS; fi

TMP=TMP/tmp-4
CACHE=TMP/CACHE
LOG=data/info
RES=data/res
RAW=data/raw
SOFTWARE=../../software

if [[ -f ${LOG}/PLog4.log ]]; then rm ${LOG}/PLog4.log; touch ${LOG}/PLog4.log; fi

### ### ### 4.1. PROCESSING OF PRS WEIGHTS AND CREATION OF PRS

echo "BEGINNING CONSTRUCTION OF POLYGENIC RISK SCORES..."
### SEE TO DO 1.1

${SOFTWARE}/plink2 --bfile ${CACHE}/RA_QC --geno-counts --out ${TMP}/FREQ &>> ${LOG}/PLog4.log
for COHORT in RA SPOS SNEG SENS VALID; do if ! [[ -f ${CACHE}/${COHORT}_PRS.txt ]]; then Rscript scripts/misc/CleanPRS.R ${COHORT} 8; fi; done    #SEE NOTES

### ### ### 4.2. TESTING OF PRS ASSOCIATION VIA REGRESSION

for COHORT in RA SPOS SNEG SENS VALID; do

	if [[ ${COHORT} =~ ^VALID$ ]]; then
		Rscript scripts/misc/CleanVar.R VALID; mv TMP/tmp-2/VALID_COV.txt ${TMP}/VALID_COV.txt
		Rscript scripts/misc/PCScree.R ${CACHE}/VALID_QC.eigenval
		read NCOMP; COVFIELDS=$((NCOMP + 4)); cut -f 1-${COVFIELDS} ${TMP}/${COHORT}_COV.txt > ${CACHE}/${COHORT}_COVSUB.txt
		rm ${CACHE}/VALID_QC.eigenval_SCREE.jpeg
	fi

	echo "PERFORMING TESTING FOR THE ${COHORT} COHORT..."
	Rscript scripts/misc/RegPRS.R ${COHORT}

done

echo "POLYGENIC RISK SCORE TESTING COMPLETED... CLEANING UP FOLDER AND EXITING BASH..."
rm ${TMP}/*

### TO DO:
# 1.1. This script simply ASSUMES that PRS weights are available to load from a file, and neither looks for these nor creates them from scratch if finding them to be missing.
#	This is due to complexities regarding the creation of weights (a time-consuming procedure that I decided to move outside of the pipeline).
#	In the future, it may be better to slot these into the pipeline for completeness.
### NOTES:
# 2.1. Testing the code, I consistently find that 8 seems to be a stable number of PRS_WEIGHTS to aggregate into one (i.e. the first 8 of the 15 iterations gave reasonable outputs). 
#	This holds true for the RA, SPOS and SNEG inputs, which is why I've simply made the script set this value ahead of time to avoid unnecessary computations. 
#	EDIT: I've now tested it for the SENS cohort, and 8 seems to suffice there too.
#	EDIT: I've now tested it for the VALID cohort, and 8 seems to suffice there too.
# 2.2. Note that a warning is printed when running this for the first time. 
#	I can't seem to figure out WHY this warning is occurring, but it only seems to appear for SPOS, SNEG and SENS when creating the .rds file from scratch, and not when running the script with that file already existing...
# 2.3. Note that there are warnings printed for SNEG. 
#	I believe I have commented this elsewhere but I will add it here too: these warnings relate to issues converging, which was unfortunately the case previously too. 
#	I decided to move on despite this warnings, since estimates look fine regardless. 
#	One could most likely test other GLM softwares (`logbin::logbin` and `glm2`, for instance) but I decided not to do this, since it works well for all the other inputs currently.
#	Another option would be to use something like modified Poisson regression, but I decided to stick with log-binomial regression.
