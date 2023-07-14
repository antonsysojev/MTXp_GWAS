# GWAS on persistence to treatment with MTX, in a Swedish population of early RA patients

This repository hosts scripts and relevant files related to the study described in Sysojev et al. (2023). It also contains various forms of cleaned GWAS summary statistic data produced in the project, made publicly available here as well as on the GWAS Catalog (https://www.ebi.ac.uk/gwas/). However, it _does not_ contain any types of raw data used within the project.

# 1. SCRIPTS

Scripts `1*.sh`, `2*.sh` and `3*.sh` performs the cleaning and quality control of raw genotyped and imputed data, genome-wide association analysis and heritability estimation for the phenotypes of persistence to treatment with MTX, at one and three years respectively, in the primary population as well as its various subpopulations. Script `4*.sh` computes RA PRS weights through LDpred2 and uses these to quantify the RA PRS within each individual in our study population, ultimately testing this PRS for an association with our persistence outcomes. Script `5*.sh` processes the output from the four preceding scripts into a more manageable and readable format, with script `6*.R` obtaining the various cohort demographics tables presented within the publication. Lastly, script `7*.sh` takes the finished data and reformats it for sharing online, creating the contents of the `data/public` folder. Scripts from the `misc` subfolder are detailed in the header of the individual scripts.

Some of the underlying scripts from the `misc` folder have not been made publicly available here. This includes verious helper scripts not included within the pipeline as well as scripts obtaining cohort demographics. These would likely only be of interest for researchers with the raw data available to them, and can then be obtained at reasonable request.

All bash scripts and their underlying R scripts were built to be ran on a Linux machine, using default command-line tools as well as PLINK1.9, PLINK2 and R (4.3.1).

# 2. DATA

Public data made available here, is located at ´data/public´. This folder contains the primary RA GWAS, the seropositive RA GWAS and the seronegative RA GWAS for the two persistence phenotypes of persistence to treatment with methotrexate at one and three years, respectively (see Sysojev et al. (2023) for further details). All GWAS files are available in tab-separated format with variables according to what is requested by the GWAS Catalog.

Associated .xlsx-files give further details on the GWAS contents in a tabular format requested by the GWAS Catalog. Files with the `.md5sum` extension contains the ouput of an `md5sum` call in Linux, the same information being reported within the directory README.

The GWAS data from the meta-analysis of the primary RA GWASs pooled with the supplementary data will be made available shortly.
