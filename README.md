# GWAS on persistence to treatment with MTX, in a Swedish population of early RA patients

This repository hosts scripts and relevant files related to the study described in [Sysojev et al. (2023)](https://pubmed-ncbi-nlm-nih-gov.proxy.kib.ki.se/37326842/). It also contains links to the various forms of cleaned GWAS summary statistic data produced in the project, made publicly available through the GWAS Catalog ([https://www.ebi.ac.uk/gwas/](https://www.ebi.ac.uk/gwas/publications/37326842)). However, it _does not_ contain any types of raw data used within the project.

# 1. SCRIPTS

Scripts `1*.sh`, `2*.sh` and `3*.sh` performs the cleaning and quality control of raw genotyped and imputed data, genome-wide association analysis and heritability estimation for the phenotypes of persistence to treatment with MTX, at one and three years respectively, in the primary population as well as its various subpopulations. Script `4*.sh` computes RA PRS weights through LDpred2 and uses these to quantify the RA PRS within each individual in our study population, ultimately testing this PRS for an association with our persistence outcomes. Script `5*.sh` processes the output from the four preceding scripts into a more manageable and readable format, with script `6*.R` obtaining the various cohort demographics tables presented within the publication. Lastly, script `7*.sh` takes the finished data and reformats it for sharing online, creating the contents of the `data/public` folder. Scripts from the `misc` subfolder are detailed in the header of the individual scripts.

Some of the underlying scripts from the `misc` folder have not been made publicly available here. This includes verious helper scripts not included within the pipeline as well as scripts obtaining cohort demographics. These would likely only be of interest for researchers with the raw data available to them, and can then be obtained at reasonable request.

All bash scripts and their underlying R scripts were built to be ran on a Linux machine, using default command-line tools as well as PLINK1.9, PLINK2 and R (4.3.1).

# 2. DATA

Due to storage limitations, no GWAS summary statistic data could be made available here, i.e. through this GitHub page. However, data hosting via the GWAS Catalog ([https://www.ebi.ac.uk/gwas/](https://www.ebi.ac.uk/gwas/publications/37326842)) allows data to be made available elsewhere, easily accessible via link. The two primary GWASs, on persistence to treatment with MTX in DMARD-monotherapy at one and three years, as well as the GWASs on the same outcomes within subcohorts of seropositive and seronegative RA patients, can be found there. However, the meta-analysis GWAS described under the `Supplementary Analysis` is not currently hosted at the same place due to formating issues, but can be made available at reasonable request. Direct links to FPT-folders hosting the various GWASs can be found below.

- [Summary statistics from the GWAS on persistence to MTX, at one year](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90281001-GCST90282000/GCST90281046/).

- [Summary statistics from the GWAS on persistence to MTX, at three years](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90281001-GCST90282000/GCST90281047/).

- [Summary statistics from the GWAS on persistence to MTX, at one year, in seropositive RA patients](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90281001-GCST90282000/GCST90281048/).

- [Summary statistics from the GWAS on persistence to MTX, at three years, in seropositive RA patients](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90281001-GCST90282000/GCST90281049/).

- [Summary statistics from the GWAS on persistence to MTX, at one year, in seronegative RA patients](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90281001-GCST90282000/GCST90281050/).

- [Summary statistics from the GWAS on persistence to MTX, at three years, in seronegative RA patients](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90281001-GCST90282000/GCST90281051/).
