# stressful-life-event-ewas

Scripts used to examine the epigenetic correlates of Stressful Life Events (SLE) in the UK population.

Forthcoming preprint title: "The epigenetic correlates of stressful life events in two large UK population samples."

# Prerequisites:

Relevant R packages must be installed:

bigmelon, BiocManager, broom.mixed, DescTools, devtools, dplyr, DunedinPACE, DunedinPoAm38, EpiDISH, EpiSmokEr, haven, gamlss, IlluminaHumanMethylationEPICanno.ilm10b2.hg19, limma, lme4, lmerTest, missMethyl, tidyr, wateRmelon

Check all scripts and amend directories and file names as appropriate for your project.

To obtain suitable DNA methylation data with linked survey data, an application must be made to ISER.

You will need to download:
 - all three supplementary table files from McCartney et al (2016) https://www.sciencedirect.com/science/article/pii/S221359601630071X
 - the ARIES middle-aged data file at http://www.mqtldb.org/db/PLINK/FOM.ALL.M.tab
 - the GoDMC data file at http://mqtldb.godmc.org.uk/downloads (choose "assoc_meta_all.csv.gz")

in order to exclude polymorphic probes, cross-reactive probes, and genetically-influenced probes.
Run the following line of code in the command line to unzip the GoDMC data and split into columns:

gunzip -c assoc_meta_all.csv.gz | awk '{split ($0,a,","); print a[1]}' > assoc_meta_all_probes.csv

The file "deriving_complex_sle_variables.R" was sent to the data team at ISER to derive complex SLE (i.e. those involving data from multiple files, spefcifically "household member died", "premature or low birth weight child", and "stillbirth or infant death"). The team then merged this with the other requested linked survey variables.

# Run scripts:

Run "derive_sles.R" to process the survey data and generate a file for merging with episcores prior to analyses: UKHLS_SLEs.rds. N.B. this script contains extra code for measuring SLE at multiple waves beyond those used in the present study.

Run "estimate_cell_proportions.R" to read in the methylation data, calculate leukocyte proportion estimates, read in survey data files (two default epigenetic covariate files containing batch numbers etc, and UKHLS_SLEs.rds), estimate smoking episcore, estimate alcohol consumption by deriving a phosphatidylethanol episcore, calculate epigenetic clocks, ensure certain variables are stored as factors or numeric, and generate three files for use in analyses: clocks.rds, USM1survey_and_covar2.rds and USM2survey_and_covar2.rds.

Run "extract_betas.R" to read in the methylation data in .gds format, extract the beta matrices, remove polymorphic and cross-reactive probes, remove genotype-influenced probes, winsorize the beta values, and generate two files for use in EWAS: USM1betamatrixWinsorized.rds and USM2betamatrixWinsorized.rds.

Run "ewas_script_usm1.R" and "ewas_script_usm2.R" to run the main EWAS, and "ewas_script_usm1_ghqsleint.R" and "ewas_script_usm2_ghqsleint.R" to run the sensitivity analyses (note: these scripts were heavily influenced by https://github.com/psychiatric-genomics-consortium/mdd-mwas). These will produce a nice table of results for the SLE main effect term (suffix "toptable.txt") and a large table of results for all model terms with fewer columns (suffix "linear.txt").

Run "grab_interaction_terms.sh" to extract the EWAS results for the SLE interaction terms, from the "linear.txt" files.

Run "ewasresults.R" to process results and run enrichment analyses.

Run "clock_and_wbc_analyses.R" to run analyses of epigenetic clocks and leukoctes.
