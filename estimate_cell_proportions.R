# The lines of code which estimate cell proportions were based on 
# code from supervisor Prof Leo Schalkwyk

# Script to be run from grandparent directory, relative to "scripts"

rm(list = ls())

# library( devtools )
# install_github( "https://github.com/sjczheng/EpiDISH.git" )
# install_github("sailalithabollepalli/EpiSmokEr")

library( EpiDISH )
library( bigmelon )
library( haven )
library( EpiSmokEr )
library( dplyr )
library( wateRmelon )
library( DunedinPoAm38 )
library( DunedinPACE )

data("age_coefficients")

# BiocManager::install("FlowSorted.Blood.450k") 

data( "centUniLIFE.m" )

if( file.exists( "unilife.rdata" )){ load( "unilife.rdata" ) }

USM1 <- openfn.gds( "ApplicationNumber.gds" )
USM2 <- openfn.gds( "ApplicationNumber_2020.gds" )
tesch <- rownames( centUniLIFE.m )
U1 <- betas( USM1 )
U2 <- betas( USM2 )
# Remove samples with no sample ID
U1 <- U1[nchar(colnames(U1)) > 0,]
U2 <- U2[nchar(colnames(U2)) > 0,]
probes <- intersect( rownames(U1), rownames(U2) )
df <- cbind( U1[ probes, ], U2[ probes, ] )
print(dim(df))
u <- df[ rownames(df) %in% tesch, ]
allests <- epidish( u, ref.m=centUniLIFE.m, maxit=500 )

save.image( "unilife.rdata" )

covar1 <- read_dta( "ApplicationNumber_epi_sendout.dta" )
covar2 <- read_dta( "ApplicationNumber_epi_2020_sendout.dta" )

estsdf <- data.frame( allests$estF )
estsdf$id <- rownames( estsdf )

covar1 <- covar1[,c("id","nsex","confage","rackbarcode","cd8t","cd4t","nk","bcell","mono","gran")]
covar2 <- covar2[,c("id","nsex","confage","rackbarcode","cd8t","cd4t","nk","bcell","mono","gran")]

colnames(covar1)[5:10] <- paste0("old_",colnames(covar1)[5:10])
colnames(covar2)[5:10] <- paste0("old_",colnames(covar2)[5:10])

covar1$usm <- 1
covar2$usm <- 2

covar <- rbind( covar1, covar2 )

covar <- merge( covar, estsdf, by="id", all.x=T, all.y=T )

surv <- readRDS( "UKHLS_SLEs.rds" )

surv$nsex <- NULL
surv$confage <- NULL
surv$cd8t <- NULL
surv$cd4t <- NULL
surv$nk <- NULL
surv$bcell <- NULL
surv$mono <- NULL
surv$gran <- NULL
surv$rackbarcode <- NULL

surv <- merge( surv, covar, by="id", all.x=T, all.y=T )

# Use EpiSmokEr to make episcores for smoking 
# https://htmlpreview.github.io/?https://github.com/sailalithabollepalli/EpiSmokEr/blob/master/vignettes/epismoker.html

surva <- na.omit( surv[surv$id %in% colnames(df), c("id","sex")] )
df <- df[,colnames(df) %in% surva$id]
rownames(surva) <- surva$id

result <- epismoker(dataset=df, samplesheet = surva, method = "SSt")
colnames( result ) <- c("id", paste0( "SSt_", colnames(result)[2:ncol(result)]))
surva <- merge( surva, result, by="id", all.x=T, all.y=T )

result <- epismoker(dataset=df, method = "SSc")
colnames( result ) <- c("id", "SSc_smokingScore")
surva <- merge( surva, result, by="id", all.x=T, all.y=T )

df_smok <- surva[,c("id","SSt_logOdds_CS", "SSt_logOdds_FS", "SSt_logOdds_NS", "SSt_probs_CS",
                            "SSt_probs_FS", "SSt_probs_NS", "SSt_PredictedSmokingStatus",
                            "SSc_smokingScore")]

surv <- merge( surv, df_smok, by="id", all.x=T, all.y=T )

# Make alcohol exposure scores based on 
# Liang, X., Justice, A.C., So-Armah, K. et al. 
# DNA methylation signature on phosphatidylethanol, not on self-reported alcohol consumption, 
# predicts hazardous alcohol consumption in two distinct populations. 
# Mol Psychiatry 26, 2238–2253 (2021). https://doi.org/10.1038/s41380-020-0668-x
# Supplementary table S2

# Pre-checks indicate that "cg20525486" and "cg04304130" are missing from the data
# so omit these

alc_probes <- c("cg06690548","cg17962756","cg13442969","cg11376147","cg18590502",
                "cg00220102", "cg26689780","cg25221975", "cg20970380")
alc_weights <- c(-6.24, -5.58, -5.40, -5.18, -5.13, -5.11, 5.08, 5.07, -4.98 )

dfalc <- df[rownames(df) %in% alc_probes,]

dfalc <- data.frame(t(dfalc))

dfalc <- dfalc[,alc_probes]

for(i in 1:9){ dfalc[,i] <- dfalc[,i] * alc_weights[i] }

dfalc$alc_exposure <- rowSums( dfalc )

dfalc$id <- rownames( dfalc )

surv <- merge( surv, dfalc[,c("id","alc_exposure")], by="id", all.x=T, all.y=T )
surv$alc_exposure <- as.numeric(scale(surv$alc_exposure))

# Make epigenetic clocks

horvath <- agep(df, ageCoefs$Horvath, verbose=T) 
df <- rbind(df, as.matrix(horvath)[,1])
rownames(df)[nrow(df)] <- "horvath"

hannum <- agep(df, ageCoefs$Hannum, verbose=T)
df <- rbind(df, as.matrix(hannum)[,1])
rownames(df)[nrow(df)] <- "hannum"

horvath_sb <- agep(df, ageCoefs$SkinBlood, verbose=T)
df <- rbind(df, as.matrix(horvath_sb)[,1])
rownames(df)[nrow(df)] <- "horvath_sb"

levine <- agep(df, ageCoefs$PhenoAge, verbose=T)
df <- rbind(df, as.matrix(levine)[,1])
rownames(df)[nrow(df)] <- "levine"

lin <- agep(df, ageCoefs$Lin, verbose=T)
df <- rbind(df, as.matrix(lin)[,1])
rownames(df)[nrow(df)] <- "lin"

poam <- DunedinPoAm38::PoAmProjector(df)
poam2 <- t(as.matrix(unlist(poam)))
colnames(poam2) <- gsub("DunedinPoAm_38.","", colnames(poam2))
df <- rbind(df,poam2)
rownames(df)[nrow(df)] <- "poam38"

pace <- DunedinPACE::PACEProjector(df)
pace2 <- t(as.matrix(unlist(pace)))
colnames(pace2) <- gsub("DunedinPACE.","", colnames(pace2))
df <- rbind(df,pace2)
rownames(df)[nrow(df)] <- "pace"

tmp <- data.frame(t(df[c("horvath","hannum","horvath_sb","levine","lin","poam38","pace"),]))
tmp$id <- rownames(tmp)

saveRDS( tmp, "clocks.rds" )

columns <- c("id","usm","wave","sex_categorical","confage",
                  "rackbarcode","blooddelay","blooddelay_2","medtyp4",
                  "whr","B","CD4T","CD8T","Mono","nRBC","Gran",
                  "NK","aCD4Tnv","aBaso","aCD4Tmem","aBmem",
                  "aBnv","aTreg","aCD8Tmem","aCD8Tnv","aEos",
                  "aNK","aNeu","aMono",
                  "family_disruption","family_intact","father_employed",
                  "PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8", 
                  "PCA9", "PCA10",
                  "blood_sle_personalloss","blood_sle_socioeconomicloss",
                  "blood_sle_caring","blood_sle_health","blood_sles_2",
                  "blood_sles_binary_2","blood_jbsat","blood_phys",
                  "blood_urban_dv","blood_imd","blood_hiqual",
                  "b_friends","b_partner","b_family","b_smoking","b_last7_drink","b_last7_maxunits",
                  "blood_ghq", "deltaghq",
                  "SSt_logOdds_CS", "SSt_logOdds_FS", "SSt_logOdds_NS", "SSt_probs_CS",
                  "SSt_probs_FS", "SSt_probs_NS", "SSt_PredictedSmokingStatus",
                  "SSc_smokingScore", "alc_exposure"
)

cols <- columns[columns %in% colnames(surv)]

surv <- surv[,cols]

surv$sex_categorical <- as.factor( surv$sex_categorical )
surv$blood_urban_dv <- as.factor( surv$blood_urban_dv )
surv$family_intact <- as.factor( surv$family_intact )
surv$father_employed <- as.factor( surv$father_employed )
surv$medtyp4 <- as.factor( surv$medtyp4 )
surv$blood_jbsat <- as.numeric( surv$blood_jbsat )
surv$wave <- as.numeric( surv$wave )

saveRDS( surv[surv$id %in% colnames(U1),], "mh_epigenetics_project/data/USM1survey_and_covar2.rds" )
saveRDS( surv[surv$id %in% colnames(U2),], "mh_epigenetics_project/data/USM2survey_and_covar2.rds" )
