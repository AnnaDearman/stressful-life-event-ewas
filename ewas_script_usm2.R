# Code based on https://github.com/psychiatric-genomics-consortium/mdd-mwas

rm(list = ls())

library( limma )
library( IlluminaHumanMethylationEPICanno.ilm10b2.hg19 )

usm <- readRDS( "USM2survey_and_covar2.rds" )

bet <- readRDS( "USM2betamatrixWinsorized.rds" )

usm <- usm[ usm$id %in% colnames( bet ) , ]

print(dim(usm))

print(dim(bet))

mod_covars <- list( 
 # Basic model
 model1 = c( "age", "age2", "sex_categorical", "blooddelay", "blooddelay_2", "B", "CD4T", "CD8T", "Mono", "nRBC", "Gran", "NK", "aCD4Tnv", "aBaso", "aCD4Tmem", "aBmem", "aBnv", "aTreg", "aCD8Tmem", "aCD8Tnv", "aEos", "aNK", "aMono" ), 
 # Adjusting for demographic confounders 
 model2 = c( "age", "age2", "sex_categorical", "blooddelay", "blooddelay_2", "B", "CD4T", "CD8T", "Mono", "nRBC", "Gran", "NK", "aCD4Tnv", "aBaso", "aCD4Tmem", "aBmem", "aBnv", "aTreg", "aCD8Tmem", "aCD8Tnv", "aEos", "aNK", "aMono", "blood_imd", "blood_urban_dv", "PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8", "PCA9", "PCA10", "family_intact", "father_employed", "blood_hiqual" ),
 # Adjusting for potential mediators
 model3 = c( "age", "age2", "sex_categorical", "blooddelay", "blooddelay_2", "B", "CD4T", "CD8T", "Mono", "nRBC", "Gran", "NK", "aCD4Tnv", "aBaso", "aCD4Tmem", "aBmem", "aBnv", "aTreg", "aCD8Tmem", "aCD8Tnv", "aEos", "aNK", "aMono", "blood_imd", "blood_urban_dv", "PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8", "PCA9", "PCA10", "family_intact", "father_employed", "blood_hiqual", "SSc_smokingScore", "whr", "blood_phys", "medtyp4", "alc_exposure"),
 # Adjusting for chronic stressors / protective factors
 model4 = c( "age", "age2", "sex_categorical", "blooddelay", "blooddelay_2", "B", "CD4T", "CD8T", "Mono", "nRBC", "Gran", "NK", "aCD4Tnv", "aBaso", "aCD4Tmem", "aBmem", "aBnv", "aTreg", "aCD8Tmem", "aCD8Tnv", "aEos", "aNK", "aMono", "blood_imd", "blood_urban_dv", "PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8", "PCA9", "PCA10", "family_intact", "father_employed", "blood_hiqual", "SSc_smokingScore", "whr", "blood_phys", "medtyp4", "alc_exposure", "b_family", "b_friends", "b_partner", "blood_jbsat" )
)

anno <- getAnnotation( IlluminaHumanMethylationEPICanno.ilm10b2.hg19 )
genes <- data.frame( ID=anno$Name, geneSymbol=anno$UCSC_RefGene_Name, CHR=anno$chr, MAPINFO=anno$pos, FEATURE=anno$UCSC_RefGene_Group, CpGISLAND=anno$Relation_to_Island )

pheno <- "blood_sles_binary_2"

for( i in 1:length(mod_covars) ){
        covars <- mod_covars[[i]]
        tmp <- na.omit( usm[,c(pheno, "id", "rackbarcode", covars)])
        batches <- as.factor(tmp$rackbarcode)
        covariates <- tmp[ , covars ]
        bets <- bet[ , colnames( bet ) %in% tmp$id ]
        corfit <- duplicateCorrelation( bets, design=model.matrix(~1, data = covariates), block=batches )
        if( i == 4 ){
            covars <- c( covars, "b_family*blood_sles_binary_2", "b_friends*blood_sles_binary_2", "b_partner*blood_sles_binary_2", "blood_jbsat*blood_sles_binary_2" )
        }
        mods <- model.matrix(
            as.formula(
                paste0(" ~ ", paste0( c( pheno, covars ), collapse=" + ") ) 
                ),
                data=tmp
        )
        fit <- lmFit( bets, mods, block=batches, correlation=corfit$consensus.correlation )
        fit$n <- nrow( tmp )
        
        efit <- eBayes( fit )
        TT <- topTable( efit, coef=2, adjust='fdr', number=length(fit$Amean))
        TT$ID <- rownames(TT)
        TT_anno <- merge( genes, TT, by="ID", all.y=TRUE )
        TT_anno$df <- fit$df.residual
        fits_n <- data.frame(ID=names(fit$N),  N=fit$N)
        pheno_betas <- data.frame(ID=row.names(fit$coefficients), beta=fit$coefficients[,2])
        TT_linear <- TT_anno

        TT_linear$beta <- fit$coefficients[as.character(TT_linear$ID),pheno]
        TT_linear$se <- TT_linear$beta / TT_linear$t
        TT_linear$N <- fit$n
        efit_betas <- plyr::adply(efit$coefficients, 1, function(x) data.frame(param=names(x), beta=x))
        names(efit_betas)[1] <- 'ID'
        efit_t <- plyr::adply(efit$t, 1, function(x) data.frame(param=names(x), t=x))
        names(efit_t)[1] <- 'ID'
        efit_p <- plyr::adply(efit$p.value, 1, function(x) data.frame(param=names(x), p=x))
        names(efit_p)[1] <- 'ID'
        efit_linear <- efit_betas
        efit_linear$t <- efit_t$t
        efit_linear$p <- efit_p$p
        efit_linear$se <- efit_linear$beta / efit_linear$t

        write.table(TT_linear, paste0("usm2_model_", i, "of4_toptable.txt"), sep="\t", row.names=FALSE, quote=F)

        write.table(efit_linear[,c('ID', 'param', 'beta', 'se', 't', 'p')], paste0("usm2_model_", i, "of4_linear.txt"), sep="\t", row.names=FALSE, quote=F)
}
