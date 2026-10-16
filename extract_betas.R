# Script to be run from grandparent directory, relative to "scripts"

library( bigmelon )
library( DescTools )
USM1 <- openfn.gds( "ApplicationNumber.gds" )
USM2 <- openfn.gds( "ApplicationNumber_2020.gds" )
U1 <- betas( USM1 )[ , ]
U2 <- betas( USM2 )[ , ]
mccartney <- lapply( list.files( "[directory name]", pattern="S221359", full.names=T ), function( x ) read.table( x, header=T ))
toexclude <- unique( c( mccartney[[1]]$IlmnID, mccartney[[2]][,1], mccartney[[3]][,1] ) )
U1 <- U1[!rownames(U1) %in% toexclude,]
U2 <- U2[!rownames(U2) %in% toexclude,]
genotypedriven_aries <- read.csv("/storage/projects/USoc_omics/genotypedriven.csv") 
genotypedriven_godmc <- read.csv("/storage/projects/USoc_omics/assoc_meta_all_probes_unique.csv") #
genotypedriven <- unique( c( genotypedriven_aries$x, 
                             genotypedriven_godmc$cpg 
                             ) )
U1 <- U1[!rownames(U1) %in% genotypedriven,]
U2 <- U2[!rownames(U2) %in% genotypedriven,]

U1w <- t(apply( U1, 1, Winsorize ))
U2w <- t(apply( U2, 1, Winsorize ))
saveRDS(U1w,"mh_epigenetics_project/data/USM1betamatrixWinsorized.rds")
saveRDS(U2w,"mh_epigenetics_project/data/USM2betamatrixWinsorized.rds")
