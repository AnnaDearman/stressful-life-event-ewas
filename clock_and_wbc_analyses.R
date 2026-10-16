library( broom.mixed )
library( dplyr )
library( lme4 )
library( lmerTest )

# Read in and combine data files

df <- readRDS( "clocks.rds" )
usm1 <- readRDS( "USM1survey_and_covar2.rds" )
usm2 <- readRDS( "USM2survey_and_covar2.rds" )

df2 <- rbind( usm1, usm2 )

df <- merge( df, df2, by="id", all.x=T, all.y=F )

class( df$nsex ) <- "integer"

sex <- read.table( "ApplicationNumber_imputation.sample",
                   header = T )[ , c( "ID1", "sex" )]
colnames( sex ) <- c( "id", "gtsex" )

df <- merge( df, sex, by = "id", all.x = T, all.y = F )

df <- subset( df, df$gtsex==df$nsex )
df <- subset( df, !is.na(sex_categorical) )

df$sex <- case_when( df$nsex == 1 ~ "male",
                     df$nsex == 2 ~ "female" )

rm( usm1, usm2, df2, sex )

# Prepare data and useful vectors

cells <- c( "aCD4Tnv", "aBaso", "aCD4Tmem", "aBmem", "aBnv", "aTreg", "aCD8Tmem", 
            "aCD8Tnv", "aEos", "aNK", "aMono", "aNeu" ) # Adult cells for use as outcomes
covar <- c( "age", "age2", "sex", "blooddelay", "blooddelay_2", "blood_ghq", 
            "blood_imd", "blood_urban_dv", "family_intact", "father_employed",
            "blood_hiqual", "whr", "blood_phys", "medtyp4", "b_family", 
            "b_friends", "b_partner", "blood_jbsat", cells, "PCA1", "PCA2", 
            "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8", "PCA9", "PCA10",
            "SSc_smokingScore", "alc_exposure" )
clocks <- c( "horvath", "hannum", "levine", "horvath_sb", "lin", "poam38",
             "pace" )
cat <- c( "blood_sles_binary_2", "rackbarcode", "sex", "blood_urban_dv", 
             "family_intact", "father_employed", "medtyp4" ) # categorical
cont <- c( clocks, "age", "age2", "blooddelay", "blood_ghq", "blooddelay_2",
           "blood_imd", "blood_hiqual", "whr", "blood_phys", "b_family", 
           "b_friends", "b_partner", "blood_jbsat", "b_smoking",
           "b_last7_drink", cells, "PCA1", "PCA2", "PCA3", "PCA4", "PCA5",
           "PCA6", "PCA7", "PCA8", "PCA9", "PCA10",
           "SSc_smokingScore", "alc_exposure" ) # continuous

# Ensure variables are of correct type

for( i in cat ){
  df[,i] <- as.factor( df[,i] )
}

for( i in cont ){
  df[,i] <- as.numeric( scale( df[,i], center = T, scale = T ))
}

##### Prepare results dataframe #####

results <- data.frame(
  effect=vector(mode="character"),
  group=vector(mode="character"),
  term=vector(mode="character"),
  estimate=vector(mode="numeric"),
  std.error=vector(mode="numeric"),
  statistic=vector(mode="numeric"),
  df=vector(mode="numeric"),
  p.value=vector(mode="numeric"),
  sample=vector(mode="character"),
  covar=vector(mode="character"),
  outcome=vector(mode="character")
)

##### Leukocyte analyses ####

# Prepare models

vars1 <- c( "age", "age2", "sex", "blooddelay", "blooddelay_2" )
vars2 <- c( vars1, "blood_imd", "blood_urban_dv", "family_intact", 
           "father_employed", "blood_hiqual",
           "PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8",
           "PCA9", "PCA10" )
vars3 <- c( vars2, "whr", "blood_phys", "medtyp4", "SSc_smokingScore",
           "alc_exposure" )
vars4 <- c( vars3, "b_family", "b_friends", "b_partner", "blood_jbsat" )

# NB rackbarcode = batch (random effect)

model_1 <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + ",
  paste0( vars1, collapse=" + ") )
model_2 <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + ",
  paste0( vars2, collapse=" + ") )
model_3 <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + ",
  paste0( vars3, collapse=" + ") )
model_4 <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + ",
  paste0( vars4, collapse=" + "),
  " + blood_sles_binary_2*b_family + blood_sles_binary_2*b_friends + blood_sles_binary_2*b_partner + blood_sles_binary_2*blood_jbsat" )

# Sensitivity analysis model (includes psychological distress,
# which is what "_pd" suffix indicates)

vars1_pd <- c( "age", "age2", "sex", "blooddelay", "blooddelay_2", "blood_ghq" )
vars2_pd <- c( vars1_pd, "blood_imd", "blood_urban_dv", "family_intact",
              "father_employed", "blood_hiqual",
              "PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8",
              "PCA9", "PCA10" )
vars3_pd <- c( vars2_pd, "whr", "blood_phys", "medtyp4", "SSc_smokingScore",
              "alc_exposure" )
vars4_pd <- c( vars3_pd, "b_family", "b_friends", "b_partner", "blood_jbsat" )

model_1_pd <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + blood_sles_binary_2*blood_ghq + ",
  paste0( vars1_pd, collapse=" + "))
model_2_pd <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + blood_sles_binary_2*blood_ghq + ",
  paste0( vars2_pd, collapse=" + "))
model_3_pd <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + blood_sles_binary_2*blood_ghq + ",
  paste0( vars3_pd, collapse=" + "))
model_4_pd <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + blood_sles_binary_2*blood_ghq + ",
  paste0( vars4_pd, collapse=" + "),
  " + blood_sles_binary_2*b_family + blood_sles_binary_2*b_friends + blood_sles_binary_2*b_partner + blood_sles_binary_2*blood_jbsat" )

# Run leukocyte analyses and append to results table

for( i in cells ){
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1 ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "full"
  tmp$covar <- model_1
  tmp$outcome <- i
  results <- rbind( results, tmp )

  tmp <- tidy( lmer( as.formula( paste( i, model_1 ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_1
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1 ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_1
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1_pd ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ),]
  tmp$sample <- "full"
  tmp$covar <- model_1_pd
  tmp$outcome <- i
  results <- rbind( results,tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1_pd ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_1_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1_pd ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_1_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  ###
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2 ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "full"
  tmp$covar <- model_2
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2 ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_2
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2 ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_2
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2_pd ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ),]
  tmp$sample <- "full"
  tmp$covar <- model_2_pd
  tmp$outcome <- i
  results <- rbind( results,tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2_pd ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_2_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2_pd ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_2_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  ###
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3 ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "full"
  tmp$covar <- model_3
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3 ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_3
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3 ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_3
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3_pd ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ),]
  tmp$sample <- "full"
  tmp$covar <- model_3_pd
  tmp$outcome <- i
  results <- rbind( results,tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3_pd ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_3_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3_pd ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_3_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  ###
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4 ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "full"
  tmp$covar <- model_4
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4 ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_4
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4 ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_4
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4_pd ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ),]
  tmp$sample <- "full"
  tmp$covar <- model_4_pd
  tmp$outcome <- i
  results <- rbind( results,tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4_pd ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_4_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4_pd ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_4_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
}

##### Clock analyses ####

# Prepare models
# NB for clock analyses, 18/19 cell types are included as covariates (including 
# estimated cord blood cells)

vars1 <- c( "age", "age2", "sex", "blooddelay", "blooddelay_2", 
           "aCD4Tnv", "aBaso", "aCD4Tmem", "aBmem", "aBnv", "aTreg", "aCD8Tmem",
           "aCD8Tnv", "aEos", "aNK", "aMono", "B", "CD4T", "CD8T", "Mono", 
           "nRBC", "Gran", "NK" )
vars2 <- c( vars1, "blood_imd", "blood_urban_dv", "family_intact", 
           "father_employed", "blood_hiqual",
           "PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8",
           "PCA9", "PCA10" )
vars3 <- c( vars2, "whr", "blood_phys", "medtyp4", "SSc_smokingScore",
           "alc_exposure" )
vars4 <- c( vars3, "b_family", "b_friends", "b_partner", "blood_jbsat" )

model_1 <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + ",
  paste0( vars1, collapse=" + "))
model_2 <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + ",
  paste0( vars2, collapse=" + "))
model_3 <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + ",
  paste0( vars3, collapse=" + "))
model_4 <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + ",
  paste0( vars4, collapse=" + "),
  " + blood_sles_binary_2*b_family + blood_sles_binary_2*b_friends + blood_sles_binary_2*b_partner + blood_sles_binary_2*blood_jbsat")

# Sensitivity analysis model (includes psychological distress,
# which is what "_pd" suffix indicates)

vars1_pd <- c( "age", "age2", "sex", "blooddelay", "blooddelay_2", "blood_ghq", 
              "aCD4Tnv", "aBaso", "aCD4Tmem", "aBmem", "aBnv", "aTreg", "aCD8Tmem", 
              "aCD8Tnv", "aEos", "aNK", "aMono", "B", "CD4T", "CD8T", "Mono", 
              "nRBC", "Gran", "NK" )
vars2_pd <- c( vars1_pd, "blood_imd", "blood_urban_dv", "family_intact",
              "father_employed", "blood_hiqual",
              "PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8",
              "PCA9", "PCA10" )
vars3_pd <- c( vars2_pd, "whr", "blood_phys", "medtyp4", "SSc_smokingScore",
              "alc_exposure" )
vars4_pd <- c( vars3_pd, "b_family", "b_friends", "b_partner", "blood_jbsat" )

model_1_pd <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + blood_sles_binary_2*blood_ghq + ",
  paste0( vars1_pd, collapse=" + "))
model_2_pd <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + blood_sles_binary_2*blood_ghq + ",
  paste0( vars2_pd, collapse=" + "))
model_3_pd <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + blood_sles_binary_2*blood_ghq + ",
  paste0( vars3_pd, collapse=" + "))
model_4_pd <- paste0( " ~ blood_sles_binary_2 + (1|rackbarcode) + blood_sles_binary_2*blood_ghq + ",
  paste0( vars4_pd, collapse=" + "),
  " + blood_sles_binary_2*b_family + blood_sles_binary_2*b_friends + blood_sles_binary_2*b_partner + blood_sles_binary_2*blood_jbsat")

# Run clock analyses and append to results table

for( i in clocks ){
  
  # NB: for horvath_sb, a warning message is produced:
  # "boundary (singular) fit: see help('isSingular')" for
  # model 2 USM 1 sensitivity analysis
  # model 3 USM 1
  # model 4 USM 1
  # model 4 USM 1 sensitivity analysis

  tmp <- tidy( lmer( as.formula( paste( i, model_1 ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "full"
  tmp$covar <- model_1
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1 ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_1
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1 ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_1
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1_pd ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ),]
  tmp$sample <- "full"
  tmp$covar <- model_1_pd
  tmp$outcome <- i
  results <- rbind( results,tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1_pd ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_1_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_1_pd ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_1_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  ###
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2 ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "full"
  tmp$covar <- model_2
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2 ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_2
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2 ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_2
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2_pd ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ),]
  tmp$sample <- "full"
  tmp$covar <- model_2_pd
  tmp$outcome <- i
  results <- rbind( results,tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2_pd ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_2_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_2_pd ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_2_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  ###
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3 ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "full"
  tmp$covar <- model_3
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3 ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_3
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3 ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_3
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3_pd ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ),]
  tmp$sample <- "full"
  tmp$covar <- model_3_pd
  tmp$outcome <- i
  results <- rbind( results,tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3_pd ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_3_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_3_pd ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_3_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  ###
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4 ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "full"
  tmp$covar <- model_4
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4 ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_4
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4 ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_4
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4_pd ) ), df ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ),]
  tmp$sample <- "full"
  tmp$covar <- model_4_pd
  tmp$outcome <- i
  results <- rbind( results,tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4_pd ) ), subset( df, usm==1 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM1"
  tmp$covar <- model_4_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
  tmp <- tidy( lmer( as.formula( paste( i, model_4_pd ) ), subset( df, usm==2 ) ) )
  tmp <- tmp[ grepl( "blood_sles_binary|ghq", tmp$term ), ]
  tmp$sample <- "USM2"
  tmp$covar <- model_4_pd
  tmp$outcome <- i
  results <- rbind( results, tmp )
  
}

write.csv( results, 
           "cell_and_clock_analysis_results.csv",
           quote = F, row.names = F)
