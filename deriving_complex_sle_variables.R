library( dplyr )
library( haven )

# Fill in the paths to the directories where the xwave and indresp files are saved, 
# and fill in the directory for saving output file

bhps_directory <- ""
ukhls_directory <- ""
output_directory <- ""

# Read in the xwave file and retain only individuals with genetic data

df <- subset( 
  read_dta( 
    paste0( ukhls_directory, "/xwavedat.dta" )
    )[,c("pidp","genetics","epigenetics")],
  genetics == 1
)

##### HOUSEHOLD MEMBER DIED #####

# Read in the xwave file again, but this time we are interested in the
# deceased, who we will subsequently identify as relatives of those with 
# genetic data (or not)

xw <- subset( 
  read_dta( 
    paste0( ukhls_directory, "/xwavedat.dta" )
  )[,c("pidp","dcsedw_dv")],
  dcsedw_dv > 0 # "Wave deceased" will be a value greater than 0 if respondent
) # is recorded as deceased

# GET HIDPs (household IDs)

# a) make vector of waves

waves <- c( paste0( "b", letters[1:18] ), # BHPS waves ba-br
            letters[1:12] # UKHLS waves a-l
)

# b) store the hidps of the genetics subsample, AND
# c) store the hidps of the deceased individuals

for( i in 1:30 ){
  if( nchar( waves[i] ) == 2 ){ # If the wave is a BHPS wave (i.e. two characters)
    tmp <- read_dta( # Read in the corresponsing indresp file (pidp and hidp only)
      paste0( bhps_directory, "/", waves[i], "_indresp.dta" ))[ , c( "pidp", paste0( waves[i], "_hidp" ) )]
    df <- merge( df, tmp, by = "pidp", all.x = T, all.y = F ) # Merge with the genetics subsample
    xw <- merge( xw, tmp, by = "pidp", all.x = T, all.y = F ) # Merge with the deceased sample
  }
  if( nchar( waves[i] ) == 1 ){ # Repeat for UKHLS waves (i.e. one character)
    tmp <- read_dta(
      paste0( ukhls_directory, "/", waves[i], "_indresp.dta" ))[ , c( "pidp", paste0( waves[i], "_hidp" ) )]
    df <- merge( df, tmp, by = "pidp", all.x = T, all.y = F )
    xw <- merge( xw, tmp, by = "pidp", all.x = T, all.y = F )
  }
}

# d) Identify the last known hidps of the deceased individuals, 
# i.e. hidp at the previous wave, e.g. if dcsedw_dv == 23 (UKHLS e), get d_hidp

xw$hidp_at_death <- case_when(
  xw$dcsedw_dv == 2 ~ xw$ba_hidp, xw$dcsedw_dv == 3 ~ xw$bb_hidp,
  xw$dcsedw_dv == 4 ~ xw$bc_hidp, xw$dcsedw_dv == 5 ~ xw$bd_hidp,
  xw$dcsedw_dv == 6 ~ xw$be_hidp, xw$dcsedw_dv == 7 ~ xw$bf_hidp,
  xw$dcsedw_dv == 8 ~ xw$bg_hidp, xw$dcsedw_dv == 9 ~ xw$bh_hidp,
  xw$dcsedw_dv == 10 ~ xw$bi_hidp, xw$dcsedw_dv == 11 ~ xw$bj_hidp,
  xw$dcsedw_dv == 12 ~ xw$bk_hidp, xw$dcsedw_dv == 13 ~ xw$bl_hidp,
  xw$dcsedw_dv == 14 ~ xw$bm_hidp, xw$dcsedw_dv == 15 ~ xw$bn_hidp,
  xw$dcsedw_dv == 16 ~ xw$bo_hidp, xw$dcsedw_dv == 17 ~ xw$bp_hidp,
  xw$dcsedw_dv == 18 ~ xw$bq_hidp, xw$dcsedw_dv == 19 ~ xw$br_hidp,
  xw$dcsedw_dv == 20 ~ xw$a_hidp, xw$dcsedw_dv == 21 ~ xw$b_hidp,
  xw$dcsedw_dv == 22 ~ xw$c_hidp, xw$dcsedw_dv == 23 ~ xw$d_hidp,
  xw$dcsedw_dv == 24 ~ xw$e_hidp, xw$dcsedw_dv == 25 ~ xw$f_hidp,
  xw$dcsedw_dv == 26 ~ xw$g_hidp, xw$dcsedw_dv == 27 ~ xw$h_hidp,
  xw$dcsedw_dv == 28 ~ xw$i_hidp, xw$dcsedw_dv == 29 ~ xw$j_hidp,
  xw$dcsedw_dv == 30 ~ xw$k_hidp, xw$dcsedw_dv == 31 ~ xw$l_hidp
)

# Identify the bereaved members of the genetics subsample using this
# information

for( i in 2:30 ){ # Going through the waves at which people died (numeric value)
  wave <- waves[i-1] # Identify the wave (letter value) at which they were last 
  # alive (i.e. previous wave)
  bereaved_households <- na.omit( 
    unique( 
      subset( xw, 
              dcsedw_dv == i )$hidp_at_death )) # For a given wave, grab the 
  # hidps of the subsequently-deceased (i.e. "bereaved households")
  
  # For the genetics subsample, identify those whose hidp indicates a "bereaved
  # household", at the correct wave
  
  # Variable name needs to contain the wave they were reported to be deceased
  # i.e. for which there is no available hidp
  
  df[ , paste0( waves[i], "_hh_mem_died" ) ] <- case_when(
    
    # The hidp relating to the event of death/bereavement needs to belong 
    # to the previous wave
    
    df[ , paste0( waves[i-1], "_hidp" ) ] %in% bereaved_households ~ 1, 
    # 1 = someone in the household died since the last wave
    TRUE ~ 0 
  )
}

##### PREGNANCY-RELATED EVENTS #####

for( i in 20:30){
  wave <- waves[i]
  tmp <- read_dta( paste0( ukhls_directory, "/", wave, "_newborn.dta" ) ) # pidp == mother
  
  # Had a baby which was premature or had a low birth weight
  
  tmp$premature_or_lowbw <- case_when(
    tmp[ , paste0( wave, "_bwtel" ) ] == 1 & # child born early
      tmp[ , paste0( wave, "_bwtwk" ) ] %in% c(4:20) ~ 1, # how many weeks early
    tmp[ , paste0( wave, "_bwt" ) ] == 1 & # birthweight units == lb and oz
      tmp[ , paste0( wave, "_bwtlb" ) ] < 5 ~ 1, # less than 5 lb
    tmp[ , paste0( wave, "_bwt" ) ] == 1 & # birthweight units == lb and oz
      tmp[ , paste0( wave, "_bwtlb" ) ] == 5 & # 5lb and...
      tmp[ , paste0( wave, "_bwtoz" ) ] <= 8 ~ 1, # ...1-8 oz
    tmp[ , paste0( wave, "_bwt" ) ] == 2 & # birthweight units == kg
      tmp[ , paste0( wave, "_bwtk" ) ] <= 2.5 ~ 1, # 2.5kg
    tmp[ , paste0( wave, "_bwtg5" ) ] == 2 ~ 1, # birthweight was not more than 5.5lb
    
    TRUE ~ 0
  
  )

  df[ , paste0( wave, "_premature_or_lowbw" ) ] <- case_when(
    df$pidp %in% subset( tmp, premature_or_lowbw == 1 )$pidp ~ 1,
    TRUE ~ 0
  )
  
  # Lost a baby to stillbirth or infant death
  
  tmp2 <- read_dta( paste0( ukhls_directory, "/", wave, "_indresp.dta" ) )
  # Keep only pregnancy outcome variables (the number of variables varies by 
  # wave, depending on the maximum number of pregnancies experienced by someone
  # during that wave)
  tmp2 <- tmp2[ , c( "pidp", colnames(tmp2)[ grep( "pregend", colnames(tmp2))])]
  
  tmp2$stillbirth_or_infant_death <- rowSums( # Add up the number of "TRUE"s 
    # returned when querying whether one of the pregnancies ended in a stillbirth
    tmp2[,c(colnames(tmp2)[grep("pregend",colnames(tmp2))])] == 2,
    na.rm=T
  )
  
  tmp$stillbirth_or_infant_death <- case_when(
    tmp[ , paste0( wave, "_lchlv" ) ] > 2 ~ 1 # "Child still lives with parent?" died (3) or stillborn (4)
  )
  
  df[ , paste0( wave, "_stillbirth_or_infant_death" ) ] <- case_when(
    df$pidp %in% subset( tmp, stillbirth_or_infant_death == 1 )$pidp ~ 1,
    df$pidp %in% subset( tmp2, stillbirth_or_infant_death >= 1 )$pidp ~ 1,
    TRUE ~ 0
  )

}

columns <- colnames(df)[grep("died|premature|stillbirth",colnames(df))]

write_dta( df[ , c("pidp", columns )], paste0( output_directory, "/derived_variables.dta" ) )
