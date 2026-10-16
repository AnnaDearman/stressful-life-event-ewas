library( dplyr ) 
library( haven ) 
library( tidyr ) 

survey_dir <- ""
genotype_dir <- ""
processed_dir <- ""
output_dir <- ""

# Read in main survey data file

df <- read_dta( paste0( survey_dir, "/ApplicationNumber_sendout.dta" ) )

# Read in and merge file containing genotype-based sex info

sex <- read.table( paste0( genotype_dir, "/ApplicationNumber_imputation.sample" ), 
                   header = T )[ , c( "ID1", "sex" )]
colnames( sex ) <- c( "id", "gtsex" )

df <- merge( df, sex, by = "id", all.x = T, all.y = F )

rm( sex ) 

# Make categorical sex variable if survey and genotype data agree

df$sex_categorical <- case_when(
  df$sex == 1 & df$gtsex == 1 ~ "male",
  df$sex == 2 & df$gtsex == 2 ~ "female" 
)

# Childhood circumstances

# Family disruption / family intact

df$family_disruption <- case_when( 
  # Convert inconsistent answers to missing first
  df$paju == 3 & df$lvag16 == 1  ~ NA, # Father died by 14 but respondent lived 
  # with father or both biological parents at 16 (error) 
  df$ynlp14 == 2 & df$lvag16 == 1  ~ NA, # Parent/s died by 14 
  # but respondent lived with both biological parents at 16 (error) 
  df$paju == 4 & df$lvag14 == 1 ~ NA, # Father recorded as "not living with"
  # in labour market status variable, but respondent living with both 
  # biological parents at 14 (error)
  df$paju == 3 & is.na( df$ynlp14 ) ~ "par_died", # Father recorded as deceased 
  # in labour market status variable, and no alternative reason given
  df$ynlp14 == 2 ~ "par_died",          # One or both parents died by age 14
  df$ynlp14 == 1 ~ "other",             # Bio parents never married/lived 
  # together by 14 so not moved/died 
  df$ynlp14 == 3 | df$paju == 4  ~ "par_moved", # At age 14, parents separated 
  # or divorced, OR father's labour market status is "not living with, so 
  # don't know" 
  df$lvag14 == 7 ~ "other",             # Living in care / foster home at 
  df$lvag16 == 1 & df$paju != 3 ~ "intact", # Living with both bio parents at 16 
  # and father's labour market status @ 14 isn't "deceased"
  df$lvag16 == 1 & is.na( df$paju ) ~ "intact", # Living with both bio parents 
  # at 16 and father's labour market status is missing
  !is.na( df$lvag16 ) ~ "other"         # Not living with both bio parents at 16
)

df$family_intact <- case_when(
  df$family_disruption == "intact" ~ 1, # Did not experience disruption 
  # from age 14-16 as far as we can tell
  !is.na( df$family_disruption ) ~ 0 # Captures all non-missing values
)

# Father's occupation: employed / unemployed

df$father_employed <- case_when( 
  # Convert inconsistent answers to missing first:
  # Father working or not working (but not deceased) at 14 but "why not living
  # with both at 14" == parent(s) died, and "who living with at 14" includes
  # "mother", i.e. father must have died, but isn't recorded as dead in paju
  df$paju %in% c( 1, 2 ) & df$ynlp14 == 2 & df$lvag14 %in% c( 3, 5 ) ~ NA, 
  df$paju == 3 & (df$lvag16 == 1 | df$lvag14 %in% c( 1, 4, 6 )) ~ NA, # Father 
  # died by 14 but respondent lived with father or both biological parents at 
  # 14 or 16 (error) 
  df$paju == 1 ~ 1,     # 1 if father known to have been employed
  !is.na( df$paju ) ~ 0 # Captures all non-missing values (including deceased
  # and don't know)
) 

# We based our measure of Stressful Life Events (SLEs) on Anderson et al (2022)
# https://doi.org/10.1177/00221465211055993
# with a few minor adjustments

# Some SLEs were derived using deriving_complex_sle_variables.R and are
# attached in the previously-imported file (now stored in R session as df).
# Other SLEs now need to be derived using this script.

# First, check whether any of the pre-derived SLEs have all-zeros

tmp <- df[ , colnames( df )[ grep( "died|premature|stillbirth", colnames( df ))] ]
lapply( tmp, table )

# Job loss

# "We included job losses attributed to being made redundant, dismissed/sacked, 
# or a temporary job ending" Anderson et al (2022)
colnames(df)[ grep( "stendreas", colnames(df) )]

# For waves where "reason for job loss" was collapsed into one variable
waves <- c( "b_", "c_", "d_", "e_", "f_" ) 
for( i in waves ){
  df[ , paste0( i, "jobloss" ) ] <- case_when(
    df[ , paste0( i, "stendreas" ) ] %in% c(3, 4, 5) ~ 1,
    df[ , paste0( i, "ivfio" ) ] %in% c(1, 2, 3) ~ 0
  )
}

# For waves where "reason for job loss" was split across variables by reason
waves <- c( "g_", "h_", "i_", "j_", "k_", "l_" ) 
for( i in waves ){
  df[ , paste0( i, "jobloss" ) ] <- case_when(
    df[ , paste0( i, "stendreas3" ) ] == 1 ~ 1,
    df[ , paste0( i, "stendreas4" ) ] == 1 ~ 1,
    df[ , paste0( i, "stendreas5" ) ] == 1 ~ 1,
    df[ , paste0( i, "ivfio" ) ] %in% c(1, 2, 3) ~ 0
  )
}

# Wave b event variables - job loss (a value of 24 indicates "threat of or
# actual")

df$b_jobloss2 <- case_when( df$b_event1 == 24 & df$b_event1s == 2 ~ 1,
                            df$b_event2 == 24 & df$b_event2s == 2 ~ 1,
                            df$b_event3 == 24 & df$b_event3s == 2 ~ 1,
                            TRUE ~ df$b_jobloss )

# Became unemployed

# "We included those who report their economic status as unemployed and did not 
# report being unemployed in the previous wave" Anderson et al (2022)
colnames(df)[ grep( "jbstat", colnames(df) )]
waves <- c( "ba_", "bb_", "bc_", "bd_", "be_", "bf_", "bg_", "bh_", "bi_",  
            "bj_", "bk_", "bl_", "bm_", "bn_", "bo_", "bp_", "bq_", "br_", 
            "a_","b_", "c_", "d_", "e_", "f_", 
            "g_", "h_", "i_", "j_", "k_", "l_" ) 

# For waves without "nmpsp_dv" variables, measure using current and previous
# wave variables
for( i in c(2:19) ){
  df[ , paste0( waves[i], "becameunemployed" ) ] <- case_when(
    # previous wave's "current" employment status wasn't unemployed
    df[ , paste0( waves[i-1], "jbstat" ) ] != 3 & 
    # but current status is unemployed
      df[ , paste0( waves[i], "jbstat" ) ] == 3 ~ 1,
    # completed both waves' interviews
    df[ , paste0( waves[i-1], "ivfio" ) ] %in% c(1, 2, 3) & 
      df[ , paste0( waves[i], "ivfio" ) ] %in% c(1, 2, 3) ~ 0
  )
}

# For waves with "nmpsp_dv" variables, can measure using current wave variables
# only
for( i in c(20:30) ){
  df[ , paste0( waves[i], "becameunemployed" ) ] <- case_when(
    # "since last interview, had at least one employment spell"
    df[ , paste0( waves[i], "nmpsp_dv" ) ] > 0 &
      df[ , paste0( waves[i], "jbstat" ) ] == 3 ~ 1,
    # No need to condition on available data at previous wave any more
      df[ , paste0( waves[i], "ivfio" ) ] %in% c(1, 2, 3) ~ 0
  )
}

# Combine jobloss and becameunemployed

waves <- c( "bb_", "bc_", "bd_", "be_", "bf_", "bg_", "bh_", "bi_",  
            "bj_", "bk_", "bl_", "bm_", "bn_", "bo_", "bp_", "bq_", "br_", 
            "a_","b_", "c_", "d_", "e_", "f_", 
            "g_", "h_", "i_", "j_", "k_", "l_" )
variables <- c("jobloss","becameunemployed")
for( i in waves ){
  cols <- colnames(df)[ colnames(df) %in% paste0(i,variables)]
  if( length(cols) == 1 ){
    df[ , paste0(i, "joblossorunemp") ] <- df[ , cols ]
  }
  if( length(cols) > 1){
    df[ , paste0(i, "joblossorunemp") ] <- rowSums(
      df[,cols], na.rm=TRUE
    )
  }
}

variables <- c("jobloss2","becameunemployed")
cols <- colnames(df)[ colnames(df) %in% paste0("b_",variables)]
df[ , "b_joblossorunemp2" ] <- rowSums(
  df[,cols], na.rm=TRUE
)

# Union dissolution

# "anyone who has at any point in the past year transitioned from being in a 
# marital or cohabiting relationship (lasting at least three months) to single 
# (i.e. not partnered, whether legally “separated,” “divorced,” or “never 
# married”)" Anderson et al (2022)

colnames(df)[ grep( "mastat", colnames(df) )]
colnames(df)[ grep( "mstatch", colnames(df) )]
colnames(df)[ grep( "cohab", colnames(df) )]
colnames(df)[ grep( "relup", colnames(df) )]

# For earlier waves, based on mastat only
waves <- c( "ba_", "bb_", "bc_", "bd_", "be_", "bf_", "bg_", "bh_", "bi_",  
            "bj_", "bk_", "bl_", "bm_", "bn_", "bo_", "bp_", "bq_", "br_", 
            "a_" ) 
for( i in c(2:19) ){
  
    # variable gets a "_dv" suffix from wave a onward
    prev <- colnames( df )[ grep( paste0( "^", waves[i-1], "mastat" ), colnames( df ))]
    curr <- colnames( df )[ grep( paste0( "^", waves[i], "mastat" ), colnames( df ))]
  
  df[ , paste0( waves[i], "uniondissolution" ) ] <- case_when(
    # Previous marital status was married, civil partnership or living as couple
    df[ , prev ] %in% c(2,3,10) &
      # Current marital status is single, separated or divorced
      df[ , curr ] %in% c(1,4,5) ~ 1,
    # Productive interview
    df[ , paste0( waves[i], "ivfio" ) ] %in% c(1, 2, 3)~ 0
  )
}

# For later waves, based on mastat plus other variables
waves <- c( "a_", "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_", "k_", "l_" ) 

for( i in c(2:12) ){
  
  df[ , paste0( waves[i], "uniondissolution" ) ] <- case_when(
    # Previous marital status was married, civil partnership or living as couple
    df[ , paste0( waves[i-1], "mastat_dv" ) ] %in% c(2,3,10) &
      # Current marital status is single, separated or divorced
      df[ , paste0( waves[i], "mastat_dv" ) ] %in% c(1,4,5) ~ 1,
    # or, marital status change since last interview (two variables) is either 
    # married or civil partnership, but 
    (df[ , paste0( waves[i], "mstatch1" ) ] %in% c(2,3) | 
      df[ , paste0( waves[i], "mstatch2" ) ] %in% c(2,3)) &
      # marital status is now single, separated or divorced
      df[ , paste0( waves[i], "mastat_dv" ) ] %in% c(1,4,5) ~ 1,
    # Cohabitation spell since last interview but now single, separated or div
    df[ , paste0( waves[i], "cohab" ) ] == 1 & 
      df[ , paste0( waves[i], "mastat_dv" ) ] %in% c(1,4,5) ~ 1,
    # Last wave partner is not resident
    df[ , paste0( waves[i], "relup" ) ] %in% c(3,4) ~ 1,
    # Productive interview
    df[ , paste0( waves[i], "ivfio" ) ] %in% c(1, 2, 3)~ 0
  )
}

# Forced move

# “evicted from rented accommodation/repossessed/other forced moves.” We 
# included those who indicated this reason." Anderson et al (2022)

colnames(df)[ grep( "movy", colnames(df) )]
colnames(df)[ grep( "housing", colnames(df) )]
waves <- c( "bc_", "bd_", "be_", "bf_", "bg_", "bh_", "bi_", "bj_", 
            "bk_", "bl_", "bm_", "bn_", "bo_", "bp_", "bq_", "br_" )

for( i in waves){
  df[ , paste0( i, "forcedmove" ) ] <- case_when(
    df[ , paste0( i, "movy1" ) ] == 15 ~ 1,
    df[ , paste0( i, "movy2" ) ] == 15 ~ 1,
    df[ , paste0( i, "ivfio" ) ] %in% c(1, 2, 3) ~ 0
  )
}

waves <- c( "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_", "k_", "l_" )

for( i in waves){
  df[ , paste0( i, "forcedmove" ) ] <- case_when(
    df[ , paste0( i, "housing" ) ] == 12 ~ 1,
    df[ , paste0( i, "ivfio" ) ] %in% c(1, 2, 3) ~ 0
  )
}

# Had child

# How many children have you had since [the previous interview]?” 
# “Did this pregnancy result in a live birth with a normal delivery or caesarean 
# section?”
colnames(df)[ grep( "nchild", colnames(df) )] # Males
colnames(df)[ grep( "pregout", colnames(df) )] # Females
colnames(df)[ grep( "nnewborn", colnames(df))] # This variable is what the original 
# authors used, but I have found that using the other variables gives 
# exactly the same values at most waves and four extra "yes" values at wave g

waves <- c( "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_", "k_", "l_" ) 

for( i in waves ){
  cols <- colnames(df)[grep( paste0( i, "pregout" ), colnames(df) )]
  if( length(cols) == 3 ){
    df[ , paste0( i, "hadchild" ) ] <- case_when(
      df[ , paste0( i, "nchild" ) ] > 0 ~ 1,
      df[ , paste0( i, "pregout1" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout2" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout3" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 4 ){
    df[ , paste0( i, "hadchild" ) ] <- case_when(
      df[ , paste0( i, "nchild" ) ] > 0 ~ 1,
      df[ , paste0( i, "pregout1" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout2" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout3" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout4" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 5 ){
    df[ , paste0( i, "hadchild" ) ] <- case_when(
      df[ , paste0( i, "nchild" ) ] > 0 ~ 1,
      df[ , paste0( i, "pregout1" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout2" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout3" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout4" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "pregout5" ) ] %in% c( 1, 2 ) ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
}

# Miscarriage, termination, ectopic pregnancy

# "We included indicated miscarriages, terminations, and ectopic or tubal 
# pregnancies in this category." Anderson et al (2022)
colnames(df)[ grep( "pregout", colnames(df) )]
colnames(df)[ grep( "pregend", colnames(df) )]
waves <- c( "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_", "k_", "l_" ) 

for( i in waves ){
  cols <- colnames(df)[grep( paste0( i, "pregend" ), colnames(df) )]
  if( length(cols) == 3 ){
    df[ , paste0( i, "misc_term_ectop" ) ] <- case_when(
      df[ , paste0( i, "pregend1" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend2" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend3" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 4 ){
    df[ , paste0( i, "misc_term_ectop" ) ] <- case_when(
      df[ , paste0( i, "pregend1" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend2" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend3" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend4" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 5 ){
    df[ , paste0( i, "misc_term_ectop" ) ] <- case_when(
      df[ , paste0( i, "pregend1" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend2" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend3" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend4" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "pregend5" ) ] %in% c( 1, 3, 4 ) ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
}

# Onset of serious health condition

# "We included new diagnoses of: 
# congestive heart failure, 3
# coronary heart disease, 4
# angina, 5
# heart attack or myocardial infarction, 6
# stroke, 7
# emphysema, 8
# hyperthyroidism, 9 
# hypothyroidism, 10
# chronic bronchitis, 11
# “any kind of liver condition,” 12
# cancer or malignancy, 13
# epilepsy 15" - based on Anderson et al (2022)

colnames(df)[ grep( "hcondn", colnames(df) )]
waves <- c( "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_" ) 
for( i in waves ){
  df[ , paste0( i, "onsethealthcond" ) ] <- case_when(
    df[ , paste0( i, "hcondn3" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn4" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn5" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn6" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn7" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn8" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn9" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn10" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn11" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn12" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn13" ) ] == 1 ~ 1,
    df[ , paste0( i, "hcondn15" ) ] == 1 ~ 1,
    df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
  )
}

# Hospital stay

# "individuals who had spent more than two days as a hospital inpatient for any 
# one new condition" Anderson et al (2022)
colnames(df)[ grep( "hospc", colnames(df) )]
colnames(df)[ grep( "hospdc", colnames(df) )]
waves <- c( "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_", "k_", "l_" ) 

for( i in waves ){
  cols <- colnames(df)[grep( paste0( i, "hospdc" ), colnames(df) )]
  print(length(cols))
  if( length(cols) == 5 ){ # 7, 8, 9, 10, 13
    df[ , paste0( i, "hospstay" ) ] <- case_when(
      df[ , paste0( i, "hospdc1" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc2" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc3" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc4" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc5" ) ] > 2 ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 7 ){ # 8, 9, 10, 13
    df[ , paste0( i, "hospstay" ) ] <- case_when(
      df[ , paste0( i, "hospdc1" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc2" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc3" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc4" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc5" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc6" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc7" ) ] > 2 ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 8 ){ # 9, 10, 13
    df[ , paste0( i, "hospstay" ) ] <- case_when(
      df[ , paste0( i, "hospdc1" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc2" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc3" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc4" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc5" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc6" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc7" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc8" ) ] > 2 ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 9 ){ # 10, 13
    df[ , paste0( i, "hospstay" ) ] <- case_when(
      df[ , paste0( i, "hospdc1" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc2" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc3" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc4" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc5" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc6" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc7" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc8" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc9" ) ] > 2 ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 10 ){ # 13
    df[ , paste0( i, "hospstay" ) ] <- case_when(
      df[ , paste0( i, "hospdc1" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc2" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc3" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc4" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc5" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc6" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc7" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc8" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc9" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc10" ) ] > 2 ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
  if( length(cols) == 13 ){ 
    df[ , paste0( i, "hospstay" ) ] <- case_when(
      df[ , paste0( i, "hospdc1" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc2" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc3" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc4" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc5" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc6" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc7" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc8" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc9" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc10" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc11" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc12" ) ] > 2 ~ 1,
      df[ , paste0( i, "hospdc13" ) ] > 2 ~ 1,
      df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
    )}
}

# Disability onset or deterioration

# "We included individuals reporting a long-standing condition and three or 
# more of these limitations (numbers 1–11) who also had not previously met these 
# criteria" Anderson et al (2022)
colnames(df)[ grep( "health", colnames(df) )] # If we condition on this being 
# "yes, I have a long standing condition" we seem to miss people at later
# waves, so we don't use it
colnames(df)[ grep( "disdif", colnames(df) )]
waves <- c( "a_", "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_", "k_", "l_" ) 

for( i in c(2:12) ){
  df[ , paste0( waves[i], "disab" ) ] <- case_when(
    df[ , paste0( waves[i], "ivfio" ) ] %in% c( 1, 2, 3 ) &
      df[ , paste0( waves[i-1], "ivfio" ) ] %in% c( 1, 2, 3 ) &
      rowSums( df[ , paste0( waves[i], "disdif", c(1:11)) ], na.rm=TRUE ) >= 3 &
      rowSums( df[ , paste0( waves[i-1], "disdif", c(1:11)) ], na.rm=TRUE ) < 3 ~ 1,
    df[ , paste0( waves[i], "ivfio" ) ] %in% c( 1, 2, 3 ) &
      df[ , paste0( waves[i-1], "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
  )
}

# Became full-time carer

# "We included individuals who were not a carer at all in the previous wave,
# but later reported spending 35 hours a week or more on such caring 
# responsibilities" Anderson et al (2022)
colnames(df)[ grep( "aidxhh", colnames(df) )]
colnames(df)[ grep( "aidhh", colnames(df) )]
colnames(df)[ grep( "aidhrs", colnames(df) )]
waves <- c( "bb_", "bc_", "bd_", "be_", "bf_", "bg_", "bh_", "bi_", "bj_", 
            "bk_", "bl_", "bm_", "bn_", "bo_", "bp_", "bq_", "br_", "a_",
            "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_", "k_", "l_" )

for( i in c(2:29) ){
  df[ , paste0( waves[i], "becameftcarer" ) ] <- case_when(
    df[ , paste0( waves[i], "aidhrs" ) ] %in% c(5, 6, 7) &
      df[ , paste0( waves[i-1], "aidhh" ) ] == 2 & 
      df[ , paste0( waves[i-1], "aidxhh" ) ] == 2 ~ 1,
    df[ , paste0( waves[i], "ivfio" ) ] %in% c( 1, 2, 3 ) &
      df[ , paste0( waves[i-1], "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
  )
}

# Household financial problems

# household reference person answered any of these items positively (i.e., 
# “Yes” to either of the first two questions 
# (Many people find it hard to keep up with their housing payments. In the last 
# twelve months, have you ever found yourself behind with your rent/mortgage?”; 
# "In the last twelve months, have you ever found yourself behind with paying 
# your Council Tax?”) 
# and/or, for the third, “Behind with some bills” or “Behind with all bills”). 
# - Anderson et al (2022)

colnames(df)[ grep( "xphsd", colnames(df) )]
waves <- c( 
  "bb_", "bc_", "bd_", "be_", "bf_", "bg_", "bh_", "bi_", "bj_",
  "bk_", "bl_", "bm_", "bn_", "bo_", "bp_", "bq_", "br_"
)

for( i in waves ){
  df[ , paste0( i, "hhfinancprob" )] <- case_when(
    df[ , paste0( i, "xphsdb" )] == 1 ~ 1, 
    df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
  )
}

waves <- c( "a_", "b_", "c_", "d_", "e_", "f_", 
            "g_", "h_", "i_", "j_", "k_", "l_" ) 
for( i in waves ){
  df[ , paste0( i, "hhfinancprob" )] <- case_when(
    df[ , paste0( i, "xphsdb" )] == 1 ~ 1, # problems paying for housing
    df[ , paste0( i, "xphsdct" )] == 1 ~ 1, # problems paying council tax
    df[ , paste0( i, "xphsdba" )] %in% c(2,3) ~ 1, # behind with some/all bills
    df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
  )
}

# Using the above SLE measures, we are missing some instances of widowhood,
# so we add those here based on variables mstatch1 and mstatch2
waves <- c( 
  "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_", "k_", "l_"
)

for( i in waves ){
  df[ , paste0( i, "0_in_hh_mem_died_but_spousedied" )] <- case_when(
    df[ , paste0( i, "mstatch1" )] == 6 ~ 1,
    df[ , paste0( i, "mstatch2" )] == 6 ~ 1,
    df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
  )
  df[ , paste0( i, "hh_mem_died_2" )] <- case_when(
    df[ , paste0( i, "0_in_hh_mem_died_but_spousedied" )] == 1 ~ 1,
    df[ , paste0( i, "hh_mem_died" )] == 1 ~1,
    df[ , paste0( i, "ivfio" ) ] %in% c( 1, 2, 3 ) ~ 0
  )
}

# Add up stressful life events

sle_variables_2 <- c( "hh_mem_died_2", "premature_or_lowbw", 
                    "stillbirth_or_infantdeath", "joblossorunemp",
                    "uniondissolution", "forcedmove", "hadchild", 
                    "misc_term_ectop", "onsethealthcond", "hospstay", "disab", 
                    "becameftcarer", "hhfinancprob" ) 
# We use the "_2" suffix because we modified the original authors'
# "household member died" approach by adding widowhood using 
# mstatch1 ==6 and mstatch2 ==6

waves <- c( "bb_", "bc_", "bd_", "be_", "bf_", "bg_", "bh_", "bi_",  
            "bj_", "bk_", "bl_", "bm_", "bn_", "bo_", "bp_", "bq_", "br_", 
            "a_","b_", "c_", "d_", "e_", "f_", 
            "g_", "h_", "i_", "j_", "k_", "l_" )
for( i in waves ){
  
  sle_cols <- colnames(df)[ colnames(df) %in% paste0(i,sle_variables_2)]
  if( length(sle_cols) == 1 ){
    df[ , paste0(i, "sles_2") ] <- df[ , sle_cols ]
  }
  if( length(sle_cols) > 1){
    df[ , paste0(i, "sles_2") ] <- rowSums( 
      df[,sle_cols], na.rm=TRUE
    )
  }
}

# N.B. the use of "blood" as a variable prefix indicates that the variable is
# as contemporaneous as possible with the blood collection. It was necessary
# to indicate this because some people had blood collected at wave 2 (b)
# and others at wave 3 (c) each of which had a separate main interview

df$blood_sles_2 <- case_when( 
  df$wave == 2 ~ df$b_sles_2,
  df$wave == 3 ~ df$c_sles_2
)

df$blood_sles_binary_2 <- case_when( # THIS IS OUR EWAS PHENOTYPE
  df$blood_sles_2 > 0 ~ 1,
  df$blood_sles_2 == 0 ~ 0
)

# Derive covariates

# Physical health

df$blood_phys <- case_when(
  df$wave == 2 ~ df$b_sf12pcs_dv,
  df$wave == 3 ~ df$c_sf12pcs_dv
)

# Area deprivation

df$blood_imd <- case_when(
  df$wave == 2 & !is.na(df$b_imd2015qe_dv) ~ df$b_imd2015qe_dv,
  df$wave == 3 & !is.na(df$c_imd2015qe_dv) ~ df$c_imd2015qe_dv,

  df$wave == 2 & !is.na(df$b_imd2014qw_dv) ~ df$b_imd2014qw_dv,
  df$wave == 3 & !is.na(df$c_imd2014qw_dv) ~ df$c_imd2014qw_dv,

  df$wave == 2 & !is.na(df$b_imd2016qs_dv) ~ df$b_imd2016qs_dv,
  df$wave == 3 & !is.na(df$c_imd2016qs_dv) ~ df$c_imd2016qs_dv,

  df$wave == 2 & !is.na(df$b_imd2017qni_dv) ~ df$b_imd2017qni_dv,
  df$wave == 3 & !is.na(df$c_imd2017qni_dv) ~ df$c_imd2017qni_dv
)

waves <- c( "a_", "b_", "c_", "d_", "e_", "f_", "g_", "h_",
            "i_", "j_", "k_", "l_" )

for( i in c(1:12) ){
  df[ , paste0( waves[i], "imd1417" ) ] <- case_when(
    !is.na( df[ , paste0( waves[i], "imd2015qe_dv" ) ] ) ~ df[ , paste0( waves[i], "imd2015qe_dv" ) ],
    !is.na( df[ , paste0( waves[i], "imd2016qs_dv" ) ] ) ~ df[ , paste0( waves[i], "imd2016qs_dv" ) ],
    !is.na( df[ , paste0( waves[i], "imd2014qw_dv" ) ] ) ~ df[ , paste0( waves[i], "imd2014qw_dv" ) ],
    !is.na( df[ , paste0( waves[i], "imd2017qni_dv" ) ] ) ~ df[ , paste0( waves[i], "imd2017qni_dv" ) ],
  )
  
  df[ , paste0( waves[i], "imd1920" ) ] <- case_when(
    !is.na( df[ , paste0( waves[i], "imd2019qe_dv" ) ] ) ~ df[ , paste0( waves[i], "imd2019qe_dv" ) ],
    !is.na( df[ , paste0( waves[i], "imd2020qs_dv" ) ] ) ~ df[ , paste0( waves[i], "imd2020qs_dv" ) ],
    !is.na( df[ , paste0( waves[i], "imd2019qw_dv" ) ] ) ~ df[ , paste0( waves[i], "imd2019qw_dv" ) ]
  )
}

# Social support

colnames( df )[ grep( "rely|undstnd|openup", colnames( df ) ) ]

# waves b, e, k only

waves <- c( "b_", "e_", "k_" )
for( i in waves ){
  df[,paste0(i,"friends")] <- as.numeric( scale( rowSums(
  df[,c(paste0(i,"scfopenup"), paste0(i,"scfrely"), paste0(i,"scfundstnd"))]
))) * -1
df[,paste0(i,"partner")] <- as.numeric( scale( rowSums(
  df[,c(paste0(i,"scpopenup"), paste0(i,"scprely"), paste0(i,"scpundstnd"))]
))) * -1
df[,paste0(i,"family")] <- as.numeric( scale( rowSums(
  df[,c(paste0(i,"scropenup"), paste0(i,"scrrely"), paste0(i,"scrundstnd"))]
))) * -1
df[,paste0(i,"socialsupport")] <- as.numeric( scale( rowSums(
  df[,c(paste0(i,"friends"), paste0(i,"partner"), paste0(i,"family"))], na.rm = T
)))
}

# Adiposity (could affect EWAS signal) - use waist-to-height ratio

df$whr <- df$wstval / df$htval

# Educational attainment

df$blood_hiqual <- as.numeric( case_when(
  df$wave == 2 ~ df$b_hiqual_dv,
  df$wave == 3 ~ df$c_hiqual_dv
) )

# Amend partner social support so that those without a partner have 
# the median value

df$b_partner <- case_when(
  !is.na( df$b_partner ) ~ df$b_partner,
  df$b_mastat_dv %in% c(0,1,4,5,6,7) ~ median( df$b_partner, na.rm=T )
)

# Job satisfaction

df$blood_jbsat <- case_when(
  df$wave == 2 & !is.na(df$b_jbsat) ~ df$b_jbsat,
  df$wave == 3 & !is.na(df$c_jbsat)~ df$c_jbsat,
  TRUE ~ 4 # neutral response, used for anyone without a job
)
df$blood_jbsat <- as.numeric(
    scale( df$blood_jbsat, center=T, scale=T )
)	

# Urban vs rural

df$blood_urban_dv <- case_when(
  df$wave == 2 ~ df$b_urban_dv,
  df$wave == 3 ~ df$c_urban_dv
)

# Health beh: smoking (only used for exploratory plotting of smoking episcore)

df$b_smoking <- case_when(
  df$b_smever == 2 ~ 0,
  df$b_smever == 1 & df$b_smnow == 2 & df$b_smcigs == 3 ~ 0,
  df$b_smever == 1 & df$b_smnow == 2 & df$b_smcigs == 2 ~ 0,
  df$b_smever == 1 & df$b_smnow == 2 & df$b_smcigs == 1 ~ 1,
  df$b_ncigs %in% c(1:10) ~ 2,
  df$b_ncigs %in% c(11:20) ~ 3,
  df$b_ncigs > 20 ~ 4
)

# Health beh: drinking (only used for exploratory plotting of episcore)

df$b_last7_drink <- case_when(
  df$b_scalcl7d == 2 ~ 0,
  df$b_scnalcl7d == 1 | df$b_scnalcl7d == 2 ~ 1,
  df$b_scnalcl7d == 3 | df$b_scnalcl7d == 4 ~ 2,
  df$b_scnalcl7d > 4 ~ 3
)

df$answeredalcquestions <- rowSums(
  is.na(df[,c("b_scnalcpint","b_scnalcshot","b_scnalcwine","b_scnalcpops")])
)
df$b_scnalcpint[ is.na(df$b_scnalcpint) ] <- 0
df$b_scnalcshot[ is.na(df$b_scnalcshot) ] <- 0
df$b_scnalcwine[ is.na(df$b_scnalcwine) ] <- 0
df$b_scnalcpops[ is.na(df$b_scnalcpops) ] <- 0
df[df$answeredalcquestions == 4, c("b_scnalcpint","b_scnalcshot","b_scnalcwine","b_scnalcpops")] <- NA

df$b_last7_maxunits <- df$b_scnalcpint*2.5 + df$b_scnalcshot*1 + 
  df$b_scnalcwine*2.1 + df$b_scnalcpops*1.5
df$b_last7_maxunits[ is.na( df$b_last7_drink ) ] <- NA

df$blood_hh_mem_died_2 <- case_when(
  df$wave==2 ~ df$b_hh_mem_died_2,
  df$wave==3 ~ df$c_hh_mem_died_2
)
df$blood_misc_term_ectop <- case_when(
  df$wave==2 ~ df$b_misc_term_ectop,
  df$wave==3 ~ df$c_misc_term_ectop
)
df$blood_stillbirth_or_infant_death <- case_when(
  df$wave==2 ~ df$b_stillbirth_or_infant_death,
  df$wave==3 ~ df$c_stillbirth_or_infant_death
)
df$blood_uniondissolution <- case_when(
  df$wave==2 ~ df$b_uniondissolution,
  df$wave==3 ~ df$c_uniondissolution
)
df$blood_joblossorunemp <- case_when(
  df$wave==2 ~ df$b_joblossorunemp,
  df$wave==3 ~ df$c_joblossorunemp
)
df$blood_hhfinancprob <- case_when(
  df$wave==2 ~ df$b_hhfinancprob,
  df$wave==3 ~ df$c_hhfinancprob
)
df$blood_forcedmove <- case_when(
  df$wave==2 ~ df$b_forcedmove,
  df$wave==3 ~ df$c_forcedmove
)
df$blood_hadchild <- case_when(
  df$wave==2 ~ df$b_hadchild,
  df$wave==3 ~ df$c_hadchild
)
df$blood_premature_or_lowbw <- case_when(
  df$wave==2 ~ df$b_premature_or_lowbw,
  df$wave==3 ~ df$c_premature_or_lowbw
)
df$blood_becameftcarer <- case_when(
  df$wave==2 ~ df$b_becameftcarer,
  df$wave==3 ~ df$c_becameftcarer
)
df$blood_onsethealthcond <- case_when(
  df$wave==2 ~ df$b_onsethealthcond,
  df$wave==3 ~ df$c_onsethealthcond
)
df$blood_hospstay <- case_when(
  df$wave==2 ~ df$b_hospstay,
  df$wave==3 ~ df$c_hospstay
)
df$blood_disab <- case_when(
  df$wave==2 ~ df$b_disab,
  df$wave==3 ~ df$c_disab
)

# For looking at subtypes of SLE

df$blood_sle_personalloss <- case_when(
  df$blood_hh_mem_died_2 == 1 ~ 1, 
  df$blood_misc_term_ectop == 1 ~ 1,  
  df$blood_stillbirth_or_infant_death == 1 ~ 1, 
  df$blood_uniondissolution == 1 ~ 1, 
  TRUE ~ 0
)
df$blood_sle_socioeconomicloss <- case_when(
  df$blood_joblossorunemp == 1 ~ 1,  
  df$blood_hhfinancprob == 1 ~ 1, 
  df$blood_forcedmove == 1 ~ 1, 
  TRUE ~ 0
)
df$blood_sle_caring <- case_when(
  df$blood_hadchild == 1 ~ 1, 
  df$blood_premature_or_lowbw == 1 ~ 1,  
  df$blood_becameftcarer == 1 ~ 1,
  TRUE ~ 0
)
df$blood_sle_health <- case_when(
  df$blood_onsethealthcond == 1 ~ 1, 
  df$blood_hospstay == 1 ~ 1,  
  df$blood_disab == 1 ~ 1,
  TRUE ~ 0
)

# calculate all interview dates

df$a_intdat <- as.POSIXct( as.Date( paste(
  df$a_intdaty_dv, df$a_intdatm_dv, df$a_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$b_intdat <- as.POSIXct( as.Date( paste(
  df$b_intdaty_dv, df$b_intdatm_dv, df$b_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$c_intdat <- as.POSIXct( as.Date( paste(
  df$c_intdaty_dv, df$c_intdatm_dv, df$c_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$d_intdat <- as.POSIXct( as.Date( paste(
  df$d_intdaty_dv, df$d_intdatm_dv, df$d_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$e_intdat <- as.POSIXct( as.Date( paste(
  df$e_intdaty_dv, df$e_intdatm_dv, df$e_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$f_intdat <- as.POSIXct( as.Date( paste(
  df$f_intdaty_dv, df$f_intdatm_dv, df$f_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$g_intdat <- as.POSIXct( as.Date( paste(
  df$g_intdaty_dv, df$g_intdatm_dv, df$g_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$h_intdat <- as.POSIXct( as.Date( paste(
  df$h_intdaty_dv, df$h_intdatm_dv, df$h_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$i_intdat <- as.POSIXct( as.Date( paste(
  df$i_intdaty_dv, df$i_intdatm_dv, df$i_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$j_intdat <- as.POSIXct( as.Date( paste(
  df$j_intdaty_dv, df$j_intdatm_dv, df$j_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$k_intdat <- as.POSIXct( as.Date( paste(
  df$k_intdaty_dv, df$k_intdatm_dv, df$k_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))
df$l_intdat <- as.POSIXct( as.Date( paste(
  df$l_intdaty_dv, df$l_intdatm_dv, df$l_intdatd_dv, sep="-" ),
  format="%Y-%m-%d" ))

# Calculate delay between blood collection date and lab processing date

df$nursedat <- as.POSIXct( as.Date( paste(
  df$nurdayy, df$nurdaym, df$nurdayd, sep="-" ), format="%Y-%m-%d" ))
df$labdat <- as.POSIXct( as.Date( paste(
  df$dory, df$dorm, df$dord, sep="-" ), format="%Y-%m-%d" ))
df$blooddelay <- as.numeric( as.Date( df$labdat ) - as.Date( df$nursedat ))
df$blooddelay[ df$blooddelay < 1 ] <- NA
df$blooddelay_2 <- df$blooddelay ^ 2

df$blood_ghq <- as.numeric( case_when(
  df$wave == 2 ~ df$b_scghq1_dv,
  df$wave == 3 ~ df$c_scghq1_dv
  ) )

df$blood_sf12 <- as.numeric( case_when(
  df$wave == 2 ~ df$b_sf12mcs_dv,
  df$wave == 3 ~ df$c_sf12mcs_dv
  ) )

df$blood_ghq <- as.numeric( scale( df$blood_ghq, center = T, scale = T ) )
df$blood_sf12 <- as.numeric( scale( df$blood_sf12, center = T, scale = T ) )
df$blood_phys <- as.numeric( scale( df$blood_phys, center = T, scale = T ) )
df$deltaghq <- case_when(
  df$wave == 2 ~ df$b_scghq1_dv - df$a_scghq1_dv,
  df$wave == 3 ~ df$c_scghq1_dv - df$b_scghq1_dv,
)

# Merge with epigenetic data covariates

tmp <- read_dta( "./rawdata/epigenetic_data/ApplicationNumber_epi_sendout.dta" )
tmp2 <- read_dta( "./rawdata/epigenetic_data/ApplicationNumber_epi_2020_sendout.dta" )

tmp <- rbind( tmp[,c("id","confage","nsex","rackbarcode","cd8t","cd4t","nk","bcell","mono","gran")], 
              tmp2[,c("id","confage","nsex","rackbarcode","cd8t","cd4t","nk","bcell","mono","gran")] )

df <- merge( df, tmp, by="id", all.x=T, all.y=T )

df$confage2 <- df$confage ^ 2
df$age <- as.numeric( scale( df$confage, scale = T, center = T ))
df$age2 <- as.numeric( scale( df$confage2, scale = T, center = T ))

saveRDS( df, paste0( processed_dir, "/UKHLS_SLEs.rds" ) )

