library( missMethyl )

input_dir <- ""
output_dir <- ""

# MAIN ANALYSES #####

# Main effects of SLE (read in files and print lambda) #####

toptables1 <- lapply( list.files(
  input_dir,
  pattern="usm1_model_[1-4]of4.*toptable", full.names = TRUE ),
  function(x) read.table( x, header=T, sep="\t", quote=""))

toptables1 <- lapply( toptables1, function(x) x[ order( x$P.Value ), ])

toptables2 <- lapply( list.files( 
  input_dir,
  pattern="usm2_model_[1-4]of4.*toptable", full.names = TRUE ),
  function(x) read.table( x, header=T, sep="\t", quote=""))

toptables2 <- lapply( toptables2, function(x) x[ order( x$P.Value ), ])

# Calculate lambda (code influenced by "Genomics for Social Scientists" course
# run by University of Michigan) https://genomicsworkshop.isr.umich.edu/

for( i in 1:length(toptables1) ){
  observed <- -log10(sort(toptables1[[i]]$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(toptables1[[i]]$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
}

for( i in 1:length(toptables2) ){
  observed <- -log10(sort(toptables2[[i]]$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(toptables2[[i]]$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
}

# Interaction terms (SLE*moderator) (read in files and print lambda) #####

toptablesintterms <- lapply( list.files(
  input_dir,
  pattern="model_4of4_interactions", full.names = TRUE ),
  function(x) read.table( x, header=F, sep="\t", quote=""))

names(toptablesintterms) <- list.files(
  input_dir,
  pattern="model_4of4_interactions", full.names = FALSE )
for( i in names(toptablesintterms) ){
  colnames(toptablesintterms[[i]]) <- c(
    "ID", "term", "beta", "se", "t", "P.Value"
  )
}

toptablesintterms$usm1_model_4of4_interactions.txt <- merge(
  toptablesintterms$usm1_model_4of4_interactions.txt,
  toptables1[[1]][,c("ID","geneSymbol","CHR","MAPINFO","FEATURE","CpGISLAND")], 
  by="ID", all.x=T, all.y=F
)
toptablesintterms$usm2_model_4of4_interactions.txt <- merge(
  toptablesintterms$usm2_model_4of4_interactions.txt,
  toptables2[[1]][,c("ID","geneSymbol","CHR","MAPINFO","FEATURE","CpGISLAND")], 
  by="ID", all.x=T, all.y=F
)

toptablesintterms <- lapply( toptablesintterms, function(x) x[ order( x$P.Value ), ])

for( i in 1:length(toptablesintterms) ){
  tmp <- subset(toptablesintterms[[i]],term=="blood_sles_binary_2:b_family")
  observed <- -log10(sort(tmp$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(tmp$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
  tmp <- subset(toptablesintterms[[i]],term=="blood_sles_binary_2:b_friends")
  observed <- -log10(sort(tmp$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(tmp$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
  tmp <- subset(toptablesintterms[[i]],term=="blood_sles_binary_2:b_partner")
  observed <- -log10(sort(tmp$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(tmp$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
  tmp <- subset(toptablesintterms[[i]],term=="blood_sles_binary_2:blood_jbsat")
  observed <- -log10(sort(tmp$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(tmp$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
}

# QQ plots #####

qqman::qq(toptables1[[1]]$P.Value, main="USM1, SLE main effect, model 1")
qqman::qq(toptables1[[2]]$P.Value, main="USM1, SLE main effect, model 2")
qqman::qq(toptables1[[3]]$P.Value, main="USM1, SLE main effect, model 3")
qqman::qq(toptables1[[4]]$P.Value, main="USM1, SLE main effect, model 4" )

qqman::qq(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                 term=="blood_sles_binary_2:b_family")$P.Value,
          main="USM1, SLE*family-support interaction")
qqman::qq(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                 term=="blood_sles_binary_2:b_friends")$P.Value,
          main="USM1, SLE*friends-support interaction")
qqman::qq(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                 term=="blood_sles_binary_2:b_partner")$P.Value,
          main="USM1, SLE*partner-support interaction")
qqman::qq(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                 term=="blood_sles_binary_2:blood_jbsat")$P.Value,
          main="USM1, SLE*job-satisfaction interaction")

qqman::qq(toptables2[[1]]$P.Value, main="USM2, SLE main effect, model 1")
qqman::qq(toptables2[[2]]$P.Value, main="USM2, SLE main effect, model 2")
qqman::qq(toptables2[[3]]$P.Value, main="USM2, SLE main effect, model 3")
qqman::qq(toptables2[[4]]$P.Value, main="USM2, SLE main effect, model 4" )

qqman::qq(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                 term=="blood_sles_binary_2:b_family")$P.Value,
          main="USM2, SLE*family-support interaction")
qqman::qq(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                 term=="blood_sles_binary_2:b_friends")$P.Value,
          main="USM2, SLE*friends-support interaction")
qqman::qq(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                 term=="blood_sles_binary_2:b_partner")$P.Value,
          main="USM2, SLE*partner-support interaction")
qqman::qq(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                 term=="blood_sles_binary_2:blood_jbsat")$P.Value,
          main="USM2, SLE*job-satisfaction interaction")

# Look at top hits ####

# Main effects of SLE

for( i in 1:length(toptables1) ){
  print(paste0("model ",i," SLE USM 1 n experiment-wide sig: ",nrow(subset(toptables1[[i]], P.Value < 1.345479e-7 )))) 
  print(paste0("model ",i," SLE USM 1 n FDR sig: ",nrow(subset(toptables1[[i]], adj.P.Val < 0.05)))) 
  write.csv( subset(toptables1[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm1_tophits0_001_model_", i, "of4.csv" ), row.names = FALSE, quote = FALSE
  )
}

for( i in 1:length(toptables2) ){
  print(paste0("model ",i," SLE USM 2 n experiment-wide sig: ",nrow(subset(toptables2[[i]], P.Value < 1.338326e-7)))) 
  print(paste0("model ",i," SLE USM 2 n FDR sig: ",nrow(subset(toptables2[[i]], adj.P.Val < 0.05))))
  write.csv( subset(toptables2[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm2_tophits0_001_model_", i, "of4.csv" ), row.names = FALSE, quote = FALSE
  )
}

# Interaction terms (SLE*moderator) in model 4

print(paste0("SLE*family USM 1 n experiment-wide sig: ",
             nrow(subset(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_family"), P.Value < 1.345479e-7))))
write.csv( subset(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_family"), P.Value < 0.001), paste0(
                           output_dir, "/usm1_SLEfamilyinteraction_tophits0_001_model_4of4.csv" ),
           row.names = FALSE, quote = FALSE
)
print(paste0("SLE*friends USM 1 n experiment-wide sig: ",
             nrow(subset(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_friends"), P.Value < 1.345479e-7))))
write.csv( subset(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_friends"), P.Value < 0.001), paste0(
                           output_dir, "/usm1_SLEfriendsinteraction_tophits0_001_model_4of4.csv" ),
           row.names = FALSE, quote = FALSE
)
print(paste0("SLE*partner USM 1 n experiment-wide sig: ",
             nrow(subset(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_partner"), P.Value < 1.345479e-7))))
write.csv( subset(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_partner"), P.Value < 0.001), paste0(
                           output_dir, "/usm1_SLEpartnerinteraction_tophits0_001_model_4of4.csv" ),
           row.names = FALSE, quote = FALSE
)
print(paste0("SLE*jbsat USM 1 n experiment-wide sig: ",
             nrow(subset(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:blood_jbsat"), P.Value < 1.345479e-7))))
write.csv( subset(subset(toptablesintterms$usm1_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:blood_jbsat"), P.Value < 0.001), paste0(
                           output_dir, "/usm1_SLEjobsatinteraction_tophits0_001_model_4of4.csv" ),
           row.names = FALSE, quote = FALSE
)
print(paste0("SLE*family USM 2 n experiment-wide sig: ",
             nrow(subset(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_family"), P.Value < 1.338326e-7))))
write.csv( subset(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_family"), P.Value < 0.001), paste0(
                           output_dir, "/usm2_SLEfamilyinteraction_tophits0_001_model_4of4.csv" ),
           row.names = FALSE, quote = FALSE
)
print(paste0("SLE*friends USM 2 n experiment-wide sig: ",
      nrow(subset(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_friends"), P.Value < 1.338326e-7))))
write.csv( subset(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_friends"), P.Value < 0.001), paste0(
                           output_dir, "/usm2_SLEfriendsinteraction_tophits0_001_model_4of4.csv" ),
           row.names = FALSE, quote = FALSE
)
print(paste0("SLE*partner USM 2 n experiment-wide sig: ",
             nrow(subset(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_partner"), P.Value < 1.338326e-7))))
write.csv( subset(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:b_partner"), P.Value < 0.001), paste0(
                           output_dir, "/usm2_SLEpartnerinteraction_tophits0_001_model_4of4.csv" ),
           row.names = FALSE, quote = FALSE
)
print(paste0("SLE*jbsat USM 2 n experiment-wide sig: ",
             nrow(subset(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:blood_jbsat"), P.Value < 1.338326e-7))))
write.csv( subset(subset(toptablesintterms$usm2_model_4of4_interactions.txt,
                         term=="blood_sles_binary_2:blood_jbsat"), P.Value < 0.001), paste0(
                           output_dir, "/usm2_SLEjobsatinteraction_tophits0_001_model_4of4.csv" ),
           row.names = FALSE, quote = FALSE
)

# Enrichment analyses for model 4 #####

# Top 50 for SLE, top 50 for each moderator, hence 250 probes (unless overlap)
# per dataset

tmp <- c(
  head( toptables1[[4]], 50 )$ID,
  head( subset( toptablesintterms[[1]],
                term=="blood_sles_binary_2:b_family" ), 50 )$ID,
  head( subset( toptablesintterms[[1]],
                term=="blood_sles_binary_2:b_friends" ), 50 )$ID,
  head( subset( toptablesintterms[[1]],
                term=="blood_sles_binary_2:b_partner" ), 50 )$ID,
  head( subset( toptablesintterms[[1]],
                term=="blood_sles_binary_2:blood_jbsat" ), 50 )$ID
)
en <- gometh( sig.cpg=tmp, 
              all.cpg=toptables1[[1]]$ID, 
              collection="KEGG", 
              array.type="EPIC", 
              plot.bias=FALSE,
              prior.prob=TRUE, 
              equiv.cpg=TRUE, 
              fract.counts=TRUE,
              genomic.features="ALL", 
              sig.genes=TRUE )
en <- en[order(en$P.DE),]
write.csv( en, paste0( output_dir, "/USM1_SLEandmods_250_KEGG.csv" ),
           quote=T, row.names = F)
en <- gometh( sig.cpg=tmp, 
              all.cpg=toptables1[[1]]$ID,  
              collection="GO", 
              array.type="EPIC", 
              plot.bias=FALSE,
              prior.prob=TRUE, 
              equiv.cpg=TRUE, 
              fract.counts=TRUE,
              genomic.features="ALL", 
              sig.genes=TRUE )
en <- en[order(en$P.DE),]
write.csv( en, paste0( output_dir, "/USM1_SLEandmods_250_GO.csv" ),
           quote=T, row.names = F)

tmp <- c(
  head( toptables2[[4]], 50 )$ID,
  head( subset( toptablesintterms[[2]],
                term=="blood_sles_binary_2:b_family" ), 50 )$ID,
  head( subset( toptablesintterms[[2]],
                term=="blood_sles_binary_2:b_friends" ), 50 )$ID,
  head( subset( toptablesintterms[[2]],
                term=="blood_sles_binary_2:b_partner" ), 50 )$ID,
  head( subset( toptablesintterms[[2]],
                term=="blood_sles_binary_2:blood_jbsat" ), 50 )$ID
)
en <- gometh( sig.cpg=tmp, 
              all.cpg=toptables2[[1]]$ID, 
              collection="KEGG", 
              array.type="EPIC", 
              plot.bias=FALSE,
              prior.prob=TRUE, 
              equiv.cpg=TRUE, 
              fract.counts=TRUE,
              genomic.features="ALL", 
              sig.genes=TRUE )
en <- en[order(en$P.DE),]
write.csv( en, paste0( output_dir, "/USM2_SLEandmods_250_KEGG.csv" ),
           quote=T, row.names = F)
en <- gometh( sig.cpg=tmp, 
              all.cpg=toptables2[[1]]$ID,  
              collection="GO", 
              array.type="EPIC", 
              plot.bias=FALSE,
              prior.prob=TRUE, 
              equiv.cpg=TRUE, 
              fract.counts=TRUE,
              genomic.features="ALL", 
              sig.genes=TRUE )
en <- en[order(en$P.DE),]
write.csv( en, paste0( output_dir, "/USM2_SLEandmods_250_GO.csv" ),
           quote=T, row.names = F)

rm( toptables1, toptables2, toptablesintterms, en )

# SENSITIVITY ANALYSES #####

# Main effects of SLE (read in files and print lambda) #####

toptables1 <- lapply( list.files(
  input_dir,
  pattern="usm1_model_sle-int-ghq.*of4_toptable", full.names = TRUE ),
  function(x) read.table( x, header=T, sep="\t", quote=""))

toptables1 <- lapply( toptables1, function(x) x[ order( x$P.Value ), ])

toptables2 <- lapply( list.files( 
  input_dir,
  pattern="usm2_model_sle-int-ghq.*of4_toptable", full.names = TRUE ),
  function(x) read.table( x, header=T, sep="\t", quote=""))

toptables2 <- lapply( toptables2, function(x) x[ order( x$P.Value ), ])

# Calculate lambda (code influenced by "Genomics for Social Scientists" course
# run by University of Michigan) https://genomicsworkshop.isr.umich.edu/

for( i in 1:length(toptables1) ){
  observed <- -log10(sort(toptables1[[i]]$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(toptables1[[i]]$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
}

for( i in 1:length(toptables2) ){
  observed <- -log10(sort(toptables2[[i]]$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(toptables2[[i]]$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
}

# Interaction terms (read in files and print lambda) #####

toptablesintterms1 <- lapply( list.files(
  input_dir,
  pattern="usm1_model_sle.*of4_interactions", full.names = TRUE ),
  function(x) read.table( x, header=F, sep="\t", quote=""))

names(toptablesintterms1) <- list.files(
  input_dir,
  pattern="usm1_model_sle.*of4_interactions", full.names = FALSE )

for( i in names(toptablesintterms1) ){
  colnames(toptablesintterms1[[i]]) <- c(
    "ID", "term", "beta", "se", "t", "P.Value"
  )
}

toptablesintterms2 <- lapply( list.files(
  input_dir,
  pattern="usm2_model_sle.*of4_interactions", full.names = TRUE ),
  function(x) read.table( x, header=F, sep="\t", quote=""))

names(toptablesintterms2) <- list.files(
  input_dir,
  pattern="usm2_model_sle.*of4_interactions", full.names = FALSE )

for( i in names(toptablesintterms2) ){
  colnames(toptablesintterms2[[i]]) <- c(
    "ID", "term", "beta", "se", "t", "P.Value"
  )
}

for( i in names(toptablesintterms1) ){
  merge(toptablesintterms1[[i]],
        toptables1[[1]][,c("ID","geneSymbol","CHR","MAPINFO","FEATURE","CpGISLAND")], 
        by="ID", all.x=T, all.y=F
  )
}
for( i in names(toptablesintterms2) ){
  merge(toptablesintterms2[[i]],
        toptables2[[1]][,c("ID","geneSymbol","CHR","MAPINFO","FEATURE","CpGISLAND")], 
        by="ID", all.x=T, all.y=F
  )
}

toptablesintterms1 <- lapply( toptablesintterms1, function(x) x[ order( x$P.Value ), ])
toptablesintterms2 <- lapply( toptablesintterms2, function(x) x[ order( x$P.Value ), ])

# Calculate lambda for SLE*GHQ interactions

for( i in 1:length(toptablesintterms1) ){
  observed <- -log10(sort(subset(toptablesintterms1[[i]],term=="blood_sles_binary_2:blood_ghq")$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(subset(toptablesintterms1[[i]],term=="blood_sles_binary_2:blood_ghq")$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
}

for( i in 1:length(toptablesintterms2) ){
  observed <- -log10(sort(subset(toptablesintterms2[[i]],term=="blood_sles_binary_2:blood_ghq")$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(subset(toptablesintterms2[[i]],term=="blood_sles_binary_2:blood_ghq")$P.Value)))
  lambda <- median(observed) / median(expected)
  print(round(lambda,3))
}

# Calculate lambda for SLE*moderator interaction terms

for( i in c("blood_sles_binary_2:b_family","blood_sles_binary_2:b_friends",
            "blood_sles_binary_2:b_partner","blood_sles_binary_2:blood_jbsat") ){
  observed <- -log10(sort(subset(toptablesintterms1$`usm1_model_sle-int-ghq_4of4_interactions.txt`,
                                 term==i)$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(subset(
    toptablesintterms1$`usm1_model_sle-int-ghq_4of4_interactions.txt`,
    term==i)$P.Value)))
  lambda <- median(observed) / median(expected)
  print(paste("USM1 model 4",i,"lambda",round(lambda,3)))
  
  observed <- -log10(sort(subset(toptablesintterms2$`usm2_model_sle-int-ghq_4of4_interactions.txt`,
                                 term==i)$P.Value, decreasing = F))
  expected <- -log10(ppoints(length(subset(
    toptablesintterms2$`usm2_model_sle-int-ghq_4of4_interactions.txt`,
    term==i)$P.Value)))
  lambda <- median(observed) / median(expected)
  print(paste("USM2 model 4",i,"lambda",round(lambda,3)))
}

# Separate lists and delete combined lists

toptablesghqintterms1 <- lapply( toptablesintterms1, function(x) subset(
  x, term=="blood_sles_binary_2:blood_ghq"))

toptablesfamilyintterms1 <- lapply( toptablesintterms1, function(x) subset(
  x, term=="blood_sles_binary_2:b_family"))

toptablesfriendsintterms1 <- lapply( toptablesintterms1, function(x) subset(
  x, term=="blood_sles_binary_2:b_friends"))

toptablespartnerintterms1 <- lapply( toptablesintterms1, function(x) subset(
  x, term=="blood_sles_binary_2:b_partner"))

toptablesjbsatintterms1 <- lapply( toptablesintterms1, function(x) subset(
  x, term=="blood_sles_binary_2:blood_jbsat"))

toptablesghqintterms2 <- lapply( toptablesintterms2, function(x) subset(
  x, term=="blood_sles_binary_2:blood_ghq"))

toptablesfamilyintterms2 <- lapply( toptablesintterms2, function(x) subset(
  x, term=="blood_sles_binary_2:b_family"))

toptablesfriendsintterms2 <- lapply( toptablesintterms2, function(x) subset(
  x, term=="blood_sles_binary_2:b_friends"))

toptablespartnerintterms2 <- lapply( toptablesintterms2, function(x) subset(
  x, term=="blood_sles_binary_2:b_partner"))

toptablesjbsatintterms2 <- lapply( toptablesintterms2, function(x) subset(
  x, term=="blood_sles_binary_2:blood_jbsat"))

rm(toptablesintterms1,toptablesintterms2)

# QQ plots #####

qqman::qq(toptables1[[1]]$P.Value, main="USM1, SLE main effect, model 1 ")
qqman::qq(toptables1[[2]]$P.Value, main="USM1, SLE main effect, model 2")
qqman::qq(toptables1[[3]]$P.Value, main="USM1, SLE main effect, model 3")
qqman::qq(toptables1[[4]]$P.Value, main="USM1, SLE main effect, model 4")

qqman::qq(toptables2[[1]]$P.Value, main="USM2, SLE main effect, model 1")
qqman::qq(toptables2[[2]]$P.Value, main="USM2, SLE main effect, model 2")
qqman::qq(toptables2[[3]]$P.Value, main="USM2, SLE main effect, model 3")
qqman::qq(toptables2[[4]]$P.Value, main="USM2, SLE main effect, model 4")

qqman::qq(toptablesghqintterms1[[1]]$P.Value,
          main="USM1, SLE*GHQ interaction term, model 1")
qqman::qq(toptablesghqintterms1[[2]]$P.Value,
          main="USM1, SLE*GHQ interaction term, model 2")
qqman::qq(toptablesghqintterms1[[3]]$P.Value,
          main="USM1, SLE*GHQ interaction term, model 3")
qqman::qq(toptablesghqintterms1[[4]]$P.Value,
          main="USM1, SLE*GHQ interaction term, model 4")

qqman::qq(toptablesghqintterms2[[1]]$P.Value,
          main="USM2, SLE*GHQ interaction term, model 1")
qqman::qq(toptablesghqintterms2[[2]]$P.Value,
          main="USM2, SLE*GHQ interaction term, model 2")
qqman::qq(toptablesghqintterms2[[3]]$P.Value,
          main="USM2, SLE*GHQ interaction term, model 3")
qqman::qq(toptablesghqintterms2[[4]]$P.Value,
          main="USM2, SLE*GHQ interaction term, model 4")

qqman::qq(toptablesfamilyintterms1[[4]]$P.Value,
          main="USM1, SLE*family-support interaction")
qqman::qq(toptablesfriendsintterms1[[4]]$P.Value,
          main="USM1, SLE*friends-support interaction")
qqman::qq(toptablespartnerintterms1[[4]]$P.Value,
          main="USM1, SLE*partner-support interaction")
qqman::qq(toptablesjbsatintterms1[[4]]$P.Value,
          main="USM1, SLE*job-satisfaction interaction")

qqman::qq(toptablesfamilyintterms2[[4]]$P.Value,
          main="USM2, SLE*family-support interaction")
qqman::qq(toptablesfriendsintterms2[[4]]$P.Value,
          main="USM2, SLE*friends-support interaction")
qqman::qq(toptablespartnerintterms2[[4]]$P.Value,
          main="USM2, SLE*partner-support interaction")
qqman::qq(toptablesjbsatintterms2[[4]]$P.Value,
          main="USM2, SLE*job-satisfaction interaction")

# Look at top hits #### AMEND FROM HERE DOWN

# Main effects of SLE

for( i in 1:length(toptables1) ){
  print(paste0("model ",i," SLE USM 1 n experiment-wide sig: ",nrow(subset(toptables1[[i]], P.Value < 1.345479e-7 )))) 
  print(paste0("model ",i," SLE USM 1 n FDR sig: ",nrow(subset(toptables1[[i]], adj.P.Val < 0.05)))) 
  write.csv( subset(toptables1[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm1_topSLEhits_ghqinttophits0_001_model_", i, "of4.csv" ), row.names = FALSE, quote = FALSE
  )
}

for( i in 1:length(toptables2) ){
  print(paste0("model ",i," SLE USM 2 n experiment-wide sig: ",nrow(subset(toptables2[[i]], P.Value < 1.338326e-7)))) 
  print(paste0("model ",i," SLE USM 2 n FDR sig: ",nrow(subset(toptables2[[i]], adj.P.Val < 0.05))))
  write.csv( subset(toptables2[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm2_topSLEhits_ghqinttophits0_001_model_", i, "of4.csv" ), row.names = FALSE, quote = FALSE
  )
}

# Interaction terms 

for( i in 1:length(toptablesghqintterms1) ){
  print(paste0("SLE*GHQ-12 USM 1 n experiment-wide sig: ",
               nrow(subset(toptablesghqintterms1[[i]], P.Value < 1.345479e-7 ))))
  write.csv( subset(toptablesghqintterms1[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm1_SLEGHQinteraction_tophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
} 

for( i in 1:length(toptablesghqintterms2) ){
  print(paste0("SLE*GHQ-12 USM 2 n experiment-wide sig: ",
               nrow(subset(toptablesghqintterms2[[i]], P.Value < 1.338326e-7 ))))
  write.csv( subset(toptablesghqintterms2[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm2_SLEGHQinteraction_tophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
}                                   

for( i in 1:length(toptablesfamilyintterms1) ){
  print(paste0("SLE*family USM 1 n experiment-wide sig (sens): ",
               nrow(subset(toptablesfamilyintterms1[[i]], P.Value < 1.345479e-7 ))))
  write.csv( subset(toptablesfamilyintterms1[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm1_SLEfamilyinteraction_ghqinttophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
} 
for( i in 1:length(toptablesfamilyintterms2) ){
  print(paste0("SLE*family USM 2 n experiment-wide sig (sens): ",
               nrow(subset(toptablesfamilyintterms2[[i]], P.Value < 1.338326e-7 ))))
  write.csv( subset(toptablesfamilyintterms2[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm2_SLEfamilyinteraction_ghqinttophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
} 

for( i in 1:length(toptablesfriendsintterms1) ){
  print(paste0("SLE*friends USM 1 n experiment-wide sig (sens): ",
               nrow(subset(toptablesfriendsintterms1[[i]], P.Value < 1.345479e-7 ))))
  write.csv( subset(toptablesfriendsintterms1[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm1_SLEfriendsinteraction_ghqinttophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
} 

for( i in 1:length(toptablesfriendsintterms2) ){
  print(paste0("SLE*friends USM 2 n experiment-wide sig (sens): ",
               nrow(subset(toptablesfriendsintterms2[[i]], P.Value < 1.338326e-7 ))))
  write.csv( subset(toptablesfriendsintterms2[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm2_SLEfriendsinteraction_ghqinttophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
} 

for( i in 1:length(toptablespartnerintterms1) ){
  print(paste0("SLE*partner USM 1 n experiment-wide sig (sens): ",
               nrow(subset(toptablespartnerintterms1[[i]], P.Value < 1.345479e-7 ))))
  write.csv( subset(toptablespartnerintterms1[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm1_SLEpartnerinteraction_ghqinttophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
} 

for( i in 1:length(toptablespartnerintterms2) ){
  print(paste0("SLE*partner USM 2 n experiment-wide sig (sens): ",
               nrow(subset(toptablespartnerintterms2[[i]], P.Value < 1.338326e-7 ))))
  write.csv( subset(toptablespartnerintterms2[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm2_SLEpartnerinteraction_ghqinttophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
} 


for( i in 1:length(toptablesjbsatintterms1) ){
  print(paste0("SLE*jbsat USM 1 n experiment-wide sig (sens): ",
               nrow(subset(toptablesjbsatintterms1[[i]], P.Value < 1.345479e-7 ))))
  write.csv( subset(toptablesjbsatintterms1[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm1_SLEjobsatinteraction_ghqinttophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
}

for( i in 1:length(toptablesjbsatintterms2) ){
  print(paste0("SLE*jbsat USM 2 n experiment-wide sig (sens): ",
               nrow(subset(toptablesjbsatintterms2[[i]], P.Value < 1.338326e-7 ))))
  write.csv( subset(toptablesjbsatintterms2[[i]], P.Value < 0.001), paste0(
    output_dir,
    "/usm2_SLEjobsatinteraction_ghqinttophits0_001_model_", i, "of4.csv" ),
    row.names = FALSE, quote = FALSE
  )
}
                                  
# Enrichment analyses for model 4 #####

# Top 50 for SLE, top 50 for SLE*GHQ-12, top 50 for each moderator, 
# hence 300 probes (unless overlap) per dataset

tmp <- c(
  head( toptables1[[4]], 50 )$ID,
  head( toptablesghqintterms1[[4]], 50 )$ID,
  head( toptablesfamilyintterms1[[4]], 50 )$ID,
  head( toptablesfriendsintterms1[[4]], 50 )$ID,
  head( toptablespartnerintterms1[[4]], 50 )$ID,
  head( toptablesjbsatintterms1[[4]], 50 )$ID
)
en <- gometh( sig.cpg=tmp,
              all.cpg=toptables1[[4]]$ID,
              collection="KEGG", array.type="EPIC", plot.bias=FALSE,
              prior.prob=TRUE, equiv.cpg=TRUE, fract.counts=TRUE,
              genomic.features="ALL", sig.genes=TRUE )
en <- en[order(en$P.DE),]
write.csv( en, paste0( output_dir, "/USM1_mod4_SLEandSLEGHQ_300_KEGG.csv" ),
           quote=T, row.names = F)
en <- gometh( sig.cpg=tmp,
              all.cpg=toptables1[[4]]$ID,
              collection="GO", array.type="EPIC", plot.bias=FALSE,
              prior.prob=TRUE, equiv.cpg=TRUE, fract.counts=TRUE,
              genomic.features="ALL", sig.genes=TRUE )
en <- en[order(en$P.DE),]
write.csv( en, paste0( output_dir, "/USM1_mod4_SLEandSLEGHQ_300_GO.csv" ),
           quote=T, row.names = F)

tmp <- c(
  head( toptables2[[4]], 50 )$ID,
  head( toptablesghqintterms2[[4]], 50 )$ID,
  head( toptablesfamilyintterms2[[4]], 50 )$ID,
  head( toptablesfriendsintterms2[[4]], 50 )$ID,
  head( toptablespartnerintterms2[[4]], 50 )$ID,
  head( toptablesjbsatintterms2[[4]], 50 )$ID
)
en <- gometh( sig.cpg=tmp,
              all.cpg=toptables2[[4]]$ID,
              collection="KEGG", array.type="EPIC", plot.bias=FALSE,
              prior.prob=TRUE, equiv.cpg=TRUE, fract.counts=TRUE,
              genomic.features="ALL", sig.genes=TRUE )
en <- en[order(en$P.DE),]
write.csv( en, paste0( output_dir, "/USM2_mod4_SLEandSLEGHQ_300_KEGG.csv" ),
           quote=T, row.names = F)
en <- gometh( sig.cpg=tmp,
              all.cpg=toptables2[[4]]$ID,
              collection="GO", array.type="EPIC", plot.bias=FALSE,
              prior.prob=TRUE, equiv.cpg=TRUE, fract.counts=TRUE,
              genomic.features="ALL", sig.genes=TRUE )
en <- en[order(en$P.DE),]
write.csv( en, paste0( output_dir, "/USM2_mod4_SLEandSLEGHQ_300_GO.csv" ),
           quote=T, row.names = F)
