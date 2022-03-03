#install.packages('tidyverse')
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)

setwd("~/data_project/assemblies/EC")

# read in samples to be analyzed
samples <- read.table('samples_EC.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_EC.csv', header = TRUE, sep = ';')
mic <- read.csv2('mic_EC.csv', header = TRUE, sep = ';')
pheno_mic <- read.csv2('pheno_mic_EC.csv', header = TRUE, sep = ';')

phenotypes[phenotypes == "esccol2" | phenotypes == "esccol3"] <- "esccol"
mic[mic == "esccol2" | mic == "esccol3"] <- "esccol"
pheno_mic[pheno_mic == "esccol2" | pheno_mic == "esccol3"] <- "esccol"
phenotypes[phenotypes == "esccolco"] <- "esccol"
mic[mic == "esccolco"] <- "esccol"
pheno_mic[pheno_mic == "esccolco"] <- "esccol"

# check sample ids; and correct them 
which(phenotypes$Sample_id %in% samples$V1 == FALSE)
phenotypes$Sample_id[7] <- "700109-6-13"
phenotypes$Sample_id[48] <- "708018-16-wh"
phenotypes$Sample_id[69] <- "804158-1-17"
phenotypes$Sample_id[73] <- "714843-17-wh"
phenotypes$Sample_id[89] <- "712814-18-wh"
phenotypes$Sample_id[91] <- "713302-18-wh"
phenotypes$Sample_id[96] <- "717578-1-18"
phenotypes$Sample_id[99] <- "719819-2-18"
phenotypes$Sample_id[147] <- "721913-18-wh"
phenotypes$Sample_id[148] <- "721915-18-wh"
phenotypes$Sample_id[168] <- "701215-7-19"
phenotypes$Sample_id[217] <- "806847-11-19"
phenotypes$Sample_id[232] <- "126031-20-wh"
phenotypes$Sample_id[248] <- "808154-21-wh"

#which(mic$Sample_id %in% samples$V1 == FALSE)
mic$Sample_id[7] <- "700109-6-13"
mic$Sample_id[48] <- "708018-16-wh"
mic$Sample_id[69] <- "804158-1-17"
mic$Sample_id[73] <- "714843-17-wh"
mic$Sample_id[89] <- "712814-18-wh"
mic$Sample_id[91] <- "713302-18-wh"
mic$Sample_id[96] <- "717578-1-18"
mic$Sample_id[99] <- "719819-2-18"
mic$Sample_id[147] <- "721913-18-wh"
mic$Sample_id[148] <- "721915-18-wh"
mic$Sample_id[168] <- "701215-7-19"
mic$Sample_id[217] <- "806847-11-19"
mic$Sample_id[232] <- "126031-20-wh"
mic$Sample_id[248] <- "808154-21-wh"

which(pheno_mic$Sample_id %in% samples$V1 == FALSE)
pheno_mic$Sample_id[40] <- "126031-20-wh"
pheno_mic$Sample_id[61] <- "700109-6-13"
pheno_mic$Sample_id[67] <- "701215-7-19"
pheno_mic$Sample_id[102] <- "708018-16-wh"
pheno_mic$Sample_id[108] <- "712814-18-wh"
pheno_mic$Sample_id[109] <- "713302-18-wh"
pheno_mic$Sample_id[115] <- "714843-17-wh"
pheno_mic$Sample_id[119] <- "717578-1-18"
pheno_mic$Sample_id[122] <- "719819-2-18"
pheno_mic$Sample_id[170] <- "721913-18-wh"
pheno_mic$Sample_id[171] <- "721915-18-wh"
pheno_mic$Sample_id[221] <- "804158-1-17"
pheno_mic$Sample_id[241] <- "806847-11-19"
pheno_mic$Sample_id[245] <- "808154-21-wh"

which(phenotypes$Sample_id == "700109-6-13")
phenotypes$LGM_BK[7] <- "esccole"
mic$MIC.LGM_BK[7] <- "esccole"

which(pheno_mic$Sample_id == "700109-6-13")
pheno_mic$LGM_BK[61] <- "esccole"

pheno_mic$Meropenem[which(pheno_mic$Sample_id == "701215-7-19")] <- "S"
pheno_mic$MIC.Meropenem[which(pheno_mic$Sample_id == "701215-7-19")] <- "<= 0.25"

pheno_mic$Meropenem[which(pheno_mic$Sample_id == "806291-19")] <- "S"
pheno_mic$MIC.Meropenem[which(pheno_mic$Sample_id == "806291-19")] <- "<= 0.25"


# combine antibiotics with several columns
pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "S"] <- "S"
pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "R"] <- "R"
pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "I"] <- "I"

pheno_mic$Ceftriaxon[pheno_mic$Ceftriaxon.ohne.Meningitis  == "S"] <- "S"
pheno_mic$Ceftriaxon[pheno_mic$Ceftriaxon.ohne.Meningitis  == "R"] <- "R"
pheno_mic$Ceftriaxon[pheno_mic$Ceftriaxon.ohne.Meningitis  == "I"] <- "I"

for (i in 1:length(pheno_mic$Sample_id)){
  if(pheno_mic$MIC.Meropenem[i] == "NULL") {
    pheno_mic$MIC.Meropenem[i] <- pheno_mic$MIC.Meropenem.ohne.Meningitis[i]
  }
}

for (i in 1:length(pheno_mic$Sample_id)){
  if(pheno_mic$MIC.Ceftriaxon.ohne.Meningitis[i] != "NULL") {
    pheno_mic$MIC.Ceftriaxon[i] <- pheno_mic$MIC.Ceftriaxon.ohne.Meningitis[i]
  }
}

#unique(as.vector(as.matrix(pheno_mic)))
pheno_mic[pheno_mic == "c(R, R)" | pheno_mic == "c(S, R)" | pheno_mic == "c(R, S)" | pheno_mic == "c(S, U, R)" |   
                      pheno_mic == "c(R, S, S)" | pheno_mic == "c(R, R, R)" | pheno_mic == "c(R, U)" | 
                      pheno_mic == "c(U, R)" | pheno_mic == "c(U, R, U, R)" | pheno_mic == "c(R, R, I)" |
                      pheno_mic == "c(R, U, U, I, I)" | pheno_mic == "c(I, U, U, R)" | pheno_mic == "c(R, U, U, I)"] <- "R"

pheno_mic[pheno_mic == "c(S, S)" | pheno_mic == "c(U, S)" | pheno_mic == "c(S, S, S)" | pheno_mic == "c(S, U)" |
                      pheno_mic == "c(S, U, S)" | pheno_mic == "c(S, S, U)" | pheno_mic == "c(S, U, U, S)" |
                      pheno_mic == "c(U, S, S)"] <- "S"

pheno_mic[pheno_mic == "c(I, S)" | pheno_mic == "c(S, I)" | pheno_mic == "c(I, S, S)" |
                      pheno_mic == "c(S, I, S, I)" | pheno_mic == "c(I, I)" | pheno_mic == "c(I, U)" | 
                      pheno_mic == "c(U, I)" | pheno_mic == "c(U, I, U, U)" | pheno_mic == "c(U, I, U, S)" |
                      pheno_mic == "c(U, S, U, I)" | pheno_mic == "c(U, I, U)" | pheno_mic == "c(I, I, U)" |
                      pheno_mic == "c(I, U, U)" | pheno_mic == "c(U, S, I)"] <- "I"

pheno_mic[pheno_mic == "U" | pheno_mic == "c(U, U)" | pheno_mic == "c(U, U, U)"] <- 0


pheno_mic[pheno_mic == "c(n.a. 0, = 0)" | pheno_mic == "n.a. 0" | pheno_mic == "c(n.a. 0, n.a. 0)" | pheno_mic == "= 0" ] <- "NA" 
pheno_mic[pheno_mic == "c(= 0.016, n.a. 0)"] <- "= 0.016"
pheno_mic[pheno_mic == "c(<= 0.25, <= 0.25)" | pheno_mic == "c(= 0.047, <= 0.25)" | pheno_mic == "c(<= 0.25, <= 0.25, <= 0.25)"] <-  "<= 0.25"
pheno_mic[pheno_mic == "c(<= 0.5, <= 0.5)" | pheno_mic == "c(= 0.38, = 0.5)" | pheno_mic == "c(= 0, = 0.5)" | 
            pheno_mic == "c(<= 0.5, <= 0.5, <= 0.5)" | pheno_mic == "c(= 0.75, <= 0.5)" | pheno_mic == "c(= 0.5, = 0.19)" | pheno_mic == "= 0.5"] <-  "<= 0.5"
pheno_mic[pheno_mic == "c(= 0, = 0.75)"] <- "= 0.75"
pheno_mic[pheno_mic == "c(n.a. 0, <= 1)" | pheno_mic == "c(<= 1, < 0.5)" | pheno_mic == "c(<= 1, = 1)" | 
           pheno_mic == "c(= 0, <= 1)" | pheno_mic == "c(<= 1, <= 1)" | pheno_mic == "c(= 1.5, <= 1)" | pheno_mic == "c(<= 1, = 6)" |
           pheno_mic == "c(<= 1, = 0, = 0, <= 1)" | pheno_mic == "c(<= 1, = 0)" | pheno_mic == "c(<= 1, <= 1, = 0)" |
          pheno_mic == "c(n.a. 0, = 1, = 1)" | pheno_mic == "c(= 0.5, = 1)" | pheno_mic == "c(= 0.5, <= 1)" | pheno_mic == "c(>= 8, = 1)" |
          pheno_mic == "c(= 1, = 0.25)" | pheno_mic == "= 1"] <- "<= 1" 
pheno_mic[pheno_mic == "c(= 0, = 1.5)" ] <- "= 1.5" 
pheno_mic[pheno_mic == "c(= 2, = 0)" | pheno_mic == "c(<= 1, = 2)" | pheno_mic == "c(n.a. 0, = 2)" | pheno_mic == "c(= 2, = 0, n.a. 0)" | 
          pheno_mic == "c(= 2, = 0.25, = 0.19)" | pheno_mic == "c(= 2, = 0.25)" | pheno_mic == "c(= 2, = 0, n.a. 0, = 12)" |
            pheno_mic == "c(= 0, <= 1, n.a. 0, = 2)" | pheno_mic == "c(= 1, = 2)" | pheno_mic == "c(= 0, = 2, = 0)" |
            pheno_mic == "c(= 2, = 2)" | pheno_mic == "c(n.a. 0, <= 1, = 2)"] <- "= 2" 
pheno_mic[pheno_mic == "c(= 4, = 4)" | pheno_mic == "c(= 4, n.a. 0)" | pheno_mic == "c(n.a. 0, = 4)" | pheno_mic == "c(= 0, = 4)" | 
           pheno_mic == "c(> 32, = 4)" | pheno_mic == "c(= 4, = 6)" | pheno_mic == "c(= 4, = 2)" | pheno_mic == "c(= 1.5, = 4)" | 
           pheno_mic == "c(= 1, = 4)" | pheno_mic == "c(n.a. 0, = 4, n.a. 0)" | pheno_mic == "c(= 0, = 4, = 0, = 1.5)" |
           pheno_mic == "c(= 2, = 4, <= 0.25, = 3)" | pheno_mic == "c(= 4, = 4, = 0)" | pheno_mic == "c(= 4, = 0)"] <- "= 4"
pheno_mic[pheno_mic == "c(= 6, = 6)" | pheno_mic == "c(= 0, = 6, = 4, = 0)"] <-  "= 6"
pheno_mic[pheno_mic == "c(n.a. 0, = 8)" | pheno_mic == "c(>= 8, = 4)" | pheno_mic == "c(= 8, = 0)" | pheno_mic == "c(= 0.75, = 8)"] <- "= 8" 
pheno_mic[pheno_mic == "c(= 16, = 24)" | pheno_mic == "c(= 1.5, >= 16)" | pheno_mic == "c(= 6, >= 16)" | pheno_mic == "c(= 16, = 4)" |
           pheno_mic == "c(n.a. 0, = 16)" | pheno_mic == "= 16" | pheno_mic == "c(n.a. 0, = 16, n.a. 0, = 16)" | 
           pheno_mic == "c(>= 64, = 16)" | pheno_mic == "c(= 16, = 16)" | pheno_mic == "c(>= 64, = 16, = 4)" | pheno_mic == "c(= 32, = 16)" |
           pheno_mic == "c(> 32, = 16)" | pheno_mic == "c(= 8, = 16)" | pheno_mic == "c(>= 32, = 16)" | pheno_mic == "c(= 16, > 32)" |
           pheno_mic == "c(= 16, = 12)" | pheno_mic == "c(= 16, = 24)" | pheno_mic == "c(>= 16, = 24)"] <- ">= 16" 
pheno_mic[pheno_mic == "= 32" | pheno_mic == "c(= 32, n.a. 0)" | pheno_mic == "c(n.a. 0, = 32)" | pheno_mic == "c(>= 8, > 32)"] <- ">= 32" 
pheno_mic[pheno_mic == "c(> 32, >= 64)" | pheno_mic == "c(>= 64, > 32)" | pheno_mic == "c(>= 64, >= 64)" 
         | pheno_mic == "c(>= 64, >= 64, >= 64)" | pheno_mic == "c(>= 64, = 32)" | pheno_mic == "c(>= 64, <= 1)" | 
           pheno_mic == "c(>= 64, n.a. 0)" | pheno_mic == "c(n.a. 0, >= 64)" | pheno_mic == "c(= 0, >= 64)" |
           pheno_mic == "c(>= 64, = 0, n.a. 0, = 2)" | pheno_mic == "c(>= 64, = 0, n.a. 0, = 2, = 2)" | pheno_mic == "c(>= 64, = 0)"] <- ">= 64" 

names(pheno_mic)[names(pheno_mic) == "Ceftriaxon"] <- "Ceftriaxone"
names(pheno_mic)[names(pheno_mic) == "MIC.Ceftriaxon"] <- "MIC.Ceftriaxone"
names(pheno_mic)[names(pheno_mic) == "Cefepim"] <- "Cefepime"
names(pheno_mic)[names(pheno_mic) == "MIC.Cefepim"] <- "MIC.Cefepime"
names(pheno_mic)[names(pheno_mic) == "Ceftazidim"] <- "Ceftazidime"
names(pheno_mic)[names(pheno_mic) == "MIC.Ceftazidim"] <- "MIC.Ceftazidime"

#-------------------------------------------------------------------------------
# get result tables
# abricate
abr_ncbi <- read.delim('EC_abricate_ncbi.tab', sep = "\t", header = T)
if (length(which(abr_ncbi$X.IDENTITY < 95)) != 0) {
  abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$X.IDENTITY < 95)),]
}
abr_ncbi$X.FILE <- gsub(".fna","", abr_ncbi$X.FILE)
names(abr_ncbi)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_ncbi$Tool <- 'abr_ncbi'

abr_card <- read.delim('EC_abricate_card.tab', sep = "\t", header = T)
if (length(which(abr_card$X.IDENTITY < 95)) != 0) {
  abr_card <- abr_card[-c(which(abr_card$X.IDENTITY < 95)),]
}
abr_card$X.FILE <- gsub(".fna","", abr_card$X.FILE)
names(abr_card)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_card$Tool <- 'abr_card'

abr_megares <- read.delim('EC_abricate_megares.tab', sep = "\t", header = T)
if (length(which(abr_megares$X.IDENTITY < 95)) != 0) {
  abr_megares <- abr_megares[-c(which(abr_megares$X.IDENTITY < 95)),]
}
abr_megares$X.FILE <- gsub(".fna","", abr_megares$X.FILE)
names(abr_megares)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_megares$Tool <- 'abr_megares'

abr_resfinder <- read.delim('EC_abricate_resfinder.tab', sep = "\t", header = T)
if (length(which(abr_resfinder$X.IDENTITY < 95)) != 0) {
  abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$X.IDENTITY < 95)),]
}
abr_resfinder$X.FILE <- gsub(".fna","", abr_resfinder$X.FILE)
names(abr_resfinder)[c(1,14)] <- c('Sample_id', 'amr_genes')
abr_resfinder$Tool <- 'abr_resfinder'

abr_argannot <- read.delim('EC_abricate_argannot.tab', sep = "\t", header = T)
if (length(which(abr_argannot$X.IDENTITY < 95)) != 0) {
  abr_argannot <- abr_argannot[-c(which(abr_argannot$X.IDENTITY < 95)),]
}
abr_argannot$X.FILE <- gsub(".fna","", abr_argannot$X.FILE)
names(abr_argannot)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_argannot$Tool <- 'abr_argannot'


#amrfinder
amrf_nuc <- read.delim('EC_amrfinder_nuc.tab', sep = "\t")
if (length(which(amrf_nuc$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_nuc)[1] <- 'Sample_id'
amrf_nuc$Sample_id <- gsub("gnl\\|USB\\|","", amrf_nuc$Contig.id)
amrf_nuc$Sample_id <- gsub('(_\\d*)',"", amrf_nuc$Sample_id)
names(amrf_nuc)[6] <- 'amr_genes'
amrf_nuc$Tool <- 'amrf_nuc'

amrf_prot <- read.delim('EC_amrfinder_prot.tab', sep = "\t")
if (length(which(amrf_prot$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_prot <- amrf_prot[-c(which(amrf_prot$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_prot)[c(1,2)] <- c('Sample_id', 'amr_genes')
amrf_prot$Sample_id <- gsub('(_\\d*)',"", amrf_prot$Sample_id)
amrf_prot$Tool <- 'amrf_prot'


# rgi
rgi <- read.delim(('EC_rgi.tab'), sep = "\t")
if (length(which(rgi$Best_Identities < 95)) != 0) {
  rgi <- rgi[-c(which(rgi$Best_Identities < 95)),]
}
names(rgi)[c(1,2)] <- c('Sample_id', 'amr_genes')
rgi$Sample_id <- gsub("gnl\\|USB\\|","", rgi$Sample_id)
rgi$Sample_id <- gsub('(_\\d*)',"", rgi$Sample_id)
rgi$Sample_id <- gsub('\\s+', '', rgi$Sample_id)
rgi$Tool <- 'rgi'

# sraX
srax_basic <- read.delim(('EC_srax_basic.tab'), sep = "\t")
if (length(which(srax_basic$Identity_p < 95)) != 0) {
  srax_basic <- srax_basic[-c(which(srax_basic$Identity_p < 95)),]
}
names(srax_basic)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_basic$Sample_id <- gsub(".fna","", srax_basic$Sample_id)
srax_basic$Tool <- 'sraX_basic'

srax_ext <- read.delim(('EC_srax_ext.tab'), sep = "\t")
if (length(which(srax_ext$Identity_p < 95)) != 0) {
  srax_ext <- srax_ext[-c(which(srax_ext$Identity_p < 95)),]
}
names(srax_ext)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_ext$Sample_id <- gsub(".fna","", srax_ext$Sample_id)
srax_ext$Tool <- 'sraX_ext'


# Deeparg
deeparg_LS <- read.delim(('EC_deeparg_LS_noPot.tab'), sep = "\t")
if (length(which(deeparg_LS$identity < 95)) != 0) {
  deeparg_LS <- deeparg_LS[-c(which(deeparg_LS$identity < 95)),]
}
if (length(which(deeparg_LS$probability < 0.8)) != 0) {
  deeparg_LS <- deeparg_LS[-c(which(deeparg_LS$probability < 0.8)),]
}
names(deeparg_LS)[4] <- 'Sample_id'
deeparg_LS$Sample_id <- gsub("gnl\\|USB\\|","", deeparg_LS$Sample_id)
deeparg_LS$Sample_id <- gsub('(_\\d*)',"", deeparg_LS$Sample_id)
for(i in 1:length(deeparg_LS[,1])) {
  genes <- unlist(strsplit(deeparg_LS[i,6], "\\|"))
  deeparg_LS[i,6] <- genes[3]
}
names(deeparg_LS)[6] <- 'amr_genes'
deeparg_LS$Tool <- 'deeparg_LS'

deeparg_SR <- read.delim(('EC_deeparg_SR_noPot.tab'), sep = "\t")
if (length(which(deeparg_SR$identity < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$identity < 95)),]
}
if (length(which(deeparg_SR$probability < 0.8)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$probability < 0.8)),]
}
deeparg_SR_unique <- deeparg_SR %>% distinct(X.ARG, best.hit, .keep_all = TRUE)
names(deeparg_SR_unique)[1] <- 'Sample_id'
deeparg_SR_unique$Sample_id <- gsub("esccol/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
deeparg_SR_unique$Sample_id <- gsub("esccole/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
deeparg_SR_unique$Sample_id <- gsub("esccolc/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
deeparg_SR_unique$Sample_id <- gsub("esccolco/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
for(i in 1:length(deeparg_SR_unique[,1])) {
  genes <- unlist(strsplit(deeparg_SR_unique[i,1], "\\_"))
  deeparg_SR_unique[i,1] <- genes[1]
}
for(i in 1:length(deeparg_SR_unique[,1])) {
  genes <- unlist(strsplit(deeparg_SR_unique[i,6], "\\|"))
  deeparg_SR_unique[i,6] <- genes[3]
}
#deeparg_SR_unique <- deeparg_SR_unique[-c(which(deeparg_SR_unique$best.hit == "undefined")), ]
names(deeparg_SR_unique)[6] <- 'amr_genes'
deeparg_SR_unique$Tool <- 'deeparg_SR'

rm(deeparg_SR)

# resfinder
resfinder_as <- read.delim(('EC_resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$Identity < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$Identity < 95)),]
}
names(resfinder_as)[c(6,1)] <- c('Sample_id', 'amr_genes')
resfinder_as$Sample_id <- gsub('gnl\\|USB\\|','', resfinder_as$Sample_id)
resfinder_as$Sample_id <- gsub('(_\\d*)','', resfinder_as$Sample_id)
resfinder_as$Tool <- 'resfinder_as'

resfinder_re <- read.delim(('EC_resfinder_reads.tab'), sep = "\t", header = F)
names(resfinder_re) <- c('Resistance gene',	'Identity',	'Alignment Length/Gene Length',	'Coverage',	'Position in reference',	'Contig',	'Position in contig',  'Phenotype',	'Accession no.', 'Sample_id')
resfinder_re <- resfinder_re[-c(1),]
resfinder_re$Identity <- as.numeric(resfinder_re$Identity )
if (length(which(resfinder_re$Identity < 95)) != 0) {
  resfinder_re <- resfinder_re[-c(which(resfinder_re$Identity < 95)),]
}
names(resfinder_re)[1] <- 'amr_genes'
resfinder_re$Tool <- 'resfinder_re'


# extract lines with bla genes / ESBL, Carba resistance 
ncbi_bla <-abr_ncbi[c(which(abr_ncbi$RESISTANCE == "CARBAPENEM" | abr_ncbi$RESISTANCE == "BETA-LACTAM" | abr_ncbi$RESISTANCE == "CEPHALOSPORIN")), ]
names(ncbi_bla)[6] <- 'bla_genes'
ncbi_bla$Tool <- 'abr_ncbi'
card_bla <-abr_card[c(grep("beta-lactam", abr_card$PRODUCT)), ]
names(card_bla)[6] <- 'bla_genes'
card_bla$Tool <- 'abr_card'
megares_bla <-abr_megares[c(grep("betalactam", abr_megares$PRODUCT)), ]
names(megares_bla)[6] <- 'bla_genes'
megares_bla$Tool <- 'abr_megares'
resfinder_bla <-abr_resfinder[c(grep("bla", abr_resfinder$amr_genes)), ]
names(resfinder_bla)[14] <- 'bla_genes'
resfinder_bla$Tool <- 'abr_resfinder'
argannot_bla <-abr_argannot[c(grep("bla", abr_argannot$PRODUCT)), ]
names(argannot_bla)[6] <- 'bla_genes'
argannot_bla$Tool <- 'abr_argannot'

amrf_nuc_bla <- amrf_nuc[c(which(amrf_nuc$Subclass == "CARBAPENEM" | amrf_nuc$Subclass == "BETA-LACTAM" | amrf_nuc$Subclass == "CEPHALOSPORIN")), ]
names(amrf_nuc_bla)[6] <- 'bla_genes'
amrf_nuc_bla$Tool <- 'amrf_nuc'
amrf_prot_bla <- amrf_prot[c(which(amrf_prot$Subclass == "CARBAPENEM" | amrf_prot$Subclass == "BETA-LACTAM" | amrf_prot$Subclass == "CEPHALOSPORIN")), ]
names(amrf_prot_bla)[2] <- 'bla_genes'
amrf_prot_bla$Tool <- 'amrf_prot'

deeparg_LS_bla <- deeparg_LS[c(which(deeparg_LS$predicted_ARG.class == "beta-lactam")), ]
names(deeparg_LS_bla)[6] <- 'bla_genes'
deeparg_LS_bla$Tool <- 'deeparg_LS'
deeparg_SR_bla <- deeparg_SR_unique[c(which(deeparg_SR_unique$predicted_ARG.class == "beta-lactam")), ]
names(deeparg_SR_bla)[6] <- 'bla_genes'
deeparg_SR_bla$Tool <- 'deeparg_SR'

resfinder_as_bla <- resfinder_as[c(grep("bla", resfinder_as$amr_genes)), ]
names(resfinder_as_bla)[1] <- 'bla_genes'
resfinder_as_bla$Tool <- 'resfinder_as'
resfinder_re_bla <- resfinder_re[c(grep("bla", resfinder_re$amr_genes)), ]
names(resfinder_re_bla)[1] <- 'bla_genes'
resfinder_re_bla$Tool <- 'resfinder_re'

rgi_bla <- rgi[c(grep("beta-lactamase", rgi$AMR.Gene.Family)), ]
names(rgi_bla)[2] <- 'bla_genes'
rgi_bla$Tool <- 'rgi'

srax_basic_bla <- srax_basic[c(grep("beta-lactamase", srax_basic$Gene_description, ignore.case = TRUE)), ]
names(srax_basic_bla)[6] <- 'bla_genes'
srax_basic_bla$Tool <- 'srax_basic'
srax_ext_bla <- srax_ext[c(grep("beta-lactamase", srax_ext$Gene_description, ignore.case = TRUE)), ]
names(srax_ext_bla)[6] <- 'bla_genes'
srax_ext_bla$Tool <- 'srax_ext'

# collect all genes found in each tool in one table
bla <- rbind(ncbi_bla[,c(1,6,16)], card_bla[,c(1,6,16)], megares_bla[,c(1,6,16)], resfinder_bla[,c(1,14,16)], argannot_bla[,c(1,6,16)],
             amrf_nuc_bla[,c(1,6,23)], amrf_prot_bla[,c(1,2,19)], deeparg_LS_bla[,c(4,6,13)], deeparg_SR_bla[,c(1,6,13)],
             resfinder_as_bla[,c(6,1,10)], resfinder_re_bla[,c(10,1,11)], rgi_bla[,c(1,2,10)], srax_basic_bla[,c(2,6,17)], srax_ext_bla[,c(2,6,17)])

# remove blank spaces in sample_ids
bla$Sample_id <- gsub('\\s+', '', bla$Sample_id)
bla <- unique(bla)


# include clinical phenotype 
for (i in 1:length(bla$Sample_id)){
  bla$Type[i] <- pheno_mic$LGM_BK[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}


for (i in 1:length(bla$Sample_id)){
  bla$MeropenemMIC[i] <- pheno_mic$MIC.Meropenem[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$ImipenemMIC[i] <- pheno_mic$MIC.Imipenem[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$ErtapenemMIC[i] <- pheno_mic$MIC.Ertapenem[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}
  
for (i in 1:length(bla$Sample_id)){
  bla$CeftriaxoneMIC[i] <- pheno_mic$MIC.Ceftriaxone[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}
  
for (i in 1:length(bla$Sample_id)){
  bla$CeftazidiemMIC[i] <- pheno_mic$MIC.Ceftazidime[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}
  
for (i in 1:length(bla$Sample_id)){
  bla$CefepimeMIC[i] <- pheno_mic$MIC.Cefepime[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}    
  

for (i in 1:length(bla$Sample_id)){
  bla$Meropenem[i] <- pheno_mic$Meropenem[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$Imipenem[i] <- pheno_mic$Imipenem[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$Ertapenem[i] <- pheno_mic$Ertapenem[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$Ceftriaxone[i] <- pheno_mic$Ceftriaxone[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$Ceftazidime[i] <- pheno_mic$Ceftazidime[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$Cefepime[i] <- pheno_mic$Cefepime[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}  


#-------------------------------------------------------------------------------
# gene cluster

bla_gene_clusters <- bla
bla_gene_clusters$bla_genes <- gsub("\\(Bla\\)", "", bla_gene_clusters$bla_genes)
bla_gene_clusters$bla_genes <- gsub("bla", "", bla_gene_clusters$bla_genes)

bla_gene_clusters[bla_gene_clusters == "KPC-1" | bla_gene_clusters == "KPC-2"] <- "KPC-1/2"
bla_gene_clusters[bla_gene_clusters == "Escherichia_coli_ampH_beta-lactamase" | bla_gene_clusters == "Escherichia coli ampH beta-lactamase" |
                  bla_gene_clusters == "Escherichia_coli_ampH" | bla_gene_clusters == "AMPH"] <- "ampH"
bla_gene_clusters[bla_gene_clusters == "Escherichia coli ampC beta-lactamase" | bla_gene_clusters == "Escherichia_coli_ampC_beta-lactamase" |
                  bla_gene_clusters == "Escherichia_coli_ampC" | bla_gene_clusters == "AMPC"] <- "ampC"
bla_gene_clusters[bla_gene_clusters == "Escherichia coli ampC1 beta-lactamase" | bla_gene_clusters == "Escherichia_coli_ampC1_beta-lactamase"] <- "ampC1"


gene_clusters <- pivot_wider(bla_gene_clusters, names_from = Tool, values_from = bla_genes)

gene_cluster_table <- pivot_longer(data = gene_clusters, 
                                   cols = -c(1:14),
                                   names_to = "ToolDB", 
                                   values_to = "gene_cluster")



gene_cluster_table <- as.data.frame(gene_cluster_table)
gene_cluster_table[gene_cluster_table == "NULL"] <- "none"

for (i in 1:length(gene_cluster_table$gene_cluster)) {
  gene_cluster_table$gene_cluster[[i]] <- sort(gene_cluster_table$gene_cluster[[i]])
}

gene_cluster_table$Tool <- ifelse(gene_cluster_table$ToolDB == "rgi", gene_cluster_table$Tool <- "RGI",
                                  ifelse(gene_cluster_table$ToolDB == "srax_basic" | gene_cluster_table$ToolDB == "srax_ext", gene_cluster_table$Tool <- "sraX",
                                         ifelse(gene_cluster_table$ToolDB == "resfinder_as" | gene_cluster_table$ToolDB == "resfinder_re", gene_cluster_table$Tool <- "ResFinder", 
                                                ifelse(gene_cluster_table$ToolDB == "deeparg_LS" | gene_cluster_table$ToolDB == "deeparg_SR", gene_cluster_table$Tool <- "DeepARG",
                                                       ifelse(gene_cluster_table$ToolDB == "amrf_nuc" | gene_cluster_table$ToolDB == "amrf_prot", gene_cluster_table$Tool <- "AMRFinder",
                                                              gene_cluster_table$Tool <- "ABRicate")))))

gene_cluster_table[gene_cluster_table == "abr_ncbi"] <- "ncbi"
gene_cluster_table[gene_cluster_table == "abr_card"] <- "card"
gene_cluster_table[gene_cluster_table == "abr_resfinder"] <- "resfinder"
gene_cluster_table[gene_cluster_table == "abr_megares"] <- "megares"
gene_cluster_table[gene_cluster_table == "abr_argannot"] <- "argannot"
gene_cluster_table[gene_cluster_table == "amrf_nuc"] <- "nuc"
gene_cluster_table[gene_cluster_table == "amrf_prot"] <- "prot"
gene_cluster_table[gene_cluster_table == "srax_basic"] <- "basic"
gene_cluster_table[gene_cluster_table == "srax_ext"] <- "ext"
gene_cluster_table[gene_cluster_table == "deeparg_LS"] <- "LS"
gene_cluster_table[gene_cluster_table == "deeparg_SR"] <- "SR"
gene_cluster_table[gene_cluster_table == "resfinder_as"] <- "assembly"
gene_cluster_table[gene_cluster_table == "resfinder_re"] <- "reads"


### clusters esccolc
cluster_esccolc <- gene_cluster_table[c(which(gene_cluster_table$Type == "esccolc")),]


cluster_esccolc$cluster_rep <- ifelse(grepl("KPC-1/2",cluster_esccolc$gene_cluster) == "TRUE", 1, 
                                      ifelse(grepl("NDM-1",cluster_esccolc$gene_cluster) == "TRUE", 2,
                                             ifelse(grepl("NDM-5",cluster_esccolc$gene_cluster) == "TRUE", 3,
                                                    ifelse(grepl("NDM",cluster_esccolc$gene_cluster) == "TRUE", 4, 
                                                           ifelse(grepl("OXA-48",cluster_esccolc$gene_cluster) == "TRUE", 5,
                                                                  ifelse(grepl("OXA-181",cluster_esccolc$gene_cluster) == "TRUE", 6,
                                                                         ifelse(grepl("OXA-244",cluster_esccolc$gene_cluster) == "TRUE", 7, 
                                                                                ifelse(grepl("KPC",cluster_esccolc$gene_cluster) == "TRUE", 8,
                                                                                       ifelse(grepl("\\bOXA\\b",cluster_esccolc$gene_cluster) == "TRUE", 11, 10)))))))))

for (i in 1:length(cluster_esccolc$Sample_id)) {
  if (cluster_esccolc$ToolDB[i] == "SR") {
    cluster_esccolc$cluster_rep[i] <- 9
  }
}

cluster_esccolc$cluster_rep <- as.character(cluster_esccolc$cluster_rep)

a <- sort(unique(cluster_esccolc$Sample_id[cluster_esccolc$ToolDB == "argannot" & cluster_esccolc$cluster_rep == "1"]))
b <- sort(unique(cluster_esccolc$Sample_id[cluster_esccolc$ToolDB == "argannot" & cluster_esccolc$cluster_rep == "2"]))
c <- sort(unique(cluster_esccolc$Sample_id[cluster_esccolc$ToolDB == "argannot" & (cluster_esccolc$cluster_rep == "3")]))
d <- sort(unique(cluster_esccolc$Sample_id[cluster_esccolc$ToolDB == "argannot" & (cluster_esccolc$cluster_rep == "5")]))
e <- sort(unique(cluster_esccolc$Sample_id[cluster_esccolc$ToolDB == "argannot" & (cluster_esccolc$cluster_rep == "6")]))
f <- sort(unique(cluster_esccolc$Sample_id[cluster_esccolc$ToolDB == "argannot" & (cluster_esccolc$cluster_rep == "7")]))
g <- sort(unique(cluster_esccolc$Sample_id[cluster_esccolc$ToolDB == "argannot" & (cluster_esccolc$cluster_rep == "10")]))

cluster_esccolc$Sample_id <- factor(cluster_esccolc$Sample_id, levels = c(a, b, c, d, e, f, g)) 

# order legend
cluster_esccolc$cluster_rep <- factor(cluster_esccolc$cluster_rep, levels = c(1,8,2,3,4,5,6,7,11,9,10)) 
colors <- c( "seagreen3", "seagreen1", "blue4", "blue", "lightsteelblue", "orangered4", "orangered", "orange", "mediumpurple1", "yellowgreen", "moccasin")

cluster_esccolc_heatmap <- ggplot(data = cluster_esccolc, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("KPC-1/2", "KPC", "NDM-1", "NDM-5",
                                                "NDM", "OXA-48", "OXA-181", "OXA-244", 
                                                "undefined OXA", "mixed", "none")) +
  labs(x = "", y = "") +
  labs(fill = "Carbapenemase genes")
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(cluster_esccolc$MeropenemMIC, levels = c("<= 0.25", "<= 0.5", "= 0.75","<= 1", "= 2",
                                                                     "= 4", "= 8", ">= 16"))), vars(Tool), scales = "free", space = "free")

cluster_esccolc_heatmap


### clusters esccole
cluster_esccole <- gene_cluster_table[c(which(gene_cluster_table$Type == "esccole")),]

ctxm1 <- c("\\bCTX-M-1\\b|\\bCTX-M-3\\b|\\bCTX-M-15\\b|\\bCTX-M-55\\b")
ctxm8 <- c("\\bCTX-M-8\\b")
ctxm9 <- c("\\bCTX-M-14\\b|\\bCTX-M-24\\b|\\bCTX-M-27\\b|\\bCTX-M-122\\b|\\bCTX-M-125\\b|\\bCTX-M-129\\b|\\bCTX-M-174\\b")
tem <- c("\\bTEM-7\\b|\\bTEM-15\\b|\\bTEM-17\\b|\\bTEM-106\\b|\\bTEM-52B\\b")
shv <- c("\\bSHV-5\\b|\\bSHV-12\\b|\\bSHV-66\\b|\\bSHV-134\\b")
cmy <- "\\bCMY-42\\b"
ampc <- c("\\bCMY-2\\b|\\bCMY-4\\b|\\bCMY-111\\b|\\bCMY-59\\b|\\bCMY-163\\b")

cluster_esccole$cluster_rep <- ifelse(grepl(cmy,cluster_esccole$gene_cluster) == "TRUE" & grepl(ctxm1,cluster_esccole$gene_cluster) == "TRUE" & grepl(ctxm9,cluster_esccole$gene_cluster) == "TRUE", 1,
                                      ifelse(grepl(ctxm1,cluster_esccole$gene_cluster) == "TRUE" & grepl(ctxm9,cluster_esccole$gene_cluster) == "TRUE" & grepl(ampc,cluster_esccole$gene_cluster) == "TRUE", 2,
                                             ifelse(grepl(ctxm9,cluster_esccole$gene_cluster) == "TRUE" & grepl(tem, cluster_esccole$gene_cluster) == "TRUE", 3,
                                                    ifelse(grepl(ctxm1,cluster_esccole$gene_cluster) == "TRUE", 4,
                                                           ifelse(grepl(ctxm8,cluster_esccole$gene_cluster) == "TRUE", 5,
                                                                  ifelse(grepl(ctxm9,cluster_esccole$gene_cluster) == "TRUE", 6,
                                                                         ifelse(grepl(cmy,cluster_esccole$gene_cluster) == "TRUE", 7,
                                                                                ifelse(grepl(shv,cluster_esccole$gene_cluster) == "TRUE", 8,
                                                                                       ifelse(grepl(tem,cluster_esccole$gene_cluster) == "TRUE", 9,
                                                                                              ifelse(grepl(ampc,cluster_esccole$gene_cluster) == "TRUE", 12, 11))))))))))


cluster_esccole$cluster_rep <- as.character(cluster_esccole$cluster_rep)

for (i in 1:length(cluster_esccole$Sample_id)) {
  if (cluster_esccole$ToolDB[i] == "SR") {
    cluster_esccole$cluster_rep[i] <- 10
  }
}

for (i in 1:length(cluster_esccole$Sample_id)) {
  if (cluster_esccole$ToolDB[i] == "megares" & grepl("\\bOXA\\b|\\bSHV\\b|\\bTEM\\b",cluster_esccolc$gene_cluster) == "TRUE") {
    cluster_esccole$cluster_rep[i] <- 13
  }
}


# order legend
cluster_esccole$cluster_rep <- factor(cluster_esccole$cluster_rep, levels = c(4,5,6,1,2,3,8,9,7,12,13,10,11)) 

h <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & cluster_esccole$cluster_rep == "4"]))
z <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & cluster_esccole$cluster_rep == "5"]))
j <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & (cluster_esccole$cluster_rep == "6")]))
k <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & (cluster_esccole$cluster_rep == "1")]))
l <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & (cluster_esccole$cluster_rep == "2")]))
m <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & (cluster_esccole$cluster_rep == "7")]))
n <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & (cluster_esccole$cluster_rep == "8")]))
o <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & (cluster_esccole$cluster_rep == "9")]))
p <- sort(unique(cluster_esccole$Sample_id[cluster_esccole$ToolDB == "argannot" & (cluster_esccole$cluster_rep == "11")]))

cluster_esccole$Sample_id <- factor(cluster_esccole$Sample_id, levels = c(h,z,j,k,l,m,n,o,p))

colors <- c( "darkblue", "orangered", "seagreen4", "salmon", "goldenrod1", "royalblue", "turquoise2", "maroon", "magenta", "black", "mediumpurple1", "yellowgreen","moccasin")

cluster_esccole_heatmap <- ggplot(data = cluster_esccole, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("CTX-M-1 Group", "CTX-M-8 Group", "CTX-M-9 Group", "CTX-M-1 Group & CTX-M-9 Group & CMY ESBL type",  
                                                "CTX-M-1 Group & CTX-M-9 Group & CMY AmpC type", "CTX-M-9 Group & TEM ESBL type", 
                                                "SHV ESBL type", "TEM ESBL type", "CMY ESBL type", "CMY AmpC type","undefined CTX/SHV/TEM","mixed", "none")) +
  labs(x = "", y = "") +
  labs(fill = "ESBL genes") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(cluster_esccole$CeftriaxoneMIC, levels = c("<= 1", "= 4", "= 8", ">= 16", ">= 32", ">= 64", "none","NA"))), vars(Tool), scales = "free", space = "free")

cluster_esccole_heatmap


### clusters esccol
cluster_esccol <- gene_cluster_table[c(which(gene_cluster_table$Type == "esccol")),]

mbl <- c("metallo-beta-lactamase")
tem <- c("\\bTEM-7\\b|\\bTEM-12\\b|\\bTEM-19\\b|\\bTEM-52B\\b")
shv <- c("\\bSHV-102\\b")
cmy <- c("\\bCMY-8\\b|\\bCMY-9\\b|\\bCMY-30\\b|\\bCMY-37\\b|\\bCMY-42\\b")
ampc <- c("\\bCMY-2\\b|\\bCMY-4\\b|\\bCMY-111\\b|\\bCMY-59\\b")

cluster_esccol$cluster_rep <- ifelse(grepl(mbl,cluster_esccol$gene_cluster) == "TRUE", 5,
                                ifelse(grepl(ampc,cluster_esccol$gene_cluster) == "TRUE",7,
                                  ifelse(grepl(tem,cluster_esccol$gene_cluster) == "TRUE", 1,
                                     ifelse(grepl(shv,cluster_esccol$gene_cluster) == "TRUE", 2, 
                                            ifelse(grepl(cmy,cluster_esccol$gene_cluster) == "TRUE", 3, 
                                                   ifelse(grepl("IMP-66",cluster_esccol$gene_cluster) == "TRUE",4,6))))))

cluster_esccol$cluster_rep <- as.character(cluster_esccol$cluster_rep)
cluster_esccol$cluster_rep <- factor(cluster_esccol$cluster_rep, levels = c(1,2,3,7,4,5,6)) 

v <- sort(unique(cluster_esccol$Sample_id[cluster_esccol$ToolDB == "SR" & cluster_esccol$cluster_rep == "7"]))
q <- sort(unique(cluster_esccol$Sample_id[cluster_esccol$ToolDB == "SR" & cluster_esccol$cluster_rep == "1"]))
r <- sort(unique(cluster_esccol$Sample_id[cluster_esccol$ToolDB == "SR" & cluster_esccol$cluster_rep == "2"]))
s <- sort(unique(cluster_esccol$Sample_id[cluster_esccol$ToolDB == "SR" & cluster_esccol$cluster_rep == "3"]))
t <- sort(unique(cluster_esccol$Sample_id[cluster_esccol$ToolDB == "SR" & cluster_esccol$cluster_rep == "5"]))
u <- sort(unique(cluster_esccol$Sample_id[cluster_esccol$ToolDB == "SR" & cluster_esccol$cluster_rep == "6"]))

cluster_esccol$Sample_id <- factor(cluster_esccol$Sample_id, levels = c(v,q,r,s,t,u))

colors <- c( "blue4", "orangered", "cyan", "magenta", "yellowgreen", "moccasin")

cluster_esccol_heatmap <- ggplot(data = cluster_esccol, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("TEM ESBL type", "SHV ESBL type", "CMY AmpC type", "IMP-66", "mixed", "none")) +
  labs(x = "", y = "") +
  labs(fill = "ESBL/Carbapenemase genes") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(cluster_esccol$CeftriaxoneMIC, levels = c("<= 1", "= 2","= 4", "= 8", ">= 16", ">= 32", ">= 64", "none","NA"))), vars(Tool), scales = "free", space = "free")

cluster_esccol_heatmap

#-------------------------------------------------------------------------------
# concordance

bla_for_heatmap <-bla_gene_clusters

# classification of the found bla genes
carba <- read.table('../EC/carba.txt')
esbl <- read.table('../EC/esbl.txt')
#ampc <- read.table('../EC/ampc.txt')


for (i in 1:length(bla_for_heatmap$bla_genes)) { 
  ifelse(bla_for_heatmap$bla_genes[i] %in% carba$V1, bla_for_heatmap$Classification[i] <- 'Carbapenem', 
         ifelse(bla_for_heatmap$bla_genes[i] %in% esbl$V1, bla_for_heatmap$Classification[i] <- 'ESBL',
                #ifelse(bla_for_heatmap$bla_genes[i] %in% ampc$V1, bla_for_heatmap$Classification[i] <- 'AmpC', 
                bla_for_heatmap$Classification[i] <- 'Beta-lactam'))
}

carb_correlation <- pivot_wider(bla_for_heatmap[,c(1,3:17)], names_from = Tool, values_from = Classification)

concordance_table <- carb_correlation

# Concordance values: esccolc/Carbapenem = 1, esccolc/ESBL = 2, esccolc/none = 3, esccole/Carbapenem = 4, esccole/ESBL = 5
#                     esccole/no Carba or ESLB = 6, esccol/Carbapenem = 7, esccol/ESBL = 8, esccol/no Carba or ESLB = 9
concordance_table$abr_ncbi <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$abr_ncbi) == "TRUE", 2,
                                            ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$abr_ncbi) == "FALSE" & grepl("ESBL", concordance_table$abr_ncbi) == "FALSE"), 3,   
                                                   ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 4, 
                                                          ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$abr_ncbi) == "TRUE", 5, 
                                                                 ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$abr_ncbi) == "FALSE" & grepl("ESBL", concordance_table$abr_ncbi) == "FALSE"), 6,
                                                                        ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 7, 
                                                                               ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$abr_ncbi) == "TRUE", 8, 9))))))))

concordance_table$abr_card <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$abr_card) == "TRUE", 2,
                                            ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$abr_card) == "FALSE" & grepl("ESBL", concordance_table$abr_card) == "FALSE"), 3,   
                                                   ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 4, 
                                                          ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$abr_card) == "TRUE", 5, 
                                                                 ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$abr_card) == "FALSE" & grepl("ESBL", concordance_table$abr_card) == "FALSE"), 6,
                                                                        ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 7, 
                                                                               ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$abr_card) == "TRUE", 8, 9))))))))

concordance_table$abr_megares <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 1,
                                        ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$abr_megares) == "TRUE", 2,
                                               ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$abr_megares) == "FALSE" & grepl("ESBL", concordance_table$abr_megares) == "FALSE"), 3,   
                                                      ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 4, 
                                                             ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$abr_megares) == "TRUE", 5, 
                                                                    ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$abr_megares) == "FALSE" & grepl("ESBL", concordance_table$abr_megares) == "FALSE"), 6,
                                                                           ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 7, 
                                                                                  ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$abr_megares) == "TRUE", 8, 9))))))))

concordance_table$abr_resfinder <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 1,
                                          ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$abr_resfinder) == "TRUE", 2,
                                                 ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$abr_resfinder) == "FALSE" & grepl("ESBL", concordance_table$abr_resfinder) == "FALSE"), 3,   
                                                        ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 4, 
                                                               ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$abr_resfinder) == "TRUE", 5, 
                                                                      ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$abr_resfinder) == "FALSE" & grepl("ESBL", concordance_table$abr_resfinder) == "FALSE"), 6,
                                                                             ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 7, 
                                                                                    ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$abr_resfinder) == "TRUE", 8, 9))))))))

concordance_table$abr_argannot <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$abr_argannot) == "TRUE", 2,
                                                ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$abr_argannot) == "FALSE" & grepl("ESBL", concordance_table$abr_argannot) == "FALSE"), 3,   
                                                       ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 4, 
                                                              ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$abr_argannot) == "TRUE", 5, 
                                                                     ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$abr_argannot) == "FALSE" & grepl("ESBL", concordance_table$abr_argannot) == "FALSE"), 6,
                                                                            ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 7, 
                                                                                   ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$abr_argannot) == "TRUE", 8, 9))))))))

concordance_table$amrf_nuc <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$amrf_nuc) == "TRUE", 2,
                                            ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$amrf_nuc) == "FALSE" & grepl("ESBL", concordance_table$amrf_nuc) == "FALSE"), 3,   
                                                   ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 4, 
                                                          ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$amrf_nuc) == "TRUE", 5, 
                                                                 ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$amrf_nuc) == "FALSE" & grepl("ESBL", concordance_table$amrf_nuc) == "FALSE"), 6,
                                                                        ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 7, 
                                                                               ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$amrf_nuc) == "TRUE", 8, 9))))))))

concordance_table$amrf_prot <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 1,
                                      ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$amrf_prot) == "TRUE", 2,
                                             ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$amrf_prot) == "FALSE" & grepl("ESBL", concordance_table$amrf_prot) == "FALSE"), 3,   
                                                    ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 4, 
                                                           ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$amrf_prot) == "TRUE", 5, 
                                                                  ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$amrf_prot) == "FALSE" & grepl("ESBL", concordance_table$amrf_prot) == "FALSE"), 6,
                                                                         ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 7, 
                                                                                ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$amrf_prot) == "TRUE", 8, 9))))))))

concordance_table$deeparg_LS <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$deeparg_LS) == "TRUE", 2,
                                              ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$deeparg_LS) == "FALSE" & grepl("ESBL", concordance_table$deeparg_LS) == "FALSE"), 3,   
                                                     ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 4, 
                                                            ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$deeparg_LS) == "TRUE", 5, 
                                                                   ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$deeparg_LS) == "FALSE" & grepl("ESBL", concordance_table$deeparg_LS) == "FALSE"), 6,
                                                                          ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 7, 
                                                                                 ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$deeparg_LS) == "TRUE", 8, 9))))))))

concordance_table$deeparg_SR <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$deeparg_SR) == "TRUE", 2,
                                              ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$deeparg_SR) == "FALSE" & grepl("ESBL", concordance_table$deeparg_SR) == "FALSE"), 3,   
                                                     ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 4, 
                                                            ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$deeparg_SR) == "TRUE", 5, 
                                                                   ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$deeparg_SR) == "FALSE" & grepl("ESBL", concordance_table$deeparg_SR) == "FALSE"), 6,
                                                                          ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 7, 
                                                                                 ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$deeparg_SR) == "TRUE", 8, 9))))))))

concordance_table$resfinder_as <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$resfinder_as) == "TRUE", 2,
                                                ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$resfinder_as) == "FALSE" & grepl("ESBL", concordance_table$resfinder_as) == "FALSE"), 3,   
                                                       ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 4, 
                                                              ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$resfinder_as) == "TRUE", 5, 
                                                                     ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$resfinder_as) == "FALSE" & grepl("ESBL", concordance_table$resfinder_as) == "FALSE"), 6,
                                                                            ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 7, 
                                                                                   ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$resfinder_as) == "TRUE", 8, 9))))))))

concordance_table$resfinder_re <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$resfinder_re) == "TRUE", 2,
                                                ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$resfinder_re) == "FALSE" & grepl("ESBL", concordance_table$resfinder_re) == "FALSE"), 3,   
                                                       ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 4, 
                                                              ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$resfinder_re) == "TRUE", 5, 
                                                                     ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$resfinder_re) == "FALSE" & grepl("ESBL", concordance_table$resfinder_re) == "FALSE"), 6,
                                                                            ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 7, 
                                                                                   ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$resfinder_re) == "TRUE", 8, 9))))))))

concordance_table$rgi <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 1,
                                ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$rgi) == "TRUE", 2,
                                       ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$rgi) == "FALSE" & grepl("ESBL", concordance_table$rgi) == "FALSE"), 3,   
                                              ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 4, 
                                                     ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$rgi) == "TRUE", 5, 
                                                            ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$rgi) == "FALSE" & grepl("ESBL", concordance_table$rgi) == "FALSE"), 6,
                                                                   ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 7, 
                                                                          ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$rgi) == "TRUE", 8, 9))))))))

concordance_table$srax_basic <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$srax_basic) == "TRUE", 2,
                                              ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$srax_basic) == "FALSE" & grepl("ESBL", concordance_table$srax_basic) == "FALSE"), 3,   
                                                     ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 4, 
                                                            ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$srax_basic) == "TRUE", 5, 
                                                                   ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$srax_basic) == "FALSE" & grepl("ESBL", concordance_table$srax_basic) == "FALSE"), 6,
                                                                          ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 7, 
                                                                                 ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$srax_basic) == "TRUE", 8, 9))))))))

concordance_table$srax_ext <- ifelse(concordance_table$Type == "esccolc" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "esccolc" & grepl("ESBL", concordance_table$srax_ext) == "TRUE", 2,
                                            ifelse(concordance_table$Type == "esccolc" & (grepl("Carbapenem", concordance_table$srax_ext) == "FALSE" & grepl("ESBL", concordance_table$srax_ext) == "FALSE"), 3,   
                                                   ifelse(concordance_table$Type == "esccole" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 4, 
                                                          ifelse(concordance_table$Type == "esccole" & grepl("ESBL", concordance_table$srax_ext) == "TRUE", 5, 
                                                                 ifelse(concordance_table$Type == "esccole" & (grepl("Carbapenem", concordance_table$srax_ext) == "FALSE" & grepl("ESBL", concordance_table$srax_ext) == "FALSE"), 6,
                                                                        ifelse(concordance_table$Type == "esccol" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 7, 
                                                                               ifelse(concordance_table$Type == "esccol" & grepl("ESBL", concordance_table$srax_ext) == "TRUE", 8, 9))))))))

# transform data for heatmap
heatmap_concordance <- pivot_longer(data = concordance_table, 
                                    cols = -c(1:14),
                                    names_to = "ToolDB", 
                                    values_to = "Concordance")

heatmap_concordance$Concordance <- as.character(heatmap_concordance$Concordance)
heatmap_concordance[heatmap_concordance == "abr_ncbi"] <- "ncbi"
heatmap_concordance[heatmap_concordance == "abr_card"] <- "card"
heatmap_concordance[heatmap_concordance == "abr_megares"] <- "megares"
heatmap_concordance[heatmap_concordance == "abr_resfinder"] <- "resfinder"
heatmap_concordance[heatmap_concordance == "abr_argannot"] <- "argannot"
heatmap_concordance[heatmap_concordance == "amrf_nuc"] <- "nuc"
heatmap_concordance[heatmap_concordance == "amrf_prot"] <- "prot"
heatmap_concordance[heatmap_concordance == "srax_basic"] <- "basic"
heatmap_concordance[heatmap_concordance == "srax_ext"] <- "ext"
heatmap_concordance[heatmap_concordance == "deeparg_LS"] <- "LS"
heatmap_concordance[heatmap_concordance == "deeparg_SR"] <- "SR"
heatmap_concordance[heatmap_concordance == "resfinder_as"] <- "assembly"
heatmap_concordance[heatmap_concordance == "resfinder_re"] <- "reads"


# include column for the different tools
heatmap_concordance$Tool <- ifelse(heatmap_concordance$ToolDB == "rgi", "RGI",
                                   ifelse(heatmap_concordance$ToolDB == "assembly" | heatmap_concordance$ToolDB == "reads", "ResFinder",
                                          ifelse(heatmap_concordance$ToolDB == "LS" | heatmap_concordance$ToolDB == "SR", "DeepARG",
                                                 ifelse(heatmap_concordance$ToolDB == "basic" | heatmap_concordance$ToolDB == "ext", "sraX",
                                                        ifelse(heatmap_concordance$ToolDB == "nuc" | heatmap_concordance$ToolDB == "prot", "AMRFinder", "ABRicate")))))

# create heatmap
colors <- c("darkgreen", "red2", "red4", "gold", "green3", "orangered", "yellow", "goldenrod1", "green")
label <- c("Carbapenemase/Carbapenemase", "Carbapenemase/ESBL", "Carbapenemase/none", "ESBL/Carbapenemase", "ESBL/ESBL", "ESBL/none", "none/Carbapenemase", "none/ESBL", "none/none")

heatmap_concordance$Sample_id <- factor(heatmap_concordance$Sample_id, levels = c(a, b, c, d, e, f, g, h,z,j,k,l,m,n,o,p, v,q,r,s,t,u)) 

concordance_heatmap <- ggplot(data = heatmap_concordance, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = label) +
  labs(x = "", y = "") +
  labs(fill = "Reported/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance$Type, levels = c("esccol", "esccole", "esccolc"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap

### esccolc
heatmap_concordance_esccolc <- heatmap_concordance[c(which(heatmap_concordance$Type == "esccolc")),]
colors_c <- c("darkgreen", "red2", "red4")
label_c <- c("Carbapenemase/Carbapenemase", "Carbapenemase/ESBL", "Carbapenemase/none")

heatmap_concordance_esccolc$Sample_id <- factor(heatmap_concordance_esccolc$Sample_id, levels = c(a, b, c, d, e, f, g)) 

concordance_heatmap_esccolc <- ggplot(data = heatmap_concordance_esccolc, mapping = aes(x = ToolDB,
                                                                                        y = Sample_id,
                                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors_c, labels = label_c) +
  labs(x = "", y = "") +
  labs(fill = "Reported/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance_esccolc$MeropenemMIC, levels = c("<= 0.25", "<= 0.5", "= 0.75","<= 1", "= 2",
                                                                  "= 4", "= 8", ">= 16"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap_esccolc


### esccole
heatmap_concordance_esccole <- heatmap_concordance[c(which(heatmap_concordance$Type == "esccole")),]
colors_e <- c("gold", "green3", "orangered")
label_e <- c("ESBL/Carbapenemase", "ESBL/ESBL", "ESBL/none")

heatmap_concordance_esccole$Sample_id <- factor(heatmap_concordance_esccole$Sample_id, levels = c(h,z,j,k,l,m,n,o,p))

concordance_heatmap_esccole <- ggplot(data = heatmap_concordance_esccole, mapping = aes(x = ToolDB,
                                                                                        y = Sample_id,
                                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors_e, labels = label_e) +
  labs(x = "", y = "") +
  labs(fill = "Reported/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance_esccole$CeftriaxoneMIC, levels = c("<= 1", "= 4", "= 8", ">= 16", ">= 32", ">= 64", "none","NA"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap_esccole


### esccol
heatmap_concordance_esccol <- heatmap_concordance[c(which(heatmap_concordance$Type == "esccol")),]

heatmap_concordance_esccol$Sample_id <- factor(heatmap_concordance_esccol$Sample_id, levels = c(v,q,r,s,t,u))

colors_n <- c("yellow", "goldenrod1", "green")
label_n <- c("none/Carbapenemase", "none/ESBL", "none/none")

concordance_heatmap_esccol <- ggplot(data = heatmap_concordance_esccol, mapping = aes(x = ToolDB,
                                                                                      y = Sample_id,
                                                                                      fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors_n, labels = label_n) +
  labs(x = "", y = "") +
  labs(fill = "Reported/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(cluster_esccol$CeftriaxoneMIC, levels = c("<= 1", "= 2","= 4", "= 8", ">= 16", ">= 32", ">= 64", "none","NA"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap_esccol


#-------------------------------------------------------------------------------
# plot mic & phenotype

all_phenotypes <- pheno_mic[,c(1,2,3,32,39,117,146,153,6,28,33,120,142,147)]

all_phenotypes_table <- pivot_longer(data = all_phenotypes, 
                               cols = -c(1:2,6:8, 12:14),
                               names_to = "Antibiotic", 
                               values_to = "Phenotype")

colors <- c("gray", "lemonchiffon", "black", "red")
all_phenotypes_table$Phenotype <- factor(all_phenotypes_table$Phenotype, levels = c("S", "I", "R", "0")) 

all_phenotypes_heatmap <- ggplot(data = all_phenotypes_table, mapping = aes(x = Antibiotic,
                                                                y = Sample_id,
                                                                fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R", "NA")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(all_phenotypes_table$MIC.Meropenem, levels = c("<= 0.25", "<= 0.5", "= 0.75","<= 1", "= 2",
                                                                  "= 4", "= 8", ">= 16"))), scales = "free", space = "free")

all_phenotypes_heatmap

all_phenotypes_heatmap <- ggplot(data = all_phenotypes_table, mapping = aes(x = Antibiotic,
                                                                            y = Sample_id,
                                                                            fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R", "NA")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(LGM_BK), scales = "free", space = "free")

all_phenotypes_heatmap

# Carbapenems
## mero, imi, erta

mero_imi <- pheno_mic[,c(1,2,3,32,39,117,146,153)]
mero_imi <- mero_imi[c(which(mero_imi$LGM_BK == "esccolc")),]

mero_imi_table <- pivot_longer(data = mero_imi, 
                               cols = -c(1:2,6:8),
                               names_to = "Antibiotic", 
                               values_to = "Phenotype")


mero_imi_table$Sample_id <- factor(mero_imi_table$Sample_id, levels = c(a,b,c,d,e,f,g))

colors <- c("gray", "lemonchiffon", "black")
mero_imi_table$Phenotype <- factor(mero_imi_table$Phenotype, levels = c("S", "I", "R")) 

mero_imi_table$Sample_id <- factor(mero_imi_table$Sample_id, levels = c(a, b, c, d, e, f, g)) 

mero_imi_heatmap <- ggplot(data = mero_imi_table, mapping = aes(x = Antibiotic,
                                                                y = Sample_id,
                                                                fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R", "NA")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("<= 0.25", "<= 0.5", "= 0.75","<= 1", "= 2",
                                                                  "= 4", "= 8", ">= 16"))), scales = "free", space = "free")

mero_imi_heatmap

mero_imi_table$MIC.Meropenem <- factor(mero_imi_table$MIC.Meropenem, levels = c("<= 0.25", "<= 0.5", "= 0.75","<= 1", "= 2", "= 4", "= 8", ">= 16"))

mero_mic <- ggplot(data = mero_imi_table, mapping = aes(x = "Meropenem MIC",
                                                  y = Sample_id,
                                                  fill = MIC.Meropenem)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("<= 0.25", "<= 0.5", "= 0.75","<= 1", "= 2",
                                                            "= 4", "= 8", ">= 16"))), scales = "free", space = "free")

mero_mic


imi_mic <- ggplot(data = mero_imi_table, mapping = aes(x = "Imipenem MIC",
                                                 y = Sample_id,
                                                 fill = MIC.Imipenem)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("<= 0.25", "<= 0.5", "= 0.75","<= 1", "= 2",
                                                                  "= 4", "= 8", ">= 16"))), scales = "free", space = "free")

imi_mic


etp_mic <- ggplot(data = mero_imi_table, mapping = aes(x = "Ertapenem MIC",
                                                       y = Sample_id,
                                                       fill = MIC.Ertapenem)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("<= 0.25", "<= 0.5", "= 0.75","<= 1", "= 2",
                                                                  "= 4", "= 8", ">= 16"))), scales = "free", space = "free")

etp_mic

#-------------------------------------------------------------------------------
## ESBL - cephalosporine
cephalo <- pheno_mic[,c(1,2,6,28,33,120,142,147)]
cephalo <- cephalo[c(which(cephalo$LGM_BK == "esccole")),]

cephalo_table <- pivot_longer(data = cephalo, 
                               cols = -c(1:2,6:8),
                               names_to = "Antibiotic", 
                               values_to = "Phenotype")

colors <- c("gray", "lemonchiffon", "black")
cephalo_table$Phenotype <- factor(cephalo_table$Phenotype, levels = c("S", "I", "R")) 

cephalo_table$Sample_id <- factor(cephalo_table$Sample_id, levels = c(h,z,j,k,l,m,n,o,p))

hm_cephalo_table <- ggplot(data = cephalo_table, mapping = aes(x = Antibiotic,
                                                                y = Sample_id,
                                                                fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table$MIC.Ceftriaxone, levels = c("<= 1", "= 4", "= 8", ">= 16", ">= 32", ">= 64", "none","NA"))), scales = "free", space = "free")

hm_cephalo_table


## Ceftriaxon

heatmap_ctx <- ggplot(data = cephalo_table, mapping = aes(x = "Ceftriaxone MIC",
                                                          y = Sample_id,
                                                          fill = MIC.Ceftriaxone)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table$MIC.Ceftriaxone, levels = c("<= 1", "= 4", "= 8", ">= 16", ">= 32", ">= 64"))), scales = "free", space = "free")

heatmap_ctx

## ceftazidim

heatmap_caz <- ggplot(data = cephalo_table, mapping = aes(x = "Ceftazidime MIC",
                                                          y = Sample_id,
                                                          fill = MIC.Ceftazidim)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table$MIC.Ceftriaxon, levels = c("<= 1", "= 4", "= 8", ">= 16", ">= 32", ">= 64"))), scales = "free", space = "free")

heatmap_caz


## cefepim

heatmap_fep <- ggplot(data = cephalo_table, mapping = aes(x = "Cefepime MIC",
                                                          y = Sample_id,
                                                          fill = MIC.Cefepim)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table$MIC.Ceftriaxon, levels = c("<= 1", "= 4", "= 8", ">= 16", ">= 32", ">= 64"))), scales = "free", space = "free")

heatmap_fep


#-------------------------------------------------------------------------------
## none - cephalosporine

cephalo2 <- pheno_mic[,c(1,2,6,28,33,120,142,147)]
cephalo2 <- cephalo2[c(which(cephalo2$LGM_BK == "esccol")),]

cephalo_table2 <- pivot_longer(data = cephalo2, 
                              cols = -c(1:2,6:8),
                              names_to = "Antibiotic", 
                              values_to = "Phenotype")

colors <- c("gray", "lemonchiffon", "black", "red")
cephalo_table2$Phenotype <- factor(cephalo_table2$Phenotype, levels = c("S", "I", "R", "0")) 

cephalo_table2$Sample_id <- factor(cephalo_table2$Sample_id, levels = c(v,q,r,s,t,u))

hm_cephalo_table2 <- ggplot(data = cephalo_table2, mapping = aes(x = Antibiotic,
                                                               y = Sample_id,
                                                               fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R", "NA")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table2$MIC.Ceftriaxone, levels = c("<= 1", "= 2","= 4", "= 8", ">= 16", ">= 32", ">= 64", "none","NA"))), scales = "free", space = "free")

hm_cephalo_table2


## Ceftriaxon

heatmap_ctx2 <- ggplot(data = cephalo_table2, mapping = aes(x = "Ceftriaxone MIC",
                                                          y = Sample_id,
                                                          fill = MIC.Ceftriaxone)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table2$MIC.Ceftriaxone, levels = c("<= 1", "= 2","= 4", "= 8", ">= 16", ">= 32", ">= 64"))), scales = "free", space = "free")

heatmap_ctx2
