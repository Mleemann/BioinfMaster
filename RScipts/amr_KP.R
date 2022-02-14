#install.packages('tidyverse')
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)

setwd("~/data_project/assemblies/KP")

# read in samples to be analyzed
samples <- read.table('samples_KP.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_KP.csv', header = TRUE, sep = ';')
mic <- read.csv2('mic_KP.csv', header = TRUE, sep = ';')
pheno_mic <- read.csv2('pheno_mic_KP.csv', header = TRUE, sep = ';')

phenotypes[phenotypes == "klepneco"] <- "klepne"
mic[mic == "klepneco"] <- "klepne"
pheno_mic[pheno_mic == "klepneco"] <- "klepne"


# check sample ids; and correct them 
which(phenotypes$Sample_id %in% samples$V1 == FALSE)
phenotypes$Sample_id[21] <- "700505-17-wh"
phenotypes$Sample_id[38] <- "804158-2-17"
phenotypes$Sample_id[58] <- "807949-2-17"
phenotypes$Sample_id[63] <- "808716-18-17"
phenotypes$Sample_id[76] <- "803993-4-18"
phenotypes$Sample_id[96] <- "612633-12-18"
phenotypes$Sample_id[144] <- "133066-1-20"
phenotypes$Sample_id[145] <- "133506-1-20"
phenotypes$Sample_id[146] <- "608151-4-20"
phenotypes$Sample_id[171] <- "801570-5-21"
phenotypes$Sample_id[181] <- "802728-2-21"
phenotypes$Sample_id[188] <- "712842-21-wh"

#which(mic$Sample_id %in% samples$V1 == FALSE)
mic$Sample_id[21] <- "700505-17-wh"
mic$Sample_id[38] <- "804158-2-17"
mic$Sample_id[58] <- "807949-2-17"
mic$Sample_id[63] <- "808716-18-17"
mic$Sample_id[76] <- "803993-4-18"
mic$Sample_id[96] <- "612633-12-18"
mic$Sample_id[144] <- "133066-1-20"
mic$Sample_id[145] <- "133506-1-20"
mic$Sample_id[146] <- "608151-4-20"
mic$Sample_id[171] <- "801570-5-21"
mic$Sample_id[181] <- "802728-2-21"
mic$Sample_id[188] <- "712842-21-wh"

which(pheno_mic$Sample_id %in% samples$V1 == FALSE)
pheno_mic$Sample_id[5] <- "133066-1-20"
pheno_mic$Sample_id[7] <- "133506-1-20"
pheno_mic$Sample_id[31] <- "608151-4-20"
pheno_mic$Sample_id[34] <- "612633-12-18"
pheno_mic$Sample_id[36] <- "700505-17-wh"
pheno_mic$Sample_id[63] <- "712842-21-wh"
pheno_mic$Sample_id[110] <- "801570-5-21"
pheno_mic$Sample_id[121] <- "802728-2-21"
pheno_mic$Sample_id[135] <- "803993-4-18"
pheno_mic$Sample_id[138] <- "804158-2-17"
pheno_mic$Sample_id[191] <- "807949-2-17"
pheno_mic$Sample_id[202] <- "808716-18-17"

#unique(as.vector(as.matrix(pheno_mic)))
pheno_mic[pheno_mic == "c(R, R)" | pheno_mic == "c(S, R)" | pheno_mic == "c(R, S)" | pheno_mic == "c(S, U, R)" |   
            pheno_mic == "c(R, S, S)" | pheno_mic == "c(R, R, R)" | pheno_mic == "c(R, U)" | pheno_mic == "c(R, R, U)" |
            pheno_mic == "c(U, R)" | pheno_mic == "c(U, R, U, R)" | pheno_mic == "c(R, R, I)" | pheno_mic == "c(R, U, R, U)" |
            pheno_mic == "c(R, U, U, I, I)" | pheno_mic == "c(I, U, U, R)" | pheno_mic == "c(R, U, U, I)" | 
            pheno_mic == "c(U, R, R)" | pheno_mic == "c(I, R)"] <- "R"

pheno_mic[pheno_mic == "c(S, S)" | pheno_mic == "c(U, S)" | pheno_mic == "c(S, S, S)" | pheno_mic == "c(S, U)" |
            pheno_mic == "c(S, U, S)" | pheno_mic == "c(S, S, U)" | pheno_mic == "c(S, U, U, S)" |
            pheno_mic == "c(U, S, S)" | pheno_mic == "c(S, U, S, U)"] <- "S"

pheno_mic[pheno_mic == "c(I, S)" | pheno_mic == "c(S, I)" | pheno_mic == "c(I, S, S)" |
            pheno_mic == "c(S, I, S, I)" | pheno_mic == "c(I, I)" | pheno_mic == "c(I, U)" | 
            pheno_mic == "c(U, I)" | pheno_mic == "c(U, I, U, U)" | pheno_mic == "c(U, I, U, S)" |
            pheno_mic == "c(U, S, U, I)" | pheno_mic == "c(U, I, U)" | pheno_mic == "c(I, I, U)" |
            pheno_mic == "c(I, U, U)" | pheno_mic == "c(U, S, I)" | pheno_mic == "c(I, U, I)"] <- "I"

pheno_mic[pheno_mic == "U" | pheno_mic == "c(U, U)" | pheno_mic == "c(U, U, U)"] <- 0


# combine antibiotics with several columns
pheno_mic$Meropenem[pheno_mic$Meropenem.bei.Meningitis  == "S"] <- "S"
pheno_mic$Meropenem[pheno_mic$Meropenem.bei.Meningitis  == "R"] <- "R"
pheno_mic$Meropenem[pheno_mic$Meropenem.bei.Meningitis  == "I"] <- "I"

pheno_mic$Ceftriaxon[pheno_mic$Ceftriaxon.bei.Meningitis  == "S"] <- "S"
pheno_mic$Ceftriaxon[pheno_mic$Ceftriaxon.bei.Meningitis  == "R"] <- "R"
pheno_mic$Ceftriaxon[pheno_mic$Ceftriaxon.bei.Meningitis  == "I"] <- "I"

for (i in 1:length(pheno_mic$Sample_id)){
  if(pheno_mic$MIC.Meropenem[i] == "NULL") {
    pheno_mic$MIC.Meropenem[i] <- pheno_mic$MIC.Meropenem.ohne.Meningitis[i]
  }
}

pheno_mic[pheno_mic == "c(n.a. 0, = 0)" | pheno_mic == "n.a. 0" | pheno_mic == "c(n.a. 0, n.a. 0)"] <- "NA" 
pheno_mic[pheno_mic == "c(= 0.016, n.a. 0)"] <- "= 0.016"
pheno_mic[pheno_mic == "c(<= 0.25, <= 0.25)" | pheno_mic == "c(= 0.047, <= 0.25)" | pheno_mic == "c(<= 0.25, <= 0.25, <= 0.25)"] <-  "<= 0.25"
pheno_mic[pheno_mic == "c(<= 0.5, <= 0.5)" | pheno_mic == "c(= 0.38, = 0.5)" | pheno_mic == "c(= 0, = 0.5)" | 
            pheno_mic == "c(<= 0.5, <= 0.5, <= 0.5)" | pheno_mic == "c(= 0.75, <= 0.5)" | pheno_mic == "c(= 0.5, = 0.19)" | pheno_mic == "= 0.5"] <-  "<= 0.5"
pheno_mic[pheno_mic == "c(= 0, = 0.75)"] <- "= 0.75"
pheno_mic[pheno_mic == "c(n.a. 0, <= 1)" | pheno_mic == "= 0" | pheno_mic == "c(<= 1, < 0.5)" | pheno_mic == "c(<= 1, = 1)" | 
            pheno_mic == "c(= 0, <= 1)" | pheno_mic == "c(<= 1, <= 1)" | pheno_mic == "c(= 1.5, <= 1)" | pheno_mic == "c(<= 1, = 6)" |
            pheno_mic == "c(<= 1, = 0, = 0, <= 1)" | pheno_mic == "c(<= 1, = 0)" | pheno_mic == "c(<= 1, <= 1, = 0)" |
            pheno_mic == "c(n.a. 0, = 1, = 1)" | pheno_mic == "c(= 0.5, = 1)" | pheno_mic == "c(= 0.5, <= 1)" | pheno_mic == "c(>= 8, = 1)" |
            pheno_mic == "c(= 1, = 0.25)" | pheno_mic == "= 1"] <- "<= 1" 
pheno_mic[pheno_mic == "c(= 0, = 1.5)" ] <- "= 1.5" 
pheno_mic[pheno_mic == "c(= 2, = 0)" | pheno_mic == "c(<= 1, = 2)" | pheno_mic == "c(n.a. 0, = 2)" | pheno_mic == "c(= 2, = 0, n.a. 0)" | 
            pheno_mic == "c(= 2, = 0.25, = 0.19)" | pheno_mic == "c(= 2, = 0.25)" | pheno_mic == "c(= 2, = 0, n.a. 0, = 12)" |
            pheno_mic == "c(= 0, <= 1, n.a. 0, = 2)" | pheno_mic == "c(= 1, = 2)" | pheno_mic == "c(= 0, = 2, = 0)" |
            pheno_mic == "c(= 2, = 2)" | pheno_mic == "c(n.a. 0, <= 1, = 2)"] <- "= 2"
pheno_mic[pheno_mic == "c(= 0, = 3)" ] <- "= 3" 
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
            pheno_mic == "c(= 16, = 12)" | pheno_mic == "c(= 16, = 24)" | pheno_mic == "c(>= 16, = 24)" | pheno_mic == "c(> 32, >= 16)"] <- ">= 16" 
pheno_mic[pheno_mic == "= 32" | pheno_mic == "c(= 32, n.a. 0)" | pheno_mic == "c(n.a. 0, = 32)" | pheno_mic == "c(>= 8, > 32)"] <- ">= 32" 
pheno_mic[pheno_mic == "c(> 32, >= 64)" | pheno_mic == "c(>= 64, > 32)" | pheno_mic == "c(>= 64, >= 64)" 
          | pheno_mic == "c(>= 64, >= 64, >= 64)" | pheno_mic == "c(>= 64, = 32)" | pheno_mic == "c(>= 64, <= 1)" | 
            pheno_mic == "c(>= 64, n.a. 0)" | pheno_mic == "c(n.a. 0, >= 64)" | pheno_mic == "c(= 0, >= 64)" |
            pheno_mic == "c(>= 64, = 0, n.a. 0, = 2)" | pheno_mic == "c(>= 64, = 0, n.a. 0, = 2, = 2)" | pheno_mic == "c(>= 64, = 0)"] <- ">= 64" 


for (i in 1:length(pheno_mic$Sample_id)){
  if(pheno_mic$MIC.Ceftriaxon[i] == "NA") {
    pheno_mic$MIC.Ceftriaxon[i] <- pheno_mic$MIC.Ceftriaxon.ohne.Meningitis[i]
  }
}

#-------------------------------------------------------------------------------
# get result tables
# abricate
abr_ncbi <- read.delim('KP_abricate_ncbi.tab', sep = "\t", header = T)
if (length(which(abr_ncbi$X.IDENTITY < 95)) != 0) {
  abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$X.IDENTITY < 95)),]
}
abr_ncbi$X.FILE <- gsub(".fna","", abr_ncbi$X.FILE)
names(abr_ncbi)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_ncbi$Tool <- 'abr_ncbi'

abr_card <- read.delim('KP_abricate_card.tab', sep = "\t", header = T)
if (length(which(abr_card$X.IDENTITY < 95)) != 0) {
  abr_card <- abr_card[-c(which(abr_card$X.IDENTITY < 95)),]
}
abr_card$X.FILE <- gsub(".fna","", abr_card$X.FILE)
names(abr_card)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_card$Tool <- 'abr_card'

abr_megares <- read.delim('KP_abricate_megares.tab', sep = "\t", header = T)
if (length(which(abr_megares$X.IDENTITY < 95)) != 0) {
  abr_megares <- abr_megares[-c(which(abr_megares$X.IDENTITY < 95)),]
}
abr_megares$X.FILE <- gsub(".fna","", abr_megares$X.FILE)
names(abr_megares)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_megares$Tool <- 'abr_megares'

abr_resfinder <- read.delim('KP_abricate_resfinder.tab', sep = "\t", header = T)
if (length(which(abr_resfinder$X.IDENTITY < 95)) != 0) {
  abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$X.IDENTITY < 95)),]
}
abr_resfinder$X.FILE <- gsub(".fna","", abr_resfinder$X.FILE)
names(abr_resfinder)[c(1,14)] <- c('Sample_id', 'amr_genes')
abr_resfinder$Tool <- 'abr_resfinder'

abr_argannot <- read.delim('KP_abricate_argannot.tab', sep = "\t", header = T)
if (length(which(abr_argannot$X.IDENTITY < 95)) != 0) {
  abr_argannot <- abr_argannot[-c(which(abr_argannot$X.IDENTITY < 95)),]
}
abr_argannot$X.FILE <- gsub(".fna","", abr_argannot$X.FILE)
names(abr_argannot)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_argannot$Tool <- 'abr_argannot'


#amrfinder
amrf_nuc <- read.delim('KP_amrfinder_nuc.tab', sep = "\t")
if (length(which(amrf_nuc$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_nuc)[1] <- 'Sample_id'
amrf_nuc$Sample_id <- gsub("gnl\\|USB\\|","", amrf_nuc$Contig.id)
amrf_nuc$Sample_id <- gsub('(_\\d*)',"", amrf_nuc$Sample_id)
names(amrf_nuc)[6] <- 'amr_genes'
amrf_nuc$Tool <- 'amrf_nuc'

amrf_prot <- read.delim('KP_amrfinder_prot.tab', sep = "\t")
if (length(which(amrf_prot$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_prot <- amrf_prot[-c(which(amrf_prot$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_prot)[c(1,2)] <- c('Sample_id', 'amr_genes')
amrf_prot$Sample_id <- gsub('(_\\d*)',"", amrf_prot$Sample_id)
amrf_prot$Tool <- 'amrf_prot'


# rgi
rgi <- read.delim(('KP_rgi.tab'), sep = "\t")
if (length(which(rgi$Best_Identities < 95)) != 0) {
  rgi <- rgi[-c(which(rgi$Best_Identities < 95)),]
}
names(rgi)[c(1,2)] <- c('Sample_id', 'amr_genes')
rgi$Sample_id <- gsub("gnl\\|USB\\|","", rgi$Sample_id)
rgi$Sample_id <- gsub('(_\\d*)',"", rgi$Sample_id)
rgi$Sample_id <- gsub('\\s+', '', rgi$Sample_id)
rgi$Tool <- 'rgi'

# sraX
srax_basic <- read.delim(('KP_srax_basic.tab'), sep = "\t")
if (length(which(srax_basic$Identity_p < 95)) != 0) {
  srax_basic <- srax_basic[-c(which(srax_basic$Identity_p < 95)),]
}
names(srax_basic)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_basic$Sample_id <- gsub(".fna","", srax_basic$Sample_id)
srax_basic$Tool <- 'sraX_basic'

srax_ext <- read.delim(('KP_srax_ext.tab'), sep = "\t")
if (length(which(srax_ext$Identity_p < 95)) != 0) {
  srax_ext <- srax_ext[-c(which(srax_ext$Identity_p < 95)),]
}
names(srax_ext)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_ext$Sample_id <- gsub(".fna","", srax_ext$Sample_id)
srax_ext$Tool <- 'sraX_ext'


# Deeparg
deeparg_LS <- read.delim(('KP_deeparg_LS_noPot.tab'), sep = "\t")
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

deeparg_SR <- read.delim(('KP_deeparg_SR_noPot.tab'), sep = "\t")
if (length(which(deeparg_SR$identity < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$identity < 95)),]
}
if (length(which(deeparg_SR$probability < 0.8)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$probability < 0.8)),]
}
deeparg_SR_unique <- deeparg_SR %>% distinct(X.ARG, best.hit, .keep_all = TRUE)
names(deeparg_SR_unique)[1] <- 'Sample_id'
deeparg_SR_unique$Sample_id <- gsub("klepne/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
deeparg_SR_unique$Sample_id <- gsub("klepnee/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
deeparg_SR_unique$Sample_id <- gsub("klepnec/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
deeparg_SR_unique$Sample_id <- gsub("klepneco/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
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
resfinder_as <- read.delim(('KP_resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$Identity < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$Identity < 95)),]
}
names(resfinder_as)[c(6,1)] <- c('Sample_id', 'amr_genes')
resfinder_as$Sample_id <- gsub('gnl\\|USB\\|','', resfinder_as$Sample_id)
resfinder_as$Sample_id <- gsub('(_\\d*)','', resfinder_as$Sample_id)
resfinder_as$Tool <- 'resfinder_as'

resfinder_re <- read.delim(('KP_resfinder_reads.tab'), sep = "\t", header = F)
resfinder_re <- resfinder_re[,1:10]
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
  bla$CeftriaxonMIC[i] <- pheno_mic$MIC.Ceftriaxon[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$CeftazidimMIC[i] <- pheno_mic$MIC.Ceftazidim[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$CefepimMIC[i] <- pheno_mic$MIC.Cefepim[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
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
  bla$Ceftriaxon[i] <- pheno_mic$Ceftriaxon[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$Ceftazidim[i] <- pheno_mic$Ceftazidim[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(bla$Sample_id)){
  bla$Cefepim[i] <- pheno_mic$Cefepim[which(bla$Sample_id[i] == pheno_mic$Sample_id)]
}


#-------------------------------------------------------------------------------
# gene cluster

bla_gene_clusters <- bla
bla_gene_clusters$bla_genes <- gsub("\\(Bla\\)", "", bla_gene_clusters$bla_genes)
bla_gene_clusters$bla_genes <- gsub("bla", "", bla_gene_clusters$bla_genes)

bla_gene_clusters[bla_gene_clusters == "KPC-1" | bla_gene_clusters == "KPC-2"] <- "KPC-1/2"


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
                                                ifelse(gene_cluster_table$ToolDB == "deeparg_LS" | gene_cluster_table$ToolDB == "deeparg_SR", gene_cluster_table$Tool <- "deepARG",
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


### clusters klepnec
cluster_klepnec <- gene_cluster_table[c(which(gene_cluster_table$Type == "klepnec")),]


cluster_klepnec$cluster_rep <- ifelse(grepl("KPC-1/2",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 1, 
                                      ifelse(grepl("KPC-3",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 2,
                                             ifelse(grepl("\\bKPC\\b",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 3,
                                                    ifelse(grepl("\\bNDM-1\\b",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 4, 
                                                           ifelse(grepl("NDM-5",cluster_klepnec$gene_cluster) == "TRUE" & grepl("OXA-232",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 5,
                                                                  ifelse(grepl("NDM-5",cluster_klepnec$gene_cluster) == "TRUE" & grepl("OXA-181",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 6,
                                                                         ifelse(grepl("NDM-7",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 7,
                                                                                ifelse(grepl("NDM-19",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 8,
                                                                                       ifelse(grepl("\\bNDM\\b",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 9, 
                                                                                              ifelse(grepl("OXA-48",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 10,
                                                                                                     ifelse(grepl("OXA-181",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 11,
                                                                                                            ifelse(grepl("OXA-232",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 12,
                                                                                                                   ifelse(grepl("SHV-38",cluster_klepnec$gene_cluster) == "TRUE", cluster_klepnec$cluster_rep <- 13, 15)))))))))))))

for (i in 1:length(cluster_klepnec$Sample_id)) {
  if (cluster_klepnec$ToolDB[i] == "SR") {
    cluster_klepnec$cluster_rep[i] <- 14
  }
}

cluster_klepnec$cluster_rep <- as.character(cluster_klepnec$cluster_rep)

a <- sort(unique(cluster_klepnec$Sample_id[cluster_klepnec$ToolDB == "argannot" & cluster_klepnec$cluster_rep == "1"]))
b <- sort(unique(cluster_klepnec$Sample_id[cluster_klepnec$ToolDB == "argannot" & cluster_klepnec$cluster_rep == "2"]))
c <- sort(unique(cluster_klepnec$Sample_id[cluster_klepnec$ToolDB == "argannot" & (cluster_klepnec$cluster_rep == "4")]))
d <- sort(unique(cluster_klepnec$Sample_id[cluster_klepnec$ToolDB == "argannot" & (cluster_klepnec$cluster_rep == "6")]))
e <- sort(unique(cluster_klepnec$Sample_id[cluster_klepnec$ToolDB == "argannot" & (cluster_klepnec$cluster_rep == "8")]))
f <- sort(unique(cluster_klepnec$Sample_id[cluster_klepnec$ToolDB == "argannot" & (cluster_klepnec$cluster_rep == "10")]))
g <- sort(unique(cluster_klepnec$Sample_id[cluster_klepnec$ToolDB == "argannot" & (cluster_klepnec$cluster_rep == "11")]))
h <- sort(unique(cluster_klepnec$Sample_id[cluster_klepnec$ToolDB == "argannot" & (cluster_klepnec$cluster_rep == "15")]))

cluster_klepnec$Sample_id <- factor(cluster_klepnec$Sample_id, levels = c(a, b, c, d, e, f, g,h)) 

# order legend
cluster_klepnec$cluster_rep <- factor(cluster_klepnec$cluster_rep, levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) 
colors <- c("seagreen4", "seagreen3", "seagreen1", "blue4", "blue", "dodgerblue3", "deepskyblue2", "cyan2", "lightsteelblue", "orangered4", "orangered", "orange", "magenta", "yellowgreen", "moccasin")

cluster_klepnec_heatmap <- ggplot(data = cluster_klepnec, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("KPC-1/2", "KPC-3", "KPC", "NDM-1", "NDM-5 & OXA-232", "NDM-5 & OXA-181",
                                                "NDM-7", "NDM-19","NDM", "OXA-48", "OXA-181", "OXA-232", "SHV-38","mixed", "none")) +
  labs(x = "", y = "") +
  labs(fill = "Carbapenemase genes") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(cluster_klepnec$MeropenemMIC, levels = c("= 2", "= 4", "= 8", ">= 16", "> 32"))), vars(Tool), scales = "free", space = "free")

cluster_klepnec_heatmap


### clusters klepnee
cluster_klepnee <- gene_cluster_table[c(which(gene_cluster_table$Type == "klepnee")),]

ctxm1 <- c("\\bCTX-M-1\\b|\\bCTX-M-3\\b|\\bCTX-M-9\\b|\\bCTX-M-15\\b|\\bCTX-M-55\\b|\\bCTX-M-107\\b|\\bCTX-M-203\\b")
ctxm8 <- c("\\bCTX-M-8\\b|\\bCTX-M-63\\b")
ctxm9 <- c("\\bCTX-M-14\\b|\\bCTX-M-14b\\b|\\bCTX-M-24\\b|\\bCTX-M-27\\b|\\bCTX-M-122\\b|\\bCTX-M-125\\b|\\bCTX-M-125\\b")
shv <- c("\\bSHV-2\\b|\\bSHV-2a\\b|\\bSHV-2A\\b|\\bSHV-5\\b|\\bSHV-12\\b|\\bSHV-13\\b|\\bSHV-27\\b|\\bSHV-31\\b|\\bSHV-40\\b|\\bSHV-41\\b|\\bSHV-45\\b|\\bSHV-55\\b|\\bSHV-57\\b|\\bSHV-65\\b|\\bSHV-66\\b|\\bSHV-70\\b|\\bSHV-105\\b|\\bSHV-106\\b|\\bSHV-120\\b|\\bSHV-129\\b|\\bSHV-134\\b")
tem <- c("\\bTEM-7\\b|\\bTEM-20\\b|\\bTEM-47\\b|\\bTEM-112\\b")
cmy <- c("\\bCMY-42\\b")

cluster_klepnee$cluster_rep <- ifelse(grepl(ctxm1,cluster_klepnee$gene_cluster) == "TRUE" & grepl("\\bSHV-38\\b",cluster_klepnee$gene_cluster) == "TRUE", 1,
                                      ifelse(grepl(ctxm1,cluster_klepnee$gene_cluster) == "TRUE" & grepl(shv,cluster_klepnee$gene_cluster) == "TRUE", 2,
                                             ifelse(grepl(ctxm1,cluster_klepnee$gene_cluster) == "TRUE" & grepl(tem, cluster_klepnee$gene_cluster) == "TRUE", 3,
                                                    ifelse(grepl(ctxm1,cluster_klepnee$gene_cluster) == "TRUE", 4,
                                                           ifelse(grepl(ctxm9,cluster_klepnee$gene_cluster) == "TRUE" & grepl("OXA-48",cluster_klepnee$gene_cluster) == "TRUE", 11,
                                                                  ifelse(grepl(ctxm9,cluster_klepnee$gene_cluster) == "TRUE" & grepl(shv,cluster_klepnee$gene_cluster) == "TRUE", 5,
                                                                         ifelse(grepl(shv,cluster_klepnee$gene_cluster) == "TRUE", 7, 
                                                                                ifelse(grepl(ctxm8,cluster_klepnee$gene_cluster) == "TRUE", 8, 
                                                                                       ifelse(grepl(ctxm9,cluster_klepnee$gene_cluster) == "TRUE", 6, 
                                                                                              ifelse(grepl("\\CTX-M\\b",cluster_klepnee$gene_cluster) == "TRUE", 12, 10))))))))))


cluster_klepnee$cluster_rep <- as.character(cluster_klepnee$cluster_rep)

for (i in 1:length(cluster_klepnee$Sample_id)) {
  if (cluster_klepnee$ToolDB[i] == "SR") {
    cluster_klepnee$cluster_rep[i] <- 9
  }
}


# order legend
cluster_klepnee$cluster_rep <- factor(cluster_klepnee$cluster_rep, levels = c(4,1,2,3,8,6,11,5,12,7,9,10)) 

j <- sort(unique(cluster_klepnee$Sample_id[cluster_klepnee$ToolDB == "argannot" & cluster_klepnee$cluster_rep == "4"]))
k <- sort(unique(cluster_klepnee$Sample_id[cluster_klepnee$ToolDB == "argannot" & (cluster_klepnee$cluster_rep == "2")]))
l <- sort(unique(cluster_klepnee$Sample_id[cluster_klepnee$ToolDB == "argannot" & (cluster_klepnee$cluster_rep == "8")]))
m <- sort(unique(cluster_klepnee$Sample_id[cluster_klepnee$ToolDB == "argannot" & (cluster_klepnee$cluster_rep == "6")]))
n <- sort(unique(cluster_klepnee$Sample_id[cluster_klepnee$ToolDB == "argannot" & (cluster_klepnee$cluster_rep == "5")]))
o <- sort(unique(cluster_klepnee$Sample_id[cluster_klepnee$ToolDB == "argannot" & (cluster_klepnee$cluster_rep == "7")]))
p <- sort(unique(cluster_klepnee$Sample_id[cluster_klepnee$ToolDB == "argannot" & (cluster_klepnee$cluster_rep == "10")]))

cluster_klepnee$Sample_id <- factor(cluster_klepnee$Sample_id, levels = c(j,k,l,m,n,o,p))

colors <- c( "blue4", "gold", "dodgerblue3", "deepskyblue2", "cyan2", "orangered", "seagreen4", "orange", "maroon", "magenta", "yellowgreen","moccasin")

cluster_klepnee_heatmap <- ggplot(data = cluster_klepnee, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("CTX-M-1 Group", "CTX-M-1 Group & SHV Carbapenemase type", "CTX-M-1 Group & SHV ESBL type",
                                                "CTX-M-1 Group & TEM ESBL type", "CTX-M-8 Group", "CTX-M-9 Group", "CTX-M-9 Group & OXA-48",
                                                "CTX-M-9 Group & SHV ESBL type", "CTX-M", "SHV ESBL type", "mixed", "none")) +
  labs(x = "", y = "") +
  labs(fill = "ESBL genes") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(cluster_klepnee$CeftriaxonMIC, levels = c("<= 1", "= 8", ">= 16", ">= 32", ">= 64"))), vars(Tool), scales = "free", space = "free")

cluster_klepnee_heatmap


### clusters klepne
cluster_klepne <- gene_cluster_table[c(which(gene_cluster_table$Type == "klepne")),]

mbl <- c("metallo-beta-lactamase")
shv <- c("\\bSHV-2\\b|\\bSHV-2a\\b|\\bSHV-2A\\b|\\bSHV-5\\b|\\bSHV-12\\b|\\bSHV-13\\b|\\bSHV-27\\b|\\bSHV-31\\b|\\bSHV-40\\b|\\bSHV-41\\b|\\bSHV-45\\b|\\bSHV-55\\b|\\bSHV-57\\b|\\bSHV-65\\b|\\bSHV-66\\b|\\bSHV-70\\b|\\bSHV-105\\b|\\bSHV-106\\b|\\bSHV-120\\b|\\bSHV-129\\b|\\bSHV-134\\b")
tem <- c("\\bTEM-7\\b|\\bTEM-20\\b|\\bTEM-47\\b|\\bTEM-112\\b")
cmy <- c("CMY-8|CMY-9|CMY-30|CMY-37|CMY-42")

cluster_klepne$cluster_rep <- ifelse(grepl(mbl,cluster_klepne$gene_cluster) == "TRUE", 4,
                                     ifelse(grepl(tem,cluster_klepne$gene_cluster) == "TRUE", 1,
                                            ifelse(grepl(shv,cluster_klepne$gene_cluster) == "TRUE", 2, 
                                                   ifelse(grepl("\\bCTX-M-9\\b",cluster_klepne$gene_cluster) == "TRUE", 3,5))))

cluster_klepne$cluster_rep <- as.character(cluster_klepne$cluster_rep)

cluster_klepne$cluster_rep <- factor(cluster_klepne$cluster_rep, levels = c(1,2,3,4,5)) 

q <- sort(unique(cluster_klepne$Sample_id[cluster_klepne$ToolDB == "argannot" & cluster_klepne$cluster_rep == "2"]))
r <- sort(unique(cluster_klepne$Sample_id[cluster_klepne$ToolDB == "argannot" & cluster_klepne$cluster_rep == "3"]))
s <- sort(unique(cluster_klepne$Sample_id[cluster_klepne$ToolDB == "argannot" & cluster_klepne$cluster_rep == "5"]))

cluster_klepne$Sample_id <- factor(cluster_klepne$Sample_id, levels = c(q,r,s))

colors <- c( "blue4", "orangered", "seagreen4", "yellowgreen", "moccasin")

cluster_klepne_heatmap <- ggplot(data = cluster_klepne, mapping = aes(x = ToolDB,
                                                                      y = Sample_id,
                                                                      fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("TEM ESBL type", "SHV ESBL type", "OXA-9 Group", "mixed", "none")) +
  labs(x = "", y = "") +
  labs(fill = "ESBL/Carbapenemase genes") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(cluster_klepne$CeftriaxonMIC, levels = c("<= 1", "= 3","= 4", "= 8", ">= 16"))), vars(Tool), scales = "free", space = "free")

cluster_klepne_heatmap

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

# Concordance values: klepnec/Carbapenem = 1, klepnec/ESBL = 2, klepnec/none = 3, klepnee/Carbapenem = 4, klepnee/ESBL = 5
#                     klepnee/no Carba or ESLB = 6, klepne/Carbapenem = 7, klepne/ESBL = 8, klepne/no Carba or ESLB = 9
concordance_table$abr_ncbi <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$abr_ncbi) == "TRUE", 2,
                                            ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$abr_ncbi) == "FALSE" & grepl("ESBL", concordance_table$abr_ncbi) == "FALSE"), 3,   
                                                   ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 4, 
                                                          ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$abr_ncbi) == "TRUE", 5, 
                                                                 ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$abr_ncbi) == "FALSE" & grepl("ESBL", concordance_table$abr_ncbi) == "FALSE"), 6,
                                                                        ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 7, 
                                                                               ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$abr_ncbi) == "TRUE", 8, 9))))))))

concordance_table$abr_card <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$abr_card) == "TRUE", 2,
                                            ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$abr_card) == "FALSE" & grepl("ESBL", concordance_table$abr_card) == "FALSE"), 3,   
                                                   ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 4, 
                                                          ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$abr_card) == "TRUE", 5, 
                                                                 ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$abr_card) == "FALSE" & grepl("ESBL", concordance_table$abr_card) == "FALSE"), 6,
                                                                        ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 7, 
                                                                               ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$abr_card) == "TRUE", 8, 9))))))))

concordance_table$abr_megares <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 1,
                                        ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$abr_megares) == "TRUE", 2,
                                               ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$abr_megares) == "FALSE" & grepl("ESBL", concordance_table$abr_megares) == "FALSE"), 3,   
                                                      ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 4, 
                                                             ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$abr_megares) == "TRUE", 5, 
                                                                    ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$abr_megares) == "FALSE" & grepl("ESBL", concordance_table$abr_megares) == "FALSE"), 6,
                                                                           ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 7, 
                                                                                  ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$abr_megares) == "TRUE", 8, 9))))))))

concordance_table$abr_resfinder <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 1,
                                          ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$abr_resfinder) == "TRUE", 2,
                                                 ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$abr_resfinder) == "FALSE" & grepl("ESBL", concordance_table$abr_resfinder) == "FALSE"), 3,   
                                                        ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 4, 
                                                               ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$abr_resfinder) == "TRUE", 5, 
                                                                      ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$abr_resfinder) == "FALSE" & grepl("ESBL", concordance_table$abr_resfinder) == "FALSE"), 6,
                                                                             ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 7, 
                                                                                    ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$abr_resfinder) == "TRUE", 8, 9))))))))

concordance_table$abr_argannot <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$abr_argannot) == "TRUE", 2,
                                                ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$abr_argannot) == "FALSE" & grepl("ESBL", concordance_table$abr_argannot) == "FALSE"), 3,   
                                                       ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 4, 
                                                              ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$abr_argannot) == "TRUE", 5, 
                                                                     ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$abr_argannot) == "FALSE" & grepl("ESBL", concordance_table$abr_argannot) == "FALSE"), 6,
                                                                            ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 7, 
                                                                                   ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$abr_argannot) == "TRUE", 8, 9))))))))

concordance_table$amrf_nuc <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$amrf_nuc) == "TRUE", 2,
                                            ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$amrf_nuc) == "FALSE" & grepl("ESBL", concordance_table$amrf_nuc) == "FALSE"), 3,   
                                                   ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 4, 
                                                          ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$amrf_nuc) == "TRUE", 5, 
                                                                 ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$amrf_nuc) == "FALSE" & grepl("ESBL", concordance_table$amrf_nuc) == "FALSE"), 6,
                                                                        ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 7, 
                                                                               ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$amrf_nuc) == "TRUE", 8, 9))))))))

concordance_table$amrf_prot <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 1,
                                      ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$amrf_prot) == "TRUE", 2,
                                             ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$amrf_prot) == "FALSE" & grepl("ESBL", concordance_table$amrf_prot) == "FALSE"), 3,   
                                                    ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 4, 
                                                           ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$amrf_prot) == "TRUE", 5, 
                                                                  ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$amrf_prot) == "FALSE" & grepl("ESBL", concordance_table$amrf_prot) == "FALSE"), 6,
                                                                         ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 7, 
                                                                                ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$amrf_prot) == "TRUE", 8, 9))))))))

concordance_table$deeparg_LS <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$deeparg_LS) == "TRUE", 2,
                                              ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$deeparg_LS) == "FALSE" & grepl("ESBL", concordance_table$deeparg_LS) == "FALSE"), 3,   
                                                     ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 4, 
                                                            ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$deeparg_LS) == "TRUE", 5, 
                                                                   ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$deeparg_LS) == "FALSE" & grepl("ESBL", concordance_table$deeparg_LS) == "FALSE"), 6,
                                                                          ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 7, 
                                                                                 ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$deeparg_LS) == "TRUE", 8, 9))))))))

concordance_table$deeparg_SR <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$deeparg_SR) == "TRUE", 2,
                                              ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$deeparg_SR) == "FALSE" & grepl("ESBL", concordance_table$deeparg_SR) == "FALSE"), 3,   
                                                     ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 4, 
                                                            ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$deeparg_SR) == "TRUE", 5, 
                                                                   ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$deeparg_SR) == "FALSE" & grepl("ESBL", concordance_table$deeparg_SR) == "FALSE"), 6,
                                                                          ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 7, 
                                                                                 ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$deeparg_SR) == "TRUE", 8, 9))))))))

concordance_table$resfinder_as <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$resfinder_as) == "TRUE", 2,
                                                ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$resfinder_as) == "FALSE" & grepl("ESBL", concordance_table$resfinder_as) == "FALSE"), 3,   
                                                       ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 4, 
                                                              ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$resfinder_as) == "TRUE", 5, 
                                                                     ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$resfinder_as) == "FALSE" & grepl("ESBL", concordance_table$resfinder_as) == "FALSE"), 6,
                                                                            ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 7, 
                                                                                   ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$resfinder_as) == "TRUE", 8, 9))))))))

concordance_table$resfinder_re <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$resfinder_re) == "TRUE", 2,
                                                ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$resfinder_re) == "FALSE" & grepl("ESBL", concordance_table$resfinder_re) == "FALSE"), 3,   
                                                       ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 4, 
                                                              ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$resfinder_re) == "TRUE", 5, 
                                                                     ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$resfinder_re) == "FALSE" & grepl("ESBL", concordance_table$resfinder_re) == "FALSE"), 6,
                                                                            ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 7, 
                                                                                   ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$resfinder_re) == "TRUE", 8, 9))))))))

concordance_table$rgi <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 1,
                                ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$rgi) == "TRUE", 2,
                                       ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$rgi) == "FALSE" & grepl("ESBL", concordance_table$rgi) == "FALSE"), 3,   
                                              ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 4, 
                                                     ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$rgi) == "TRUE", 5, 
                                                            ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$rgi) == "FALSE" & grepl("ESBL", concordance_table$rgi) == "FALSE"), 6,
                                                                   ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 7, 
                                                                          ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$rgi) == "TRUE", 8, 9))))))))

concordance_table$srax_basic <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$srax_basic) == "TRUE", 2,
                                              ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$srax_basic) == "FALSE" & grepl("ESBL", concordance_table$srax_basic) == "FALSE"), 3,   
                                                     ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 4, 
                                                            ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$srax_basic) == "TRUE", 5, 
                                                                   ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$srax_basic) == "FALSE" & grepl("ESBL", concordance_table$srax_basic) == "FALSE"), 6,
                                                                          ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 7, 
                                                                                 ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$srax_basic) == "TRUE", 8, 9))))))))

concordance_table$srax_ext <- ifelse(concordance_table$Type == "klepnec" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "klepnec" & grepl("ESBL", concordance_table$srax_ext) == "TRUE", 2,
                                            ifelse(concordance_table$Type == "klepnec" & (grepl("Carbapenem", concordance_table$srax_ext) == "FALSE" & grepl("ESBL", concordance_table$srax_ext) == "FALSE"), 3,   
                                                   ifelse(concordance_table$Type == "klepnee" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 4, 
                                                          ifelse(concordance_table$Type == "klepnee" & grepl("ESBL", concordance_table$srax_ext) == "TRUE", 5, 
                                                                 ifelse(concordance_table$Type == "klepnee" & (grepl("Carbapenem", concordance_table$srax_ext) == "FALSE" & grepl("ESBL", concordance_table$srax_ext) == "FALSE"), 6,
                                                                        ifelse(concordance_table$Type == "klepne" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 7, 
                                                                               ifelse(concordance_table$Type == "klepne" & grepl("ESBL", concordance_table$srax_ext) == "TRUE", 8, 9))))))))

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
                                          ifelse(heatmap_concordance$ToolDB == "LS" | heatmap_concordance$ToolDB == "SR", "deepARG",
                                                 ifelse(heatmap_concordance$ToolDB == "basic" | heatmap_concordance$ToolDB == "ext", "sraX",
                                                        ifelse(heatmap_concordance$ToolDB == "nuc" | heatmap_concordance$ToolDB == "prot", "AMRFinder", "ABRicate")))))

# create heatmap
colors <- c("darkgreen", "red4", "gold", "green3", "orangered", "yellow", "goldenrod1", "green")
label <- c("Carbapenemase/Carbapenemase", "Carbapenemase/none", "ESBL/Carbapenemase", "ESBL/ESBL", "ESBL/none", "none/Carbapenemase", "none/ESBL", "none/none")

heatmap_concordance$Sample_id <- factor(heatmap_concordance$Sample_id, levels = c(a,b,c,d,e,f,g,h,j,k,l,m,n,o,p,q,r,s)) 

concordance_heatmap <- ggplot(data = heatmap_concordance, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = label) +
  labs(x = "", y = "") +
  labs(fill = "Reported/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance$Type, levels = c("klepne", "klepnee", "klepnec"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap

### klepnec
heatmap_concordance_klepnec <- heatmap_concordance[c(which(heatmap_concordance$Type == "klepnec")),]
colors_c <- c("darkgreen", "red4")
label_c <- c("Carbapenemase/Carbapenemase", "Carbapenemase/none")

heatmap_concordance_klepnec$Sample_id <- factor(heatmap_concordance_klepnec$Sample_id, levels = c(a, b, c, d, e, f, g, h)) 

concordance_heatmap_klepnec <- ggplot(data = heatmap_concordance_klepnec, mapping = aes(x = ToolDB,
                                                                                        y = Sample_id,
                                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors_c, labels = label_c) +
  labs(x = "", y = "") +
  labs(fill = "Reported/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance_klepnec$MeropenemMIC, levels = c("= 2", "= 4", "= 8", ">= 16", "> 32"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap_klepnec


### klepnee
heatmap_concordance_klepnee <- heatmap_concordance[c(which(heatmap_concordance$Type == "klepnee")),]
colors_e <- c("gold", "green3", "orangered")
label_e <- c("ESBL/Carbapenemase", "ESBL/ESBL", "ESBL/none")

heatmap_concordance_klepnee$Sample_id <- factor(heatmap_concordance_klepnee$Sample_id, levels = c(j,k,l,m,n,o,p))

concordance_heatmap_klepnee <- ggplot(data = heatmap_concordance_klepnee, mapping = aes(x = ToolDB,
                                                                                        y = Sample_id,
                                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors_e, labels = label_e) +
  labs(x = "", y = "") +
  labs(fill = "Reported/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance_klepnee$CeftriaxonMIC, levels = c("<= 1", "= 8", ">= 16", ">= 32", ">= 64"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap_klepnee


### klepne
heatmap_concordance_klepne <- heatmap_concordance[c(which(heatmap_concordance$Type == "klepne")),]

heatmap_concordance_klepne$Sample_id <- factor(heatmap_concordance_klepne$Sample_id, levels = c(q,r,s))

colors_n <- c("yellow", "goldenrod1", "green")
label_n <- c("none/Carbapenemase", "none/ESBL", "none/none")

concordance_heatmap_klepne <- ggplot(data = heatmap_concordance_klepne, mapping = aes(x = ToolDB,
                                                                                      y = Sample_id,
                                                                                      fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors_n, labels = label_n) +
  labs(x = "", y = "") +
  labs(fill = "Reported/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(cluster_klepne$CeftriaxonMIC, levels = c("<= 1", "= 3","= 4", "= 8", ">= 16"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap_klepne


#-------------------------------------------------------------------------------
# plot mic & phenotype

# Carbapenems
## mero, imi, erta

mero_imi <- pheno_mic[,c(1,2,5,29,40,112,136,147)]
mero_imi <- mero_imi[c(which(mero_imi$LGM_BK == "klepnec")),]

mero_imi_table <- pivot_longer(data = mero_imi, 
                               cols = -c(1:2,6:8),
                               names_to = "Antibiotic", 
                               values_to = "Phenotype")


mero_imi_table$Sample_id <- factor(mero_imi_table$Sample_id, levels = c(a,b,c,d,e,f,g,h))

colors <- c("gray", "lemonchiffon", "black", "red")
mero_imi_table$Phenotype <- factor(mero_imi_table$Phenotype, levels = c("S", "I", "R", "0")) 

mero_imi_heatmap <- ggplot(data = mero_imi_table, mapping = aes(x = Antibiotic,
                                                                y = Sample_id,
                                                                fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R", "NA")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("= 2", "= 4", "= 8", ">= 16", "> 32"))), scales = "free", space = "free")

mero_imi_heatmap

mero_imi_table$MIC.Meropenem <- factor(mero_imi_table$MIC.Meropenem, levels = c("= 2", "= 4", "= 8", ">= 16", "> 32"))

mero_mic <- ggplot(data = mero_imi_table, mapping = aes(x = "Meropenem MIC",
                                                        y = Sample_id,
                                                        fill = MIC.Meropenem)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("= 2", "= 4", "= 8", ">= 16", "> 32"))), scales = "free", space = "free")

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
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("<= 0.5", "= 0.75","<= 1", "= 2",
                                                                  "= 4", "= 8", ">= 16", "none","NA"))), scales = "free", space = "free")

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
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("<= 0.5", "= 0.75","<= 1", "= 2",
                                                                  "= 4", "= 8", ">= 16", "none","NA"))), scales = "free", space = "free")

etp_mic

#-------------------------------------------------------------------------------
## ESBL - cephalosporine
cephalo <- pheno_mic[,c(1,2,20,27,31,127,134,138)]
cephalo <- cephalo[c(which(cephalo$LGM_BK == "klepnee")),]

cephalo_table <- pivot_longer(data = cephalo, 
                              cols = -c(1:2,6:8),
                              names_to = "Antibiotic", 
                              values_to = "Phenotype")

colors <- c("gray", "lemonchiffon", "black", "red")
cephalo_table$Phenotype <- factor(cephalo_table$Phenotype, levels = c("S", "I", "R", "0")) 

cephalo_table$Sample_id <- factor(cephalo_table$Sample_id, levels = c(j,k,l,m,n,o,p))

hm_cephalo_table <- ggplot(data = cephalo_table, mapping = aes(x = Antibiotic,
                                                               y = Sample_id,
                                                               fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R", "NA")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table$MIC.Ceftriaxon, levels = c("<= 1", "= 8", ">= 16", ">= 32", ">= 64"))), scales = "free", space = "free")

hm_cephalo_table


## Ceftriaxon

heatmap_ctx <- ggplot(data = cephalo_table, mapping = aes(x = "Ceftriaxon MIC",
                                                          y = Sample_id,
                                                          fill = MIC.Ceftriaxon)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table$MIC.Ceftriaxon, levels = c("<= 1", "= 8", ">= 16", ">= 32", ">= 64"))), scales = "free", space = "free")

heatmap_ctx

## ceftazidim

heatmap_caz <- ggplot(data = cephalo_table, mapping = aes(x = "Ceftazidim MIC",
                                                          y = Sample_id,
                                                          fill = MIC.Ceftazidim)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table$MIC.Ceftriaxon, levels = c("<= 1", "= 8", ">= 16", ">= 32", ">= 64"))), scales = "free", space = "free")

heatmap_caz


## cefepim

heatmap_fep <- ggplot(data = cephalo_table, mapping = aes(x = "Cefepim MIC",
                                                          y = Sample_id,
                                                          fill = MIC.Cefepim)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table$MIC.Ceftriaxon, levels = c("<= 1", "= 8", ">= 16", ">= 32", ">= 64"))), scales = "free", space = "free")

heatmap_fep


#-------------------------------------------------------------------------------
## none - cephalosporine

cephalo2 <- pheno_mic[,c(1,2,20,27,31,127,134,138)]
cephalo2 <- cephalo2[c(which(cephalo2$LGM_BK == "klepne")),]

cephalo_table2 <- pivot_longer(data = cephalo2, 
                               cols = -c(1:2,6:8),
                               names_to = "Antibiotic", 
                               values_to = "Phenotype")

colors <- c("gray", "lemonchiffon", "black", "red")
cephalo_table2$Phenotype <- factor(cephalo_table2$Phenotype, levels = c("S", "I", "R", "0")) 

cephalo_table2$Sample_id <- factor(cephalo_table2$Sample_id, levels = c(q,r,s))

hm_cephalo_table2 <- ggplot(data = cephalo_table2, mapping = aes(x = Antibiotic,
                                                                 y = Sample_id,
                                                                 fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R", "NA")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table2$MIC.Ceftriaxon, levels = c("<= 1", "= 3","= 4", "= 8", ">= 16"))), scales = "free", space = "free")

hm_cephalo_table2


## Ceftriaxon

heatmap_ctx2 <- ggplot(data = cephalo_table2, mapping = aes(x = "Ceftriaxon MIC",
                                                            y = Sample_id,
                                                            fill = MIC.Ceftriaxon)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(cephalo_table2$MIC.Ceftriaxon, levels = c("<= 1", "= 3","= 4", "= 8", ">= 16"))), scales = "free", space = "free")

heatmap_ctx2
