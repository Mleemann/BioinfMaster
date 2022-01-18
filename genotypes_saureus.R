library(ggplot2)
library(tidyr)
install.packages('tidyverse')
library(tidyverse)

# read in samples to be analyzed
samples <- read.table('samples_SA.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_SA.csv', header = TRUE, sep = ';')
mic <- read.csv2('mic_SA.csv', header = TRUE, sep = ';')
pheno_mic <- read.csv2('pheno_mic_SA.csv', header = TRUE, sep = ';')


# check sample ids; and correct them 
#which(phenotypes$Sample_id %in% samples$V1 == FALSE)
phenotypes$Sample_id[90] <- "504423-19-wh"
phenotypes$Sample_id[161] <- "503081-1-20"
phenotypes$Sample_id[162] <- "402988-1-20"
phenotypes$Sample_id[163] <- "608029-1-20"
phenotypes$Sample_id[164] <- "503156-1-20"
phenotypes$Sample_id[165] <- "805048-1-20-wh"
phenotypes$Sample_id[179] <- "806818-20-wh"
phenotypes$Sample_id[219] <- "606179-21-wh"
phenotypes$Sample_id[233] <- "804519-01"
phenotypes$Sample_id[242] <- "608085-21-wh"
phenotypes$Sample_id[256] <- "205044-21-wh"
phenotypes$Sample_id[263] <- "808152-21-wh"
phenotypes$Sample_id[267] <- "404692-21-wh"
phenotypes$Sample_id[268] <- "808201-21-wh"

#which(mic$Sample_id %in% samples$V1 == FALSE)
mic$Sample_id[90] <- "504423-19-wh"
mic$Sample_id[161] <- "503081-1-20"
mic$Sample_id[162] <- "402988-1-20"
mic$Sample_id[163] <- "608029-1-20"
mic$Sample_id[164] <- "503156-1-20"
mic$Sample_id[165] <- "805048-1-20-wh"
mic$Sample_id[179] <- "806818-20-wh"
mic$Sample_id[219] <- "606179-21-wh"
mic$Sample_id[233] <- "804519-01"
mic$Sample_id[242] <- "608085-21-wh"
mic$Sample_id[256] <- "205044-21-wh"
mic$Sample_id[263] <- "808152-21-wh"
mic$Sample_id[267] <- "404692-21-wh"
mic$Sample_id[268] <- "808201-21-wh"

#which(pheno_mic$Sample_id %in% samples$V1 == FALSE)
pheno_mic$Sample_id[102] <- "504423-19-wh"
pheno_mic$Sample_id[86] <- "503081-1-20"
pheno_mic$Sample_id[50] <- "402988-1-20"
pheno_mic$Sample_id[134] <- "608029-1-20"
pheno_mic$Sample_id[87] <- "503156-1-20"
pheno_mic$Sample_id[226] <- "805048-1-20-wh"
pheno_mic$Sample_id[248] <- "806818-20-wh"
pheno_mic$Sample_id[123] <- "606179-21-wh"
pheno_mic$Sample_id[223] <- "804519-01"
pheno_mic$Sample_id[136] <- "608085-21-wh"
pheno_mic$Sample_id[27] <- "205044-21-wh"
pheno_mic$Sample_id[256] <- "808152-21-wh"
pheno_mic$Sample_id[69] <- "404692-21-wh"
pheno_mic$Sample_id[258] <- "808201-21-wh"


# get result tables
# abricate
abr_ncbi <- read.delim('abricate_ncbi.tab', sep = "\t", header = F)
if (length(which(abr_ncbi$V11 < 95)) != 0) {
  abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$V11 < 95)),]
}
abr_ncbi$V1 <- gsub(".fna","", abr_ncbi$V1)
names(abr_ncbi)[1] <- c('Sample_id')

abr_card <- read.delim('abricate_card.tab', sep = "\t", header = F)
if (length(which(abr_card$V11 < 95)) != 0) {
  abr_card <- abr_card[-c(which(abr_card$V11 < 95)),]
}
abr_card$V1 <- gsub(".fna","", abr_card$V1)
names(abr_card)[1] <- c('Sample_id')

abr_megares <- read.delim('abricate_megares.tab', sep = "\t", header = F)
if (length(which(abr_megares$V11 < 95)) != 0) {
  abr_megares <- abr_megares[-c(which(abr_megares$V11 < 95)),]
}
abr_megares$V1 <- gsub(".fna","", abr_megares$V1)
names(abr_megares)[1] <- c('Sample_id')

abr_resfinder <- read.delim('abricate_resfinder.tab', sep = "\t", header = F)
if (length(which(abr_resfinder$V11 < 95)) != 0) {
  abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$V11 < 95)),]
}
abr_resfinder$V1 <- gsub(".fna","", abr_resfinder$V1)
names(abr_resfinder)[1] <- c('Sample_id')

abr_argannot <- read.delim('abricate_argannot.tab', sep = "\t", header = F)
if (length(which(abr_argannot$V11 < 95)) != 0) {
  abr_argannot <- abr_argannot[-c(which(abr_argannot$V11 < 95)),]
}
abr_argannot$V1 <- gsub(".fna","", abr_argannot$V1)
names(abr_argannot)[1] <- c('Sample_id')

#amrfinder
amrf_nuc <- read.delim('amrfinder_nuc.tab', sep = "\t")
if (length(which(amrf_nuc$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_nuc)[1] <- 'Sample_id'
amrf_nuc$Sample_id <- gsub("gnl\\|USB\\|","", amrf_nuc$Contig.id)
amrf_nuc$Sample_id <- gsub('(_\\d*)',"", amrf_nuc$Sample_id)

amrf_prot <- read.delim('amrfinder_prot.tab', sep = "\t")
if (length(which(amrf_prot$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_prot <- amrf_prot[-c(which(amrf_prot$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_prot)[1] <- 'Sample_id'
amrf_prot$Sample_id <- gsub('(_\\d*)',"", amrf_prot$Sample_id)


# rgi
rgi <- read.delim(('rgi.tab'), sep = "\t")
if (length(which(rgi$Best_Identities < 95)) != 0) {
  rgi <- rgi[-c(which(rgi$Best_Identities < 95)),]
}
names(rgi)[1] <- 'Sample_id'
rgi$Sample_id <- gsub("gnl\\|USB\\|","", rgi$Sample_id)
rgi$Sample_id <- gsub('(_\\d*)',"", rgi$Sample_id)


# sraX
srax_basic <- read.delim(('srax_basic.tab'), sep = "\t")
if (length(which(srax_basic$Identity_p < 95)) != 0) {
  srax_basic <- srax_basic[-c(which(srax_basic$Identity_p < 95)),]
}
names(srax_basic)[2] <- 'Sample_id'
srax_basic$Sample_id <- gsub(".fna","", srax_basic$Sample_id)

srax_ext <- read.delim(('srax_ext.tab'), sep = "\t")
if (length(which(srax_ext$Identity_p < 95)) != 0) {
  srax_ext <- srax_ext[-c(which(srax_ext$Identity_p < 95)),]
}
names(srax_ext)[2] <- 'Sample_id'
srax_ext$Sample_id <- gsub(".fna","", srax_ext$Sample_id)


# Deeparg
deeparg_LS <- read.delim(('deeparg_LS.tab'), sep = "\t")
if (length(which(deeparg_LS$identity < 95)) != 0) {
  deeparg_LS <- deeparg_LS[-c(which(deeparg_LS$identity < 95)),]
}
names(deeparg_LS)[4] <- 'Sample_id'
deeparg_LS$Sample_id <- gsub("gnl\\|USB\\|","", deeparg_LS$Sample_id)
deeparg_LS$Sample_id <- gsub('(_\\d*)',"", deeparg_LS$Sample_id)
for(i in 1:length(deeparg_LS[,1])) {
  genes <- unlist(strsplit(deeparg_LS[i,6], "\\|"))
  deeparg_LS[i,6] <- genes[length(genes)]
}

deeparg_SR <- read.delim(('deeparg_SR.tab'), sep = "\t")
if (length(which(deeparg_SR$identity < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$identity < 95)),]
}
deeparg_SR_unique <- deeparg_SR %>% distinct(X.ARG, best.hit, .keep_all = TRUE)
names(deeparg_SR_unique)[1] <- 'Sample_id'
deeparg_SR_unique$Sample_id <- gsub("mrsa/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
deeparg_SR_unique$Sample_id <- gsub("staaur/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
for(i in 1:length(deeparg_SR_unique[,1])) {
  genes <- unlist(strsplit(deeparg_SR_unique[i,1], "\\_"))
  deeparg_SR_unique[i,1] <- genes[1]
}
for(i in 1:length(deeparg_SR_unique[,1])) {
  genes <- unlist(strsplit(deeparg_SR_unique[i,6], "\\|"))
  deeparg_SR_unique[i,6] <- genes[length(genes)]
}
deeparg_SR_unique <- deeparg_SR_unique[-c(which(deeparg_SR_unique$best.hit == "undefined")), ]


# resfinder
resfinder_as <- read.delim(('resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$Identity < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$Identity < 95)),]
}
names(resfinder_as)[6] <- 'Sample_id'
resfinder_as$Sample_id <- gsub('gnl\\|USB\\|','', resfinder_as$Sample_id)
resfinder_as$Sample_id <- gsub('(_\\d*)','', resfinder_as$Sample_id)

resfinder_re <- read.delim(('resfinder_reads.tab'), sep = "\t")
if (length(which(resfinder_re$Identity < 95)) != 0) {
  resfinder_re <- resfinder_re[-c(which(resfinder_re$Identity < 95)),]
}


# extract lines with methicillin resistance 
ncbi_meth <-abr_ncbi[c(which(abr_ncbi$V15 == "METHICILLIN")), ]
names(ncbi_meth)[6] <- 'abr_ncbi'
card_meth <-abr_card[c(which(abr_card$V6 == "mecA" | abr_card$V6 == "mecC")), ]
names(card_meth)[6] <- 'abr_card'
megares_meth <-abr_megares[c(which(abr_megares$V6 == "MECA" | abr_megares$V6 == "MECI" | abr_megares$V6 == "PBP2")), ]
names(megares_meth)[6] <- 'abr_megares'
resfinder_meth <-abr_resfinder[c(which(abr_resfinder$V14 == "mecA" | abr_resfinder$V14 == "mecC")), ]
names(resfinder_meth)[14] <- 'abr_resfinder'
argannot_meth <-abr_argannot[c(which(abr_argannot$V6 == "(Bla)mecA" | abr_argannot$V6 == "(Bla)mecC")), ]
names(argannot_meth)[6] <- 'abr_argannot'
amrf_nuc_meth <- amrf_nuc[c(which(amrf_nuc$Subclass == "METHICILLIN")), ]
names(amrf_nuc_meth)[6] <- 'amrf_nuc'
amrf_prot_meth <- amrf_prot[c(which(amrf_prot$Subclass == "METHICILLIN")), ]
names(amrf_prot_meth)[2] <- 'amrf_prot'
deeparg_LS_meth <- deeparg_LS[c(which(deeparg_LS$best.hit == "mecA")), ]
names(deeparg_LS_meth)[6] <- 'deeparg_LS'
deeparg_SR_meth <- deeparg_SR_unique[c(which(deeparg_SR_unique$best.hit == "mecA" | deeparg_SR_unique$best.hit == "mecR1" | deeparg_SR_unique$best.hit == "mecI", deeparg_SR_unique$best.hit == "mecC")), ]
names(deeparg_SR_meth)[6] <- 'deeparg_SR'
resfinder_as_meth <- resfinder_as[c(which(resfinder_as$Resistance.gene == "mecA" | resfinder_as$Resistance.gene == "mecC")), ]
names(resfinder_as_meth)[1] <- 'resfinder_as'
resfinder_re_meth <- resfinder_re[c(which(resfinder_re$Resistance.gene == "mecA" | resfinder_re$Resistance.gene == "mecC")), ]
names(resfinder_re_meth)[1] <- 'resfinder_re'
rgi_meth <- rgi[c(which(rgi$Drug.Class == "penam" & rgi$Best_Hit_ARO != "PC1 beta-lactamase (blaZ)" & rgi$Best_Hit_ARO != "mecC-type BlaZ")), ]
names(rgi_meth)[2] <- 'rgi'
srax_basic_meth <- srax_basic[c(which(srax_basic$AMR_Class == 'Penam' & srax_basic$AMR_gene != 'PC1_beta-lactamase_(blaZ)')), ]
names(srax_basic_meth)[6] <- 'sraX_basic'
srax_ext_meth <- srax_ext[c(which(srax_ext$AMR_Class == 'Penam' & srax_ext$AMR_gene != 'PC1_beta-lactamase_(blaZ)')), ]
names(srax_ext_meth)[6] <- 'sraX_ext'

meth_dfs <- list(ncbi_meth[,c(1,6)], card_meth[,c(1,6)], megares_meth[,c(1,6)], resfinder_meth[,c(1,14)], argannot_meth[,c(1,6)],
                 amrf_nuc_meth[,c(1,6)], amrf_prot_meth[,c(1,2)], deeparg_LS_meth[,c(4,6)], deeparg_SR_meth[,c(1,6)],
                 resfinder_as_meth[,c(6,1)], resfinder_re_meth[,c(10,1)], rgi_meth[,c(1,2)], srax_basic_meth[,c(2,6)], srax_ext_meth[,c(2,6)])
meth <- meth_dfs %>% reduce(full_join, by = 'Sample_id')

# extract oxacillin pattern
SA_ox <- pheno_mic[,c(1,2,27,93)]

# check if methicillin resistance found for sample_id and include it into the table --> TRUE present
SA_ox$meth_ncbi <- SA_ox$Sample_id %in% ncbi_meth$Sample_id
SA_ox$meth_card <- SA_ox$Sample_id %in% card_meth$Sample_id
SA_ox$meth_megares <- SA_ox$Sample_id %in% megares_meth$Sample_id
SA_ox$meth_resfinder <- SA_ox$Sample_id %in% resfinder_meth$Sample_id
SA_ox$meth_argannot <- SA_ox$Sample_id %in% argannot_meth$Sample_id
SA_ox$meth_amrf_nuc <- SA_ox$Sample_id %in% amrf_nuc_meth$Sample_id
SA_ox$meth_amrf_prot <- SA_ox$Sample_id %in% amrf_prot_meth$Sample_id
SA_ox$meth_deeparg_LS <- SA_ox$Sample_id %in% deeparg_LS_meth$Sample_id
SA_ox$meth_deeparg_SR <- SA_ox$Sample_id %in% deeparg_SR_meth$Sample_id
SA_ox$meth_resfinder_as <- SA_ox$Sample_id %in% resfinder_as_meth$Sample_id
SA_ox$meth_resfinder_re <- SA_ox$Sample_id %in% resfinder_re_meth$Sample_id
SA_ox$meth_rgi <- SA_ox$Sample_id %in% rgi_meth$Sample_id
SA_ox$meth_srax_basic <- SA_ox$Sample_id %in% srax_basic_meth$Sample_id
SA_ox$meth_srax_ext <- SA_ox$Sample_id %in% srax_ext_meth$Sample_id


#-------------------------------------------------------------------------------
# creating a heatmap 
for_heatmap_ox <- SA_ox
for_heatmap_ox[for_heatmap_ox == "c(R, R)"] <- "R"
for_heatmap_ox[for_heatmap_ox == "c(S, S)"] <- "S"
for_heatmap_ox[for_heatmap_ox == "c(U, R)"] <- "R"
for_heatmap_ox[for_heatmap_ox == "c(R, U)"] <- "R"
for_heatmap_ox[for_heatmap_ox == "c(R, S)"] <- "R"
for_heatmap_ox[for_heatmap_ox == "c(S, R)"] <- "R"

for_heatmap_ox <- for_heatmap_ox[-c(which(for_heatmap_ox$Oxacillin == 'U')), ]

# Concordance values: R/R = 1, S/S = 2, S/R = 3, R/S = 4
for_heatmap_ox$meth_ncbi <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_ncbi == "FALSE", 2,
                               ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_ncbi == "TRUE", 3,
                                      ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_ncbi == "TRUE", 1, 4)))

for_heatmap_ox$meth_card <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_card == "FALSE", 2,
                               ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_card == "TRUE", 3,
                                      ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_card == "TRUE", 1, 4)))

for_heatmap_ox$meth_megares <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_megares == "FALSE", 2,
                                   ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_megares == "TRUE", 3,
                                          ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_megares == "TRUE", 1, 4)))

for_heatmap_ox$meth_resfinder <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_resfinder == "FALSE", 2,
                                   ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_resfinder == "TRUE", 3,
                                          ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_resfinder == "TRUE", 1, 4)))

for_heatmap_ox$meth_argannot <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_argannot == "FALSE", 2,
                                   ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_argannot == "TRUE", 3,
                                          ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_argannot == "TRUE", 1, 4)))

for_heatmap_ox$meth_amrf_nuc <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_amrf_nuc == "FALSE", 2,
                                   ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_amrf_nuc == "TRUE", 3,
                                          ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_amrf_nuc == "TRUE", 1, 4)))

for_heatmap_ox$meth_amrf_prot <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_amrf_prot == "FALSE", 2,
                                    ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_amrf_prot == "TRUE", 3,
                                           ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_amrf_prot == "TRUE", 1, 4)))

for_heatmap_ox$meth_rgi <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_rgi == "FALSE", 2,
                              ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_rgi == "TRUE", 3,
                                     ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_rgi == "TRUE", 1, 4)))

for_heatmap_ox$meth_srax_basic <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_srax_basic == "FALSE", 2,
                                 ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_srax_basic == "TRUE", 3,
                                        ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_srax_basic == "TRUE", 1, 4)))

for_heatmap_ox$meth_srax_ext <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_srax_ext == "FALSE", 2,
                                 ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_srax_ext == "TRUE", 3,
                                        ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_srax_ext == "TRUE", 1, 4)))

for_heatmap_ox$meth_deeparg_LS <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_deeparg_LS == "FALSE", 2,
                                     ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_deeparg_LS == "TRUE", 3,
                                            ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_deeparg_LS == "TRUE", 1, 4)))

for_heatmap_ox$meth_deeparg_SR <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_deeparg_SR == "FALSE", 2,
                                            ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_deeparg_SR == "TRUE", 3,
                                                   ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_deeparg_SR == "TRUE", 1, 4)))

for_heatmap_ox$meth_resfinder_as <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_resfinder_as == "FALSE", 2,
                                       ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_resfinder_as == "TRUE", 3,
                                              ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_resfinder_as == "TRUE", 1, 4)))

for_heatmap_ox$meth_resfinder_re <- ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_resfinder_re == "FALSE", 2,
                                       ifelse(for_heatmap_ox$Oxacillin == "S" & for_heatmap_ox$meth_resfinder_re == "TRUE", 3,
                                              ifelse(for_heatmap_ox$Oxacillin == "R" & for_heatmap_ox$meth_resfinder_re == "TRUE", 1, 4)))

# transform data for heatmap
heatmap_ox <- pivot_longer(data = for_heatmap_ox, 
                              cols = -c(1:4),
                              names_to = "ToolDB", 
                              values_to = "Concordance")

heatmap_ox$Concordance <- as.character(heatmap_ox$Concordance)
heatmap_ox[heatmap_ox == "meth_ncbi"] <- "ncbi"
heatmap_ox[heatmap_ox == "meth_card"] <- "card"
heatmap_ox[heatmap_ox == "meth_megares"] <- "megares"
heatmap_ox[heatmap_ox == "meth_resfinder"] <- "resfinder"
heatmap_ox[heatmap_ox == "meth_argannot"] <- "argannot"
heatmap_ox[heatmap_ox == "meth_amrf_nuc"] <- "nuc"
heatmap_ox[heatmap_ox == "meth_amrf_prot"] <- "prot"
heatmap_ox[heatmap_ox == "meth_srax_basic"] <- "basic"
heatmap_ox[heatmap_ox == "meth_srax_ext"] <- "ext"
heatmap_ox[heatmap_ox == "meth_deeparg_LS"] <- "LS"
heatmap_ox[heatmap_ox == "meth_deeparg_SR"] <- "SR"
heatmap_ox[heatmap_ox == "meth_resfinder_as"] <- "assembly"
heatmap_ox[heatmap_ox == "meth_resfinder_re"] <- "reads"

# include column for the different tools
heatmap_ox$Tool <- ifelse(heatmap_ox$ToolDB == "meth_rgi", "rgi",
                             ifelse(heatmap_ox$ToolDB == "assembly" | heatmap_ox$ToolDB == "reads", "resfinder",
                                    ifelse(heatmap_ox$ToolDB == "LS" | heatmap_ox$ToolDB == "SR", "deeparg",
                                           ifelse(heatmap_ox$ToolDB == "basic" | heatmap_ox$ToolDB == "ext", "sraX",
                                                  ifelse(heatmap_ox$ToolDB == "nuc" | heatmap_ox$ToolDB == "prot", "amrfinder", "abricate")))))

# do it after step before otherwise card of rgi is counted as abricate
heatmap_ox[heatmap_ox == "meth_rgi"] <- "card"

# create heatmap
colors <- c("darkgreen", "green", "gold", "red3")

OX_heatmap <- ggplot(data = heatmap_ox, mapping = aes(x = ToolDB,
                                                         y = Sample_id,
                                                         fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "S/R", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype/Genotype") +
  theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 6, angle = 90)) +
  facet_grid(vars(LGM_BK), vars(Tool), scales = "free", space = "free")

OX_heatmap


#-------------------------------------------------------------------------------
# visualization of gene clusters per strain abricate ncbi
genes_abr_ncbi <- abr_ncbi[, c(1,6)]
genes_abr_ncbi$V3 <- 1
for(i in 1:length(genes_abr_ncbi[,1])) {
  genes_abr_ncbi[i,4] <- pheno_mic[which(genes_abr_ncbi[i,1] == pheno_mic$Sample_id), 2]
}

#gen_ord <- c('vanA', 'vanH-A', 'vanR-A', 'vanS-A', 'vanX-A', 'vanY-A', 'vanZ-A', 'vanB', 'vanH-B', 'vanR-B', 'vanS-B', 'vanW-B', 'vanX-B', 'vanY-B')
#for(i in gen_ord) {
#  van_genes$V3[van_genes$V6 == i] <- which(gen_ord == i)
#}
genes_abr_ncbi <- genes_abr_ncbi[-c(which(is.na(genes_abr_ncbi$V6))), ]
genes_abr_ncbi$V3 <- as.character(genes_abr_ncbi$V3)

#gene_order <- factor(van_genes$V6, level = gen_ord)

colors2 <- c("blue4", "red2", "red", "orangered", "sienna1", "darkorange", "blue3", "blue1", "dodgerblue2", "dodgerblue", "deepskyblue2", "deepskyblue", "red4", "red3")
ncbi_genes_heatmap <- ggplot(data = genes_abr_ncbi, mapping = aes(x = V6,
                                                            y = Sample_id,
                                                            fill = V3)) +
  geom_tile(colour="white", size=1) +
  labs(x = "", y = "") +
  scale_fill_manual(values = colors2) +
  theme(legend.position = 'none', axis.text.x = element_text(size = 8, angle = 90), axis.text.y = element_text(size = 6)) + 
  facet_grid(V4 ~ ., scales = "free_y", space = "free")  +
  ggtitle("Abricate - NCBI")

ncbi_genes_heatmap


#-------------------------------------------------------------------------------
# methicillin resistance genes 

# transform data for heatmap
heatmap_meth_genes_big <- pivot_longer(data = unique(meth), 
                           cols = -c(1),
                           names_to = "ToolDB", 
                           values_to = "Genes")

heatmap_meth_genes <- unique(heatmap_meth_genes_big)
heatmap_meth_genes <- heatmap_meth_genes[-c(which(is.na(heatmap_meth_genes$Genes))), ]
heatmap_meth_genes$Sample_id <- gsub('\\s+', '', heatmap_meth_genes$Sample_id)
heatmap_meth_genes$V4 <- 1

heatmap_meth_genes[heatmap_meth_genes == "(Bla)mecA"] <- "mecA"
heatmap_meth_genes[heatmap_meth_genes == "(Bla)mecC"] <- "mecC"
heatmap_meth_genes[heatmap_meth_genes == "MECA"] <- "mecA"
heatmap_meth_genes[heatmap_meth_genes == "MECI"] <- "mecI"
heatmap_meth_genes[heatmap_meth_genes == "MecC-type methicillin resistance repressor MecI"] <- "mecI_of_mecC"

for (i in 1:length(heatmap_meth_genes$Sample_id)) { 
  ifelse (heatmap_meth_genes$Genes[i] == 'mecA', heatmap_meth_genes$V4[i] <- length(which(heatmap_meth_genes$Sample_id == heatmap_meth_genes$Sample_id[i] & heatmap_meth_genes$Genes == 'mecA')),
          ifelse(heatmap_meth_genes$Genes[i] == 'mecC', heatmap_meth_genes$V4[i] <- length(which(heatmap_meth_genes$Sample_id == heatmap_meth_genes$Sample_id[i] & heatmap_meth_genes$Genes == 'mecC')),
                 ifelse(heatmap_meth_genes$Genes[i] == 'mecI', heatmap_meth_genes$V4[i] <- length(which(heatmap_meth_genes$Sample_id == heatmap_meth_genes$Sample_id[i] & heatmap_meth_genes$Genes == 'mecI')), 
                        ifelse(heatmap_meth_genes$Genes[i] == 'mecR1', heatmap_meth_genes$V4[i] <- length(which(heatmap_meth_genes$Sample_id == heatmap_meth_genes$Sample_id[i] & heatmap_meth_genes$Genes == 'mecR1')),
                               ifelse(heatmap_meth_genes$Genes[i] == 'mecI_of_mecA', heatmap_meth_genes$V4[i] <- length(which(heatmap_meth_genes$Sample_id == heatmap_meth_genes$Sample_id[i] & heatmap_meth_genes$Genes == 'mecI_of_mecA')), 
                                      ifelse(heatmap_meth_genes$Genes[i] == 'mecI_of_mecC', heatmap_meth_genes$V4[i] <- length(which(heatmap_meth_genes$Sample_id == heatmap_meth_genes$Sample_id[i] & heatmap_meth_genes$Genes == 'mecI_of_mecC')),
                                             heatmap_meth_genes$V4[i] <- length(which(heatmap_meth_genes$Sample_id == heatmap_meth_genes$Sample_id[i] & heatmap_meth_genes$Genes == 'PBP2'))))))))
  }

heatmap_meth_genes$V4 <- as.character(heatmap_meth_genes$V4)
colors3 <- c('green4', 'orange', 'yellow', 'blue', 'violet', 'cyan', 'black', 'red')

heatmap_meth_genes$V4 <- factor(heatmap_meth_genes$V4, level = c(14,11,7,5,4,3,2,1))

meth_genes_heatmap <- ggplot(data = heatmap_meth_genes[1:2245,], mapping = aes(x = Genes,
                                                                  y = Sample_id,
                                                                  fill = V4)) +
  geom_tile(colour="white", size=1) +
  labs(x = "", y = "") +
  labs(fill = 'in # of tools found') +
  scale_fill_manual(values = colors3) +
  coord_fixed(ratio = 0.5) +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust=1, vjust=0.5), axis.text.y = element_text(size = 4)) +  
  #facet_grid(. ~ ToolDB, scales = "free", space = "free")  +
  ggtitle("Meth genes")

meth_genes_heatmap

#-------------------------------------------------------------------------------
# plot phenotypes
phenotypes_filtered <- phenotypes
phenotypes_filtered$Dummy <- NULL
phenotypes_filtered[phenotypes_filtered == "NULL" | phenotypes_filtered == "U"] <- 0
phenotypes_filtered <- Filter(function(x) !(all(x==0)), phenotypes_filtered)

phenotypes_filtered[phenotypes_filtered == "c(R, R)"] <- "R"
phenotypes_filtered[phenotypes_filtered == "c(S, R)"] <- "R"
phenotypes_filtered[phenotypes_filtered == "c(R, S)"] <- "R"
phenotypes_filtered[phenotypes_filtered == "c(U, R)"] <- "R"
phenotypes_filtered[phenotypes_filtered == "c(R, U)"] <- "R"
phenotypes_filtered[phenotypes_filtered == "c(R, I)"] <- "R"
phenotypes_filtered[phenotypes_filtered == "c(I, R)"] <- "R"
phenotypes_filtered[phenotypes_filtered == "c(S, S)"] <- "S"
phenotypes_filtered[phenotypes_filtered == "c(U, S, S)"] <- "S"
phenotypes_filtered[phenotypes_filtered == "c(U, S)"] <- "S"
phenotypes_filtered[phenotypes_filtered == "c(S, U)"] <- "S"
phenotypes_filtered[phenotypes_filtered == "c(S, S, S)"] <- "S"
phenotypes_filtered[phenotypes_filtered == "c(I, I)"] <- "I"
phenotypes_filtered[phenotypes_filtered == "c(S, I)"] <- "I"
phenotypes_filtered[phenotypes_filtered == "c(I, S)"] <- "I"
phenotypes_filtered[phenotypes_filtered == "c(U, U)"] <- 0

heatmap_table_phenotypes <- pivot_longer(data = phenotypes_filtered, 
                              cols = -c(1:2),
                              names_to = "Antibiotic", 
                              values_to = "Phenotype")

heatmap_table_phenotypes <- heatmap_table_phenotypes[-c(which(heatmap_table_phenotypes$Phenotype == 0)), ]
  
heatmap_table_phenotypes$Phenotype <- factor(heatmap_table_phenotypes$Phenotype, levels = c("S", "I", "R"))
 

# create heatmap
colors <- c("darkgreen", "orange", "red")

heatmap_phenotypes <- ggplot(data = heatmap_table_phenotypes, mapping = aes(x = Antibiotic,
                                                         y = Sample_id,
                                                         fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") + 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust=1), axis.text.y = element_text(size = 6)) +
  #coord_fixed(ratio = 0.2) +
  #theme(axis.text.y = element_text(size = 6)) +
  facet_grid(LGM_BK ~ ., scales = "free_y", space = "free")

heatmap_phenotypes

# heatmap for selected antibiotics
heatmap_table_phenotypes_filtered <- pivot_longer(data = phenotypes_filtered[,c(1,2,7,8,12,15,16,18,26)], 
                                         cols = -c(1:2),
                                         names_to = "Antibiotic", 
                                         values_to = "Phenotype")

heatmap_table_phenotypes_filtered <- heatmap_table_phenotypes_filtered[-c(which(heatmap_table_phenotypes_filtered$Phenotype == 0)), ]
heatmap_table_phenotypes_filtered$Phenotype <- factor(heatmap_table_phenotypes_filtered$Phenotype, levels = c("S", "I", "R"))

heatmap_phenotypes_filtered <- ggplot(data = heatmap_table_phenotypes_filtered, mapping = aes(x = Antibiotic,
                                                                            y = Sample_id,
                                                                            fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") + 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust=1), axis.text.y = element_text(size = 6)) +
  #coord_fixed(ratio = 0.2) +
  #theme(axis.text.y = element_text(size = 6)) +
  facet_grid(LGM_BK ~ ., scales = "free_y", space = "free")

heatmap_phenotypes_filtered
