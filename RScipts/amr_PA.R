#install.packages('tidyverse')
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)

setwd("~/data_project/assemblies/PA")

# read in samples to be analyzed
samples <- read.table('samples_PA.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_PA.csv', header = TRUE, sep = ';')
mic <- read.csv2('mic_PA.csv', header = TRUE, sep = ';')
pheno_mic <- read.csv2('pheno_mic_PA.csv', header = TRUE, sep = ';')

# check sample ids; and correct them 
which(phenotypes$Sample_id %in% samples$V1 == FALSE)
phenotypes$Sample_id[5] <- "503076-2-16"
phenotypes$Sample_id[18] <- "703687-17-wh"
phenotypes$Sample_id[22] <- "804523-10-18"
phenotypes$Sample_id[23] <- "712331-3-18"
phenotypes$Sample_id[26] <- "807875-6-18"
phenotypes$Sample_id[53] <- "716436-19-wh"
phenotypes$Sample_id[58] <- "613038-19-wh"
phenotypes$Sample_id[68] <- "713135-2-20"
phenotypes$Sample_id[74] <- "504345-2-20"
phenotypes$Sample_id[84] <- "503825-21-wh2"

#which(mic$Sample_id %in% samples$V1 == FALSE)
mic$Sample_id[5] <- "503076-2-16"
mic$Sample_id[18] <- "703687-17-wh"
mic$Sample_id[22] <- "804523-10-18"
mic$Sample_id[23] <- "712331-3-18"
mic$Sample_id[26] <- "807875-6-18"
mic$Sample_id[53] <- "716436-19-wh"
mic$Sample_id[58] <- "613038-19-wh"
mic$Sample_id[68] <- "713135-2-20"
mic$Sample_id[74] <- "504345-2-20"
mic$Sample_id[84] <- "503825-21-wh2"

which(pheno_mic$Sample_id %in% samples$V1 == FALSE)
pheno_mic$Sample_id[49] <- "503076-2-16"
pheno_mic$Sample_id[56] <- "503825-21-wh2"
pheno_mic$Sample_id[59] <- "504345-2-20"
pheno_mic$Sample_id[72] <- "613038-19-wh"
pheno_mic$Sample_id[76] <- "703687-17-wh"
pheno_mic$Sample_id[77] <- "712331-3-18"
pheno_mic$Sample_id[78] <- "713135-2-20"
pheno_mic$Sample_id[81] <- "716436-19-wh"
pheno_mic$Sample_id[84] <- "804523-10-18"
pheno_mic$Sample_id[85] <- "807875-6-18"


pheno_mic[pheno_mic == "c(= 0, > 32)" | pheno_mic == "c(> 32, > 32)" | pheno_mic == "c(> 32, >= 16)"] <- "> 32"
pheno_mic[pheno_mic == "c(>= 16, >= 16)" | pheno_mic == "c(>= 16, = 16)" | pheno_mic == "c(> 32, >= 16, > 32, >= 16)"] <- ">= 16"
pheno_mic[pheno_mic == "c(= 16, = 8)"] <- "= 8"
pheno_mic[pheno_mic == "c(= 1, = 1)" | pheno_mic == "c(= 1, = 1, = 1)"] <- "= 1"
pheno_mic[pheno_mic == "c(<= 0.25, <= 0.25, <= 0.25)" | pheno_mic == "c(= 0.125, <= 0.25)" | pheno_mic == "c(<= 0.25, <= 0.25)"] <- "<= 0.25"
pheno_mic[pheno_mic == "c(= 0.5, n.a. 0)" | pheno_mic == "= 0.5"] <- "<= 0.5"
pheno_mic[pheno_mic == "n.a. 10016"] <- "NA"

for (i in 1:length(pheno_mic$Sample_id)){
  if(pheno_mic$MIC.Meropenem[i] == "NULL") {
    pheno_mic$MIC.Meropenem[i] <- pheno_mic$MIC.Meropenem.ohne.Meningitis[i]
  }
}


# combine antibiotics with several columns
pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "S"] <- "S"
pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "R"] <- "R"
pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "I"] <- "I"


pheno_mic[pheno_mic == "c(R, R)" | pheno_mic == "c(S, R)" | pheno_mic == "c(R, S)" | pheno_mic == "c(S, U, R)" |   
                      pheno_mic == "c(R, S, S)" | pheno_mic == "c(R, R, R)" | pheno_mic == "c(R, U)" | 
                      pheno_mic == "c(U, R)" | pheno_mic == "c(U, R, U, R)" | pheno_mic == "c(R, R, I)" |
                      pheno_mic == "c(R, U, U, I, I)" | pheno_mic == "c(I, U, U, R)" | pheno_mic == "c(R, U, U, I)" |
                      pheno_mic == "c(R, I, R, I)" | pheno_mic == "c(S, S, R)" | pheno_mic == "c(R, I, R, R)" |
                      pheno_mic == "c(U, R, U)" | pheno_mic == "c(U, R, R)" | pheno_mic == "c(U, R, U, U, R)" |
                      pheno_mic == "c(U, U, R, S)" | pheno_mic == "c(U, R, U, I)" | pheno_mic == "c(U, U, R)" |
                      pheno_mic == "c(U, R, R, U, R, R)" |pheno_mic == "c(R, U, R)" | pheno_mic == "c(R, R, R, R)" |
                      pheno_mic == "c(R, U, S)"] <- "R"

pheno_mic[pheno_mic == "c(S, S)" | pheno_mic == "c(U, S)" | pheno_mic == "c(S, S, S)" | pheno_mic == "c(S, U)" |
                      pheno_mic == "c(S, U, S)" | pheno_mic == "c(S, S, U)" | pheno_mic == "c(S, U, U, S)" |
                      pheno_mic == "c(U, S, S)"] <- "S"

pheno_mic[pheno_mic == "c(I, S)" | pheno_mic == "c(S, I)" | pheno_mic == "c(I, S, S)" |
                      pheno_mic == "c(S, I, S, I)" | pheno_mic == "c(I, I)" | pheno_mic == "c(I, U)" | 
                      pheno_mic == "c(U, I)" | pheno_mic == "c(U, I, U, U)" | pheno_mic == "c(U, I, U, S)" |
                      pheno_mic == "c(U, S, U, I)" | pheno_mic == "c(U, I, U)" | pheno_mic == "c(I, I, U)" |
                      pheno_mic == "c(I, U, U)" | pheno_mic == "c(U, S, I)" | pheno_mic == "c(R, I)" |
                      pheno_mic == "c(I, S, I)"] <- "I"



#-------------------------------------------------------------------------------
# get result tables
# abricate
abr_ncbi <- read.delim('PA_abricate_ncbi.tab', sep = "\t", header = T)
if (length(which(abr_ncbi$X.IDENTITY < 95)) != 0) {
  abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$X.IDENTITY < 95)),]
}
abr_ncbi$X.FILE <- gsub(".fna","", abr_ncbi$X.FILE)
names(abr_ncbi)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_ncbi$Tool <- 'abr_ncbi'

abr_card <- read.delim('PA_abricate_card.tab', sep = "\t", header = T)
if (length(which(abr_card$X.IDENTITY < 95)) != 0) {
  abr_card <- abr_card[-c(which(abr_card$X.IDENTITY < 95)),]
}
abr_card$X.FILE <- gsub(".fna","", abr_card$X.FILE)
names(abr_card)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_card$Tool <- 'abr_card'

abr_megares <- read.delim('PA_abricate_megares.tab', sep = "\t", header = T)
if (length(which(abr_megares$X.IDENTITY < 95)) != 0) {
  abr_megares <- abr_megares[-c(which(abr_megares$X.IDENTITY < 95)),]
}
abr_megares$X.FILE <- gsub(".fna","", abr_megares$X.FILE)
names(abr_megares)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_megares$Tool <- 'abr_megares'

abr_resfinder <- read.delim('PA_abricate_resfinder.tab', sep = "\t", header = T)
if (length(which(abr_resfinder$X.IDENTITY < 95)) != 0) {
  abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$X.IDENTITY < 95)),]
}
abr_resfinder$X.FILE <- gsub(".fna","", abr_resfinder$X.FILE)
names(abr_resfinder)[c(1,14)] <- c('Sample_id', 'amr_genes')
abr_resfinder$Tool <- 'abr_resfinder'

abr_argannot <- read.delim('PA_abricate_argannot.tab', sep = "\t", header = T)
if (length(which(abr_argannot$X.IDENTITY < 95)) != 0) {
  abr_argannot <- abr_argannot[-c(which(abr_argannot$X.IDENTITY < 95)),]
}
abr_argannot$X.FILE <- gsub(".fna","", abr_argannot$X.FILE)
names(abr_argannot)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_argannot$Tool <- 'abr_argannot'


#amrfinder
amrf_nuc <- read.delim('PA_amrfinder_nuc.tab', sep = "\t")
if (length(which(amrf_nuc$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_nuc)[1] <- 'Sample_id'
amrf_nuc$Sample_id <- gsub("gnl\\|USB\\|","", amrf_nuc$Contig.id)
amrf_nuc$Sample_id <- gsub('(_\\d*)',"", amrf_nuc$Sample_id)
names(amrf_nuc)[6] <- 'amr_genes'
amrf_nuc$Tool <- 'amrf_nuc'

amrf_prot <- read.delim('PA_amrfinder_prot.tab', sep = "\t")
if (length(which(amrf_prot$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_prot <- amrf_prot[-c(which(amrf_prot$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_prot)[c(1,2)] <- c('Sample_id', 'amr_genes')
amrf_prot$Sample_id <- gsub('(_\\d*)',"", amrf_prot$Sample_id)
amrf_prot$Tool <- 'amrf_prot'


# rgi
rgi <- read.delim(('PA_rgi.tab'), sep = "\t")
if (length(which(rgi$Best_Identities < 95)) != 0) {
  rgi <- rgi[-c(which(rgi$Best_Identities < 95)),]
}
names(rgi)[c(1,2)] <- c('Sample_id', 'amr_genes')
rgi$Sample_id <- gsub("gnl\\|USB\\|","", rgi$Sample_id)
rgi$Sample_id <- gsub('(_\\d*)',"", rgi$Sample_id)
rgi$Sample_id <- gsub('\\s+', '', rgi$Sample_id)
rgi$Tool <- 'rgi'

# sraX
srax_basic <- read.delim(('PA_srax_basic.tab'), sep = "\t")
if (length(which(srax_basic$Identity_p < 95)) != 0) {
  srax_basic <- srax_basic[-c(which(srax_basic$Identity_p < 95)),]
}
names(srax_basic)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_basic$Sample_id <- gsub(".fna","", srax_basic$Sample_id)
srax_basic$Tool <- 'sraX_basic'

srax_ext <- read.delim(('PA_srax_ext.tab'), sep = "\t")
if (length(which(srax_ext$Identity_p < 95)) != 0) {
  srax_ext <- srax_ext[-c(which(srax_ext$Identity_p < 95)),]
}
names(srax_ext)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_ext$Sample_id <- gsub(".fna","", srax_ext$Sample_id)
srax_ext$Tool <- 'sraX_ext'


# Deeparg
deeparg_LS <- read.delim(('PA_deeparg_LS_noPot.tab'), sep = "\t")
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

deeparg_SR <- read.delim(('PA_deeparg_SR_noPot.tab'), sep = "\t")
if (length(which(deeparg_SR$identity < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$identity < 95)),]
}
if (length(which(deeparg_SR$probability < 0.8)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$probability < 0.8)),]
}
deeparg_SR_unique <- deeparg_SR %>% distinct(X.ARG, best.hit, .keep_all = TRUE)
names(deeparg_SR_unique)[1] <- 'Sample_id'
deeparg_SR_unique$Sample_id <- gsub("pseaer/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)

for(i in 1:length(deeparg_SR_unique[,1])) {
  genes <- unlist(strsplit(deeparg_SR_unique[i,1], "\\_"))
  deeparg_SR_unique[i,1] <- genes[1]
}
for(i in 1:length(deeparg_SR_unique[,1])) {
  genes <- unlist(strsplit(deeparg_SR_unique[i,6], "\\|"))
  deeparg_SR_unique[i,6] <- genes[3]
}

names(deeparg_SR_unique)[6] <- 'amr_genes'
deeparg_SR_unique$Tool <- 'deeparg_SR'

rm(deeparg_SR)

# resfinder
resfinder_as <- read.delim(('PA_resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$Identity < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$Identity < 95)),]
}
names(resfinder_as)[c(6,1)] <- c('Sample_id', 'amr_genes')
resfinder_as$Sample_id <- gsub('gnl\\|USB\\|','', resfinder_as$Sample_id)
resfinder_as$Sample_id <- gsub('(_\\d*)','', resfinder_as$Sample_id)
resfinder_as$Tool <- 'resfinder_as'

resfinder_re <- read.delim(('PA_resfinder_reads.tab'), sep = "\t", header = F)
names(resfinder_re) <- c('Resistance gene',	'Identity',	'Alignment Length/Gene Length',	'Coverage',	'Position in reference',	'Contig',	'Position in contig',  'Phenotype',	'Accession no.', 'Sample_id')
resfinder_re <- resfinder_re[-c(1),]
resfinder_re$Identity <- as.numeric(resfinder_re$Identity )
if (length(which(resfinder_re$Identity < 95)) != 0) {
  resfinder_re <- resfinder_re[-c(which(resfinder_re$Identity < 95)),]
}
names(resfinder_re)[1] <- 'amr_genes'
resfinder_re$Tool <- 'resfinder_re'


#-------------------------------------------------------------------------------
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

srax_basic_bla <- srax_basic[c(grep("beta-lactamase|OprD", srax_basic$Gene_description, ignore.case = TRUE)), ]
names(srax_basic_bla)[6] <- 'bla_genes'
srax_basic_bla$Tool <- 'srax_basic'
srax_ext_bla <- srax_ext[c(grep("beta-lactamase|OprD", srax_ext$Gene_description, ignore.case = TRUE)), ]
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


#-------------------------------------------------------------------------------
# gene cluster

bla_gene_clusters <- bla
bla_gene_clusters$bla_genes <- gsub("\\(Bla\\)", "", bla_gene_clusters$bla_genes)
bla_gene_clusters$bla_genes <- gsub("bla", "", bla_gene_clusters$bla_genes)


gene_clusters <- pivot_wider(bla_gene_clusters, names_from = Tool, values_from = bla_genes)

gene_cluster_table <- pivot_longer(data = gene_clusters, 
                                   cols = -c(1:4),
                                   names_to = "ToolDB", 
                                   values_to = "gene_cluster")

gene_cluster_table <- as.data.frame(gene_cluster_table)
gene_cluster_table[gene_cluster_table == "NULL"] <- "none"

for (i in 1:length(gene_cluster_table$gene_cluster)) {
  gene_cluster_table$gene_cluster[[i]] <- sort(gene_cluster_table$gene_cluster[[i]])
}


gene_cluster_table$cluster_rep <- ifelse(grepl("\\bVIM-1\\b",gene_cluster_table$gene_cluster) == "TRUE", 1, 
                                         ifelse(grepl("\\bVIM-2\\b",gene_cluster_table$gene_cluster) == "TRUE", 2,
                                                ifelse(grepl("\\bVIM\\b",gene_cluster_table$gene_cluster) == "TRUE", 3,
                                                       ifelse(grepl("OXA-573",gene_cluster_table$gene_cluster) == "TRUE", 4, 
                                                              ifelse(grepl("metallo-beta-lactamase", gene_cluster_table$gene_cluster) == "TRUE", 5,
                                                                     ifelse(grepl("OXA-418",gene_cluster_table$gene_cluster) == "TRUE", 5, 
                                                                            ifelse(grepl("OXA-202",gene_cluster_table$gene_cluster) == "TRUE", 5, 
                                                                                   ifelse(grepl("oprD",gene_cluster_table$gene_cluster) == "TRUE", 7,6))))))))

for (i in 1:length(gene_cluster_table$Sample_id)) {
  if (gene_cluster_table$ToolDB[i] == "deeparg_SR" & gene_cluster_table$Type[i] == "pseaerc") {
    gene_cluster_table$cluster_rep[i] <- 5
  }
}

part_oprd_srax_b <- srax_basic$Sample_id[which(srax_basic$amr_genes == "oprD" & srax_basic$Status_hit == "Partial gene, no gaps")]
part_oprd_srax_e <- srax_ext$Sample_id[which(srax_ext$amr_genes == "oprD" & srax_ext$Status_hit == "Partial gene, no gaps")]

for (i in 1:length(gene_cluster_table$Sample_id)) {
  if ((gene_cluster_table$Sample_id[i] %in% part_oprd_srax_b & gene_cluster_table$ToolDB[i] == "srax_basic") | 
      (gene_cluster_table$Sample_id[i] %in% part_oprd_srax_e & gene_cluster_table$ToolDB[i] == "srax_ext")) {
    gene_cluster_table$cluster_rep[i] <- 8
  }
}  


gene_cluster_table$cluster_rep <- as.character(gene_cluster_table$cluster_rep)

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


# order legend
gene_cluster_table$cluster_rep <- factor(gene_cluster_table$cluster_rep, levels = c(1,2,3,5,4,7,8,6)) 
colors <- c("blue4", "blue", "lightsteelblue", "yellowgreen", "orangered", "magenta", "mediumpurple1", "moccasin")

a <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "1"]))
b <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "2"]))
c <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & (gene_cluster_table$cluster_rep == "6")]))

gene_cluster_table$Sample_id <- factor(gene_cluster_table$Sample_id, levels = c(a, b, c)) 


gene_cluster_heatmap <- ggplot(data = gene_cluster_table, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("VIM-1", "VIM-2", "VIM", "mixed", "OXA-573", "oprD", "partial oprD","none")) +
  labs(x = "", y = "") +
  labs(fill = "Carbapenemase genes") +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, angle = 90, hjust=1, vjust=0.5), 
        strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6), panel.background = element_rect(color = "white", linetype = "solid")) +
  facet_grid(vars(factor(gene_cluster_table$Type, levels = c("pseaerc", "pseaer", "pseaermu"))), vars(Tool), scales = "free", space = "free")

gene_cluster_heatmap

gene_cluster_heatmap <- ggplot(data = gene_cluster_table, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("VIM-1", "VIM-2", "VIM", "mixed", "OXA-573", "oprD", "partial oprD", "none")) +
  labs(x = "", y = "") +
  labs(fill = "Carbapenemase genes") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(gene_cluster_table$MeropenemMIC, levels = c("<= 0.25", "","<= 0.5", "= 0.75","= 1", "= 2",
                                                                      "= 3", "= 4", "= 6", "= 8", ">= 16", "> 32", "NA"))), vars(Tool), scales = "free", space = "free")

gene_cluster_heatmap


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


carb_correlation <- pivot_wider(bla_for_heatmap[,c(1,3:7)], names_from = Tool, values_from = Classification)

concordance_table <- carb_correlation

# Concordance values: R/R = 1, S/S = 2, S/R = 3, R/S = 4
concordance_table$abr_ncbi <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 1,
                                     ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 3,
                                            ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$abr_ncbi) == "FALSE", 2, 4)))

concordance_table$abr_card <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 1,
                                     ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 3,
                                            ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$abr_card) == "FALSE", 2, 4)))

concordance_table$abr_megares <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 1,
                                        ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 3,
                                               ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$abr_megares) == "FALSE", 2, 4)))

concordance_table$abr_resfinder <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 1,
                                          ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 3,
                                                 ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$abr_resfinder) == "FALSE", 2, 4)))

concordance_table$abr_argannot <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 1,
                                         ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 3,
                                                ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$abr_argannot) == "FALSE", 2, 4)))

concordance_table$amrf_nuc <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 1,
                                     ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 3,
                                            ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$amrf_nuc) == "FALSE", 2, 4)))

concordance_table$amrf_prot <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 1,
                                      ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 3,
                                             ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$amrf_prot) == "FALSE", 2, 4)))

concordance_table$deeparg_LS <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 1,
                                       ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 3,
                                              ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$deeparg_LS) == "FALSE", 2, 4)))

concordance_table$deeparg_SR <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 1,
                                       ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 3,
                                              ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$deeparg_SR) == "FALSE", 2, 4)))

concordance_table$resfinder_as <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 1,
                                         ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 3,
                                                ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$resfinder_as) == "FALSE", 2, 4)))

concordance_table$resfinder_re <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 1,
                                         ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 3,
                                                ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$resfinder_re) == "FALSE", 2, 4)))

concordance_table$rgi <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 1,
                                ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 3,
                                       ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$rgi) == "FALSE", 2, 4)))

concordance_table$srax_basic <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 1,
                                       ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 3,
                                              ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$srax_basic) == "FALSE", 2, 4)))


concordance_table$srax_ext <- ifelse(concordance_table$Type == "pseaerc" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 1,
                                     ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu")  & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 3,
                                            ifelse((concordance_table$Type == "pseaer" | concordance_table$Type == "pseaermu") & grepl("Carbapenem", concordance_table$srax_ext) == "FALSE", 2, 4)))


# transform data for heatmap
heatmap_concordance <- pivot_longer(data = concordance_table, 
                                    cols = -c(1:4),
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

# do it after step before otherwise card of rgi is counted as abricate
#heatmap_concordance[heatmap_concordance == "rgi"] <- "card"

# create heatmap
colors <- c("darkgreen", "green", "gold", "red3")

heatmap_concordance$Sample_id <- factor(heatmap_concordance$Sample_id, levels = c(a, b, c)) 

concordance_heatmap <- ggplot(data = heatmap_concordance, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "S/R", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance$MeropenemMIC, levels = c("<= 0.25", "","<= 0.5", "= 0.75","= 1", "= 2",
                                                                       "= 3", "= 4", "= 6", "= 8", ">= 16", "> 32", "NA"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap


#-------------------------------------------------------------------------------
# plot mero imi phenotypes

mero_imi <- pheno_mic[,c(1,2,18,32,103,117)]

mero_imi_table <- pivot_longer(data = mero_imi, 
                               cols = -c(1:2,5:6),
                               names_to = "Antibiotic", 
                               values_to = "Phenotype")

mero_imi_table$Sample_id <- factor(mero_imi_table$Sample_id, levels = c(a,b,c))

colors <- c("gray", "lemonchiffon", "black", "coral")
mero_imi_table$Phenotype <- factor(mero_imi_table$Phenotype, levels = c("S", "I", "R")) 

mero_imi_heatmap <- ggplot(data = mero_imi_table, mapping = aes(x = Antibiotic,
                                                                y = Sample_id,
                                                                fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(LGM_BK), scales = "free", space = "free")

mero_imi_heatmap

mero_imi_heatmap <- ggplot(data = mero_imi_table, mapping = aes(x = Antibiotic,
                                                          y = Sample_id,
                                                          fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi_table$MIC.Meropenem, levels = c("<= 0.25", "","<= 0.5", "= 0.75","= 1", "= 2",
                                                            "= 3", "= 4", "= 6", "= 8", ">= 16", "> 32", "NA"))), scales = "free", space = "free")

mero_imi_heatmap


mero_imi$Sample_id <- factor(mero_imi$Sample_id, levels = c(a,b,c))

mero_mic <- ggplot(data = mero_imi, mapping = aes(x = "Meropenem MIC",
                                                  y = Sample_id,
                                                  fill = MIC.Meropenem)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi$MIC.Meropenem, levels = c("<= 0.25", "","<= 0.5", "= 0.75","= 1", "= 2",
                                                            "= 3", "= 4", "= 6", "= 8", ">= 16", "> 32", "NA"))), scales = "free", space = "free")

mero_mic


#-------------------------------------------------------------------------------
# plot porin

porin_data <- read.csv2('porin.csv', header = FALSE, sep = ',')

porin <- heatmap_concordance

for (i in 1:length(porin$Sample_id)){
  porin$oprD[i] <- porin_data$V2[which(porin$Sample_id[i] == porin_data$V1)]
}

porin$Sample_id <- factor(porin$Sample_id, levels = c(a,b,c))
colors <- c("darkturquoise", "darkmagenta")

porin_heatmap <- ggplot(data = porin, mapping = aes(x = "oprD",
                                                          y = Sample_id,
                                                          fill = oprD)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("intact", "mutation")) +
  labs(x = "", y = "") +
  labs(fill = "oprD") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(porin$MeropenemMIC, levels = c("<= 0.25", "","<= 0.5", "= 0.75","= 1", "= 2",
                                                                      "= 3", "= 4", "= 6", "= 8", ">= 16", "> 32", "NA"))), scales = "free", space = "free")

porin_heatmap

