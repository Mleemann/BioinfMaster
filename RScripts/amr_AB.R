#install.packages('tidyverse')
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)

setwd("~/data_project/assemblies/AB")

# read in samples to be analyzed
samples <- read.table('samples_AB.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_AB.csv', header = TRUE, sep = ';')
mic <- read.csv2('mic_AB.csv', header = TRUE, sep = ';')
pheno_mic <- read.csv2('pheno_mic_AB.csv', header = TRUE, sep = ';')

for (i in 1:length(pheno_mic$Sample_id)) {
  if (pheno_mic$MIC.Meropenem.bei.Meningitis[i] != "NULL") {
    pheno_mic$MIC.Meropenem[i] <- pheno_mic$MIC.Meropenem.bei.Meningitis[i]
  }
}

pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "S"] <- "S"
pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "R"] <- "R"
pheno_mic$Meropenem[pheno_mic$Meropenem.ohne.Meningitis  == "I"] <- "I"


#unique(as.vector(as.matrix(pheno_mic)))
pheno_mic[pheno_mic == "c(R, R)" | pheno_mic == "c(S, R)" | pheno_mic == "c(R, S)" | pheno_mic == "c(R, I)" |   
                      pheno_mic == "c(R, U)" | pheno_mic == "c(U, R)"] <- "R"
pheno_mic[pheno_mic == "c(S, S)" | pheno_mic == "c(S, U)" | pheno_mic == "c(S, S, U)"] <- "S"
pheno_mic[pheno_mic == "c(U, U)" | pheno_mic == "c(U, U, U)" | pheno_mic == "K"] <- 0

pheno_mic[pheno_mic == "c(>= 16, >= 16)" | pheno_mic == "c(> 32, >= 16)" | pheno_mic == "c(>= 16, > 32)" | pheno_mic == "c(>= 16, = 0)"] <- ">= 16"
pheno_mic[pheno_mic == "c(> 32, = 8)"] <- "= 8"
pheno_mic[pheno_mic == "c(= 16, = 2)"] <- "= 2"




# check sample ids; and correct them 
which(phenotypes$Sample_id %in% samples$V1 == FALSE)
phenotypes$Sample_id[5] <- "808937-18-16"
phenotypes$Sample_id[14] <- "806847-8-19"
phenotypes$Sample_id[21] <- "805022-16-20"
phenotypes$Sample_id[22] <- "805143-3-20"
phenotypes$Sample_id[26] <- "503504-21-wh"

#which(mic$Sample_id %in% samples$V1 == FALSE)
mic$Sample_id[5] <- "808937-18-16"
mic$Sample_id[14] <- "806847-8-19"
mic$Sample_id[21] <- "805022-16-20"
mic$Sample_id[22] <- "805143-3-20"
mic$Sample_id[26] <- "503504-21-wh"

which(pheno_mic$Sample_id %in% samples$V1 == FALSE)
pheno_mic$Sample_id[30] <- "808937-18-16"
pheno_mic$Sample_id[29] <- "806847-8-19"
pheno_mic$Sample_id[23] <- "805022-16-20"
pheno_mic$Sample_id[24] <- "805143-3-20"
pheno_mic$Sample_id[10] <- "503504-21-wh"


#-------------------------------------------------------------------------------
# get result tables
# abricate
abr_ncbi <- read.delim('AB_abricate_ncbi.tab', sep = "\t", header = T)
if (length(which(abr_ncbi$X.IDENTITY < 95)) != 0) {
  abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$X.IDENTITY < 95)),]
}
abr_ncbi$X.FILE <- gsub(".fna","", abr_ncbi$X.FILE)
names(abr_ncbi)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_ncbi$Tool <- 'abr_ncbi'

abr_card <- read.delim('AB_abricate_card.tab', sep = "\t", header = T)
if (length(which(abr_card$X.IDENTITY < 95)) != 0) {
  abr_card <- abr_card[-c(which(abr_card$X.IDENTITY < 95)),]
}
abr_card$X.FILE <- gsub(".fna","", abr_card$X.FILE)
names(abr_card)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_card$Tool <- 'abr_card'

abr_megares <- read.delim('AB_abricate_megares.tab', sep = "\t", header = T)
if (length(which(abr_megares$X.IDENTITY < 95)) != 0) {
  abr_megares <- abr_megares[-c(which(abr_megares$X.IDENTITY < 95)),]
}
abr_megares$X.FILE <- gsub(".fna","", abr_megares$X.FILE)
names(abr_megares)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_megares$Tool <- 'abr_megares'

abr_resfinder <- read.delim('AB_abricate_resfinder.tab', sep = "\t", header = T)
if (length(which(abr_resfinder$X.IDENTITY < 95)) != 0) {
  abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$X.IDENTITY < 95)),]
}
abr_resfinder$X.FILE <- gsub(".fna","", abr_resfinder$X.FILE)
names(abr_resfinder)[c(1,14)] <- c('Sample_id', 'amr_genes')
abr_resfinder$Tool <- 'abr_resfinder'

abr_argannot <- read.delim('AB_abricate_argannot.tab', sep = "\t", header = T)
if (length(which(abr_argannot$X.IDENTITY < 95)) != 0) {
  abr_argannot <- abr_argannot[-c(which(abr_argannot$X.IDENTITY < 95)),]
}
abr_argannot$X.FILE <- gsub(".fna","", abr_argannot$X.FILE)
names(abr_argannot)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_argannot$Tool <- 'abr_argannot'


#amrfinder
amrf_nuc <- read.delim('AB_amrfinder_nuc.tab', sep = "\t")
if (length(which(amrf_nuc$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_nuc)[1] <- 'Sample_id'
amrf_nuc$Sample_id <- gsub("gnl\\|USB\\|","", amrf_nuc$Contig.id)
amrf_nuc$Sample_id <- gsub('(_\\d*)',"", amrf_nuc$Sample_id)
names(amrf_nuc)[6] <- 'amr_genes'
amrf_nuc$Tool <- 'amrf_nuc'

amrf_prot <- read.delim('AB_amrfinder_prot.tab', sep = "\t")
if (length(which(amrf_prot$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_prot <- amrf_prot[-c(which(amrf_prot$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_prot)[c(1,2)] <- c('Sample_id', 'amr_genes')
amrf_prot$Sample_id <- gsub('(_\\d*)',"", amrf_prot$Sample_id)
amrf_prot$Tool <- 'amrf_prot'


# rgi
rgi <- read.delim(('AB_rgi.tab'), sep = "\t")
if (length(which(rgi$Best_Identities < 95)) != 0) {
  rgi <- rgi[-c(which(rgi$Best_Identities < 95)),]
}
names(rgi)[c(1,2)] <- c('Sample_id', 'amr_genes')
rgi$Sample_id <- gsub("gnl\\|USB\\|","", rgi$Sample_id)
rgi$Sample_id <- gsub('(_\\d*)',"", rgi$Sample_id)
rgi$Sample_id <- gsub('\\s+', '', rgi$Sample_id)
rgi$Tool <- 'rgi'

# sraX
srax_basic <- read.delim(('AB_srax_basic.tab'), sep = "\t")
if (length(which(srax_basic$Identity_p < 95)) != 0) {
  srax_basic <- srax_basic[-c(which(srax_basic$Identity_p < 95)),]
}
names(srax_basic)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_basic$Sample_id <- gsub(".fna","", srax_basic$Sample_id)
srax_basic$Tool <- 'sraX_basic'

srax_ext <- read.delim(('AB_srax_ext.tab'), sep = "\t")
if (length(which(srax_ext$Identity_p < 95)) != 0) {
  srax_ext <- srax_ext[-c(which(srax_ext$Identity_p < 95)),]
}
names(srax_ext)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_ext$Sample_id <- gsub(".fna","", srax_ext$Sample_id)
srax_ext$Tool <- 'sraX_ext'


# Deeparg
deeparg_LS <- read.delim(('AB_deeparg_LS_noPot.tab'), sep = "\t")
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

deeparg_SR <- read.delim(('AB_deeparg_SR_noPot.tab'), sep = "\t")
if (length(which(deeparg_SR$identity < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$identity < 95)),]
}
if (length(which(deeparg_SR$probability < 0.8)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$probability < 0.8)),]
}
deeparg_SR_unique <- deeparg_SR %>% distinct(X.ARG, best.hit, .keep_all = TRUE)
names(deeparg_SR_unique)[1] <- 'Sample_id'
deeparg_SR_unique$Sample_id <- gsub("acibau/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)
deeparg_SR_unique$Sample_id <- gsub("acibauc/resistance/deeparg_SR/","", deeparg_SR_unique$Sample_id)

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
resfinder_as <- read.delim(('AB_resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$Identity < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$Identity < 95)),]
}
names(resfinder_as)[c(6,1)] <- c('Sample_id', 'amr_genes')
resfinder_as$Sample_id <- gsub('gnl\\|USB\\|','', resfinder_as$Sample_id)
resfinder_as$Sample_id <- gsub('(_\\d*)','', resfinder_as$Sample_id)
resfinder_as$Tool <- 'resfinder_as'

resfinder_re <- read.delim(('AB_resfinder_reads.tab'), sep = "\t", header = F)
resfinder_re <- resfinder_re[1:10]
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
bla[bla == "802700-15_R" | bla == "802700-15R"] <- "802700-15"
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


# 1: OXA-23 & OXA-66  2: OXA-23 & OXA-58  3: OXA-23 & OXA-66 & OXA-58   4: OXA-23 & OXA-68
# 5: OXA-24           6: OXA-72 & OXA-66  7: NDM-1 & OXA-72 & OXA-66    8: OXA-90
# 9: OXA-72 & OXA-92  10: OXA-66

oxa51 <- "OXA-51|OXA-66|OXA-68|OXA-69|OXA-70|OXA-90|OXA-91|OXA-92|OXA-94|OXA-120|OXA-317|OXA-430|OXA-562|OXA-685|OXA-714|OXA-853"

#gene_cluster_table$cluster_rep <- ifelse(grepl("OXA-23",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-58",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-66",gene_cluster_table$gene_cluster) == "TRUE", 3, 
 #                                        ifelse(grepl("NDM-1",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-72",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-66",gene_cluster_table$gene_cluster) == "TRUE", 7,
  #                                              ifelse(grepl("OXA-23",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-66",gene_cluster_table$gene_cluster) == "TRUE", 1,
   #                                                    ifelse(grepl("OXA-23",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-58",gene_cluster_table$gene_cluster) == "TRUE", 2,
    #                                                          ifelse(grepl("OXA-23",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-68",gene_cluster_table$gene_cluster) == "TRUE", 4,
     #                                                                ifelse(grepl("OXA-23",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-91",gene_cluster_table$gene_cluster) == "TRUE", 12,      
      #                                                                      ifelse(grepl("OXA-24",gene_cluster_table$gene_cluster) == "TRUE", 5,
       #                                                                            ifelse(grepl("OXA-72",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-66",gene_cluster_table$gene_cluster) == "TRUE", 6, 
        #                                                                                  ifelse(grepl("OXA-90",gene_cluster_table$gene_cluster) == "TRUE", 8,
         #                                                                                        ifelse(grepl("OXA-72",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-92",gene_cluster_table$gene_cluster) == "TRUE", 9,
          #                                                                                              ifelse(grepl("OXA-66",gene_cluster_table$gene_cluster) == "TRUE", gene_cluster_table$cluster_rep <- 10,
           #                                                                                                    ifelse(grepl("OXA-51",gene_cluster_table$gene_cluster) == "TRUE", gene_cluster_table$cluster_rep <- 13,
            #                                                                                                          ifelse(grepl("OXA-92",gene_cluster_table$gene_cluster) == "TRUE", gene_cluster_table$cluster_rep <- 14,
             #                                                                                                                ifelse(grepl(oxa51,gene_cluster_table$gene_cluster) == "TRUE", gene_cluster_table$cluster_rep <- 11,0))))))))))))))


gene_cluster_table$cluster_rep <- ifelse(grepl("OXA-23",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-58",gene_cluster_table$gene_cluster) == "TRUE" & grepl(oxa51,gene_cluster_table$gene_cluster) == "TRUE", 4, 
                                        ifelse(grepl("NDM-1",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-72",gene_cluster_table$gene_cluster) == "TRUE" & grepl(oxa51,gene_cluster_table$gene_cluster) == "TRUE", 1,
                                              ifelse(grepl("OXA-23",gene_cluster_table$gene_cluster) == "TRUE" & grepl(oxa51,gene_cluster_table$gene_cluster) == "TRUE", 2,
                                                    ifelse(grepl("OXA-23",gene_cluster_table$gene_cluster) == "TRUE" & grepl("OXA-58",gene_cluster_table$gene_cluster) == "TRUE", 3,
                                                           ifelse(grepl("OXA-24",gene_cluster_table$gene_cluster) == "TRUE" & grepl(oxa51,gene_cluster_table$gene_cluster) == "TRUE", 5,
                                                                  ifelse(grepl("OXA-72",gene_cluster_table$gene_cluster) == "TRUE" & grepl(oxa51,gene_cluster_table$gene_cluster) == "TRUE", 6, 
                                                                         ifelse(grepl(oxa51,gene_cluster_table$gene_cluster) == "TRUE", gene_cluster_table$cluster_rep <- 7,
                                                                                ifelse(grepl("\\bOXA\\b",gene_cluster_table$gene_cluster) == "TRUE", 8, 10))))))))


for (i in 1:length(gene_cluster_table$Sample_id)) {
  if (gene_cluster_table$ToolDB[i] == "deeparg_SR") {
    gene_cluster_table$cluster_rep[i] <- 9
  }
} 

gene_cluster_table$cluster_rep <- as.character(gene_cluster_table$cluster_rep)

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


# order legend
#gene_cluster_table$cluster_rep <- factor(gene_cluster_table$cluster_rep, levels = c(1,2,3,4,12,5,6,9,7,10,8,13,14,11,15,0)) 
#c("OXA-23 & OXA-66", "OXA-23 & OXA-58", "OXA-23 & OXA-66 & OXA-58", 
#  "OXA-23 & OXA-68", "OXA-23 & OXA-91", "OXA-24", "OXA-72 & OXA-66", "OXA-72 & OXA-92", 
#  "NDM-1 & OXA-72 & OXA-66", "OXA-66", "OXA-90", "OXA-51",
#  "OXA-92","OXA-51 type", "mixed", "none"))

colors <- c("black", "blue", "deepskyblue3", "green3", "seagreen2", "orange", "mediumpurple1", "yellowgreen", "moccasin" )

a <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "1"]))
b <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "2"]))
c <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "3"]))
d <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "4"]))
e <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "5"]))
f <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "6"]))
g <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "7"]))

gene_cluster_table$Sample_id <- factor(gene_cluster_table$Sample_id, levels = c(a,b,c,d,e,f,g))
gene_cluster_table$cluster_rep <- factor(gene_cluster_table$cluster_rep, levels = c(1,2,3,4,5,6,7,8,9,10)) 

gene_cluster_heatmap <- ggplot(data = gene_cluster_table, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("NDM-1 & OXA-72 & OXA-51 family", "OXA-23 & OXA-51 family", "OXA-23 & OXA-58 & OXA-51 family",  
                                                "OXA-24 & OXA-51 family" , "OXA-72 & OXA-51 family",  "OXA-51 family", "undefined OXA", "mixed", "none")) +
  labs(x = "", y = "") +
  labs(fill = "Carbapenemase genes") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(gene_cluster_table$MeropenemMIC, levels = c("<= 0.25", ">= 16"))), vars(Tool), scales = "free", space = "free")

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
concordance_table$abr_ncbi <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_ncbi) == "TRUE", 3,
                                            ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_ncbi) == "FALSE", 2, 4)))

concordance_table$abr_card <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_card) == "TRUE", 3,
                                            ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_card) == "FALSE", 1, 4)))

concordance_table$abr_megares <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 1,
                                        ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_megares) == "TRUE", 3,
                                               ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_megares) == "FALSE", 2, 4)))

concordance_table$abr_resfinder <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 1,
                                          ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_resfinder) == "TRUE", 3,
                                                 ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_resfinder) == "FALSE", 2, 4)))

concordance_table$abr_argannot <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_argannot) == "TRUE", 3,
                                                ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$abr_argannot) == "FALSE", 2, 4)))

concordance_table$amrf_nuc <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$amrf_nuc) == "TRUE", 3,
                                            ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$amrf_nuc) == "FALSE", 2, 4)))

concordance_table$amrf_prot <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 1,
                                      ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$amrf_prot) == "TRUE", 3,
                                             ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$amrf_prot) == "FALSE", 2, 4)))

concordance_table$deeparg_LS <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$deeparg_LS) == "TRUE", 3,
                                              ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$deeparg_LS) == "FALSE", 2, 4)))

concordance_table$deeparg_SR <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$deeparg_SR) == "TRUE", 3,
                                              ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$deeparg_SR) == "FALSE", 2, 4)))

concordance_table$resfinder_as <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$resfinder_as) == "TRUE", 3,
                                                ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$resfinder_as) == "FALSE", 2, 4)))

concordance_table$resfinder_re <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 1,
                                         ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$resfinder_re) == "TRUE", 3,
                                                ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$resfinder_re) == "FALSE", 2, 4)))

concordance_table$rgi <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 1,
                                ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$rgi) == "TRUE", 3,
                                       ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$rgi) == "FALSE", 2, 4)))

concordance_table$srax_basic <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 1,
                                       ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$srax_basic) == "TRUE", 3,
                                              ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$srax_basic) == "FALSE", 2, 4)))

concordance_table$srax_ext <- ifelse(concordance_table$Type == "acibauc" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 1,
                                     ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$srax_ext) == "TRUE", 3,
                                            ifelse(concordance_table$Type == "acibau" & grepl("Carbapenem", concordance_table$srax_ext) == "FALSE", 2, 4)))


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
                                          ifelse(heatmap_concordance$ToolDB == "LS" | heatmap_concordance$ToolDB == "SR", "DeepARG",
                                                 ifelse(heatmap_concordance$ToolDB == "basic" | heatmap_concordance$ToolDB == "ext", "sraX",
                                                        ifelse(heatmap_concordance$ToolDB == "nuc" | heatmap_concordance$ToolDB == "prot", "AMRFinder", "ABRicate")))))

# do it after step before otherwise card of rgi is counted as abricate
#heatmap_concordance[heatmap_concordance == "rgi"] <- "card"

# create heatmap
colors <- c("darkgreen", "green", "gold", "red3")

heatmap_concordance$Sample_id <- factor(heatmap_concordance$Sample_id, levels = c(a,b,c,d,e,f,g))

concordance_heatmap <- ggplot(data = heatmap_concordance, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "S/R", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance$MeropenemMIC, levels = c("<= 0.25", ">= 16"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap


heatmap_concordance2 <- heatmap_concordance

heatmap_concordance2 <- cbind(heatmap_concordance2, gene_cluster_table[7])
heatmap_concordance2$Concordance[which(heatmap_concordance2$Concordance == 3 & heatmap_concordance2$cluster_rep == 7)] <- 2
heatmap_concordance2$Concordance[which(heatmap_concordance2$Concordance == 1 & heatmap_concordance2$cluster_rep == 7)] <- 4

concordance_heatmap2 <- ggplot(data = heatmap_concordance2, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "S/R", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_concordance2$MeropenemMIC, levels = c("<= 0.25", ">= 16"))), vars(Tool), scales = "free", space = "free")

concordance_heatmap2


#-------------------------------------------------------------------------------
# plot mero imi phenotypes

mero_imi <- pheno_mic[,c(1,2,13,32,75,94)]

mero_imi_table <- pivot_longer(data = mero_imi, 
                                    cols = -c(1:2,5:6),
                                    names_to = "Antibiotic", 
                                    values_to = "Phenotype")
res <- which(mero_imi_table$Phenotype == "NULL")
mero_imi_table$Phenotype[res] <- "R"

mero_imi_table$Sample_id <- factor(mero_imi_table$Sample_id, levels = c(a,b,c,d,e,f,g))

colors <- c("gray", "lemonchiffon", "black")
mero_imi_table$Phenotype <- factor(mero_imi_table$Phenotype, levels = c("S", "I", "R")) 

mero_imi_heatmap <- ggplot(data = mero_imi_table, mapping = aes(x = Antibiotic,
                                                          y = Sample_id,
                                                          fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("S", "I", "R")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(MIC.Meropenem), scales = "free", space = "free")

mero_imi_heatmap


mero_imi$Sample_id <- factor(mero_imi$Sample_id, levels = c(a,b,c,d,e,f,g))

colors <- c("gray", "black")

mero_mic <- ggplot(data = mero_imi, mapping = aes(x = "Meropenem MIC",
                                                      y = Sample_id,
                                                      fill = MIC.Meropenem)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi$MIC.Meropenem, levels = c("<= 0.25", ">= 16"))), scales = "free", space = "free")

mero_mic


imi_mic <- ggplot(data = mero_imi, mapping = aes(x = "Imipenem MIC",
                                                      y = Sample_id,
                                                      fill = MIC.Imipenem)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(mero_imi$MIC.Meropenem, levels = c("<= 0.25", ">= 16"))), scales = "free", space = "free")

imi_mic
