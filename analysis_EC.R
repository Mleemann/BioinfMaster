library(ggplot2)
library(tidyr)
install.packages('tidyverse')
library(tidyverse)


# read in samples to be analyzed
samples <- read.table('samples_EC.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_EC.csv', header = TRUE, sep = ';')
mic <- read.csv2('mic_EC.csv', header = TRUE, sep = ';')
pheno_mic <- read.csv2('pheno_mic_EC.csv', header = TRUE, sep = ';')

phenotypes[phenotypes == "esccol2" | phenotypes == "esccol3"] <- "esccol"
mic[mic == "esccol2" | mic == "esccol3"] <- "esccol"
pheno_mic[pheno_mic == "esccol2" | pheno_mic == "esccol3"] <- "esccol"

# check sample ids; and correct them 
#which(phenotypes$Sample_id %in% samples$V1 == FALSE)
phenotypes$Sample_id[7] <- "700109-6-13"
phenotypes$Sample_id[52] <- "708018-16-wh"
phenotypes$Sample_id[74] <- "804158-1-17"
phenotypes$Sample_id[78] <- "714843-17-wh"
phenotypes$Sample_id[95] <- "712814-18-wh"
phenotypes$Sample_id[97] <- "713302-18-wh"
phenotypes$Sample_id[102] <- "717578-1-18"
phenotypes$Sample_id[105] <- "719819-2-18"
phenotypes$Sample_id[106] <- "720780-18-wh"
phenotypes$Sample_id[155] <- "721913-18-wh"
phenotypes$Sample_id[156] <- "721915-18-wh"
phenotypes$Sample_id[176] <- "701215-7-19"
phenotypes$Sample_id[236] <- "806847-11-19"
phenotypes$Sample_id[252] <- "126031-20-wh"
phenotypes$Sample_id[254] <- "804813-1-20"
phenotypes$Sample_id[256] <- "804895-1-20"
phenotypes$Sample_id[257] <- "133707-1-20"
phenotypes$Sample_id[271] <- "808154-21-wh"

#which(mic$Sample_id %in% samples$V1 == FALSE)
mic$Sample_id[7] <- "700109-6-13"
mic$Sample_id[52] <- "708018-16-wh"
mic$Sample_id[74] <- "804158-1-17"
mic$Sample_id[78] <- "714843-17-wh"
mic$Sample_id[95] <- "712814-18-wh"
mic$Sample_id[97] <- "713302-18-wh"
mic$Sample_id[102] <- "717578-1-18"
mic$Sample_id[105] <- "719819-2-18"
mic$Sample_id[106] <- "720780-18-wh"
mic$Sample_id[155] <- "721913-18-wh"
mic$Sample_id[156] <- "721915-18-wh"
mic$Sample_id[176] <- "701215-7-19"
mic$Sample_id[236] <- "806847-11-19"
mic$Sample_id[252] <- "126031-20-wh"
mic$Sample_id[254] <- "804813-1-20"
mic$Sample_id[256] <- "804895-1-20"
mic$Sample_id[257] <- "133707-1-20"
mic$Sample_id[271] <- "808154-21-wh"

#which(pheno_mic$Sample_id %in% samples$V1 == FALSE)
pheno_mic$Sample_id[74] <- "700109-6-13"
pheno_mic$Sample_id[117] <- "708018-16-wh"
pheno_mic$Sample_id[240] <- "804158-1-17"
pheno_mic$Sample_id[131] <- "714843-17-wh"
pheno_mic$Sample_id[123] <- "712814-18-wh"
pheno_mic$Sample_id[124] <- "713302-18-wh"
pheno_mic$Sample_id[135] <- "717578-1-18"
pheno_mic$Sample_id[139] <- "719819-2-18"
pheno_mic$Sample_id[144] <- "720780-18-wh"
pheno_mic$Sample_id[189] <- "721913-18-wh"
pheno_mic$Sample_id[190] <- "721915-18-wh"
pheno_mic$Sample_id[80] <- "701215-7-19"
pheno_mic$Sample_id[263] <- "806847-11-19"
pheno_mic$Sample_id[50] <- "126031-20-wh"
pheno_mic$Sample_id[246] <- "804813-1-20"
pheno_mic$Sample_id[248] <- "804895-1-20"
pheno_mic$Sample_id[56] <- "133707-1-20"
pheno_mic$Sample_id[267] <- "808154-21-wh"


# get result tables
# abricate
abr_ncbi <- read.delim('EC_abricate_ncbi.tab', sep = "\t", header = T)
if (length(which(abr_ncbi$X.IDENTITY < 95)) != 0) {
  abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$X.IDENTITY < 95)),]
}
abr_ncbi$X.FILE <- gsub(".fna","", abr_ncbi$X.FILE)
names(abr_ncbi)[1] <- c('Sample_id')

abr_card <- read.delim('EC_abricate_card.tab', sep = "\t", header = T)
if (length(which(abr_card$X.IDENTITY < 95)) != 0) {
  abr_card <- abr_card[-c(which(abr_card$X.IDENTITY < 95)),]
}
abr_card$X.FILE <- gsub(".fna","", abr_card$X.FILE)
names(abr_card)[1] <- c('Sample_id')

abr_megares <- read.delim('EC_abricate_megares.tab', sep = "\t", header = T)
if (length(which(abr_megares$X.IDENTITY < 95)) != 0) {
  abr_megares <- abr_megares[-c(which(abr_megares$X.IDENTITY < 95)),]
}
abr_megares$X.FILE <- gsub(".fna","", abr_megares$X.FILE)
names(abr_megares)[1] <- c('Sample_id')

abr_resfinder <- read.delim('EC_abricate_resfinder.tab', sep = "\t", header = T)
if (length(which(abr_resfinder$X.IDENTITY < 95)) != 0) {
  abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$X.IDENTITY < 95)),]
}
abr_resfinder$X.FILE <- gsub(".fna","", abr_resfinder$X.FILE)
names(abr_resfinder)[1] <- c('Sample_id')

abr_argannot <- read.delim('EC_abricate_argannot.tab', sep = "\t", header = T)
if (length(which(abr_argannot$X.IDENTITY < 95)) != 0) {
  abr_argannot <- abr_argannot[-c(which(abr_argannot$X.IDENTITY < 95)),]
}
abr_argannot$X.FILE <- gsub(".fna","", abr_argannot$X.FILE)
names(abr_argannot)[1] <- c('Sample_id')

#amrfinder
amrf_nuc <- read.delim('EC_amrfinder_nuc.tab', sep = "\t")
if (length(which(amrf_nuc$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_nuc)[1] <- 'Sample_id'
amrf_nuc$Sample_id <- gsub("gnl\\|USB\\|","", amrf_nuc$Contig.id)
amrf_nuc$Sample_id <- gsub('(_\\d*)',"", amrf_nuc$Sample_id)

amrf_prot <- read.delim('EC_amrfinder_prot.tab', sep = "\t")
if (length(which(amrf_prot$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_prot <- amrf_prot[-c(which(amrf_prot$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_prot)[1] <- 'Sample_id'
amrf_prot$Sample_id <- gsub('(_\\d*)',"", amrf_prot$Sample_id)


# rgi
rgi <- read.delim(('EC_rgi.tab'), sep = "\t")
if (length(which(rgi$Best_Identities < 95)) != 0) {
  rgi <- rgi[-c(which(rgi$Best_Identities < 95)),]
}
names(rgi)[1] <- 'Sample_id'
rgi$Sample_id <- gsub("gnl\\|USB\\|","", rgi$Sample_id)
rgi$Sample_id <- gsub('(_\\d*)',"", rgi$Sample_id)


# sraX
srax_basic <- read.delim(('EC_srax_basic.tab'), sep = "\t")
if (length(which(srax_basic$Identity_p < 95)) != 0) {
  srax_basic <- srax_basic[-c(which(srax_basic$Identity_p < 95)),]
}
names(srax_basic)[2] <- 'Sample_id'
srax_basic$Sample_id <- gsub(".fna","", srax_basic$Sample_id)

srax_ext <- read.delim(('EC_srax_ext.tab'), sep = "\t")
if (length(which(srax_ext$Identity_p < 95)) != 0) {
  srax_ext <- srax_ext[-c(which(srax_ext$Identity_p < 95)),]
}
names(srax_ext)[2] <- 'Sample_id'
srax_ext$Sample_id <- gsub(".fna","", srax_ext$Sample_id)


# Deeparg
deeparg_LS <- read.delim(('EC_deeparg_LS.tab'), sep = "\t")
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
  deeparg_LS[i,6] <- genes[length(genes)]
}

deeparg_SR <- read.delim(('EC_deeparg_SR.tab'), sep = "\t")
if (length(which(deeparg_SR$identity < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$identity < 95)),]
}
if (length(which(deeparg_SR$probability < 0.8)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$probability < 0.8)),]
}
deeparg_SR_unique <- deeparg_SR %>% distinct(X.ARG, best.hit, .keep_all = TRUE)
names(deeparg_SR_unique)[1] <- 'Sample_id'
rm(deeparg_SR)
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
  deeparg_SR_unique[i,6] <- genes[length(genes)]
}
deeparg_SR_unique <- deeparg_SR_unique[-c(which(deeparg_SR_unique$best.hit == "undefined")), ]


# resfinder
resfinder_as <- read.delim(('EC_resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$Identity < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$Identity < 95)),]
}
names(resfinder_as)[6] <- 'Sample_id'
resfinder_as$Sample_id <- gsub('gnl\\|USB\\|','', resfinder_as$Sample_id)
resfinder_as$Sample_id <- gsub('(_\\d*)','', resfinder_as$Sample_id)

resfinder_re <- read.delim(('EC_resfinder_reads.tab'), sep = "\t", header = F)
names(resfinder_re) <- c('Resistance gene',	'Identity',	'Alignment Length/Gene Length',	'Coverage',	'Position in reference',	'Contig',	'Position in contig',  'Phenotype',	'Accession no.', 'Sample_id')
resfinder_re <- resfinder_re[-c(1),]
if (length(which(resfinder_re$Identity < 95)) != 0) {
  resfinder_re <- resfinder_re[-c(which(resfinder_re$Identity < 95)),]
}


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
resfinder_bla <-abr_resfinder[c(grep("bla", abr_resfinder$PRODUCT)), ]
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

resfinder_as_bla <- resfinder_as[c(grep("bla", resfinder_as$Resistance.gene)), ]
names(resfinder_as_bla)[1] <- 'bla_genes'
resfinder_as_bla$Tool <- 'resfinder_as'
resfinder_re_bla <- resfinder_re[c(grep("bla", resfinder_re$`Resistance gene`)), ]
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

# include clinical phenotype 
for (i in 1:length(bla$Sample_id)){
  bla$Type[i] <- phenotypes$LGM_BK[grep(bla$Sample_id[i], phenotypes$Sample_id)]
  }

bla <- unique(bla)

#-------------------------------------------------------------------------------
# plot phenotypes of all antibiotics
phenotypes_filtered <- phenotypes
phenotypes_filtered$Dummy <- NULL
phenotypes_filtered[phenotypes_filtered == "NULL" | phenotypes_filtered == "U"] <- 0
phenotypes_filtered <- Filter(function(x) !(all(x==0)), phenotypes_filtered)

#unique(as.vector(as.matrix(phenotypes_filtered)))
phenotypes_filtered[phenotypes_filtered == "c(R, R)" | phenotypes_filtered == "c(S, R)" | phenotypes_filtered == "c(R, S)" | phenotypes_filtered == "c(S, U, R)" |   
                      phenotypes_filtered == "c(R, S, S)" | phenotypes_filtered == "c(R, R, R)" | phenotypes_filtered == "c(R, U)" | 
                      phenotypes_filtered == "c(U, R)" | phenotypes_filtered == "c(U, R, U, R)" | phenotypes_filtered == "c(R, R, I)" |
                      phenotypes_filtered == "c(R, U, U, I, I)" | phenotypes_filtered == "c(I, U, U, R)" | phenotypes_filtered == "c(R, U, U, I)"] <- "R"

phenotypes_filtered[phenotypes_filtered == "c(S, S)" | phenotypes_filtered == "c(U, S)" | phenotypes_filtered == "c(S, S, S)" | phenotypes_filtered == "c(S, U)" |
                      phenotypes_filtered == "c(S, U, S)" | phenotypes_filtered == "c(S, S, U)" | phenotypes_filtered == "c(S, U, U, S)" |
                      phenotypes_filtered == "c(U, S, S)"] <- "S"

phenotypes_filtered[phenotypes_filtered == "c(I, S)" | phenotypes_filtered == "c(S, I)" | phenotypes_filtered == "c(I, S, S)" |
                      phenotypes_filtered == "c(S, I, S, I)" | phenotypes_filtered == "c(I, I)" | phenotypes_filtered == "c(I, U)" | 
                      phenotypes_filtered == "c(U, I)" | phenotypes_filtered == "c(U, I, U, U)" | phenotypes_filtered == "c(U, I, U, S)" |
                      phenotypes_filtered == "c(U, S, U, I)" | phenotypes_filtered == "c(U, I, U)" | phenotypes_filtered == "c(I, I, U)" |
                      phenotypes_filtered == "c(I, U, U)" | phenotypes_filtered == "c(U, S, I)"] <- "I"

phenotypes_filtered[phenotypes_filtered == "c(U, U)" | phenotypes_filtered == "c(U, U, U)"] <- 0

# combine antibiotics with several columns
phenotypes_filtered$Meropenem[phenotypes_filtered$Meropenem.bei.Meningitis  == "S"] <- "S"
phenotypes_filtered$Meropenem[phenotypes_filtered$Meropenem.bei.Meningitis  == "R"] <- "R"
phenotypes_filtered$Meropenem[phenotypes_filtered$Meropenem.bei.Meningitis  == "I"] <- "I"

phenotypes_filtered$Cefuroxim[phenotypes_filtered$Cefuroxim.bei.Meningitis  == "S"] <- "S"
phenotypes_filtered$Cefuroxim[phenotypes_filtered$Cefuroxim.bei.Meningitis  == "R"] <- "R"
phenotypes_filtered$Cefuroxim[phenotypes_filtered$Cefuroxim.bei.Meningitis  == "I"] <- "I"

phenotypes_filtered <- phenotypes_filtered[,-c(16, 27, 42, 44, 45, 46, 47, 48)]

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


#-------------------------------------------------------------------------------
# plot bla genes

bla_for_heatmap <- bla
bla_for_heatmap$bla_genes <- gsub("\\(Bla\\)", "", bla_for_heatmap$bla_genes)
bla_for_heatmap$bla_genes <- gsub("bla", "", bla_for_heatmap$bla_genes, ignore.case=T)

bla_for_heatmap[bla_for_heatmap == "KPC-1" | bla_for_heatmap == "KPC-2"] <- "KPC-1/2"

bla_for_heatmap$Occurance <- 1
bla_for_heatmap$Occurance <- as.character(bla_for_heatmap$Occurance)

# classification of the found bla genes
carba <- read.table('carba.txt')
esbl <- read.table('esbl.txt')
ampc <- read.table('ampc.txt')

for (i in 1:length(bla_for_heatmap$bla_genes)) { 
  ifelse(bla_for_heatmap$bla_genes[i] %in% carba$V1, bla_for_heatmap$Classification[i] <- 'Carbapenem', 
         ifelse(bla_for_heatmap$bla_genes[i] %in% esbl$V1, bla_for_heatmap$Classification[i] <- 'ESBL',
                ifelse(bla_for_heatmap$bla_genes[i] %in% ampc$V1, bla_for_heatmap$Classification[i] <- 'AmpC', 
                       bla_for_heatmap$Classification[i] <- 'Beta-lactam')))
}

bla_heatmap <- ggplot(data = bla_for_heatmap, mapping = aes(x = bla_genes,
                                                            y = Sample_id,
                                                            fill = Occurance)) +
  geom_tile(colour="white", size=1) +
  labs(x = "", y = "") +
  #scale_fill_manual(values = colors2) +
  theme(legend.position = 'none', axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 6)) + 
  facet_grid(vars(Type), vars(Classification), scales = "free", space = "free")
  #ggtitle()

bla_heatmap

# show how many tools have found each gene
for (i in 1:length(bla_for_heatmap$Sample_id)) { 
  j <- bla_for_heatmap$bla_genes[i]
  bla_for_heatmap$Occurance[i] <- length(which(bla_for_heatmap$Sample_id == bla_for_heatmap$Sample_id[i] & bla_for_heatmap$bla_genes == j))
}
bla_for_heatmap$Occurance <- as.character(bla_for_heatmap$Occurance)

bla_heatmap2 <- ggplot(data = bla_for_heatmap[1:3000,], mapping = aes(x = bla_genes,
                                                            y = Sample_id,
                                                            fill = Occurance)) +
  geom_tile(colour="white", size=1) +
  labs(x = "", y = "") +
  #scale_fill_manual(values = colors2) +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 6)) + 
  facet_grid(vars(Type), vars(Classification), scales = "free", space = "free")
#ggtitle()

bla_heatmap2


# test abricate tools
abricate_bla <- bla_for_heatmap[grep('abr', bla_for_heatmap$Tool),]
abricate_bla <- abricate_bla[-c(grep('ncbi', abricate_bla$Tool)),]

for (i in 1:length(abricate_bla$Sample_id)) { 
  j <- abricate_bla$bla_genes[i]
  abricate_bla$Occurance[i] <- length(which(abricate_bla$Sample_id == abricate_bla$Sample_id[i] & abricate_bla$bla_genes == j))
}
abricate_bla$Occurance <- as.character(abricate_bla$Occurance)

abricate_bla_heatmap <- ggplot(data = abricate_bla, mapping = aes(x = bla_genes,
                                                                      y = Sample_id,
                                                                      fill = Occurance)) +
  geom_tile(colour="white", size=1) +
  labs(x = "", y = "") +
  #scale_fill_manual(values = colors2) +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 6)) + 
  facet_grid(vars(Type), vars(Classification), scales = "free", space = "free")
#ggtitle()

abricate_bla_heatmap
