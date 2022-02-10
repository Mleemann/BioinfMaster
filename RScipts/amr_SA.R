#install.packages('tidyverse')
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)

setwd("~/data_project/assemblies/SA")

# read in samples to be analyzed
samples <- read.table('samples_SA.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_SA.csv', header = TRUE, sep = ';')
mic <- read.csv2('mic_SA.csv', header = TRUE, sep = ';')
pheno_mic <- read.csv2('pheno_mic_SA.csv', header = TRUE, sep = ';')

# check sample ids; and correct them 
which(phenotypes$Sample_id %in% samples$V1 == FALSE)
phenotypes$Sample_id[89] <- "504423-19-wh"
phenotypes$Sample_id[158] <- "503081-1-20"
phenotypes$Sample_id[159] <- "402988-1-20"
phenotypes$Sample_id[160] <- "608029-1-20"
phenotypes$Sample_id[161] <- "503156-1-20"
phenotypes$Sample_id[162] <- "805048-1-20-wh"
phenotypes$Sample_id[176] <- "806818-20-wh"
phenotypes$Sample_id[215] <- "606179-21-wh"
phenotypes$Sample_id[229] <- "804519-01"
phenotypes$Sample_id[238] <- "608085-21-wh"
phenotypes$Sample_id[251] <- "205044-21-wh"
phenotypes$Sample_id[258] <- "808152-21-wh"
phenotypes$Sample_id[262] <- "404692-21-wh"
phenotypes$Sample_id[263] <- "808201-21-wh"

phenotypes <- phenotypes[-c(which(phenotypes$Sample_id == "808339-21")),]

#which(mic$Sample_id %in% samples$V1 == FALSE)
mic$Sample_id[89] <- "504423-19-wh"
mic$Sample_id[158] <- "503081-1-20"
mic$Sample_id[159] <- "402988-1-20"
mic$Sample_id[160] <- "608029-1-20"
mic$Sample_id[161] <- "503156-1-20"
mic$Sample_id[162] <- "805048-1-20-wh"
mic$Sample_id[176] <- "806818-20-wh"
mic$Sample_id[215] <- "606179-21-wh"
mic$Sample_id[229] <- "804519-01"
mic$Sample_id[238] <- "608085-21-wh"
mic$Sample_id[251] <- "205044-21-wh"
mic$Sample_id[258] <- "808152-21-wh"
mic$Sample_id[262] <- "404692-21-wh"
mic$Sample_id[263] <- "808201-21-wh"

mic <- mic[-c(which(mic$Sample_id == "808339-21")),]

which(pheno_mic$Sample_id %in% samples$V1 == FALSE)
pheno_mic$Sample_id[27] <- "205044-21-wh"
pheno_mic$Sample_id[49] <- "402988-1-20"
pheno_mic$Sample_id[68] <- "404692-21-wh"
pheno_mic$Sample_id[85] <- "503081-1-20"
pheno_mic$Sample_id[86] <- "503156-1-20"
pheno_mic$Sample_id[101] <- "504423-19-wh"
pheno_mic$Sample_id[121] <- "606179-21-wh"
pheno_mic$Sample_id[132] <- "608029-1-20"
pheno_mic$Sample_id[134] <- "608085-21-wh"
pheno_mic$Sample_id[219] <- "804519-01"
pheno_mic$Sample_id[222] <- "805048-1-20-wh"
pheno_mic$Sample_id[243] <- "806818-20-wh"
pheno_mic$Sample_id[251] <- "808152-21-wh"
pheno_mic$Sample_id[253] <- "808201-21-wh"

pheno_mic <- pheno_mic[-c(which(pheno_mic$Sample_id == "808339-21")),]


pheno_mic[pheno_mic == "c(R, R)"] <- "R"
pheno_mic[pheno_mic == "c(S, S)"] <- "S"
pheno_mic[pheno_mic == "c(U, R)"] <- "R"
pheno_mic[pheno_mic == "c(R, U)"] <- "R"
pheno_mic[pheno_mic == "c(R, S)"] <- "R"
pheno_mic[pheno_mic == "c(S, R)"] <- "R"

pheno_mic[pheno_mic == "c(<= 0.25, = 0.38)" | pheno_mic == "c(<= 0.25, <= 0.25)" |
             pheno_mic == "c(<= 0.25, = 0.5)" | pheno_mic == "c(<= 0.25, = 0.25)"] <- "<= 0.25"
pheno_mic[pheno_mic == "c(= 1, = 1)"] <- "= 1"


# get result tables
# abricate
abr_ncbi <- read.delim('abricate_ncbi.tab', sep = "\t", header = F)
if (length(which(abr_ncbi$V11 < 95)) != 0) {
  abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$V11 < 95)),]
}
abr_ncbi$V1 <- gsub(".fna","", abr_ncbi$V1)
names(abr_ncbi)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_ncbi$Tool <- 'abr_ncbi'

abr_card <- read.delim('abricate_card.tab', sep = "\t", header = F)
if (length(which(abr_card$V11 < 95)) != 0) {
  abr_card <- abr_card[-c(which(abr_card$V11 < 95)),]
}
abr_card$V1 <- gsub(".fna","", abr_card$V1)
names(abr_card)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_card$Tool <- 'abr_card'

abr_megares <- read.delim('abricate_megares.tab', sep = "\t", header = F)
if (length(which(abr_megares$V11 < 95)) != 0) {
  abr_megares <- abr_megares[-c(which(abr_megares$V11 < 95)),]
}
abr_megares$V1 <- gsub(".fna","", abr_megares$V1)
names(abr_megares)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_megares$Tool <- 'abr_megares'

abr_resfinder <- read.delim('abricate_resfinder.tab', sep = "\t", header = F)
if (length(which(abr_resfinder$V11 < 95)) != 0) {
  abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$V11 < 95)),]
}
abr_resfinder$V1 <- gsub(".fna","", abr_resfinder$V1)
names(abr_resfinder)[c(1,14)] <- c('Sample_id', 'amr_genes')
abr_resfinder$Tool <- 'abr_resfinder'

abr_argannot <- read.delim('abricate_argannot.tab', sep = "\t", header = F)
if (length(which(abr_argannot$V11 < 95)) != 0) {
  abr_argannot <- abr_argannot[-c(which(abr_argannot$V11 < 95)),]
}
abr_argannot$V1 <- gsub(".fna","", abr_argannot$V1)
names(abr_argannot)[c(1,6)] <- c('Sample_id', 'amr_genes')
abr_argannot$Tool <- 'abr_argannot'


#amrfinder
amrf_nuc <- read.delim('amrfinder_nuc.tab', sep = "\t")
if (length(which(amrf_nuc$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_nuc)[1] <- 'Sample_id'
amrf_nuc$Sample_id <- gsub("gnl\\|USB\\|","", amrf_nuc$Contig.id)
amrf_nuc$Sample_id <- gsub('(_\\d*)',"", amrf_nuc$Sample_id)
names(amrf_nuc)[6] <- 'amr_genes'
amrf_nuc$Tool <- 'amrf_nuc'

amrf_prot <- read.delim('amrfinder_prot.tab', sep = "\t")
if (length(which(amrf_prot$X..Identity.to.reference.sequence < 95)) != 0) {
  amrf_prot <- amrf_prot[-c(which(amrf_prot$X..Identity.to.reference.sequence < 95)),]
}
names(amrf_prot)[c(1,2)] <- c('Sample_id', 'amr_genes')
amrf_prot$Sample_id <- gsub('(_\\d*)',"", amrf_prot$Sample_id)
amrf_prot$Tool <- 'amrf_prot'


# rgi
rgi <- read.delim(('rgi.tab'), sep = "\t")
if (length(which(rgi$Best_Identities < 95)) != 0) {
  rgi <- rgi[-c(which(rgi$Best_Identities < 95)),]
}
names(rgi)[c(1,2)] <- c('Sample_id', 'amr_genes')
rgi$Sample_id <- gsub("gnl\\|USB\\|","", rgi$Sample_id)
rgi$Sample_id <- gsub('(_\\d*)',"", rgi$Sample_id)
rgi$Sample_id <- gsub('\\s+', '', rgi$Sample_id)
rgi$Tool <- 'rgi'


# sraX
srax_basic <- read.delim(('srax_basic.tab'), sep = "\t")
if (length(which(srax_basic$Identity_p < 95)) != 0) {
  srax_basic <- srax_basic[-c(which(srax_basic$Identity_p < 95)),]
}
names(srax_basic)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_basic$Sample_id <- gsub(".fna","", srax_basic$Sample_id)
srax_basic$Tool <- 'sraX_basic'

srax_ext <- read.delim(('srax_ext.tab'), sep = "\t")
if (length(which(srax_ext$Identity_p < 95)) != 0) {
  srax_ext <- srax_ext[-c(which(srax_ext$Identity_p < 95)),]
}
names(srax_ext)[c(2,6)] <- c('Sample_id', 'amr_genes')
srax_ext$Sample_id <- gsub(".fna","", srax_ext$Sample_id)
srax_ext$Tool <- 'sraX_ext'


# Deeparg
deeparg_LS <- read.delim(('deeparg_LS_noPot.tab'), sep = "\t")
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
names(deeparg_LS)[6] <- 'amr_genes'
deeparg_LS$Tool <- 'deeparg_LS'


deeparg_SR <- read.delim(('deeparg_SR_noPot.tab'), sep = "\t")
if (length(which(deeparg_SR$identity < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$identity < 95)),]
}
if (length(which(deeparg_SR$probability < 0.8)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$probability < 0.8)),]
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
  deeparg_SR_unique[i,6] <- genes[3]
}
# deeparg_SR_unique <- deeparg_SR_unique[-c(which(deeparg_SR_unique$best.hit == "undefined")), ]
names(deeparg_SR_unique)[6] <- 'amr_genes'
deeparg_SR_unique$Tool <- 'deeparg_SR'

rm(deeparg_SR)

# resfinder
resfinder_as <- read.delim(('resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$Identity < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$Identity < 95)),]
}
names(resfinder_as)[c(6,1)] <- c('Sample_id', 'amr_genes')
resfinder_as$Sample_id <- gsub('gnl\\|USB\\|','', resfinder_as$Sample_id)
resfinder_as$Sample_id <- gsub('(_\\d*)','', resfinder_as$Sample_id)
resfinder_as$Tool <- 'resfinder_as'

resfinder_re <- read.delim(('resfinder_reads.tab'), sep = "\t")
if (length(which(resfinder_re$Identity < 95)) != 0) {
  resfinder_re <- resfinder_re[-c(which(resfinder_re$Identity < 95)),]
}
names(resfinder_re)[1] <- 'amr_genes'
resfinder_re$Tool <- 'resfinder_re'


#-------------------------------------------------------------------------------
# combine all data into one table
all_genes <- rbind(abr_ncbi[,c(1,6,16)], abr_card[,c(1,6,16)], abr_megares[,c(1,6,15)], abr_resfinder[,c(1,14,16)], abr_argannot[,c(1,6,15)],
                   amrf_nuc[,c(1,6,23)], amrf_prot[,c(1,2,19)], deeparg_LS[,c(4,6,13)], deeparg_SR_unique[,c(1,6,13)],
                   resfinder_as[,c(6,1,10)], resfinder_re[,c(10,1,11)], rgi[,c(1,2,10)], srax_basic[,c(2,6,17)], srax_ext[,c(2,6,17)])

# remove samples that have not passed quality control
all_genes <- all_genes[-c(which(all_genes$Sample_id == "402449-20" | all_genes$Sample_id == "610751-19" |
                                  all_genes$Sample_id == "600352-20" | all_genes$Sample_id == "719192-20" | 
                                  all_genes$Sample_id == "806052-21" | all_genes$Sample_id == "808339-21")),]


#-------------------------------------------------------------------------------
# methicillin resistance genes 

meth_all <- all_genes
meth_all[meth_all == "(Bla)mecA" | meth_all == "MECA"] <- "mecA"
meth_all[meth_all == "(Bla)mecC"] <- "mecC"
meth_all[meth_all == "MECI"] <- "mecI"
meth_all[meth_all == "MecC-type methicillin resistance repressor MecI"] <- "mecI_of_mecC"

meth_all <- unique(meth_all)
meth_all <- meth_all[c(which(meth_all$amr_genes == "mecA" | meth_all$amr_genes == "mecC" |
                               meth_all$amr_genes == "mecI" | meth_all$amr_genes == "mecI_of_mecC" |
                               meth_all$amr_genes == "mecI_of_mecA" | meth_all$amr_genes == "mecR1" |
                               meth_all$amr_genes == "PBP2")), ]

colnames(meth_all) <- c("Sample_id", "meth_genes", "ToolDB")  

gene_clusters <- pivot_wider(meth_all, names_from = ToolDB, values_from = meth_genes)

for (i in samples$V1[which(samples$V1 %in% gene_clusters$Sample_id == "FALSE")]) {
  j <- length(gene_clusters$Sample_id) + 1
  gene_clusters[j,1] <- i
}

# sort gene clusters
for (i in 1:length(gene_clusters$abr_ncbi)) {
  gene_clusters$abr_ncbi[[i]] <- sort(gene_clusters$abr_ncbi[[i]])
}

for (i in 1:length(gene_clusters$abr_card)) {
  gene_clusters$abr_card[[i]] <- sort(gene_clusters$abr_card[[i]])
}

for (i in 1:length(gene_clusters$abr_megares)) {
  gene_clusters$abr_megares[[i]] <- sort(gene_clusters$abr_megares[[i]])
}

for (i in 1:length(gene_clusters$abr_resfinder)) {
  gene_clusters$abr_resfinder[[i]] <- sort(gene_clusters$abr_resfinder[[i]])
}

for (i in 1:length(gene_clusters$abr_argannot)) {
  gene_clusters$abr_argannot[[i]] <- sort(gene_clusters$abr_argannot[[i]])
}

for (i in 1:length(gene_clusters$amrf_nuc)) {
  gene_clusters$amrf_nuc[[i]] <- sort(gene_clusters$amrf_nuc[[i]])
}

for (i in 1:length(gene_clusters$amrf_prot)) {
  gene_clusters$amrf_prot[[i]] <- sort(gene_clusters$amrf_prot[[i]])
}

for (i in 1:length(gene_clusters$deeparg_LS)) {
  gene_clusters$deeparg_LS[[i]] <- sort(gene_clusters$deeparg_LS[[i]])
}

for (i in 1:length(gene_clusters$deeparg_SR)) {
  gene_clusters$deeparg_SR[[i]] <- sort(gene_clusters$deeparg_SR[[i]])
}

for (i in 1:length(gene_clusters$resfinder_as)) {
  gene_clusters$resfinder_as[[i]] <- sort(gene_clusters$resfinder_as[[i]])
}

for (i in 1:length(gene_clusters$resfinder_re)) {
  gene_clusters$resfinder_re[[i]] <- sort(gene_clusters$resfinder_re[[i]])
}

for (i in 1:length(gene_clusters$rgi)) {
  gene_clusters$rgi[[i]] <- sort(gene_clusters$rgi[[i]])
}

for (i in 1:length(gene_clusters$sraX_basic)) {
  gene_clusters$sraX_basic[[i]] <- sort(gene_clusters$sraX_basic[[i]])
}

for (i in 1:length(gene_clusters$sraX_ext)) {
  gene_clusters$sraX_ext[[i]] <- sort(gene_clusters$sraX_ext[[i]])
}


gene_cluster_table <- pivot_longer(data = gene_clusters, 
                                   cols = -c(1),
                                   names_to = "ToolDB", 
                                   values_to = "gene_cluster")

# include clinical phenotype 
for (i in 1:length(gene_cluster_table$Sample_id)){
  gene_cluster_table$LGM_BK[i] <- pheno_mic$LGM_BK[which(gene_cluster_table$Sample_id[i] == pheno_mic$Sample_id)]
}

for (i in 1:length(gene_cluster_table$Sample_id)){
  gene_cluster_table$MIC.Oxacillin[i] <- pheno_mic$MIC.Oxacillin[which(gene_cluster_table$Sample_id[i] == pheno_mic$Sample_id)]
}


gene_cluster_table$cluster_nr <- 0
ind <- unique(gene_cluster_table$gene_cluster)

for (i in 1:length(gene_cluster_table$gene_cluster)) {
  for (j in 1:length(ind)) {
    if(identical(ind[[j]],gene_cluster_table$gene_cluster[[i]])) {
      gene_cluster_table$cluster_nr[i] <- j
    }
  }
}

gene_cluster_table$cluster_nr <- as.character(gene_cluster_table$cluster_nr)

# mecA: cluster_rep = 1, mecA + regulator cluster_rep = 2
# mecC: cluster_rep = 3, mecC + regulator cluster_rep = 4
# PBP2 = 5; mecA and mecC: 6, only reg:7, no real determinant = 0 (this includes singel vanS, vanR)

rep_one <- c(1)
rep_two <- c(2,3,4,5)
rep_three <- c(7)
rep_four <- c(6,9,12)
rep_five <- c(8)
rep_six <- c(11)
rep_seven <- c(13)

gene_cluster_table$cluster_rep <- ifelse(gene_cluster_table$cluster_nr %in% rep_one, gene_cluster_table$cluster_rep <- 1,
                                         ifelse(gene_cluster_table$cluster_nr %in% rep_two, gene_cluster_table$cluster_rep <- 2,
                                                ifelse(gene_cluster_table$cluster_nr %in% rep_three, gene_cluster_table$cluster_rep <- 3, 
                                                       ifelse(gene_cluster_table$cluster_nr %in% rep_four, gene_cluster_table$cluster_rep <- 4,
                                                              ifelse(gene_cluster_table$cluster_nr %in% rep_five, gene_cluster_table$cluster_rep <- 5, 
                                                                     ifelse(gene_cluster_table$cluster_nr %in% rep_six, gene_cluster_table$cluster_rep <- 6, 
                                                                            ifelse(gene_cluster_table$cluster_nr %in% rep_seven, gene_cluster_table$cluster_rep <- 7, gene_cluster_table$cluster_rep <- 0))))))) 

gene_cluster_table$cluster_rep <- as.character(gene_cluster_table$cluster_rep)

gene_cluster_table$Tool <- ifelse(gene_cluster_table$ToolDB == "rgi", gene_cluster_table$Tool <- "RGI",
                                  ifelse(gene_cluster_table$ToolDB == "sraX_basic" | gene_cluster_table$ToolDB == "sraX_ext", gene_cluster_table$Tool <- "sraX",
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
gene_cluster_table[gene_cluster_table == "sraX_basic"] <- "basic"
gene_cluster_table[gene_cluster_table == "sraX_ext"] <- "ext"
gene_cluster_table[gene_cluster_table == "deeparg_LS"] <- "LS"
gene_cluster_table[gene_cluster_table == "deeparg_SR"] <- "SR"
gene_cluster_table[gene_cluster_table == "resfinder_as"] <- "assembly"
gene_cluster_table[gene_cluster_table == "resfinder_re"] <- "reads"

# order legend
gene_cluster_table$cluster_rep <- factor(gene_cluster_table$cluster_rep, levels = c(1,2,3,4,5,6,7,0)) 

a <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & (gene_cluster_table$cluster_rep == "1" | gene_cluster_table$cluster_rep == "2")]))
b <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & (gene_cluster_table$cluster_rep == "3" | gene_cluster_table$cluster_rep == "4")]))
c <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & (gene_cluster_table$cluster_rep == "0")]))

gene_cluster_table$Sample_id <- factor(gene_cluster_table$Sample_id, levels = c(a, b, c)) 
colors <- c("blue4", "blue1", "red3", "sienna1", "turquoise", "magenta", "yellowgreen","moccasin")

gene_cluster_heatmap <- ggplot(data = gene_cluster_table, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("mecA", "mecA + regulator", "mecC", "mecC + regulator", "PBP2", "mecA & mecC", "only regulator","none")) +
  labs(x = "", y = "") +
  labs(fill = "Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(MIC.Oxacillin), vars(Tool), scales = "free", space = "free")

gene_cluster_heatmap


#-------------------------------------------------------------------------------
## analyse oxacillin resistance

# extract oxacillin pattern
SA_ox <- pheno_mic[,c(1,2,27,93)]

# get all samples with either mecA or mecC
meth_genes <- all_genes
meth_genes[meth_genes == "(Bla)mecA"] <- "mecA"
meth_genes[meth_genes == "(Bla)mecC"] <- "mecC"
meth_genes[meth_genes == "MECA"] <- "mecA"

meth_genes <- meth_genes[c(which(meth_genes$amr_genes == 'mecA' | meth_genes$amr_genes == 'mecC')),]
meth_genes <- unique(meth_genes)

SA_ox$abr_ncbi <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'abr_ncbi')]
SA_ox$abr_card <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'abr_card')]
SA_ox$abr_megares <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'abr_megares')]
SA_ox$abr_resfinder <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'abr_resfinder')]
SA_ox$abr_argannot <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'abr_argannot')]
SA_ox$amrf_nuc <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'amrf_nuc')]
SA_ox$amrf_prot <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'amrf_prot')]
SA_ox$deeparg_LS <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'deeparg_LS')]
SA_ox$deeparg_SR <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'deeparg_SR')]
SA_ox$resfinder_as <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'resfinder_as')]
SA_ox$resfinder_re <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'resfinder_re')]
SA_ox$rgi <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'rgi')]
SA_ox$sraX_basic <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'sraX_basic')]
SA_ox$sraX_ext <- SA_ox$Sample_id %in% meth_genes$Sample_id[which(meth_genes$Tool == 'sraX_ext')]


# Concordance values: R/R = 1, S/S = 2, S/R = 3, R/S = 4
SA_ox$abr_ncbi <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_ncbi == "FALSE", 2,
                         ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_ncbi == "TRUE", 3,
                                ifelse(SA_ox$Oxacillin == "R" & SA_ox$abr_ncbi == "TRUE", 1, 4)))

SA_ox$abr_card <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_card == "FALSE", 2,
                         ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_card == "TRUE", 3,
                                ifelse(SA_ox$Oxacillin == "R" & SA_ox$abr_card == "TRUE", 1, 4)))

SA_ox$abr_megares <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_megares == "FALSE", 2,
                            ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_megares == "TRUE", 3,
                                   ifelse(SA_ox$Oxacillin == "R" & SA_ox$abr_megares == "TRUE", 1, 4)))

SA_ox$abr_resfinder <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_resfinder == "FALSE", 2,
                              ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_resfinder == "TRUE", 3,
                                     ifelse(SA_ox$Oxacillin == "R" & SA_ox$abr_resfinder == "TRUE", 1, 4)))

SA_ox$abr_argannot <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_argannot == "FALSE", 2,
                             ifelse(SA_ox$Oxacillin == "S" & SA_ox$abr_argannot == "TRUE", 3,
                                    ifelse(SA_ox$Oxacillin == "R" & SA_ox$abr_argannot == "TRUE", 1, 4)))

SA_ox$amrf_nuc <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$amrf_nuc == "FALSE", 2,
                         ifelse(SA_ox$Oxacillin == "S" & SA_ox$amrf_nuc == "TRUE", 3,
                                ifelse(SA_ox$Oxacillin == "R" & SA_ox$amrf_nuc == "TRUE", 1, 4)))

SA_ox$amrf_prot <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$amrf_prot == "FALSE", 2,
                          ifelse(SA_ox$Oxacillin == "S" & SA_ox$amrf_prot == "TRUE", 3,
                                 ifelse(SA_ox$Oxacillin == "R" & SA_ox$amrf_prot == "TRUE", 1, 4)))

SA_ox$rgi <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$rgi == "FALSE", 2,
                    ifelse(SA_ox$Oxacillin == "S" & SA_ox$rgi == "TRUE", 3,
                           ifelse(SA_ox$Oxacillin == "R" & SA_ox$rgi == "TRUE", 1, 4)))

SA_ox$sraX_basic <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$sraX_basic == "FALSE", 2,
                           ifelse(SA_ox$Oxacillin == "S" & SA_ox$sraX_basic == "TRUE", 3,
                                  ifelse(SA_ox$Oxacillin == "R" & SA_ox$sraX_basic == "TRUE", 1, 4)))

SA_ox$sraX_ext <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$sraX_ext == "FALSE", 2,
                         ifelse(SA_ox$Oxacillin == "S" & SA_ox$sraX_ext == "TRUE", 3,
                                ifelse(SA_ox$Oxacillin == "R" & SA_ox$sraX_ext == "TRUE", 1, 4)))

SA_ox$deeparg_LS <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$deeparg_LS == "FALSE", 2,
                           ifelse(SA_ox$Oxacillin == "S" & SA_ox$deeparg_LS == "TRUE", 3,
                                  ifelse(SA_ox$Oxacillin == "R" & SA_ox$deeparg_LS == "TRUE", 1, 4)))

SA_ox$deeparg_SR <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$deeparg_SR == "FALSE", 2,
                           ifelse(SA_ox$Oxacillin == "S" & SA_ox$deeparg_SR == "TRUE", 3,
                                  ifelse(SA_ox$Oxacillin == "R" & SA_ox$deeparg_SR == "TRUE", 1, 4)))

SA_ox$resfinder_as <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$resfinder_as == "FALSE", 2,
                             ifelse(SA_ox$Oxacillin == "S" & SA_ox$resfinder_as == "TRUE", 3,
                                    ifelse(SA_ox$Oxacillin == "R" & SA_ox$resfinder_as == "TRUE", 1, 4)))

SA_ox$resfinder_re <- ifelse(SA_ox$Oxacillin == "S" & SA_ox$resfinder_re == "FALSE", 2,
                             ifelse(SA_ox$Oxacillin == "S" & SA_ox$resfinder_re == "TRUE", 3,
                                    ifelse(SA_ox$Oxacillin == "R" & SA_ox$resfinder_re == "TRUE", 1, 4)))

# transform data for heatmap
heatmap_ox <- pivot_longer(data = SA_ox, 
                           cols = -c(1:4),
                           names_to = "ToolDB", 
                           values_to = "Concordance")

heatmap_ox$Concordance <- as.character(heatmap_ox$Concordance)
heatmap_ox[heatmap_ox == "abr_ncbi"] <- "ncbi"
heatmap_ox[heatmap_ox == "abr_card"] <- "card"
heatmap_ox[heatmap_ox == "abr_megares"] <- "megares"
heatmap_ox[heatmap_ox == "abr_resfinder"] <- "resfinder"
heatmap_ox[heatmap_ox == "abr_argannot"] <- "argannot"
heatmap_ox[heatmap_ox == "amrf_nuc"] <- "nuc"
heatmap_ox[heatmap_ox == "amrf_prot"] <- "prot"
heatmap_ox[heatmap_ox == "sraX_basic"] <- "basic"
heatmap_ox[heatmap_ox == "sraX_ext"] <- "ext"
heatmap_ox[heatmap_ox == "deeparg_LS"] <- "LS"
heatmap_ox[heatmap_ox == "deeparg_SR"] <- "SR"
heatmap_ox[heatmap_ox == "resfinder_as"] <- "assembly"
heatmap_ox[heatmap_ox == "resfinder_re"] <- "reads"


# include column for the different tools
heatmap_ox$Tool <- ifelse(heatmap_ox$ToolDB == "rgi", "RGI",
                          ifelse(heatmap_ox$ToolDB == "assembly" | heatmap_ox$ToolDB == "reads", "ResFinder",
                                 ifelse(heatmap_ox$ToolDB == "LS" | heatmap_ox$ToolDB == "SR", "deepARG",
                                        ifelse(heatmap_ox$ToolDB == "basic" | heatmap_ox$ToolDB == "ext", "sraX",
                                               ifelse(heatmap_ox$ToolDB == "nuc" | heatmap_ox$ToolDB == "prot", "AMRFinder", "ABRicate")))))

# do it after step before otherwise card of rgi is counted as abricate
# heatmap_ox$ToolDB[heatmap_ox$ToolDB == "rgi"] <- "card"

# create heatmap
colors <- c("darkgreen", "green", "gold", "red3")

OX_heatmap <- ggplot(data = heatmap_ox, mapping = aes(x = ToolDB,
                                                      y = Sample_id,
                                                      fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "S/R", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_ox$LGM_BK, levels = c('staaur', 'mrsa'))), vars(Tool), scales = "free", space = "free")

OX_heatmap

# MIC as facet

heatmap_ox$Sample_id <- factor(heatmap_ox$Sample_id, levels = c(a, b, c)) 

OX_heatmap2 <- ggplot(data = heatmap_ox, mapping = aes(x = ToolDB,
                                                       y = Sample_id,
                                                       fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "S/R", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Resistance intepretation/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(MIC.Oxacillin), vars(Tool), scales = "free", space = "free")

OX_heatmap2

SA_ox %>% count(abr_ncbi)
SA_ox %>% count(abr_megares)
SA_ox %>% count(deeparg_LS)
SA_ox %>% count(deeparg_SR)


#-------------------------------------------------------------------------------
# plot ox peni phenotypes
pheno_ox_peni <- pheno_mic[,c(1,2,27,28,93,94)]

pheno_ox_peni_table <- pivot_longer(data = pheno_ox_peni, 
                                    cols = -c(1:2,5:6),
                                    names_to = "Antibiotic", 
                                    values_to = "Phenotype")

pheno_ox_peni_table$Phenotype <- factor(pheno_ox_peni_table$Phenotype, levels = c("S", "I", "R"))

# create heatmap
colors <- c("grey", "black")

# remove peni
pheno_ox_peni_table <-  pheno_ox_peni_table[c(which(pheno_ox_peni_table$Antibiotic == "Oxacillin")),]

pheno_ox_peni_table$Sample_id <- factor(pheno_ox_peni_table$Sample_id, levels = c(a, b, c)) 

heatmap_ox_peni <- ggplot(data = pheno_ox_peni_table, mapping = aes(x = Antibiotic,
                                                                    y = Sample_id,
                                                                    fill = Phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "Resistance interpretation") + 
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(MIC.Oxacillin), scales = "free_y", space = "free")

heatmap_ox_peni

heatmap_ox_peni <- ggplot(data = pheno_ox_peni_table, mapping = aes(x = "Oxacillin MIC",
                                                                    y = Sample_id,
                                                                    fill = MIC.Oxacillin)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") + 
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(MIC.Oxacillin), scales = "free_y", space = "free", switch = "y")

heatmap_ox_peni
