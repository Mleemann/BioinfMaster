#install.packages('tidyverse')
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)

setwd("~/data_project/assemblies/EF")

# read in samples to be analyzed
samples <- read.table('samples_EF.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_van_tei_EF.csv', header = TRUE, sep = ';')
vre <- which(phenotypes$Sample_id == "705423-17")
phenotypes$LGM_BK[vre] <- "entfaemv"

# check sample ids; and correct them 
which(phenotypes$Sample_id %in% samples$V1 == FALSE)

phenotypes$Sample_id[9] <- "110167-1-11"
phenotypes$Sample_id[72] <- "302293-18-wh"
phenotypes$Sample_id[129] <- "607982-2-20"
phenotypes$Sample_id[134] <- "705943-2-21"
phenotypes$Sample_id[142] <- "802725-20-wh"

mhk <- which(phenotypes$Sample_id == "302613-17")
phenotypes$MIC.Vancomycin[mhk] <- "= 24"

# get rid of multiple values
phenotypes[phenotypes == "c(R, R)"] <- "R"
phenotypes[phenotypes == "c(S, S)"] <- "S"
phenotypes[phenotypes == "c(U, R)"] <- "R"
phenotypes[phenotypes == "c(R, U)"] <- "R"
phenotypes[phenotypes == "c(R, S)"] <- "R"
phenotypes[phenotypes == "c(S, R)"] <- "R"
phenotypes[phenotypes == "c(U, S)"] <- "S"

phenotypes[phenotypes == "c(>= 32, = 0)" | phenotypes == "c(= 4, >= 32)" | phenotypes == "c(= 0, >= 32)" | 
                phenotypes == "c(>= 32, = 4)" | phenotypes == "c(= 32, = 32)" | phenotypes == "c(> 256, >= 32)" |
                phenotypes == "c(> 256, > 32)"] <- ">= 32"
phenotypes[phenotypes == "c(> 256, > 256)"] <- ">= 32"
phenotypes[phenotypes == "c(= 0, = 1.5)"] <- "= 1.5"
phenotypes[phenotypes == "c(= 0, = 12)"] <- "= 12"
phenotypes[phenotypes == "c(= 8, = 4)" | phenotypes == "c(= 1.5, = 8)"] <- "= 8"
phenotypes[phenotypes == "c(= 2, = 1)" | phenotypes == "<= 2"] <- "= 2"
phenotypes[phenotypes == "c(<= 0.5, = 1)" | phenotypes == "c(<= 0.5, <= 0.5)"] <- "<= 0.5"
phenotypes[phenotypes == "n.a. 20002"] <- "= 0"


#-------------------------------------------------------------------------------
setwd("~/data_project/assemblies/EF/amr_results/van")

# get result tables
# ncbi
abr_ncbi <- read.table('van_abricate_ncbi.tab', sep = "\t")
if (length(which(abr_ncbi$V11 < 95)) != 0) {
  abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$V11 < 95)),]
}
abr_ncbi$V1 <- gsub(".fna","", abr_ncbi$V1)
names(abr_ncbi)[1] <- c('Sample_id')

abr_card <- read.table('van_abricate_card.tab', sep = "\t")
if (length(which(abr_card$V11 < 95)) != 0) {
  abr_card <- abr_card[-c(which(abr_card$V11 < 95)),]
}
abr_card$V1 <- gsub(".fna","", abr_card$V1)
names(abr_card)[1] <- c('Sample_id')

abr_resfinder <- read.table('van_abricate_resfinder.tab', sep = "\t")
if (length(which(abr_resfinder$V11 < 95)) != 0) {
  abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$V11 < 95)),]
}
abr_resfinder$V1 <- gsub(".fna","", abr_resfinder$V1)
names(abr_resfinder)[1] <- c('Sample_id')

abr_megares <- read.table('van_abricate_megares.tab', sep = "\t")
if (length(which(abr_megares$V11 < 95)) != 0) {
  abr_megares <- abr_megares[-c(which(abr_megares$V11 < 95)),]
}
abr_megares$V1 <- gsub(".fna","", abr_megares$V1)
names(abr_megares)[1] <- c('Sample_id')

abr_argannot <- read.table('van_abricate_argannot.tab', sep = "\t")
if (length(which(abr_argannot$V11 < 95)) != 0) {
  abr_argannot <- abr_argannot[-c(which(abr_argannot$V11 < 95)),]
}
abr_argannot$V1 <- gsub(".fna","", abr_argannot$V1)
names(abr_argannot)[1] <- c('Sample_id')

#amrfinder
amrf_nuc <- read.table('van_amrfinder_nuc.tab', sep = "\t")
if (length(which(amrf_nuc$V11 < 95)) != 0) {
  amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$V17 < 95)),]
}
amrf_nuc$V1 <- gsub("gnl\\|USB\\|","", amrf_nuc$V2)
amrf_nuc$V1 <- gsub('(_\\d*)',"", amrf_nuc$V1)
names(amrf_nuc)[1] <- 'Sample_id'

amrf_prot <- read.table('van_amrfinder_prot.tab', sep = "\t")
if (length(which(amrf_prot$V11 < 95)) != 0) {
  amrf_prot <- amrf_prot[-c(which(amrf_prot$V13 < 95)),]
}
amrf_prot$V1 <- gsub('(_\\d*)',"", amrf_prot$V1)
names(amrf_prot)[1] <- 'Sample_id'


# rgi
rgi <- read.table(('van_rgi.tab'), sep = "\t")
if (length(which(rgi$V11 < 95)) != 0) {
  rgi <- rgi[-c(which(rgi$V11 < 95)),]
}
rgi$V1 <- gsub("gnl\\|USB\\|","", rgi$V1)
rgi$V1 <- gsub('(_\\d*)',"", rgi$V1)
names(rgi)[1] <- 'Sample_id'


# sraX
srax_basic <- read.table(('van_srax_basic.tab'), sep = "\t")
if (length(which(srax_basic$V11 < 95)) != 0) {
  srax_basic <- srax_basic[-c(which(srax_basic$V11 < 95)),]
}
srax_basic$V2 <- gsub(".fna","", srax_basic$V2)
names(srax_basic)[2] <- 'Sample_id'

srax_ext <- read.table(('van_srax_ext.tab'), sep = "\t")
if (length(which(srax_ext$V11 < 95)) != 0) {
  srax_ext <- srax_ext[-c(which(srax_ext$V11 < 95)),]
}
srax_ext$V2 <- gsub(".fna","", srax_ext$V2)
names(srax_ext)[2] <- 'Sample_id'


# Deeparg
deeparg_LS <- read.table(('van_deepargLS_noPot.tab'), sep = "\t")
if (length(which(deeparg_LS$V8 < 95)) != 0) {
  deeparg_LS <- deeparg_LS[-c(which(deeparg_LS$V8 < 95)),]
}
deeparg_LS$V4 <- gsub("gnl\\|USB\\|","", deeparg_LS$V4)
deeparg_LS$V4 <- gsub('(_\\d*)',"", deeparg_LS$V4)
for(i in 1:length(deeparg_LS[,1])) {
  genes <- unlist(strsplit(deeparg_LS[i,6], "\\|"))
  deeparg_LS[i,6] <- genes[length(genes)]
}
names(deeparg_LS)[4] <- 'Sample_id'

deeparg_SR <- read.table(('van_deepargSR_noPot.tab'), sep = "\t")
if (length(which(deeparg_SR$V8 < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$V8 < 95)),]
}
if (length(which(deeparg_SR$V7 < 0.8)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$V7 < 0.8)),]
}
deeparg_SR_unique <- deeparg_SR %>% distinct(V1, V6, .keep_all = TRUE)
deeparg_SR_unique$V1 <- gsub("entfaem/resistance/deeparg_SR/","", deeparg_SR_unique$V1)
deeparg_SR_unique$V1 <- gsub("entfaemv/resistance/deeparg_SR/","", deeparg_SR_unique$V1)
for(i in 1:length(deeparg_SR_unique[,1])) {
  genes <- unlist(strsplit(deeparg_SR_unique[i,1], "\\_"))
  deeparg_SR_unique[i,1] <- genes[1]
}
for(i in 1:length(deeparg_SR_unique[,1])) {
  genes <- unlist(strsplit(deeparg_SR_unique[i,6], "\\|"))
  deeparg_SR_unique[i,6] <- genes[length(genes)]
}
#deeparg_SR_unique <- deeparg_SR_unique[-c(which(deeparg_SR_unique$V6 == "undefined")), ]
names(deeparg_SR_unique)[1] <- 'Sample_id'
rm(deeparg_SR)


# resfinder
resfinder_as <- read.table(('van_resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$V2 < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$V2 < 95)),]
}
resfinder_as$V6 <- gsub('gnl\\|USB\\|','', resfinder_as$V6)
resfinder_as$V6 <- gsub('(_\\d*)','', resfinder_as$V6)
names(resfinder_as)[6] <- 'Sample_id'

resfinder_re <- read.table(('van_resfinder_reads.tab'), sep = "\t")
if (length(which(resfinder_re$V2 < 95)) != 0) {
  resfinder_re <- resfinder_re[-c(which(resfinder_re$V2 < 95)),]
}
names(resfinder_re)[10] <- 'Sample_id'


#-------------------------------------------------------------------------------
# combine all van genes found in one table

names(abr_ncbi)[6] <- 'van_genes'
abr_ncbi$Tool <- 'abr_ncbi'
names(abr_card)[6] <- 'van_genes'
abr_card$Tool <- 'abr_card'
names(abr_megares)[6] <- 'van_genes'
abr_megares$Tool <- 'abr_megares'
names(abr_resfinder)[14] <- 'van_genes'
abr_resfinder$Tool <- 'abr_resfinder'
names(abr_argannot)[6] <- 'van_genes'
abr_argannot$Tool <- 'abr_argannot'
names(amrf_nuc)[6] <- 'van_genes'
amrf_nuc$Tool <- 'amrf_nuc'
names(amrf_prot)[2] <- 'van_genes'
amrf_prot$Tool <- 'amrf_prot'
names(deeparg_LS)[6] <- 'van_genes'
deeparg_LS$Tool <- 'deeparg_LS'
names(deeparg_SR_unique)[6] <- 'van_genes'
deeparg_SR_unique$Tool <- 'deeparg_SR'
names(resfinder_as)[1] <- 'van_genes'
resfinder_as$Tool <- 'resfinder_as'
names(resfinder_re)[1] <- 'van_genes'
resfinder_re$Tool <- 'resfinder_re'
names(rgi)[2] <- 'van_genes'
rgi$Tool <- 'rgi'
names(srax_basic)[6] <- 'van_genes'
srax_basic$Tool <- 'srax_basic'
names(srax_ext)[6] <- 'van_genes'
srax_ext$Tool <- 'srax_ext'

all_van <- rbind(abr_ncbi[,c(1,6,16)], abr_card[,c(1,6,16)], abr_megares[,c(1,6,16)], abr_resfinder[,c(1,14,16)], abr_argannot[,c(1,6,16)],
                 amrf_nuc[,c(1,6,23)], amrf_prot[,c(1,2,19)], deeparg_LS[,c(4,6,13)], deeparg_SR_unique[,c(1,6,13)],
                 resfinder_as[,c(6,1,10)], resfinder_re[,c(10,1,11)], rgi[,c(1,2,9)], srax_basic[,c(2,6,16)], srax_ext[,c(2,6,16)])

all_van <- unique(all_van)

all_van$van_genes <- gsub("\\(Gly\\)", "", all_van$van_genes)
all_van$van_genes <- gsub("\\-", "", all_van$van_genes)


all_van[all_van == "VANRA"] <- "vanRA"
all_van[all_van == "VANHA"] <- "vanHA"
all_van[all_van == "VANRB"] <- "vanRB"
all_van[all_van == "VANSB"] <- "vanSB"
all_van[all_van == "VANYB"] <- "vanYB"
all_van[all_van == "VANWB"] <- "vanWB"
all_van[all_van == "VANHB"] <- "vanHB"
all_van[all_van == "vanAA"] <- "vanA"
all_van[all_van == "vanAB"] <- "vanB"


#-------------------------------------------------------------------------------
# plot van genes found in all tools

gene_clusters <- pivot_wider(all_van[,c(1,2,3)], names_from = Tool, values_from = van_genes)

# add samples that have no van genes
for (i in samples$V1[which(samples$V1 %in% gene_clusters$Sample_id == "FALSE")]) {
  j <- length(gene_clusters$Sample_id) + 1
  gene_clusters[j,1] <- i
}

gene_clusters <- gene_clusters[-c(which(gene_clusters$Sample_id =="304084-1-20" | gene_clusters$Sample_id =="305162-21" | gene_clusters$Sample_id =="716692-17")),]


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

for (i in 1:length(gene_clusters$srax_basic)) {
  gene_clusters$srax_basic[[i]] <- sort(gene_clusters$srax_basic[[i]])
}

for (i in 1:length(gene_clusters$srax_ext)) {
  gene_clusters$srax_ext[[i]] <- sort(gene_clusters$srax_ext[[i]])
}



gene_cluster_table <- pivot_longer(data = gene_clusters, 
                                   cols = -c(1),
                                   names_to = "ToolDB", 
                                   values_to = "gene_cluster")

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

# include clinical phenotype 
gene_cluster_table$Type <- 0
for (i in 1:length(gene_cluster_table$Sample_id)){
  gene_cluster_table$Type[i] <- phenotypes$LGM_BK[which(gene_cluster_table$Sample_id[i] == phenotypes$Sample_id)]
}

# complete vanA operon: cluster_rep = 1, partially cluster_rep = 2, partially with ZA cluster_rep = 6
# complete vanB operon: cluster_rep = 3, partially cluster_rep = 4
# undefined = 5; no real determinant = 0 (this includes singel vanS, vanR)

rep_one <- c(14,23,25,27,35,42,49,59)
rep_two <- c(15,16,19,20,22,24,36,38)
rep_six <- c(21)
rep_three <- c(1,6,13,32,40,52,53)
rep_four <- c(2,3,4,7,8,9,10,11,12,28,29,30,31,33,34,43,44)
rep_five <- c(18,26,37,39,41,45,46,47,48,50,51,54,55,56,57,58)
rep_seven <- c(61,17,5)

gene_cluster_table$cluster_rep <- ifelse(gene_cluster_table$cluster_nr %in% rep_one, gene_cluster_table$cluster_rep <- 1,
                                         ifelse(gene_cluster_table$cluster_nr %in% rep_two, gene_cluster_table$cluster_rep <- 2,
                                                ifelse(gene_cluster_table$cluster_nr %in% rep_three, gene_cluster_table$cluster_rep <- 3, 
                                                       ifelse(gene_cluster_table$cluster_nr %in% rep_four, gene_cluster_table$cluster_rep <- 4,
                                                              ifelse(gene_cluster_table$cluster_nr %in% rep_five, gene_cluster_table$cluster_rep <- 5, 
                                                                     ifelse(gene_cluster_table$cluster_nr %in% rep_six, gene_cluster_table$cluster_rep <- 6, 
                                                                            ifelse(gene_cluster_table$cluster_nr %in% rep_seven, gene_cluster_table$cluster_rep <- 7, gene_cluster_table$cluster_rep <- 0))))))) 

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


for (i in 1:length(gene_cluster_table$Sample_id)){
  gene_cluster_table$MIC[i] <- phenotypes$MIC.Vancomycin[which(gene_cluster_table$Sample_id[i] == phenotypes$Sample_id)]
}



# order legend
gene_cluster_table$cluster_rep <- factor(gene_cluster_table$cluster_rep, levels = c(1,6,2,3,4,5,7,0)) 

a <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "1"])) 
b <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "2"]))
c <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "6"]))
d <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "3"]))
e <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "4"]))
f <- sort(unique(gene_cluster_table$Sample_id[gene_cluster_table$ToolDB == "argannot" & gene_cluster_table$cluster_rep == "0"]))

gene_cluster_table$Sample_id <- factor(gene_cluster_table$Sample_id, levels = c(a, b, c, d, e, f)) 
colors <- c("blue4", "blue1", "dodgerblue2", "orangered", "sienna1", "plum1", "yellowgreen", "moccasin")

gene_cluster_heatmap <- ggplot(data = gene_cluster_table, mapping = aes(x = ToolDB,
                                                                        y = Sample_id,
                                                                        fill = cluster_rep)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("complete vanA operon", "partial vanA operon incl. vanZA", "partial vanA operon", "complete vanB operon", "partial vanB operon", "mixed", "only regulator" , "none")) +
  labs(x = "", y = "") +
  labs(fill = "Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(gene_cluster_table$MIC, levels = c("= 0", "<= 0.5", "= 1", "= 2",  "= 8", "= 12", "= 24",">= 32", "> 256" ))), vars(Tool), scales = "free", space = "free")

gene_cluster_heatmap



#-------------------------------------------------------------------------------
# create concordance table

EF_van_all <- phenotypes

# check if any van gene found for sample_id and include it into the table --> TRUE van present
EF_van_all$van_ncbi <- EF_van_all$Sample_id %in% abr_ncbi$Sample_id
EF_van_all$van_card <- EF_van_all$Sample_id %in% abr_card$Sample_id
EF_van_all$van_resfinder <- EF_van_all$Sample_id %in% abr_resfinder$Sample_id
EF_van_all$van_megares <- EF_van_all$Sample_id %in% abr_megares$Sample_id
EF_van_all$van_argannot <- EF_van_all$Sample_id %in% abr_argannot$Sample_id
EF_van_all$van_amrf_nuc <- EF_van_all$Sample_id %in% amrf_nuc$Sample_id
EF_van_all$van_amrf_prot <- EF_van_all$Sample_id %in% amrf_prot$Sample_id
EF_van_all$van_rgi <- EF_van_all$Sample_id %in% rgi$Sample_id 
EF_van_all$van_srax_b <- EF_van_all$Sample_id %in% srax_basic$Sample_id
EF_van_all$van_srax_e <- EF_van_all$Sample_id %in% srax_ext$Sample_id
EF_van_all$van_deeparg_LS <- EF_van_all$Sample_id %in% deeparg_LS$Sample_id
EF_van_all$van_deeparg_SR_unique <- EF_van_all$Sample_id %in% deeparg_SR_unique$Sample_id
EF_van_all$van_resfinder_as <- EF_van_all$Sample_id %in% resfinder_as$Sample_id
EF_van_all$van_resfinder_re <- EF_van_all$Sample_id %in% resfinder_re$Sample_id



EF_van_single <- phenotypes

# check if vanA or vanB found for sample_id and include it into the table --> TRUE van present
EF_van_single$van_ncbi <- EF_van_single$Sample_id %in% abr_ncbi$Sample_id[c(which(abr_ncbi$van_genes == "vanA" | abr_ncbi$van_genes == "vanB"))]
EF_van_single$van_card <- EF_van_single$Sample_id %in% abr_card$Sample_id[c(which(abr_card$van_genes == "vanA" | abr_card$van_genes == "vanB"))]
EF_van_single$van_resfinder <- EF_van_single$Sample_id %in% abr_resfinder$Sample_id[c(which(abr_resfinder$van_genes == "VanHAX" | abr_resfinder$van_genes == "VanHBX"))]
EF_van_single$van_megares <- EF_van_single$Sample_id %in% abr_megares$Sample_id[c(which(abr_megares$van_genes == "VANHA" | abr_megares$van_genes == "VANHB"))]
EF_van_single$van_argannot <- EF_van_single$Sample_id %in% abr_argannot$Sample_id[c(which(abr_argannot$van_genes == "(Gly)vanA-A" | abr_argannot$van_genes == "(Gly)vanB" | abr_argannot$van_genes == "(Gly)vanA-B"))]
EF_van_single$van_amrf_nuc <- EF_van_single$Sample_id %in% amrf_nuc$Sample_id[c(which(amrf_nuc$van_genes == "vanA" | amrf_nuc$van_genes == "vanB"))]
EF_van_single$van_amrf_prot <- EF_van_single$Sample_id %in% amrf_prot$Sample_id[c(which(amrf_prot$van_genes == "vanA" | amrf_prot$van_genes == "vanB"))]
EF_van_single$van_rgi <- EF_van_single$Sample_id %in% rgi$Sample_id [c(which(rgi$van_genes == "vanA" | rgi$van_genes == "vanB"))]
EF_van_single$van_srax_b <- EF_van_single$Sample_id %in% srax_basic$Sample_id[c(which(srax_basic$van_genes == "vanA" | srax_basic$van_genes == "vanB"))]
EF_van_single$van_srax_e <- EF_van_single$Sample_id %in% srax_ext$Sample_id[c(which(srax_ext$van_genes == "vanA" | srax_ext$van_genes == "vanB"))]
EF_van_single$van_deeparg_LS <- EF_van_single$Sample_id %in% deeparg_LS$Sample_id[c(which(deeparg_LS$van_genes == "vanA" | deeparg_LS$van_genes == "vanB"))]
EF_van_single$van_deeparg_SR_unique <- EF_van_single$Sample_id %in% deeparg_SR_unique$Sample_id[c(which(deeparg_SR_unique$van_genes == "vanA" | deeparg_SR_unique$van_genes == "vanB"))]
EF_van_single$van_resfinder_as <- EF_van_single$Sample_id %in% resfinder_as$Sample_id[c(which(resfinder_as$van_genes == "VanHAX" | resfinder_as$van_genes == "VanHBX"))]
EF_van_single$van_resfinder_re <- EF_van_single$Sample_id %in% resfinder_re$Sample_id[c(which(resfinder_re$van_genes == "VanHAX" | resfinder_re$van_genes == "VanHBX"))]


# include also other van types
EF_all_VRE <- EF_van_single

EF_all_VRE$van_deeparg_SR_unique <- EF_all_VRE$Sample_id %in% deeparg_SR_unique$Sample_id[c(which(deeparg_SR_unique$van_genes == "vanA" | deeparg_SR_unique$van_genes == "vanB" | 
                                                                                                    deeparg_SR_unique$van_genes == "vanC" | deeparg_SR_unique$van_genes == "vanD" | 
                                                                                                    deeparg_SR_unique$van_genes == "vanE" | deeparg_SR_unique$van_genes == "vanG" | 
                                                                                                    deeparg_SR_unique$van_genes == "vanL" | deeparg_SR_unique$van_genes == "vanM" | 
                                                                                                    deeparg_SR_unique$van_genes == "vanN"))]



for_heatmap <- EF_all_VRE
for_heatmap <- for_heatmap[-c(which(for_heatmap$Sample_id =="304084-1-20" | for_heatmap$Sample_id =="305162-21" | for_heatmap$Sample_id =="716692-17")),]


# Concordance values: R/R = 1, S/S = 2, S/R = 3, R/S = 4
for_heatmap$van_ncbi <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_ncbi == "FALSE", 2,
                               ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_ncbi == "TRUE", 3,
                                      ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_ncbi == "TRUE", 1, 4)))

for_heatmap$van_card <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_card == "FALSE", 2,
                               ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_card == "TRUE", 3,
                                      ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_card == "TRUE", 1, 4)))

for_heatmap$van_resfinder <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_resfinder == "FALSE", 2,
                                    ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_resfinder == "TRUE", 3,
                                           ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_resfinder == "TRUE", 1, 4)))

for_heatmap$van_megares <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_megares == "FALSE", 2,
                                  ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_megares == "TRUE", 3,
                                         ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_megares == "TRUE", 1, 4)))

for_heatmap$van_argannot <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_argannot == "FALSE", 2,
                                   ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_argannot == "TRUE", 3,
                                          ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_argannot == "TRUE", 1, 4)))

for_heatmap$van_amrf_nuc <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_amrf_nuc == "FALSE", 2,
                                   ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_amrf_nuc == "TRUE", 3,
                                          ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_amrf_nuc == "TRUE", 1, 4)))

for_heatmap$van_amrf_prot <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_amrf_prot == "FALSE", 2,
                                    ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_amrf_prot == "TRUE", 3,
                                           ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_amrf_prot == "TRUE", 1, 4)))

for_heatmap$van_rgi <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_rgi == "FALSE", 2,
                              ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_rgi == "TRUE", 3,
                                     ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_rgi == "TRUE", 1, 4)))

for_heatmap$van_srax_b <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_srax_b == "FALSE", 2,
                                 ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_srax_b == "TRUE", 3,
                                        ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_srax_b == "TRUE", 1, 4)))

for_heatmap$van_srax_e <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_srax_e == "FALSE", 2,
                                 ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_srax_e == "TRUE", 3,
                                        ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_srax_e == "TRUE", 1, 4)))

for_heatmap$van_deeparg_LS <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_deeparg_LS == "FALSE", 2,
                                     ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_deeparg_LS == "TRUE", 3,
                                            ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_deeparg_LS == "TRUE", 1, 4)))

for_heatmap$van_deeparg_SR_unique <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_deeparg_SR_unique == "FALSE", 2,
                                            ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_deeparg_SR_unique == "TRUE", 3,
                                                   ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_deeparg_SR_unique == "TRUE", 1, 4)))

for_heatmap$van_resfinder_as <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_resfinder_as == "FALSE", 2,
                                       ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_resfinder_as == "TRUE", 3,
                                              ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_resfinder_as == "TRUE", 1, 4)))

for_heatmap$van_resfinder_re <- ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_resfinder_re == "FALSE", 2,
                                       ifelse(for_heatmap$Vancomycin == "S" & for_heatmap$van_resfinder_re == "TRUE", 3,
                                              ifelse(for_heatmap$Vancomycin == "R" & for_heatmap$van_resfinder_re == "TRUE", 1, 4)))


# transform data for heatmap
heatmap_table <- pivot_longer(data = for_heatmap, 
                              cols = -c(1:6),
                              names_to = "ToolDB", 
                              values_to = "Concordance")

heatmap_table$Concordance <- as.character(heatmap_table$Concordance)
heatmap_table[heatmap_table == "van_ncbi"] <- "ncbi"
heatmap_table[heatmap_table == "van_card"] <- "card"
heatmap_table[heatmap_table == "van_resfinder"] <- "resfinder"
heatmap_table[heatmap_table == "van_megares"] <- "megares"
heatmap_table[heatmap_table == "van_argannot"] <- "argannot"
heatmap_table[heatmap_table == "van_amrf_nuc"] <- "nuc"
heatmap_table[heatmap_table == "van_amrf_prot"] <- "prot"
heatmap_table[heatmap_table == "van_srax_b"] <- "basic"
heatmap_table[heatmap_table == "van_srax_e"] <- "ext"
heatmap_table[heatmap_table == "van_deeparg_LS"] <- "LS"
heatmap_table[heatmap_table == "van_deeparg_SR_unique"] <- "SR"
heatmap_table[heatmap_table == "van_resfinder_as"] <- "assembly"
heatmap_table[heatmap_table == "van_resfinder_re"] <- "reads"

# include column for the different tools
heatmap_table$Tool <- ifelse(heatmap_table$ToolDB == "van_rgi", "RGI",
                             ifelse(heatmap_table$ToolDB == "assembly" | heatmap_table$ToolDB == "reads", "ResFinder",
                                    ifelse(heatmap_table$ToolDB == "LS" | heatmap_table$ToolDB == "SR", "DeepARG",
                                           ifelse(heatmap_table$ToolDB == "basic" | heatmap_table$ToolDB == "ext", "sraX",
                                                  ifelse(heatmap_table$ToolDB == "nuc" | heatmap_table$ToolDB == "prot", "AMRFinder", "ABRicate")))))

# do it after step before otherwise card of rgi is counted as abricate
#heatmap_table[heatmap_table == "van_rgi"] <- "card"


# create heatmap
colors <- c("darkgreen", "green", "red")

EF_heatmap <- ggplot(data = heatmap_table, mapping = aes(x = ToolDB,
                                                         y = Sample_id,
                                                         fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype/Genotype") +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, angle = 90, hjust=1, vjust=0.5),
        strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6), panel.background = element_rect(color = "white", linetype = "solid")) +
  facet_grid(vars(LGM_BK), vars(Tool), scales = "free", space = "free")

EF_heatmap


heatmap_table$Sample_id <- factor(heatmap_table$Sample_id, levels = c(a, b, c, d, e, f)) 

EF_heatmap <- ggplot(data = heatmap_table, mapping = aes(x = ToolDB,
                                                         y = Sample_id,
                                                         fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype/Genotype") +
  theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 10, face = "bold"), axis.ticks.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"), plot.background = element_rect(fill = 'white')) +
  facet_grid(vars(factor(heatmap_table$MIC.Vancomycin, levels = c("= 0", "<= 0.5", "= 1", "= 2",  "= 8", "= 12", "= 24", ">= 32", "> 256" ))), 
             vars(Tool), scales = "free", space = "free")

EF_heatmap


#-------------------------------------------------------------------------------
# vanco teico phenotype

van_tei_table <- pivot_longer(data = heatmap_table, 
                              cols = -c(1:2, 5:9),
                              names_to = "Antibiotic", 
                              values_to = "phenotype")

colors <- c("black", "grey")

van_tei_table$Sample_id <- factor(van_tei_table$Sample_id, levels = c(a, b, c, d, e, f)) 

van_tei_heatmap <- ggplot(data = van_tei_table, mapping = aes(x = Antibiotic,
                                                              y = Sample_id,
                                                              fill = phenotype)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(van_tei_table$MIC.Vancomycin, levels = c("= 0", "<= 0.5", "= 1", "= 2",  "= 8", "= 12", "= 24", ">= 32", "> 256" ))), scales = "free", space = "free")

van_tei_heatmap

van_tei_table$MIC.Vancomycin <- factor(van_tei_table$MIC.Vancomycin, levels = c("= 0", "<= 0.5", "= 1", "= 2",  "= 8", "= 12", "= 24", ">= 32", "> 256" ))


van_tei_heatmap2 <- ggplot(data = van_tei_table, mapping = aes(x = "Vancomycin MIC",
                                                               y = Sample_id,
                                                               fill = MIC.Vancomycin)) +
  geom_tile(colour="white", size=0.5) +
  #scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  labs(fill = "MIC") +
  theme(strip.text.y = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust=1, vjust=0.5), 
        strip.background = element_blank()) +
  facet_grid(vars(factor(van_tei_table$MIC.Vancomycin, levels = c("= 0", "<= 0.5", "= 1", "= 2",  "= 8", "= 12", "= 24", ">= 32", "> 256" ))), scales = "free", space = "free")

van_tei_heatmap2

