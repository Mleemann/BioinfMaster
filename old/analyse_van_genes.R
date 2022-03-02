library(ggplot2)
library("tidyr")

# read in samples to be analyzed
samples <- read.table('samples_EF.txt')

# get phenotype table
phenotypes <- read.csv2('phenotypes_vanco_EF.csv', header = TRUE, sep = ';')

# check sample ids; and correct them 
which(phenotypes$Sample_id %in% samples$V1 == FALSE)
# phenotypes$Sample_id[131] <- "correct value"
phenotypes$Sample_id[9] <- "110167-1-11"
phenotypes$Sample_id[72] <- "302293-18-wh"
phenotypes$Sample_id[99] <- "304084-1-20"
phenotypes$Sample_id[131] <- "607982-2-20"
phenotypes$Sample_id[136] <- "705943-2-21"
phenotypes$Sample_id[144] <- "802725-20-wh"

# include phenotype data into tabel
EF_abricate <- merge(phenotypes, samples, by.x = 'Sample_id', by.y = 'V1')

# get result tables
abr_ncbi <- read.table('van_abricate_ncbi.tab', sep = "\t")
abr_ncbi$V1 <- gsub(".fna","", abr_ncbi$V1)
if (length(which(abr_ncbi$V11 < 95)) != 0) {
abr_ncbi <- abr_ncbi[-c(which(abr_ncbi$V11 < 95)),]
}

abr_card <- read.table('van_abricate_card.tab', sep = "\t")
abr_card$V1 <- gsub(".fna","", abr_card$V1)
if (length(which(abr_card$V11 < 95)) != 0) {
abr_card <- abr_card[-c(which(abr_card$V11 < 95)),]
}

abr_resfinder <- read.table('van_abricate_resfinder.tab', sep = "\t")
abr_resfinder$V1 <- gsub(".fna","", abr_resfinder$V1)
if (length(which(abr_resfinder$V11 < 95)) != 0) {
abr_resfinder <- abr_resfinder[-c(which(abr_resfinder$V11 < 95)),]
}

abr_megares <- read.table('van_abricate_megares.tab', sep = "\t")
abr_megares$V1 <- gsub(".fna","", abr_megares$V1)
if (length(which(abr_megares$V11 < 95)) != 0) {
abr_megares <- abr_megares[-c(which(abr_megares$V11 < 95)),]
}

abr_argannot <- read.table('van_abricate_argannot.tab', sep = "\t")
abr_argannot$V1 <- gsub(".fna","", abr_argannot$V1)
if (length(which(abr_argannot$V11 < 95)) != 0) {
abr_argannot <- abr_argannot[-c(which(abr_argannot$V11 < 95)),]
}

amrf_nuc <- read.table('van_amrfinder_nuc.tab', sep = "\t")
amrf_nuc$V1 <- gsub("gnl\\|USB\\|","", amrf_nuc$V2)
amrf_nuc$V1 <- gsub('(_\\d*)',"", amrf_nuc$V1)
if (length(which(amrf_nuc$V11 < 95)) != 0) {
amrf_nuc <- amrf_nuc[-c(which(amrf_nuc$V17 < 95)),]
}

amrf_prot <- read.table('van_amrfinder_prot.tab', sep = "\t")
amrf_prot$V1 <- gsub('(_\\d*)',"", amrf_prot$V1)
if (length(which(amrf_prot$V11 < 95)) != 0) {
amrf_prot <- amrf_prot[-c(which(amrf_prot$V13 < 95)),]
}

rgi <- read.table(('van_rgi.tab'), sep = "\t")
rgi$V1 <- gsub("gnl\\|USB\\|","", rgi$V1)
rgi$V1 <- gsub('(_\\d*)',"", rgi$V1)
if (length(which(rgi$V11 < 95)) != 0) {
rgi <- rgi[-c(which(rgi$V11 < 95)),]
}

srax_basic <- read.table(('van_srax_basic.tab'), sep = "\t")
srax_basic$V2 <- gsub(".fna","", srax_basic$V2)
if (length(which(srax_basic$V11 < 95)) != 0) {
srax_basic <- srax_basic[-c(which(srax_basic$V11 < 95)),]
}

srax_ext <- read.table(('van_srax_ext.tab'), sep = "\t")
srax_ext$V2 <- gsub(".fna","", srax_ext$V2)
if (length(which(srax_ext$V11 < 95)) != 0) {
srax_ext <- srax_ext[-c(which(srax_ext$V11 < 95)),]
}

deeparg_LS <- read.table(('van_deeparg_LS.tab'), sep = "\t")
deeparg_LS$V4 <- gsub("gnl\\|USB\\|","", deeparg_LS$V4)
deeparg_LS$V4 <- gsub('(_\\d*)',"", deeparg_LS$V4)
for(i in 1:length(deeparg_LS[,1])) {
  genes <- unlist(strsplit(deeparg_LS[i,6], "\\|"))
  deeparg_LS[i,6] <- genes[length(genes)]
}
if (length(which(deeparg_LS$V8 < 95)) != 0) {
  deeparg_LS <- deeparg_LS[-c(which(deeparg_LS$V8 < 95)),]
}

deeparg_SR <- read.table(('van_deeparg_SR.tab'), sep = "\t")
if (length(which(deeparg_SR$V8 < 95)) != 0) {
  deeparg_SR <- deeparg_SR[-c(which(deeparg_SR$V8 < 95)),]
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
deeparg_SR_unique <- deeparg_SR_unique[-c(which(deeparg_SR_unique$V6 == "undefined")), ]

resfinder_as <- read.table(('van_resfinder_assembly.tab'), sep = "\t")
if (length(which(resfinder_as$V2 < 95)) != 0) {
  resfinder_as <- resfinder_as[-c(which(resfinder_as$V2 < 95)),]
}
resfinder_as$V6 <- gsub('gnl\\|USB\\|','', resfinder_as$V6)
resfinder_as$V6 <- gsub('(_\\d*)','', resfinder_as$V6)

resfinder_re <- read.table(('van_resfinder_reads.tab'), sep = "\t")
if (length(which(resfinder_re$V2 < 95)) != 0) {
  resfinder_re <- resfinder_re[-c(which(resfinder_re$V2 < 95)),]
}


# check if van found for sample_id and include it into the table --> TRUE van present
EF_abricate$van_ncbi <- EF_abricate$Sample_id %in% abr_ncbi$V1
EF_abricate$van_card <- EF_abricate$Sample_id %in% abr_card$V1
EF_abricate$van_resfinder <- EF_abricate$Sample_id %in% abr_resfinder$V1
EF_abricate$van_megares <- EF_abricate$Sample_id %in% abr_megares$V1
EF_abricate$van_argannot <- EF_abricate$Sample_id %in% abr_argannot$V1
EF_abricate$van_amrf_nuc <- EF_abricate$Sample_id %in% amrf_nuc$V1
EF_abricate$van_amrf_prot <- EF_abricate$Sample_id %in% amrf_prot$V1
EF_abricate$van_rgi <- EF_abricate$Sample_id %in% rgi$V1 
EF_abricate$van_srax_b <- EF_abricate$Sample_id %in% srax_basic$V2
EF_abricate$van_srax_e <- EF_abricate$Sample_id %in% srax_ext$V2
EF_abricate$van_deeparg_LS <- EF_abricate$Sample_id %in% deeparg_LS$V4
EF_abricate$van_deeparg_SR_unique <- EF_abricate$Sample_id %in% deeparg_SR_unique$V1
EF_abricate$van_resfinder_as <- EF_abricate$Sample_id %in% resfinder_as$V6
EF_abricate$van_resfinder_re <- EF_abricate$Sample_id %in% resfinder_re$V10

# creating a heatmap 
for_heatmap <- EF_abricate
for_heatmap[for_heatmap == "c(R, R)"] <- "R"
for_heatmap[for_heatmap == "c(S, S)"] <- "S"
for_heatmap[for_heatmap == "c(U, R)"] <- "R"
for_heatmap[for_heatmap == "c(R, U)"] <- "R"
for_heatmap[for_heatmap == "c(R, S)"] <- "R"
for_heatmap[for_heatmap == "c(S, R)"] <- "R"

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
                          cols = -c(1:3),
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

heatmap_table$Tool <- ifelse(heatmap_table$ToolDB == "van_rgi", "rgi",
                             ifelse(heatmap_table$ToolDB == "assembly" | heatmap_table$ToolDB == "reads", "resfinder",
                                ifelse(heatmap_table$ToolDB == "LS" | heatmap_table$ToolDB == "SR", "deeparg",
                                  ifelse(heatmap_table$ToolDB == "basic" | heatmap_table$ToolDB == "ext", "sraX",
                                      ifelse(heatmap_table$ToolDB == "nuc" | heatmap_table$ToolDB == "prot", "amrfinder", "abricate")))))

heatmap_table[heatmap_table == "van_rgi"] <- "card"

# create heatmap
colors <- c("darkgreen", "green", "orange", "red")

EF_heatmap <- ggplot(data = heatmap_table, mapping = aes(x = ToolDB,
                                                       y = Sample_id,
                                                       fill = Concordance)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_manual(values = colors, labels = c("R/R", "S/S", "S/R", "R/S")) +
  labs(x = "", y = "") +
  labs(fill = "Phenotype/Genotype") +
  theme(axis.text.y = element_text(size = 6)) +
  #coord_fixed(ratio = 0.18) +

  facet_grid(~ Tool, scales = "free_x", space = "free_x")
  
EF_heatmap


test_genes <- abr_ncbi[, c(1,6)]
gen_ord <- c('vanA', 'vanH-A', 'vanR-A', 'vanS-A', 'vanX-A', 'vanY-A', 'vanZ-A', 'vanB', 'vanH-B', 'vanR-B', 'vanS-B', 'vanW-B', 'vanX-B', 'vanY-B')
for(i in gen_ord) {
  test_genes$V3[test_genes$V6 == i] <- which(gen_ord == i)
}

gene_order <- factor(test_genes$V6, level = c('vanA', 'vanH-A', 'vanR-A', 'vanS-A', 'vanX-A', 'vanY-A', 'vanZ-A', 'vanB', 'vanH-B', 'vanR-B', 'vanS-B', 'vanW-B', 'vanX-B', 'vanY-B'))

colors2 <- c("blue4", "red2", "red", "orangered", "sienna1", "darkorange", "blue3", "blue1", "dodgerblue2", "dodgerblue", "deepskyblue2", "deepskyblue", "red4", "red3")
test_genes_heatmap <- ggplot(data = test_genes, mapping = aes(x = gene_order,
                                                              y = V1,
                                                              fill = V3)) +
  geom_tile(colour="white", size=1) +
  labs(x = "", y = "") +
  scale_fill_manual(values = colors2) +
  theme(legend.position = 'none', axis.text.y = element_text(size = 6)) + 
  ggtitle("Abricate - NCBI")
test_genes_heatmap

