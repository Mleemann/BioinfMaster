library(dplyr)

# get abricate ncbi data for EC
EC_abricate_ncbi <- read.delim('~/data_project/assemblies/EC/abricate_ncbi.tab', header = FALSE, sep = "\t")
EC_abricate_ncbi$V1 <- gsub(".fna","", EC_abricate_ncbi$V1)
if (length(which(EC_abricate_ncbi$V11 < 95)) != 0) {
  EC_abricate_ncbi <- EC_abricate_ncbi[-c(which(EC_abricate_ncbi$V11 < 95)),]
}

# extract genes with resistance CARBAPENEM, BETA-LACTAM, or CEPHALOSPORIN
EC_ncbi_cephBlaCarb <-EC_abricate_ncbi[c(which(EC_abricate_ncbi$V15 == "CARBAPENEM" | EC_abricate_ncbi$V15 == "BETA-LACTAM" | EC_abricate_ncbi$V15 == "CEPHALOSPORIN")), ]

# get names of genes and its resistance 
resistance_EC_ncbi <- EC_ncbi_cephBlaCarb[,c(6,15)] %>% distinct(V6,.keep_all=TRUE)

# count occurence of the genes 
EC_ncbi_occ <- EC_ncbi_cephBlaCarb %>% count(V6)

# include resistance into occurence table and sort genes according to decreasing count
EC_ncbi_occ <- merge(resistance_EC_ncbi, EC_ncbi_occ, by.x = "V6", by.y = "V6")
EC_ncbi_occ <- EC_ncbi_occ[order(EC_ncbi_occ$n,decreasing = TRUE), ]


# same for KP
KP_abricate_ncbi <- read.delim('~/data_project/assemblies/KP/abricate_ncbi.tab', header = FALSE, sep = "\t")
KP_abricate_ncbi$V1 <- gsub(".fna","", KP_abricate_ncbi$V1)
if (length(which(KP_abricate_ncbi$V11 < 95)) != 0) {
  KP_abricate_ncbi <- KP_abricate_ncbi[-c(which(KP_abricate_ncbi$V11 < 95)),]
}

KP_ncbi_cephBlaCarb <-KP_abricate_ncbi[c(which(KP_abricate_ncbi$V15 == "CARBAPENEM" | KP_abricate_ncbi$V15 == "BETA-LACTAM" | KP_abricate_ncbi$V15 == "CEPHALOSPORIN")), ]
resistance_KP_ncbi <- KP_ncbi_cephBlaCarb[,c(6,15)] %>% distinct(V6,.keep_all=TRUE)
KP_ncbi_occ <- KP_ncbi_cephBlaCarb %>% count(V6)
KP_ncbi_occ <- merge(resistance_KP_ncbi, KP_ncbi_occ, by.x = "V6", by.y = "V6")
KP_ncbi_occ <- KP_ncbi_occ[order(KP_ncbi_occ$n,decreasing = TRUE), ]


# same for PA
PA_abricate_ncbi <- read.delim('~/data_project/assemblies/PA/abricate_ncbi.tab', header = FALSE, sep = "\t")
PA_abricate_ncbi$V1 <- gsub(".fna","", PA_abricate_ncbi$V1)
if (length(which(PA_abricate_ncbi$V11 < 95)) != 0) {
  PA_abricate_ncbi <- PA_abricate_ncbi[-c(which(PA_abricate_ncbi$V11 < 95)),]
}

PA_ncbi_cephBlaCarb <-PA_abricate_ncbi[c(which(PA_abricate_ncbi$V15 == "CARBAPENEM" | PA_abricate_ncbi$V15 == "BETA-LACTAM" | PA_abricate_ncbi$V15 == "CEPHALOSPORIN")), ]
resistance_PA_ncbi <- PA_ncbi_cephBlaCarb[,c(6,15)] %>% distinct(V6,.keep_all=TRUE)
PA_ncbi_occ <- PA_ncbi_cephBlaCarb %>% count(V6)
PA_ncbi_occ <- merge(resistance_PA_ncbi, PA_ncbi_occ, by.x = "V6", by.y = "V6")
PA_ncbi_occ <- PA_ncbi_occ[order(PA_ncbi_occ$n,decreasing = TRUE), ]


# same for AB
AB_abricate_ncbi <- read.delim('~/data_project/assemblies/AB/abricate_ncbi.tab', header = FALSE, sep = "\t")
AB_abricate_ncbi$V1 <- gsub(".fna","", AB_abricate_ncbi$V1)
if (length(which(AB_abricate_ncbi$V11 < 95)) != 0) {
  AB_abricate_ncbi <- AB_abricate_ncbi[-c(which(AB_abricate_ncbi$V11 < 95)),]
}

AB_ncbi_cephBlaCarb <-AB_abricate_ncbi[c(which(AB_abricate_ncbi$V15 == "CARBAPENEM" | AB_abricate_ncbi$V15 == "BETA-LACTAM" | AB_abricate_ncbi$V15 == "CEPHALOSPORIN")), ]
resistance_AB_ncbi <- AB_ncbi_cephBlaCarb[,c(6,15)] %>% distinct(V6,.keep_all=TRUE)
AB_ncbi_occ <- AB_ncbi_cephBlaCarb %>% count(V6)
AB_ncbi_occ <- merge(resistance_AB_ncbi, AB_ncbi_occ, by.x = "V6", by.y = "V6")
AB_ncbi_occ <- AB_ncbi_occ[order(AB_ncbi_occ$n,decreasing = TRUE), ]

# test for co-occurrence of genes
unique(AB_ncbi_cephBlaCarb$V1)
which(duplicated(AB_ncbi_cephBlaCarb$V1) == "FALSE")

a <- AB_ncbi_cephBlaCarb[,c(1,6)]
for(j in (unique(AB_ncbi_cephBlaCarb$V1))) {
  for(i in which(AB_ncbi_cephBlaCarb$V1 == j)) {
    a$j <- AB_ncbi_cephBlaCarb[i, 6]
  }
}  
data_frame_mod <- a[a$V1 %in% c("403606-21"),]


test_genes <- AB_ncbi_cephBlaCarb[, c(1,6)]
test_genes$V3 <- 1
test_heatmap <- ggplot(data = test_genes, mapping = aes(x = V6,
                                                        y = V1,
                                                        fill = V3)) +
  geom_tile(colour="white", size=1) +
  labs(x = "", y = "") +
  #scale_fill_manual(values = colors2) +
  theme(legend.position = 'none', axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90)) + 
  ggtitle("Abricate - NCBI - AB") + 
  coord_fixed(ratio = 0.25)
test_heatmap

#--------------------------------------------------------------------------------
EC_abricate_card <- read.delim('~/data_project/assemblies/EC/abricate_card.tab', header = FALSE, sep = "\t")
EC_abricate_card$V1 <- gsub(".fna","", EC_abricate_card$V1)
if (length(which(EC_abricate_card$V11 < 95)) != 0) {
  EC_abricate_card <- EC_abricate_card[-c(which(EC_abricate_card$V11 < 95)),]
}

patterns <- c("carbapenem", "cephalosporin")
card_cephBlaCarb <- filter(EC_abricate_card, grepl(paste(patterns, collapse = "|"), EC_abricate_card$V15))

genes_count_card <- card_cephBlaCarb %>% count(V6)
genes_count_card <- genes_count_card[order(genes_count_card$n,decreasing = TRUE), ]
