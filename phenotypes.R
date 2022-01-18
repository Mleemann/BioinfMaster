#install.packages("tidyverse")
library(tidyr)
library(dplyr)

# E. faecium

# load file
metadata <- 'KEIMIDENTIFIKATION_ANTIBIOGRAM_RESULTATE_EnterococcusFaecium.csv'
all_metadata <- read.csv(metadata, sep=';')

# remove first part of LSM_BK
all_metadata$LSM_BK <- gsub('(\\d*\\|)', '', all_metadata$LSM_BK)

# get year for the sample_id
all_metadata$LSM_SAMPLE_ENTRY_DATE_TS <- gsub('(^\\d{4})(.*)','\\1', all_metadata$LSM_SAMPLE_ENTRY_DATE_TS)
all_metadata$LSM_SAMPLE_ENTRY_DATE_TS <- gsub('(^\\d{2})','', all_metadata$LSM_SAMPLE_ENTRY_DATE_TS)

# combine sample-nr. with the year to sample_id
all_metadata$LSM_BK <- paste0(all_metadata$LSM_BK, '-',all_metadata$LSM_SAMPLE_ENTRY_DATE_TS)
names(all_metadata)[2] <- c('Sample_id')

# get LSM_BK, Antibiotic, Result
ab_results <- all_metadata[,c(2,29,30)]
mic <- all_metadata[,c(2,29,35,36)]

# combine prefix and value of mic
mic$LARM_MIC <- paste0(mic$LARM_MIC_PREFIX, ' ', mic$LARM_MIC)
mic <- mic[ ,c(1,2,4)]

# remove rows with empty fields
ab_results <- ab_results[!(ab_results$LARM_ANTIBIOTIC == ''), ]
mic <- mic[!(mic$LARM_ANTIBIOTIC == ''), ]

# transform antibiotics from rows to columns
all_phenotypes <- pivot_wider(ab_results, names_from = LARM_ANTIBIOTIC, values_from = LARM_INTERNAL_RESULT_ID)
all_mic <- pivot_wider(mic, names_from = LARM_ANTIBIOTIC, values_from = LARM_MIC)

# import sample_ids to be analyzed
samples_to_analyze <- read.table("samples_EF.txt")

# subset phenotype dataframe with the samples to be analyzed
phenotypes <- all_phenotypes[all_phenotypes$Sample_id %in% samples_to_analyze$V1, ]
mic_data <- all_mic[all_mic$Sample_id %in% samples_to_analyze$V1, ]

# create phenotype and mic table for only Vancomycin
vanco_pheno <- phenotypes[,c(1,31)]
vanco_mic <- mic_data[,c(1,31)]
names(vanco_mic)[2] <- c('MIC-Vancomycin')

pheno_mic <- merge(vanco_pheno, vanco_mic, by.x = 'Sample_id')

# save phenotype and mic table
pheno_mic <- apply(pheno_mic,2,as.character)
write.csv2(pheno_mic, file = 'phenotypes_vanco_EF.csv', row.names = F, quote = F)


# save phenotype and mic table
#phenotypes <- apply(phenotypes,2,as.character)
#write.csv2(phenotypes, file = 'phenotypes_EF.csv', row.names = F, quote = F)

#mic_data <- apply(mic_data,2,as.character)
#write.csv2(mic_data, file = 'mic_EF.csv', row.names = F, quote = F)


#-------------------------------------------------------------------------------
# S. aureus

# load file
metadata <- 'KEIMIDENTIFIKATION_ANTIBIOGRAM_RESULTATE_StaphAureus.csv'
all_metadata <- read.csv(metadata, sep=';')

# remove first part of LSM_BK
all_metadata$LSM_BK <- gsub('(\\d*\\|)', '', all_metadata$LSM_BK)

# get year for the sample_id
all_metadata$LSM_SAMPLE_ENTRY_DATE_TS <- gsub('(^\\d{4})(.*)','\\1', all_metadata$LSM_SAMPLE_ENTRY_DATE_TS)
all_metadata$LSM_SAMPLE_ENTRY_DATE_TS <- gsub('(^\\d{2})','', all_metadata$LSM_SAMPLE_ENTRY_DATE_TS)

# combine sample-nr. with the year to sample_id
all_metadata$LSM_BK <- paste0(all_metadata$LSM_BK, '-',all_metadata$LSM_SAMPLE_ENTRY_DATE_TS)
names(all_metadata)[2] <- c('Sample_id')

# get LSM_BK, Antibiotic, Result
ab_results <- all_metadata[,c(2,19,29,30)]
mic <- all_metadata[,c(2,19,29,35,36)]

# combine prefix and value of mic
mic$LARM_MIC <- paste0(mic$LARM_MIC_PREFIX, ' ', mic$LARM_MIC)
mic <- mic[ ,c(1,2,3,5)]

# remove rows with empty fields
ab_results <- ab_results[!(ab_results$LARM_ANTIBIOTIC == ''), ]
mic <- mic[!(mic$LARM_ANTIBIOTIC == ''), ]

# transform antibiotics from rows to columns
all_phenotypes <- pivot_wider(ab_results, names_from = LARM_ANTIBIOTIC, values_from = LARM_INTERNAL_RESULT_ID)
all_mic <- pivot_wider(mic, names_from = LARM_ANTIBIOTIC, values_from = LARM_MIC)

# import sample_ids to be analyzed
samples_to_analyze <- read.table("samples_SA.txt")

# subset phenotype dataframe with the samples to be analyzed
phenotypes <- all_phenotypes[all_phenotypes$Sample_id %in% samples_to_analyze$V1, ]
phenotypes <- phenotypes[-c(182),]
mic_data <- all_mic[all_mic$Sample_id %in% samples_to_analyze$V1, ]
colnames(mic_data) <- paste("MIC", colnames(mic_data), sep = "-")
names(mic_data)[1] <- c('Sample_id')
mic_data <- mic_data[-c(182),]

pheno_mic <- merge(phenotypes, mic_data, by.x = 'Sample_id', by.y = 'Sample_id')

# save phenotype and mic table
pheno_mic <- apply(pheno_mic,2,as.character)
write.csv2(pheno_mic, file = 'pheno_mic_SA.csv', row.names = F, quote = F)
phenotypes <- apply(phenotypes,2,as.character)
write.csv2(phenotypes, file = 'phenotypes_SA.csv', row.names = F, quote = F)
mic_data <- apply(mic_data,2,as.character)
write.csv2(mic_data, file = 'mic_SA.csv', row.names = F, quote = F)


