library(tidyverse)
library(data.table)
library (dplyr)

args = commandArgs(TRUE)

METAANALYSIS = args[1]
LOCI = args[2]
ALLELEFREQ_EUR = args[3]
ALLELEFREQ_AMR = args[4]
ALLELEFREQ_AFR = args[5]
LOGISTIC_EUR = args[6]
LOGISTIC_AMR = args[7]
LOGISTIC_AFR = args[8]
FREQUENCY_OFALL = args[9]
OUT_DIR = args[10]


metal_results <- read.delim (METAANALYSIS, sep ="\t", header = T) %>% rename(SNP=MarkerName)
alleleEUR <- read.delim (ALLELEFREQ_EUR, sep ="", header = T) %>% select(SNP, MAF_A, MAF_U, NCHROBS_A, NCHROBS_U ) %>% mutate(NCase_EUR=NCHROBS_A/2,NCntrl_U_EUR=NCHROBS_U/2 ) %>% select(-NCHROBS_A,-NCHROBS_U)
alleleAFR <- read.delim (ALLELEFREQ_AFR, sep ="", header = T) %>% select(SNP, MAF_A, MAF_U, NCHROBS_A, NCHROBS_U ) %>% mutate(NCase_AFR=NCHROBS_A/2,NCntrl_U_AFR=NCHROBS_U/2 ) %>% select(-NCHROBS_A,-NCHROBS_U)
alleleAMR <- read.delim (ALLELEFREQ_AMR, sep ="", header = T) %>% select(SNP, MAF_A, MAF_U, NCHROBS_A, NCHROBS_U ) %>%mutate(NCase_AMR=NCHROBS_A/2,NCntrl_U_AMR=NCHROBS_U/2 ) %>% select(-NCHROBS_A,-NCHROBS_U)


if (file.exists(LOGISTIC_EUR)){
logisticEUR <- read.delim (LOGISTIC_EUR, sep ="\t", header = T) %>% select(SNP, OR, P)}

if (file.exists(LOGISTIC_AFR)){
logisticAFR <- read.delim (LOGISTIC_AFR, sep ="\t", header = T) %>% select(SNP, OR, P)}

if (file.exists(LOGISTIC_AMR)){
logisticAMR <- read.delim (LOGISTIC_AMR, sep ="\t", header = T) %>% select(SNP, OR, P)}



combined1<- left_join(metal_results,alleleEUR, by="SNP") %>% rename(MAFinCase_EUR=MAF_A, MAFinCntrl_EUR=MAF_U)

if (file.exists(LOGISTIC_EUR)){
  combined2<-left_join(combined1,logisticEUR,by="SNP") %>% rename(OR_EUR=OR, P_EUR=P)}

if (file.exists(LOGISTIC_EUR)){
  combined3<- left_join(combined2,alleleAFR, by="SNP") %>% rename(MAFinCase_AFR=MAF_A, MAFinCntrl_AFR=MAF_U)} else{
  combined3<- left_join(combined1,alleleAFR, by="SNP") %>% rename(MAFinCase_AFR=MAF_A, MAFinCntrl_AFR=MAF_U)}

if (file.exists(LOGISTIC_AFR)){
  combined4<-left_join(combined3,logisticAFR, by="SNP") %>% rename(OR_AFR=OR, P_AFR=P)}

if (file.exists(LOGISTIC_AFR)){
  combined5<- left_join(combined4,alleleAMR, by="SNP") %>% rename(MAFinCase_AMR=MAF_A, MAFinCntrl_AMR=MAF_U)} else{
  combined5<- left_join(combined3,alleleAMR, by="SNP") %>% rename(MAFinCase_AMR=MAF_A, MAFinCntrl_AMR=MAF_U)}

if (file.exists(LOGISTIC_AMR)){
  combinedFinal<-left_join(combined5,logisticAMR,by="SNP") %>% rename(OR_AMR=OR, P_AMR=P)} else{
  combinedFinal<-combined5}

frequency_of_all <- read.delim (FREQUENCY_OFALL, sep ="", header = T) %>% select(SNP, MAF) %>% rename(allele1FreqInALL=MAF)

final<- left_join(combinedFinal,frequency_of_all,by="SNP")

setwd(OUT_DIR)
write.table (final, paste0("METAANALYSIS1.", LOCI,".annotated.txt"), sep ="\t",row.names =F , quote = F)

cat ("The end\n")

