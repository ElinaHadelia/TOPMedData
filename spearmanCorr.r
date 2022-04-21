library (ggplot2)
library (dplyr)
library(tidyverse)
library(stats)
library(scales)
library(gridExtra)
library(data.table)

args = commandArgs(TRUE)

META_ANALYSIS = args[1]
GWAS = args[2]
OUTFILE_NAME = args[3]
OUT_DIR = args[4]
GWAS_INFO = args[5]
GENE=args[6]

#metalResults results
metalResults_TOPMed_final = read.delim (META_ANALYSIS,  sep ="\t", header = T) %>% select( SNP, P.value, allele1FreqInALL)%>%
  separate(SNP, c("SNP", "REF","ALT"), sep = "_")

#filtered TOPmed results
metalResults_TOPMed_filtered<-filter(metalResults_TOPMed_final,allele1FreqInALL > 0.001, P.value<=1e-05) %>% mutate(rankTOPMed = dense_rank(P.value))
numbOfSnpsTOPMed001 = count (metalResults_TOPMed_filtered) %>% as.numeric

#filtering TOPmed results to know the number of SNPs when using MAF005
forSNPnumber<-filter(metalResults_TOPMed_final,allele1FreqInALL > 0.005) 
numbOfSnpsTOPMed005 = count (forSNPnumber) %>% as.numeric

#for GWAS studies
GWASstudy = read.delim (GWAS,  sep ="\t", header =F) %>% select(V1,V2,V5,V6,V10,V13,V15) %>%
  rename(gwasStudy=V1, chr.x=V2, REF=V5, ALT=V6, gene =V10, posBegin_hg38=V15, LP_GWAS=V13) %>%
  mutate(posBegin_hg38_1base= posBegin_hg38+1) %>%
  mutate(SNP=paste(chr.x,":",posBegin_hg38_1base, sep="")) %>%
  mutate(SNP_full=paste(chr.x,":",posBegin_hg38_1base,"_", REF,"_",ALT,sep=""))  %>% 
  mutate(GWAS.P=10**(-LP_GWAS))
#need to add these 2 lines of code for CSTB, NOTCH2, and SETD6 instead of the last mutate
#GWASstudy$new_LP_GWAS<- as.numeric(as.character(GWASstudy$LP_GWAS)) 
#GWASstudy<- GWASstudy%>% mutate(GWAS.P=10**(-new_LP_GWAS))

setwd(OUT_DIR)

final <-GWASstudy %>% filter(GWAS.P<=1e-03) %>%group_by(gwasStudy) %>%  mutate(rankGWAS= dense_rank(GWAS.P)) %>% data.frame
final_forSpearman = final %>% left_join(metalResults_TOPMed_filtered, by="SNP") %>% distinct %>% select(gwasStudy,gene,SNP,GWAS.P,P.value,rankTOPMed,rankGWAS)

#to check number of snps overlapping aka count number of snsps per study in TOPmed when using MAF 001
toCheck001 = GWASstudy %>% inner_join(metalResults_TOPMed_filtered, by="SNP") %>% select(gwasStudy,SNP_full) %>% group_by(gwasStudy) %>%
  summarize (nSNPsinTOPmedPerStudyMAF001 =length(unique(SNP_full)))

#to check number of snps overlapping aka count number of snsps per study in TOPmed when using MAF 005
toCheck005 = GWASstudy %>% inner_join(forSNPnumber, by="SNP") %>% select(gwasStudy,SNP_full) %>% group_by(gwasStudy) %>%
  summarize (nSNPsinTOPmedPerStudyMAF005 =length(unique(SNP_full)))

#to make stats to add later
stats <- GWASstudy %>% group_by(gwasStudy) %>% summarize (nsnpsGWASForThisGene_1MbTotal = length(unique(SNP_full)))
stats_p03 <- final %>% group_by(gwasStudy) %>% summarize (nsnpsGWASafterFilteringP03 = n(), minPvalGWASafterFilteringP03= min(GWAS.P, na.rm=T))

numberOfRows = count (final_forSpearman) %>% as.numeric

if (file.exists(OUTFILE_NAME)) { unlink (OUTFILE_NAME) }
for (i in 1:numberOfRows) {
  # GWASstudyName=final_forSpearman[i] %>% distinct
  GWASstudyName<- unique(final_forSpearman$gwasStudy) %>% as.character 
  
  y <- final_forSpearman %>% filter (gwasStudy ==GWASstudyName[i]) 
  leftover<- y  %>% filter(P.value!= "NA" | GWAS.P== "NA")
  min = y %>% summarize (minPGWASP= min(GWAS.P, na.rm=T)) %>% as.numeric 
  
  if(nrow(leftover) > 2)
  {
    tobj = try (resultSpearman <- cor.test(leftover$P.value, leftover$GWAS.P, method="spearman", exact = FALSE),silent=TRUE)
    if(is(tobj,"try-error"))
    {
      SpearmanR <- NA
      nObsForSpear <- NA
      SpearmanP <- NA
    } else {
      SpearmanR <- resultSpearman$estimate 
      nObsForSpear <- nrow (leftover)
      SpearmanP <- resultSpearman$p.value
    }
  } else {
    SpearmanR <- NA
    nObsForSpear <- NA
    SpearmanP <- NA
  }
  
  results <- c(GWASstudyName[i],SpearmanR, nObsForSpear,SpearmanP,min,GENE); rm ( SpearmanR, nObsForSpear, SpearmanP, gene)
  names(results) =c("gwasStudy","SpearmanR","nObsForSpear", "SpearmanP","minPGWASForThisGene","gene")
  results <- as.data.frame(t(results))
  cat("... Saving data ....\n")
  if (i==1) {
    write.table(results,OUTFILE_NAME ,append = T,quote = F,sep = "\t", row.names = F,col.names = T)
    rm (results)
  } else {
    write.table(results,OUTFILE_NAME,append = T,quote = F,sep = "\t", row.names = F,col.names = F)
    rm (results)
  }
}

cat ("Annotate final file\n")
spearman_results <- fread (OUTFILE_NAME, sep ="\t", check.names = F, header = T) %>% filter(SpearmanR!= "NA" & SpearmanP!=0)
gwasStudyInfo <- read.delim (GWAS_INFO, sep ="\t", check.names = F, header = T) %>% select (-batch) %>%
  rename("gwasStudy"="id", "GWAStrait"= "trait", "populationGWAS" = "population","sampleSize_totalGWAS"="nsnp","TotalVariants_totalGWAS"="nsnp", "ncase_totalGWAS"="ncase","ncontrol_totalGWAS"="ncontrol")

finalToWrite <- spearman_results %>% inner_join(stats, by = "gwasStudy") %>% inner_join(stats_p03, by = "gwasStudy")  %>% 
  inner_join (gwasStudyInfo, by = "gwasStudy")%>% 
  inner_join(toCheck001, by = "gwasStudy")  %>% mutate(numberOfSnpsInTOPMedFilteredMAF001=numbOfSnpsTOPMed001) %>%
  inner_join(toCheck005, by = "gwasStudy")  %>%  mutate(numberOfSnpsInTOPMedFilteredMAF005=numbOfSnpsTOPMed005)
setwd(OUT_DIR)
write.table (finalToWrite, OUTFILE_NAME, quote = F,sep = "\t", row.names = F)

cat ("The end\n")

