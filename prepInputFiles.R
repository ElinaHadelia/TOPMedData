#preparing GWAS input file
library (dplyr)
library(stats)
library(tidyr)

args = commandArgs(TRUE)

GWAS = args[1]
GENE=args[2]
GWAS_STYDYTYPE=args[3]
GWAS_NUMBER=args[4]
OUTFILE_NAME = args[5]
DIRECTORY =args[6]

GWAS_summStats = read.delim (GWAS,  sep ="\t", header = F) %>% select(V1, V2, V10, V14, V15, V17)%>%
  rename(SNPID=V1, id=V2, mafVal=V10, logP=V14, CHR=V15, POS=V17) %>%
  separate(mafVal, c("delete","MAF"),sep="=") %>% select(-delete) %>%
  mutate(PVAL=10**(-logP))

GWASinfo=read.delim(GWAS_STYDYTYPE, sep="\t", header=F) %>% rename(id=V1, typeofStudy=V2) %>% 
  mutate(type=ifelse(typeofStudy=="CaseControl", "cc", "quant"))

GWASnumber = read.delim(GWAS_NUMBER, sep="\t", header=T) %>% select(id, sample_size,ncase) %>% rename(N=sample_size, Ncases=ncase)  

GWAS_combined= GWAS_summStats %>% left_join(GWASinfo,by="id") %>% left_join(GWASnumber,by="id") %>%select(PVAL, type, MAF, CHR, POS, N, Ncases, SNPID, id) %>%
  rename(ProbeID=id)

setwd(DIRECTORY)
write.table (GWAS_combined, OUTFILE_NAME, quote = F,sep = "\t", row.names = F)

cat ("The end\n")






#preparing my input file (aka eql)
library (dplyr)
library(stats)
library(tidyr)

args = commandArgs(TRUE)

INPUT_MINE = args[1]
GENE=args[2]
DIRECTORY=args[3]
OUTFILE_NAME = args[4]

my_results = read.delim (INPUT_MINE,  sep ="\t", header = T) %>% select(SNP, Weight, P.value, allele1FreqInALL)%>%
    rename(SNPID=SNP, N=Weight, PVAL=P.value, MAF=allele1FreqInALL) %>% mutate(ProbeID=GENE) %>%
    separate(SNPID, c("chrom","snpName"),sep=":", remove=FALSE) %>% separate(snpName, c("POS","Ref","Alt"),sep="_")
  
my_results$chrom <- sub("^","chr",my_results$chrom)
my_results<- my_results %>%  rename(CHR= chrom) %>% select(-Ref, -Alt) %>% select(PVAL, ProbeID, CHR, POS, MAF, N, SNPID)

setwd(DIRECTORY)
write.table (my_results, OUTFILE_NAME, quote = F,sep = "\t", row.names = F)

cat ("The end\n")
