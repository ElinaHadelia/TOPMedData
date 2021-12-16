library (ggplot2)
library (dplyr)
library(tidyverse)
library(stats)
library(scales)

args = commandArgs(trailingOnly=TRUE)

META_ANALYSIS = args[1]
OUTFILE_NAME = args[2]
OUT_DIR = args[3]
GENE = args[4]
LOCIBP= args[5]

#the file to plot needs to have base-pair location and p-value of each SNP
metalResults = read.delim (META_ANALYSIS,  sep ="\t", header = T) %>% select(SNP, P.value)

final<- metalResults %>% separate(col=SNP, into=c("a", "b", "c"), sep="_") %>%
  separate(col=a, into=c("a", "SNP_bp"), sep=":") %>% select(-a, -b) %>% data.frame %>% mutate(SNP_kb= as.numeric(SNP_bp)/1000) %>%
  mutate(FDR=p.adjust(P.value, method="fdr"), bonf=p.adjust(P.value, method="bonferroni"))

nSNPs = final %>% select(SNP_bp)  %>% distinct %>% nrow                    

loci_start = strsplit(LOCIBP, '-') [[1]] [1] %>% as.numeric 
loci_start<- loci_start/1000
loci_end = strsplit(LOCIBP, '-') [[1]] [2] %>%as.numeric
loci_end<- loci_end/1000
 
g <- ggplot(final, aes(x=SNP_kb, y=-log10(FDR))) +
  geom_point(alpha=0.8, size=1.3, color =  "grey") +    # Show all points
  geom_point(data = subset ( final, -log10(FDR) >5, alpha=0.8, size=1.3, color =  "black")) +  
  scale_y_continuous(expand = c(0, 0.3)) +     # remove space between plot area and x axis
  scale_x_continuous(breaks = pretty_breaks(),labels = comma) +
  theme_bw() +
  geom_vline(xintercept = c(loci_start, loci_end), linetype='dotted', col = 'red') +
  labs (x = "Genomic position(kb)",
        y ="FDR(-log10)",
        title =paste0(GENE, ":", LOCIBP, "± 1 Mb"),
        subtitle = paste0("Number of SNPs: ", format(nSNPs, big.mark   = ","))) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(face="bold",size = 16,colour = "black",hjust = 0.5),
    plot.subtitle = element_text(size = 12,colour= "black",hjust = 0.5),
    axis.title = element_text (size = 12,face= "bold"),
    axis.text.x= element_text(vjust = 0.5,hjust =0.5,color = "black", size = 12,angle= 0), #,face="bold" vjust centered to tick
    axis.text.y= element_text(colour = "black", size = 12,vjust = 0.5,hjust = 1)) 

png(OUTFILE_NAME, width = 3000, height = 2400, res = 300, pointsize = 12)
setwd(OUT_DIR)
print(g)
dev.off()



