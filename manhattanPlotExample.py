library (ggplot2)
library (dplyr)
library(tidyverse)
library(stats)
library(scales)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

META_ANALYSIS_WITH_GWAS = args[1]
OUTFILE_NAME = args[2]
OUT_DIR = args[3]
GENE = args[4]
LOCIBP= args[5]

#the file to plot needs to have base-pair location and p-value of each SNP
metalResultsWithGWAS = read.delim (META_ANALYSIS_WITH_GWAS,  sep ="\t", header = F) %>% select(V2, V4, V9, V27, V32) %>% 
  rename(SNP_bp=V2, SNP=V4, metal.P=V9, allele1FreqInALL=V27, GWAS.P=V32) %>%
  data.frame %>%  mutate(SNP_kb= as.numeric(SNP_bp)/1000) 

nSNPs = metalResultsWithGWAS %>% select(SNP_bp)  %>% distinct %>% nrow                    

loci_start = strsplit(LOCIBP, '-') [[1]] [1] %>% as.numeric 
loci_start<- loci_start/1000
loci_end = strsplit(LOCIBP, '-') [[1]] [2] %>%as.numeric
loci_end<- loci_end/1000

BREAKS <- c(0, 0.001, 0.01, 0.1, 1)
toplot <- metalResultsWithGWAS%>% 
  mutate (Bins_MAF= cut(allele1FreqInALL,breaks= BREAKS,right = FALSE,include.lowest = TRUE),
          highlight=ifelse(-log10(metal.P)>= 5,"yes","no")) %>% data.frame

THEME= theme( 
  legend.position="bottom",
  legend.direction="horizontal",
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  plot.title = element_text(face="bold",size = 16,colour = "black",hjust = 0.5),
  plot.subtitle = element_text(size = 12,colour= "black",hjust = 0.5),
  axis.title = element_text (size = 12,face= "bold"),
  axis.text.x= element_text(vjust = 0.5,hjust =0.5,color = "black", size = 12,angle= 0), #,face="bold" vjust centered to tick
  axis.text.y= element_text(colour = "black", size = 12,vjust = 0.5,hjust = 1))
 
g <- ggplot(toplot, aes(x=SNP_kb, y=-log10(metal.P))) +
  geom_point(alpha=0.5, size=1.3, color =  "grey") +
  geom_point(data = subset (toplot, highlight=="yes"), size=1.5, aes(color =  Bins_MAF)) +
  scale_color_manual(name ="MAF", values = c("[0,0.001)"="#80e27e", "[0.001,0.01)"="#00600f", "[0.01,0.1)"="#FF7566",  "[0.1,1]"="#EFA6FF")) +
  scale_y_continuous(expand = c(0, 0.6)) +     # remove space between plot area and x axis
  scale_x_continuous(breaks = pretty_breaks(),labels = comma) +
  theme_bw() +
  geom_vline(xintercept = c(loci_start, loci_end), linetype='dotted', col = 'red') +
  labs (x = "Genomic position(kb)",
        y ="P(-log10) from meta-analysis",
        title =paste0(GENE, ":", LOCIBP, " ± 1 Mb"),
        subtitle = paste0("Number of SNPs: ", format(nSNPs, big.mark   = ","))) +
  THEME

toPlotGWAS<- toplot
#use the following instead for studies with missing GWAS.P values
toPlotGWAS<- toplot%>% filter(GWAS.P!=-1.000e+00)
#use the following instead for heelBone study
toPlotGWAS<- toplot%>% filter(GWAS.P!=".")
toPlotGWAS <- transform(toPlotGWAS, GWAS.P = as.numeric(GWAS.P))

nSNPs_GWAS = toPlotGWAS %>% select(SNP_bp)  %>% distinct %>% nrow        
h <- ggplot(toPlotGWAS, aes(x=SNP_kb, y=-log10(GWAS.P))) +
  geom_point(alpha=0.5, size=1.3, color =  "grey") +
  geom_point(data = subset (toplot, highlight=="yes"), size=1.5, aes(color =  Bins_MAF)) + #change toplot to toPlotGWAS for heelBone study
  scale_color_manual(name ="MAF", values = c("[0,0.001)"="#80e27e", "[0.001,0.01)"="#00600f", "[0.01,0.1)"="#FF7566",  "[0.1,1]"="#EFA6FF")) +
  scale_y_continuous(expand = c(0, 0.6)) +     # remove space between plot area and x axis
  scale_x_continuous(breaks = pretty_breaks(),labels = comma) +
  theme_bw() +
  geom_vline(xintercept = c(loci_start, loci_end), linetype='dotted', col = 'red') +
  labs (x = "Genomic position(kb)",
        y ="P(-log10) from GWAS",
        title =paste0(GENE, ":", LOCIBP, " ± 1 Mb. Results for heelBone from GWAS"),
        #insert for which trait it is
        subtitle = paste0("Number of SNPs: ", format(nSNPs_GWAS, big.mark   = ","))) +
  THEME 

png(OUTFILE_NAME, width = 3000, height = 2600, res = 300, pointsize = 12); grid.arrange (g,h, nrow = 2)
setwd(OUT_DIR)
dev.off()



