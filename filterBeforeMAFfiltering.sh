# first, create the list of files to be merged
cd /YOURDIRECTORY
ls -d */| sed 's/.$//'| while read INPUT; do
	cd /YOURDIRECTORY/$INPUT/
	ls -d */ | sed 's/.$//'| while read ANCESTRY; do
	cd /YOURDIRECTORY/$INPUT/$ANCESTRY/
	ls *.bed| grep -v log|sed 's/.bed//g'| grep -v MESA|sort -V |grep -v $ANCESTRYlist.list.txt > $INPUT.$ANCESTRY.list.txt
	done
done



###Filtering -  After merging, create a txt file of SNPs to keep(only keep SNPs with 1 nucleotide in reference and reported)
#cat combined.bim | awk 'length($5) ==1 && length($6) ==1'|cut -f2 > SNPsToKeep.txt
cd /YOURDIRECTORY
ls -d */| sed 's/.$//'|while read INPUT; do
	cd /YOURDIRECTORY/$INPUT
	ls -d */ |sed 's/.$//'| while read ANCESTRY; do
		cd /YOURDIRECTORY/$INPUT/$ANCESTRY
		cat $INPUT.$ANCESTRY.combined.bim| awk 'length($5) !=1 || length($6) !=1'|cut -f2 > $INPUT.SNPsToRemove.txt
	done
done

 

# To calculate the fraction of missingness per SNV, generate .lmiss[stats per SNP](also generates .imiss file[stats per sample]) file using:
#plink --bfile newCombined --missing --out newCombinedTable
#cat newCombinedTable.lmiss| cut -f5|sort -V|head

cd /YOURDIRECTORY
Exp=keepingOnlySingleNucleotideVariants
ls -d */| sed 's/.$//'|while read INPUT; do
	cd /YOURDIRECTORY/$INPUT
	ls -d */ |sed 's/.$//'|while read ANCESTRY; do
		cd /YOURDIRECTORY/$INPUT/$ANCESTRY
		ls *combined.SNVs.bed|sed 's/.bed//g' |while read IN; do 
			plink --bfile $IN --missing --out $IN; rm -f $IN.nosex
		done
	done
done


		
# to check the F_MISS (Proportion of sample missing for this SNP) . To see if it's other than 0
cd /YOURDIRECTORY
Exp=keepingOnlySingleNucleotideVariants
ls -d */| sed 's/.$//'|while read INPUT; do
	cd /YOURDIRECTORY/$INPUT
	ls -d */ |sed 's/.$//'|while read ANCESTRY; do
		cd /YOURDIRECTORY/$INPUT/$ANCESTRY
		ls *.lmiss|while read IN; do 
		cat $IN| cut -f5|sort -V| awk '{if($5!=0 && $5!="F_MISS") print $1,$2,$3,$4,$5}'
		done
	done
done



