#obtaining from logistic file what is needed for METAL

cd /YOURDIRECTORY
Exp=afterLogForMetal
ls -d */ | sed 's/.$//'| grep -v phenotypes| while read INPUT; do
	cd /YOURDIRECTORY/$INPUT
	ls -d */ |  sed 's/.$//'| while read ANCESTRY; do
		cd /YOURDIRECTORY/$INPUT/$ANCESTRY
		rm $ANCESTRY.logistic.txt
		ls *.assoc.logistic|while read IN; do	 			
    		cat $IN|awk '$5=="ADD" && $7!="NA"'|awk -F'_' '{print $0,$3}' OFS="\t"|awk 'BEGIN {print"SNP\tSNP1\tSNP2\tNMISS\tOR\tP"}{print $2,$4,$10,$6,$7,$9}' OFS="\t" > $ANCESTRY.logistic.txt
		done
	done	
done
