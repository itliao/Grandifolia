#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=PicardAus_3
#SBATCH -p rausherlab,common,scavenger 
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

##update 201008

while read i
do

echo "Processing sample $i ..." 

java -jar /dscrhome/itl4/Programs/picard/build/libs/picard.jar FixMateInformation \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/$i\_L003Aligned.sortedByCoord.out.bam

java -jar /dscrhome/itl4/Programs/picard/build/libs/picard.jar CleanSam \
	VALIDATION_STRINGENCY=LENIENT \
	INPUT=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/$i\_L003Aligned.sortedByCoord.out.bam \
	OUTPUT=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/tmp_clean_L3_$i.bam

java -jar /dscrhome/itl4/Programs/picard/build/libs/picard.jar SortSam \
	INPUT=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/tmp_clean_L3_$i.bam \
	OUTPUT=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/tmp_clean_sorted_L3_$i.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/tmp
	
java -jar /dscrhome/itl4/Programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/tmp_clean_sorted_L3_$i.bam \
	O=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/mergeready_L3_$i.bam  \
	RGID=3$i \
	RGLB=lane3 \
	RGPL=illumina \
	RGPU=unit3$i \
	RGSM=sample$i 
	  
	    
done < /dscrhome/itl4/RI_Austinii/IpoSp_3.txt
