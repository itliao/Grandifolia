#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=ValidateBAM
#SBATCH -p rausherlab,common,scavenger 
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL


while read i
do

echo "Processing sample $i ..." 

java -jar /dscrhome/itl4/Programs/picard/build/libs/picard.jar ValidateSamFile \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/merged_$i.bam \
	O=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/out_$i.out \
	MODE=SUMMARY 

done < /dscrhome/itl4/RI_Austinii/IpoSp_merge.txt


while read i
do

echo "Processing sample $i ..." 

java -jar /dscrhome/itl4/Programs/picard/build/libs/picard.jar ValidateSamFile \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/mergeready_$i.bam \
	O=/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/out_$i.out \
	MODE=SUMMARY 

done < /dscrhome/itl4/RI_Austinii/IpoSp_AC_clean.txt

