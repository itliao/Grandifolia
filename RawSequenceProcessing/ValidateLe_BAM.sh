#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=ValidateBAM
#SBATCH -p rausherlab,common,scavenger 
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

#updated 201127
while read i
do

echo "Processing sample $i ..." 

/hpc/group/biodept/itl4/Programs/jdk1.8.0_211/bin/java -jar \
	/hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar ValidateSamFile \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/All_merged/$i.bam \
	O=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/All_merged/out_$i.out \
	MODE=SUMMARY
	
done < /hpc/group/biodept/itl4/WGS/WGS_Ind_loop.txt 