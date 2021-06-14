#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=Trim
#SBATCH -p common,scavenger,rausherlab
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

while read i
do

echo "Processing sample $i ..." 

java -jar /dscrhome/itl4/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	/work/itl4/Austinii/Raw_sequence/$i\_R1_001.fastq.gz \
	/work/itl4/Austinii/Raw_sequence/$i\_R2_001.fastq.gz \
	/work/itl4/Austinii/Trim/$i\_R1_paired.fastq.gz \
	/work/itl4/Austinii/Trim/$i\_R1_unpaired.fastq.gz \
	/work/itl4/Austinii/Trim/$i\_R2_paired.fastq.gz \
	/work/itl4/Austinii/Trim/$i\_R2_unpaired.fastq.gz \
	ILLUMINACLIP:/dscrhome/itl4/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
	LEADING:3 TRAILING:3 MINLEN:36

done < /dscrhome/itl4/RI_Austinii/IpoInd.txt

# for Ipomoea lacunosa and Ipomoea cordatotriloba, use IpoSpAC.txt