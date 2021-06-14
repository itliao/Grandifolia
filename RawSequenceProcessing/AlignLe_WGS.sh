#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=AlignWGS_A
#SBATCH -p rausherlab,common,scavenger
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

##############
# for lane 1 #
##############
while read i

do
echo "Processing sample $i ..." 

/hpc/group/biodept/itl4/Programs/NextGenMap-0.5.5/bin/ngm-0.5.5/ngm \
        -r /hpc/group/biodept/itl4/LacGenome/Ipomoea_Lep-Anchor_final.fasta \
        -b -1 /work/itl4/WGS/RAUSHER_5099_180928A1/$i\_L003_R1_001.fastq.gz -2 /work/itl4/WGS/RAUSHER_5099_180928A1/$i\_L003_R2_001.fastq.gz \
        -o /work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL1_$i.bam \

done < /hpc/group/biodept/itl4/WGS/WGS_Ind_A.txt 

##############
# for lane 2 #
##############
while read i
do

echo "Processing sample $i ..." 

/hpc/group/biodept/itl4/Programs/NextGenMap-0.5.5/bin/ngm-0.5.5/ngm \
	-r /hpc/group/biodept/itl4/LacGenome/Ipomoea_Lep-Anchor_final.fasta \
	-b -1 /work/itl4/WGS/RAUSHER_5099_181005A1/$i\_L007_R1_001.fastq.gz -2 /work/itl4/WGS/RAUSHER_5099_181005A1/$i\_L007_R2_001.fastq.gz \
	-o /work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL2_$i.bam \

done < /hpc/group/biodept/itl4/WGS/WGS_Ind_B.txt 

##############
# for lane 3 #
##############

while read i
do

echo "Processing sample $i ..." 

/hpc/group/biodept/itl4/Programs/NextGenMap-0.5.5/bin/ngm-0.5.5/ngm \
	-r /hpc/group/biodept/itl4/LacGenome/Ipomoea_Lep-Anchor_final.fasta \
	-b -1 /work/itl4/WGS/RAUSHER_5099_181005A1/$i\_L008_R1_001.fastq.gz -2 /work/itl4/WGS/RAUSHER_5099_181005A1/$i\_L008_R2_001.fastq.gz \
	-o /work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL3_$i.bam \

done < /hpc/group/biodept/itl4/WGS/WGS_Ind_B.txt 