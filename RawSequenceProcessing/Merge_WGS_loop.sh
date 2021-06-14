#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=mergeWGS
#SBATCH -p rausherlab,common,scavenger
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL


#updated 201127
while read i

do
echo "Processing sample $i ..."

/hpc/group/biodept/itl4/Programs/samtools-1.9/samtools merge \
	/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/All_merged/$i.bam \
	/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/mergeready_L1_$i.bam \
	/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/mergeready_L2_$i.bam \
	/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/mergeready_L3_$i.bam

done < /hpc/group/biodept/itl4/WGS/WGS_Ind_loop.txt

