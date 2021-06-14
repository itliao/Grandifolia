#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=MergeLeaf
#SBATCH -p rausherlab,common,scavenger 
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

while read i

do
echo "Processing sample $i ..."

/dscrhome/itl4/Programs/samtools-1.9/samtools merge \
	/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/merged_$i.bam \
	/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/mergeready_L1_$i.bam\
	/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/mergeready_L2_$i.bam 
	
done < /dscrhome/itl4/RI_Austinii/IpoSp_merge.txt