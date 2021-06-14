#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=VCall_Leaf
#SBATCH -p rausherlab 
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

/hpc/group/biodept/itl4/Programs/bcftools-1.10.2/bcftools mpileup \
        -d 1000 --threads 20 \
        -f /hpc/group/biodept/itl4/LacGenome/Ipomoea_Lep-Anchor_final.fasta \
        -b /hpc/group/biodept/itl4/RI_Austinii/bamlist_IpoLeaf_2Le.txt | \
/hpc/group/biodept/itl4/Programs/bcftools-1.10.2/bcftools call \
         --threads 20 -f gq -m -Oz \
         -o /work/itl4/Austinii/raw_calls_LeafIpo_2Le_2020.vcf.gz
