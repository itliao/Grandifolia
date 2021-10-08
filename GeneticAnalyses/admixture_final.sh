#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=5
#SBATCH --job-name=admixture_test
#SBATCH -p rausherlab,common,scavenger
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

#test admixture for K =2-10 
for K in {2..10}; do
/hpc/group/biodept/itl4/Programs/dist/admixture_linux-1.3.0/admixture \
	--cv -s 100 /work/itl4/Austinii/admixture/plink6_210110.bed $K -j10 
done