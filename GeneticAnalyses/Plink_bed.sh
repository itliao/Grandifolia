#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Plink_test
#SBATCH -p rausherlab,common,scavenger
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

#bed file for both admixture and faststructure
/hpc/group/biodept/itl4/Programs/plink1.9/plink \
	--vcf /work/itl4/Austinii/filtered/filter6_201218.vcf.gz \
	--double-id \
	--maf 0.05 \
	--make-bed \
	--out /work/itl4/Austinii/admixture/plink6_210110
