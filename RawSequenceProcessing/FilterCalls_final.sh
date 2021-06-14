#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=FilterRI_snp
#SBATCH -p rausherlab,common,scavenger
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

#update 201218
/hpc/group/biodept/itl4/Programs/bcftools-1.10.2/bcftools view -m2 -M2 -v snps -S \
	/hpc/group/biodept/itl4/RI_Austinii/filtersamplesLe.txt \
	/work/itl4/Austinii/raw_calls_LeafIpo_2Le_2020.vcf.gz | \
/hpc/group/biodept/itl4/Programs/bcftools-1.10.2/bcftools filter -Oz -i  \
	'QUAL >= 20 & AVG(FMT/GQ)>10 & INFO/DP>=900 & INFO/DP<=2700 & TYPE="snp" & F_MISSING<0.1' \
	>/work/itl4/Austinii/filtered/filter6_201218.vcf.gz \
	2>/work/itl4/Austinii/filtered/filter6_201218.err

#test filter 201028 - biallelic

/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools view -m2 -M2 -v snps -S /dscrhome/itl4/RI_Austinii/filtersamples.txt \
/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/raw_calls_LeafIpo_201009.bcf | \
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools filter -Oz -i  'QUAL >= 20 & AVG(FMT/GQ)>10 & INFO/DP>=900 & INFO/DP<=2700 & TYPE="snp"' \
>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter4_201028.vcf.gz 2>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter4_201028.err

/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools view -m2 -M2 -v snps -S /dscrhome/itl4/RI_Austinii/filtersamples.txt \
        /work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/raw_calls_LeafIpo_201009.bcf | \
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools filter -Oz -i  'QUAL >= 20 & AVG(FMT/GQ)>10 & INFO/DP>=900 & TYPE="snp"' \
        >/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter5_201028.vcf.gz 2>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter5_201028.err

/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools view -m2 -M2 -v snps -S /dscrhome/itl4/RI_Austinii/filtersamples.txt \
        /work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/raw_calls_LeafIpo_201009.bcf | \
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools filter -Oz -i  'QUAL >= 20 & AVG(FMT/GQ)>10 & INFO/DP>=900 & INFO/DP<=2700 & TYPE="snp" & F_MISSING<0.1' \
        >/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter6_201028.vcf.gz 2>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter6_201028.err
        
#Filters_applied_201019
#90 samples total

#filter 1 - common one I've used before 
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools view -S /dscrhome/itl4/RI_Austinii/filtersamples.txt \
/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/raw_calls_LeafIpo_201009.bcf | \
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools filter -Oz -i  'QUAL >= 20 & AVG(FMT/GQ)>10 & INFO/DP>=900 & INFO/DP<=2700 & TYPE="snp"' \
>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter1_201019.vcf.gz 2>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter1_201019.err

#filter 2 - a bit less stringent 
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools view -S /dscrhome/itl4/RI_Austinii/filtersamples.txt \
        /work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/raw_calls_LeafIpo_201009.bcf | \
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools filter -Oz -i  'QUAL >= 20 & AVG(FMT/GQ)>10 & INFO/DP>=900 & TYPE="snp"' \
>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter2_201019.vcf.gz 2>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter2_201019.err

#filter 3 - a bit more stringent - add filter missing
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools view -S /dscrhome/itl4/RI_Austinii/filtersamples.txt \
        /work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/raw_calls_LeafIpo_201009.bcf | \
/dscrhome/itl4/Programs/bcftools-1.10.2/bcftools filter -Oz -i  'QUAL >= 20 & AVG(FMT/GQ)>10 & INFO/DP>=900 & INFO/DP<=2700 & TYPE="snp" & F_MISSING<0.1' \
>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter3_201019.vcf.gz 2>/work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/All_merged/Calls/filter3_201019.err
