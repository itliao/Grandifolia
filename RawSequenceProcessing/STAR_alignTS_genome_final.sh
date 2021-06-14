#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=STAR_leaf
#SBATCH -p common,scavenger,rausherlab
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL


while read i
do

echo "Processing sample $i ..." 

/dscrhome/itl4/Programs/STAR-2.7.5c/bin/Linux_x86_64/STAR \
        --runThreadN 24 \
        –-runMode alignReads \
        --genomeDir /dscrhome/itl4/LacGenome/STAR_final_GFF \
        --readFilesCommand zcat \
        --readFilesIn /work/itl4/Austinii/Trim/$i\_R1_paired.fastq.gz /work/itl4/Austinii/Trim/$i\_R2_paired.fastq.gz \
 		--sjdbGTFfile /dscrhome/itl4/LacGenome/annotation/maker.4.2.ipr.sort.gff \
		--sjdbGTFtagExonParentTranscript Parent \
        --sjdbOverhang 149 \
        --quantMode TranscriptomeSAM \
        --twopassMode Basic \
        --outMultimapperOrder Random \
        --outFileNamePrefix /work/itl4/Alignments/STAR/2020/Austinii/Genome/TransCord/GFF/$i \
        --outSAMtype BAM SortedByCoordinate \

done < /dscrhome/itl4/RI_Austinii/IpoSp.txt


