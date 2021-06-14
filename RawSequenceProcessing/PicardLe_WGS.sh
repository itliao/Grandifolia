#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=10
#SBATCH --job-name=PicardWGS_A
#SBATCH -p rausherlab,common,scavenger 
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

##########
# lane 1 #
##########
while read i j
do

echo "Processing sample $i ..." 

java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar FixMateInformation \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL1_$i$j.bam
 
java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar CleanSam \
	VALIDATION_STRINGENCY=LENIENT \
	INPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL1_$i$j.bam \
	OUTPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_L1_$i.bam

java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar SortSam \
	INPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_L1_$i.bam \
	OUTPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_sorted_L1_$i.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp 
 
java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
        VALIDATION_STRINGENCY=LENIENT \
        I=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_sorted_L1_$i.bam  \
        O=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/mergeready_L1_$i.bam  \
        RGID=1$i \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=HWJMMBBXX.3$i \
        RGSM=sample$i
	    
done < /hpc/group/biodept/itl4/WGS/WGS_Atab.txt 

##########
# lane 2 #
##########

#updated 201127
while read i j
do

echo "Processing sample $i ..." 

java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar FixMateInformation \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL2_$i$j.bam
 
java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar CleanSam \
	VALIDATION_STRINGENCY=LENIENT \
	INPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL2_$i$j.bam \
	OUTPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_L2_$i.bam

java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar SortSam \
	INPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_L2_$i.bam \
	OUTPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_sorted_L2_$i.bam \
	SORT_ORDER=coordinate \ 
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp 
 
java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_sorted_L2_$i.bam  \
	O=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/mergeready_L2_$i.bam  \
	RGID=2$i \
    RGLB=lib2 \
    RGPL=illumina \
    RGPU=HWJVGBBXX.7$i \
    RGSM=sample$i 
	    
done < /hpc/group/biodept/itl4/WGS/WGS_Btab.txt 

##########
# lane 3 #
##########

while read i j

do
echo "Processing sample $i ..." 

java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar FixMateInformation \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL3_$i$j.bam

java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar CleanSam \
	VALIDATION_STRINGENCY=LENIENT \
	INPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/rawL3_$i$j.bam \
	OUTPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_L3_$i.bam

java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar SortSam \
	INPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_L3_$i.bam \
	OUTPUT=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_sorted_L3_$i.bam \
	SORT_ORDER=coordinate \ 
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp 
 
java -jar /hpc/group/biodept/itl4/Programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
	VALIDATION_STRINGENCY=LENIENT \
	I=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/tmp_clean_sorted_L3_$i.bam  \
	O=/work/itl4/Alignments/NGM/Ipomoea/2020/WGS/Aligned/mergeready_L3_$i.bam  \
	RGID=3$i \
	RGLB=lib3 \
	RGPL=illumina \
    RGPU=HWJVGBBXX.8$i \
    RGSM=sample$i
	    
done < /hpc/group/biodept/itl4/WGS/WGS_Btab.txt 