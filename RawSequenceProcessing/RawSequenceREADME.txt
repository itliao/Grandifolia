Raw sequence read processing

Taking raw sequence reads from leaf transcriptomes from various Ipomoea species and calling SNP variants by aligning sequences to the lab Ipomoea lacunosa draft genome with a draft annotation.

For 2 Ipomoea leucantha samples, leaf transcriptomes were not sequenced, but we have whole genome sequences from another project. These were processed slightly differently (see below)

Sequence files - NCBI SRA accession numbers:
(1) PRJNA735523: Ipomoea samples ("austinii", cynanchifolia, grandifolia, lacunosa (purple), leucantha, ramosissima, splendor-sylvae, tenuissima)
(2) XXX: Ipomoea lacunosa and Ipomoea cordatotriloba from Rifkin et al., 2019 Molecular Ecology
(3) PRJNA732922: Ipomoea leucantha samples

###########################
# For leaf transcriptomes #
###########################
Step 1: Trim with Trimmomatic v. 0.39
IpoInd.txt 
IpoSpAC.txt 
Trim_final.sh

Step 2: Align with STAR v.2.7.5c with draft genome and annotation in gff
For draft genome and gff annotation, contact Mark D. Rausher (mrausher@duke.edu)
IpoSp.txt 
STAR_alignTS_genome.sh

Step 3: Clean sequences with Picard Tools 2.19.2-1-g0d1e881
IpoSp_1.txt 
IpoSp_2.txt
IpoSp_3.txt 

Picard1.sh
Picard2.sh
Picard3.sh

Step 4: Merge files (if more than 1 lane of the same sample was sequenced; (samtools v.1.9)
IpoSp_merge.txt
MergeClean.sh

Step 5: Validate BAM before variant calling
IpoSp_merge.txt
IpoSp_AC_clean.txt
Validate_BAM_RNALeaf.sh

Step 6: Variant calling with mpileup (includes Ipomoea leucantha WGS)
bamlist_IpoLeaf_2Le.txt
VariantCall_Ipo_final.sh

Step 7: Filtering variants with only SNPs
FilterCalls_final.sh

Final output: filter6_201218.vcf.gz

#############################################
# For whole genome sequences (I. leucantha) #
#############################################
Step 1: Align with NextGenMap v.0.5.5 with draft genome
WGS_Ind_A.txt 
WGS_Ind_B.txt 
AlignLe_WGS.sh

Step 2: Clean sequences with Picard Tools 2.19.2-1-g0d1e881
WGS_Atab.txt 
WGS_Btab.txt 
PicardLe_WGS.sh

Step 3: Merge files (samtools v.1.9)
WGS_Ind_loop.txt
Merge_WGS_loop.sh

Step 4: Validate BAM before variant calling
ValidateLe_BAM.sh
