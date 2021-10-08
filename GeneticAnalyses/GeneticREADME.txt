Genetic differentiation analyses

Used three complementary analyses to test for genetic differentiation between Ipomoea "austinii" and Ipomoea grandifolia

(1) fastSTRUCTURE and ADMIXTURE (slight different algorithms, but similar results)
(2) PCA
(3) Neighbor joining tree
(4) Euclidean distance tree

Filtered VCF file for all analyses: filter6_201218.vcf.gz

########################################################
# Convert VCF to files for fastSTRUCTURE and ADMIXTURE #
########################################################
use Plink v.1.9
Plink_bed.sh

#################
# fastSTRUCTURE #
#################
requires python 2.7, so used compatible virtual environment
randomNumbers_fastStructure.txt 
faststructure_final.sh

for best K, summary of model complexity output:
chooseK_complexity.txt

#############
# ADMIXTURE #
#############
requires python 3+, so used compatible virtual environment
randomNumber_admixture.txt - random number seed for each run; create a separate folder for each run and run the script within the specific folder
admixture_final.sh

for best K, summary of cross validation output:
admixtureCV.txt

##################################################
# Visualizing fastSTRUCTURE and ADMIXTURE output #
##################################################
use pophelper v2.3.1 in R (Figures 5, S2)
for species list: species4.txt
for fastSTRUCTURE files: (admixQ)
for ADMIXTURE files: (meanQ)
R file: pophelper_final_210813.R

#######
# PCA #
#######
use SNPRelate in R (Figures 5, S3, S5)
convert to gds file: newfilter6_redo.gds
population/species codes: species3.txt
R file: SNPRelate_final_210813.R

#########################
# Neighbor-joining tree #
#########################
use phangorn in R (Figure 5)
convert to DNAbin file to convert to phyDat file from VCF: filter6_DNAbin.fasta
R file: NJTree_final_210813.R

##############################
# Euclidean distance cluster #
##############################
use poppr v.2.8.7 in R (Figure S4)
R file: poppr_EuclideanClustering_final_210813.R
