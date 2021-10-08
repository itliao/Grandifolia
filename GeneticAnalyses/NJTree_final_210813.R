work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/RI_Austinii"
setwd (work_dir)

# building neighbor-joining tree
library(ape)
library(ade4)
library(adegenet)
library(pegas)
library(vcfR)
library(phangorn)

################################
# Convert VCF into DNAbin file #
################################
# these steps can take a while
# don't need to redo if already have DNAbin file

# read in VCF file
JG1v <- read.vcfR(paste0(work_dir,"/VCF_Filter/filter6_201218.vcf.gz"), verbose = FALSE)
# convert VCF
JG1Dbin <- vcfR2DNAbin(JG1v, consensus = TRUE, extract.haps = FALSE, verbose = FALSE)
JG1Dbin
# write DNAbin file (to skip these steps)
write.dna(JG1Dbin, file = 'filter6_DNAbin.fasta', format = 'fasta')

##########################
# NJ tree using phangorn #
##########################

# read in DNAbin file
testDNA <- fasta2DNAbin("filter6_DNAbin.fasta")
# changing data class from DNAbin to phyDat
JG1pD <- as.phyDat(testDNA)

# calculate distance and build tree
dm <- dist.logDet(JG1pD) 
tree <- NJ(dm)

###########################
# test with 10 bootstraps #
###########################
set.seed(123)
NJtrees <- bootstrap.phyDat(JG1pD, FUN=function(x)NJ(dist.logDet(x)), bs=10)
treeNJ <- plotBS(tree, NJtrees, "phylogram")

# choose outgroup to root tree(s)
tree.2 <- root(NJtrees ,outgroup=c("sampleIpoSpl_03_S63"),resolve.root=TRUE)
tree.3 <- root(tree ,outgroup=c("sampleIpoSpl_03_S63"),resolve.root=TRUE)

# plot trees
treeNJroot <- plotBS(tree.3, NJtrees, "phylogram")
densiPhy1 <- densiTree(tree.2, type="phylogram", col="blue")
#write.tree(tree.2, "phangorn_NJ_b10_root_201218.tre")

##################################
# final run with 1000 bootstraps #
##################################

NJtrees1000 <- bootstrap.phyDat(JG1pD, FUN=function(x)NJ(dist.logDet(x)), bs=1000)
treeNJ <- plotBS(tree, NJtrees1000, "phylogram")
#write.tree(NJtrees1000, "NJtrees1000_201219.tre")

tree2.1000 <- root(NJtrees1000 ,outgroup=c("sampleIpoSpl_03_S63"),resolve.root=TRUE)
tree3.1000 <- root(tree ,outgroup=c("sampleIpoSpl_03_S63"),resolve.root=TRUE)
treeNJroot.1000 <- plotBS(tree3.1000, NJtrees1000, "phylogram")

densiTree(tree2.1000, optim=TRUE, type="phylogram", col="blue")
#write.tree(tree2.1000, "phangorn_NJ_b1000_root_201218.tre")