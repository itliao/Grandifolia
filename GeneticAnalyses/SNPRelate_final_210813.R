work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/RI_Austinii"
setwd (work_dir)

# Using SNPRelate to plot PCA and identity-by-state dendrogram
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(plotly)

###################################################
# Format conversion from VCF files using SeqArray #
###################################################
# Good for large-scale whole-exome and whole-genome sequencing variant data
# Don't need to rerun if already converted

# the VCF file
vcf.fn.1 <- paste0(work_dir,"/VCF_Filter/filter6_201218.vcf.gz")
# convert, save in "tmp.gds" with the default lzma compression algorithm
seqVCF2GDS(vcf.fn.1, "newfilter6_redo.gds")

################
# PCA analysis #
################
# open  GDS file
genofile <- seqOpen("newfilter6_redo.gds")

# take out sample id
head(samp.id <- seqGetData(genofile, "sample.id"))

# get what population/species each sample belongs to
pop_code <- scan("species3.txt", what=character())

########################
# LD-based SNP pruning #
########################
set.seed(1000)

# Rifkin et al., 2019 MolEco: ld.threshold=0.3

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.3)
snpset.id <- unlist(unname(snpset))

####################
# Run PCA analysis #
####################
#The functions in SNPRelate for PCA include:
# 1) calculating the genetic covariance matrix from genotypes
# 2) computing the correlation coefficients between sample loadings and genotypes for each SNP
# 3) calculating SNP eigenvectors (loadings)
# 4) estimating the sample loadings of a new dataset from specified SNP eigenvectors.

pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2) 

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# the order of sample IDs is as the same as population codes
sample.id <- pca$sample.id
head(cbind(samp.id, pop_code))

#make a dataframe for first 2 eigenvectors
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

# write.table(tab, file = "filter6_EV1_EV2_201217.txt", sep = "\t")

# Plot PC1 and PC2
plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), xlab="eigenvector 1", ylab="eigenvector 2")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))

# for use in plot_ly
fig <- plot_ly(tab, x = ~EV1, y= ~EV2, color = ~pop)
fig

# Plot the principal component pairs for the first four PCs:
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)

###########################################
# Identity-By-State Analysis & Dendrogram #
###########################################
# cluster analysis on the n×n matrix of genome-wide IBS pairwise distances
# determine the groups by a permutation score

set.seed(100)
ipsp <- as.factor(pop_code)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2))

# determine groups of individuals by population information
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))

# plot dendrogram
plot(rv2$dendrogram, leaflab="perpendicular", main="Dendrogram_Ipomoea")
legend("topright", legend=levels(ipsp), col=1:nlevels(ipsp), pch=19, ncol=4)
