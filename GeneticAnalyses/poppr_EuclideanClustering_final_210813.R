work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/RI_Austinii"
setwd (work_dir)

# Euclidean distance and clustering
library(vcfR)
library(poppr)
library(ape)
library(adegenet)
library(RColorBrewer)
library(ggplot2)
library(dartR)
library(rrBLUP)

###########################################
# read VCF and convert to genlight object #
###########################################

# path and read in VCF
vcf_path <- paste0(work_dir,"/VCF_Filter/filter6_201218.vcf.gz")
ipo.VCF <- read.vcfR(vcf_path)

#file with the VCF sample names and the species associated with the sample
pop.data <- read.table("species5.txt", sep = "\t", header = TRUE)

#check that all samples are in the VCF
all(colnames(ipo.VCF@gt)[-1] == pop.data$sample.id)

#convert VCF into genlight object
gl.ipo <- vcfR2genlight(ipo.VCF)

#####################################
# Euclidean distance and clustering #
#####################################
#calculate Euclidean distance, number of different alleles
eudist <- poppr::bitwise.dist(gl.ipo, euclidean = TRUE)
eudist.matrix <- as.matrix(eudist)
eudist.df <- as.data.frame(eudist.matrix)
# write_csv(eudist.df, "poppr_adegenet/EuclideanDistance_210126.csv")

# clustering methods:
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA)
# "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

eucluster <- hclust(eudist, method = "average")
plot(eucluster,hang = -1, main = "Cluster Dendrogram", ylab = "Height")

# save pdf of clustering output (test several types)
pdf("poppr_adegenet/SNPdistance_cluster_210126.pdf")
d <- eudist
hcorr <- hclust(d, method="complete")
plot(hcorr, hang = -1, main="Euclidean")
hcorr <- hclust(d, method="ward.D")
plot(hcorr, hang = -1, main="Euclidean")
hcorr <- hclust(d, method="ward.D2")
plot(hcorr, hang = -1, main="Euclidean")
hcorr <- hclust(d, method="average")
plot(hcorr, hang = -1, main="Euclidean")
hcorr <- hclust(d, method="mcquitty")
plot(hcorr, hang = -1, main="Euclidean")
dev.off()

# heatmap of Euclidean distances
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(eudist.matrix, col = col, scale = c("row"))
heatmap(dist.default.matrix, col = col)
