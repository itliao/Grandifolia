work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/RI_Austinii"
setwd (work_dir)

#visualizing admixture and faststructure output using pophelper
library(pophelper)
library(ggplot2)
library(gridExtra)

path <- work_dir
###############
# group names #
###############
# in order of 
grp <- read.delim("species4.txt", header=F)

##############
# color list #
##############
clist <- list(
  "shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
  "strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
  "oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
  "keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
  "vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
  "muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
  "teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
  "merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
  "funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
  "retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
  "cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
  "cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
  "morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
  "wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
  "krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"))


###################
# admixture files #
###################
fileAll <- paste0(path,"/admixture/admixQ/")
afiles <- list.files(path=fileAll, full.names=T)
alist <- readQ(files=afiles)
#alignK by colors first
alist <- alignK(alist)
length(alist)
#merge by K
length(mergeQ(alist))
alist_merge <- mergeQ(alist)

#######################
# faststructure files #
#######################
fileFS <- paste0(path,"/fastStructure/meanQ/")
fsfiles <- list.files(path=fileFS, full.names=T)
fslist <- readQ(files=fsfiles)
#alignK by colors first
fslist <- alignK(fslist)
length(fslist)
#merge by K
length(mergeQ(fslist))
fslist_merge <- mergeQ(fslist)

###################
# admixture plots #
###################

# merged admixture K2-10 all on one file
p1 <- plotQ(alignK(sortQ(alist_merge[c(1:9)])),imgoutput="join",returnplot=T, exportplot=F,
            showindlab=T,basesize=11, ordergrp = T,grplab=data.frame(lab1=grp$V1),
            sharedindlab=F, clustercol=clist$cb_set3, splab=paste0("K=",sapply(alist_merge[c(1:9)],ncol)))
grid.arrange(p1$plot[[1]])

# merged K=5
# rearranged, bars do not correspond to same individuals as previous plot

p5 <- plotQ(alist_merge[4],returnplot=T, exportplot=F,
            showindlab=T,basesize=11, ordergrp = T, grplab=data.frame(lab1=grp$V1), sortind="all",
            sharedindlab=F, clustercol=clist$cb_set3, splab=paste0("K=",sapply(alist_merge[4],ncol)))
grid.arrange(p5$plot[[1]])

#######################
# faststructure plots #
#######################

# merged faststructure K2-10 all on one file
p2 <- plotQ(alignK(sortQ(fslist_merge[c(1:9)])),imgoutput="join",returnplot=T, exportplot=F,
            showindlab=T,basesize=11, ordergrp = T,grplab=data.frame(lab1=grp$V1),
            sharedindlab=F, clustercol=clist$cb_set3, splab=paste0("K=",sapply(fslist_merge[c(1:9)],ncol)))

grid.arrange(p2$plot[[1]])

#plot merged K=7 (6)
p5 <- plotQ(fslist_merge[6],returnplot=T, exportplot=F,
            showindlab=T,basesize=11, ordergrp = T, grplab=data.frame(lab1=grp$V1), sortind="all",
            sharedindlab=F, clustercol=clist$cb_set3, splab=paste0("K=",sapply(fslist_merge[6],ncol)))
grid.arrange(p5$plot[[1]])

p7 <- plotQ(fslist_merge[5],returnplot=T, exportplot=F,
            showindlab=T,basesize=11, ordergrp = T, grplab=data.frame(lab1=grp$V1), sortind="all",
            sharedindlab=F, clustercol=clist$cb_set3, splab=paste0("K=",sapply(fslist_merge[5],ncol)))
grid.arrange(p7$plot[[1]])
