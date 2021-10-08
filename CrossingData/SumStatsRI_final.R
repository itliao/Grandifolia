work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/ReproductiveIsolation/"
setwd (work_dir)

library(lme4) # for linear mixed-effect models
library(glmmTMB) # for models including zero-inflation 
library(DHARMa) # for model fit
library(emmeans) # for model comparisons
library(lattice)
library(car)
library(blme)

library(plyr) #for summarizing data
library(tidyverse)
library(gplots) 
library(scales) # for transparency
library(ggpubr)
library(RColorBrewer)


# define a function that runs multiple model fitting tests from DHARMa
DHARMa_test_model_fit <- function(model) {
  print(summary(model))
  hist(resid(model))
  sim_res <- simulateResiduals(model)
  testDispersion(sim_res) # simulation test for over/under dispersion
  print(testOutliers(sim_res))
  plotQQunif(sim_res)
  plotResiduals(sim_res)
  print(testZeroInflation(sim_res))
  print(Anova(model))
}

#input seed file
path <- work_dir

cross_file <- "RIseeds_final.txt"
CL_KLO_file <- "Ipo_RI_survey.csv"
#file modified to change "S" into the associated species to combine the within-individual
#with the within species crosses
#cross_file <- "RIseeds_final2.txt"

crossData <- read.table(cross_file, sep="\t", header=TRUE, stringsAsFactor=FALSE)
CLdata <- read.csv(CL_KLO_file, stringsAsFactors = TRUE, na.strings = c("NA", ""))

#############################################################
# Summary of fruited and seeds/fruit for each species cross #
# Tables 2 and S5                                           #
#############################################################

# use ddply from plyr 
# Table S5
# summarize seed/fruit by 2 groups, FemaleSp and MaleSp, unweighted average
seedData <- ddply(crossData, c("FemaleSp", "MaleSp"), summarise,
               N    = length(Number_of_seeds),
               mean = mean(Number_of_seeds), 
               sd   = sd(Number_of_seeds),
               se   = sd / sqrt(N)
)
seedData

seedDF <- as.data.frame(seedData)
# write.table(seedDF, file = "RIseed_meanSummaryCombo_201211.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

#fruited data, summarized by two groups, FemaleSp and MaleSp, unweighed average
fruitData <- ddply(crossData, c("FemaleSp", "MaleSp"), summarise,
                  N    = length(Fruited),
                  mean = mean(Fruited), 
                  sd   = sd(Fruited),
                  se   = sd / sqrt(N)
)
fruitData

fruitDF <- as.data.frame(fruitData)
# write.table(fruitDF, file = "RIfruit_meanSummaryCombo_201211.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# seed data summarized by the following variables 
# find mean seeds per fruit for each unique population pair
# used for statistical analyses below
seedFamily <- ddply(crossData, c("Female_parent", "Male_parent", "FemaleSp", "MaleSp", "Cross", "CrossFamily"), summarise,
                  N    = length(Number_of_seeds),
                  mean = mean(Number_of_seeds), 
                  sd   = sd(Number_of_seeds),
                  se   = sd / sqrt(N)
)
seedFamily

seedFamilyDF <- as.data.frame(seedFamily)
#write.table(seedFamilyDF, file = "RIseedFamily_meanSummaryCombo_201211.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# then find the average of the species level crosses; weighted average of the seeds
# Table 2
seedFamilyAvg <- ddply(seedFamily, c("FemaleSp", "MaleSp", "Cross"), summarise,
                    N    = length(mean),
                    sd   = sd(mean),
                    se   = sd / sqrt(N),
                    mean = mean(mean)
)
seedFamilyAvg

seedFamilyAvg <- as.data.frame(seedFamilyAvg)
#write.table(seedFamilyAvg, file = "RIseedFamilyAvg_meanSummaryCombo_201228.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# Kate Ostevik's CL data to calculate the unweighted average of seeds per fruit
# Table S5
CLseedData <- ddply(CLdata, c("mom.sp", "dad.sp"), summarise,
                  N    = length(n.seeds),
                  mean = mean(n.seeds), 
                  sd   = sd(n.seeds),
                  se   = sd / sqrt(N)
)
CLseedData

# seed per fruit data, average of the individual crosses

CLseed <- ddply(CLdata, c("mom", "dad", "mom.sp", "dad.sp", "cross.ind", "r.sp"), summarise,
                    N    = length(n.seeds),
                    mean = mean(n.seeds), 
                    sd   = sd(n.seeds),
                    se   = sd / sqrt(N)
)
CLseed

CLseedDF <- as.data.frame(CLseed)
#write.table(CLseedDF, file = "RIseedCLKate_meanSummary_201207.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# seed per fruit data, average of the individual crosses
# find mean seeds per fruit for each unique population pair
# Kate's CL data to calculate the weighted average of seeds per fruit
CLseedAvg <- ddply(CLseedDF, c("mom.sp", "dad.sp", "r.sp"), summarise,
                N    = length(mean),
                sd   = sd(mean),
                se   = sd / sqrt(N),
                mean = mean(mean) 
)
CLseedAvg

CLseedAvgDF <- as.data.frame(CLseedAvg)
# write.table(CLseedAvgDF, file = "RIseedCLKate_meanWeightedSummary_210125.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

#########################################################################
# Subsetting data for downstream pairwise comparisons, plots, and tests #
#########################################################################

# subset mean seed per fruit data by cross type 
AA <- subset(seedFamilyDF, FemaleSp=="Au" & MaleSp=="Au")
AS <- subset(seedFamilyDF, FemaleSp=="Au" & MaleSp=="S")
AX <- subset(seedFamilyDF, FemaleSp=="Au" & MaleSp=="X")
GG <- subset(seedFamilyDF, FemaleSp=="G" & MaleSp=="G")
GS <- subset(seedFamilyDF, FemaleSp=="G" & MaleSp=="S")
GX <- subset(seedFamilyDF, FemaleSp=="G" & MaleSp=="X")
LeLe <- subset(seedFamilyDF, FemaleSp=="Le" & MaleSp=="Le")
LeS <- subset(seedFamilyDF, FemaleSp=="Le" & MaleSp=="S")
LeX <- subset(seedFamilyDF, FemaleSp=="Le" & MaleSp=="X")
LaLa <- subset(seedFamilyDF, FemaleSp=="La" & MaleSp=="La")
LaS <- subset(seedFamilyDF, FemaleSp=="La" & MaleSp=="S")
LaX <- subset(seedFamilyDF, FemaleSp=="La" & MaleSp=="X")
CC <- subset(seedFamilyDF, FemaleSp=="C" & MaleSp=="C")
CS <- subset(seedFamilyDF, FemaleSp=="C" & MaleSp=="S")
CX <- subset(seedFamilyDF, FemaleSp=="C" & MaleSp=="X")

AG <- subset(seedFamilyDF, FemaleSp=="Au" & MaleSp=="G")
GA <- subset(seedFamilyDF, FemaleSp=="G" & MaleSp=="Au")
ALe <- subset(seedFamilyDF, FemaleSp=="Au" & MaleSp=="Le")
LeA <- subset(seedFamilyDF, FemaleSp=="Le" & MaleSp=="Au")
ALa <- subset(seedFamilyDF, FemaleSp=="Au" & MaleSp=="La")
LaA <- subset(seedFamilyDF, FemaleSp=="La" & MaleSp=="Au")
AC <- subset(seedFamilyDF, FemaleSp=="Au" & MaleSp=="C")
CA <- subset(seedFamilyDF, FemaleSp=="C" & MaleSp=="Au")
GLe <- subset(seedFamilyDF, FemaleSp=="G" & MaleSp=="Le")
LeG <- subset(seedFamilyDF, FemaleSp=="Le" & MaleSp=="G")
GLa <- subset(seedFamilyDF, FemaleSp=="G" & MaleSp=="La")
LaG <- subset(seedFamilyDF, FemaleSp=="La" & MaleSp=="G")
GC <- subset(seedFamilyDF, FemaleSp=="G" & MaleSp=="C")
CG <- subset(seedFamilyDF, FemaleSp=="C" & MaleSp=="G")
LeLa <- subset(seedFamilyDF, FemaleSp=="Le" & MaleSp=="La")
LaLe <- subset(seedFamilyDF, FemaleSp=="La" & MaleSp=="Le")
LeC <- subset(seedFamilyDF, FemaleSp=="Le" & MaleSp=="C")
CLe <- subset(seedFamilyDF, FemaleSp=="C" & MaleSp=="Le")

#subset original crossing data by cross type
AAc <- subset(crossData, FemaleSp=="Au" & MaleSp=="Au")
ASc <- subset(crossData, FemaleSp=="Au" & MaleSp=="S")
AXc <- subset(crossData, FemaleSp=="Au" & MaleSp=="X")
GGc <- subset(crossData, FemaleSp=="G" & MaleSp=="G")
GSc <- subset(crossData, FemaleSp=="G" & MaleSp=="S")
GXc <- subset(crossData, FemaleSp=="G" & MaleSp=="X")
LeLec <- subset(crossData, FemaleSp=="Le" & MaleSp=="Le")
LeSc <- subset(crossData, FemaleSp=="Le" & MaleSp=="S")
LeXc <- subset(crossData, FemaleSp=="Le" & MaleSp=="X")
LaLac <- subset(crossData, FemaleSp=="La" & MaleSp=="La")
LaSc <- subset(crossData, FemaleSp=="La" & MaleSp=="S")
LaXc <- subset(crossData, FemaleSp=="La" & MaleSp=="X")
CCc <- subset(crossData, FemaleSp=="C" & MaleSp=="C")
CSc <- subset(crossData, FemaleSp=="C" & MaleSp=="S")
CXc <- subset(crossData, FemaleSp=="C" & MaleSp=="X")

AGc <- subset(crossData, FemaleSp=="Au" & MaleSp=="G")
GAc <- subset(crossData, FemaleSp=="G" & MaleSp=="Au")
ALec <- subset(crossData, FemaleSp=="Au" & MaleSp=="Le")
LeAc <- subset(crossData, FemaleSp=="Le" & MaleSp=="Au")
ALac <- subset(crossData, FemaleSp=="Au" & MaleSp=="La")
LaAc <- subset(crossData, FemaleSp=="La" & MaleSp=="Au")
ACc <- subset(crossData, FemaleSp=="Au" & MaleSp=="C")
CAc <- subset(crossData, FemaleSp=="C" & MaleSp=="Au")
GLec <- subset(crossData, FemaleSp=="G" & MaleSp=="Le")
LeGc <- subset(crossData, FemaleSp=="Le" & MaleSp=="G")
GLac <- subset(crossData, FemaleSp=="G" & MaleSp=="La")
LaGc <- subset(crossData, FemaleSp=="La" & MaleSp=="G")
GCc <- subset(crossData, FemaleSp=="G" & MaleSp=="C")
CGc <- subset(crossData, FemaleSp=="C" & MaleSp=="G")
LeLac <- subset(crossData, FemaleSp=="Le" & MaleSp=="La")
LaLec <- subset(crossData, FemaleSp=="La" & MaleSp=="Le")
LeCc <- subset(crossData, FemaleSp=="Le" & MaleSp=="C")
CLec <- subset(crossData, FemaleSp=="C" & MaleSp=="Le")

###################################
# Plot for all the data (Fig. S1) #
###################################

#for ALL the crosses
#stripchart(mean ~ Cross, data=seedFamilyDF, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
#           pch = 16, cex = 1.5, lty = 1)
mosaicplot(table(crossData$Cross, crossData$Number_of_seeds), las=2, main=NULL, ylab = "Number of seeds per fruit", xlab = "Cross type")

###############################
# remove crosses with 0 seeds #
###############################
# first remove crosses with 0 seeds - for testing if Au and G interchangable
oneSeedCross <- subset(crossData, Number_of_seeds != 0)

# then find average seed per fruit for unique population pairs with at least 1 seed per fruit
seed1Family <- ddply(oneSeedCross, c("Female_parent", "Male_parent", "FemaleSp", "MaleSp", "Cross", "CrossFamily"), summarise,
                    N    = length(Number_of_seeds),
                    mean = mean(Number_of_seeds), 
                    sd   = sd(Number_of_seeds),
                    se   = sd / sqrt(N)
)
seed1Family

seed1FamilyDF <- as.data.frame(seed1Family)
#write.table(seed1FamilyDF, file = "RIseedFamily_meanSummary_1seed_Combo_210225.txt")

# subset average seed per fruit for all fruited
AA1 <- subset(seed1Family, FemaleSp=="Au" & MaleSp=="Au")
GG1 <- subset(seed1Family, FemaleSp=="G" & MaleSp=="G")
AG1 <- subset(seed1Family, FemaleSp=="Au" & MaleSp=="G")
GA1 <- subset(seed1Family, FemaleSp=="G" & MaleSp=="Au")

ALe1 <- subset(seed1Family, FemaleSp=="Au" & MaleSp=="Le")
LeA1 <- subset(seed1Family, FemaleSp=="Le" & MaleSp=="Au")
ALa1 <- subset(seed1Family, FemaleSp=="Au" & MaleSp=="La")
LaA1 <- subset(seed1Family, FemaleSp=="La" & MaleSp=="Au")
AC1 <- subset(seed1Family, FemaleSp=="Au" & MaleSp=="C")
CA1 <- subset(seed1Family, FemaleSp=="C" & MaleSp=="Au")
GLe1 <- subset(seed1Family, FemaleSp=="G" & MaleSp=="Le")
LeG1 <- subset(seed1Family, FemaleSp=="Le" & MaleSp=="G")
GLa1 <- subset(seed1Family, FemaleSp=="G" & MaleSp=="La")
LaG1 <- subset(seed1Family, FemaleSp=="La" & MaleSp=="G")
GC1 <- subset(seed1Family, FemaleSp=="G" & MaleSp=="C")
CG1 <- subset(seed1Family, FemaleSp=="C" & MaleSp=="G")

LeLe1 <- subset(seed1Family, FemaleSp=="Le" & MaleSp=="Le")
LeLa1 <- subset(seed1Family, FemaleSp=="Le" & MaleSp=="La")
LaLe1 <- subset(seed1Family, FemaleSp=="La" & MaleSp=="Le")
LaS1 <- subset(seed1Family, FemaleSp=="La" & MaleSp=="S")

#########################################
# austinii and grandifolia only crosses #
# Table 1, Fig. 2B-D                    #
#########################################
#not useful, but needed for good merging of files
AA["Au_G"] <- AA$FemaleSp
AA["Sp_Male_Female"] <- "Male"
GG["Au_G"] <- GG$MaleSp
GG["Sp_Male_Female"] <- "Female"

AG["Au_G"] <- AG$FemaleSp
AG["Sp_Male_Female"] <- "Male"
GA["Au_G"] <- GA$MaleSp
GA["Sp_Male_Female"] <- "Female"
#subsetted data
AGseed <- rbind(AA,AG,GA,GG)
AGcross <- rbind(AAc,AGc,GAc,GGc)

AGseed$Cross <- factor(AGseed$Cross, levels = c("AuxAu", "AuxG", "GxAu", "GxG"))
AGcross$Cross <- factor(AGcross$Cross, levels = c("AuxAu", "AuxG", "GxAu", "GxG"))

# plot data
#stripchart(mean ~ Cross, data=AGseed, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
#           pch = 16, cex = 1.5, lty = 1)

# Fig. 2B
mosaicplot(table(AGcross$Cross, AGcross$Number_of_seeds), las=2, main=NULL, ylab = "Number of seeds per fruit", xlab = "Cross type")

# test models of mean seed number for Au and G
# with zero-inflation, model chosen: AG_S_model2 <- glmmTMB(mean ~ Cross + (1 | Female_parent), ziformula=~Cross, data = AGseed)
# without zero-inflation, model chosen: AG_simple <- lm(mean ~ Cross, data = AGseed)
# Table 1, all crosses
AG_simple <- lm(mean ~ Cross, data = AGseed)
AG_S_model1 <- glmmTMB(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), ziformula=~Cross, data = AGseed)
AG_S_model2 <- glmmTMB(mean ~ Cross + (1 | Female_parent), ziformula=~Cross, data = AGseed)
AG_S_model3 <- glmmTMB(mean ~ (1 | Female_parent), ziformula=~Cross, data = AGseed)
AG_S_model4 <- glmmTMB(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), ziformula=~Cross, family = tweedie, data = AGseed)
AG_S_model5 <- glmmTMB(mean ~ Cross + (1 | Female_parent), ziformula=~Cross, family = tweedie, data = AGseed)
AG_S_model6 <- glmmTMB(mean ~ (1 | Female_parent), ziformula=~Cross, family = tweedie, data = AGseed)
AG_S_model7 <- lmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = AGseed)
AG_S_model8 <- lmer(mean ~ Cross + (1 | Female_parent), data = AGseed)
AG_S_model9 <- lmer(mean ~ (1 | Female_parent), data = AGseed)
anova(AG_S_model2, AG_S_model3, AG_S_model4, AG_S_model5, AG_S_model6)
anova(AG_S_model8, AG_S_model9, AG_simple)

# zero-inflation model
DHARMa_test_model_fit(AG_S_model2)
emmeans(AG_S_model2, pairwise ~ Cross)
emmeans(AG_S_model2, pairwise ~ Cross, adjust="none")
plot(emmeans(AG_S_model2, pairwise ~ Cross))
contrast(emmeans(AG_S_model2, "Cross"), list(AvG=c(1,0,0,-1), H1vH2=c(0,1,-1,0), AGvH = c(1,-1,-1,1)))
contrast(emmeans(AG_S_model2, "Cross"), list(AvG=c(1,0,0,-1), H1vH2=c(0,1,-1,0), AGvH = c(1,-1,-1,1)), adjust = "fdr")

# no zero-inflation model
DHARMa_test_model_fit(AG_simple)
emmeans(AG_simple, pairwise ~ Cross)
emmeans(AG_simple, pairwise ~ Cross, adjust="none")
plot(emmeans(AG_simple, pairwise ~ Cross))
contrast(emmeans(AG_simple, "Cross"), list(AvG=c(1,0,0,-1), H1vH2=c(0,1,-1,0), AGvH = c(1,-1,-1,1)))
contrast(emmeans(AG_simple, "Cross"), list(AvG=c(1,0,0,-1), H1vH2=c(0,1,-1,0), AGvH = c(1,-1,-1,1)), adjust = "fdr")

# went with zero-inflation model, but both give similar results
estimatesAG_all <- emmeans(AG_S_model2, pairwise ~Cross, type="response")[[1]] %>% data.frame()

# Fig. 2D
AGall_comp <- ggplot(AGseed, aes(x=Cross, y=mean)) + ylim(0, 4) +
  geom_jitter(aes(fill = Cross, color = Cross), size = 2, alpha = 0.7,
              position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("#A8DDB5","#7BCCC4","#43A2CA","#0868AC")) +
  scale_color_manual(values = c("#A8DDB5","#7BCCC4","#43A2CA","#0868AC")) +
  geom_pointrange(aes(x = Cross, y = emmean, ymin = lower.CL, ymax = upper.CL, color=Cross), 
                  data = estimatesAG_all,
                  position = position_nudge(x=0.25)) +
  labs(x = "Cross", y = "Mean number of seeds/fruit") +
  theme_bw()
AGall_comp


#~~~~~~~~~~ crosses >= 1 seed; 0 seeds removed ~~~~~~~~~~#

AA_GG1 <- rbind(AA1,GA1,AG1,GG1)

AA_GG1$Cross <- factor(AA_GG1$Cross, levels = c("AuxAu", "AuxG", "GxAu", "GxG"))

# plot data
stripchart(mean ~ Cross, data=AA_GG1, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AA_GG1_simple <- lm(mean ~ Cross, data = AA_GG1)
AA_GG1_simple <- lm(mean ~ Cross, data = AA_GG1)
AA_GG1_S_model1 <- lmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = AA_GG1)
AA_GG1_S_model2 <- lmer(mean ~ Cross + (1 | Female_parent), data = AA_GG1)
AA_GG1_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AA_GG1)
AA_GG1_S_model4 <- blmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = AA_GG1)
AA_GG1_S_model5 <- blmer(mean ~ Cross + (1 | Female_parent), data = AA_GG1)
AA_GG1_S_model6 <- blmer(mean ~ (1 | Female_parent), data = AA_GG1)
anova(AA_GG1_S_model4, AA_GG1_S_model5, AA_GG1_S_model6, AA_GG1_simple)
Anova(AA_GG1_S_model5)
DHARMa_test_model_fit(AA_GG1_simple)
emmeans(AA_GG1_simple, pairwise ~ Cross)
emmeans(AA_GG1_simple, pairwise ~ Cross, adjust="none")
plot(emmeans(AA_GG1_simple, pairwise ~ Cross))
contrast(emmeans(AA_GG1_simple, "Cross"), list(AvG=c(1,0,0,-1), AvH=c(2,-1,-1,0), GvH=c(0,-1,-1,2), H1vH2=c(0,1,-1,0), AGvH = c(1,-1,-1,1)))

estimatesAG <- emmeans(AA_GG1_simple, pairwise ~Cross, type="response")[[1]] %>% data.frame()

# Fig. 2C
AAGG_comp <- ggplot(AA_GG1, aes(x=Cross, y=mean)) + ylim(0, 4.1) +
  geom_jitter(aes(fill = Cross, color = Cross), size = 2, alpha = 0.7,
              position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("#A8DDB5","#7BCCC4","#43A2CA","#0868AC")) +
  scale_color_manual(values = c("#A8DDB5","#7BCCC4","#43A2CA","#0868AC")) +
  geom_pointrange(aes(x = Cross, y = emmean, ymin = lower.CL, ymax = upper.CL, color=Cross), 
                  data = estimatesAG,
                  position = position_nudge(x=0.25)) +
  labs(x = "Cross", y = "Mean number of seeds/fruit") +
  theme_bw()
AAGG_comp


############################################
# pie chart and bar chart of fruit success #
# Fig. 2A and S2                           #
############################################
pdf("PieChartSummary_ComboIntraSp_201211.pdf")

t <- table(crossData$Cross, crossData$Fruited) %>% as.data.frame.matrix()
t <- tibble::rownames_to_column(t)
names(t) <- c("cross","fail", "success")
t
t_fruitSuccess <- t %>% mutate(sum = fail + success, percent = (success/sum)*100)
par(mfrow=c(2,4))
par(mar=c(2,2,2,2))
pie(c(t[1,2],t[1,3]), label = c("fail","success"), main="AuxAu")
pie(c(t[2,2],t[2,3]), label = c("fail","success"), main="AuxC") 
pie(c(t[3,2],t[3,3]), label = c("fail","success"), main="AuxG")
pie(c(t[4,2],t[4,3]), label = c("fail","success"), main="AuxLa")
pie(c(t[5,2],t[5,3]), label = c("fail","success"), main="AuxLe")
pie(c(t[6,2],t[6,3]), label = c("fail","success"), main="AuxS")
pie(c(t[7,2],t[7,3]), label = c("fail","success"), main="AuxX")
pie(c(t[8,2],t[8,3]), label = c("fail","success"), main="CxAu")
par(mfrow=c(1,1))

par(mfrow=c(2,4))
par(mar=c(2,2,2,2))
pie(c(t[9,2],t[9,3]), label = c("fail","success"), main="CxG")
pie(c(t[10,2],t[10,3]), label = c("fail","success"), main="CxLe")
pie(c(t[11,2],t[11,3]), label = c("fail","success"), main="CxS")
pie(c(t[12,2],t[12,3]), label = c("fail","success"), main="CxX")
pie(c(t[13,2],t[13,3]), label = c("fail","success"), main="GxAu") 
pie(c(t[14,2],t[14,3]), label = c("fail","success"), main="GxC")
pie(c(t[15,2],t[15,3]), label = c("fail","success"), main="GxG")
pie(c(t[16,2],t[16,3]), label = c("fail","success"), main="GxLa")

par(mfrow=c(2,4))
par(mar=c(2,2,2,2))
pie(c(t[17,2],t[17,3]), label = c("fail","success"), main="GxLe")
pie(c(t[18,2],t[18,3]), label = c("fail","success"), main="GxS") 
pie(c(t[19,2],t[19,3]), label = c("fail","success"), main="GxX")
pie(c(t[20,2],t[20,3]), label = c("fail","success"), main="LaxAu")
pie(c(t[21,2],t[21,3]), label = c("fail","success"), main="LaxG")
pie(c(t[22,2],t[22,3]), label = c("fail","success"), main="LaxLe")
pie(c(t[23,2],t[23,3]), label = c("fail","success"), main="LaxS")
pie(c(t[24,2],t[24,3]), label = c("fail","success"), main="LaxX")

par(mfrow=c(2,4))
par(mar=c(2,2,2,2))
pie(c(t[25,2],t[25,3]), label = c("fail","success"), main="LexAu")
pie(c(t[26,2],t[26,3]), label = c("fail","success"), main="LexC") 
pie(c(t[27,2],t[27,3]), label = c("fail","success"), main="LexG")
pie(c(t[28,2],t[28,3]), label = c("fail","success"), main="LexLa")
pie(c(t[29,2],t[29,3]), label = c("fail","success"), main="LexLe")
pie(c(t[30,2],t[30,3]), label = c("fail","success"), main="LexS") 
pie(c(t[31,2],t[31,3]), label = c("fail","success"), main="LexX")

dev.off()

#Fruit success of several crosses (Fig. SX)
subset_t <- t_fruitSuccess[c(2:5,8:9,13:14,16:17,20:21,25,27),c(1:5)]
subset_t$cross <- factor(subset_t$cross, 
                         levels = c("AuxG", "GxAu",
                                    "AuxLa", "LaxAu", "GxLa", "LaxG",
                                    "AuxLe", "LexAu", "GxLe","LexG",
                                    "AuxC", "CxAu", "GxC","CxG"))
barFruit <-ggplot(subset_t, aes(x=cross, y=percent)) +
  geom_bar(stat="identity", fill="black") +
  geom_text(aes(label=sprintf("%.1f", percent)), vjust=1.6, color="white", size=3.5) +
  xlab("cross") + ylab("percent fruit success") +
  theme_bw()
barFruit

pdf("barchart_fruitSuccess_210822.pdf")
barFruit
dev.off()

#############################
# for AxG and A/GxSpecies A #
# Figure 3, Table 3         #
#############################

#~~~~~~~~~ crosses with at least 1 seed (0 seeds removed) ~~~~~~~~~#
AG_CLaLe1Seed <- rbind(AG1,GA1,AC1,GC1,CA1,CG1,ALa1,GLa1,LaA1,LaG1,ALe1,GLe1,LeA1,LeG1)
AG_CLaLe1Seed$Cross <- factor(AG_CLaLe1Seed$Cross, 
                             levels = c("AuxG", "GxAu",
                                        "AuxLa", "LaxAu", "GxLa", "LaxG",
                                        "AuxLe", "LexAu", "GxLe","LexG",
                                        "AuxC", "CxAu", "GxC","CxG"))

# plot data
stripchart(mean ~ Cross, data=AG_CLaLe1Seed, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AG_CLaLe1Seed_simple <- lm(mean ~ Cross, data = AG_CLaLe1Seed)
AG_CLaLe1Seed_simple <- lm(mean ~ Cross, data = AG_CLaLe1Seed)
AG_CLaLe1Seed_S_model1 <- lmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = AG_CLaLe1Seed)
AG_CLaLe1Seed_S_model2 <- lmer(mean ~ Cross + (1 | Female_parent), data = AG_CLaLe1Seed)
AG_CLaLe1Seed_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AG_CLaLe1Seed)
AG_CLaLe1Seed_S_model4 <- blmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = AG_CLaLe1Seed)
AG_CLaLe1Seed_S_model5 <- blmer(mean ~ Cross + (1 | Female_parent), data = AG_CLaLe1Seed)
AG_CLaLe1Seed_S_model6 <- blmer(mean ~ (1 | Female_parent), data = AG_CLaLe1Seed)
anova(AG_CLaLe1Seed_S_model3,AG_CLaLe1Seed_S_model4, AG_CLaLe1Seed_S_model5, AG_CLaLe1Seed_S_model6, AG_CLaLe1Seed_simple)
Anova(AG_CLaLe1Seed_simple)
DHARMa_test_model_fit(AG_CLaLe1Seed_simple)
emmeans(AG_CLaLe1Seed_simple, pairwise ~ Cross)
emmeans(AG_CLaLe1Seed_simple, pairwise ~ Cross, adjust="fdr")
plot(emmeans(AG_CLaLe1Seed_simple, pairwise ~ Cross))

#contrast for each reciprocal cross, Table 3A
contrast(emmeans(AG_CLaLe1Seed_simple, "Cross", type="response"), list(ALa=c(0,0,1,-1,0,0,0,0,0,0,0,0,0,0),
                                                                       GLa=c(0,0,0,0,1,-1,0,0,0,0,0,0,0,0),
                                                                       ALe=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0),
                                                                       GLe=c(0,0,0,0,0,0,0,0,1,-1,0,0,0,0),
                                                                       AC=c(0,0,0,0,0,0,0,0,0,0,1,-1,0,0),
                                                                       GC=c(0,0,0,0,0,0,0,0,0,0,0,0,1,-1)))
contrast(emmeans(AG_CLaLe1Seed_simple, "Cross", type="response"), list(ALa=c(0,0,1,-1,0,0,0,0,0,0,0,0,0,0),
                                                                       GLa=c(0,0,0,0,1,-1,0,0,0,0,0,0,0,0),
                                                                       ALe=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0),
                                                                       GLe=c(0,0,0,0,0,0,0,0,1,-1,0,0,0,0),
                                                                       AC=c(0,0,0,0,0,0,0,0,0,0,1,-1,0,0),
                                                                       GC=c(0,0,0,0,0,0,0,0,0,0,0,0,1,-1)), adjust="fdr")

#contrasts between A/G X other species, Table 3B
contrast(emmeans(AG_CLaLe1Seed_simple, "Cross", type="response"), list(AGvALa=c(1,1,-2,0,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvLaA=c(1,1,0,-2,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvGLa=c(1,1,0,0,-2,0,0,0,0,0,0,0,0,0),
                                                                       AGvLaG=c(1,1,0,0,0,-2,0,0,0,0,0,0,0,0),
                                                                       AGvALe=c(1,1,0,0,0,0,-2,0,0,0,0,0,0,0),
                                                                       AGvLeA=c(1,1,0,0,0,0,0,-2,0,0,0,0,0,0),
                                                                       AGvGLe=c(1,1,0,0,0,0,0,0,-2,0,0,0,0,0),
                                                                       AGvLeG=c(1,1,0,0,0,0,0,0,0,-2,0,0,0,0),
                                                                       AGvAC=c(1,1,0,0,0,0,0,0,0,0,-2,0,0,0),
                                                                       AGvCA=c(1,1,0,0,0,0,0,0,0,0,0,-2,0,0),
                                                                       AGvGC=c(1,1,0,0,0,0,0,0,0,0,0,0,-2,0),
                                                                       AGvCG=c(1,1,0,0,0,0,0,0,0,0,0,0,0,-2)))
contrast(emmeans(AG_CLaLe1Seed_simple, "Cross", type="response"), list(AGvALa=c(1,1,-2,0,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvLaA=c(1,1,0,-2,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvGLa=c(1,1,0,0,-2,0,0,0,0,0,0,0,0,0),
                                                                       AGvLaG=c(1,1,0,0,0,-2,0,0,0,0,0,0,0,0),
                                                                       AGvALe=c(1,1,0,0,0,0,-2,0,0,0,0,0,0,0),
                                                                       AGvLeA=c(1,1,0,0,0,0,0,-2,0,0,0,0,0,0),
                                                                       AGvGLe=c(1,1,0,0,0,0,0,0,-2,0,0,0,0,0),
                                                                       AGvLeG=c(1,1,0,0,0,0,0,0,0,-2,0,0,0,0),
                                                                       AGvAC=c(1,1,0,0,0,0,0,0,0,0,-2,0,0,0),
                                                                       AGvCA=c(1,1,0,0,0,0,0,0,0,0,0,-2,0,0),
                                                                       AGvGC=c(1,1,0,0,0,0,0,0,0,0,0,0,-2,0),
                                                                       AGvCG=c(1,1,0,0,0,0,0,0,0,0,0,0,0,-2)), adjust="fdr")

#contrast with each pair, Table 3C
contrast(emmeans(AG_CLaLe1Seed_simple, "Cross", type="response"), list(AGvALa=c(1,1,-1,-1,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvGLa=c(1,1,0,0,-1,-1,0,0,0,0,0,0,0,0),
                                                                       AGvALe=c(1,1,0,0,0,0,-1,-1,0,0,0,0,0,0),
                                                                       AGvGLe=c(1,1,0,0,0,0,0,0,-1,-1,0,0,0,0),
                                                                       AGvAC=c(1,1,0,0,0,0,0,0,0,0,-1,-1,0,0),
                                                                       AGvGC=c(1,1,0,0,0,0,0,0,0,0,0,0,-1,-1)))
contrast(emmeans(AG_CLaLe1Seed_simple, "Cross", type="response"), list(AGvALa=c(1,1,-1,-1,0,0,0,0,0,0,0,0,0,0),
                                                                    AGvGLa=c(1,1,0,0,-1,-1,0,0,0,0,0,0,0,0),
                                                                    AGvALe=c(1,1,0,0,0,0,-1,-1,0,0,0,0,0,0),
                                                                    AGvGLe=c(1,1,0,0,0,0,0,0,-1,-1,0,0,0,0),
                                                                    AGvAC=c(1,1,0,0,0,0,0,0,0,0,-1,-1,0,0),
                                                                    AGvGC=c(1,1,0,0,0,0,0,0,0,0,0,0,-1,-1)), adjust="fdr")

# Figure 3A
estimatesAGCLaLe1 <- emmeans(AG_CLaLe1Seed_simple, pairwise ~Cross, type="response")[[1]] %>% data.frame()

AGCLaLe1_comp <- ggplot(AG_CLaLe1Seed, aes(x=Cross, y=mean)) + ylim(0, 4.1) +
  geom_jitter(aes(fill = Cross, color = Cross), size = 2, alpha = 0.7,
              position = position_jitter(width = 0.2)) +
  geom_pointrange(aes(x = Cross, y = emmean, ymin = lower.CL, ymax = upper.CL, color=Cross), 
                  data = estimatesAGCLaLe1,
                  position = position_nudge(x=0.25)) +
  labs(x = "Cross", y = "Mean number of seeds/fruit") +
  theme_bw()
AGCLaLe1_comp

#~~~~~~~~~~~~~~~~~~~~~ all crosses (including 0 seeds) ~~~~~~~~~~~~~~~~~~~~~# 
AG_CLaLeSeed <- rbind(AG,GA,AC,GC,CA,CG,ALa,GLa,LaA,LaG,ALe,GLe,LeA,LeG)
AG_CLaLeSeed$Cross <- factor(AG_CLaLeSeed$Cross, 
                              levels = c("AuxG", "GxAu",
                                         "AuxLa", "LaxAu", "GxLa", "LaxG",
                                         "AuxLe", "LexAu", "GxLe","LexG",
                                         "AuxC", "CxAu", "GxC","CxG"))

# plot data
stripchart(mean ~ Cross, data=AG_CLaLeSeed, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AG_CLaLeSeed_S_model4 <- glmmTMB(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), family=tweedie, ziformula=~Cross, data = AG_CLaLeSeed)
AG_CLaLeSeed_simple <- lm(mean ~ Cross, data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model1 <- lmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model2 <- lmer(mean ~ Cross + (1 | Female_parent), data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model4 <- glmmTMB(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), family=tweedie, ziformula=~Cross, data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model5 <- glmmTMB(mean ~ Cross + (1 | Female_parent), family=tweedie, ziformula=~Cross,data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model6 <- glmmTMB(mean ~ (1 | Female_parent), family=tweedie, ziformula=~Cross, data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model7 <- glmmTMB(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), ziformula=~Cross, data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model8 <- glmmTMB(mean ~ Cross + (1 | Female_parent), ziformula=~Cross,data = AG_CLaLeSeed)
AG_CLaLeSeed_S_model9 <- glmmTMB(mean ~ (1 | Female_parent), ziformula=~Cross, data = AG_CLaLeSeed)
anova(AG_CLaLeSeed_S_model2,AG_CLaLeSeed_S_model3, AG_CLaLeSeed_simple)
anova(AG_CLaLeSeed_S_model4, AG_CLaLeSeed_S_model9)
Anova(AG_CLaLeSeed_simple)
Anova(AG_CLaLeSeed_S_model4)
DHARMa_test_model_fit(AG_CLaLeSeed_simple) #zero-inflated
DHARMa_test_model_fit(AG_CLaLeSeed_S_model4)
emmeans(AG_CLaLeSeed_S_model4, pairwise ~ Cross)
emmeans(AG_CLaLeSeed_S_model4, pairwise ~ Cross, adjust="fdr")
plot(emmeans(AG_CLaLeSeed_S_model4, pairwise ~ Cross))

#contrast for each reciprocal cross, Table 3A
contrast(emmeans(AG_CLaLeSeed_S_model4, "Cross", type="response"), list(ALa=c(0,0,1,-1,0,0,0,0,0,0,0,0,0,0),
                                                                        GLa=c(0,0,0,0,1,-1,0,0,0,0,0,0,0,0),
                                                                        ALe=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0),
                                                                        GLe=c(0,0,0,0,0,0,0,0,1,-1,0,0,0,0),
                                                                        AC=c(0,0,0,0,0,0,0,0,0,0,1,-1,0,0),
                                                                        GC=c(0,0,0,0,0,0,0,0,0,0,0,0,1,-1)))
contrast(emmeans(AG_CLaLeSeed_S_model4, "Cross", type="response"), list(ALa=c(0,0,1,-1,0,0,0,0,0,0,0,0,0,0),
                                                                        GLa=c(0,0,0,0,1,-1,0,0,0,0,0,0,0,0),
                                                                        ALe=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0),
                                                                        GLe=c(0,0,0,0,0,0,0,0,1,-1,0,0,0,0),
                                                                        AC=c(0,0,0,0,0,0,0,0,0,0,1,-1,0,0),
                                                                        GC=c(0,0,0,0,0,0,0,0,0,0,0,0,1,-1)), adjust="fdr")

#contrasts between A/G X other species, Table 3B
contrast(emmeans(AG_CLaLeSeed_S_model4, "Cross", type="response"), list(AGvALa=c(1,1,-2,0,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvLaA=c(1,1,0,-2,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvGLa=c(1,1,0,0,-2,0,0,0,0,0,0,0,0,0),
                                                                       AGvLaG=c(1,1,0,0,0,-2,0,0,0,0,0,0,0,0),
                                                                       AGvALe=c(1,1,0,0,0,0,-2,0,0,0,0,0,0,0),
                                                                       AGvLeA=c(1,1,0,0,0,0,0,-2,0,0,0,0,0,0),
                                                                       AGvGLe=c(1,1,0,0,0,0,0,0,-2,0,0,0,0,0),
                                                                       AGvLeG=c(1,1,0,0,0,0,0,0,0,-2,0,0,0,0),
                                                                       AGvAC=c(1,1,0,0,0,0,0,0,0,0,-2,0,0,0),
                                                                       AGvCA=c(1,1,0,0,0,0,0,0,0,0,0,-2,0,0),
                                                                       AGvGC=c(1,1,0,0,0,0,0,0,0,0,0,0,-2,0),
                                                                       AGvCG=c(1,1,0,0,0,0,0,0,0,0,0,0,0,-2)))
contrast(emmeans(AG_CLaLeSeed_S_model4, "Cross", type="response"), list(AGvALa=c(1,1,-2,0,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvLaA=c(1,1,0,-2,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvGLa=c(1,1,0,0,-2,0,0,0,0,0,0,0,0,0),
                                                                       AGvLaG=c(1,1,0,0,0,-2,0,0,0,0,0,0,0,0),
                                                                       AGvALe=c(1,1,0,0,0,0,-2,0,0,0,0,0,0,0),
                                                                       AGvLeA=c(1,1,0,0,0,0,0,-2,0,0,0,0,0,0),
                                                                       AGvGLe=c(1,1,0,0,0,0,0,0,-2,0,0,0,0,0),
                                                                       AGvLeG=c(1,1,0,0,0,0,0,0,0,-2,0,0,0,0),
                                                                       AGvAC=c(1,1,0,0,0,0,0,0,0,0,-2,0,0,0),
                                                                       AGvCA=c(1,1,0,0,0,0,0,0,0,0,0,-2,0,0),
                                                                       AGvGC=c(1,1,0,0,0,0,0,0,0,0,0,0,-2,0),
                                                                       AGvCG=c(1,1,0,0,0,0,0,0,0,0,0,0,0,-2)), adjust="fdr")


#contrast with each pair, Table 3C
contrast(emmeans(AG_CLaLeSeed_S_model4, "Cross", type="response"), list(AGvALa=c(1,1,-1,-1,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvGLa=c(1,1,0,0,-1,-1,0,0,0,0,0,0,0,0),
                                                                       AGvALe=c(1,1,0,0,0,0,-1,-1,0,0,0,0,0,0),
                                                                       AGvGLe=c(1,1,0,0,0,0,0,0,-1,-1,0,0,0,0),
                                                                       AGvAC=c(1,1,0,0,0,0,0,0,0,0,-1,-1,0,0),
                                                                       AGvGC=c(1,1,0,0,0,0,0,0,0,0,0,0,-1,-1)))
contrast(emmeans(AG_CLaLeSeed_S_model4, "Cross", type="response"), list(AGvALa=c(1,1,-1,-1,0,0,0,0,0,0,0,0,0,0),
                                                                       AGvGLa=c(1,1,0,0,-1,-1,0,0,0,0,0,0,0,0),
                                                                       AGvALe=c(1,1,0,0,0,0,-1,-1,0,0,0,0,0,0),
                                                                       AGvGLe=c(1,1,0,0,0,0,0,0,-1,-1,0,0,0,0),
                                                                       AGvAC=c(1,1,0,0,0,0,0,0,0,0,-1,-1,0,0),
                                                                       AGvGC=c(1,1,0,0,0,0,0,0,0,0,0,0,-1,-1)), adjust="fdr")

# Figure 3B
estimatesAGCLaLe <- emmeans(AG_CLaLeSeed_S_model4, pairwise ~Cross, type="response")[[1]] %>% data.frame()

AGCLaLe_comp <- ggplot(AG_CLaLeSeed, aes(x=Cross, y=mean)) + ylim(0, 4) +
  geom_jitter(aes(fill = Cross, color = Cross), size = 2, alpha = 0.7,
              position = position_jitter(width = 0.2)) +
  geom_pointrange(aes(x = Cross, y = response, ymin = lower.CL, ymax = upper.CL, color=Cross), 
                  data = estimatesAGCLaLe,
                  position = position_nudge(x=0.25)) +
  labs(x = "Cross", y = "Mean number of seeds/fruit") +
  theme_bw()
AGCLaLe_comp

# combine A and B for Figure 3
combAGCLaLe <- ggarrange(AGCLaLe1_comp, AGCLaLe_comp,
                         labels = c("A","B"), common.legend = TRUE,
                         nrow = 2)
combAGCLaLe


########################################
# data tables for 2 way ANOVA analyses #
########################################
# create and add to new columns - each column a factor
# Au_G - whether cross with Au (austinii) or G (grandifolia)
# Sp_Male_Female - whether species A (C, La, Le) was the male or female parent in the cross

ALe1["Au_G"] <- ALe1$FemaleSp
ALe1["Sp_Male_Female"] <- "Male"
LeA1["Au_G"] <- LeA1$MaleSp
LeA1["Sp_Male_Female"] <- "Female"

GLe1["Au_G"] <- GLe1$FemaleSp
GLe1["Sp_Male_Female"] <- "Male"
LeG1["Au_G"] <- LeG1$MaleSp
LeG1["Sp_Male_Female"] <- "Female"

ALa1["Au_G"] <- ALa1$FemaleSp
ALa1["Sp_Male_Female"] <- "Male"
LaA1["Au_G"] <- LaA1$MaleSp
LaA1["Sp_Male_Female"] <- "Female"

GLa1["Au_G"] <- GLa1$FemaleSp
GLa1["Sp_Male_Female"] <- "Male"
LaG1["Au_G"] <- LaG1$MaleSp
LaG1["Sp_Male_Female"] <- "Female"

AC1["Au_G"] <- AC1$FemaleSp
AC1["Sp_Male_Female"] <- "Male"
CA1["Au_G"] <- CA1$MaleSp
CA1["Sp_Male_Female"] <- "Female"

GC1["Au_G"] <- GC1$FemaleSp
GC1["Sp_Male_Female"] <- "Male"
CG1["Au_G"] <- CG1$MaleSp
CG1["Sp_Male_Female"] <- "Female"

#not useful, but needed for good merging of files
AA1["Au_G"] <- AA1$FemaleSp
AA1["Sp_Male_Female"] <- "Male"
GG1["Au_G"] <- GG1$MaleSp
GG1["Sp_Male_Female"] <- "Female"

AG1["Au_G"] <- AG1$FemaleSp
AG1["Sp_Male_Female"] <- "Male"
GA1["Au_G"] <- GA1$MaleSp
GA1["Sp_Male_Female"] <- "Female"

###############################################
# 2-way ANOVA: crosses with Lacunosa          #
# Table 4, Figure 4A                          #
###############################################

#~~~~~~~~~~ crosses with at least 1 seed (0 seeds removed) ~~~~~~~~~~#
#subsetted data
AGLa1 <- rbind(ALa1,LaA1,GLa1,LaG1)

AGLa1$Cross <- factor(AGLa1$Cross, levels = c("LaxAu", "LaxG", "AuxLa", "GxLa"))

# plot data
stripchart(mean ~ Cross, data=AGLa1, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AGLa1_simple <- lm(mean ~ Au_G + Sp_Male_Female, data = AGLa1)
AGLa1_simple <- lm(mean ~ Au_G + Sp_Male_Female, data = AGLa1)
AGLa1_S_model1 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGLa1)
AGLa1_S_model2 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGLa1)
AGLa1_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AGLa1)
AGLa1_S_model4 <- blmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGLa1)
AGLa1_S_model5 <- blmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGLa1)
AGLa1_S_model6 <- blmer(mean ~ (1 | Female_parent), data = AGLa1)
anova(AGLa1_S_model3, AGLa1_S_model4, AGLa1_S_model6, AGLa1_simple)
Anova(AGLa1_simple)
DHARMa_test_model_fit(AGLa1_simple)
emmeans(AGLa1_simple, pairwise ~Au_G|Sp_Male_Female)
emmeans(AGLa1_simple, pairwise ~Au_G|Sp_Male_Female, adjust="none")
plot(emmeans(AGLa1_simple, pairwise ~Au_G|Sp_Male_Female))

estimates2La <- emmeans(AGLa1_simple, pairwise ~Au_G|Sp_Male_Female, type="response")[[1]] %>% data.frame()

La_comp <- ggplot(AGLa1, aes(x=Sp_Male_Female, y=mean)) + ylim(0, 4) +
  geom_jitter(aes(fill = Au_G, color = Au_G), size = 2, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  scale_fill_manual(values = c("#A8DDB5", "#0868AC")) +
  scale_color_manual(values = c("#A8DDB5", "#0868AC")) +
  geom_pointrange(aes(x = Sp_Male_Female, y = emmean, ymin = lower.CL, ymax = upper.CL, color=Au_G), 
                  data = estimates2La,
                  position = position_dodge(width=.5)) +
  labs(x = "I. lacunosa", y = "Mean number of seeds/fruit") +
  theme_bw()

La_comp

#~~~~~~~~~~~~~~~~~~~~~~~~~~ all crosses ~~~~~~~~~~~~~~~~~~~~~~~~~~#

# assign categories
ALa["Au_G"] <- ALa$FemaleSp
ALa["Sp_Male_Female"] <- "Male"
LaA["Au_G"] <- LaA$MaleSp
LaA["Sp_Male_Female"] <- "Female"

GLa["Au_G"] <- GLa$FemaleSp
GLa["Sp_Male_Female"] <- "Male"
LaG["Au_G"] <- LaG$MaleSp
LaG["Sp_Male_Female"] <- "Female"

#subsetted data
AGLa <- rbind(ALa,LaA,GLa,LaG)

AGLa$Cross <- factor(AGLa$Cross, levels = c("LaxAu", "LaxG", "AuxLa", "GxLa"))

# plot data
stripchart(mean ~ Cross, data=AGLa, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AGLa_S_model5 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), family=tweedie, ziformula=~Au_G, data = AGLa)
AGLa_simple <- lm(mean ~ Au_G + Sp_Male_Female, data = AGLa)
AGLa_S_model1 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGLa)
AGLa_S_model2 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGLa)
AGLa_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AGLa)
AGLa_S_model4 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), family=tweedie, ziformula=~Au_G, data = AGLa)
AGLa_S_model5 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), family=tweedie, ziformula=~Au_G, data = AGLa)
AGLa_S_model6 <- glmmTMB(mean ~ (1 | Female_parent), family=tweedie, ziformula=~Au_G, data = AGLa)
AGLa_S_model7 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), ziformula=~Au_G, data = AGLa)
AGLa_S_model8 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), ziformula=~Au_G, data = AGLa)
AGLa_S_model9 <- glmmTMB(mean ~ (1 | Female_parent),  ziformula=~Au_G, data = AGLa)
anova(AGLa_S_model2, AGLa_S_model3, AGLa_simple)
anova(AGLa_S_model4,AGLa_S_model5, AGLa_S_model6, AGLa_S_model8, AGLa_S_model9)
Anova(AGLa_S_model2)
Anova(AGLa_S_model8)
Anova(AGLa_S_model5)
DHARMa_test_model_fit(AGLa_S_model2) #zero-inflated
DHARMa_test_model_fit(AGLa_S_model8) #zero-inflated
DHARMa_test_model_fit(AGLa_S_model5)
emmeans(AGLa_S_model5, pairwise ~Au_G|Sp_Male_Female)
emmeans(AGLa_S_model5, pairwise ~Au_G|Sp_Male_Female, adjust="none")
plot(emmeans(AGLa_S_model5, pairwise ~Au_G|Sp_Male_Female))

estimates2La_all <- emmeans(AGLa_S_model5, pairwise ~Au_G|Sp_Male_Female, type="response")[[1]] %>% data.frame()

La_comp_all <- ggplot(AGLa, aes(x=Sp_Male_Female, y=mean)) + ylim(0, 4) +
  geom_jitter(aes(fill = Au_G, color = Au_G), size = 2, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  scale_fill_manual(values = c("#A8DDB5", "#0868AC")) +
  scale_color_manual(values = c("#A8DDB5", "#0868AC")) +
  geom_pointrange(aes(x = Sp_Male_Female, y = response, ymin = lower.CL, ymax = upper.CL, color=Au_G), 
                  data = estimates2La_all,
                  position = position_dodge(width=.5)) +
  labs(x = "I. lacunosa", y = "Mean number of seeds/fruit") +
  theme_bw()

La_comp_all

###############################################
# 2-way ANOVA: crosses with Leucantha         #
# Table 4, Figure 4B                          #
###############################################

#~~~~~~~~~~ crosses with at least 1 seed (0 seeds removed) ~~~~~~~~~~#
#subsetted data
AGLe1 <- rbind(ALe1,LeA1,GLe1,LeG1)

AGLe1$Cross <- factor(AGLe1$Cross, levels = c("LexAu", "LexG", "AuxLe", "GxLe"))

# plot data
stripchart(mean ~ Cross, data=AGLe1, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AGLe1_simple <- lm(mean ~ Au_G + Sp_Male_Female, data = AGLe1)
AGLe1_simple <- lm(mean ~ Au_G + Sp_Male_Female, data = AGLe1)
AGLe1_S_model1 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGLe1)
AGLe1_S_model2 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGLe1)
AGLe1_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AGLe1)
AGLe1_S_model4 <- blmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGLe1)
AGLe1_S_model5 <- blmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGLe1)
AGLe1_S_model6 <- blmer(mean ~ (1 | Female_parent), data = AGLe1)
anova(AGLe1_S_model4, AGLe1_S_model5, AGLe1_S_model6, AGLe1_simple)
Anova(AGLe1_simple)
DHARMa_test_model_fit(AGLe1_simple)
emmeans(AGLe1_simple, pairwise ~Au_G|Sp_Male_Female)
emmeans(AGLe1_simple, pairwise ~Au_G|Sp_Male_Female, adjust="none")
plot(emmeans(AGLe1_simple, pairwise ~Au_G|Sp_Male_Female))


estimates2Le <- emmeans(AGLe1_simple, pairwise ~Au_G|Sp_Male_Female, type="response")[[1]] %>% data.frame()
pvalue2Le <- emmeans(AGLe1_simple, pairwise ~Au_G|Sp_Male_Female, type="response")[[2]] %>% data.frame()

brewer.pal(n = 9, name = "GnBu")

Le_comp <- ggplot(AGLe1, aes(x=Sp_Male_Female, y=mean)) + ylim(0, 4) +
  geom_jitter(aes(fill = Au_G, color = Au_G), size = 2, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  scale_fill_manual(values = c("#A8DDB5", "#0868AC")) +
  scale_color_manual(values = c("#A8DDB5", "#0868AC")) +
  geom_pointrange(aes(x = Sp_Male_Female, y = emmean, ymin = lower.CL, ymax = upper.CL, color=Au_G), 
                  data = estimates2Le,
                  position = position_dodge(width=.5)) +
  labs(x = "I. leucantha", y = "Mean number of seeds/fruit") +
  theme_bw()

Le_comp

#~~~~~~~~~~~~~~~~~~~~~~~~~~ all crosses ~~~~~~~~~~~~~~~~~~~~~~~~~~#

# assign categories
ALe["Au_G"] <- ALe$FemaleSp
ALe["Sp_Male_Female"] <- "Male"
LeA["Au_G"] <- LeA$MaleSp
LeA["Sp_Male_Female"] <- "Female"

GLe["Au_G"] <- GLe$FemaleSp
GLe["Sp_Male_Female"] <- "Male"
LeG["Au_G"] <- LeG$MaleSp
LeG["Sp_Male_Female"] <- "Female"

#subsetted data
AGLe <- rbind(ALe,LeA,GLe,LeG)

AGLe$Cross <- factor(AGLe$Cross, levels = c("LexAu", "LexG", "AuxLe", "GxLe"))

# plot data
stripchart(mean ~ Cross, data=AGLe, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AGLe_S_model5 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent),  family=tweedie, ziformula=~Au_G, data = AGLe)
AGLe_simple <- lm(mean ~ Au_G + Sp_Male_Female, data = AGLe)
AGLe_S_model1 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGLe)
AGLe_S_model2 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGLe)
AGLe_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AGLe)
AGLe_S_model4 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), family=tweedie, ziformula=~Au_G, data = AGLe)
AGLe_S_model5 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent),  family=tweedie, ziformula=~Au_G, data = AGLe)
AGLe_S_model6 <- glmmTMB(mean ~ (1 | Female_parent),  family=tweedie, ziformula=~Au_G, data = AGLe)
AGLe_S_model7 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), ziformula=~Au_G, data = AGLe)
AGLe_S_model8 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), ziformula=~Au_G, data = AGLe)
AGLe_S_model9 <- glmmTMB(mean ~ (1 | Female_parent), ziformula=~Au_G, data = AGLe)
anova(AGLe_S_model1, AGLe_S_model2, AGLe_S_model3, AGLe_simple)
anova(AGLe_S_model4, AGLe_S_model5, AGLe_S_model6, AGLe_S_model8, AGLe_S_model9)
Anova(AGLe_simple)
Anova(AGLe_S_model8)
Anova(AGLe_S_model5)
DHARMa_test_model_fit(AGLe_simple) #zero-inflated
DHARMa_test_model_fit(AGLe_S_model8) #zero-inflated
DHARMa_test_model_fit(AGLe_S_model5)
emmeans(AGLe_S_model5, pairwise ~Au_G|Sp_Male_Female)
emmeans(AGLe_S_model5, pairwise ~Au_G|Sp_Male_Female, adjust="none")
plot(emmeans(AGLe_S_model5, pairwise ~Au_G|Sp_Male_Female))

estimates2Le_all <- emmeans(AGLe_S_model5, pairwise ~Au_G|Sp_Male_Female, type="response")[[1]] %>% data.frame()

Le_comp_all <- ggplot(AGLe, aes(x=Sp_Male_Female, y=mean)) + ylim(0, 4) +
  geom_jitter(aes(fill = Au_G, color = Au_G), size = 2, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  scale_fill_manual(values = c("#A8DDB5", "#0868AC")) +
  scale_color_manual(values = c("#A8DDB5", "#0868AC")) +
  geom_pointrange(aes(x = Sp_Male_Female, y = response, ymin = lower.CL, ymax = upper.CL, color=Au_G), 
                  data = estimates2Le_all,
                  position = position_dodge(width=.5)) +
  labs(x = "I. leucantha", y = "Mean number of seeds/fruit") +
  theme_bw()

Le_comp_all

###############################################
# 2-way ANOVA: crosses with Cordatotriloba    #
# Table 4, Figure 4C                          #
###############################################

#~~~~~~~~~~ crosses with at least 1 seed (0 seeds removed) ~~~~~~~~~~#
#subsetted data
AGC1 <- rbind(AC1,CA1,GC1,CG1)

AGC1$Cross <- factor(AGC1$Cross, levels = c("CxAu", "CxG", "AuxC", "GxC"))

# plot data
stripchart(mean ~ Cross, data=AGC1, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AGC1_simple <- lm(mean ~ Au_G + Sp_Male_Female, data = AGC1)
AGC1_simple<- lm(mean ~ Au_G + Sp_Male_Female, data = AGC1)
AGC1_S_model1 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGC1)
AGC1_S_model2 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGC1)
AGC1_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AGC1)
AGC1_S_model4 <- blmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGC1)
AGC1_S_model5 <- blmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGC1)
AGC1_S_model6 <- blmer(mean ~ (1 | Female_parent), data = AGC1)
anova(AGC1_S_model1, AGC1_S_model2, AGC1_S_model3, AGC1_simple)
Anova(AGC1_simple)
DHARMa_test_model_fit(AGC1_simple)
emmeans(AGC1_simple, pairwise ~Au_G|Sp_Male_Female)
emmeans(AGC1_simple, pairwise ~Au_G|Sp_Male_Female, adjust="none")
plot(emmeans(AGC1_simple, pairwise ~Au_G|Sp_Male_Female))


estimates2C <- emmeans(AGC1_simple, pairwise ~Au_G|Sp_Male_Female, type="response")[[1]] %>% data.frame()

C_comp <- ggplot(AGC1, aes(x=Sp_Male_Female, y=mean)) + ylim(0, 4) +
  geom_jitter(aes(fill = Au_G, color = Au_G), size = 2, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  scale_fill_manual(values = c("#A8DDB5", "#0868AC")) +
  scale_color_manual(values = c("#A8DDB5", "#0868AC")) +
  geom_pointrange(aes(x = Sp_Male_Female, y = emmean, ymin = lower.CL, ymax = upper.CL, color=Au_G), 
                  data = estimates2C,
                  position = position_dodge(width=.5)) +
  labs(x = "I. cordatotriloba", y = "Mean number of seeds/fruit") +
  theme_bw()

C_comp

#~~~~~~~~~~~~~~~~~~~~~~~~~~ all crosses ~~~~~~~~~~~~~~~~~~~~~~~~~~#

# assign categories
AC["Au_G"] <- AC$FemaleSp
AC["Sp_Male_Female"] <- "Male"
CA["Au_G"] <- CA$MaleSp
CA["Sp_Male_Female"] <- "Female"
GC["Au_G"] <- GC$FemaleSp
GC["Sp_Male_Female"] <- "Male"
CG["Au_G"] <- CG$MaleSp
CG["Sp_Male_Female"] <- "Female"

#subsetted data
AGC <- rbind(AC,CA,GC,CG)

AGC$Cross <- factor(AGC$Cross, levels = c("CxAu", "CxG", "AuxC", "GxC"))

# plot data
stripchart(mean ~ Cross, data=AGC, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: AGC_S_model4 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), family=tweedie, ziformula=~Au_G, data = AGC)
AGC_simple<- lm(mean ~ Au_G + Sp_Male_Female, data = AGC)
AGC_S_model1 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), data = AGC)
AGC_S_model2 <- lmer(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), data = AGC)
AGC_S_model3 <- lmer(mean ~ (1 | Female_parent), data = AGC)
AGC_S_model4 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), family=tweedie, ziformula=~Au_G, data = AGC)
AGC_S_model5 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent),family=tweedie, ziformula=~Au_G, data = AGC)
AGC_S_model6 <- glmmTMB(mean ~ (1 | Female_parent), family=tweedie, ziformula=~Au_G, data = AGC)
AGC_S_model7 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent) + (1 | Male_parent), ziformula=~Au_G, data = AGC)
AGC_S_model8 <- glmmTMB(mean ~ Au_G + Sp_Male_Female + (1 | Female_parent), ziformula=~Au_G, data = AGC)
AGC_S_model9 <- glmmTMB(mean ~ (1 | Female_parent), ziformula=~Au_G, data = AGC)
anova(AGC_S_model2, AGC_S_model3, AGC_simple)
anova(AGC_S_model4, AGC_S_model5, AGC_S_model6, AGC_S_model7, AGC_S_model8, AGC_S_model9)
Anova(AGC_S_model2)
Anova(AGC_S_model4)
DHARMa_test_model_fit(AGC_S_model2) #zero-inflated, deviation in QQ plot residuals
DHARMa_test_model_fit(AGC_S_model4)
emmeans(AGC_S_model4, pairwise ~Au_G|Sp_Male_Female)
emmeans(AGC_S_model4, pairwise ~Au_G|Sp_Male_Female, adjust="none")
plot(emmeans(AGC_S_model4, pairwise ~Au_G|Sp_Male_Female))

# Figure 4F
estimates2C_all <- emmeans(AGC_S_model4, pairwise ~Au_G|Sp_Male_Female, type="response")[[1]] %>% data.frame()

C_comp_all <- ggplot(AGC, aes(x=Sp_Male_Female, y=mean)) + ylim(0, 4) +
  geom_jitter(aes(fill = Au_G, color = Au_G), size = 2, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  scale_fill_manual(values = c("#A8DDB5", "#0868AC")) +
  scale_color_manual(values = c("#A8DDB5", "#0868AC")) +
  geom_pointrange(aes(x = Sp_Male_Female, y = response, ymin = lower.CL, ymax = upper.CL, color=Au_G), 
                  data = estimates2C_all,
                  position = position_dodge(width=.5)) +
  labs(x = "I. cordatotriloba", y = "Mean number of seeds/fruit") +
  theme_bw()

C_comp_all

############
# Figure 4 #
############
# combine 2-way ANOVA 0 seeds removed into 1 row - Figure 4 A,B,C
comb2way2 <- ggarrange(La_comp, Le_comp, C_comp,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1, common.legend = TRUE)
comb2way2

# combine 2-way ANOVA all crosses into 1 row - Figure 4 D,E,F
comb2way2b <- ggarrange(La_comp_all, Le_comp_all, C_comp_all,
                        labels = c("D", "E", "F"),
                        ncol = 3, nrow = 1, common.legend = TRUE)
comb2way2b


###################################
# for leucantha and lacunosa only #
# Figure S6, Table S7             #
###################################
#~~~~~~~~~~ crosses >= 1 seed; 0 seeds removed ~~~~~~~~~~#
LeLeLaS_1 <- rbind(LeLe1,LeLa1,LaLe1,LaS1)

LeLeLaS_1$Cross <- factor(LeLeLaS_1$Cross, levels = c("LexLe", "LexLa", "LaxLe", "LaxS"))

# plot data
stripchart(mean ~ Cross, data=LeLeLaS_1, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per fruit", 
           pch = 16, cex = 1.5, lty = 1)

# test models of mean seed number for Au and G
# model chosen: LeLeLaS_1_simple <- lm(mean ~ Cross, data = LeLeLaS_1)
LeLeLaS_1_simple <- lm(mean ~ Cross, data = LeLeLaS_1)
LeLeLaS_1_S_model1 <- lmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = LeLeLaS_1)
LeLeLaS_1_S_model2 <- lmer(mean ~ Cross + (1 | Female_parent), data = LeLeLaS_1)
LeLeLaS_1_S_model3 <- lmer(mean ~ (1 | Female_parent), data = LeLeLaS_1)
LeLeLaS_1_S_model4 <- blmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = LeLeLaS_1)
LeLeLaS_1_S_model5 <- blmer(mean ~ Cross + (1 | Female_parent), data = LeLeLaS_1)
LeLeLaS_1_S_model6 <- blmer(mean ~ (1 | Female_parent), data = LeLeLaS_1)
anova(LeLeLaS_1_S_model3, LeLeLaS_1_S_model5, LeLeLaS_1_S_model6, LeLeLaS_1_simple)
Anova(LeLeLaS_1_simple)
DHARMa_test_model_fit(LeLeLaS_1_simple)
emmeans(LeLeLaS_1_simple, pairwise ~ Cross)
emmeans(LeLeLaS_1_simple, pairwise ~ Cross, adjust="none")
plot(emmeans(LeLeLaS_1_simple, pairwise ~ Cross))
contrast(emmeans(LeLeLaS_1_simple, "Cross"), list(LevLa=c(1,0,0,-1), LavH=c(0,-1,-1,2), H1vH2=c(0,1,-1,0), LeLavH = c(1,-1,-1,1)))
contrast(emmeans(LeLeLaS_1_simple, "Cross"), list(LevLa=c(1,0,0,-1), LavH=c(0,-1,-1,2), H1vH2=c(0,1,-1,0), LeLavH = c(1,-1,-1,1)), adjust = "fdr")

estimatesLeLa1 <- emmeans(LeLeLaS_1_simple, pairwise ~Cross, type="response")[[1]] %>% data.frame()

# Figure S6C
LeLeLaS1_comp <- ggplot(LeLeLaS_1, aes(x=Cross, y=mean)) + ylim(-0.1, 4.1) +
  geom_jitter(aes(fill = Cross, color = Cross), size = 2, alpha = 0.7,
              position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("#A8DDB5","#7BCCC4","#43A2CA","#0868AC")) +
  scale_color_manual(values = c("#A8DDB5","#7BCCC4","#43A2CA","#0868AC")) +
  geom_pointrange(aes(x = Cross, y = emmean, ymin = lower.CL, ymax = upper.CL, color=Cross), 
                  data = estimatesLeLa1,
                  position = position_nudge(x=0.25)) +
  labs(x = "Cross", y = "Mean number of seeds/fruit") +
  theme_bw()
LeLeLaS1_comp

#~~~~~~~~~~~~ all crosses ~~~~~~~~~~~~#
#subsetted data
LeLaSeed <- rbind(LeLe,LeLa,LaLe,LaS)
LeLaCross <- rbind(LeLec,LeLac,LaLec,LaSc)

LeLaSeed$Cross <- factor(LeLaSeed$Cross, levels = c("LexLe", "LexLa", "LaxLe", "LaxS"))
LeLaCross$Cross <- factor(LeLaCross$Cross, levels = c("LexLe", "LexLa", "LaxLe", "LaxS"))

# plot data
stripchart(mean ~ Cross, data=LeLaSeed, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per pod", 
           pch = 16, cex = 1.5, lty = 1)
# Figure S6B
mosaicplot(table(LeLaCross$Cross, LeLaCross$Number_of_seeds), las=2, main=NULL, ylab = "Number of seeds per fruit", xlab = "Cross type")

# test models of mean seed number for Le and La
# model chosen: LeLa_S_simple <- lm(mean ~ Cross, data = LeLaSeed)
LeLa_S_simple <- lm(mean ~ Cross, data = LeLaSeed)
LeLa_S_model1 <- glmmTMB(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), family = tweedie, ziformula=~Cross, data = LeLaSeed)
LeLa_S_model2 <- glmmTMB(mean ~ Cross + (1 | Female_parent), ziformula=~Cross, family = tweedie, data = LeLaSeed)
LeLa_S_model3 <- glmmTMB(mean ~ (1 | Female_parent), ziformula=~Cross, family = tweedie, data = LeLaSeed)
LeLa_S_model4 <- glmmTMB(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), ziformula=~Cross, data = LeLaSeed)
LeLa_S_model5 <- glmmTMB(mean ~ Cross + (1 | Female_parent), ziformula=~Cross, data = LeLaSeed)
LeLa_S_model6 <- glmmTMB(mean ~ (1 | Female_parent), ziformula=~Cross, data = LeLaSeed)
LeLa_S_model7 <- lmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = LeLaSeed)
LeLa_S_model8 <- lmer(mean ~ Cross + (1 | Female_parent), data = LeLaSeed)
LeLa_S_model9 <- lmer(mean ~ (1 | Female_parent), data = LeLaSeed)
LeLa_S_model10 <- blmer(mean ~ Cross + (1 | Female_parent) + (1 | Male_parent), data = LeLaSeed)
LeLa_S_model11 <- blmer(mean ~ Cross + (1 | Female_parent), data = LeLaSeed)
LeLa_S_model12 <- blmer(mean ~ (1 | Female_parent), data = LeLaSeed)
anova(LeLa_S_model9, LeLa_S_model11, LeLa_S_model12, LeLa_S_simple)
anova(LeLa_S_model3, LeLa_S_model4, LeLa_S_model5)
Anova(LeLa_S_simple)
Anova(LeLa_S_model5)
DHARMa_test_model_fit(LeLa_S_simple) #zero-inflated
DHARMa_test_model_fit(LeLa_S_model5) #zero-inflated
emmeans(LeLa_S_simple, pairwise ~ Cross)
plot(emmeans(LeLa_S_simple, pairwise ~ Cross))
contrast(emmeans(LeLa_S_simple, "Cross"), list(LevLa=c(1,0,0,-1), LavH=c(0,-1,-1,2), H1vH2=c(0,1,-1,0), LeLavH = c(1,-1,-1,1)))
contrast(emmeans(LeLa_S_simple, "Cross"), list(LevLa=c(1,0,0,-1), LavH=c(0,-1,-1,2), H1vH2=c(0,1,-1,0), LeLavH = c(1,-1,-1,1)), adjust = "fdr")

estimatesLeLa_all <- emmeans(LeLa_S_simple, pairwise ~Cross, type="response")[[1]] %>% data.frame()

#Figure S6D
LeLeLaS_comp <- ggplot(LeLaSeed, aes(x=Cross, y=mean)) + ylim(-0.1, 4.1) +
  geom_jitter(aes(fill = Cross, color = Cross), size = 2, alpha = 0.7,
              position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("#A8DDB5","#7BCCC4","#43A2CA","#0868AC")) +
  scale_color_manual(values = c("#A8DDB5","#7BCCC4","#43A2CA","#0868AC")) +
  geom_pointrange(aes(x = Cross, y = emmean, ymin = lower.CL, ymax = upper.CL, color=Cross), 
                  data = estimatesLeLa_all,
                  position = position_nudge(x=0.25)) +
  labs(x = "Cross", y = "Mean number of seeds/fruit") +
  theme_bw()
LeLeLaS_comp

#Combine Figure S6C and S6D (delete top row later)
combS6 <- ggarrange(LeLeLaS1_comp, LeLeLaS_comp,LeLeLaS1_comp, LeLeLaS_comp,
                       labels = c("A", "B", "C", "D"),
                       ncol = 2, nrow = 2, common.legend = TRUE)
combS6

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
pie(c(t[29,2],t[29,3]), label = c("fail","success"), main="LexLe")
pie(c(t[28,2],t[28,3]), label = c("fail","success"), main="LexLa")
pie(c(t[22,2],t[22,3]), label = c("fail","success"), main="LaxLe")
pie(c(t[23,2],t[23,3]), label = c("fail","success"), main="LaxS")

pdf("FigureS6_210827.pdf", useDingbats=FALSE)
combS6
mosaicplot(table(LeLaCross$Cross, LeLaCross$Number_of_seeds), las=2, main=NULL, ylab = "Number of seeds per fruit", xlab = "Cross type")
par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
pie(c(t[29,2],t[29,3]), label = c("fail","success"), main="LexLe")
pie(c(t[28,2],t[28,3]), label = c("fail","success"), main="LexLa")
pie(c(t[22,2],t[22,3]), label = c("fail","success"), main="LaxLe")
pie(c(t[23,2],t[23,3]), label = c("fail","success"), main="LaxS")
dev.off()

################################
# G tests for fruit set        #
# Stats for Tables 1, 3, 4, S6 #
################################
library(DescTools)
library(RVAideMemoire)

seedSF <- read.table("SeedPercent_Gtest.txt", sep="\t", header=TRUE, row.names = 1)

seedSF <- as.matrix(seedSF)
G.test(seedSF)
allNoCorr <- pairwise.G.test(seedSF,p.method = "none") 
allCorr <- pairwise.G.test(seedSF,p.method = "fdr") 

#write.table(allNoCorr$p.value, "SeedPairwiseGtest_noCorr_210228.txt")
#write.table(allCorr$p.value, "SeedPairwiseGtest_CorrFDF_210228.txt")

sink("Gtest_attempt_210228.txt")
# A: Table 1
AAGG_g <- seedSF[row.names(seedSF) %in% c("AuxAu","GxG"),]
G.test(AAGG_g)
AG_g <- seedSF[row.names(seedSF) %in% c("AuxG","GxAu"),]
G.test(AG_g)
AAGG_AG_g <- seedSF[row.names(seedSF) %in% c("AAGG","(AuxG)"),]
G.test(AAGG_AG_g)

# B: Table 3A, reciprocal crosses
AuLa_g <- seedSF[row.names(seedSF) %in% c("AuxLa","LaxAu"),]
G.test(AuLa_g)
GLa_g <- seedSF[row.names(seedSF) %in% c("GxLa","LaxG"),]
G.test(GLa_g)
AuLe_g <- seedSF[row.names(seedSF) %in% c("AuxLe","LexAu"),]
G.test(AuLe_g)
GLe_g <- seedSF[row.names(seedSF) %in% c("GxLe","LexG"),]
G.test(GLe_g)
AuC_g <- seedSF[row.names(seedSF) %in% c("AuxC","CxAu"),]
G.test(AuC_g)
GC_g <- seedSF[row.names(seedSF) %in% c("GxC","CxG"),]
G.test(GC_g)

# B: Table 3B, (AuxG) vs. individual cross
AuG_AuLa_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","AuxLa"),]
G.test(AuG_AuLa_g)
AuG_LaAu_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","LaxAu"),]
G.test(AuG_LaAu_g)
AuG_GLa_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","GxLa"),]
G.test(AuG_GLa_g)
AuG_LaG_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","LaxG"),]
G.test(AuG_LaG_g)

AuG_AuLe_g <- seedSF[row.names(seedSF) %in% c("(AuxG)", "AuxLe"),]
G.test(AuG_AuLe_g)
AuG_LeAu_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","LexAu"),]
G.test(AuG_LeAu_g)
AuG_GLe_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","GxLe"),]
G.test(AuG_GLe_g)
AuG_LeG_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","LexG"),]
G.test(AuG_LeG_g)

AuG_AuC_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","AuxC"),]
G.test(AuG_AuC_g)
AuG_CAu_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","CxAu"),]
G.test(AuG_CAu_g)
AuG_GC_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","GxC"),]
G.test(AuG_GC_g)
AuG_GC_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","CxG"),]
G.test(AuG_GC_g)

# B: Table 3C: (AuxG) vs. reciprocal crosses pooled
AuG_AuLa2_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","(AuxLa)"),]
G.test(AuG_AuLa2_g)
AuG_GLa2_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","(GxLa)"),]
G.test(AuG_GLa2_g)
AuG_AuLe2_g <- seedSF[row.names(seedSF) %in% c("(AuxG)", "(AuxLe)"),]
G.test(AuG_AuLe2_g)
AuG_GLe2_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","(GxLe)"),]
G.test(AuG_GLe2_g)
AuG_AuC2_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","(AuxC)"),]
G.test(AuG_AuC2_g)
AuG_GC2_g <- seedSF[row.names(seedSF) %in% c("(AuxG)","(GxC)"),]
G.test(AuG_GC2_g)

# C: Table 4
AuG_La1_g <- seedSF[row.names(seedSF) %in% c("(AuxLa)","(GxLa)"),]
G.test(AuG_La1_g)
AuG_La2_g <- seedSF[row.names(seedSF) %in% c("(AuGxLa)","(LaxAuG)"),]
G.test(AuG_La2_g)
AuG_Le1_g <- seedSF[row.names(seedSF) %in% c("(AuxLe)","(GxLe)"),]
G.test(AuG_Le1_g)
AuG_Le2_g <- seedSF[row.names(seedSF) %in% c("(AuGxLe)","(LexAuG)"),]
G.test(AuG_Le2_g)
AuG_C1_g <- seedSF[row.names(seedSF) %in% c("(AuxC)","(GxC)"),]
G.test(AuG_C1_g)
AuG_C2_g <- seedSF[row.names(seedSF) %in% c("(AuGxC)","(CxAuG)"),]
G.test(AuG_C2_g)
sink()

# D: Table S6
LeLeLaS_g <- seedSF[row.names(seedSF) %in% c("LexLe","LaxS"),]
G.test(LeLeLaS_g)
LeLa_g <- seedSF[row.names(seedSF) %in% c("LexLa","LaxLe"),]
G.test(LeLa_g)
LaS_LeLa_g <- seedSF[row.names(seedSF) %in% c("LaxS","(LexLa)"),]
G.test(LaS_LeLa_g)
LeLeLaS_LeLa_g <- seedSF[row.names(seedSF) %in% c("LeLeLaS","(LexLa)"),]
G.test(LeLeLaS_LeLa_g)

#correct for pvalues for each section
# A: Table 1
Apvalue <- c(0.4975,0.6987,0.0002664)
p.adjust(Apvalue, method = "fdr")
# B: Table 3
# Table 3A
Brecip <- c(0.2154,0.1263,0.851,0.3307,0.003585,0.0893)
p.adjust(Brecip, method="fdr")
# Table 3B
Bind <- c(0.04017,0.000486,0.002738,1.29E-06,0.0001402,0.0002741,6.38E-05,
          3.46E-07,1.24E-09,2.20E-16,1.52E-13,2.20E-16)
p.adjust(Bind, method="fdr")
# Table 3C
Bpool <- c(0.0006999,1.69E-06,6.04E-06,4.12E-08,2.20E-16,2.20E-16)
p.adjust(Bpool, method="fdr")
# C: Table 4G
Cpvalue <- c(0.1877,0.04685,0.3899,0.5715,0.2781,9.58E-04)
p.adjust(Cpvalue, method = "fdr")
C0seed <- c(0.436, 1.51E-04, 0.277, 3.82E-05, 0.209, 0.046)
p.adjust(C0seed, method = "fdr")
C_all_p <- c(0.167, 0.002, 0.139, 0.023, 0.128, 0.271)
p.adjust(C_all_p, method = "fdr")
# D: Table S6
Dpvalue <- c(3.54E-04,0.148,0.8395,0.0155)
p.adjust(Dpvalue, method = "fdr")


