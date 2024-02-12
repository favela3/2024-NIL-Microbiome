### Alonso Favela ====
### Date: May 17th 2019
###qPCR rough Analysis


#This is the location of these genes
setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL Maize Microbiome/Data/qPCR")
setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL MM2018/Data/qPCR")

#Merging Files
library(dplyr)
fluidigm <- read.csv("All.Samples.NIL.qPCR.csv", sep=",", header=TRUE)
samplenames <- read.csv("qPCR.Concentration.csv", sep=",", header=TRUE)

merged <- left_join(fluidigm, samplenames[c('Sample.ID', 'ng.ul')], by='Sample.ID')

merged$copies <- (merged$Calibrated.rConc)/(merged$DNAconc)

write.csv(merged, "RawNILqPCR.csv")

#Now I'm going to combine these files and pull out the data
library(igraph)
library(vegan)
library(Hmisc)

#This is the protocal to subet the data
#Using a master file and subseting instead of making excel sheets
master = read.csv(file ="RawNILqPCR.csv", head = T)
#Subseting data This method only runs sum on passed data
#Now I am going to go through each of the genes 
#I'm only doing my data processing with those samples that passed. 
passed=subset(master,Call== "Pass")

###NifH
nifH = subset(passed, Name == "nifH")
nifH.results.m= aggregate(nifH[, 8:17], list(nifH$Sample.ID), mean)
write.csv(nifH.results.m, file = "SubsettedFiles/nifHqPCR.csv",row.names=FALSE)


#Bacterial amoA
BamoA = subset(passed, Name == "Bacterial amoA")
BamoA.results.m= aggregate(BamoA[, 8:17], list(BamoA$Sample.ID), mean)
write.csv(BamoA.results.m, file = "SubsettedFiles/BamoAqPCR.csv",row.names=FALSE)

#Arch amoA
AamoA = subset(passed, Name == "Archaeal amoA")
AamoA.results.m= aggregate(AamoA[, 8:17], list(AamoA$Sample.ID), mean)
write.csv(AamoA.results.m, file = "SubsettedFiles/AamoAqPCR.csv",row.names=FALSE)

#nirK
nirK = subset(passed, Name == "nirK")
nirK.results.m= aggregate(nirK[, 8:17], list(nirK$Sample.ID), mean)
write.csv(nirK.results.m, file = "SubsettedFiles/nirKqPCR.csv",row.names=FALSE)

#nirS
nirS = subset(passed, Name == "nirS")
nirS.results.m= aggregate(nirS[, 8:17], list(nirS$Sample.ID), mean)
write.csv(nirS.results.m, file = "SubsettedFiles/nirSqPCR.csv",row.names=FALSE)

#nrfA | This data should be taken with a grain of salt- Seems like it didnt work 
#Seems like all my nrfA samples failed
#Perhaps I will use a different functional gene to look at these samples 
nrfA = subset(master, Name == "nrfA")
nrfA.results.m= aggregate(nrfA[, 8:17], list(nrfA$Sample.ID), mean)
write.csv(nrfA.results.m, file = "SubsettedFiles/nrfAqPCR.csv",row.names=FALSE)

#nosZI
nosZI = subset(passed, Name == "nosZI")
nosZI.results.m= aggregate(nosZI[, 8:17], list(nosZI$Sample.ID), mean)
write.csv(nosZI.results.m, file = "SubsettedFiles/nosZIqPCR.csv",row.names=FALSE)


#Outside of R I divided the Calibrated.rConcentration by the measured qubit values. 
#Then I combined them all into a master qPCR file. 
#Now I will use dplyer to merge them with the meta info

library(dplyr)
Master.qPCR <- read.csv("Master.NIL.qPCR.csv", sep=",", header=TRUE)
Metainformation <- read.csv("NIL2018MappingDNA.csv", sep=",", header=TRUE)


merged <- full_join(Master.qPCR, Metainformation, by = c("Samples"="Samples"))
 
boxplot(merged$Corrected.Value ~ merged$Genes)

merged=(subset(merged,Time== "T2"))

par(cex.axis=.60)
par(mfrow=c(2,3))

nifH=subset(merged,Genes== "nifH")
boxplot(log(nifH$Corrected.Value) ~ nifH$Pedigree,las=2,xlab =NA,
        ylab="Gene Copy Number (ng/ul)",main="nifH in NIL Rhizosphere")

amoA=subset(merged,Genes== "AamoA")
boxplot(log(amoA$Corrected.Value) ~ amoA$Pedigree,las=2, xlab=NA, 
        ylab="Gene Copy Number (ng/ul)", main= "archeal amoA in NIL Rhizopshere")

BamoA=subset(merged,Genes== "BamoA")
boxplot(log(BamoA$Corrected.Value) ~ BamoA$Pedigree,las=2, xlab=NA, 
        ylab="Gene Copy Number (ng/ul)", main="bacterial amoA in NIL Rhizosphere")

nirK=subset(merged,Genes== "nirK ")
boxplot(log(nirK$Corrected.Value) ~ nirK$Pedigree,las=2, xlab=NA, 
        ylab="Gene Copy Number (ng/ul)", main="nirK in NIL Rhizosphere")

nirS=subset(merged,Genes== "nirS")
boxplot(log(nirS$Corrected.Value) ~ nirS$Pedigree,las=2, xlab=NA,
        ylab="Gene Copy Number (ng/ul)", main="nirS in NIL Rhizosphere")

#It may be interesting to compare the relationship between the different data to see if they have a solid relationshiop
#Lets make a simple liner model to look at the data.

#Actually I'm goign to look at all the data for a sec, and see how varaince 



# These models are returning different values 

model<-lm(nifH$Corrected.Value ~nifH$Pedigree)
summary(model) 

model<-lm(nifH$Corrected.Value ~nifH$Pedigree)
summary(model) 

# 
n.dataframe <-aov(merged$Corrected.Value ~(merged$Pedigree+merged$Range
                                           +merged$Row+merged$Genes))
summary(n.dataframe)
TukeyHSD(n.dataframe)

#There are differencesin the data set for sure
# Df    Sum Sq   Mean Sq F value  Pr(>F)    
# merged$Pedigree   45 1.233e+09 2.740e+07   2.429 7.3e-07 ***
#   merged$Range       8 4.442e+09 5.553e+08  49.216 < 2e-16 ***
#   merged$Row        21 4.306e+08 2.051e+07   1.817  0.0134 *  
#   merged$Genes       6 8.000e+09 1.333e+09 118.176 < 2e-16 ***
#   Residuals       1112 1.255e+10 1.128e+07                    



n.dataframe <-aov(nifH$Corrected.Value ~(nifH$Pedigree+nifH$Range+nifH$Row))
summary(n.dataframe)
TukeyHSD(n.dataframe)

# > summary(n.dataframe)
# Df    Sum Sq  Mean Sq F value Pr(>F)    
# nifH$Pedigree 45  47544835  1056552   1.494 0.0532 .  
# nifH$Range     8 124939336 15617417  22.081 <2e-16 ***
#   nifH$Row      21  17861532   850549   1.203 0.2681    
# Residuals     92  65070089   707284                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

n.dataframe <-aov(amoA$Corrected.Value ~(amoA$Pedigree+amoA$Range+amoA$Row))
summary(n.dataframe)
TukeyHSD(n.dataframe)


n.dataframe <-aov(BamoA$Corrected.Value ~(BamoA$Pedigree+BamoA$Range+BamoA$Row))
summary(n.dataframe)
TukeyHSD(n.dataframe)


n.dataframe <-aov(nirK$Corrected.Value ~(nirK$Pedigree+nirK$Range+nirK$Row))
summary(n.dataframe)
TukeyHSD(n.dataframe)


n.dataframe <-aov(nirS$Corrected.Value ~(nirS$Pedigree+nirS$Range+nirS$Row))
summary(n.dataframe)
TukeyHSD(n.dataframe)

#Most of the genes do not show a clear pattern.
# Iam interested in looking at the predictive relationship between the two differet factors
