#Date: 9/20/18
#Name: Alonso Favela ====

#Purpose: To give a rough analysis of the NIL nitfication Dataset to see if we can identfiy
# NILs that would be useful for further analysis and study 

#Ive decided to also add the DEA data to this analysis 

#This directory is where I will store all scripts related to the NIL project
setwd("~/Documents/Thesis /Chapter 4/Science /Rstats/Old R-Analysis /R-Scripts/NIL CSV/Functional Data")


library(igraph)
library(vegan)
library(ggplot2)
library(forcats)

##Master Data sheet 3.2.21
Master = read.csv(file ="MasterFunction.DataSheet.csv", head = T)
Master = read.csv(file ="MasterFunction.DataSheet.Corrected.3.14.21.csv", head = T)
Master = read.csv(file ="MasterFunction.DataSheet.Corrected.3.14.21.csv", head = T)


#Read CSV. This csv in particular has all of the nitrfication data from t1 and t2
#these are little tester data sets to get a better look at the data (Still has typos)
NILNitf = read.csv(file ="NIL18.All.Nitrfication.DATA.csv", head = T)
#SingleOutlier Removed
NILNitf = read.csv(file ="NIL18.All.Nitrfication.DATA.OUTlIE.csv", head = T)
#Without teosinte
NILNitf = read.csv(file ="NIL18.All.Nitrfication.DATA.without.maize.csv", head = T)
#Corrected some typos
# NILNitf = read.csv(file ="NIL18.All.Nitrfication.DATA.Check.csv", head = T)
# NILNitf = read.csv(file ="NIL18.All.Nitrfication.DATA.Check.2.csv", head = T)
# NILNitf = read.csv(file ="NIL18.All.Nitrfication.DATA.Check.3.csv", head = T)
NILAngelaFigure = read.csv(file ="NIL18.All.Nitrfication.DATA.Check.Angela.csv", head = T)
Selectnitf=NILNitf
#There are a couple of outliers in the data that are altering the found results. 

#Simple linar model will all data. I ran two models one with pedigree as a factor and one as time. 
#Pedigree was not sig
#Time was very sig
#Results below
#Row had no effect/ range slightly does
n.dataframe <-aov(NILNitf$Nitrifcation.ngN.gh ~(NILNitf$Pedigree))
summary(n.dataframe)


n.dataframe <-aov(NILNitf$Log ~(NILNitf$Pedigree))
summary(n.dataframe)

# Df    Sum Sq   Mean Sq F value Pr(>F)
# NILNitf$Pedigree  47 3.827e+12 8.143e+10   0.896  0.668
# Residuals        309 2.808e+13 9.088e+10   
#-------------------------------------------------------------
# Df    Sum Sq   Mean Sq F value   Pr(>F)    
# NILNitf$Time   1 1.127e+12 1.127e+12   12.99 0.000357 ***
#   Residuals    355 3.078e+13 8.671e+10
#-------------------------------------------------------------
# Df    Sum Sq   Mean Sq F value Pr(>F) 
# NILNitf$Range   8 1.232e+12 1.540e+11   1.748 0.0865 .
# Residuals     348 3.068e+13 8.815e+10                 


#Now I want to visualize it 
boxplot(NILNitf$Nitrifcation.ngN.gh ~(NILNitf$Pedigree))
boxplot(NILNitf$Nitrifcation.ngN.gh ~(NILNitf$Time))
boxplot(NILNitf$Nitrifcation.ngN.gh ~(NILNitf$Range))


#Now Lets subset in time
NitrT1=subset(NILNitf, Time== "T1")
NitrT2=subset(NILNitf, Time== "T2")
NitrT2Angela=subset(NILAngelaFigure, Time== "T2")

##Genotype selection =====
# Nitrfication: NIL 0578, NIL 0057, NIL 0071, B73, teosinte, NIL 0047, NIL 0021
# NitrT2Angela=subset(NILAngelaFigure, Time== "T2")


Selectnitf=subset(NitrT2, Pedigree== "Z031E0028"| Pedigree== "Z031E0059"| Pedigree== "B73xPI384071"|
                    Pedigree== "Z031E0047"| Pedigree== "Z031E0021"| Pedigree== "B73"
                  | Pedigree== "PI384071" )

Selectnitf=subset(Selectnitf, Sample.ID!="NIL18.T2.152")
Selectnitf=subset(Selectnitf, Sample.ID!="NIL18.T2.009")
Selectnitf=subset(Selectnitf, Sample.ID!="NIL18.T2.083")
	



# T1 overview=====
#What is key here is understnign how things vary in time.
boxplot(NitrT1$Nitrifcation.ngN.gh ~(NitrT1$Pedigree))
boxplot(NitrT1$Log ~(NitrT1$Pedigree))
ggplot(NitrT1, aes(x = Nitrifcation.ngN.gh, fill = Pedigree)) + geom_density(alpha = 0.5)
ggplot(NitrT1, aes(x = fct_reorder(Pedigree, Nitrifcation.ngN.gh, .desc=TRUE), y=Nitrifcation.ngN.gh, fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


boxplot(NitrT1$Nitrifcation.ngN.gh ~(NitrT1$Time))
boxplot(NitrT1$Nitrifcation.ngN.gh ~(NitrT1$Range))


n.dataframe <-aov(NitrT1$Nitrifcation.ngN.gh ~(NitrT1$Pedigree))
summary(n.dataframe)


boxplot(NitrT2$Nitrifcation.ngN.gh ~(NitrT2$Pedigree))
boxplot(NitrT2$Log ~(NitrT2$Pedigree))
ggplot(NitrT2, aes(x = Nitrifcation.ngN.gh, fill = Pedigree)) + geom_density(alpha = 0.5)
ggplot(NitrT2, aes(x = fct_reorder(Pedigree, Nitrifcation.ngN.gh, .fun = median, .desc=TRUE), y=Nitrifcation.ngN.gh, fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

boxplot(NitrT2$Nitrifcation.ngN.gh ~(NitrT2$Time))
boxplot(NitrT2$Nitrifcation.ngN.gh ~(NitrT2$Range))

n.dataframe <-aov(NitrT2$Nitrifcation.ngN.gh ~(NitrT2$Pedigree))
summary(n.dataframe)
TukeyHSD(n.dataframe)

n.dataframe <-aov(NitrT2$Log ~(NitrT2$Pedigree))
summary(n.dataframe)
TukeyHSD(n.dataframe)
###Code For makeing the Specifc figure for Angela======

ggplot(Selectnitf, aes(x = fct_reorder(Pedigree, Nitrifcation.ngN.gh, .fun = median, .desc=TRUE), y=Nitrifcation.ngN.gh, fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+
  theme_classic()+scale_fill_manual(values=c("blue", "limegreen", "lightblue","lightblue", "lightblue", "lightblue", "lightblue"))


#Seems like there are a couple of outliers in the data
#this data set apears to have to much intrinc varability 
##Fig3. Nitrfication Reduced Figures=====
Master = read.csv(file ="MasterFunction.DataSheet.Corrected.3.14.21.csv", head = T)

Master=subset(Master, X.SampleID!="NIL.T1.110")
Master=subset(Master, X.SampleID!="NIL.T2.009")
Master=subset(Master, X.SampleID!="NIL.T2.152")


ggplot(Master, aes(x = fct_reorder(Pedigree, Nitrification.Rate..ugN.ghr., .fun = median, .desc=TRUE), y=Nitrification.Rate..ugN.ghr., fill=Pedigree))+
  geom_boxplot()+labs(y="Potential nitrification rate (ngN/gh)", title="Potential Nitrification Rate", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")

NitrT1=subset(Master, Time== "T1")
NitrT2=subset(Master, Time== "T2")

#Timepoint1 Fig.3----
#Outlier in Z31E507
NitrT1=subset(NitrT1, X.SampleID!="NIL.T1.110")

ggplot(NitrT1, aes(x = fct_reorder(Pedigree, Nitrification.Rate..ugN.ghr., .fun = median, .desc=TRUE), y=Nitrification.Rate..ugN.ghr., fill=Pedigree))+
  geom_boxplot()+labs(y="Potential nitrification rate (ngN/gh)", title="Timepoint 1 - Potential Nitrification Rate", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")

SelectT1=subset(NitrT1, Pedigree== "Z031E0068"| Pedigree== "Z031E0507"| Pedigree== "B73xPI384071"|
                    Pedigree== "Z031E0058"| Pedigree== "Z031E0559"| Pedigree== "B73"
                  | Pedigree== "PI384071" )

ggplot(SelectT1, aes(x = fct_reorder(Pedigree, Nitrification.Rate..ugN.ghr., .fun = median, .desc=TRUE), y=Nitrification.Rate..ugN.ghr., fill=Pedigree))+
  geom_boxplot()+labs(y="Potential nitrification rate (ngN/gh)", title="T1- Select Potential Nitrification Rate", x="Genotypes")+
  theme_classic()+scale_fill_manual(values=c("blue", "limegreen", "yellow","lightblue", "lightblue", "lightblue", "lightblue"))+
  geom_jitter( width=0.20)

##Subsetting T2 Lower


#Teosinte Outlierss
NitrT2=subset(NitrT2, Sample.ID!="NIL.T2.009")
#Outlier Z31E0021
NitrT2=subset(NitrT2, Sample.ID!="NIL.T2.152")
NitrT2=subset(NitrT2, Sample.ID!="NIL.T2.152")
#B73 Outlier
NitrT2=subset(NitrT2, Sample.ID!="NIL.T2.070")


ggplot(NitrT2, aes(x = fct_reorder(Pedigree, Nitrification.Rate..ugN.ghr., .fun = median, .desc=TRUE), y=Nitrification.Rate..ugN.ghr., fill=Pedigree))+
  geom_boxplot()+labs(y="Potential nitrification rate (ngN/gh)", title="Timepoint 2 - Potential Nitrification Rate", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")

#Subset figure NILS

Selectnitf=subset(NitrT2, Pedigree== "Z031E0028"| Pedigree== "Z031E0059"| Pedigree== "B73xPI384071"|
                    Pedigree== "Z031E0047"| Pedigree== "Z031E0021"| Pedigree== "B73"
                  | Pedigree== "PI384071" )

Selectnitf$Nitrification.Rate..ugN.ghr.[8]<-0

Selectnitf$Nitrification.Rate..ugN.ghr.[8]<-.1
Selectnitf$Nitrification.Rate..ngN.ghr.[8]<-.1

NitrT2$Nitrification.Rate..ugN.ghr.[54]<-1

NitrT2$Nitrification.Rate..ngN.ghr.[54]<-.1


#Tukey Test with Selected NILs for Manuscript|| ANVOA VALUES USED in manuscript Fig 2 =====
n.dataframe <-aov(log10(Selectnitf$Nitrification.Rate..ugN.ghr.) ~(Selectnitf$Pedi2gree))
n.dataframe <-aov(log10(Selectnitf$Nitrification.Rate..ugN.ghr.) ~(Selectnitf$Pedigree))
n.dataframe <-aov((Selectnitf$Nitrification.Rate..ugN.ghr.) ~(Selectnitf$Pedigree))

summary(n.dataframe)

TukeyHSD(n.dataframe)

colnames(NitrT2)


Selectnitf=subset(NitrT2, Pedigree== "B73"
                  | Pedigree== "PI384071" )

n.dataframe <-aov(log1p(Selectnitf$Nitrification.Rate..ugN.ghr.) ~(Selectnitf$Pedigree))

summary(n.dataframe)


Selectnitf=subset(NitrT2,  Pedigree== "Z031E0021"| Pedigree== "B73")
# 0.8210   10.27 0.0327 *
n.dataframe <-aov((Selectnitf$Nitrification.Rate..ugN.ghr.) ~(Selectnitf$Pedigree))
# 14.96    11.6 0.0271 *
summary(n.dataframe)


Selectnitf=subset(NitrT2,  Pedigree== "Z031E0047"| Pedigree== "B73")

n.dataframe <-aov((Selectnitf$Nitrification.Rate..ugN.ghr.) ~(Selectnitf$Pedigree))

summary(n.dataframe)
# 23.957   28.54 0.00308 **

Selectnitf=subset(NitrT2,  Pedigree== "B73xPI384071"| Pedigree== "B73")

n.dataframe <-aov((Selectnitf$Nitrification.Rate..ugN.ghr.) ~(Selectnitf$Pedigree))

summary(n.dataframe)

# 7.090   6.198 0.0552 .
#4x4
##Manuscript Figure=====
ggplot(Selectnitf, aes(x = fct_reorder(Pedigree, Nitrification.Rate..ugN.ghr., .fun = median, .desc=TRUE), y=Nitrification.Rate..ugN.ghr., fill=Pedigree))+
  geom_boxplot()+labs(y="Potential nitrification rate (ngN/gh)", x="Genotypes")+
  theme_classic()+
  scale_fill_manual(values=c("#FE6100", "#DC267F", "#785EF0","#FFB000", "#FFB000", "#FFB000", "#FFB000"))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



ggplot(Selectnitf, aes(x = fct_reorder(Pedigree, log10(Nitrification.Rate..ugN.ghr.), .fun = median, .desc=TRUE), y=Nitrification.Rate..ugN.ghr., fill=Pedigree))+
  geom_boxplot()+labs(y="Potential nitrification rate (ngN/gh)", title="T2- Select Potential Nitrification Rate", x="Genotypes")+
  theme_classic()+
  scale_fill_manual(values=c("#FE6100", "#DC267F", "#785EF0","#FFB000", "#FFB000", "#FFB000", "#FFB000"))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



ggplot(Selectnitf, aes(x = fct_reorder(Pedigree, log1p(Nitrification.Rate..ugN.ghr.), .fun = median, .desc=TRUE), y=Nitrification.Rate..ugN.ghr., fill=Pedigree))+
  geom_boxplot()+labs(y="Potential nitrification rate (ngN/gh)", title="T2- Select Potential Nitrification Rate", x="Genotypes")+
  theme_classic()+
  scale_fill_manual(values=c("#FE6100", "#DC267F", "#785EF0","#FFB000", "#FFB000", "#FFB000", "#FFB000"))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


#Denitrfication Overall Figure 3====

Master = read.csv(file ="MasterFunction.DataSheet.Corrected.3.14.21.csv", head = T)



ggplot(Master, aes(x = fct_reorder(Pedigree, HelTotal.N2O.ng.g.hr, .fun = median, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Pedigree))+
  geom_boxplot()+labs(y="Overall Denitrfication(ngN/ghr)", title=" Overall Denitrfication (N2O+N2)", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")

NitrT1=subset(Master, Time== "T1")
NitrT2=subset(Master, Time== "T2")

NitrT2=subset(NitrT2, Sample.ID!="NIL.T2.173")



n.dataframe <-aov(log10(NitrT2$HelTotal.N2O.ng.g.hr) ~(NitrT2$Pedigree))
summary(n.dataframe)
TukeyHSD(n.dataframe)

n.dataframe <-aov(log(NitrT2$HelTotal.N2O.ng.g.hr) ~(NitrT2$Pedigree))
summary(n.dataframe)
TukeyHSD(n.dataframe)



colnames(NitrT1)
ggplot(NitrT1, aes(x = fct_reorder(Pedigree, HelTotal.N2O.ng.g.hr, .fun = median, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Pedigree))+
  geom_boxplot()+labs(y="Overall Denitrfication(ngN/ghr)", title="Timepoint 1 - Overall Denitrfication (N2O+N2)", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")


SelectT1=subset(NitrT1, Pedigree== "Z031E0035"| Pedigree== "Z031E0578"| Pedigree== "B73xPI384071"|
                    Pedigree== "Z031E0028"| Pedigree== "Z031E0021"| Pedigree== "B73"
                  | Pedigree== "PI384071" )


colnames(SelectT1)
ggplot(SelectT1, aes(x = fct_reorder(Pedigree, HelTotal.N2O.ng.g.hr, .fun = median, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Pedigree))+
  geom_boxplot()+labs(y="Potential Overall Denitrfication(ngN/ghr)", title="T1- Select Potential Overall Denitrfication(ngN/ghr)", x="Genotypes")+
  theme_classic()+scale_fill_manual(values=c("blue", "limegreen", "yellow","lightblue", "lightblue", "lightblue", "lightblue"))+
  geom_jitter( width=0.20)

##Timepoint 2

colnames(NitrT2)
ggplot(NitrT2, aes(x = fct_reorder(Pedigree, HelTotal.N2O.ng.g.hr, .fun = median, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Pedigree))+
  geom_boxplot()+labs(y="Overall Denitrfication(ngN/ghr)", title="Timepoint 2 - Overall Denitrfication (N2O+N2)", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")


SelectT2=subset(NitrT2, Pedigree== "Z031E0591"| Pedigree== "Z031E0012"| Pedigree== "B73xPI384071"|
                  Pedigree== "Z031E0047"| Pedigree== "Z031E0042"| Pedigree== "B73"
                | Pedigree== "PI384071" )


colnames(SelectT2)

###Mansucript Figure=====
ggplot(SelectT2, aes(x = fct_reorder(Pedigree, log1p(HelTotal.N2O.ng.g.hr), .fun = median, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Pedigree))+
  geom_boxplot()+labs(y="Potential Overall Denitrfication(ngN/ghr)", x="Genotypes")+
  theme_classic()+scale_fill_manual(values=c("#FE6100", "#DC267F", "#785EF0","#FFB000", "#FFB000", "#FFB000", "#FFB000"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


n.dataframe <-aov((SelectT2$HelTotal.N2O.ng.g.hr) ~(SelectT2$Pedigree))

summary(n.dataframe)

TukeyHSD(n.dataframe)

colnames(NitrT2)


Selectnitf=subset(NitrT2, Pedigree== "B73"
                  | Pedigree== "PI384071" )
n.dataframe <-aov((Selectnitf$HelTotal.N2O.ng.g.hr) ~(Selectnitf$Pedigree))
summary(n.dataframe)
#1.078  0.339

Selectnitf=subset(NitrT2,  Pedigree== "Z031E0012"| Pedigree== "B73")
n.dataframe <-aov((Selectnitf$HelTotal.N2O.ng.g.hr) ~(Selectnitf$Pedigree))
summary(n.dataframe)
# 56.54   3.053  0.131
Selectnitf=subset(NitrT2,  Pedigree== "Z031E0591"| Pedigree== "B73")
n.dataframe <-aov((Selectnitf$HelTotal.N2O.ng.g.hr) ~(Selectnitf$Pedigree))
summary(n.dataframe)

# 56.54   3.053  0.131


Selectnitf=subset(NitrT2, Pedigree== "B73"
                  | Pedigree== "B73xPI384071" )
n.dataframe <-aov(log1p(Selectnitf$HelTotal.N2O.ng.g.hr) ~(Selectnitf$Pedigree))
summary(n.dataframe)
# 0.02802   0.187   0.68

###Incomplete Denitrfication =====
Master = read.csv(file ="MasterFunction.DataSheet.Corrected.3.14.21.csv", head = T)

ggplot(Master, aes(x = fct_reorder(Pedigree, AceTotal.N2O...N2, .fun = median, .desc=TRUE), y=AceTotal.N2O...N2, fill=Pedigree))+
  geom_boxplot()+labs(y="Incomplete Denitrfication(ngN/ghr)", title=" Incomplete Denitrfication (N2O)", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")


NitrT1=subset(Master, Time== "T1")
NitrT2=subset(Master, Time== "T2")

colnames(NitrT1)
ggplot(NitrT1, aes(x = fct_reorder(Pedigree, AceTotal.N2O...N2, .fun = median, .desc=TRUE), y=AceTotal.N2O...N2, fill=Pedigree))+
  geom_boxplot()+labs(y="Incomplete Denitrfication(ngN/ghr)", title="Timepoint 1 - Incomplete Denitrfication (N2O)", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")


SelectT1=subset(NitrT1, Pedigree== "Z031E0068"| Pedigree== "Z031E0578"| Pedigree== "B73xPI384071"|
                  Pedigree== "Z031E0536"| Pedigree== "Z031E0560"| Pedigree== "B73"
                | Pedigree== "PI384071" )


colnames(SelectT1)
ggplot(SelectT1, aes(x = fct_reorder(Pedigree, AceTotal.N2O...N2, .fun = median, .desc=TRUE), y=AceTotal.N2O...N2, fill=Pedigree))+
  geom_boxplot()+labs(y="Potential Incomplete Denitrfication(ngN/ghr)", title="T1- Select Potential Incomplete Denitrfication(ngN/ghr)", x="Genotypes")+
  theme_classic()+scale_fill_manual(values=c("blue", "limegreen", "yellow","lightblue", "lightblue", "lightblue", "lightblue"))+
  geom_jitter( width=0.20)

##Timepoint 2

colnames(NitrT2)
ggplot(NitrT2, aes(x = fct_reorder(Pedigree, AceTotal.N2O...N2, .fun = median, .desc=TRUE), y=AceTotal.N2O...N2, fill=Pedigree))+
  geom_boxplot()+labs(y="Potential Incomplete Denitrfication(ngN/ghr)", title="Timepoint 2 - Incomplete Denitrfication (N2O)", x="Genotypes")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")


SelectT2=subset(NitrT2, Pedigree== "Z031E0523"| Pedigree== "Z031E0012"| Pedigree== "B73xPI384071"|
                  Pedigree== "Z031E0021"| Pedigree== "Z031E0556"| Pedigree== "B73"
                | Pedigree== "PI384071" )


colnames(SelectT2)

###Manuscript Figure=====

NIL.T2.085

SelectT2=subset(SelectT2, Sample.ID!="NIL.T2.085")

ggplot(SelectT2, aes(x = fct_reorder(Pedigree, AceTotal.N2O...N2, .fun = median, .desc=TRUE), y=AceTotal.N2O...N2, fill=Pedigree))+
  geom_boxplot()+labs(y="Potential Incomplete Denitrfication(ngN/ghr)", x="Genotypes")+
  theme_classic()+
  scale_fill_manual(values=c("#FE6100", "#DC267F", "#785EF0","#FFB000", "#FFB000", "#FFB000", "#FFB000"))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


Selectnitf=subset(SelectT2,  Pedigree== "PI384071"| Pedigree== "B73")
n.dataframe <-aov((Selectnitf$AceTotal.N2O...N2) ~(Selectnitf$Pedigree))
summary(n.dataframe)
#  0.003  0.957
Selectnitf=subset(SelectT2,  Pedigree== "B73xPI384071"| Pedigree== "B73")
n.dataframe <-aov((Selectnitf$AceTotal.N2O...N2) ~(Selectnitf$Pedigree))
summary(n.dataframe)

#  5.10   0.079  0.788
Selectnitf=subset(SelectT2,  Pedigree== "Z031E0012"| Pedigree== "B73")
n.dataframe <-aov((Selectnitf$AceTotal.N2O...N2) ~(Selectnitf$Pedigree))
summary(n.dataframe)
# 44.39   1.001  0.356
Selectnitf=subset(SelectT2,  Pedigree== "Z031E0523"| Pedigree== "B73")
n.dataframe <-aov((Selectnitf$AceTotal.N2O...N2) ~(Selectnitf$Pedigree))
summary(n.dataframe)
#115.90    2.52  0.173
R

Selectnitf=subset(SelectT2, Pedigree== "B73"
                  | Pedigree== "B73xPI384071" )
n.dataframe <-aov(log1p(Selectnitf$AceTotal.N2O...N2) ~(Selectnitf$Pedigree))
summary(n.dataframe)

# 0.03753   0.162  0.701



#Core analysis on B73, Hybrid and teo
CoreNitf = read.csv(file ="NIL18.All.Nitrfication.DATA.OUTlE.Core2.csv", head = T)

n.dataframe <-aov(CoreNitf$Log ~(CoreNitf$Pedigree))
summary(n.dataframe)
TukeyHSD(n.dataframe)

boxplot(CoreNitf$Nitrifcation.ngN.gh ~(CoreNitf$Pedigree))
plot(density(CoreNitf$Nitrifcation.ngN.gh))
lines(density(CoreNitf$Pedigree))

#Tukeys range test doesn't seem to useful for this 

#Plot.
ggplot(CoreNitf, aes(x = Nitrifcation.ngN.gh, fill = Pedigree)) + geom_density(alpha = 0.5)


#Aging per mean
Pedigree.mean.results= aggregate(NitrT2[, 6], list(NitrT2$Pedigree), mean)



#Subset and rerun
NitrCoreT1=subset(CoreNitf, Time== "T1")
NitrCoreT2=subset(CoreNitf, Time== "T2")

boxplot(NitrCoreT1$Nitrifcation.ngN.gh ~(NitrCoreT1$Pedigree))
ggplot(NitrCoreT1, aes(x = Nitrifcation.ngN.gh, fill = Pedigree)) + geom_density(alpha = 0.5)

n.dataframe <-aov(NitrCoreT1$Nitrifcation.ngN.gh ~(NitrCoreT1$Pedigree))
summary(n.dataframe)
TukeyHSD(n.dataframe)


boxplot(NitrCoreT2$Nitrifcation.ngN.gh ~(NitrCoreT2$Pedigree))
ggplot(NitrCoreT2, aes(x = Nitrifcation.ngN.gh, fill = Pedigree)) + geom_density(alpha = 0.5)

n.dataframe <-aov(NitrCoreT2$Nitrifcation.ngN.gh ~(NitrCoreT2$Pedigree))
summary(n.dataframe)
TukeyHSD(n.dataframe)



####Denitrification Data =====


#DEA Results 
NILDEA = read.csv(file ="NIL18.T2.DEA.csv", head = T)
NILDEA = read.csv(file ="NIL18.T1.DEA.csv", head = T)
NILDEA = read.csv(file ="NIL18.T2.DEA.ALt.csv", head = T)


#Running a rough analysis 
# Denitrification: NIL 0047, NIL 0021, NIL 0545, B73, teosinte, NIL 0578, NIL 0057, NIL 0071
Selectdea=subset(NILDEA, Genotype== "Z031E0578"| Genotype== "Z031E0057"| Genotype== "Z031E0071"|
                   Genotype== "Z031E0047"| Genotype== "Z031E0021"| Genotype== "B73"
                 | Genotype== "PI384071" |Genotype== "Z031E0545" )

Selectdea=subset(Selectnitf, Sample.ID!="NIL18.T2.152")


#Hel
n.dataframe <-aov(NILDEA$HelTotal.N2O.ng.g.hr ~(NILDEA$Genotype))
summary(n.dataframe)

#Ace
n.dataframe <-aov(NILDEA$AceTotal.N2O...N2 ~(NILDEA$Genotype))
summary(n.dataframe)
TukeyHSD(n.dataframe)

##Significant diff

#Difference 
n.dataframe <-aov(NILDEA$Evolved.N2.ng.g.hr ~(NILDEA$Genotype))
summary(n.dataframe)

#Now I want to look at the data 
boxplot(NILDEA$HelTotal.N2O.ng.g.hr ~(NILDEA$Genotype))
boxplot(NILDEA$AceTotal.N2O...N2 ~(NILDEA$Genotype))
boxplot(NILDEA$Evolved.N2.ng.g.hr ~(NILDEA$Genotype))


#Isolating the means values
Hel.mean.results= aggregate(NILDEA[, 6], list(NILDEA$Genotype), mean)
N2.mean.results= aggregate(NILDEA[, 7], list(NILDEA$Genotype), mean)
Ace.mean.results= aggregate(NILDEA[, 8], list(NILDEA$Genotype), mean)

Hel.mean.results= aggregate(NILDEA[, 6], list(NILDEA$Genotype), var)
N2.mean.results= aggregate(NILDEA[, 7], list(NILDEA$Genotype), mean)

#Making a Nice Figure for DEA ======

ggplot(NILDEA, aes(x = HelTotal.N2O.ng.g.hr, fill = Genotype)) + geom_density(alpha = 0.5)

ggplot(NILDEA, aes(x = fct_reorder(Genotype, HelTotal.N2O.ng.g.hr, .fun = median, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Genotype))+
  geom_boxplot()+labs(y="ngN/gh", title="Potential Denitrification Rate", x="Genotypes")+
  theme_classic()
ggplot(NILDEA, aes(x = fct_reorder(Genotype, HelTotal.N2O.ng.g.hr, .fun = median, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Genotype))+
  geom_boxplot()+labs(y="ngN/gh", title="Potential Denitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(Selectdea, aes(x = fct_reorder(Genotype, HelTotal.N2O.ng.g.hr, .fun = median, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Genotype))+
  geom_boxplot()+labs(y="ngN/gh", title="Potential Denitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()+scale_fill_manual(values=c("blue", "limegreen", "lightblue","lightblue", "lightblue", "lightblue", "lightblue","lightblue"))


 ggplot(NitrCoreT1, aes(x = Nitrifcation.ngN.gh, fill = Pedigree)) + geom_density(alpha = 0.5)
 ggplot(NitrCoreT2, aes(x = Nitrifcation.ngN.gh, fill = Pedigree)) + geom_density(alpha = 0.5)
 
 ###Caluating mean and start devation sfor the new set of maps outliers removed. =====
 

NILNitf = read.csv(file ="NIL18.All.Nitrfication.DATA.OUTlIE.csv", head = T)
 
 
Means<-aggregate(NILNitf[, 6:7], list(NILNitf$Pedigree, NILNitf$Time), mean)
# colnames(Means) #Code to understand header names
#Change col names
names(Means)[names(Means) == "Nitrifcation.ngN.gh"] <- "MeanNitrfication"
names(Means)[names(Means) == "Log"] <- "MeanLogNit"

SD<-aggregate(NILNitf[, 6:7], list(NILNitf$Pedigree, NILNitf$Time), sd)
colnames(SD)
names(SD)[names(SD) == "Nitrifcation.ngN.gh"] <- "SD.Nitrfication"
names(SD)[names(SD) == "Log"] <- "SD.LogNit"


Allruns=rbind(Means, SD)
Head<-merge(Means, SD) 
write.csv(Head,"Nitrfication.GenotypeMean.SD.csv", row.names = FALSE)

Time2Genotyp=subset(Head, Group.2== "T2")
write.csv(Time2Genotyp,"Nitrfication.Genotype.Time2.Mean.SD.csv", row.names = FALSE)

###Plotting the full scope of function

Master=subset(Master, Sample.ID!="NIL18.T2.152")
Master=subset(Master, Sample.ID!="NIL18.T2.009")
Master=subset(Master, Sample.ID!="NIL18.T2.083")


Master.t1<-subset(Master, Time="T1")
Master.t2<-subset(Master, Time="T2")

ggplot(Master, aes(x = fct_reorder(Pedigree, Nitrifcation.ngN.gh, .fun = median, .desc=TRUE), y=Nitrifcation.ngN.gh, fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="BNI", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(Master.t1, aes(x = fct_reorder(Pedigree, Nitrifcation.ngN.gh, .fun = median, .desc=TRUE), y=Nitrifcation.ngN.gh, fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="Time 1 Nitrfication Rates", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(Master.t2, aes(x = fct_reorder(Pedigree, Nitrifcation.ngN.gh, .fun = median, .desc=TRUE), y=Nitrifcation.ngN.gh, fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="Time 2 Nitrfication Rates", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

##N20 Production
ggplot(Master, aes(x = fct_reorder(Pedigree, HelTotal.N2O.ng.g.hr, .fun = mean, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="HelTotal.N2O.ng.g.hr", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(Master.t1, aes(x = fct_reorder(Pedigree, HelTotal.N2O.ng.g.hr, .fun = mean, .desc=TRUE), y=HelTotal.N2O.ng.g.hr, fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="Time 1 HelTotal.N2O.ng.g.hr", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(Master.t2, aes(x = fct_reorder(Pedigree, (HelTotal.N2O.ng.g.hr), .fun = median, .desc=TRUE), y= (HelTotal.N2O.ng.g.hr), fill=Pedigree))+
  geom_boxplot()+labs(y="ngN/gh", title="Time 2 HelTotal.N2O.ng.g.hr", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))




##ASREML-R Code ====


library(asreml) 
library(asremlPlus)
library(tidyverse)

Master = read.csv(file ="MasterFunction.DataSheet.csv", head = T)
Master = read.csv(file ="MasterFunction.DataSheet.Corrected.3.14.21.csv", head = T)

head(Master)
Function.Data<-Master%>%unite(Block,c(Range,Row), remove=FALSE)

Function.Data$Block <- as.factor(Function.Data$Block)
Function.Data$Time <- as.factor(Function.Data$Time)
Function.Data$Range <- as.factor(Function.Data$Range)
Function.Data$Row <- as.factor(Function.Data$Row)
Function.Data$Pedigree <- as.factor(Function.Data$Pedigree)

Function.Data2<-subset(Function.Data, Pedigree!="Bulk")
Function.Data2<-subset(Function.Data2, Pedigree!="Bulk ")

Function.Data2.nit<-Function.Data2[!is.na(Function.Data2$Nitrification.Rate..ugN.ghr.),]
Function.Data2.hel<-Function.Data2[!is.na(Function.Data2$HelTotal.N2O.ng.g.hr),]
Function.Data2.ace<-Function.Data2[!is.na(Function.Data2$AceTotal.N2O...N2),]


Function.Data2<-subset(Function.Data, Pedigree=="B73"|Pedigree=="Z031E0578"|Pedigree=="Z031E0028")
# 
# NitrificationModel<-asreml(data=Function.Data2,
#                            fixed=(Nitrifcation.ngN.gh)~1,
#                            random= ~~Pedigree+ar1(Range):ar1(Row),
#                            na.action = na.method(x = 'fail', y = 'omit'))
# 
# wald(NitrificationModel)

colnames(Function.Data2)
##Here is the model I used in the 2017 Data set|  But I'm coming across some singularity issues
#Which seems to mean Im overspecified the data
NitrificationModel<-asreml(data=Function.Data2,
                           fixed=log1p(Nitrification.Rate..ugN.ghr.)~Pedigree,
                           random= ~1+Time+Row+Range+Block,
                           na.action = na.method(x = 'fail', y = 'omit'))

wald(NitrificationModel)


#Function.Data2T<-na.omit(Function.Data2$Nitrifcation.ngN.gh)

##IT has something to do with the data strucutre be
NitrificationModel<-asreml(data=Function.Data2,
                           fixed=log1p(Nitrification.Rate..ugN.ghr.)~1,
                           random= ~1+Pedigree+Time+Row+Range+Block,
                           na.action = na.method(x = 'fail', y = 'omit'))

wald(NitrificationModel)

plot(NitrificationModel)
wald(NitrificationModel)
summary(NitrificationModel)
vpredict(NitrificationModel, h2~(V4)/(V1+V2+V3+V4))

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    903.43        2690.64    <2e-16 ***
#   Pedigree      43     13.91          41.42    0.5402    
# residual (MS)         0.34                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

colnames(Function.Data2.nit)
NitrificationModel<-asreml(data=Function.Data2.nit,
                           fixed=Nitrification.Rate..ugN.ghr.~Pedigree,
                           random= ~1+(Time)+(Range)+(Row))

wald(NitrificationModel)


#Predict Means
predict(NitrificationModel, classify="Pedigree",sed=T)

#This Data frame has the LS means and the standard error
LSMeans<-predict(NitrificationModel, classify="Pedigree")
LSMeans<-LSMeans$pvals
ggplot(LSMeans, aes(x=reorder(Pedigree,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error*1.96, ymax=predicted.value+std.error*1.96), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    380391         18.252 1.935e-05 ***
#   Pedigree      43    991643         47.581    0.2916    
# residual (MS)        20841                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##DEA Models: Issue occured with Nitrfication
N.asr.DEA<-asreml(data=Function.Data2.hel,
                           fixed=log1p(HelTotal.N2O.ng.g.hr)~1,
                           random= ~Pedigree+(Range)+(Row)+Block,
                           na.action = na.method(x = 'fail', y = 'omit'),
                           ai.loadings=0)

plot(N.asr.DEA)
summary(N.asr.DEA)
summary.asreml(N.asr.DEA,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.asr.DEA)
wald(N.asr.DEA)


wald(N.asr.DEA)
summary(N.asr.DEA)
vpredict(N.asr.DEA, h2~(V1)/(V1+V2+V3+V4))

##Model could not reach singlularity 
N.asr.DEA$loglik

##Looking at the fixed and random coffeftne
round(N.asr.DEA$coefficients$fixed,2)
round(head(N.asr.DEA$coefficients$random, 10),2)

#Variance compenets of error
N.asr.DEA$sigma2
N.asr.DEA$vparameters

# We can get the variance components estimates by multiplying the gamma estimates by the residual error variance:
N.asr.DEA$sigma2*N.asr.DEA$vparameters


summary(N.asr.DEA)$varcomp

#Wal Test 
wald(N.asr.DEA, denDF = "numeric")
plot(varioGram(N.asr.DEA))


N.asr.DEA<-asreml(data=Function.Data2.hel,
                  fixed=log1p(HelTotal.N2O.ng.g.hr)~Pedigree,
                  random= ~(Time)+ar1(Range):ar1(Row),
                  na.action = na.method(x = 'fail', y = 'omit'))


plot(N.asr.DEA)
summary(N.asr.DEA)
summary.asreml(N.asr.DEA,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.asr.DEA)
wald(N.asr.DEA)

#Predict Means
predict(N.asr.DEA, classify="Pedigree",sed=T)

#This Data frame has the LS means and the standard error
LSMeans<-predict(N.asr.DEA, classify="Pedigree")
LSMeans<-LSMeans$pvals

ggplot(LSMeans, aes(x=reorder(Pedigree,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error*1.96, ymax=predicted.value+std.error*1.96), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


###Preparing the data to look at introgressions===
Function.Data<-subset(Function.Data, Pedigree!="Bulk")
Function.Data<-subset(Function.Data, Pedigree!="B73xPI384071")
Function.Data<-subset(Function.Data, Pedigree!="PI384071")

Function.Data.t1<-subset(Function.Data2.hel, Time=="T1")
Function.Data.t2<-subset(Function.Data2.hel, Time=="T2")



##Time 1 DEA HEL====
Function.Data<-subset(Function.Data, Pedigree!="Bulk")
Function.Data<-subset(Function.Data, Pedigree!="Bulk")
Function.Data<-subset(Function.Data, Pedigree!="Bulk")

Function.Data.t1<-subset(Function.Data, Time=="T1")
Function.Data.t2<-subset(Function.Data, Time=="T2")

N.DEA.t1<-asreml(data=Function.Data.t1,
                  fixed=log1p(HelTotal.N2O.ng.g.hr)~Pedigree,
                  random= ~+ar1(Range):ar1(Row),
                  na.action = na.method(x = 'fail', y = 'omit'))


plot(N.DEA.t1)
summary(N.DEA.t1)
summary.asreml(N.DEA.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.DEA.t1)
wald(N.DEA.t1)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   13.3095         63.896 1.332e-15 ***
#   Pedigree      43    8.9015         42.735    0.4827    
# residual (MS)       0.2083                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Predict Means
predict(N.DEA.t1, classify="Pedigree",sed=T)

#This Data frame has the LS means and the standard error
LSMeans<-predict(N.DEA.t1, classify="Pedigree")
LSMeans<-LSMeans$pvals
write.csv(LSMeans,"Lsmeans/HEL.DEA.T1.csv", row.names = FALSE)

#dir.create("Lsmeans")

ggplot(LSMeans, aes(x=reorder(Pedigree,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error*1.96, ymax=predicted.value+std.error*1.96), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



##Time 2 DEA HEL====
Function.Data<-subset(Function.Data, Pedigree!="Bulk")

Function.Data.t1<-subset(Function.Data, Time=="T1")
Function.Data.t2<-subset(Function.Data, Time=="T2")

N.DEA.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(HelTotal.N2O.ng.g.hr)~Pedigree,
                 random= ~+ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))

N.DEA.t2<-asreml(data=Function.Data.t2,
                 fixed=(HelTotal.N2O.ng.g.hr)~Pedigree,
                 random= ~+ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))

plot(N.DEA.t2)
summary(N.DEA.t2)
summary.asreml(N.DEA.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.DEA.t2)
wald(N.DEA.t2)


# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    46.466         336.75 < 2.2e-16 ***
#   Pedigree      43    10.095          73.16  0.002786 ** 
#   residual (MS)        0.138                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Predict Means
predict(N.DEA.t2, classify="Pedigree",sed=T)



#This Data frame has the LS means and the standard error
LSMeans<-predict(N.DEA.t2, classify="Pedigree")
LSMeans<-LSMeans$pvals
write.csv(LSMeans,"Lsmeans/HEL.DEA.T2.csv", row.names = FALSE)


ggplot(LSMeans, aes(x=reorder(Pedigree,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error*1.96, ymax=predicted.value+std.error*1.96), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##Low Is Z031E071,E0591, E070, 

##Time 1 DEA ACE====
Function.Data<-subset(Function.Data, Pedigree!="Bulk")

Function.Data.t1<-subset(Function.Data2.ace, Time=="T1")
Function.Data.t2<-subset(Function.Data2.ace, Time=="T2")

N.DEA.t1<-asreml(data=Function.Data.t1,
                 fixed=(AceTotal.N2O...N2)~Pedigree,
                 random= ~+ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.DEA.t1)
summary(N.DEA.t1)
summary.asreml(N.DEA.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.DEA.t1)
wald(N.DEA.t1)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   1.67241         121.82 < 2.2e-16 ***
#   Pedigree      43   1.09625          79.85 0.0005432 ***
#   residual (MS)      0.01373                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Predict Means
predict(N.DEA.t1, classify="Pedigree",sed=T)

#This Data frame has the LS means and the standard error
LSMeans<-predict(N.DEA.t1, classify="Pedigree")
LSMeans<-LSMeans$pvals

write.csv(LSMeans,"Lsmeans/ACE.DEA.T1.csv", row.names = FALSE)


ggplot(LSMeans, aes(x=reorder(Pedigree,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error*1.96, ymax=predicted.value+std.error*1.96), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



##Time 2 DEA ACE====
Function.Data<-subset(Function.Data, Pedigree!="Bulk")

Function.Data.t1<-subset(Function.Data, Time=="T1")
Function.Data.t2<-subset(Function.Data, Time=="T2")

N.DEA.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(AceTotal.N2O...N2)~Pedigree,
                 random= ~+ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))
N.DEA.t2<-asreml(data=Function.Data.t2,
                 fixed=(AceTotal.N2O...N2)~Pedigree,
                 random= ~+ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.DEA.t2)
summary(N.DEA.t2)
summary.asreml(N.DEA.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.DEA.t2)
wald(N.DEA.t2)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    98.818        1028.85    <2e-16 ***
#   Pedigree      43     5.280          54.98    0.1041    
# residual (MS)        0.096                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

predict(N.DEA.t2, classify="Pedigree",sed=T)



#This Data frame has the LS means and the standard error
LSMeans<-predict(N.DEA.t2, classify="Pedigree")
LSMeans<-LSMeans$pvals
write.csv(LSMeans,"Lsmeans/ACE.DEA.T2.csv", row.names = FALSE)


ggplot(LSMeans, aes(x=reorder(Pedigree,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error, ymax=predicted.value+std.error), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





##Time 1 NIT====
Function.Data<-subset(Function.Data, Pedigree!="Bulk")


Function.Data.t1<-subset(Function.Data, Time=="T1")
Function.Data.t2<-subset(Function.Data, Time=="T2")


Function.Data.t1<-subset(Function.Data2, Time=="T1")
Function.Data.t2<-subset(Function.Data2, Time=="T2")


N.NIT.t1<-asreml(data=Function.Data.t1,
                 fixed=log1p(Nitrifcation.ngN.gh)~Pedigree,
                 random= ~+ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))

N.NIT.t1<-asreml(data=Function.Data.t1,
                 fixed=(Nitrification.Rate..ngN.ghr./1000)~Pedigree,
                 random= ~+ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t1)
summary(N.NIT.t1)
summary.asreml(N.NIT.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t1)
wald(N.NIT.t1)
##THis wa the log1p model
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    4897.0        27807.7    <2e-16 ***
#   Pedigree      43       6.4           36.4    0.7509    
# residual (MS)          0.2                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Nit/1000 model
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   2009801        273.482    <2e-16 ***
#   Pedigree      43    314687         42.821     0.479    
# residual (MS)         7349                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Predict Means
predict(N.NIT.t1, classify="Pedigree",sed=T)

#This Data frame has the LS means and the standard error
LSMeans<-predict(N.NIT.t1, classify="Pedigree")
LSMeans<-LSMeans$pvals
#write.csv(LSMeans,"Lsmeans/NIT.T1.csv", row.names = FALSE)



ggplot(LSMeans, aes(x=reorder(Pedigree,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error*1.96, ymax=predicted.value+std.error*1.96), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



##Time 2 NIT====
Function.Data<-subset(Function.Data, Pedigree!="Bulk")

Function.Data.t1<-subset(Function.Data, Time=="T1")
Function.Data.t2<-subset(Function.Data, Time=="T2")

#Outliers||| Removing 1 from Teo and 1 from Z31E0031
Function.Data.t2=subset(Function.Data.t2, X.SampleID!="NIL.T2.152")
Function.Data.t2=subset(Function.Data.t2, X.SampleID!="NIL.T2.009")



N.NIT.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(Nitrification.Rate..ngN.ghr.)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


N.NIT.t2<-asreml(data=Function.Data.t2,
                 fixed=(Nitrification.Rate..ngN.ghr./1000)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))
plot(N.NIT.t2)
summary(N.NIT.t2)
summary.asreml(N.NIT.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t2)
wald(N.NIT.t2)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   26826.1          79479 < 2.2e-16 ***
#   Pedigree      43      23.9             71  0.004753 ** 
#   residual (MS)          0.3                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Predict Means
predict(N.NIT.t2, classify="Pedigree",sed=T)

#This Data frame has the LS means and the standard error
LSMeans<-predict(N.NIT.t2, classify="Pedigree")
LSMeans<-LSMeans$pvals
# write.csv(LSMeans,"Lsmeans/NIT.T2.nooutlier.csv", row.names = FALSE)

ggplot(LSMeans, aes(x=reorder(Pedigree,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error, ymax=predicted.value+std.error), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Z031E0047 Z031E0071 low Z031E0021
ggplot(Function.Data.t2, aes(x=reorder(Pedigree,Nitrification.Rate..ugN.ghr.,median), y=(Nitrification.Rate..ugN.ghr.+10)))+ geom_boxplot()+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(Function.Data.t2, aes(x=reorder(Pedigree,Nitrification.Rate..ugN.ghr.,median), y=log(Nitrification.Rate..ugN.ghr.+10)))+ geom_boxplot()+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##Great at least this replicates our model results 
##There seems to be some pattern switching between the two sample sets.
##At T1 B73 seems to have some sus traits but this pattern switchs later in the data set

##Functional Gene Models =====

##Timepoint 1=
N.NIT.t1<-asreml(data=Function.Data.t1,
                 fixed=log1p(amoA..copies.ul.)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t1)
summary(N.NIT.t1)
summary.asreml(N.NIT.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t1)
wald(N.NIT.t1)


# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    7729.9         7850.2    <2e-16 ***
#   Pedigree      43      39.5           40.1     0.597    
# residual (MS)          1.0                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

N.NIT.t1<-asreml(data=Function.Data.t1,
                 fixed=log1p(bamoA..copies.ul.)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t1)
summary(N.NIT.t1)
summary.asreml(N.NIT.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t1)
wald(N.NIT.t1)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   1684.65        2039.26    <2e-16 ***
#   Pedigree      43     30.91          37.42    0.7113    
# residual (MS)         0.83                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


N.NIT.t1<-asreml(data=Function.Data.t1,
                 fixed=log1p(X.4)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t1)
summary(N.NIT.t1)
summary.asreml(N.NIT.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t1)
wald(N.NIT.t1)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    463.20         3406.9    <2e-16 ***
#   Pedigree      43      6.66           49.0    0.2444    
# residual (MS)         0.14                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


N.NIT.t1<-asreml(data=Function.Data.t1,
                 fixed=log1p(X.7)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t1)
summary(N.NIT.t1)
summary.asreml(N.NIT.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t1)
wald(N.NIT.t1)
#singularity not reached

N.NIT.t1<-asreml(data=Function.Data.t1,
                 fixed=log1p(X.10)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t1)
summary(N.NIT.t1)
summary.asreml(N.NIT.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t1)
wald(N.NIT.t1)
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   1232.35        1969.53    <2e-16 ***
#   Pedigree      43     27.00          43.16    0.4646    
# residual (MS)         0.63                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

N.NIT.t1<-asreml(data=Function.Data.t1,
                 fixed=log1p(X.13)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t1)
summary(N.NIT.t1)
summary.asreml(N.NIT.t1,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t1)
wald(N.NIT.t1)
#singulaity not reached


##Timepoint 2=
N.NIT.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(amoA..copies.ul.)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t2)
summary(N.NIT.t2)
summary.asreml(N.NIT.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t2)
wald(N.NIT.t2)
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1     58.09         57.541 3.308e-14 ***
#   Pedigree      43     25.73         25.487    0.9845    
# residual (MS)         1.01                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

N.NIT.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(bamoA..copies.ul.)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t2)
summary(N.NIT.t2)
summary.asreml(N.NIT.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t2)
wald(N.NIT.t2)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    30.233         33.305 7.877e-09 ***
#   Pedigree      43    38.708         42.641    0.4868    
# residual (MS)        0.908                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


N.NIT.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(X.4)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t2)
summary(N.NIT.t2)
summary.asreml(N.NIT.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t2)
wald(N.NIT.t2)

# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1   18.7638         21.055 4.463e-06 ***
#   Pedigree      43   27.6259         30.999    0.9139    
# residual (MS)       0.8912                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


N.NIT.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(X.7)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t2)
summary(N.NIT.t2)
summary.asreml(N.NIT.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t2)
wald(N.NIT.t2)
#singularity not reached
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    42.044         35.930 2.046e-09 ***
#   Pedigree      43    56.006         47.862     0.282    
# residual (MS)        1.170                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > 
N.NIT.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(X.10)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t2)
summary(N.NIT.t2)
summary.asreml(N.NIT.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t2)
wald(N.NIT.t2)
# Df Sum of Sq Wald statistic Pr(Chisq)  
# (Intercept)    1     8.818          5.488   0.01914 *
#   Pedigree      43    77.917         48.497   0.26093  
# residual (MS)        1.607                           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

N.NIT.t2<-asreml(data=Function.Data.t2,
                 fixed=log1p(X.13)~Pedigree,
                 random= ~ar1(Range):ar1(Row),
                 na.action = na.method(x = 'fail', y = 'omit'))


plot(N.NIT.t2)
summary(N.NIT.t2)
summary.asreml(N.NIT.t2,param = c("sigma", "gamma"),vparameters = TRUE)

names(N.NIT.t2)
wald(N.NIT.t2)
# Df Sum of Sq Wald statistic Pr(Chisq)    
# (Intercept)    1    36.419         42.415 7.381e-11 ***
#   Pedigree      43    47.068         54.817    0.1068    
# residual (MS)        0.859                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Non of the qPCR data worked 
###Denitification means====
NILDEA = read.csv(file ="NIL18.T2.DEA.csv", head = T)


Means<-aggregate(NILDEA[, 6:8], list(NILDEA$Genotype, NILDEA$Time), mean)
 colnames(Means) #Code to understand header names
#Change col names
names(Means)[names(Means) == "HelTotal.N2O.ng.g.hr"] <- "HelTotal.N2O.ng.g.hr.M"
names(Means)[names(Means) == "Evolved.N2.ng.g.hr"] <- "Evolved.N2.ng.g.hr.M"
names(Means)[names(Means) == "AceTotal.N2O...N2"] <- "AceTotal.N2O...N2.M"

SD<-aggregate(NILDEA[, 6:8], list(NILDEA$Genotype, NILDEA$Time), sd)
colnames(SD)
names(SD)[names(SD) == "HelTotal.N2O.ng.g.hr"] <- "SD.HelTotal.N2O.ng.g.hr"
names(SD)[names(SD) == "Evolved.N2.ng.g.hr"] <- "SD.Evolved.N2.ng.g.hr"
names(SD)[names(SD) == "AceTotal.N2O...N2"] <- "SD.AceTotal.N2O...N2"


Head<-merge(Means, SD) 
write.csv(Head,"DEA.Mean.SD.T2.csv", row.names = FALSE)


Time2Genotyp=subset(Head, Group.2== "T2")
write.csv(Time2Genotyp,"Nitrfication.Genotype.Time2.Mean.SD.csv", row.names = FALSE)


###N20 emission in Space====
setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL Maize Microbiome /R-Scripts/NIL CSV/Functional Data")
setwd("~/Documents/Thesis /Chapter 4/Science /Rstats/Old R-Analysis /R-Scripts/NIL CSV")
N20emission<- read.csv(file="NIL.N20.Emissions.csv",head=T)
N20emission<- read.csv(file="MasterN20Emissions.csv",head=T)
N20emission<- read.csv(file="MasterN20Emissions2.csv",head=T)

N20emission[,9:14]<-exp(N20emission[,9:14])
logData<-log1p(N20emission[,9:14])


hist(N20emission$Flux.min)


shapiro.test(N20emission$N2O.Concentation)
shapiro.test(N20emission$Flux.min)

qqnorm(N20emission$N2O.Concentation)
qqnorm(N20emission$Flux.min)

#Lets look at the data very coursly
#Boxplot: N20 Per treatment Type====
boxplot(N20emission$N2O.Concentation~ N20emission$Type)
boxplot(N20emission$Slope~ N20emission$Type)
boxplot(N20emission$Flux.min~ N20emission$Type)
boxplot(N20emission$Flux.hr~ N20emission$Type)
boxplot(N20emission$Headspace.Concentration..ppm.~ N20emission$Type)

#Boxplot: N20 Per treatment Time| Important Factor aside from just rate====
boxplot(N20emission$N2O.Concentation~ N20emission$Time)
boxplot(N20emission$Slope~ N20emission$Time)
boxplot(N20emission$Flux.min~ N20emission$Time)
boxplot(N20emission$Flux.hr~ N20emission$Time)
boxplot(N20emission$Headspace.Concentration..ppm.~ N20emission$Time)
boxplot(N20emission$Headspace.Concentration..ppm.~ N20emission$X)

#The Headspace.Concentration..ppm. is the T0 N20 concentration in the sample| I.E. The ambient N20
#This is consistent with the lit. Now what is more suppreising is the rate.
#The data is not super normal, but there is a clear trend in time
model<-lm((N20emission$Headspace.Concentration..ppm.)~ N20emission$X)
summary(model)
plot(model)
plot((N20emission$Headspace.Concentration..ppm.)~ N20emission$X)
abline(model)


model<-lm((N20emission$Flux.min)~ N20emission$X)
summary(model)
plot((N20emission$Flux.min)~ N20emission$X)
abline(model)

#It may be best to use a fitted model to these data. There is a clear relationship,
#but I think the stand curves should be corrected for. There's some week 2 wonky-ness and some patterns driven by the standard curve. 
library(ggplot2)
library(tidyverse)

ggplot(N20emission, aes(X, Flux.min)) +
  geom_point() +
  geom_smooth(method="loess", se = TRUE)+theme_bw()+
  ylab("log of N20 per min")+ggtitle("N20 flux in Time:Gen 11/19/19 ")+xlab("Week")

#Lets color the points by type. This figure is a little more complicated I need to figure out how to plot every week

ggplot(N20emission, aes(X, Flux.min)) +
  geom_point(aes(color=factor(Type), shape=factor(Space))) +
  geom_smooth(method="loess", se = TRUE)+theme_bw()+
  ylab("log of N20 per min")+ggtitle("N20 flux in Time:Gen 11/19/19 ")+xlab("Time")

#Here I want to look at the factor factor differences 
library(tidyverse)

#lets look at how the centeral space point shapes the n20 emisison.===
Far<-N20emission %>%
  subset(Space=="F")%>%
  ggplot(aes(X, Flux.min)) +
  geom_point(aes(color=factor(Type), shape=factor(Space))) +
  geom_smooth(method="loess", se = TRUE)+theme_bw()+
  ylab("log of N20 per min")+ggtitle("Farthest form Plant")+xlab("Time")

Middle<-N20emission %>%
  subset(Space=="M")%>%
  ggplot(aes(X, Flux.min)) +
  geom_point(aes(color=factor(Type), shape=factor(Space))) +
  geom_smooth(method="loess", se = TRUE)+theme_bw()+
  ylab("log of N20 per min")+ggtitle("Middle to Plant")+xlab("Time")

Center<-N20emission %>%
  subset(Space=="C")%>%
  ggplot(aes(X, Flux.min)) +
  geom_point(aes(color=factor(Type), shape=factor(Space))) +
  geom_smooth(method="loess", se = TRUE)+theme_bw()+
  ylab("log of N20 per min")+ggtitle("Closest to Plant")+xlab("Time")



library(gridExtra)
library(grid)
grid.arrange(Center, Middle, Far,nrow = 1)
#Saved at 12x4 inches
# It may be useful to calculate the varaince for each of the of the different treatments 

#
ggplot(N20emission, aes(X, Headspace.Concentration..ppm.)) +
  geom_point() +
  geom_smooth(method="loess", se = TRUE)+theme_bw()+
  ylab("log of N20 per min")+ggtitle("ambient N20 in Time:Gen 11/19/19 ")+xlab("Week")

#This function allows to plot the data in space: After making this plot I'm realzing. 
#This plot makes me realzie that my locational data is a little functed up.
#Tidyverse code I could deff plot these values in a different way====
library(tidyverse)

#this is soo cool. If I pipe out the last bit of the code we can easily manipulate the plot 
Week1<-
  N20emission %>%
  subset(Time== "Week1")%>%
  subset(Type!="Bulk ")%>%
  subset(Type!="Bulk")%>%
  ggplot(aes(RowPlus, Range, size=Flux.hr))+theme_classic()+
  labs(title="N20 Emissions in Space", subtitle = "Week1",ylab="Row")+
  geom_point(aes(color=factor(Type)))+scale_size_area()#position = position_jitter(width =.5,height = 0))#+theme_void()
#I need to figure out how to plot all of these intervals together
library(grid)
library(gridExtra)
#multiplot(Week1, Week2,Week3,Week5,Week5,Week6,Week7,Week8,Week9, nrow=3)
#This figure has scaling issues
grid.arrange(Week1, Week2,Week3,Week5,Week5,Week6,Week7,Week8,Week9,nrow = 3)


# I may want to creat a set of figures for all 9
#This jitter position is closert to my needs. I think setting the jitter may be a little bit more useful

#The different time points, Can probably be prepresened as asome sort of jitter
#This is soooo cool! This allows me to quickly make files and edit them on the fly. 
#I may want to come back to this to see how the values are actually changing in time.

#Note: most of the N20 emissions are occuring early in the growing season and reduce dramatically.

#Boxplot: N20 Per treatment Space====
boxplot(N20emission$N2O.Concentation~ N20emission$Space)
boxplot(N20emission$Slope~ N20emission$Space)
boxplot(N20emission$Flux.min~ N20emission$Space)
boxplot(N20emission$Flux.hr~ N20emission$Space)
boxplot(N20emission$Headspace.Concentration..ppm.~ N20emission$Space)

#Boxplot: N20 Per treatment Sample.name====
boxplot(N20emission$N2O.Concentation~ N20emission$Sample.name)
boxplot(N20emission$Slope~ N20emission$Sample.name)
boxplot(N20emission$Flux.min~ N20emission$Sample.name)
boxplot(N20emission$Flux.hr~ N20emission$Sample.name)
boxplot(N20emission$HeadSample.name.Concentration..ppm.~ N20emission$Sample.name)

#Lets make a liner mixed effects model to see how all of the different factors ====

library(lme4)
library(Matrix)

model<-lmer(Slope~Space, data = N20emission)
#model<-aov(N20emission$Slope~N20emission$Space)
#I like this model of interaction and range and row
model<-aov(Flux.min~Time*Space*Type+Range+Row, data = N20emission)

model<-aov(Flux.min~Time*Space*Type+Range+Row, data = N20emission)
summary(model)
plot(model)
# Df    Sum Sq   Mean Sq F value   Pr(>F)    
# Time                 7 1.906e-08 2.723e-09 2957.33 3.51e-10 ***
#   Space                3 2.680e-10 8.930e-11   97.01 1.79e-05 ***
#   Type                 2 1.890e-10 9.450e-11  102.58 2.29e-05 ***
#   Location            17 1.721e-09 1.012e-10  109.94 4.55e-06 ***
#   Time:Space          20 5.680e-10 2.840e-11   30.83 0.000184 ***
#   Time:Type           14 1.695e-09 1.211e-10  131.48 2.84e-06 ***
#   Space:Type           4 1.120e-10 2.790e-11   30.33 0.000404 ***
#   Time:Location       37 2.472e-09 6.680e-11   72.57 1.33e-05 ***
#   Space:Location      23 9.830e-10 4.270e-11   46.42 5.40e-05 ***
#   Time:Space:Type     28 1.074e-09 3.840e-11   41.66 7.14e-05 ***
#   Time:Space:Location 68 1.261e-09 1.860e-11   20.14 0.000532 ***
#   Residuals            6 6.000e-12 9.000e-13                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#This model say they are all signifcant and all important. So if everything is interacting then everything is important?


model<-aov(Slope~Time+Space+Type+Range+Row, data = N20emission)
summary(model)
# Df    Sum Sq   Mean Sq F value  Pr(>F)    
# Time          7 1.906e-08 2.723e-09  66.655 < 2e-16 ***
# Space         3 2.680e-10 8.930e-11   2.186 0.09085 .  
# Type          2 1.890e-10 9.450e-11   2.312 0.10169    
# Location     17 1.721e-09 1.012e-10   2.478 0.00146 ** 
#   Residuals   200 8.171e-09 4.090e-11                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Type seems to be important in some models but not in other
model<-aov(Slope~Time*Space*Type, data = N20emission)
summary(model)

model<-aov(Flux.hr~Time*Space*Type, data = N20emission)
summary(model)
#Double check to make sure the head space volumes are correct for week 3
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Time              7 0.9829 0.14042  38.787  < 2e-16 ***
#   Space             3 0.0135 0.00449   1.240 0.297238    
# Type              2 0.0030 0.00150   0.415 0.661132    
# Time:Space       20 0.0556 0.00278   0.768 0.748683    
# Time:Type        14 0.1520 0.01086   3.000 0.000438 ***
#   Space:Type        4 0.0079 0.00197   0.545 0.702994    
# Time:Space:Type  28 0.0747 0.00267   0.737 0.826647    
# Residuals       151 0.5467 0.00362                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

model<-lmer(Slope~Time*Space*Type+(1|Range)+(1|Row), data = N20emission)
summary(model)
#Lets look at the varaince across the different treatments====

library(data.table)
setDT(N20emission)[, list(GroupVariance=var(rep(Slope, Flux.min))), keyby = Space]
# Space GroupVariance
# 1: Bulk 1  3.012640e-10
# 2: Bulk 2  1.949342e-10
# 3:      C  1.228240e-10
# 4:      F  1.450779e-10
# 5:      M  1.141068e-10
setDT(N20emission)[, list(GroupVariance=var(rep(Slope, Flux.min))), keyby = Type]
#Apperently hybrids have a bunch of varaince. 
# Type GroupVariance
# 1:      B73  1.123864e-10
# 2:     Bulk  5.561329e-10
# 3:    Bulk   1.895747e-10
# 4:   Hybrid  7.469320e-11
# 5: Teosinte  1.945172e-10

#It may be useful to look at the combined factor between the two datasets. 



#Lets remove all of the different time factors.
Week1=subset(N20emission, Time== "Week1")

boxplot(Week1$N2O.Concentation~ Week1$Type)
boxplot(Week1$Slope~ Week1$Type)
boxplot(Week1$Flux.min~ Week1$Type)
boxplot(Week1$Flux.hr~ Week1$Type)
boxplot(Week1$Headspace.Concentration..ppm.~ Week1$Type)

Week2=subset(N20emission, Time== "Week2")

boxplot(Week2$N2O.Concentation~ Week2$Type)
boxplot(Week2$Slope~ Week2$Type)
boxplot(Week2$Flux.min~ Week2$Type)
boxplot(Week2$Flux.hr~ Week2$Type)
boxplot(Week2$Headspace.Concentration..ppm.~ Week2$Type)


Week3=subset(N20emission, Time== "Week3")

boxplot(Week3$N2O.Concentation~ Week3$Type)
boxplot(Week3$Slope~ Week3$Type)
boxplot(Week3$Flux.min~ Week3$Type)
boxplot(Week3$Flux.hr~ Week3$Type)
boxplot(Week3$Headspace.Concentration..ppm.~ Week3$Type)

Week5=subset(N20emission, Time== "Week5")

boxplot(Week5$N2O.Concentation~ Week5$Type)
boxplot(Week5$Slope~ Week5$Type)
boxplot(Week5$Flux.min~ Week5$Type)
boxplot(Week5$Flux.hr~ Week5$Type)
boxplot(Week5$Headspace.Concentration..ppm.~ Week5$Type)

Week5=subset(N20emission, Time== "Week5")

boxplot(Week5$N2O.Concentation~ Week5$Type)
boxplot(Week5$Slope~ Week5$Type)
boxplot(Week5$Flux.min~ Week5$Type)
boxplot(Week5$Flux.hr~ Week5$Type)
boxplot(Week5$Headspace.Concentration..ppm.~ Week5$Type)


Week6=subset(N20emission, Time== "Week6")

boxplot(Week6$N2O.Concentation~ Week6$Type)
boxplot(Week6$Slope~ Week6$Type)
boxplot(Week6$Flux.min~ Week6$Type)
boxplot(Week6$Flux.hr~ Week6$Type)
boxplot(Week6$Headspace.Concentration..ppm.~ Week6$Type)

#The relationship amoong the different genotypes seems to be chnaging. 
#Perhaps I will attempt to look into bootstrapping the data
# Bootstrap 95% CI for R-Squared
library(boot)

# function to obtain R-Squared from the data
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
}
# bootstrapping with 1000 replications
results <- boot(data=mtcars, statistic=rsq,
                R=10000, formula=mpg~wt+disp)

results <- boot(data=NILNitf, statistic=rsq,
                R=1000, formula=Average.rate.ug.NO2.N.g.hr~Genotypes)
# view results
results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")

##I'm not super sure what bootstrapping is doing



###NIL 2020 GreenHouse analysis=====

library(igraph)
library(vegan)
library(ggplot2)
library(forcats)
library(plyr)

setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL follow- Exp 2019/Data /")
#New Folder Set -- 10/26/23
setwd("~/Documents/Thesis /Chapter 4/Science /Rstats/NIL follow- Exp 2019/Data /")

  NILNitf = read.csv(file ="NIL2020.Nit.csv", head = T)
  
  Means <- ddply(NILNitf, "Genotypes", summarise, Nit.mean=mean(Average.rate.ug.NO2.N.g.hr))
  Varaince <- ddply(NILNitf, "Genotypes", summarise, Nit.mean=var(Average.rate.ug.NO2.N.g.hr))
  Means <- ddply(NILNitf, "Genotypes", summarise, Nit.mean=mean(Ratio.BM))
  Means <- ddply(NILNitf, "Genotypes", summarise, Nit.mean=mean(Ratio.BM.2))
  
  #NILNitf=subset(NILNitf, Rep!="Rep5")
  
  #  NILNitf=subset(NILNitf, Sample.ID!="NIL2-161")
  #  NILNitf=subset(NILNitf, Sample.ID!="NIL2-008")
  # 
  NILNitf=subset(NILNitf, Genotypes!= "Bulk Soil Mix")
  
  
  NILNitf=subset(NILNitf, Genotypes!= "Z030E1023")
  NILNitf=subset(NILNitf, Genotypes!= "Z038E0074")
  
  Chromosome5=subset(NILNitf, Chromosome=="N chr 5"| Chromosome=="Control")
  Chromosome9=subset(NILNitf, Chromosome=="N chr 9"| Chromosome=="Control")
  
  ggplot(NILNitf, aes(x =Average.rate.ug.NO2.N.g.hr)) + geom_density(alpha = 0.5)
  ggplot(NILNitf, aes(x =diff.rate.ng.NO2.N.g.hr)) + geom_density(alpha = 0.5)
  
#Genotypic view || This first methods was the one used by Mark and Sandy
  ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Average.rate.ug.NO2.N.g.hr, .fun = median, .desc=TRUE), y=Average.rate.ug.NO2.N.g.hr, fill=Genotypes))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 229, linetype="dashed", color = "red")
  
  ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Average.rate.ug.NO2.N.g.hr, .fun = median, .desc=TRUE), y=Average.rate.ug.NO2.N.g.hr/10, fill=Genotypes))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 229/100, linetype="dashed", color = "red")

  #This is the orginal method
  ggplot(NILNitf, aes(x = fct_reorder(Genotypes, diff.rate.ng.NO2.N.g.hr, .fun = median, .desc=TRUE), y=diff.rate.ng.NO2.N.g.hr, fill=Genotypes))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 229, linetype="dashed", color = "red")
  
  par(mfrow=c(1,2))
  
  #Genotype does not seem to be signficant--- This may be caused by high levels fo chance in the dataset
  model<-aov((Average.rate.ug.NO2.N.g.hr)~Genotypes+Plant.BM+Error(Rep), data = NILNitf)
  
  model<-aov((diff.rate.ng.NO2.N.g.hr)~Genotypes+Plant.BM+Error(Rep), data = NILNitf)
  #TukeyHSD(model)
  #Very close to signicance|| If one removes the highly varaibel samples we will see signficanvce in the samples
  summary(model)
  ##P value is slightly smaller when you take the averges earelier 
  
#Chromosomeal view of the data
  ggplot(NILNitf, aes(x = fct_reorder(Chromosome, Average.rate.ug.NO2.N.g.hr, .fun = mean, .desc=TRUE), y=Average.rate.ug.NO2.N.g.hr, fill=Chromosome))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_hline(yintercept=260, linetype="dashed", color = "red")
  
  model<-aov(((Average.rate.ug.NO2.N.g.hr))~Chromosome+Error(Rep), data = NILNitf)
  summary(model)
  TukeyHSD(model)

#Weird stuff happened across the replicates  
  ggplot(NILNitf, aes(x = fct_reorder(Rep, Average.rate.ug.NO2.N.g.hr, .fun = mean, .desc=TRUE), y=Average.rate.ug.NO2.N.g.hr, fill=Rep))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_point()
  
  #Rep seems to be a signficiant factor in the model 
  model<-aov((Average.rate.ug.NO2.N.g.hr)~Rep, data = NILNitf)
  TukeyHSD(model)
  summary(model)
  
  
#This is the un-averaged data - From the two runs
  #Run 1
  ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Rep.1.C.adj, .fun = median, .desc=TRUE), y=Rep.1.C.adj, fill=Genotypes))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  model<-aov((Rep.1.C.adj)~Genotypes+Error(Rep), data = NILNitf)
  #TukeyHSD(model)
  summary(model)
  
  #Run 2
  ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Rep.2.C.adj, .fun = median, .desc=TRUE), y=Rep.2.C.adj, fill=Genotypes))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  model<-aov((Rep.2.C.adj)~Genotypes+Error(Rep), data = NILNitf)
  TukeyHSD(model)
  summary(model)
#The means of the runs are significantly different and this should be noted in the analysis process. 
# Although most of the patterns are the same across the runs. 
  
  
  ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Abs.N.Rep1, .fun = median, .desc=TRUE), y=Abs.N.Rep1, fill=Genotypes))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Abs.N.Rep.2, .fun = median, .desc=TRUE), y=Abs.N.Rep.2, fill=Genotypes))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
   ggplot(NILNitf, aes(x = fct_reorder(Genotypes, ((Rep.1.C.adj+Rep.2.C.adj)/2), .fun = median, .desc=TRUE), y=((Rep.1.C.adj+Rep.2.C.adj)/2), fill=Genotypes))+
    geom_boxplot()+labs(y="ngN/gh", title="Potential Nitrification Rate", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
   
   model<-aov((Abs.N.Rep1)~Genotypes+Error(Rep), data = NILNitf)
   TukeyHSD(model)
   summary(model)
   
   
   Testcomapr=subset(NILNitf, Genotypes=="Z031E0021"| Genotypes=="B73")
   
   t.test(Average.rate.ug.NO2.N.g.hr~Genotypes, data=Testcomapr)
   
   
   
   ###Lets look at these from the chromosome perspective
  
  #Chromosome 5 Anlaysis  
   ggplot(Chromosome5, aes(x = fct_reorder(Genotypes, Average.rate.ug.NO2.N.g.hr, .fun = median, .desc=TRUE), y=Average.rate.ug.NO2.N.g.hr, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="Chr5 PNR -Median", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     geom_hline(yintercept = 229, linetype="dashed", color = "red")
   
   model<-aov((Average.rate.ug.NO2.N.g.hr)~Genotypes+Error(Rep), data = Chromosome5)
   summary(model)
   plot(summary(model))
   
   
   
   
   #Chromosome 9 Analysis 
   ggplot(Chromosome9, aes(x = fct_reorder(Genotypes, Average.rate.ug.NO2.N.g.hr, .fun = median, .desc=TRUE), y=Average.rate.ug.NO2.N.g.hr/10, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="Chr9 PNR -Median", x="Genotypes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     geom_hline(yintercept = 22.9, linetype="dashed", color = "red")
   
   model<-aov((Average.rate.ug.NO2.N.g.hr)~Genotypes+Error(Rep), data = Chromosome9)
   summary(model)
   
  #Correcting for plant biomass is in the pot setting really helps the results and conclusions drawn.
   #THe major issue with the green house expirements are that they seem be confounded with plant size|
   #Since this nitrfication spression is a function of plant size and time I want to make a correction for growth .
#Plant Biomass
   #Genotypic view 
   ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Plant.BM, .fun = median, .desc=TRUE), y=Plant.BM, fill=Genotypes))+
     geom_boxplot()+labs(y="grams", title="Aboveground Plant Biomass", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #Here I attmpt to contorol for biomass
   ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Ratio.BM, .fun = median, .desc=TRUE), y=Ratio.BM, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="PNR/PlantBiomass", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

   ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Ratio.BM.2, .fun = median, .desc=TRUE), y=Ratio.BM, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="PNR/PlantBiomass", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
   
   ##Chromosome plots=== Manuscprit figure ## Add a line for the mediaum =====
   NIL2-035
   Chr5.outlier=subset(Chromosome5, Sample.ID!="NIL2-035")
   
   
   ggplot(Chr5.outlier, aes(x = fct_reorder(Genotypes, Ratio.BM.2, .fun = median, .desc=TRUE), y=(Ratio.BM/10), fill=Genotypes))+
     geom_boxplot()+labs(y="Potential Nitrification Rate/Plant Biomass ((ngN/gh)/g)", title="PNR/PlantBiomass", x="Genotypes")+
     theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     geom_hline(yintercept = 4.125, linetype="dashed", color = "red")
   
   #Noted: Used Rep.1.Rates, excluded rep2. Large replicate effect in these data. excluding an outlier sample NIL2-035
   ggplot(Chr5.outlier, aes(x = fct_reorder(Genotypes, Ratio.BM.2, .fun = median, .desc=TRUE), y=(Rep.1.C.adj/Plant.BM)/10, fill=Genotypes))+
     geom_boxplot()+labs(y="Potential Nitrification Rate/Plant Biomass ((ngN/gh)/g)", title="Chromosome 5 NILs", x="Genotypes")+
     theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     geom_hline(yintercept = 7.35, linetype="dashed", color = "red")
   
   T.test.Ratio.BM=subset(Chr5.outlier, Genotypes=="Z031E0021"| Genotypes=="B73")
   t.test((Rep.1.C.adj/Plant.BM)~Genotypes, data=T.test.Ratio.BM)
   t.test(Ratio.BM.2~Genotypes, data=T.test.Ratio.BM)
   
   # data:  (Rep.1.C.adj/Plant.BM) by Genotypes
   # t = 3.404, df = 2.8182, p-value = 0.04652
   # T.test.Ratio.BM=subset(Chromosome5, Genotypes=="Z031E0021"| Genotypes=="B73")
   # t.test(Ratio.BM.2~Genotypes, data=T.test.Ratio.BM)
   # 
   
   model<-aov(log(Rep.1.C.adj/Plant.BM)~Genotypes+Error(Rep), data = Chromosome5)
   summary(model)
              
              
   ggplot(T.test.Ratio.BM, aes(x = fct_reorder(Genotypes, Ratio.BM.2, .fun = median, .desc=TRUE), y=(Rep.1.C.adj/Plant.BM)/10, fill=Genotypes))+
     geom_boxplot()+labs(y="Potential Nitrification Rate/Plant Biomass ((ngN/gh)/g)", title="Chromosome 5 NILs", x="Genotypes")+
     theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     geom_hline(yintercept = 7.35, linetype="dashed", color = "red")
   
   
   T.test.Ratio.BM=subset(Chr5.outlier, Genotypes=="Z031E0021"| Genotypes=="B73")
   t.test((Rep.1.C.adj/Plant.BM)~Genotypes, data=T.test.Ratio.BM)
   
   
   T.test.Ratio.BM=subset(Chr5.outlier, Genotypes=="Z035E1047"| Genotypes=="B73")
   t.test((Rep.1.C.adj/Plant.BM)~Genotypes, data=T.test.Ratio.BM)
   
   T.test.Ratio.BM=subset(Chr5.outlier, Genotypes=="Z034E1047"| Genotypes=="B73")
   t.test((Rep.1.C.adj/Plant.BM)~Genotypes, data=T.test.Ratio.BM)

  
   
   
   ggplot(Chromosome5, aes(x = fct_reorder(Genotypes, Rep.1.C.adj, .fun = median, .desc=TRUE), y=(Rep.1.C.adj/100), fill=Genotypes))+
     geom_boxplot()+labs(y="Potential Nitrification Rate/Plant Biomass ((ngN/gh)/g)", title="PNR/PlantBiomass", x="Genotypes")+
     theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     geom_hline(yintercept = 4.125, linetype="dashed", color = "red")
   
   
   ###----Chromosome 9 Figure-----
   Chr9.outlier=subset(Chromosome9, Sample.ID!="NIL2-161"& Sample.ID!="NIL2-055")
   
   NIL2-161
   
   ggplot(Chr9.outlier, aes(x = fct_reorder(Genotypes, Ratio.BM.2, .fun = median, .desc=TRUE), y=(Rep.1.C.adj/Plant.BM)/10, fill=Genotypes))+
     geom_boxplot()+labs(y="Potential Nitrification Rate/Plant Biomass ((ngN/gh)/g)", title="Chromosome 9 NILs ", x="Genotypes")+
     theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     geom_hline(yintercept = 7.35, linetype="dashed", color = "red")
   
   
   T.test.Ratio.BM=subset(Chr9.outlier, Genotypes=="Z031E0047"| Genotypes=="B73")
   t.test((Rep.1.C.adj/Plant.BM)~Genotypes, data=T.test.Ratio.BM)
   # data:  (Rep.1.C.adj/Plant.BM) by Genotypes
   # t = 2.4625, df = 3.394, p-value = 0.0808
   
   #The broader trending in means are more important for this confirmation 
   
   Chromosome5
   
   T.test.Ratio.BM=subset(Chromosome9, Genotypes=="Z031E0047"| Genotypes=="B73")
   t.test((Rep.1.C.adj/Plant.BM)~Genotypes, data=T.test.Ratio.BM)
   
   
   
   
   
   ggplot(Chromosome9, aes(x = fct_reorder(Genotypes, Ratio.BM, .fun = median, .desc=TRUE), y=Ratio.BM/10, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="PNR/PlantBiomass", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
   
   ggplot(Chromosome5, aes(x = fct_reorder(Genotypes, Ratio.BM.2, .fun = median, .desc=TRUE), y=Ratio.BM.2, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="PNR/PlantBiomass", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
   ggplot(Chromosome9, aes(x = fct_reorder(Genotypes, Ratio.BM.2, .fun = median, .desc=TRUE), y=Ratio.BM.2, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="PNR/PlantBiomass", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
   
   
   
   #Nitrfication/biomass is stastically significant
   model<-aov((Ratio.BM.2)~Genotypes+Error(Rep), data = NILNitf)
   TukeyHSD(model)
   summary(model)
   
   #Plant biomass is signficantly different across treatment 
   model<-aov((Plant.BM)~Genotypes+Error(Rep), data = NILNitf)
   summary(model)
   
   
   #Very close to signicance|| If one removes the highly varaibel samples we will see signficanvce in the samples
   summary(model)
   #Slight Negative corrlation with Biomass|| i.e. nitrfication increases, biomass decreases
   cor(NILNitf$Average.rate.ug.NO2.N.g.hr,NILNitf$Plant.BM)  
   #Very, very, very slight negative correlation
   ggplot(NILNitf,aes(x=Average.rate.ug.NO2.N.g.hr, y=Plant.BM))+geom_point()+stat_smooth(method=lm)
   
   
   ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Rep.1.C.adj/Plant.BM, .fun = median, .desc=TRUE), y=Rep.1.C.adj/Plant.BM, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="PNR/PlantBiomass", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
   ggplot(NILNitf, aes(x = fct_reorder(Genotypes, Rep.2.C.adj/Plant.BM, .fun = median, .desc=TRUE), y=Rep.2.C.adj/Plant.BM, fill=Genotypes))+
     geom_boxplot()+labs(y="ngN/gh", title="PNR/PlantBiomass", x="Genotypes")+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
   
   
#Varaince | need to figure it out |
   var(NILNitf$Average.rate.ug.NO2.N.g.hr)
   var(NILNitf$Rep.1.C.adj)
   var(NILNitf$Rep.2.C.adj)
   var(NILNitf$Plant.BM)
   
   
   
   
 