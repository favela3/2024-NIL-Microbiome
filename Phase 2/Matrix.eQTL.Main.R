#Alonso Favela
#Matrix eQTL-Toy example

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');
# base.dir = '.';

#Code below is to test all gene pairs
######
## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");

# Output file name
output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 3000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 3000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)



#####
library(MatrixEQTL)

# This option is so it can deal with the sample load 
options(max.print=1000000)
#Now I'm going to try something with my samples
######
## Settings


setwd("~/Documents/Thesis /Chapter 4/Science /Rstats/Old R-Analysis /R-Scripts/Matrix_eQTL")

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

#Here I need to edit the things. #Here I need to use dplyer to remove teostine from the files
Snp<-read.table("snps.txt")
Gene<-read.table("NILOTU.txt")
Map<-read.table("NILMap.txt")
Genes.Trans<-read_csv("~/Documents/Thesis /Chapter 4/Science /Rstats/ASV.Tansformed.csv")
Genes.ASV<-as.table(Genes.Trans[,-97])

Gene<-read.table("Full.ASV.Genes.txt")

##This is important when reformatting file=====
#write.table(Genes.ASV, "ASV.Genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

##NIL.T1.97 is missing from the datasets

#Specifcally I need to remove teosinte and hybrids.
Snp2 <- Snp[,-1]
rownames(Snp2) <-  Snp[,1]

Snp2 <- Snp2[-1,]
colnames(Snp2) <- Snp[1,]

snpfree<-Snp2[rowSums(Snp2[])>0,]

# Genotype file name
SNP_file_name = paste( "snps.txt", sep="");

# Gene expression file name
#expression_file_name = paste("NILOTU.Ordered.txt", sep="");
expression_file_name = paste("Full.ASV.Genes.txt", sep="")

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste("NILMap.txt", sep="");

# Output file name
output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-6;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 500;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  #This is an important part of code
  pvOutputThreshold = 1e-9,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)
pvalue.hist = "qqplot"

#Method for exporting the Data 
qtls<-(as.data.frame(me$all$eqtls))
#write.csv(qtls, file = "SignificantQTLs.csv")

### I want to plot the density
d<-density(qtls$pvalue)
plot(d)

d<-density(qtls$statistic)
plot(d)

d<-density(qtls$beta)

plot(d, main="Density vs Beta")
#polygon(d, col="Black", border="Black")

hist(qtls$beta, col="blue")

###Looking at the QTL for teosinte-Maize Genetic map =====
#These were the results I made for the presentation
# PastData<-read.csv("SignificantQTLs.csv")
# PastSortedCounts<-PastData%>%count(snps, sort=TRUE)
# PastLocANDCount<-G.Location%>%full_join(PastSortedCounts,by = c("snp"="snps"))
# PastLocANDCount[,6][is.na(PastLocANDCount[,6])]<-0
# ggplot(PastLocANDCount, aes((chr+(dist.cM./100)), n))+geom_bar(position="dodge", stat="identity", color="black")+
#   theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
#   ggtitle("Past Results- Blue Map")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")
# 


#This old run had about 5000 or so association
#With the filtered ASV we see about 800 or so correlations out of 6000ASV

library(tidyverse)
#qtls%>%group_by(snps)%>%tally()

SortedCounts<-qtls%>%count(snps, sort=TRUE)
#I think it would be usful to add a key to the loci information

G.Location<-read.csv("GenomicLocation.csv")

#full_join(G.Location, qtls, by = c("snps"="snp"))

LocANDCount<-G.Location%>%full_join(SortedCounts,by = c("snp"="snps"))

#This removes the Na files and turns them into zeros
LocANDCount[,6][is.na(LocANDCount[,6])]<-0

#Genetic Map Graphic====
#This is pretty cool. I've been looking at the distribution 
#This view broadly allows us to get an idea. 
ggplot(LocANDCount, aes((chr+dist.cM.), n))+geom_bar(position="dodge", stat="identity")
#+
ggplot(LocANDCount, aes((chr+(dist.cM./100)), n))+geom_bar(position="dodge", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 w teosinte")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")

# I want to set a robust threshold analysis on the micorbes being altered. Currently I'm attempting 1SD away from the mean

ggplot(LocANDCount, aes((chr+(dist.cM./100)), n))+geom_bar(position="dodge", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-4 w teosinte")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")


#write.table(LocANDCount, file = "LocusOTU-Count-WithTeosinte.csv", sep = ",",
            #row.names = FALSE,quote=FALSE,col.names = TRUE)
            
#How many of our QTL are above 1SD of the mean?
#This is the code I wanted... This will give you your QTLs with an association above the Mean+1SD====
KeyLoci<-LocANDCount%>%filter(n > mean(n,na.rm=TRUE)+sd(n))

KeyLociMean<-LocANDCount%>%filter(n > mean(n,na.rm=TRUE))

#The list of regions that I should be interesed in looking at
write.table(KeyLoci, file = "MAG-1e-06-M+SD.csv", sep = ",",
            row.names = TRUE,quote=FALSE,col.names = TRUE)

# LocANDCount[,3]<-as.numeric(LocANDCount[,3])
# write.table(LocANDCount, file = "MAG-1e-06.csv", sep = ",",
#             row.names = FALSE,quote=FALSE,col.names = TRUE)

#IMPORTANT Threshold filter was at 1e-10=====
#Threshhold filter was at 10


#These results went tdown to 831. There seems to be a clear genetic relationship. and th
#4000 or so were discovered when removing all other NILs. Its weired tha that the NIls 
#Statistic is the effect size of the interaction. when the attional 

library(ggplot2)
ggplot(qtls, aes(x = qtls$beta)) + geom_density(alpha = 0.5)

# Perform the same analysis recording information for 
# a Q-Q plot
meq = Matrix_eQTL_engine(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt, 
  output_file_name = output_file_name,
  pvOutputThreshold = 1e-4, 
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = TRUE,
  pvalue.hist = "qqplot");
unlink( output_file_name );
# png(filename = "QQplot.png", width = 650, height = 650)
plot(meq, pch = 16, cex = 0.7)
 #dev.off();


#####Removing the Teosinte and Hybrids from the data using tidyr----
#These files are labeled as New Gene, Map and Snp

library(tidyr)
library(tidyverse)
setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL Maize Microbiome /R-Scripts/Matrix_eQTL")


Snp<-read.table("snps.txt")
Gene<-read.table("NILOTU.txt")
Map<-read.table("NILMap.txt")

#Need to remove these two file types and re-run
#PI384071, B73xPI384071, "B73-maybe:

#Transpose then grep out the needed information
tMap<-data.frame(t(Map))

#This grep command grabed both the teosinte lines and the 
grep("PI384071",tMap$X4)
#[1]  10  35  47  56 119 120 139 160 185 210 222 231 294 295 314 335

#This code will give you all of the needed infromation and all other covarates
#Sample ID are important as they are the ones I need to remove from the other files
tMap[grep("PI384071",tMap$X4),]
#By placing the second number in the file we will the the sample ID list that we will need to remove from the data frames
tMap[grep("PI384071",tMap$X4),1]

# X1 X2 X3           X4 X5
# V10  NIL.T1.009  2 41     PI384071 T1
# V35  NIL.T1.034  3 49 B73xPI384071 T1
# V47  NIL.T1.046  4 13 B73xPI384071 T1
# V56  NIL.T1.055  4 49     PI384071 T1
# V119 NIL.T1.119  7 61 B73xPI384071 T1
# V120 NIL.T1.120  7 57     PI384071 T1
# V139 NIL.T1.139  8 33     PI384071 T1
# V160 NIL.T1.160  9 73 B73xPI384071 T1
# V185 NIL.T2.009  2 41     PI384071 T2
# V210 NIL.T2.034  3 49 B73xPI384071 T2
# V222 NIL.T2.046  4 13 B73xPI384071 T2
# V231 NIL.T2.055  4 49     PI384071 T2
# V294 NIL.T2.119  7 61 B73xPI384071 T2
# V295 NIL.T2.120  7 57     PI384071 T2
# V314 NIL.T2.139  8 33     PI384071 T2
# V335 NIL.T2.160  9 73 B73xPI384071 T2

#Now I need to remove the colums with these sample ID present in them! Then we can run our analysis 

#This code from the tidyverse does the same thing. Here the carrot makes the search more specifc. Regardless of method what is important is that you obtain a list of sample IDs to remove
tMap %>% filter(str_detect(X4, '^PI384071'))
tMap %>% filter(str_detect(X4, 'PI384071'))
# X1 X2 X3           X4 X5
# 1  NIL.T1.009  2 41     PI384071 T1
# 2  NIL.T1.034  3 49 B73xPI384071 T1
# 3  NIL.T1.046  4 13 B73xPI384071 T1
# 4  NIL.T1.055  4 49     PI384071 T1
# 5  NIL.T1.119  7 61 B73xPI384071 T1
# 6  NIL.T1.120  7 57     PI384071 T1
# 7  NIL.T1.139  8 33     PI384071 T1
# 8  NIL.T1.160  9 73 B73xPI384071 T1
# 9  NIL.T2.009  2 41     PI384071 T2
# 10 NIL.T2.034  3 49 B73xPI384071 T2
# 11 NIL.T2.046  4 13 B73xPI384071 T2
# 12 NIL.T2.055  4 49     PI384071 T2
# 13 NIL.T2.119  7 61 B73xPI384071 T2
# 14 NIL.T2.120  7 57     PI384071 T2
# 15 NIL.T2.139  8 33     PI384071 T2
# 16 NIL.T2.160  9 73 B73xPI384071 T2

#! Mark will mean to remove from file. This bit of code for example will remove all with PI tMap %>% filter(!str_detect(X4, 'PI384071'))
#THis combination of code also works which is pretty cool
#data.frame(t(Map))%>%filter(str_detect(X4, '^PI384071'))

#This seemingly worked. We removed #16 files from 351 and note have 336. THis math checks out
NewMap<-data.frame(t(Map))%>%filter(!str_detect(X4, 'PI384071'))%>%
  rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

#lets try this with another set of data
#Lets try this on snp data first
Removelist<-tMap[grep("PI384071",tMap$X4),1]
#Example of Filtering from List: This is how you look for muliple lists=====
NewSnp<-data.frame(t(Snp))%>%filter(!str_detect(X1, paste(Removelist, collapse=("|"))))%>%
  rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

#Now I want to do this to the gene tabel 

NewGene<-data.frame(t(Gene))%>%filter(!str_detect(X1, paste(Removelist, collapse=("|"))))%>%
  rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

#Removing the variable key the tidyverse give the infromation to keep track of everyhting.
NewGene<-NewGene[,-1]
NewSnp<-NewSnp[,-1]
NewMap<-NewMap[,-1]


#Writing table 
# write.table(NewGene, file = "NewGene.txt", sep = "\t",
#             row.names = FALSE,quote=FALSE,col.names = FALSE)
# write.table(NewSnp, file = "NewSnp.txt", sep = "\t",
#             row.names = FALSE,quote=FALSE,col.names = FALSE)
# write.table(NewMap, file = "NewMap.txt", sep = "\t",
#             row.names = FALSE,quote=FALSE,col.names = FALSE)

##Rerunning Mapping results without Teosinte, and hybrid NILs-----
##Alonso Start Here 2/18/20======
setwd("~/Documents/Projects/Research/Maize Microbiome/NIL MM2018/R-Scripts/Matrix_eQTL")

library(MatrixEQTL)

# Genotype file name
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS


##IVE RUN THIS WITH BOTH ALL AND T2
SNP_file_name = paste( "NewSnp.txt", sep="");

# Gene expression file name
expression_file_name = paste("NewGene.txt", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste("Map.txt", sep="");

useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Output file name
output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-6;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 50;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)
pvalue.hist = "qqplot"

#Method for exporting the Data 
qtls<-(as.data.frame(me$all$eqtls))
#write.csv(qtls, file = "SignificantQTLs-NoTeoHy-T2.csv")


### I want to plot the density
d<-density(qtls$pvalue)
d<-density(qtls$statistic)
d<-density(qtls$beta)

plot(d, main="Density vs Beta",xlim=c(0, 5))
polygon(d, col="Black", border="Black")

hist(qtls$beta, col="blue",xlim=c(0, 50))

#The threshold effects the amount of mapped QTLs in the analysis =====
#For example changing the threshold from pvOutputThreshold, 1e-6 to 1e-8 will cause a huge reduction in the determined loci
#About 20,000 to 600... this about a 33% fold decrease. What I wasn't expecting was the the remove of the Teosinte, and hybrid caused the lowered finding of NILs.
#There are 22306 matches. Seems like the additon of the teosinte results in the removal of false positives. 
#I'm a little shocked by these results! 

#This is a super weird distrbution.
ggplot(qtls, aes(x = qtls$pvalue)) + geom_density(alpha = 0.5)
ggplot(qtls, aes(x = qtls$statistic)) + geom_density(alpha = 0.5)

ggplot(qtls, aes(x = qtls$b)) + 
  geom_density(alpha = 0.5, fill="black")+ 
  geom_vline(aes(xintercept=mean(qtls$b)),color="red", linetype="dashed", size=1)


#Where are the QTLs found to be important=====
###Now I want to count the number of quniquw loci
#This is a method to cound the number of files 

#qtls%>%group_by(snps)%>%tally()

SortedCounts<-qtls%>%count(snps, sort=TRUE)
#I think it would be usful to add a key to the loci information

G.Location<-read.csv("GenomicLocation.csv")


full_join(G.Location, qtls, by = c("snps"="snp"))

LocANDCount<-G.Location%>%full_join(SortedCounts,by = c("snp"="snps"))

LocANDCount[,6][is.na(LocANDCount[,6])]<-0

#This is pretty cool. I've been looking at the distribution 
#This view broadly allows us to get an idea. 
ggplot(LocANDCount, aes((chr+dist.cM.), n))+geom_bar(position="dodge", stat="identity")
#
ggplot(LocANDCount, aes((chr+(dist.cM./210)), n))+geom_bar(position="dodge", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 No TeoHy")+xlab("Chromosome+cM") + ylab("Number of OTU")+geom_rug(sides = "b")


ggplot(LocANDCount, aes((chr+(dist.cM./203)), n))+geom_bar(position="dodge", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 No TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")

ggplot(LocANDCount, aes(AGP.v1.position, n))+geom_bar(position="dodge", stat="identity", color="black")+theme_classic()+ggtitle("Micorbe-Assocation across Genome | 1e-6 No TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")




ggplot(LocANDCount, aes((chr+(dist.cM./100)), n))+geom_bar(position="dodge", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")

#What does this look like on a chromosomal level| These are stacked
ggplot(LocANDCount, aes(chr, n))+geom_bar(position="stack", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")


# write.table(LocANDCount, file = "MAG-1e-06-LocCout-NOTEO-T2.csv", sep = ",",
# row.names = TRUE,quote=FALSE,col.names = TRUE)

KeyLoci<-LocANDCount%>%filter(n > mean(n,na.rm=TRUE)+sd(n))
KeyLoci2SD<-LocANDCount%>%filter(n > mean(n,na.rm=TRUE)+sd(n)*2)

#Distribution afeter filtering 
ggplot(KeyLoci, aes(chr, n))+geom_bar(position="stack", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")


#Analysis with two standard di
ggplot(KeyLoci, aes((chr+(dist.cM./203)), n))+geom_bar(position="dodge", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("MAG|ChrView|2SD| 1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")


ggplot(KeyLoci2SD, aes(chr, n))+geom_bar(position="stack", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("MAG|ChrView| 2SD| 1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")


#KeyLociMean<-LocANDCount%>%filter(n > mean(n,na.rm=TRUE))

#The list of regions that I should be interesed in looking at
# write.table(KeyLoci, file = "MAG-1e-06-MSD-NOTEO-T2.csv", sep = ",",
# row.names = TRUE,quote=FALSE,col.names = TRUE)
#About 22 percent of the maize genome maps

###lets get a clearier idea of what kind of taxa are being mapped to these regions.
OTU.Taxa<-read.csv("OTU.Taxa.csv")
QTLs<-as.data.frame(qtls)


QTLs<-qtls[,c(2,1,3,4,5,6)]
OTU.Taxa[,1]<-as.numeric(OTU.Taxa[,1])
TaxaSplit<-OTU.Taxa%>%separate(col=taxonomy,into=c("k","p","c","o","f","g",'s'),sep=";")

QTLs[,1]<-as.numeric(QTLs[,1])

QTOTUs<-TaxaSplit%>%full_join(QTLs, by = c("X.OTU.ID"="gene"))
QTOTUs<-OTU.Taxa%>%full_join(QTLs, by = c("X.OTU.ID"="gene"))
QTOTUs<-na.omit(QTOTUs)

#QTOTUs%>%filter(str_detect(snps))%>%count(k, sort=TRUE)
#Alows to subdivide

ASV.taxa<-read.csv("ASV.taxa.csv")
QTLs<-as.data.frame(qtls)
QTOTUs<-ASV.taxa%>%full_join(QTLs, by = c("Sample.ID"="gene"))
QTOTUs<- QTOTUs%>%filter(beta!="NA")


#Be aware that running other packages could disrupt this line of code. 
CountbyTaxa<-QTOTUs%>%group_by(snps)%>%count(p, sort=TRUE)
#Allows full taxa
CountbyTaxa<-QTOTUs%>%group_by(snps)%>%count(taxonomy, sort=TRUE)

Full<-CountbyTaxa%>%full_join(KeyLoci, by = c("snps"="snp"))
Full<-na.omit(Full)

Full<-CountbyTaxa%>%full_join(LocANDCount, by = c("snps"="snp"))
Full<-na.omit(Full)

###This code is for writng out all of the infromation
 Full[,6]<-as.numeric(unlist(Full[,6]))
 
 # write.table(Full, file = "DataToMakeStackPlotGenome-NEW.3.15.csv", sep = ",",
 #             row.names = FALSE,quote=FALSE,col.names = TRUE)

#I'm trying to make stacked plots at the different regions. I need to make the colors deeper.
#I'm not sure why its not doing this... D:
ggplot(Full, aes((chr+(dist.cM./203)), y=n.x, fill=p))+geom_bar(position="stack", stat="identity",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
  theme(legend.position = "none")


ggplot(Full, aes((chr+(dist.cM./203)), y=(n.x/2), fill=p))+geom_bar(position="stack", stat="identity",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+
  theme(legend.position = "none")

ggplot(Full, aes((chr+(dist.cM./203)), y=(n.y/n.x), fill=p))+geom_bar(position="stack", stat="identity",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+
  theme(legend.position = "none")

ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=taxonomy))+geom_bar(position="stack", stat="identity",width = 0.05)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | Stacked |1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+theme(legend.position = "none")

##Esa figure
ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=p))+geom_bar(position="stack", stat="identity",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
  theme(legend.position = "none")


#looking at assocations across the chromomsomes
ggplot(Full, aes((chr), y=n.x, fill=p))+geom_bar(position="stack", stat="identity")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | Stacked |1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+theme(legend.position = "none")



#WithBoxoutline
ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=p))+geom_bar(position="stack", stat="identity",width = 0.05, order=Full$n.x)+
     theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
     ggtitle("Micorbe-Assocation across Genome | Stacked|  1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+
     geom_rug(sides = "b")

#I'm pretty happy with this work flow and it gives very interesting and generally applealing figues. 



multiplot(Mapwith,Mapwith2)
# ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=taxonomy))+geom_bar(position="stack", stat="identity",width = 0.01)+
#   theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
#   ggtitle("Micorbe-Assocation across Genome | 1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
#   #geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")
  
  
  
  
ggplot(Full, aes((chr+(dist.cM./100)), y=n.y, color="black"))+geom_bar(position="stack", stat="identity")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")

#Doing this type of analysis keys are results tight and to a point. It allows us to see exactly what the genetic regions we are interested in are doing.


###Lets look at the hypo mixture
# m <- 10; n <- 7; k <- 8
# x <- 0:(k+1)
# rbind(phyper(x, m, n, k), dhyper(x, m, n, k))
# all(phyper(x, m, n, k) == cumsum(dhyper(x, m, n, k)))  # FALSE
# ## but error is very small:
# signif(phyper(x, m, n, k) - cumsum(dhyper(x, m, n, k)), digits = 3)


#These are the genomic regions at 1e-06 and filtered for the mean+SD: This has teosinte in it too
# chr1:34723974..46254743
# chr4: 20064775..177976304
# chr5:7086156..145120913
# chr7:168990870..169017086
# chr8: 23796614..161460056
# chr9:110522266..148712750
# chr10: 16264665..72565617


###Lets see if we can make figure using the abunance of those micorbial groips for the samples in question
##Reworking microbial community data into analysis=====
#I basically want to filter the data to see if we can find alterations in abundance across our key OTUs for Just our key otus
#Full file should act as my key for analysis
#NOOOOOTTTTEEE: I kinda lost track of what I wanted to do here. I can now pull out the OTUs from the analysis, 
# the NILs for a region if I want to do a tighter comparison
OTU<-read.table("NewGene.txt",header = TRUE)
MAP<-read.table("NewMap.txt",header = TRUE)
SNP<-read.table("NewSnp.txt",header = TRUE)

#I will use this to get information about exact OTUs| Currently this script pull our the Loci, region,etc and tells us whose there
#This bit of code pulls out the List of OTUs assocaited with a region
#Chr5OTU<-QTOTUs%>%filter(str_detect(snps, '^PZA00148-3'))%>%select(1)

#Grep works pretty well for this tooo. Basically I'm telling it look for x, in snps and keep column 1
Removelist<-QTOTUs[grep("PZA00148",QTOTUs$snps),1]

#Here I added a character for the paste code to work
Removelist<-prepend(Removelist, "X")

#The paste part is simply manipulating the list to add brakes
Chr5OTU<-data.frame(OTU)%>%filter(str_detect(id, paste(Removelist, collapse=("$|^"))))

#This almost works!!, just need to figure out how to do this before my list


Test<-SNP%>%filter(str_detect(id, '^PZA00148-3'))%>%slice(1)%>%
  rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

NILReg<-Test[grep("1",Test$"1"),1]

# names(Test) <-  unlist(Test[1,])
# Test<-Test[-1,]
# names(Test)
# 
# Test[grep("1",Test$"PZA00148-3"),1]
# 
# Test[grep("1",Test$"1"),1]


#This is how you look for muliple lists 


  # #This seemingly worked. We removed #16 files from 351 and note have 336. THis math checks out
  # NewMap<-data.frame(t(Map))%>%filter(!str_detect(X4, 'PI384071'))%>%
  # rownames_to_column %>%
  # gather(variable, value, -rowname) %>% 
  # spread(rowname, value)





##Let run the FastExQTL analysis with the ASV across the Genetic map across the entire dataset===



# Genotype file name
SNP_file_name = paste( "snps.txt", sep="");

# Gene expression file name
#expression_file_name = paste("NILOTU.Ordered.txt", sep="");
expression_file_name = paste("Full.ASV.Genes.txt", sep="")

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste("NILMap.txt", sep="");

# Output file name
output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-6;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 500;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  #This is an important part of code
  pvOutputThreshold = 1e-9,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)
pvalue.hist = "qqplot"

#Method for exporting the Data 
qtls<-(as.data.frame(me$all$eqtls))
#write.csv(qtls, file = "SignificantQTLs.csv")

### I want to plot the density
d<-density(qtls$pvalue)
plot(d)

d<-density(qtls$statistic)
plot(d)

d<-density(qtls$beta)

plot(d, main="Density vs Beta")
polygon(d, col="Black", border="Black")

ggplot()
hist(qtls$beta, col="blue")
##Here we are showing the small impacts, but some with very large


betafilter<-qtls%>%filter(beta>2)
Supmean<-qtls%>%filter(beta>mean(qtls$beta))
Submean<-qtls%>%filter(beta<mean(qtls$beta))

mean(qtls$beta)
median(qtls$beta)
 
plot(density(betafilter$beta), main="Density vs Beta")

ggplot(qtls, aes(x = qtls$b)) + 
  geom_density(alpha = 0.5, fill="black")+ 
  geom_vline(aes(xintercept=mean(qtls$b)),color="red", linetype="dashed", size=1)+ggtitle("Beta Distribution")+labs(x="Beta", y="Density")



library(tidyverse)
#qtls%>%group_by(snps)%>%tally()
#qtls<-QTOTU
##Counting the Number ASV per SNP
SortedCounts<-qtls%>%count(snps, sort=TRUE)

#SortedCounts<-Supmean%>%count(snps, sort=TRUE)

#I think it would be usful to add a key to the loci information

G.Location<-read.csv("GenomicLocation.RemovedOverlap.csv")

#full_join(G.Location, qtls, by = c("snps"="snp"))

#Binding SNP Location to 
LocANDCount<-G.Location%>%full_join(SortedCounts,by = c("snp"="snps"))

#This removes the Na files and turns them into zeros
LocANDCount[,6][is.na(LocANDCount[,6])]<-0

Loci<-LocANDCount%>%filter(n!=0)
mean(Loci$n)
AboveMean<-Loci%>%filter(beta<10)

#Here we are doing a compariosn 
281/504
##0.57%
##Greater then the
100/504
#0.20

##Some how the mitrocondrial DNA made it through to this stage Probabily best to just remove them. 
##WE see that there is still am impact of the MC

#Genetic Map Graphic====
#This is pretty cool. I've been looking at the distribution 
#This view broadly allows us to get an idea. 
ggplot(LocANDCount, aes((chr+dist.cM.), n))+geom_bar(position="dodge", stat="identity")
#+
ggplot(LocANDCount, aes((chr+(dist.cM./100)), n))+geom_bar(position="dodge", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 w teosinte")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")


##Adding Taxaonmoy

###Trying to replicate Figure with ASV=====

ASV.taxa<-read.csv("ASV.taxa.csv")
QTLs<-as.data.frame(qtls)

#Binding Taxanomy to OTU
QTOTUs<-ASV.taxa%>%full_join(QTLs, by = c("Sample.ID"="gene"))
QTOTUs<- QTOTUs%>%filter(beta!="NA")

QTOTU<-QTOTUs%>%filter(Family!="Mitochondria")
##This counting scripts are really important for the manuscript=======
#Be aware that running other packages could disrupt this line of code.
##This returns the number of time that ASV of said phylum has appeared at that SNP====
CountbyTaxa<-QTOTUs%>%group_by(snps)%>%count(Phylum, sort=TRUE)

##Here we can see the location of the gene, its correlation, and the SNP introgression loation-
#highlighting the linkage
CountbyTaxaLoca<-G.Location%>%full_join(QTOTUs,by = c("snp"="snps"))

CountbyTaxaLoca<-CountbyTaxaLoca%>%filter(beta!="NA")%>%filter(Family!="Mitochondria")
#Allows full taxa
#CountbyTaxa<-QTOTUs%>%group_by(snps)%>%count(taxonomy, sort=TRUE)

# Full<-CountbyTaxa%>%full_join(KeyLoci, by = c("snps"="snp"))
# Full<-na.omit(Full)

##Here we are joining the number of phlyum by the total number of taxa per snp
Full<-CountbyTaxa%>%full_join(LocANDCount, by = c("snps"="snp"))
Full<-na.omit(Full)

FullCountbyTaxaLoca<-G.Location%>%full_join(Full,by = c("snp"="snps"))

###This code is for writng out all of the infromation
#Full[,6]<-as.numeric(unlist(Full[,6]))

# write.table(Full, file = "DataToMakeStackPlotGenome-NEW.3.15.csv", sep = ",",
#             row.names = FALSE,quote=FALSE,col.names = TRUE)

Full%>%group_by(chr)%>%count(dist.cM., sort=TRUE)
#The reason we are seeing distinct SNPs in those locations is that  there are two SNPS in the same chr

#I'm trying to make stacked plots at the different regions. I need to make the colors deeper.
#I'm not sure why its not doing this... D:
ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-9 | w Teo")+xlab("Chromosome") + ylab("Number of ASV")+
  geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+2*sd(LocANDCount[,6])), linetype="dashed", color="red")+
theme(legend.position = "none")

  

  ImportantSNPS<-Full%>%filter(n.y>mean(LocANDCount[,6]+sd(LocANDCount[,6])))%>%count(snps)

  ImportantSNPS<-G.Location%>%full_join(ImportantSNPS,by = c("snp"="snps"))
  ImportantSNPS<-na.omit(ImportantSNPS)
  

  Nosure<-CountbyTaxaLoca%>%count(Sample.ID)

ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="dodge", stat="identity",width = .1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | w Teo")+xlab("Chromosome") + ylab("Number of ASV")+
  geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")#+
  theme(legend.position = "none")

# ggplot(Full, aes((chr+(dist.cM./100)), y=n.y, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
#   theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
#   ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | w Teo")+xlab("Chromosome") + ylab("Number of ASV")+
#   geom_rug(sides = "b")+
#   geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
#   theme(legend.position = "none")


ggplot(LocANDCount, aes((chr+(dist.cM./100)), n))+geom_bar(position="dodge", stat="identity", color="black",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 w teosinte")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")





ggplot(Full, aes((chr+(dist.cM./203)), y=(n.x/2), fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+
  theme(legend.position = "none")

##Count by Phylum/ total number of counts
ggplot(Full, aes((chr+(dist.cM./100)), y=(n.y/n.x), fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")#+
  theme(legend.position = "none")

ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.05)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | Stacked |1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")#+theme(legend.position = "none")

##Esa figure
ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
  theme(legend.position = "none")


#looking at assocations across the chromomsomes
ggplot(Full, aes((chr), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | Stacked |1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")+theme(legend.position = "none")



#WithBoxoutline
ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.05)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | Stacked|  1e-6 NO TeoHy")+xlab("Chromosome") + ylab("Number of OTU")+
  geom_rug(sides = "b")



####Lets remove Teosinte======
setwd("~/Documents/Thesis /Chapter 4/Science /Rstats/Old R-Analysis /R-Scripts/Matrix_eQTL")

library(tidyr)
library(tidyverse)


Snp<-read.table("snps.txt")
Gene<-read.table("Full.ASV.Genes.txt")
Map<-read.table("NILMap.txt")

#Need to remove these two file types and re-run
#PI384071, B73xPI384071, "B73-maybe:

#Transpose then grep out the needed information
tMap<-data.frame(t(Map))

#This grep command grabed both the teosinte lines and the 
grep("PI384071",tMap$X4)
#[1]  10  35  47  56 119 120 139 160 185 210 222 231 294 295 314 335

#This code will give you all of the needed infromation and all other covarates
#Sample ID are important as they are the ones I need to remove from the other files
tMap[grep("PI384071",tMap$X4),]
#By placing the second number in the file we will the the sample ID list that we will need to remove from the data frames
tMap[grep("PI384071",tMap$X4),1]



#Now I need to remove the colums with these sample ID present in them! Then we can run our analysis 

#This code from the tidyverse does the same thing. Here the carrot makes the search more specifc. Regardless of method what is important is that you obtain a list of sample IDs to remove
tMap %>% filter(str_detect(X4, '^PI384071'))
tMap %>% filter(str_detect(X4, 'PI384071'))

#! Mark will mean to remove from file. This bit of code for example will remove all with PI tMap %>% filter(!str_detect(X4, 'PI384071'))
#THis combination of code also works which is pretty cool
#data.frame(t(Map))%>%filter(str_detect(X4, '^PI384071'))

#This seemingly worked. We removed #16 files from 351 and note have 336. THis math checks out
NewMap<-data.frame(t(Map))%>%filter(!str_detect(X4, 'PI384071'))%>%
  rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

#lets try this with another set of data
#Lets try this on snp data first
Removelist<-tMap[grep("PI384071",tMap$X4),1]
#Example of Filtering from List: This is how you look for muliple lists=====
NewSnp<-data.frame(t(Snp))%>%filter(!str_detect(X1, paste(Removelist, collapse=("|"))))%>%
  rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

#Now I want to do this to the gene tabel 

NewGene<-data.frame(t(Gene))%>%filter(!str_detect(X1, paste(Removelist, collapse=("|"))))%>%
  rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

#Removing the variable key the tidyverse give the infromation to keep track of everyhting.
NewGene<-NewGene[,-1]
NewSnp<-NewSnp[,-1]
NewMap<-NewMap[,-1]


#Writing table 
# write.table(NewGene, file = "ASV.NIL.Gene.txt", sep = "\t",
#             row.names = FALSE,quote=FALSE,col.names = FALSE)
# 
# write.table(NewSnp, file = "ASV.NIL.Snp.txt", sep = "\t",
#             row.names = FALSE,quote=FALSE,col.names = FALSE)
# write.table(NewMap, file = "ASV.NIL.Map.txt", sep = "\t",
#             row.names = FALSE,quote=FALSE,col.names = FALSE)


##Rerun without Teosinte and Hybrids ASV=====
library(MatrixEQTL)

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste( "ASV.NIL.Snp.txt", sep="");

# Gene expression file name
expression_file_name = paste("ASV.NIL.Gene.txt", sep="")

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste("ASV.NIL.Map.txt", sep="");

# Output file name
output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-6;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 500;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  #This is an important part of code
  pvOutputThreshold = 1e-10,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name);

##Im gping to use the without teosint ewith a threshold of 1e-10
## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');

show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)
pvalue.hist = "qqplot"

#Method for exporting the Data 
qtls<-(as.data.frame(me$all$eqtls))
#write.csv(qtls, file = "SignificantQTLs.csv")

### I want to plot the density
d<-density(qtls$pvalue)
plot(d)

d<-density(qtls$statistic)
plot(d)

d<-density(qtls$beta)

plot(d, main="Density vs Beta")
#polygon(d, col="Black", border="Black")

hist(qtls$beta, col="blue")



#Adding Taxaonmoy

###Trying to replicate Figure with ASV=====

ASV.taxa<-read.csv("ASV.taxa.csv")
QTLs<-as.data.frame(qtls)
G.Location<-read.csv("GenomicLocation.RemovedOverlap.csv")

#Binding Taxanomy to OTU
QTOTUs<-ASV.taxa%>%full_join(QTLs, by = c("Sample.ID"="gene"))
QTOTUs<- QTOTUs%>%filter(beta!="NA")

QTOTUs<- QTOTUs%>%filter(Family!="Mitochondria")

SortedCounts<-QTOTUs%>%count(snps, sort=TRUE)

##We see 376 Loci having a correlation
##Out of 506

##Counting the Number ASV per SNP
#SortedCounts<-qtls%>%count(snps, sort=TRUE)
#I think it would be usful to add a key to the loci information

#G.Location<-read.csv("GenomicLocation.RemovedOverlap.csv")

#full_join(G.Location, qtls, by = c("snps"="snp"))

#Binding SNP Location to 
LocANDCount<-G.Location%>%full_join(SortedCounts,by = c("snp"="snps"))

LocANDCount<-na.omit(LocANDCount)


#Be aware that running other packages could disrupt this line of code. 
CountbyTaxa<-QTOTUs%>%group_by(snps)%>%count(Phylum, sort=TRUE)

#CountbyTaxa<-QTOTUs%>%group_by(Sample.ID)%>%count(Phylum, sort=TRUE)


Full<-CountbyTaxa%>%full_join(LocANDCount, by = c("snps"="snp"))

Full<-na.omit(Full)

###This code is for writng out all of the infromation
#Full[,6]<-as.numeric(unlist(Full[,6]))

##This counting scripts are really important for the manuscript=======
#Be aware that running other packages could disrupt this line of code.
##This returns the number of time that ASV of said phylum has appeared at that SNP====
CountbyTaxa<-QTOTUs%>%group_by(snps)%>%count(Phylum, sort=TRUE)

##Here we can see the location of the gene, its correlation, and the SNP introgression loation-
#highlighting the linkage
CountbyTaxaLoca<-G.Location%>%full_join(QTOTUs,by = c("snp"="snps"))


#FullCountbyTaxaLoca<-G.Location%>%full_join(Full,by = c("snp"="snps"))

###This code is for writng out all of the infromation
#Full[,6]<-as.numeric(unlist(Full[,6]))

# write.table(Full, file = "DataToMakeStackPlotGenome-NEW.3.15.csv", sep = ",",
#             row.names = FALSE,quote=FALSE,col.names = TRUE)

Full%>%group_by(chr)%>%count(dist.cM., sort=TRUE)
#The reason we are seeing distinct SNPs in those locations is that  there are two SNPS in the same chr


ggplot(LocANDCount, aes((chr+(dist.cM./100)), n))+geom_bar(position="dodge", stat="identity", color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-6 w teosinte")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")
# 
# ggplot(FullCountbyTaxaLoca, aes((chr+(dist.cM./100)), beta))+geom_bar(position="dodge", stat="identity", color="black")+
#   theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
#   ggtitle("Micorbe-Assocation across Genome | 1e-6 w teosinte")+xlab("Chromosome") + ylab("Number of OTU")+geom_rug(sides = "b")+
#   geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")

New.FUll<-G.Location%>%full_join(QTOTUs,by = c("snp"="snps"))

##Large Effect QTLs and there Taxa being Shapped
New.FUll%>%ggplot( aes((chr+(dist.cM./100)), beta, fill=Class))+geom_bar(position="dodge", stat="identity",width =1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome | 1e-10 w teosinte")+xlab("Chromosome") + 
  ylab("Beta (effect size)")+geom_rug(sides = "b")+theme(legend.position = "none")

ggplot(New.FUll, aes((chr+(dist.cM./180)), y=Order, fill=Class)) + 
  geom_dotplot(binaxis='y', stackdir='center')+scale_x_continuous(breaks=seq(1,10,by=1))+
  theme(legend.position = "none")


#I'm trying to make stacked plots at the different regions. I need to make the colors deeper.
#I'm not sure why its not doing this... D:
##Chaning up the 
ggplot(Full, aes((chr+(dist.cM./180)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.05, color="black")+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Regions Shaping Large Number of OTUs")+xlab("Chromosome") + ylab("Number of ASV")+
  geom_rug(sides = "b")+
  #geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
  theme(legend.position = "none")

##This is going to be the manuscript figure. ##Interesintly inters a sharp dropoff

ggplot(Full, aes((chr+(dist.cM./200)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity", width=0.05)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-10 | w/o Teo")+xlab("Chromosome") + ylab("Number of ASV")+
  geom_rug(sides = "b")#+
  geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
  theme(legend.position = "none")


ggplot(Full, aes((chr), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 1)+
  theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
  ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | w Teo")+xlab("Chromosome") + ylab("Number of ASV")+
  geom_rug(sides = "b")#+
  #geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
  theme(legend.position = "none")



# ggplot(Full, aes((chr+(dist.cM./100)), y=n.x/n.y, fill=Phylum))+geom_bar(position="dodge", stat="identity",width = 0.1)+
#   theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
#   ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
#   geom_rug(sides = "b")#+
#   geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")
#   theme(legend.position = "none")

  ggplot(Full, aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="dodge", stat="identity",width = 1)+
    theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
    ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
    geom_rug(sides = "b")+theme(legend.position = "none")
    geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
  theme(legend.position = "none")
    
    
    Full%>%filter(Phylum=="Actinobacteriota")%>%ggplot( aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
      theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
      ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
      geom_rug(sides = "b")+theme(legend.position = "none")
    geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
      theme(legend.position = "none")
    
    Full%>%filter(Phylum=="Proteobacteria")%>%ggplot( aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
      theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
      ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
      geom_rug(sides = "b")+theme(legend.position = "none")
    geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
      theme(legend.position = "none")
    
    Full%>%filter(Phylum=="Chloroflexi")%>%ggplot( aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
      theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
      ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
      geom_rug(sides = "b")+theme(legend.position = "none")
    geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
      theme(legend.position = "none")
    
    Full%>%filter(Phylum=="Proteobacteria"|Phylum=="Actinobacteriota")%>%ggplot( aes((chr+(dist.cM./100)), y=n.x, fill=Phylum))+geom_bar(position="stack", stat="identity",width = 0.1)+
      theme_classic()+scale_x_continuous(breaks=seq(1,10,by=1))+
      ggtitle("Micorbe-Assocation across Genome |Stacked| 1e-6 | NOTeo")+xlab("Chromosome") + ylab("Number of OTU")+
      geom_rug(sides = "b")#+theme(legend.position = "none")
    geom_hline(yintercept=mean(LocANDCount[,6]+sd(LocANDCount[,6])), linetype="dashed", color="red")+
      theme(legend.position = "none")
    
    ##Here we will do some number collecting 
    
    ImportantSNPS<-Full%>%filter(n.y>mean(LocANDCount[,6]+sd(LocANDCount[,6])))%>%count(snps)
    LocANDCount%>%filter(n>500)
    
    ImportantSNPS<-G.Location%>%full_join(ImportantSNPS,by = c("snp"="snps"))
    ImportantSNPS<-na.omit(ImportantSNPS)
    
    ##504 total number of SNPs
    375/504
    
    
    Nosure<-CountbyTaxaLoca%>%count(Sample.ID)
    
    betafilter<-QTOTUs%>%filter(beta>11)
    mean(QTOTUs$beta)+2*sd(QTOTUs$beta)
    median(QTOTUs$beta)
    Supmean<-QTOTUs%>%filter(beta>mean(QTOTUs$beta)+3*sd(QTOTUs$beta))
    Submean<-QTOTUs%>%filter(beta<mean(QTOTUs$beta)+3*sd(QTOTUs$beta))
    
    Supmean<-QTOTUs%>%filter(beta>mean(QTOTUs$beta))
    Submean<-QTOTUs%>%filter(beta<mean(QTOTUs$beta))
    
    SortedCounts<-betafilter%>%count(snps, sort=TRUE)
    LocANDCount<-G.Location%>%full_join(SortedCounts,by = c("snp"="snps"))
    
    LocANDCount<-na.omit(LocANDCount)
    
    mean(qtls$beta)
    median(qtls$beta)
    
    
    sd(Full$n.x)
    
    ggplot(QTOTUs, aes(x = beta)) + 
      geom_density(alpha = 0.5, fill="black")+ 
      geom_vline(aes(xintercept=mean(beta)),color="red", linetype="dashed", size=1)+
      ggtitle("Beta Distribution| 1e-10 w/oTeo")+labs(x="Beta (effect size)", y="Density")+xlim(0,100)
    
    
    ggplot(QTOTUs, aes(x = beta)) +geom_histogram(bins = 500)+scale_x_binned()
    
    
    ##Pulling out the important NILs
    LocANDCount%>%filter(n>20)%>%filter(chr==5)#%>%pull(snp)
      
    greaterthen100<-LocANDCount%>%filter(n>100)#%>%filter(chr==5)#%>%pull(snp)
    
    write_csv(greaterthen100, "GreaterThan100ASValtered.csv")
    
    ##This is really really coool we can 
    LocANDCount%>%filter(n>90)%>%pull(snp)%>%as.data.frame()
      NewSnp%>%filter(str_detect(V1, 'PZB00054-3'))
    
    NewSnp%>%filter(V1=="PZA00112-5"| V1=="id")%>%t()%>%as.data.frame()%>%filter(V2!=0)%>%pull(V1)#%>%str_replace_all("/s", ",")
    
    as.data.frame(t(NewMap))%>%filter(str_detect(V1, "NIL.T2.091"))
    
    Full
    
    NewSnp$V1
    
    
    Full%>%filter()
    
    

