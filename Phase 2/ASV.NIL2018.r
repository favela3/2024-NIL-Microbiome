##Alonso Favela
##Decemeber 2nd 2022
##NIL Re-Analysis for publication| Basically, we want to determine how The microbial communites are chaning in time


##DADA2 OUTPUT ANALYSIS =====

library(ggplot2)
library(tidyverse)
library(phyloseq)
library(vegan)

setwd("~/Documents/Thesis /Chapter 4/Science /Rstats")

seqtab<- readRDS("NIL2018.ASV/seqtab_final.RDS")
summaryDada<- readRDS("NIL2018.ASV/summary-of-dada2.RDS")
taxa<- readRDS("NIL2018.ASV/taxonomy_final.RDS")
samdf<-read_csv("NIL2018.ASV/NIL2018MappingDNA.csv")


##Cleaning up Data
rownames(seqtab)<-rownames(seqtab)%>%str_replace("V4_515F_New_V4_806R_New-", "")%>%str_extract(".*(?=\\_)")

samdf<-as.data.frame(samdf)
rownames(samdf)<-samdf$`Sample ID`


##Need to make metadata for this seed stuff.
##The sample data thing is here
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps@sam_data

taxa<-data.frame(tax_table(ps))

##Rename barcode as ASV
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

taxaASV<-(ps@tax_table)


#We see a clear partitioning of the seeds spray microbiome
plot_richness(ps, measures=c("Shannon", "Observed"), x="Time")

cleanps<-subset_taxa(ps, Order!="Chloroplast")
cleanps<-subset_taxa(cleanps, Family!="Mitochondria")

##Saving the Phyloseq object
#saveRDS(ps, file = "NILcleanASV.OTU")

cleanps<-subset_samples(cleanps, Sample.ID!="NIL_T2_089")
cleanps<-subset_samples(cleanps, Sample.ID!="NIL_T2_069")
cleanps<-subset_samples(cleanps, Sample.ID!="NIL_T2_090")
cleanps<-subset_samples(cleanps, Sample.ID!="NIL_T2_086")
cleanps<-subset_samples(cleanps, Sample.ID!="NIL_T1_010")

CleanNIL<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="PI384071"|Pedigree=="Z031E0047"|Pedigree=="Z031E0021"|Pedigree=="Z031E0523"|Pedigree=="Z031E0591")
CleanNIL<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="Z031E0047"|Pedigree=="Z031E0021"|Pedigree=="Z031E0523"|Pedigree=="Z031E0591")

##Here we can pull out the ASV_Table 
asv_table <- rownames_to_column(data.frame(otu_table(cleanps)), var = "Sample.ID")


GP.ord <- ordinate(CleanNIL, "NMDS", "bray",)
p1 = plot_ordination(cleanps, GP.ord, color="Time",  title="NIL Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Sample.ID), size = 5, vjust = 1.5)
print(p1)

CleanNIL%>%subset_samples(Time!="T1")%>%plot_ordination( GP.ord, color="Pedigree",  title="NIL Ordination")+
  geom_point(size=3)+
  #geom_text(mapping = aes(label = Sample.ID), size = 5, vjust = 1.5)+
  stat_ellipse()

print(p1)

metadata <- as(sample_data(cleanps), "data.frame")

adonis2( phyloseq::distance(cleanps, method="bray") ~ Pedigree+Time+Rep+Range+Row+Block,
         data = metadata)
#                Df SumOfSqs      R2      F Pr(>F)    
#   Pedigree       44    8.173 0.13211 1.1197  0.001 ***
#   Time            1    0.843 0.01362 5.0785  0.001 ***
#   Rep             1    0.137 0.00222 0.8286  0.928    
#   Range          10    3.620 0.05852 2.1822  0.001 ***
#   Row            21    4.415 0.07137 1.2673  0.001 ***
#   Pedigree:Time  43    6.689 0.10811 0.9376  0.997    
#   Residual      229   37.992 0.61407                  
#   Total         349   61.870 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##This is the interaction model and is consistent with our other results! This is still pretty exciting! 
##Strated block 
# adonis2(formula = phyloseq::distance(cleanps, method = "bray") ~ Pedigree * Time + Rep + Block, data = metadata, strata = metadata$Block)
# Df SumOfSqs      R2      F Pr(>F)    
# Pedigree       44    8.173 0.13211 1.1791  0.001 ***
#   Time            1    0.843 0.01362 5.3481  0.001 ***
#   Rep             1    0.137 0.00222 0.8726  0.793    
# Block         135   26.317 0.42536 1.2374  0.001 ***
#   Pedigree:Time  43    6.707 0.10841 0.9901  0.728    
# Residual      125   19.693 0.31829                  
# Total         349   61.870 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Lets remove the bulk 
NILsV4<-subset_samples(cleanps, Pedigree!="BULK")

##Filtering OTUS greater then 0.005%
FSr  = transform_sample_counts(NILsV4, function(x) x / sum(x) )
FSfr = filter_taxa(FSr, function(x) sum(x) > .005, TRUE)

NILsV4.05<-FSfr

##I'm going to clean up the data by removing the EXTREMELY low abundance outs from the data
adonis2( phyloseq::distance(NILsV4.05, method="bray") ~ Pedigree*Time+Rep+Block,
         data = metadata, strata = metadata$Block)

# adonis2(formula = phyloseq::distance(NILsV4.05, method = "bray") ~ Pedigree * Time + Rep + Block, data = metadata, strata = metadata$Block)
#                  Df SumOfSqs      R2      F Pr(>F)    
#   Pedigree       44    7.211 0.13328 1.2138  0.003 ** 
#   Time            1    0.832 0.01538 6.1625  0.001 ***
#   Rep             1    0.114 0.00211 0.8448  0.813    
#   Block         135   23.306 0.43074 1.2786  0.003 ** 
#   Pedigree:Time  43    5.766 0.10656 0.9931  0.650    
#   Residual      125   16.878 0.31193                  
#   Total         349   54.107 1.00000                  
# ---


Nitri<-subset_taxa(cleanps, Order=="Nitrospirales"|Family=="Nitrosomonadaceae")

Nitri<-subset_taxa(cleanps, Order=="Nitrososphaerales")#| Order=="Nitrospirales"|Family=="Nitrosomonadaceae")

Nitri<-subset_samples(Nitri, Pedigree=="B73"| Pedigree=="PI384071"| Pedigree=="B73xPI384071"|Pedigree=="Z031E0047"|Pedigree=="Z031E0021")

##This is pretty intersting for time2
plot_bar(Nitri, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")+
  facet_wrap(~Time,scales="free_x")

##Actinobacteria 
actino<-subset_taxa(cleanps, Phylum=="Actinobacteriota")
##This is pretty intersting for time2
plot_bar(actino, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")+
  facet_wrap(~Time,scales="free_x")

actinoT2<-subset_samples(actino, Time=="T2")
actinoT1<-subset_samples(actino, Time=="T1")

##Interestinyl we see the patter in arche
plot_bar(actinoT2, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")#+
  facet_wrap(~Time,scales="free_x")

  plot_bar(actinoT1, fill="Phylum",x="Pedigree")+
    geom_bar(stat="identity")#+
  facet_wrap(~Time,scales="free_x")
###Lets so a little DESEQ2
library(DESeq2)
TeoMaiz<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="PI384071")

T1<-subset_samples(TeoMaiz, Time=="T1")

diagdds = phyloseq_to_deseq2(T1, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="mean") ## This step isnt working for some reason,


##This is really, really annoying|| Dont know what to do. 
res = results(diagdds, cooksCutoff = FALSE)
alpha = .05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Subset)[rownames(sigtab), ], "matrix"))
head(sigtab)
rownames(sigtab)
dim(sigtab)

Teosinte.B73.T2.DESEQ<-sigtab
#saveRDS(Teosinte.B73.T2.DESEQ, file = "RObjects/Teosinte.B73.T2.DESEQ.rds")

#T1 [1] 109   6 OTUS
#T2 [1] 142   6

TvM<-subset_taxa(T1, rownames(tax_table(T1)) %in% c(rownames(sigtab)))

plot_bar(TvM, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")#+
#  facet_wrap(~Concentration,scales="free_x")
##This is pretty great as we previously didnt find any difference between the two

Deseqactinos<-subset_taxa(TvM, Phylum=="Actinobacteriota")
plot_bar(Deseqactinos, fill="Family",x="Pedigree")+
  geom_bar(stat="identity")

#Solirubrobacteracease


TeoMaiz<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="Z031E0021")

T1<-subset_samples(TeoMaiz, Time=="T2")

diagdds = phyloseq_to_deseq2(T1, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="mean") ## This step isnt working for some reason,
#Its because of all the zeros in the sata

##This is really, really annoying|| Dont know what to do. 
res = results(diagdds, cooksCutoff = FALSE)
alpha = .05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Subset)[rownames(sigtab), ], "matrix"))
head(sigtab)
rownames(sigtab)
dim(sigtab)

NIL21.B73.T2.DESEQ<-sigtab
#saveRDS(NIL21.B73.T2.DESEQ, file = "RObjects/NIL21.B73.T2.DESEQ.rds")


TvM<-subset_taxa(T1, rownames(tax_table(T1)) %in% c(rownames(sigtab)))

plot_bar(TvM, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")#+


TeoMaiz<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="Z031E0047")

T1<-subset_samples(TeoMaiz, Time=="T2")

diagdds = phyloseq_to_deseq2(T1, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="mean") ## This step isnt working for some reason,
#Its because of all the zeros in the sata
#Comparing control to the AZ Seeds

##This is really, really annoying|| Dont know what to do. 
res = results(diagdds, cooksCutoff = FALSE)
alpha = .05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Subset)[rownames(sigtab), ], "matrix"))
head(sigtab)
rownames(sigtab)
dim(sigtab)

NIL47.B73.T2.DESEQ<-sigtab
#saveRDS(NIL47.B73.T2.DESEQ, file = "RObjects/NIL47.B73.T2.DESEQ.rds")


TvM<-subset_taxa(T1, rownames(tax_table(T1)) %in% c(rownames(sigtab)))

plot_bar(TvM, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")#+
##We see alot more actinobacteria, Bdellovibrionota, and Armatinonadota

##Lets compare these lists: =====
##This is awesome, theres acutally OTUs present across these that are shared. 
NIL21.Inter<-intersect(rownames(NIL21.B73.T2.DESEQ),rownames(Teosinte.B73.T2.DESEQ))


DESEQCompar<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="Z031E0021"|Pedigree=="PI384071")

DESEQCompar<-subset_samples(DESEQCompar, Time=="T2")

TvM<-subset_taxa(DESEQCompar, rownames(tax_table(DESEQCompar)) %in% c(intersect(rownames(NIL21.B73.T2.DESEQ),rownames(Teosinte.B73.T2.DESEQ))))

plot_bar(TvM, fill="Order",x="Pedigree")+
  geom_bar(stat="identity")#+



NIL47.Inter<-intersect(rownames(NIL47.B73.T2.DESEQ),rownames(Teosinte.B73.T2.DESEQ)) 

DESEQCompar<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="Z031E0047"|Pedigree=="PI384071")

DESEQCompar<-subset_samples(DESEQCompar, Time=="T2")

TvM<-subset_taxa(DESEQCompar, rownames(tax_table(DESEQCompar)) %in% c(intersect(rownames(NIL47.B73.T2.DESEQ),rownames(Teosinte.B73.T2.DESEQ))))

plot_bar(TvM, fill="Order",x="Pedigree")+
  geom_bar(stat="identity")#+



intersect(NIL47.Inter, NIL21.Inter)

DESEQCompar<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="Z031E0047"|Pedigree=="PI384071"| Pedigree=="Z031E0021")

DESEQCompar<-subset_samples(DESEQCompar, Time=="T2")

TvM<-subset_taxa(DESEQCompar, rownames(tax_table(DESEQCompar)) %in% c(intersect(NIL47.Inter, NIL21.Inter)))

plot_bar(TvM, fill="Order",x="Pedigree")+
  geom_bar(stat="identity")#+


##When comparing the list it seems like there are taxa that are both ptrsent in the rhizosphere 

##These results are making alot more sense
##ASV results are the same. 
#Overall the H2 Analysis was useful to determing the taxa that are important
##Lets try to calcualte heretibility =========

##Caluclate the Most Hereitble OTUs to use as a filter 

# heritability analysis


library("tidyverse")
library("phyloseq")
library("lme4")
library("ggpubr")
library("ggrepel")


# calculate mean microbe abundance for each maize genotype
asv_table <- rownames_to_column(data.frame(otu_table(ps)), var = "Sample.ID")
NILsV4.05.1<-NILsV4.05
NILsV4.05<-subset_samples(NILsV4.05, Time=="T2")

### filter asv table for 2019 data without checks
asv_table <- rownames_to_column(data.frame(otu_table(NILsV4.05)), var = "Sample.ID")
counts <- left_join(data.frame(sample_data(NILsV4.05)), asv_table)


### calculate means of ASV counts across all subsamples by group
mean_counts <- counts %>%
  group_by(Pedigree, Time, Block) %>%
  summarize_each(funs(mean), -Sample.ID, -Range, -Row)


### calculate log relative abundance
asvtab <- mean_counts[, 5:ncol(mean_counts)]
logrel <- t(apply(asvtab, 1, function(x) log(x/sum(x) + 0.001)))
mean_counts_logrel <- cbind(data.frame(mean_counts[, 1:4]), logrel)




### split stdN, lowN samples, retain only genotypes with at least 2 reps (both blocks)
# (note that genotypes 38-11, A214N, and CI90C were planted twice and have 4 reps)

h2dat_stdN <- mean_counts_logrel %>%
  filter(Pedigree != "Bulk") %>%
  group_by(Pedigree) %>%
  add_tally(name="count") %>%
  filter(count >= 2) %>%
  dplyr::select(-count)

## calculate H2 for +N
df <- h2dat_stdN
traits <- colnames(df)[5:ncol(df)]

H2_stdN <- data.frame()

for(trait in traits){
  #trait <- "asv_000013"
  f <- formula(paste0(trait, ' ~ (1|Pedigree) + (1|Block)'))
  fit <- lmer(f, data=df)
  v <- as.data.frame(VarCorr(fit))
  Vg <- v$vcov[v$grp == "Pedigree"]
  Ve <- v$vcov[v$grp == "Residual"]
  H2 <- round(Vg/(Vg + Ve/6), 6)
  H2_stdN <- rbind(H2_stdN, data.frame("ASV"=trait, "H2_stdN_19"=H2))
}

H2_Genotype<-H2_stdN

#write_csv(H2_Genotype, file = "H2_ASV.csv")

H2Microbes<-H2_Genotype%>%filter(H2_stdN_19 > 0.3)

##This may be useful
hist(H2_Genotype$H2_stdN_19, breaks=20)
hist(H2Microbes$H2_stdN_19)

##Can use this to filter the dataset 
##I think then I will focus on those highly heretible OTUs i.e. OTUs that change with plant genotype
##HEre we are running an ANOVA on the heretible otus 
HeretibleASV<-subset_taxa(cleanps, rownames(tax_table(cleanps)) %in% c((H2Microbes$ASV)))


adonis2( phyloseq::distance(HeretibleASV, method="bray") ~ Pedigree*Time+Rep+Block,
         data = metadata, strata = metadata$Block)

##R2 improves 
# adonis2(formula = phyloseq::distance(HeretibleASV, method = "bray") ~ Pedigree * Time + Rep + Block, data = metadata, strata = metadata$Block)
# Df SumOfSqs      R2      F Pr(>F)    
# Pedigree       44   12.192 0.17834 1.6754  0.018 *  3
#   Time            1    0.694 0.01016 4.1980  0.001 ***
#   Rep             1    0.154 0.00225 0.9294  0.530    
# Block         135   27.508 0.40239 1.2320  0.019 *  
#   Pedigree:Time  43    7.140 0.10445 1.0040  0.517    
# Residual      125   20.673 0.30241                  
# Total         349   68.361 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

HeretibleASV<-subset_samples(HeretibleASV, Pedigree!="Bulk")

plot_bar(HeretibleASV, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")+
  facet_wrap(~Time,scales="free_x")


#Does this work? Should probably do this in a different part of the pipeline
RareHer<-rarefy_even_depth(HeretibleASV, replace=FALSE)
#RareHer = rarefy_even_depth(HeretibleASV, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)

plot_bar(RareHer, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")+
  facet_wrap(~Time,scales="free_x")

##This is cool it may be useful to understand these values 
##SO it we can determine which taxa are hereeibtl-

# It may be interesting to join these H2 values with the OTU table to determine the taxanomic compositon


##https://github.com/jyanglab/Maize_Rhizobiome_2022/blob/main/analysis/heritability_analysis.R
#Has code for P values using permucations

# Function to calculate H2
h2dat_stdN

getH2 <- function(h2dat) {
  #h2dat <- h2dat_stdN
  traits <- colnames(h2dat)[6:ncol(h2dat)]
  #print(traits)
  H2_df <- data.frame()
  
  for(trait in traits){
    #trait <- "asv_000013"
    f <- formula(paste0(trait, ' ~ (1|Pedigree) + (1|Block)'))
    fit <- suppressMessages(lmer(f, data=h2dat))
    v <- as.data.frame(VarCorr(fit))
    Vg <- v$vcov[v$grp == "Pedigree"]
    Ve <- v$vcov[v$grp == "Residual"]
    H2 <- round(Vg/(Vg + Ve/6), 6)
    H2_df <- rbind(H2_df, data.frame("ASV"=trait, "H2"=H2))
  }
  
  return(H2_df)
  
}

#h2dat <- h2dat_stdN



# function to shuffle abundances
shuffle <- function(h2dat) {
  h2dat$Pedigree <- sample(h2dat$Pedigree)
  return(h2dat)
}



#### calculate permutation p_values stdN

set.seed(2021)
H2_permutations_stdN <-  getH2(h2dat_stdN)
colnames(H2_permutations_stdN)[2] <- "obs"

for( i in c(1:40)){
  print(paste("permutation", i))
  perm <-  getH2(shuffle(h2dat_stdN))
  colnames(perm)[2] <- paste0("p",i)
  H2_permutations_stdN <- suppressMessages(left_join(H2_permutations_stdN, perm))
}

p_values_stdN <- H2_permutations_stdN %>%
  rowwise() %>%
  mutate(perm_p_stdN = (sum(c_across(3:length(H2_permutations_stdN))> obs)+1)/(length(H2_permutations_stdN)-1)) %>%
  dplyr::select(ASV, perm_p_stdN)


#Here we get a p vlaue form mulitpe test. 

test<-left_join(
  p_values_stdN,
  H2_stdN,
  by = "ASV",

)

taxatable<-as.data.frame(NILsV4@tax_table)

taxatable[,8]<-rownames(taxatable)
colnames(taxatable)[8]<-"ASV"


test2<-left_join(
  test,
  taxatable,
  by = "ASV",
  
)

##This is def an interesting methods
##I think what we are seeing here is that the archeal nitrifer may be the one that is being suppressed

##We are seeing high H2 of 
##Pull out Taxanomic Data can compare to the Values 




######


####Creating the Files for Fast X QTL======
setwd("~/Documents/Thesis /Chapter 4/Science /Rstats")
library(tidyverse)
library(phyloseq)
#saveRDS(ps, file = "NILcleanASV.OTU")
##This file has been cleaned, IE removed mitocondria, and chloroplast
ASVTable<- readRDS("NILcleanASV.OTU")
ASVtab<-subset_samples(ASVTable, Pedigree!="Bulk")
asv_table <- rownames_to_column(data.frame(otu_table(ASVtab)), var = "Sample.ID")
#Need to Transpose this Dataframe 

taxatab<-rownames_to_column(data.frame(tax_table(ASVtab)), var = "Sample.ID")

##Filtering OTUS greater then 0.005%
FSr  = transform_sample_counts(ASVTable, function(x) x / sum(x) )
FSfr = filter_taxa(FSr, function(x) sum(x) > .005, TRUE)
Filtered.Transformed.ASV <- rownames_to_column(data.frame(otu_table(FSfr)), var = "Sample.ID")

##Removing Bulk
NA.Soils.ASV<-asv_table[1:351,]
Filtered.Soils.ASV<-Filtered.Transformed.ASV[1:351,]

##Transposing 
Filtered.Soils.ASV.T<-t(Filtered.Soils.ASV)
df<-as.data.frame(t(asv_table))

df<-as.data.frame(Filtered.Soils.ASV.T)

#Correcting sample Names
testdf<-df %>% 
  mutate_all(funs(str_replace_all(., "_", ".")))

testdf<-(testdf[,-265])
#The writing table instructions are VERY VERY important, if they are slighly off they will be unable to be loaded into the file
write.table(testdf, "~/Documents/Thesis /Chapter 4/Science /Rstats/Old R-Analysis /R-Scripts/Matrix_eQTL/Full.ASV.Genes.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote=FALSE)


write.csv(taxatab,"~/Documents/Thesis /Chapter 4/Science /Rstats/Old R-Analysis /R-Scripts/Matrix_eQTL/ASV.taxa.csv", row.names = FALSE )
write.csv(testdf, "ASV.Tansformed.csv")

Trans<-read_csv("~/Documents/Thesis /Chapter 4/Science /Rstats/ASV.Tansformed.csv")


###Lets loook at our i
TeoMaiz<-subset_samples(cleanps, Pedigree=="B73"|Pedigree=="PI384071")




####Fungal ITS ------

##DADA2 OUTPUT ANALYSIS =====

library(ggplot2)
library(tidyverse)
library(phyloseq)
library(vegan)

setwd("~/Documents/Thesis /Chapter 4/Science /Rstats")

seqtab<- readRDS("NIL2018.ASV/ITS_seqtab_final.RDS")
summaryDada<- readRDS("NIL2018.ASV/ITS_summary-of-dada2.RDS")
taxa<- readRDS("NIL2018.ASV/ITS_taxonomy_final.RDS")
samdf<-read_csv("NIL2018.ASV/NIL2018MappingDNA.csv")


##Cleaning up Data
rownames(seqtab)<-rownames(seqtab)%>%str_replace("ITS3_ITS4-", "")%>%str_extract(".*(?=\\_)")

samdf<-as.data.frame(samdf)
rownames(samdf)<-samdf$`Sample ID`


##Need to make metadata for this seed stuff.
##The sample data thing is here
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps@sam_data

sample_names(ps)
##Rename barcode as ASV
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


#We see a clear partitioning of the seeds spray microbiome
plot_richness(ps, measures=c("Shannon", "Simpson"), x="Time")

fungiits<-ps

fungiits<-subset_samples(fungiits, Sample.ID!="NIL_T1_039")


GP.ord <- ordinate(fungiits, "NMDS", "bray",)
p1 = plot_ordination(fungiits, GP.ord, color="Time",  title="NIL Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Sample.ID), size = 5, vjust = 1.5)
print(p1)

fungiits%>%plot_ordination( GP.ord, color="Pedigree",  title="NIL Ordination")+
  geom_point(size=3)+
  #geom_text(mapping = aes(label = Sample.ID), size = 5, vjust = 1.5)+
  stat_ellipse()

print(p1)

fungiits<-subset_samples(fungiits, Sample.ID!="NIL_T1_039")
fungiits<-subset_samples(fungiits, Pedigree!="Bulk")


metadata <- as(sample_data(fungiits), "data.frame")
adonis2( phyloseq::distance(fungiits, method="bray") ~ Pedigree*Time+Rep+Block,
         data = metadata)


#The fungal ITS  is close to signifcance but not very well

ITS.T1<-subset_samples(fungiits, Time!="T2")

metadata <- as(sample_data(ITS.T1), "data.frame")
adonis2( phyloseq::distance(ITS.T1, method="bray") ~ Pedigree+Rep+Range+Row+Block,
         data = metadata)


ITS.T2<-subset_samples(fungiits, Time!="T1")
metadata <- as(sample_data(ITS.T2), "data.frame")
adonis2( phyloseq::distance(ITS.T2, method="bray") ~ Pedigree+Rep+Range+Row+Block,
         data = metadata)

##LEts filter out the OTU

##Lets look at the AMF and EMF groups
AMF<-subset_taxa(fungiits, Order=="o__Glomerales")
head(otu_table(AMF), 15) 
 
plot_richness(AMF, measures=c("Shannon" ), x="Pedigree")

plot_bar(AMF, fill="Family",x="Pedigree")+
  geom_bar(stat="identity")#+
  facet_wrap(~Time,scales="free_x")

plot_bar(AMF, fill="Genus",x="Pedigree")+
  geom_bar(stat="identity")+
  facet_wrap(~Time,scales="free_x")
  
  
AMFPedigree<-subset_samples(AMF, Pedigree=="Z031E0009"|Pedigree=="B73" |Pedigree=="PI384071"| Pedigree=="B73xPI384071"| 
                              Pedigree=="Z031E0061")
  
AMFPedigree<-subset_samples(AMF, Pedigree=="B73" |Pedigree=="PI384071"| Pedigree=="B73xPI384071")

plot_bar(AMFPedigree, fill="Genus",x="Pedigree", facet_grid=~Genus)+
  geom_bar(stat="identity")


plot_bar(AMFPedigree, fill="Genus",x="Pedigree")+
  geom_bar(stat="identity")

plot_bar(AMF, fill="Genus",x="Pedigree")+
  geom_bar(stat="identity")

Fun<-subset_taxa(fungiits, Genus=="g__Funneliformis")

Micro<-subset_taxa(fungiits, Genus=="g__Microdominikia")


plot_bar(Fun, fill="Genus",x="Pedigree")+
  geom_bar(stat="identity")

plot_bar(Micro, fill="Genus",x="Pedigree")+
  geom_bar(stat="identity")

AMFPedigree<-subset_samples(AMF, Pedigree=="B73" |Pedigree=="PI384071"| Pedigree=="B73xPI384071")

AMFT2<-subset_samples(AMFPedigree, Time=="T2")
AMFT1<-subset_samples(AMFPedigree, Time=="T1")

plot_bar(AMFT2, fill="Order",x="Pedigree")+
  geom_bar(stat="identity")

plot_bar(AMFT2, fill="Genus",x="Pedigree")+
  geom_bar(stat="identity")


plot_bar(AMFT1, fill="Genus",x="Pedigree")+
  geom_bar(stat="identity")

AMFPedigree<-subset_samples(AMF, Pedigree!="B73")
AMFPedigree<-subset_samples(AMFPedigree, Pedigree!="PI384071")
AMFPedigree<-subset_samples(AMFPedigree, Pedigree!="B73xPI384071")

AMFT2<-subset_samples(AMFPedigree, Time=="T2")
plot_bar(AMFT2, fill="Order",x="Pedigree")+
  geom_bar(stat="identity")

###Stats for fungi ====
metadata <- as(sample_data(AMF), "data.frame")
adonis2( phyloseq::distance(AMF, method="bray") ~ Pedigree+Time+Rep+Block,
         data = metadata)


x2 = transform_sample_counts(AMF, rank)
head(otu_table(x2), 15)


adonis2( phyloseq::distance(x2, method="bray") ~ Pedigree+Time+Rep+Block,
         data = metadata)


##Filtering OTUS greater then 0.005%
FSr  = transform_sample_counts(fungiits, function(x) x / sum(x) )
FSfr = filter_taxa(FSr, function(x) sum(x) > .005, TRUE)

Fungi.05<-FSfr

##I'm going to clean up the data by removing the EXTREMELY low abundance outs from the data
metadata <- as(sample_data(Fungi.05), "data.frame")
adonis2( phyloseq::distance(Fungi.05, method="bray") ~ Pedigree*Time+Rep+Block,
         data = metadata, strata = metadata$Block)



# Df SumOfSqs      R2      F Pr(>F)    
# Pedigree       43   17.993 0.12464 1.0173  0.003 ** 
#   Time            1    1.616 0.01119 3.9282  0.001 ***
#   Rep             1    0.384 0.00266 0.9343  0.377    
# Block         132   56.580 0.39193 1.0421  0.004 ** 
#   Pedigree:Time  43   17.605 0.12195 0.9954  0.536    
# Residual      122   50.183 0.34762                  
# Total         342  144.362 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##I'm going to clean up the data by removing the EXTREMELY low abundance outs from the data
plot_bar(Fungi.05, fill="Phylum",x="Pedigree")+
  geom_bar(stat="identity")#+
facet_wrap(~Time,scales="free_x")
           
           
AMF<-subset_taxa(Fungi.05, Order=="o__Glomerales")

plot_bar(AMF, fill="Family",x="Pedigree")+
  geom_bar(stat="identity")#+
facet_wrap(~Time,scales="free_x")

NIL.05<-subset_samples(Fungi.05, Pedigree!="B73"| Pedigree!="B73xP1384071"|Pedigree!="P1384071")

metadata <- as(sample_data(NIL.05), "data.frame")
adonis2( phyloseq::distance(NIL.05, method="bray") ~ Pedigree*Time+Rep+Block,
         data = metadata, strata = metadata$Block)
# Df SumOfSqs      R2      F Pr(>F)    
# Pedigree       43   17.993 0.12464 1.0173  0.003 ** 
#   Time            1    1.616 0.01119 3.9282  0.001 ***
#   Rep             1    0.384 0.00266 0.9343  0.418    
# Block         132   56.580 0.39193 1.0421  0.003 ** 
#   Pedigree:Time  43   17.605 0.12195 0.9954  0.522    
# Residual      122   50.183 0.34762                  
# Total         342  144.362 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##Here we are testing just the introgession

