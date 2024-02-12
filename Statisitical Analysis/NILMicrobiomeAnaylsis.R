#Near-Isogenic Line Microbiome #First look through analysis and plan
##Alonso Favela 1/17/18
## I'm going to start by just looking at the data visually and then I will examine the data satisitcally. 
### To do this im going to run them in phyloseq and close out with perfroming a vegan anova etc. 

library(phyloseq)
library(ggplot2)
library(vegan)
#I'm a little worried that my data is being effect by the level of rarefraction. 
#Seems like I might need to use the raw data to make it the statisics work better.
#

####16S rRNA ======
#Set working directory

setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL Maize Microbiome /R-Scripts/NIL CSV")
setwd("~/Documents/Projects/Research/Maize Microbiome/NIL MM2018/R-Scripts/NIL CSV")


###Phyloseq Code
##########
#Importing Biom File
biomOTU = import_biom('Phyloseq/otu_table_16S_30152.biom',treefilename = 'Phyloseq/rep_seq_aligned_pfiltered.tre', 
                      parseFunction = parse_taxonomy_greengenes)
#Importaning Mapping
NILMap = import_qiime_sample_data("Phyloseq/NIL2018MappingDNA.txt")
#Functional Bins are present
NILMap = import_qiime_sample_data("Sample_map.txt")
#NILMap = import_qiime_sample_data("Sample_map_ActualHighLow.txt")

#Setting themes 
theme_set(theme_bw())

##Loading and Editing Data ====
#Combining Biom and Mapping file to Make the phyloseq file
NILOTU = merge_phyloseq(biomOTU, NILMap)

NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.090")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.069")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.086")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T1.039")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.050")
NILOTU <- subset_samples(NILOTU, Pedigree != "Bulk")

#Removing Teosinte&Hybrid
NILOTU <- subset_samples(NILOTU, Pedigree != "PI384071")
NILOTU <- subset_samples(NILOTU, Pedigree != "B73xPI384071")

#Subet Data Timepoints
NILT1 <- subset_samples(NILOTU, Time != "T2")
NILT2 <- subset_samples(NILOTU, Time != "T1")
#THESE NEED TO BE REMOVED
# NIL.T2.090
# NIL.T2.069
# NIL.T2.086
#OUTLIERS
# NIL.T1.039
# NIL.T2.050


#Here I will be doing my comparison in the functional bins
#figure results are based on NILOTU- so all of the results
FunctionNiBin <- subset_samples(NILOTU, Nitrification == "Low"|Nitrification =="High")
FunctionDeNBin <- subset_samples(NILOTU, Denitrification == "Low"|Denitrification =="High")

#Time point 2 results make more sense to present and basically echo the same result
FunctionNiBin <- subset_samples(NILT2, Nitrification == "Low"|Nitrification =="High")
FunctionDeNBin <- subset_samples(NILT2, Denitrification == "Low"|Denitrification =="High")

# 
# library("DESeq2")
# ###Nitrifcation BIN
# diagdds = phyloseq_to_deseq2(FunctionNiBin, ~ Nitrification)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
# 
# res = results(diagdds, cooksCutoff = FALSE)
# alpha = 0.10
# sigtab = res[which(res$padj < alpha), ]
# sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(FunctionNiBin)[rownames(sigtab), ], "matrix"))
# head(sigtab)
# 
# dim(sigtab)

#Lookinf at the Denitrfication

# diagdds = phyloseq_to_deseq2(FunctionDeNBin, ~ Denitrification)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
# 
# res = results(diagdds, cooksCutoff = FALSE)
# alpha = 0.05
# sigtab = res[which(res$padj < alpha), ]
# sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(FunctionNiBin)[rownames(sigtab), ], "matrix"))
# head(sigtab)
# 
# dim(sigtab)
# 
# NILOTU = merge_phyloseq(biomOTU, NILMap)


#Looking at just the nitrfiying bacteria 
#make sure to transform sample counts before subsetting taxa 
# I need to alter a part of this subsetting| Let me double check with angelas version
#Plot figures 6x8 
OTU_merge = merge_samples(FunctionNiBin, "Nitrification", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTU_merge)$Nitrification <- levels(sample_data(FunctionNiBin)$Nitrification)
OTU_rel_merge = transform_sample_counts(OTU_merge, function(x) 100 * x/sum(x) )
NitrifersOTU=subset_taxa(OTU_rel_merge, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")
plot_bar(NitrifersOTU, "Nitrification", fill = "Order", facet_grid = ~Class)+labs(title = "nitrfier taxa T2 Ave-Trans")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()#+geom_bar(stat="identity")
#I'm going to try transforming, averagign, subsetting
OTU_merge = transform_sample_counts(FunctionNiBin, function(x) 100 * x/sum(x) )
OTU_rel_merge = merge_samples(OTU_merge, "Nitrification", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTU_rel_merge)$Nitrification <- levels(sample_data(FunctionNiBin)$Nitrification)
NitrifersOTU=subset_taxa(OTU_rel_merge, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")
plot_bar(NitrifersOTU, "Nitrification", fill = "Order", facet_grid = ~Class)+labs(title = "nitrifier taxa T2 Trans-Ave")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()#+geom_bar(stat="identity")

#This file here shows the clearest difference in treatment
#OTU_rel_merge= transform_sample_counts(FunctionNiBin, function(x) 100* x/sum(x) )
NitrifersOTU=subset_taxa(OTU_rel_merge, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")
NitrifersOTU=subset_taxa(OTU_merge, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")


# FunctionNiBin <- subset_samples(NitrifersOTU, Nitrification == "Low"|Nitrification =="High")
# FunctionNitrosoBin <- subset_samples(NitrosoOTU, Nitrification == "Low"|Nitrification =="High")

plot_bar(NitrifersOTU, "Nitrification", fill = "Order", facet_grid = ~Class)+geom_bar(stat="identity")+labs(title = "nitrifier taxa")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()
plot_bar(NitrosoOTU, "Nitrification", fill = "Order", facet_grid = ~Class)+geom_bar(stat="identity")+labs(title = "nitrifier taxa")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))

#Trying something funky
# FunctionNiBin2=subset_taxa(FunctionNiBin, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")
# OTU_merge = merge_samples(FunctionNiBin2, "Nitrification", fun = mean) #averages each OTU in all samples belonging to each habitat class
# sample_data(OTU_merge)$Nitrification <- levels(sample_data(FunctionNiBin)$Nitrification)
# OTU_rel_merge = transform_sample_counts(OTU_merge, function(x)   x/sum(x) )
# plot_bar(OTU_rel_merge, "Nitrification", fill = "Order", facet_grid = ~Class)+geom_bar(stat="identity")+labs(title = "nitrifier taxa")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()


NitroT1 <- subset_samples(FunctionNiBin, Time != "T2")
NitroT2 <- subset_samples(FunctionNiBin, Time != "T1")

GP.ord <- ordinate(FunctionNiBin, "NMDS", "bray",)
p1 = plot_ordination(FunctionNiBin, GP.ord, type="Samples", color="Nitrification", title="Nitrifiers Ordination")+geom_point(size=3)+
print(p1)

GP.ord <- ordinate(FunctionNitrosoBin, "NMDS", "bray",)
p1 = plot_ordination(FunctionNitrosoBin, GP.ord, type="Samples", color="Nitrification", title="Nitrifiers Ordination")+geom_point(size=3)+
  print(p1)

TopNOTUs = names(sort(taxa_sums(FunctionNiBin), TRUE) [1:30])
OTU_merge = merge_samples(FunctionNiBin, "Nitrification", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTU_merge)$Nitrification <- levels(sample_data(FunctionNiBin)$Nitrification)
OTU_rel_merge = transform_sample_counts(OTU_merge, function(x) 100 * x/sum(x) ) 
Project_10_rel_merge = prune_taxa(TopNOTUs, OTU_rel_merge)
Project_10_rel_merge #confirm that the phyloseq object has been subsetted
taxa_names(Project_10_rel_merge) 


plot_bar(FunctionNiBin, "Nitrification", fill = "Order", facet_grid = ~Class)+geom_bar(stat="identity")+labs(title = "Nitrifying Micorbes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))
plot_bar(NitrifersOTU, "Nitrification", fill = "Order", facet_grid = ~Class)+geom_bar(stat="identity")+labs(title = "nitrifier taxa")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))

plot_bar(Project_10_rel_merge, "Nitrification", fill = "OTU", facet_grid = ~Family)+geom_bar(stat="identity")+labs(title = "10 Most Abundant 16S OTUs")+xlab("Treatment")
plot_richness(FunctionNiBin, x = "Nitrification", measures = "Shannon") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")
plot_richness(FunctionNiBin, x = "Nitrification", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")


quartz()
title = "nitrifier taxa"
plot_bar(NitrifersOTU, "Nitrification", fill = "Phylum",
         title=title, facet_grid = "Class~.")+
  theme(text=element_text(family = "Helvetica",size=14))




#Alpha Diversity Plots
#c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
plot_richness(NILOTU, x = "Pedigree", measures = c("Observed", "Chao1", "Simpson", "Shannon"), sortby = "Shannon") + geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 90))+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA Diversity")
#It would be interesting to take into account the amount of enriched diversity and the lowered diversity to see if there are speciffc genetic regions that are enriching the C and some that are loweing it. 
plot_richness(NILOTU, x = "Pedigree", measures = "Shannon", sortby = "Shannon") + geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 90))+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA Diversity")

plot_richness(NILT1, x = "Pedigree", measures = "Shannon", sortby = "Shannon") +
  geom_boxplot()+theme_classic()+
  theme(strip.text.x = element_blank(),axis.text.x = element_text(angle = 90))+
  labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

plot_richness(NILT2, x = "Pedigree", measures = "Chao1", sortby="Chao1") + geom_boxplot()+theme_classic()+
  theme(strip.text.x = element_blank(),axis.text.x = element_text(angle = 90))+
  ggtitle("16S rRNA")
plot_richness(NILT2, x = "Pedigree", measures = "InvSimpson", sortby="InvSimpson") + geom_boxplot()+theme_classic()+
  theme(strip.text.x = element_blank(),axis.text.x = element_text(angle = 90))+
  ggtitle("16S rRNA")

plot_richness(NILT2, x = "Pedigree", measures = "Fisher", sortby="Fisher") + geom_boxplot()+theme_classic()+
  theme(strip.text.x = element_blank(),axis.text.x = element_text(angle = 90))+
  ggtitle("16S rRNA")


##Genotype Averages=====
GenotypeMean = merge_samples(NILOTU, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(GenotypeMean)$Pedigree <- levels(sample_data(NILOTU)$Pedigree)

T1Mean = merge_samples(NILT1, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(T1Mean)$Pedigree <- levels(sample_data(NILT1)$Pedigree)

T2Mean = merge_samples(NILT2, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(T2Mean)$Pedigree <- levels(sample_data(NILT2)$Pedigree)

GxT = merge_samples(NILOTU, "GxT", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(GxT)$GxT <- levels(sample_data(NILOTU)$GxT)
##This is not working at all. Its alining to different lines and is garbage
sample_data(GxT)$Pedigree <- as.factor(sample_data(GxT)$Pedigree)
sample_data(GxT)$Time <- as.factor(sample_data(GxT)$Time)

###Finding the mean values and ploting diversity
plot_richness(GenotypeMean, x = "Pedigree", measures = "Chao1", sortby = "Chao1") + geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 90))+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA Diversity")

plot_richness(T1Mean, x = "Pedigree", measures = "Chao1", sortby = "Chao1") + geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 90))+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA Diversity")

plot_richness(T2Mean, x = "Pedigree", measures = "Chao1", sortby = "Chao1") + geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 90))+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA Diversity")

plot_richness(GxT, x = "GxT", measures = "Chao1", sortby = "Chao1") + geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 90))+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA Diversity")

#Beta Diversity 
GP.ord <- ordinate(NILOTU, "RDA", "bray",)
p1 = plot_ordination(NILOTU, GP.ord, type="Samples", color="Pedigree", title="Bacterial 16S Ordination")+geom_point(size=3)
print(p1)
  
GP.ord <- ordinate(NILOTU, "RDA", "bray",)
p1 = plot_ordination(NILOTU, GP.ord, type="Samples", color="Time",  title="Bacterial 16S Ordination")+geom_point(size=3)+stat_ellipse()
print(p1)


GP.ord <- ordinate(NILT1, "RDA", "bray",)
p1 = plot_ordination(NILT1, GP.ord, type="Samples", color="Pedigree",  title="Bacterial 16S Ordination")+geom_point(size=3)#+stat_ellipse()
print(p1)

GP.ord <- ordinate(NILT2, "RDA", "bray",)
p1 = plot_ordination(NILT2, GP.ord, type="Samples", color="Pedigree",  title="Bacterial 16S Ordination")+geom_point(size=3)+  geom_text(mapping = aes(label = Pedigree))
labs(color="Treatment")
print(p1)

#This doesnt take a true mean, its actually a sum of the values 
GP.ord <- ordinate(GenotypeMean, "RDA", "bray",)
p1 = plot_ordination(GenotypeMean, GP.ord, type="Samples", color="Pedigree",  title="Bacterial 16S Ordination")+geom_point(size=3)+
  geom_text(mapping = aes(label = Pedigree))
print(p1)
##Genotype by Time effects
GP.ord <- ordinate(GxT, "RDA", "bray",)
p1 = plot_ordination(GxT, GP.ord, type="Samples", color="Time", title="Bacterial 16S Ordination")+geom_point(size=3)+
   geom_line(aes(fill=Time))+  geom_text(mapping = aes(label = GxT))
print(p1)

p1 = plot_ordination(GxT, GP.ord, type="Samples", color="Pedigree", title="Bacterial 16S Ordination")+geom_point(size=3)+
  geom_text(mapping = aes(label = GxT))
print(p1)


GP.ord <- ordinate(T1Mean, "RDA", "bray",)
p1 = plot_ordination(T1Mean, GP.ord, type="Samples", color="Pedigree",  title="Bacterial 16S Ordination")+geom_point(size=3)+
  geom_text(mapping = aes(label = Pedigree))
print(p1)

GP.ord <- ordinate(T2Mean, "RDA", "bray",)
p1 = plot_ordination(T2Mean, GP.ord, type="Samples", color="Pedigree",  title="Bacterial 16S Ordination")+geom_point(size=3)+
  geom_text(mapping = aes(label = Pedigree))
print(p1)


GP.ord <- ordinate(NILT2, "RDA", "bray",)
p2 = plot_ordination(NILT2, GP.ord, type="Samples", color="Pedigree",  title="Bacterial 16S Ordination")+geom_point(size=3)+  geom_text(mapping = aes(label = Pedigree))
labs(color="Treatment")
P1<-p1+p2
print(P3)


#Top OTU All samples 
TopNOTUs = names(sort(taxa_sums(GenotypeMean), TRUE) [1:10])
OTU_merge = merge_samples(GenotypeMean, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTU_merge)$Pedigree <- levels(sample_data(NILOTU)$Pedigree)
OTU_rel_merge = transform_sample_counts(OTU_merge, function(x) 100 * x/sum(x) ) 


#use the top 10 list from above to select the same OTUs from the averaged and relativeized data
Project_10_rel_merge = prune_taxa(TopNOTUs, OTU_rel_merge)
Project_10_rel_merge #confirm that the phyloseq object has been subsetted
taxa_names(Project_10_rel_merge) 

plot_bar(Project_10_rel_merge, "Pedigree", fill = "Family", facet_grid = ~OTU)+geom_bar(stat="identity")+labs(title = "10 Most Abundant 16S OTUs")+xlab("Treatment")

GP.ord <- ordinate(OTU_merge, "NMDS", "bray",)
p1 = plot_ordination(OTU_merge, GP.ord, type="Samples", color="Pedigree",  title="Bacterial 16S Ordination")+geom_point(size=3)+stat_ellipse()+
  labs(color="Treatment")
print(p1)

#These micorbial communites seem to be really similar.


##Looking to see if methanotrophs follow the same patter as the nitrifing bacteria
MethyOTU=subset_taxa(FunctionNiBin, Family == "Methylobacteriaceae"|Family=="Methylocystaceae"|Class=="[Methylacidiphilae]")
Methy2OTU=subset_taxa(FunctionNiBin, Family=="Methylocystaceae")

plot_bar(Methy2OTU, "Nitrification", fill = "OTU", facet_grid = ~Genus)+geom_bar(stat="identity")+labs(title = "Methy taxa")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))



#Im going to do some basic stats to see how the treatments boor out.
## TO reduce the varaince in the samples I'm going to alter stuff. I'm going to look at the MC at the genus level 

#Importing Biom File
biomOTU = import_biom('Phyloseq/Sample_otu_table_L6.biom')
#Importaning Mapping
NILMap = import_qiime_sample_data("Phyloseq/Sample_map.txt")

#Combining Biom and Mapping file to Make the phyloseq file
NILOTU = merge_phyloseq(biomOTU, NILMap)
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.090")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.069")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.086")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T1.039")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.050")
NILOTU <- subset_samples(NILOTU, Pedigree != "Bulk")

GP.ord <- ordinate(NILOTU, "RDA", "bray",)
p1 = plot_ordination(NILOTU, GP.ord, type="Samples", color="Pedigree",  title="Bacterial 16S Ordination")+geom_point(size=3)+stat_ellipse()+
  labs(color="Treatment")
print(p1)



#Here I just wanted to look at how the micorbail communiy supression lines differed
SupressionLines <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0021"|Pedigree == "Z031E0047"|Pedigree == "B73xPI384071"|Pedigree == "PI384071")


GP.ord <- ordinate(SupressionLines, "NMDS", "bray",)
p1 = plot_ordination(SupressionLines, GP.ord, type="Samples", color="Pedigree",  title="16S rRNA - Supression NILs")+
  geom_point(size=3)+stat_ellipse(geom = "polygon", alpha=1/100)+
  labs(color="Treatment")+geom_text(mapping = aes(label = Time), size=3,hjust=0, vjust=0)
print(p1)

##Here I would like to look at the maize-teosinte-hybrid comparison =====

OTUF1TeoB73 <- subset_samples(NILOTU, Pedigree == "PI384071"|Pedigree == "B73xPI384071"|Pedigree == "B73")

OTUF1TeoB73T1 <- subset_samples(OTUF1TeoB73, Time=="T1")
OTUF1TeoB73T2 <- subset_samples(OTUF1TeoB73, Time=="T2")

plot_richness(OTUF1TeoB73, x = "Pedigree", measures = "Shannon") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Genotype", y="Alpha diversity (chao1)")+ggtitle("16S rRNA")
plot_richness(OTUF1TeoB73T1, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Genotyoe", y="Alpha diversity (chao1)")+ggtitle("Diversity Time 1")
#In this third time point there does seem to be some speration
plot_richness(OTUF1TeoB73T2, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Genotype", y="Alpha diversity (chao1)")+ggtitle("Diversity Time 2")


##ordination
GP.ord <- ordinate(OTUF1TeoB73, "NMDS", "bray",)
p1 = plot_ordination(OTUF1TeoB73, GP.ord, type="Samples", color="Pedigree",  title="16S Ordination")+geom_point(size=3)+
  labs(color="Treatment")+theme_classic()+scale_color_manual(values=c("Grey", "Black", "White"))+geom_point(shape = 1,size = 3,colour = "black")+guides(fill="legend")
print(p1)

GP.ord <- ordinate(OTUF1TeoB73T1, "NMDS", "bray",)
p1 = plot_ordination(OTUF1TeoB73T1, GP.ord, type="Samples", color="Pedigree",  title="16S Ordination")+geom_point(size=3)+
  labs(color="Treatment")+theme_classic()+scale_color_manual(values=c("Grey", "Black", "White"))+geom_point(shape = 1,size = 3,colour = "black")+guides(fill="legend")
print(p1)

GP.ord <- ordinate(OTUF1TeoB73T2, "NMDS", "bray",)
p1 = plot_ordination(OTUF1TeoB73T2, GP.ord, type="Samples", color="Pedigree",  title="16S Ordination")+geom_point(size=3)+
  labs(color="Treatment")+theme_classic()+scale_color_manual(values=c("Grey", "Black", "White"))+geom_point(shape = 1,size = 3,colour = "black")+guides(fill="legend")+stat_ellipse()
print(p1)


###Lets look at top OTUs
TopNOTUs = names(sort(taxa_sums(OTUF1TeoB73T2), TRUE) [1:40])
OTU_merge = merge_samples(OTUF1TeoB73, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTU_merge)$Pedigree <- levels(sample_data(OTUF1TeoB73)$Pedigree)
OTU_rel_merge = transform_sample_counts(OTU_merge, function(x) 100 * x/sum(x) ) 


#use the top 10 list from above to select the same OTUs from the averaged and relativeized data
Project_10_rel_merge = prune_taxa(TopNOTUs, OTU_rel_merge)
Project_10_rel_merge #confirm that the phyloseq object has been subsetted
taxa_names(Project_10_rel_merge) 

plot_bar(Project_10_rel_merge, "Pedigree", fill = "Order", facet_grid = ~Phylum)+geom_bar(stat="identity")+labs(title = "10 Most Abundant 16S OTUs")+xlab("Treatment")

plot_bar(OTUF1TeoB73T2, "Pedigree", fill = "Phylum", facet_grid = ~Phylum)+geom_bar(stat="identity")+labs(title = "10 Most Abundant 16S OTUs")+xlab("Treatment")


###########
######Functional Gene Phyoseq =======

setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL Maize Microbiome/R-Scripts/NIL CSV")

###nirK Code====
##########
#Importing Biom File
nirKOTU = import_biom('Phyloseq/otu_table_nirK_5004.biom')
NILMap = import_qiime_sample_data("Sample_map.txt")
NILOTU = merge_phyloseq(nirKOTU, NILMap)
FunctionDeNBin <- subset_samples(NILOTU, Denitrification == "Low"|Denitrification =="High")




plot_richness(FunctionDeNBin, x = "Denitrification", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("nirK")


GP.ord <- ordinate(FunctionDeNBin, "NMDS", "bray",)
p1 = plot_ordination(FunctionDeNBin, GP.ord, type="Samples", color="Denitrification", title="nirK genes")+geom_point(size=3)
print(p1)


plot_bar(FunctionDeNBin, "Denitrification", fill = "Rank4", facet_grid = ~Rank3)+geom_bar(stat="identity")+labs(title = "DEA Micorbes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))
plot_bar(FunctionDeNBin, "Denitrification", fill = "OTU", facet_grid = ~Rank2)+geom_bar(stat="identity")+labs(title = "DEA Micorbes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme(legend.position = "none")


#Nirk Differencs in the abundance
FunctionDeNBin10 = names(sort(taxa_sums(FunctionDeNBin), TRUE) [1:10])
FunctionDeNBin100 = prune_taxa(FunctionDeNBin10, FunctionDeNBin)
plot_bar(FunctionDeNBin100, "Denitrification", fill = "OTU", facet_grid = ~Rank2)+geom_bar(stat="identity")+labs(title = "NirK Top Genes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))

plot_bar(FunctionDeNBin100, "Denitrification", fill = "Rank3", facet_grid = ~Rank5)+geom_bar(stat="identity")+labs(title = "NirK Top Genes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))

###nirS Code====
##########
#Importing Biom File
#Seems like nirS is higher in diversirty compared to the other genes
nirSOTU = import_biom('Phyloseq/otu_table_nirS_707.biom')
NILMap = import_qiime_sample_data("Sample_map.txt")
NILOTU = merge_phyloseq(nirSOTU, NILMap)
FunctionDeNBin <- subset_samples(NILOTU, Denitrification == "Low"|Denitrification =="High")

plot_richness(FunctionDeNBin, x = "Denitrification", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("nirS")


GP.ord <- ordinate(FunctionDeNBin, "NMDS", "bray",)
p1 = plot_ordination(FunctionDeNBin, GP.ord, type="Samples", color="Denitrification", title="nirK genes")+geom_point(size=3)
print(p1)

plot_bar(FunctionDeNBin, "Denitrification", fill = "Rank4", facet_grid = ~Rank3)+geom_bar(stat="identity")+labs(title = "DEA Micorbes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme(legend.position = "none")
plot_bar(FunctionDeNBin, "Denitrification", fill = "OTU", facet_grid = ~Rank2)+geom_bar(stat="identity")+labs(title = "DEA Micorbes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme(legend.position = "none")

#NirS Top 10
FunctionDeNBin10 = names(sort(taxa_sums(FunctionDeNBin), TRUE) [1:10])
FunctionDeNBin100 = prune_taxa(FunctionDeNBin10, FunctionDeNBin)
plot_bar(FunctionDeNBin100, "Denitrification", fill = "OTU", facet_grid = ~Rank5)+geom_bar(stat="identity")+labs(title = "NirS Top genes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))

plot_bar(FunctionDeNBin100, "Denitrification", fill = "OTU", facet_grid = ~Rank2)+geom_bar(stat="identity")+labs(title = "NirS Top genes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))


###nosZ Code====
##########
#Importing Biom File
nosZOTU = import_biom('Phyloseq/otu_table_nosZ_294.biom')
NILMap = import_qiime_sample_data("Sample_map.txt")
NILOTU = merge_phyloseq(nosZOTU, NILMap)
FunctionDeNBin <- subset_samples(NILOTU, Denitrification == "Low"|Denitrification =="High")

plot_richness(FunctionDeNBin, x = "Denitrification", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("nosZ")


GP.ord <- ordinate(FunctionDeNBin, "NMDS", "bray",)
p1 = plot_ordination(FunctionDeNBin, GP.ord, type="Samples", color="Denitrification", title="nosZ genes")+geom_point(size=3)
print(p1)



plot_bar(FunctionDeNBin, "Denitrification", fill = "OTU", facet_grid = ~Rank2)+geom_bar(stat="identity")+labs(title = "DEA Micorbes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme(legend.position = "none")

#Top 10 NosZ
FunctionDeNBin10 = names(sort(taxa_sums(FunctionDeNBin), TRUE) [1:10])
FunctionDeNBin100 = prune_taxa(FunctionDeNBin10, FunctionDeNBin)
plot_bar(FunctionDeNBin100, "Denitrification", fill = "OTU", facet_grid = ~Rank2)+geom_bar(stat="identity")+labs(title = "NosZ Top genes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))

plot_bar(FunctionDeNBin100, "Denitrification", fill = "OTU", facet_grid = ~Rank5)+geom_bar(stat="identity")+labs(title = "NosZ Top genes")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))


######Vegan Code
##### Make sure to remove outlier samples 
# NIL.T2.090
# NIL.T2.069
# NIL.T2.086
# NIL.T1.039
# NIL.T2.050
##########

#####ITS Fungal ======
#Importing Biom File
ITSOTU = import_biom('Phyloseq/otu_table_ITS_2302.biom')
#Importaning Mapping
NILMap = import_qiime_sample_data("NIL2018MappingDNA.txt")
#Functional Bins are present
NILMap = import_qiime_sample_data("Sample_map.txt")

ITS.NIL.OTU = merge_phyloseq(ITSOTU, NILMap)

#c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
plot_richness(ITS.NIL.OTU, x = "Pedigree", measures = "Observed", sortby = "Observed") +
  geom_boxplot()+theme_classic()+
  theme(strip.text.x = element_blank(),axis.text.x = element_text(angle = 90))+
  labs( x="Treatment", y="Alpha Diversity")+ggtitle("ITS rRNA")

# c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
GP.ord <- ordinate(ITS.NIL.OTU, "RDA", "bray",)
p1 = plot_ordination(ITS.NIL.OTU, GP.ord, type="Samples", color="Pedigree",  title="ITS Ordination")+geom_point(size=3)+stat_ellipse()+
  labs(color="Treatment")
print(p1)

NILOTU <- subset_samples(NILOTU, Pedigree != "Bulk")

##Here we will be looking at the genotypics effects on pedigree
ITS.Genotype.Mean = merge_samples(ITS.NIL.OTU, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(ITS.Genotype.Mean)$Pedigree <- levels(sample_data(ITS.NIL.OTU)$Pedigree)

#Interesting it seems like This data behaved much more during amplitifaction 
OTU_merge

##Looking at AMF By pulling them out
get_taxa_unique(ITS.NIL.OTU, "Rank2")
AMF=subset_taxa(ITS.NIL.OTU, Rank2 == "p__Glomeromycota")

plot_bar(AMF, "Pedigree", fill = "Rank4", "Abundance")


TopNOTUs = names(sort(taxa_sums(ITS.NIL.OTU), TRUE) [1:30])
OTU_merge = merge_samples(ITS.NIL.OTU, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTU_merge)$Pedigree <- levels(sample_data(ITS.NIL.OTU)$Pedigree)
OTU_rel_merge = transform_sample_counts(OTU_merge, function(x) 100 * x/sum(x) )

plot_bar(OTU_rel_merge, "Pedigree", fill = "Rank2", "Abundance")

AMF <- subset_samples(AMF, Pedigree != "Bulk")

GP.ord <- ordinate(AMF, "NMDS", "bray",)
p1 = plot_ordination(AMF, GP.ord, type="Samples", color="Pedigree",  title="AMF Ordination")+geom_point(size=3)+
  labs(color="Treatment")+geom_text(mapping = aes(label = Pedigree), size = 10, vjust = 1.5)
print(p1)


GP.ord <- ordinate(OTU_rel_merge, "NMDS", "bray",)
p1 = plot_ordination(OTU_rel_merge, GP.ord, type="Samples", color="Pedigree",  title="AMF Ordination")+geom_point(size=3)+
  labs(color="Treatment")+geom_text(mapping = aes(label = Pedigree), size = 10, vjust = 1.5)
print(p1)

##Fungal DESEQ
NILOTUE0PI384071 <- subset_samples(ITS.NIL.OTU, Pedigree == "B73"|Pedigree=="PI384071")


plot_richness(NILOTUE0PI384071, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

##Very Nice Seperation
GP.ord <- ordinate(NILOTUE0PI384071, "NMDS", "bray",)
p1 = plot_ordination(NILOTUE0PI384071, GP.ord, type="Samples", color="Pedigree",  title="Fungal ITS Ordination")+geom_point(size=3)+
  labs(color="Treatment")+theme_classic()+scale_color_manual(values=c("Grey", "Black", "White"))+geom_point(shape = 1,size = 3,colour = "black")+guides(fill="legend")
print(p1)
#Comparing
#NoN Signifcant results
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0PI384071, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE009)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)


###DeSeq2 Analysis ------
######


#Deseq is a method by which to compare two single samples to each other. Here I will perfrom Deseq2 between all NILs and b73.
#Hopefully this will give us more clarity 
# First I need to load the data 
#16S rRNA 
#Set working directory
setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL Maize Microbiome/R-Scripts/NIL CSV/")

setwd("~/Documents/Projects/Research/Maize Microbiome/NIL MM2018/R-Scripts/NIL CSV")
###Phyloseq Code
##########
#Importing Biom File
biomOTU = import_biom('Phyloseq/otu_table_16S_30152.biom',treefilename = 'Phyloseq/rep_seq_aligned_pfiltered.tre', 
                      parseFunction = parse_taxonomy_greengenes)
biomOTU = import_biom('Phyloseq/otu_table.biom',treefilename = 'Phyloseq/rep_seq_aligned_pfiltered.tre', 
                      parseFunction = parse_taxonomy_greengenes)


#Importaning Mapping
NILMap = import_qiime_sample_data("NIL2018MappingDNA.txt")

NILOTU = merge_phyloseq(biomOTU, NILMap)


NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.090")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.069")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.086")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T1.039")
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.T2.050")
##Fungi Coding

#Importing Biom File
ITSOTU = import_biom('Phyloseq/otu_table_ITS_2302.biom')
#Importaning Mapping
NILMap = import_qiime_sample_data("NIL2018MappingDNA.txt")

NILOTU = merge_phyloseq(ITSOTU, NILMap)

NILOTU  = transform_sample_counts(NILOTU, function(x) x+1 )

#Seperate by time
NILOTUT1 <- subset_samples(NILOTU,  Time == "T1")
NILOTUT2 <- subset_samples(NILOTU, Time == "T2")

##Here to run T2 comparison
NILOTU<-NILOTUT2
#Here run averages
NILOTU<-T2Mean


#This Code Removes treatments 
#Starting with E009 comparsion =====
NILOTUE009 <- subset_samples(NILOTU, Pedigree == "Z031E0009" |Pedigree == "B73" )


plot_richness(NILOTUE009, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

##Very Nice Seperation
GP.ord <- ordinate(NILOTUE009, "NMDS", "bray",)
p1 = plot_ordination(NILOTUE009, GP.ord, type="Samples", color="Pedigree",  title="Fungal ITS Ordination")+geom_point(size=3)+
  labs(color="Treatment")+theme_classic()+scale_color_manual(values=c("Grey", "Black", "White"))+geom_point(shape = 1,size = 3,colour = "black")+guides(fill="legend")
print(p1)
#Comparing B73x009
#NoN Signifcant results
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE009, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE009)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
# T2 1 19====
#0 ====
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))



##B73 xE0011 comparsion Dseq=====
#First I will start by removing all fo the files to isolate B73 and E011
NILOTUE0011 <- subset_samples(NILOTU,Pedigree== "B73"| Pedigree=="Z031E0011")


plot_richness(NILOTUE0011, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

GP.ord <- ordinate(NILOTUE0011, "NMDS", "bray",)
p1 = plot_ordination(NILOTUE0011, GP.ord, type="Samples", color="Pedigree",  title="Fungal ITS Ordination")+geom_point(size=3)+
  labs(color="Treatment")+theme_classic()+scale_color_manual(values=c("Grey", "Black", "White"))+geom_point(shape = 1,size = 3,colour = "black")+guides(fill="legend")
print(p1)
#Comparing B73x0011
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0011, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0011)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#T2 0 6====
#0 6====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))





##B73 xE0012 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E012
NILOTUE0012 <- subset_samples(NILOTU,Pedigree== "B73"| Pedigree=="Z031E0012")


plot_richness(NILOTUE0012, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0012
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0012, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric", sfType="iterate")


res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0012)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#T2  9 19 ====
#4 19====
myTaxa = rownames(sigtab) #now we have a list of just the OTU names
sig_OTU = prune_taxa(myTaxa, NILOTUE0012) #retain just the OTU list, but now from the relativized version for plotting
sig_OTU #confirm that the phyloseq object has been subsetted
taxa_names(sig_OTU) #confirm that sig_OTU contains just the significant taxa (but now as a phyloseq object)
#[1] "410" "103" "109" "157"


title = "significant bacterial taxa"
p = plot_bar(sig_OTU, "Pedigree", fill = "Family", "Abundance", 
             title=title, facet_grid = ~Family)
p + geom_bar(aes(fill=Genus), stat="identity", position="stack")



theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))








##B73 xE0016 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E016
NILOTUE0016 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0016" )


plot_richness(NILOTUE0016, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0016 DESEQ
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0016, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0016)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#T2 2 19====
#0====


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))


##B73 xE0021 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E016
NILOTUE0021 <- subset_samples(NILOTU,Pedigree== "B73"| Pedigree=="Z031E0021")



plot_richness(NILOTUE0021, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0021
#Make sure that the diasflide changes
library("DESeq2")


diagdds = phyloseq_to_deseq2(NILOTUE0021, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0021)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#T2 2 19====
#[1] 18 19=====
#Fungi 26 13
myTaxa = rownames(sigtab) #now we have a list of just the OTU names
sig_OTU = prune_taxa(myTaxa, NILOTUE0021) #retain just the OTU list, but now from the relativized version for plotting
sig_OTU #confirm that the phyloseq object has been subsetted
taxa_names(sig_OTU) #confirm that sig_OTU contains just the significant taxa (but now as a phyloseq object)


title = "significant bacterial taxa"
p = plot_bar(sig_OTU, "Pedigree", fill = "Order", "Abundance", 
             title=title, facet_grid = ~Order)
p + geom_bar(aes(fill=Order), stat="identity", position="stack")



theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

NILOTUE0022 <- subset_samples(NILOTU, Pedigree == "Bulk"| Pedigree=="PI384071")


##B73 xE0022 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E016
NILOTUE0022 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0022")

plot_richness(NILOTUE0022, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0022
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0022, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0022)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#T2 0 6===
#0 6 ====


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))






##B73 xE0028 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E028
NILOTUE0028 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0028")



plot_richness(NILOTUE0028, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0022====
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0028, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0028)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#T2 3 19 ====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))









##B73 xE0031 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E031
NILOTUE0031 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0031")


plot_richness(NILOTUE0031, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0031, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0031)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#[1]  5 19=====
#T2 2 19-----

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))


##B73 xE0035 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E035
NILOTUE0035 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0035")



plot_richness(NILOTUE0035, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031====
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0035, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0035)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#[1]  7 19====
#T2 1 19=====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))







##B73 xE0038 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E038
NILOTUE0038 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0038")


plot_richness(NILOTUE0038, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0038, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0038)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#[1] 0 6====
#T2 1 19====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))





##B73 xE0040 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0040 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0040")


plot_richness(NILOTUE0040, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0040, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0040)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6====
#T2 1 19====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))







##B73 xE0042 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E042
NILOTUE0042 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0042")



plot_richness(NILOTUE0042, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0042, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0042)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
# 0 6====
#T2 1 19 ======
#This tells us how many OTUs are differnt

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))






##B73 xE0047 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E047
NILOTUE0047 <- subset_samples(NILOTU, Pedigree =="B73"|Pedigree == "Z031E0047")
NILOTUE0047 <- subset_samples(NILOTUT2, Pedigree =="B73"|Pedigree == "Z031E0047")
NILOTUE0047 <- subset_samples(NILOTU, Pedigree =="Z031E0059"|Pedigree == "Z031E0047")


plot_richness(NILOTUE0047, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0047=====
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0047, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0047)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6====
#1 19====
#Fungi 24 13
myTaxa = rownames(sigtab) #now we have a list of just the OTU names
sig_OTU = prune_taxa(myTaxa, NILOTUE0021) #retain just the OTU list, but now from the relativized version for plotting
sig_OTU #confirm that the phyloseq object has been subsetted
taxa_names(sig_OTU) #confirm that sig_OTU contains just the significant taxa (but now as a phyloseq object)


title = "significant bacterial taxa"
p = plot_bar(sig_OTU, "Pedigree", fill = "Order", "Abundance", 
             title=title, facet_grid = ~Genus)
p + geom_bar(aes(fill=Order), stat="identity", position="stack")




theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))






##B73 xE0050 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E050
NILOTUE0050 <- subset_samples(NILOTU, Pedigree =="B73"| Pedigree=="Z031E0050")


plot_richness(NILOTUE0050, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0050
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0050, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0050)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6=====
#T2 0 6 ====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0052 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E052
NILOTUE0052 <- subset_samples(NILOTU, Pedigree == "B73"| Pedigree=="Z031E0052")

plot_richness(NILOTUE0052, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0052, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0052)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#2 19====
#13 19=====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))







##B73 xE0054 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E054
NILOTUE0054 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree== "Z031E0054")


plot_richness(NILOTUE0054, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0054
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0054, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0054)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6=====
#T2 0 6-----

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))






##B73 xE0057 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E054
NILOTUE0057 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0057")


plot_richness(NILOTUE0057, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0057
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0057, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0057)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
# 0 6=====
#T2 0 6====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0040 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0058 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree == "Z031E0058")


plot_richness(NILOTUE0058, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0058, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0058)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6=====
#T2 06 ====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))





##B73 xE0059 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0059 <- subset_samples(NILOTU, Pedigree == "B73"| Pedigree=="Z031E0059")

plot_richness(NILOTUE0059, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0059, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0059)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0====
#T2 9 19=====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0061 comparsion Dseq==== 
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0061 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0061")


plot_richness(NILOTUE0061, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0061, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0061)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6====
#T2 1 19====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))





##B73 xE0040 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0067 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0067")


plot_richness(NILOTUE0067, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0067, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0067)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
# 4 19====
#T2 0 6----

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))





##B73 xE0068 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E068
NILOTUE0068 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0068")

plot_richness(NILOTUE0068, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes

diagdds = phyloseq_to_deseq2(NILOTUE0068, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0068)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
## 0 6====
## T2 1 10=====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))






##B73 xE0070 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E070
NILOTUE0070 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0070")


plot_richness(NILOTUE0070, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0070, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0070)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6====
#0 6====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))






##B73 xE0071 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E071
NILOTUE0071 <- subset_samples(NILOTU, Pedigree == "B73" | Pedigree=="Z031E0071")

plot_richness(NILOTUE0071, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0071, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0071)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#8 19=====
#1 19=====
#Fungal [1] 26 13==
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0074 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E074
NILOTUE0074 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0074")


plot_richness(NILOTUE0074, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0031
#Make sure that the diasflide changes
#library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0074, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0074)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6====

#1 19-----
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))






##B73 xE0507 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0507 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0507")

NILOTUE0507 <- subset_samples(NILOTUT2, Pedigree == "B73"|Pedigree=="Z031E0507")

plot_richness(NILOTUE0507, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0507
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0507, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0507)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#1 19====
#T2 1 19====
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0523 comparsion Des====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0523 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0523")


plot_richness(NILOTUE0523, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0523
#Make sure that the diasflide changes
#library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0523, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0523)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6=====
#0 6-----

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))



##B73 xE0526 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0526 <- subset_samples(NILOTU, Pedigree == "Z031E0526"|Pedigree=="B73")


plot_richness(NILOTUE0526, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0526
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0526, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0526)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#28 19=====
#1 19=====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0536 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0536 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0536")


plot_richness(NILOTUE0536, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0536
#Make sure that the diasflide changes

diagdds = phyloseq_to_deseq2(NILOTUE0536, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0536)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0====
#0 6----

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0537 comparsion Dseq =======
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0537 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0537")

plot_richness(NILOTUE0537, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0537
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0537, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0537)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0=====
#0 6====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))





##B73 xE0545 comparsion Dseq ======
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0545 <- subset_samples(NILOTU, Pedigree == "B73"| Pedigree=="Z031E0545")


plot_richness(NILOTUE0545, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0545
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0545, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0545)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6=====
#T2 1 19 ------

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))



##B73 xE0559 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0559 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0559")
NILOTUE0559 <- subset_samples(NILT1, Pedigree == "B73"|Pedigree=="Z031E0559")
NILOTUE0559 <- subset_samples(NILT2, Pedigree == "B73"|Pedigree=="Z031E0559")


plot_richness(NILOTUE0559, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0559
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0559, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0559)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#47 19======
#T2 2 19-----
myTaxa = rownames(sigtab) #now we have a list of just the OTU names
sig_OTU = prune_taxa(myTaxa, NILOTUE0559) #retain just the OTU list, but now from the relativized version for plotting
sig_OTU #confirm that the phyloseq object has been subsetted
taxa_names(sig_OTU) #confirm that sig_OTU contains just the significant taxa (but now as a phyloseq object)


title = "significant bacterial taxa"
p = plot_bar(sig_OTU, "Pedigree", fill = "Order", "Abundance", 
             title=title, facet_grid = ~Order)
p + geom_bar(aes(fill=Order), stat="identity", position="stack")


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0560 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0560 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0560")


plot_richness(NILOTUE0560, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0560
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0560, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0560)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#8 19=====
#T2 9 19----

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0566 comparsion Dseq ======
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0566 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0566")



plot_richness(NILOTUE0566, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0566
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0566, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0566)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6=====
#T2 7 19 ====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))



##B73 xE0577 comparsion Dseq ======
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0577 <- subset_samples(NILOTU, Pedigree == "B73"| Pedigree=="Z031E0577")


plot_richness(NILOTUE0577, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0577
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0577, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0577)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
# 0 6=====
#T2 2 19====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))



##B73 xE0580 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0580 <- subset_samples(NILOTU, Pedigree == "B73"| Pedigree=="Z031E0580")


plot_richness(NILOTUE0580, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0580
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0580, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0580)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#5 19=====
#0 6 -----

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))




##B73 xE0585 comparsion Dseq =====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0585 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0585")


plot_richness(NILOTUE0585, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0585
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0585, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0585)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6======
#T2 0 6====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))



##B73 xE0591 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0591 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0591")


plot_richness(NILOTUE0591, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0591
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0591, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0591)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6====
# T2 0 6 ====

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))



##B73 xE0594 comparsion Dseq ====
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0594 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="Z031E0594")


plot_richness(NILOTUE0594, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")

#Comparing B73x0594
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0594, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0594)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6=====
#T2 2 19 ----

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))



##B73 xE0PI384071 comparsion Dseq ======
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0PI384071 <- subset_samples(NILOTU, Pedigree == "B73"|Pedigree=="PI384071")


plot_richness(NILOTUE0PI384071, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")


GP.ord <- ordinate(NILOTUE0PI384071, "NMDS", "bray",)
p1 = plot_ordination(NILOTUE0PI384071, GP.ord, type="Samples", color="Pedigree",  title="16S Ordination")+geom_point(size=3)+
  labs(color="Treatment")+theme_classic()+scale_color_manual(values=c("Grey", "Black", "White"))+geom_point(shape = 1,size = 3,colour = "black")+guides(fill="legend")
print(p1)
#Comparing B73x0PI384071
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0PI384071, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.20
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0PI384071)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#
#0 ====
#T2 1 19----
##This is super weird and makes me feel like I nned an additional model to test for difference. 

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))








##B73 xE0B73xPI384071 comparsion Dseq 
#First I will start by removing all fo the files to isolate B73 and E040
NILOTUE0B73xPI384071 <- subset_samples(NILOTU, Pedigree == "B73xPI384071"| Pedigree=="B73")

plot_richness(NILOTUE0B73xPI384071, x = "Pedigree", measures = "Chao1") + geom_boxplot()+theme_classic()+theme(strip.text.x = element_blank())+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA")


GP.ord <- ordinate(NILOTUE0B73xPI384071, "NMDS", "bray",)
p1 = plot_ordination(NILOTUE0B73xPI384071, GP.ord, type="Samples", color="Pedigree",  title="Fungal ITS Ordination")+geom_point(size=3)+
  labs(color="Treatment")+theme_classic()+scale_color_manual(values=c("Grey", "Black", "White"))+geom_point(shape = 1,size = 3,colour = "black")+guides(fill="legend")
print(p1)
#Comparing B73x0B73xPI384071
#Make sure that the diasflide changes
library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE0B73xPI384071, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE0B73xPI384071)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#0 6====


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
######
###END of DESEQ2 Comparisons======

###Beginning of perANOVA analysis====
######
setwd("~/Documents/Projects/Projects/Maize Microbiome/NIL Maize Microbiome /R-Scripts/NIL CSV")
setwd("~/Documents/Projects/Research/Maize Microbiome/NIL MM2018/R-Scripts/NIL CSV")
#Taxalist location. This is probably good to have when understanding my results. 
  OTU.Taxa<-read.csv("../Matrix_eQTL/OTU.Taxa.csv")
  
  library(vegan)
  library(ggplot2)
  library(reshape2)
  library(fdrtool)
  #Somewhere in here I need to write a loop to test out to also test the FDR. It may be useful to also look at how 
#I think it would be usful to creat a taxa onkect

  NIL.16S <- read.csv("NIL16S.csv")
#NIL.16S <- read.csv("Sample_otu_table_L6.csv")

  NIL.16S<-NIL.16S[-c(39,226,176,350,351,352,353),]

#Removing the Teosinte and B73 treatments
  NIL.16S<-subset(NIL.16S, Pedigree!="PI384071"|Pedigree!="B73xPI384071"	)

  NIL.16S.spec <- NIL.16S[,8:14031]
  NIL.16S.treatment <- NIL.16S[,1:7]

### Density Plot of all of the OTUs in the NIL study
  #Melt tansforms data to two axis| log stransfrom so they are on the same scale. 
  x<-melt(log(NIL.16S.spec[,1:500]))
  ggplot(data=x, aes(x=value,stat(count),fill=variable))+geom_density( position="stack")+stat_density()+theme(legend.position = "none")

#Anova Models 
  adonis(NIL.16S.spec ~ Time, NIL.16S.treatment)
  adonis(NIL.16S.spec ~ Pedigree, NIL.16S.treatment)
  adonis(NIL.16S.spec ~ Range, NIL.16S.treatment)
  adonis(NIL.16S.spec ~ Row, NIL.16S.treatment)
  
  adonis(formula= NIL.16S.spec~ Time+Pedigree,NIL.16S.treatment)
  adonis(formula= NIL.16S.spec~ Pedigree+Time+RangeRow,strata=NIL.16S.treatment$RangeRow,NIL.16S.treatment)
  adonis(formula= NIL.16S.spec~ Pedigree+Time+Range,strata=NIL.16S.treatment$Range,NIL.16S.treatment)
  adonis(formula= NIL.16S.spec~ Pedigree+Time+Row+Range,strata=NIL.16S.treatment$Row,NIL.16S.treatment)
  adonis(formula= NIL.16S.spec~ Pedigree+Time+Row+Range,strata=NIL.16S.treatment$Range,NIL.16S.treatment)
    
#Here I'm running the model and time and genotype have an interaction. 
#I'm kinda confused by whats dropped by the model 
adonis(formula= NIL.16S.spec~ Pedigree*Time+Row+Range+RangeRow,strata=NIL.16S.treatment$Range,NIL.16S.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Time        1     0.615 0.61482  7.6916 0.01885  0.001 ***
#   Pedigree   45     4.401 0.09781  1.2236 0.13496  0.001 ***
#   RangeRow  134    13.847 0.10333  1.2927 0.42459  0.001 ***
#   Residuals 172    13.749 0.07993         0.42159           
# Total     352    32.612                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Here we see no digficant effect of the genotype time intereaction
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Pedigree       43    4.0424 0.09401  1.1893 0.12835  0.001 ***
#   Time            1    0.6292 0.62916  7.9594 0.01998  0.001 ***
#   Row            21    2.2841 0.10876  1.3760 0.07252  0.001 ***
#   Range           7    1.7288 0.24697  3.1244 0.05489  0.001 ***
#   RangeRow      104    9.4227 0.09060  1.1462 0.29917  0.001 ***
#   Pedigree:Time  43    3.4297 0.07976  1.0090 0.10889  0.391    
#   Residuals     126    9.9597 0.07905         0.31622           
# Total         345   31.4966                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Loop runs indivdual OTUs to everyfactor in the model====
NIL.16S.spec <- NIL.16S[,8:14031]

NIL.16S.spec <- NIL.16S[,8:20]

pValuesOTU<-list()
for (numbers in c(1:length(NIL.16S.spec))){
  model<-aov(NIL.16S.spec[,numbers] ~ Pedigree*Time+Range+Row+RangeRow,NIL.16S.treatment)
  print("ANOVA for Vx")
  print(numbers)
  Ps<-summary(model)
  pValuesOTU[[paste0("", numbers)]]<-Ps}



#587 was a blank failed model. Need to remove| This code below worked once 

failedRuns<-array()

for (i in c(1:length(pValuesOTU))){
  if (pValuesOTU[[paste0(i)]][[1]][["Mean Sq"]][1]<=0 & pValuesOTU[[paste0(i)]][[1]][["Mean Sq"]][4]<=0)
  {{print(paste0(i))
    failedRuns[i]<- paste0(i,"")
   }}}

failedRuns<-na.omit(failedRuns)
failedRuns<-as.numeric(failedRuns)

#Now Im going to use this array to remove otus

#Now we will rerun it and confirm that we are removing everything. 
NIL.16S.spec<-NIL.16S.spec[,-c(failedRuns)]
pValuesOTU<-list()
for (numbers in c(1:length(NIL.16S.spec))){
  model<-aov(NIL.16S.spec[,numbers] ~ Pedigree*Time+Range+Row+RangeRow,NIL.16S.treatment)
  print("ANOVA for Vx")
  print(numbers)
  Ps<-summary(model)
  pValuesOTU[[paste0("", numbers)]]<-Ps}
##I Can probably add tukey to this loop to get additional infromation === 


#This checks the table
for (i in c(1:length(pValuesOTU))){
  if (pValuesOTU[[paste0(i)]][[1]][["Mean Sq"]][1]<=0 & pValuesOTU[[paste0(i)]][[1]][["Mean Sq"]][4]<=0)
  {(print(paste0(i)))}
}

##Lop tells me which models had signifcant effets||| NEED TO ADD A FDR Correction. 
#Here instead I'm just going to do a bonforni correction. The BFC just increases the signficance threshold to make correlations hardwere.
#The BFC is alpha/number to tests to get a new alpha. 
sigOTU<-array()
#With BFC Power loss is signifcant between the two tests. 
for (i in c(1:length(pValuesOTU))){
  if (pValuesOTU[[paste0(i)]][[1]][["Pr(>F)"]][1]<(0.05/length(pValuesOTU)) || pValuesOTU[[paste0(i)]][[1]][["Pr(>F)"]][6] < (0.05/length(pValuesOTU)) )
    {print(paste0(i))
    sigOTU[i]<- paste0(i,"")
  } 

}


#Without BFC
for (i in c(1:length(pValuesOTU))){
  if (pValuesOTU[[paste0(i)]][[1]][["Pr(>F)"]][1]<0.05 || pValuesOTU[[paste0(i)]][[1]][["Pr(>F)"]][6] < 0.05 )
  {print(paste0(i))
    sigOTU[i]<- paste0(i,"")
  } 
  
}
#This removes all of the non-signifcant treatment OTUs. This allows for proper decomposion of the MC
sigOTU<-na.omit(sigOTU)
sigOTU<-as.numeric(sigOTU)
sigOTU ## 107 OTUs That are infleucned by hose genetics with BFC | 1653 without BFC
pValuesOTUSig<-pValuesOTU[sigOTU]

#This prints the taxa that are significant 
OTU.Taxa[sigOTU,]
#This makes a dataframe just ot be tested for signifcance
NIL.sig.16S<-NIL.16S.spec[,c(sigOTU)]
sigA<-as.data.frame(NIL.sig.16S)

#Percentage of the MC the is reponding to treatment factor
#12% of the total number of OTUs ||| 0.007% of th MC is Responding to G
length(pValuesOTUSig)/length(pValuesOTU)
#Gives modest levels of diversity in these sentios|| I think what is most important from this method is to just pull out ant test the OTUs. 
#It also tells us how many OTUs may be signfincantly correlated. wit the model 
adonis(formula= NIL.sig.16S~ Pedigree*Time+Row+Range+RangeRow,strata=NIL.16S.treatment$RangeRow,NIL.16S.treatment)
#Large amounts of this diversity doesnt seem to be clear. 

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree       43    2.6841 0.06242  1.7557 0.15269  0.001 ***
#   Time            1    0.4577 0.45771 12.8740 0.02604  0.001 ***
#   Row            21    1.4061 0.06696  1.8832 0.07998  0.001 ***
#   Range           7    1.1301 0.16145  4.5410 0.06429  0.001 ***
#   RangeRow      104    5.5158 0.05304  1.4918 0.31377  0.001 ***
#   Pedigree:Time  43    1.9057 0.04432  1.2465 0.10840  0.001 ***
#   Residuals     126    4.4797 0.03555         0.25483           
# Total         345   17.5792                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Results from running the Model on all OTUS

####This is the model to just pull out signficant values from a individual row
RangeRowOTU<-array()
for (i in c(1:length(pValuesOTU))){
  if (pValuesOTU[[paste0(i)]][[1]][["Pr(>F)"]][5]<0.05) 
  {print(paste0(i))
    RangeRowOTU[i]<- paste0(i,"")
  } 
  else {print("boo!")}
}
RangeRowOTU<-na.omit(RangeRowOTU)
RangeRowOTU<-as.numeric(RangeRowOTU)
pRangeRowOTUSig<-pValuesOTU[RangeRowOTU]


##Here I want to add a TUKEY addition to this model 
#
library(tidyverse)

##This look test every OTU in the MC using the Model that I developed to be the best. 
##Once the model is done then it runs a Tukey HSD on the model and collects the The information about which OTU is causing that
pValuesOTU<-list()
TukeyValue<-list()
for (numbers in c(1:length(NIL.16S.spec))){
  model<-aov(NIL.16S.spec[,numbers] ~ Pedigree*Time+Range+Row+RangeRow,NIL.16S.treatment)
  print("ANOVA for Vx")
  print(numbers)
  Ps<-summary(model)
  Ts<-TukeyHSD(model,which="Pedigree",ordered = TRUE)
  result_table <- Ts$Pedigree %>% as_tibble(rownames = "Pedigree")
  pValuesOTU[[paste0("", numbers)]]<-Ps
  TukeyValue[[paste0("", numbers)]]<-result_table %>% filter(`p adj` < 0.05)%>%filter(str_detect(Pedigree, "B73$"))}

library(rlist)
##We run the model and then use the data
SigTukey<-TukeyValue[sigOTU]
SigTukeyrmNA<-list.filter(SigTukey, Pedigree!='0')

TEST<-as.data.frame(unlist(SigTukeyrmNA, recursive = FALSE))

#write.csv(TEST,"SigTurkeyLoopResults.csv", row.names = FALSE)

#Run an R package that looks at the interaction effect. 
#
#Independent Component analysis is a computatioal method for separating a multivaraite signal inot addive subcomponents. 
#This is done by assuming that the subcomponents are non-Gaussian signals and that they are statisitcally independent from each other. #This assumes that signal is non-normal-- this is probably what I would assume MC to be 
#Statisitcal indpendence is a fundumental notion in probablity theory. Two events are stastically indpendent or stochastically indepdent if the occurance of one does not affect the proabblity of occurance of the other. 
#ICA is a type of blind source seperation and is a method of decomposing multivarate signals. 
#ICA was developed to solve the famous "cocktail party" problem... and is used to decode all of the voice you hear from a cocktail party. 

#Here I am using ICA to decode the signals/abundances of the micorbial community to 1. minimize mutual infromaiton across different OTUS 2. Mazimizae non-gaussian features. 
#What interesting about the ICA model is that you can set of the number of indepdendent components it needs to find. 

###ICA Analysis of the MC=====
#this calls and runs all of LibICA commands
source('~/Downloads/LibICA.r', echo=TRUE)
library(fastICA)


#Remove garbage data points| This filter removes anything not present in 5% of samples
filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.10 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
NIL.16S.spec2 = filter(NIL.16S.spec)

#Here I want to just test the second time point
NIL.16S.spec2 = filter(NILT2.spec)
NIL.16S.treatment<-NILT2.treatment


tGenes=t(log1p(NIL.16S.spec2))
tGenes=t((NIL.16S.spec))

#Here I am testing how alterations to the mixing Matrix impacts the results. 
#Im concered that the code I was useing was currint
#Between the two methods it seems as if A and W are the same thing. 
#NILICAResults<-runICA(tGenes, ncomp=2)

#Setting the number of number of independent axis|| 200 indepedent compontents may be where I'm computationally limited at this computer. 
Dems<-50

NILICAResults<-fastICA(tGenes, n.comp=Dems, alg.typ = "deflation", fun = "logcosh", alpha = 1,
        method = "C", row.norm = FALSE, maxit = 200,
        tol = 0.0001, verbose = TRUE)

#ICAGenes<-getGenesICA(NILICAResults)

#Pulling out units in the model
X<-as.data.frame(t(NILICAResults$X))
S<-as.data.frame((NILICAResults$S))
#M<-as.data.frame(abs(t(NILICAResults$M)))
#PERNOVA only run on abs
A<-as.data.frame(abs(t(NILICAResults$A)))
#A<-as.data.frame((t(NILICAResults$A)))

rownames(A)<-t(NIL.16S.treatment[,1])
heatmap(as.matrix(A),main="A Martix",xlab="IC",ylab="Samples")

#Testing signficantce
adonis(M~ Pedigree+Time+Range+Row, strata=NIL.16S.treatment$RangeRow, NIL.16S.treatment)
adonis(M~ Pedigree+Time+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)
adonis(A~ Pedigree+Time+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)
adonis(A~ Pedigree*Time+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)


heatmap(NILICAResults$S,main="S Martix",xlab="IC",ylab="OTUs")
heatmap(NILICAResults$A,main="A Martix",xlab="Samples",ylab="IC")

#Here I'm looking at how the ICA componetns vary in explaing the stats 
##Testing the signifcance the ICAs and tells us the stastical realtionship with the ICA
pValues<-list()
for (numbers in c(1:Dems)){
model<-aov(A[,numbers] ~ Pedigree*Time+Range+Row+RangeRow,NIL.16S.treatment)
print("ANOVA for Vx")
print(numbers)
Ps<-summary(model)
pValues[[paste0("", numbers)]]<-Ps}

##This loop is for a single time point
for (numbers in c(1:Dems)){
  model<-aov(A[,numbers] ~ Pedigree+Range+Row,NIL.16S.treatment)
  print("ANOVA for Vx")
  print(numbers)
  Ps<-summary(model)
  pValues[[paste0("", numbers)]]<-Ps}
#for (i in c(1:50)){print(i)
  #print(pValues[[paste0(i)]][[1]][["Pr(>F)"]])}

#Loop to determine signifcant.Trues are sig not trues are nto
# #for (i in c(1:50)){
#   if (pValues[[paste0(i)]][[1]][["Pr(>F)"]][1]<.05)
#   {print(paste0(i,"-",TRUE))}
#   else{print(paste0(FALSE,i,"-"))}}

#Loop to determine signifcant.Trues are sig not true

#Here we are just printing the list
# This filtering step focuses on the MC that we are important to the genotype time effect
# Here & means both must be true while || means ANDalso

sigI<-array()
for (i in c(1:Dems)){
  if (pValues[[paste0(i)]][[1]][["Pr(>F)"]][1]<0.05 || pValues[[paste0(i)]][[1]][["Pr(>F)"]][6]<0.05)
  {print(paste0(i))
    sigI[i]<- paste0(i,"")
  }
}


##Loop for a single timepoint
sigI<-array()
for (i in c(1:Dems)){
  if (pValues[[paste0(i)]][[1]][["Pr(>F)"]][1]<0.05 )
  {print(paste0(i))
    sigI[i]<- paste0(i,"")
  }
}
#Timepoint 2 Deminsions.
# [1] "2"
# [1] "16"
# [1] "23"
# [1] "42"


#This filters the data by those that respond to just the location
#This looks for arrays that are only signifcant for the area space
# sigI<-array()
# for (i in c(1:Dems)){
#        if (pValues[[paste0(i)]][[1]][["Pr(>F)"]][5]<0.05 & pValues[[paste0(i)]][[1]][["Pr(>F)"]][1]>0.10)
#          {print(paste0(i))
#          sigI[i]<- paste0(i,"")
#        }
#   }

#Altering the [#]<.0.5 will tell you which factor should be checked for signifcance
sigI<-na.omit(sigI)
sigI<-as.numeric(sigI)

#This makes a dataframe just ot be tested for signifcance
sigA<-A[,c(sigI)]
sigA<-as.data.frame(sigA)
#This combination of steps allows for the denoising of the data

##Here we used the fucntion developed by others to atch the NIL results
##By doing this it will be nicce to see how the NIL preform. 
OTUsICs<-getGenesICA(NILICAResults, alpha=0.005)
OTUsig<-OTUsICs[sigI]

#Using my taxa informaiton i would like to add the additional meta infromation to make the values easier to understand and interperute

#This is a quick and dirty was of paritioning the ICs in there respective OTUs
#By looking at this data it also seems as if there are varation differences across OTU.
#Maybe incresing the stringiscy may be a better way at looking at the MC. 

#From here I can make a porpotional figure based on counts of OTUs and repeats to better understand how the data looks. Probably something from tidyverse. 
library(tidyverse)
OTUjoined<-OTUsICs[["ic40"]][["neg"]]%>%full_join(OTU.Taxa, by = c("genes"="X.OTU.ID"))
OTUjoined<-na.omit(OTUjoined)

#This is how you print all of the values
# sigA<-data.frame()
# for(i in 1:length(sigI)){
#   sigA<-rbind(c(A[,(sigI[i])]))}

#When we break up the or
#For now I'll do this manually 1,3,9,14,15
#I want to automate this part.. basically I need to create a data frame with just the sig values

#Here I make a list I'm not sure how to quickly search it. Heres how I look at a singel value 
#pValues[["1"]][[1]][["Pr(>F)"]]

#This could probably be made better, but I dont have time to do this rn
#This loop tests compoenets of the model and tells me if they are sigifcant 
#I run ICA with 50 components to determie if we could fine the axis that where most important- The noise reduction method helps when looking at the whole community at once
#The filtering has it limitations but results in modest increases
adonis(sigA~ Pedigree+Time+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)
#This is the model that I think best chararctize the data
adonis(sigA~ Pedigree*Time+Range+Row+RangeRow, strata=NIL.16S.treatment$RangeRow,NIL.16S.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Pedigree       43    2.4719 0.05749  2.9273 0.17524  0.001 ***
#   Time            1    0.4844 0.48436 24.6647 0.03434  0.001 ***
#   Range           7    1.1784 0.16834  8.5721 0.08354  0.001 ***
#   Row            21    1.6553 0.07882  4.0138 0.11735  0.001 ***
#   RangeRow      104    4.8059 0.04621  2.3532 0.34071  0.001 ***
#   Pedigree:Time  43    1.0356 0.02408  1.2264 0.07341  0.043 *  
#   Residuals     126    2.4744 0.01964         0.17541           
# Total         345   14.1058                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Single Timepoint
adonis(sigA~ Pedigree+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# Pedigree   41    4.5376 0.110672  1.7125 0.32524  0.003 **
# Range       7    1.7310 0.247291  3.8264 0.12408  0.147   
# Row        21    1.6077 0.076557  1.1846 0.11524  0.184   
# Residuals  94    6.0750 0.064628         0.43544          
# Total     163   13.9513                  1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(A[,c(3,  4,  9,  11, 18, 20, 28, 32, 41, 48, 49)]~ Pedigree+Time+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)


adonis(A[,c(3,  4,  9,  11, 18, 20, 28, 32, 41, 48, 49)]~ Pedigree+Time+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)

###This is superid 
adonis(A[,c( 9,  22, 26, 27, 28, 29, 30 )]~ Pedigree*Time+Range+Row+RangeRow, strata=NIL.16S.treatment$RangeRow,NIL.16S.treatment)
#About 45% of variance is explained by space. ~17% from the  and 2% change in time. The ridsiduals are about 25%
#the block ID is the strongest factor~ 30% i.e. Lots of the varaince is not changing in time. 


#This helps alot when we are looking at the full micorbila communitty
#adonis(A[,c(1,  3,  8,  15)]~ Pedigree+Time+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)

#I'd like to make a simple look to look at the signficant IC factors and put them directly into the anvoa model
adonis(A[,c( 1,3,9,14,15)]~ Pedigree*Time+Range+Row, strata=NIL.16S.treatment$Range,NIL.16S.treatment)
#There doesnt seem to be an interaction effect regardless of the type of model that is being run. 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Pedigree   44    1.6635 0.03781  1.7029 0.17137  0.001 ***
#   Time        1    0.4000 0.39998 18.0157 0.04121  0.001 ***
#   Range      10    0.8449 0.08449  3.8058 0.08705  0.026 *  
#   Row        21    0.6709 0.03195  1.4390 0.06912  0.007 ** 
#   Residuals 276    6.1276 0.02220         0.63126           
#   Total     352    9.7070                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###This is how you get the corrlated compoenets from the values 
OTUsICs<-getGenesICA(NILICAResults)
#
##the introgression is important
##These are the independ compoenets that aline with the NIL treatment. 
#THis DIDNT work
heatmap(NILICAResults$S[,c(3,5,6,8,11,12,14,15)],main="S Martix",xlab="IC",ylab="OTUs")

##getGenes didnt work because lack of name is colom

colnames(NILICAResults$S)<-(NIL.16S.treatment[,1])

#At an ICA of 6 we can see the the signifcance of pedigress is the highest. There does seems to be some funck stuff going on with the perMANVOA
#These  predictions and results can vey its kinda weird. THey can go up and down dending on 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#   Pedigree   44    0.9523 0.021642  1.6895 0.16008  0.001 ***
#   Time        1    0.1080 0.107965  8.4282 0.01815  0.001 ***
#   Range      10    0.7657 0.076570  5.9773 0.12871  0.001 ***
#   Row        21    0.5873 0.027968  2.1833 0.09873  0.001 ***
#   Residuals 276    3.5356 0.012810         0.59433           
# Total     352    5.9488                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#After testing a number of components 6 seemed to be the most important factor 

plot(V10~V15, data=S, col=c(rainbow(1000)[NIL.16S.treatment$RangeRow]))


#Run an R package that looks at the interaction effect. 

##Simper Analysis-----
#Filtering Fuction allows you to run analysis 
NIL.16S <- read.csv("NIL16S.csv")

# NIL.16S.1 <- subset(NIL.16S, Pedigree == "B73"|Pedigree=="Z031E0057"|Pedigree=="Z031E0021")

NIL.16S.spec <- NIL.16S[,8:14031]
NIL.16S.treatment <- NIL.16S[,1:7]

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.10 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.16S.spec)


rrna.rel <- decostand(X.filtered, method="hellinger")

object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
test<-filter(object,names(object)$"^B73")
filter(object, str_detect(names(object), "B73_Z031E0021"))
bro<-data.frame(unclass(object))
write.csv(bro, file = "SimperResults.CSV")
#since time isnt a factor I will look at it a little different
##Here Im going to do simper across every B73-introgression comparison ====

##Z031E009====
NILOTU<- subset(NIL.16S, Pedigree == "Z031E0009" |Pedigree == "B73" )

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.05 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )

bro<-data.frame(unclass(object))
bro<-dplyr::filter(bro, B73_Z031E0009.p<0.05)
write.csv(bro, file = "B73_Z031E0009.CSV")
##Z031E0011====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0011")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.05 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))
bro<-dplyr::filter(bro, B73_Z031E0011.p<0.05)
write.csv(bro, file = "B73_Z031E0011.CSV")
##Z031E0012====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0012")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))
bro<-dplyr::filter(bro, B73_Z031E0012.p<0.05)
write.csv(bro, file = "B73_Z031E0012.CSV")
##Z031E0016====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0016")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))
bro<-dplyr::filter(bro, B73_Z031E0016.p<0.05)
write.csv(bro, file = "B73_Z031E0016.CSV")
##Z031E0021====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0021")

##Fungal coding
NIL.ITS <- read.csv("otu_table_ITS_2302.csv")
###Below is to remove specifc samples from dataset
NILOTU <- subset(NIL.ITS,Pedigree== "B73"| Pedigree=="Z031E0021")
NIL.ITS.spec <- NILOTU[,8:2253]
NIL.ITS.treatment <- NILOTU[,1:7]

##Prokaytoic coding
NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
object<-simper(rrna.rel, NIL.ITS.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )

## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0021.p<0.05)
write.csv(bro, file = "B73_Z031E0021.CSV")

##Z031E0022====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0022")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0022.p<0.05)
write.csv(bro, file = "B73_Z031E0022.CSV")

##Z031E0028====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0028")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0028.p<0.05)
write.csv(bro, file = "B73_Z031E0028.CSV")

##Z031E0031====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0031")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, Z031E0031_B73.p<0.05)
write.csv(bro, file = "B73_Z031E0031.CSV")

##Z031E0035====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0035")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0035.p<0.05)
write.csv(bro, file = "B73_Z031E0035.CSV")
##Z031E0038====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0038")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0038.p<0.05)
write.csv(bro, file = "B73_Z031E0038.CSV")

##Z031E0040====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0040")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]


X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0040.p<0.05)
write.csv(bro, file = "B73_Z031E0040.CSV")

##Z031E0042====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0042")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]


X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0042.p<0.05)
write.csv(bro, file = "B73_Z031E0042.CSV")

##Z031E0047====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0047")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]


X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0047.p<0.05)
write.csv(bro, file = "B73_Z031E0047.CSV")


##Z031E0050====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0050")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0050.p<0.05)
write.csv(bro, file = "B73_Z031E0050.CSV")

##Z031E0052====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0052")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0052.p<0.05)
write.csv(bro, file = "B73_Z031E0052.CSV")

##Z031E0054====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0054")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0054.p<0.05)
write.csv(bro, file = "B73_Z031E0054.CSV")

##Z031E0057====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0057")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0057.p<0.05)
write.csv(bro, file = "B73_Z031E0057.CSV")

##Z031E0058====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0058")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0058.p<0.05)
write.csv(bro, file = "B73_Z031E0058.CSV")

##Z031E0059====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0059")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0059.p<0.05)
write.csv(bro, file = "B73_Z031E0059.CSV")

##Z031E0061====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0061")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0061.p<0.05)
write.csv(bro, file = "B73_Z031E0061.CSV")


##Z031E0067====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0067")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0067.p<0.05)
write.csv(bro, file = "B73_Z031E0067.CSV")

##Z031E0068====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0068")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, Z031E0068_B73.p<0.05)
write.csv(bro, file = "B73_Z031E0068.CSV")


##Z031E0070====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0070")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0070.p<0.05)
write.csv(bro, file = "B73_Z031E0070.CSV")


##Z031E0071====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0071")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0071.p<0.05)
write.csv(bro, file = "B73_Z031E0071.CSV")

   
##Z031E0074====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0074")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, Z031E0074_B73.p<0.05)
write.csv(bro, file = "B73_Z031E0074.CSV")

##Z031E0507====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0507")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0507.p<0.05)
write.csv(bro, file = "B73_Z031E0507.CSV")

##Z031E0523====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0523")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, Z031E0523_B73.p<0.05)
write.csv(bro, file = "B73_Z031E0523.CSV")

##Z031E0526====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0526")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0526.p<0.05)
write.csv(bro, file = "B73_Z031E0526.CSV")

##Z031E0536====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0536")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0536.p<0.05)
write.csv(bro, file = "B73_Z031E0536.CSV")

##Z031E0537====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0537")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0537.p<0.05)
write.csv(bro, file = "B73_Z031E0537.CSV")

##Z031E0545====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0545")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0545.p<0.05)
write.csv(bro, file = "B73_Z031E0545.CSV")

##Z031E0556====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0556")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0556.p<0.05)
write.csv(bro, file = "B73_Z031E0556.CSV")


##Z031E0559====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0559")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0559.p<0.05)
write.csv(bro, file = "B73_Z031E0559.CSV")

##Z031E0560====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0550")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0560.p<0.05)
write.csv(bro, file = "B73_Z031E0560.CSV")

##Z031E0566====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0556")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0566.p<0.05)
write.csv(bro, file = "B73_Z031E0566.CSV")

##Z031E0577====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0577")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0577.p<0.05)
write.csv(bro, file = "B73_Z031E0577.CSV")

##Z031E0578====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0578")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0578.p<0.05)
write.csv(bro, file = "B73_Z031E0578.CSV")

##Z031E0580====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0580")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0580.p<0.05)
write.csv(bro, file = "B73_Z031E0580.CSV")

##Z031E0585====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0585")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0585.p<0.05)
write.csv(bro, file = "B73_Z031E0585.CSV")

##Z031E0591====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0591")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0591.p<0.05)
write.csv(bro, file = "B73_Z031E0591.CSV")

##Z031E0594====
NILOTU <- subset(NIL.16S,Pedigree== "B73"| Pedigree=="Z031E0594")

NIL.16S.spec <- NILOTU[,8:14031]
NIL.16S.treatment <- NILOTU[,1:7]

X.filtered = filter(NIL.16S.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
#Simper
object<-simper(rrna.rel, NIL.16S.treatment$Pedigree, permutations = 999, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(object))

bro<-dplyr::filter(bro, B73_Z031E0594.p<0.05)
write.csv(bro, file = "B73_Z031E0594.CSV")
#Im goign to break up the data into the two different time points because time has an effect on the community compostion
NIL.T1=subset(NIL.16S,Time== "T1")
NIL.T2=subset(NIL.16S,Time== "T2")

NIL.16S.spec <- NIL.T1[,7:14030]
NIL.16S.treatment <- NIL.T1[,1:6]

adonis(NIL.16S.spec ~ Pedigree, NIL.16S.treatment)
anosim(NIL.16S.spec, NIL.16S.treatment$Pedigree)

NIL.16S.spec <- NIL.T2[,7:14030]
NIL.16S.treatment <- NIL.T2[,1:6]

adonis(NIL.16S.spec ~ Pedigree, NIL.16S.treatment)


###Subsetting Teosinte, B73 and Hybrid and soil 
#########

###Subsetting code
#######

##load Data
All = read.csv(file ="NIL16S.csv")
row.names(All) <- All$X
X=rbind(
  subset(All, Pedigree=="B73xPI384071"),
  subset(All, Pedigree=="B73"),
  subset(All, Pedigree=="PI384071"))

NIL.Set.spec <- X[,8:14031]
NIL.Set.treatment <- X[,1:7]

adonis(NIL.Set.spec~ Pedigree,NIL.Set.treatment)
adonis(formula= NIL.Set.spec~ Pedigree+Time+RangeRow,strata=NIL.Set.treatment$Time,NIL.Set.treatment)
anosim(NIL.Set.spec, NIL.Set.treatment$Pedigree)

COL.treat <- c("purple","red","Black","blue","blue")
MDS <- metaMDS(NIL.Set.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
points(MDS$points,  
       col=COL.treat[NIL.Set.treatment$Pedigree],
       pch=18, cex=2)


object<-simper(NIL.Set.spec, NIL.Set.treatment$Pedigree, permutations = 200, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
bro<-summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )
bro<-data.frame(unclass(bro))
write.csv(bro, file = "SimperResultsTeosinteMaize.CSV")



##Phyloseq Fungi Stuff

FungibiomOTU = import_biom('Phyloseq/otu_table_ITS_2302.biom',treefilename = 'Phyloseq/rep_seq_aligned_pfiltered.tre', 
                           parseFunction = parse_taxonomy_greengenes)

NILMap = import_qiime_sample_data("NIL2018MappingDNA.txt")

#Combining Biom and Mapping file to Make the phyloseq file
FunNILOTU = merge_phyloseq(FungibiomOTU, NILMap)

#Average by treatment
OTU_merge = merge_samples(FunNILOTU, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTU_merge)$Pedigree <- levels(sample_data(FunNILOTU)$Pedigree)
OTU_rel_merge = transform_sample_counts(OTU_merge, function(x) 100 * x/sum(x) )


#Attemping to subset to 100 fungi #I think this worked
TopFungiOTUs = names(sort(taxa_sums(FunNILOTU), TRUE) [1:100])
TopFungi = prune_taxa(TopFungiOTUs, OTU_rel_merge)

#This function is just to make sure the subsetting worked
GP.ord <- ordinate(TopFungi, "NMDS", "bray",)
p1 = plot_ordination(TopFungi, GP.ord, type="Samples", color="Pedigree", title=" Ordination")+
  geom_point(size=3)+geom_text(mapping = aes(label = Pedigree), size = 10, vjust = 1.5) 
  print(p1)

#Test Example
NILOTUE009 <- subset_samples(NILOTU,Pedigree== "B73"| Pedigree=="Z031E009")
  

library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUE009, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUE009)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#Based on these results it feels as if the fungi should be rerun in the Loci assignment matrix to see if we can find an effect. 

##########




######Vegan Code: Stastical Analysis ======
#It's clear that a portion of the microbial communtiy is infact driven by the genetic introfgression. 
#This can be seen in our model 
#These models need to be rerun without the teosinte and bulk soil factors this will allow us to definaty say that there is an effect 
# I want to use this infromation and grow upon it. We know that the genetic introgressions have an efffet
#This effect is ofcourse on the single loci level. I.E. altering genetic introgession within a culitvar can alter the micorbiome. 
library(vegan)

NIL.16S <- read.csv("NIL16S.csv")
#NIL.16S <- read.csv("Sample_otu_table_L6.csv")
#Remove outliers
NIL.16S<-NIL.16S[-c(39,226,176,350,351,352,353),]
NIL.16S<-subset(NIL.16S, Pedigree!="Bulk")
NIL.16S<-subset(NIL.16S, Pedigree!="B73xPI384071")
NIL.16S<-subset(NIL.16S, Pedigree!="PI384071")
NIL.16S<-subset(NIL.16S,Pedigree!="Bulk")
NIL.16S<-subset(NIL.16S,Pedigree!="Bulk ")

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.005 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}


NILT1<-subset(NIL.16S, Time=="T1")
NILT2<-subset(NIL.16S, Time=="T2")
NILF1TeoB73<-subset(NIL.16S, Pedigree=="B73"|Pedigree=="B73xPI384071"|Pedigree=="PI384071")
  
NIL.16S<-subset(NIL.16S, Pedigree!="B73xPI384071"|Pedigree!="PI384071")
NIL.16S.spec <- NIL.16S[,8:14031]
NIL.16S.treatment <- NIL.16S[,1:7]

##Flitered and Transfromed Data
NIL.16S<-filter(NIL.16S)
NIL.16S.spec <- NIL.16S[,8:11195]
NIL.16S.treatment <- NIL.16S[,1:7]
NIL.16S.spec <- decostand(NIL.16S.spec, method="hellinger")


adonis(NIL.16S.spec ~ Time, NIL.16S.treatment)
adonis(NIL.16S.spec ~ Pedigree, NIL.16S.treatment)
adonis(NIL.16S.spec ~ Range, NIL.16S.treatment)
adonis(NIL.16S.spec ~ Row, NIL.16S.treatment)

adonis(formula= NIL.16S.spec~ Time+Pedigree,NIL.16S.treatment)
adonis(formula= NIL.16S.spec~ Pedigree+Time+RangeRow,strata=NIL.16S.treatment$RangeRow,NIL.16S.treatment)
##12 to 13 % of the MC is related to changes in the MC
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   44     4.250 0.09660  1.2085 0.13033  0.001 ***
#   Time        1     0.621 0.62058  7.7636 0.01903  0.001 ***
#   RangeRow  135    13.992 0.10364  1.2966 0.42905  0.001 ***
#   Residuals 172    13.749 0.07993         0.42159           
# Total     352    32.612                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
####
adonis(formula= NIL.16S.spec~ Pedigree+Time+Range+Row+RangeRow,strata=NIL.16S.treatment$RangeRow,NIL.16S.treatment)

###This model below is just testing the effects of the thr introgression on the microbiome 

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   43    2.2624 0.052615  1.1841 0.12890  0.001 ***
#   Time        1    0.2646 0.264557  5.9541 0.01507  0.001 ***
#   Range       7    0.8701 0.124296  2.7974 0.04957  0.001 ***
#   Row        21    1.3147 0.062604  1.4089 0.07490  0.001 ***
#   RangeRow  104    5.3311 0.051260  1.1537 0.30373  0.001 ***
#   Residuals 169    7.5092 0.044433         0.42782           
# Total     345   17.5520                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



adonis(formula= NIL.16S.spec~ Pedigree+Time+Range,strata=NIL.16S.treatment$Range,NIL.16S.treatment)
adonis(formula= NIL.16S.spec~ Pedigree+Time+Row+Range,strata=NIL.16S.treatment$Row,NIL.16S.treatment)
adonis(formula= NIL.16S.spec~ Pedigree+Time+Row+Range,strata=NIL.16S.treatment$Range,NIL.16S.treatment)
#adonis(formula= NIL.16S.spec~ Pedigree*Time*Row*Range,NIL.16S.treatment)
#Stasitcal analysis showing that 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   44     4.250 0.09660  1.1459 0.13033  0.006 ** 
#   Time        1     0.621 0.62058  7.3616 0.01903  0.001 ***
#   Row        24     2.743 0.11428  1.3557 0.08410  0.001 ***
# Range       7     1.731 0.24735  2.9341 0.05309  0.003 **
#   Residuals 276    23.267 0.08430         0.71344           
# Total     352    32.612                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Large amount of the variance is being produced by the field variable 

#
#I want to run these seqeunces as ASVs, but the more I think about it the more unnessary it seems. At least from the QTL analysis there are OTUs that are correlated to the regions
#Might need to incorrporate tidyverse code here for doing this kind of stuff
adonis(formula= NIL.16S.spec~ Pedigree+Time+Row+Range+RangeRow,strata=NIL.16S.treatment$RangeRow,NIL.16S.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   41    3.3834 0.08252  1.1609 0.12836  0.001 ***
#   Time        1    0.3749 0.37490  5.2740 0.01422  0.001 ***
#   RangeRow  126   11.1567 0.08855  1.2456 0.42325  0.001 ***
#   Residuals 161   11.4445 0.07108         0.43417           
# Total     329   26.3595                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(formula= NIL.16S.spec~ Pedigree*Time+Row+Range+RangeRow,strata=NIL.16S.treatment$RangeRow,NIL.16S.treatment)

# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree       41    3.3834 0.08252  1.1597 0.12836  0.001 ***
#   Time            1    0.3749 0.37490  5.2687 0.01422  0.001 ***
#   Row            21    1.9677 0.09370  1.3168 0.07465  0.001 ***
#   Range           7    1.2873 0.18389  2.5844 0.04883  0.001 ***
#   RangeRow       98    7.9018 0.08063  1.1332 0.29977  0.001 ***
#   Pedigree:Time  41    2.9058 0.07087  0.9960 0.11024  0.593    
# Residuals     120    8.5387 0.07116         0.32393           
# Total         329   26.3595                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##Here we can see that introgression and genetic play a role in shaping the microbiome. THe intro
##The introgression status seemingly influcnes the OTUs
### Looking at composition across timepoint 3.2.21=====

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.05 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
#Time 1
NILT1<-subset(NIL.16S, Time=="T1")
NILT2<-subset(NIL.16S, Time=="T2")

NILT1.spec <- NILT1[,8:14031]
NILT1.treatment <- NILT1[,1:7]
filter.x<-filter(NILT1.spec)
rrna.rel <- decostand(filter.x, method="hellinger")
##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= rrna.rel~ Pedigree+Row*Range,NILT1.treatment, strata=NILT1.treatment$Range)

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   41    3.2530 0.079341  1.0097 0.24020  0.582  
# Row        21    1.7928 0.085371  1.0865 0.13238  0.014 *
#   Range       7    0.9534 0.136206  1.7334 0.07040  0.492  
# Residuals  96    7.5434 0.078577         0.55701         
# Total     165   13.5426                  1.00000         
# ---
##Time 2
NILT2.spec <- NILT2[,8:14031]
NILT2.treatment <- NILT2[,1:7]
filter.x<-filter(NILT2.spec)
rrna.rel <- decostand(filter.x, method="hellinger")
adonis(formula= rrna.rel~ Pedigree+Row+Range,NILT2.treatment,strata=NILT2.treatment$Row)
adonis(formula= NILT2.spec~ Pedigree+Row+Range,NILT2.treatment,strata=NILT2.treatment$Row)

##R2 is greater but not signficant when seperated in time. 

##Here Im going to looka at the teosine Maize compar
#NILF1TeoB73
##Flitered and Transfromed Data
NILF1TeoB73<-filter(NILF1TeoB73)
NIL.16S.spec <- NILF1TeoB73[,8:6836]
NIL.16S.treatment <- NILF1TeoB73[,1:7]
NIL.16S.spec <- decostand(NIL.16S.spec, method="hellinger")

NILF1TeoB73<-subset(NIL.16S, Pedigree=="B73"|Pedigree=="B73xPI384071"|Pedigree=="PI384071")
adonis(formula= NIL.16S.spec~ Pedigree,NIL.16S.treatment)

NILF1TeoB73T2<-subset(NILF1TeoB73, Time=="T1")
NIL.16S.spec <- NILF1TeoB73T2[,8:6836]
NIL.16S.treatment <- NILF1TeoB73T2[,1:7]

adonis(formula= NIL.16S.spec~ Pedigree,NIL.16S.treatment)

##Here I, going to see how the status look when I find the genotype average ===
##Not sure if this would actually be usefull in this setting 3.2.21
NIL.Means<-aggregate(NIL.16S[, (8:14031)], list(RangeRow=NIL.16S$RangeRow,
                                                                   Pedigree=NIL.16S$Pedigree, 
                                                                   Range=NIL.16S$Range,
                                                                  Row=NIL.16S$Row ), mean)

NIL.Means<-subset(NIL.Means, Pedigree!="Bulk")

#This is to pull out the data 
NIL.16S.treatment.mean <- NIL.Means[,1:4]
NIL.16S.spec.mean <- NIL.Means[,5:14028]

filter.x<-filter(NIL.16S.spec.mean)
NIL.16S.spec.mean <- decostand(filter.x, method="hellinger")

adonis(formula= NIL.16S.spec.mean~ Pedigree+Row+Range,strata=NIL.16S.treatment.mean$Range,NIL.16S.treatment.mean)
##Averaging by the block didnt really imporve the model differences. Which is run I think

##The data is close to being sig
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   44    2.6278 0.059722  1.0810 0.23609  0.310    
# Row        24    1.7715 0.073811  1.3360 0.15916  0.056 .  
# Range       7    0.9852 0.140745  2.5475 0.08852  0.001 ***
#   Residuals 104    5.7459 0.055249         0.51624           
# Total     179   11.1303                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

NIL.16S.spec <- NIL.16S[,8:14031]
NIL.16S.treatment <- NIL.16S[,1:7]
filter.x<-filter(NIL.16S.spec)
rrna.rel <- decostand(filter.x, method="hellinger")

#It seems like I need to incorrporate the RangeRow facotrs into the model and control for it. Positioning seems to play an extremely important role
adonis(formula = rrna.rel ~ Pedigree + Time + Row + Range + RangeRow,      data = NIL.16S.treatment, strata = NIL.16S.treatment$RangeRow) 
adonis(formula = rrna.rel ~ Pedigree * Time + Row + Range + RangeRow,      data = NIL.16S.treatment, strata = NIL.16S.treatment$RangeRow) 

## intorgression model=====
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   41    3.6378 0.08873  1.1545 0.12826  0.001 ***
#   Time        1    0.3922 0.39217  5.1029 0.01383  0.001 ***
#   Row        21    2.1062 0.10029  1.3050 0.07426  0.001 ***
#   Range       7    1.3530 0.19329  2.5151 0.04771  0.001 ***
#   RangeRow   98    8.4999 0.08673  1.1286 0.29969  0.001 ***
#   Residuals 161   12.3733 0.07685         0.43626           
# Total     329   28.3624                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Here we can see that genotype*Time interaction was not sigificant. 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree       41    3.6378 0.08873  1.1531 0.12826  0.001 ***
#   Time            1    0.3922 0.39217  5.0964 0.01383  0.001 ***
#   Row            21    2.1062 0.10029  1.3034 0.07426  0.001 ***
#   Range           7    1.3530 0.19329  2.5119 0.04771  0.001 ***
#   RangeRow       98    8.4999 0.08673  1.1271 0.29969  0.001 ***
#   Pedigree:Time  41    3.1393 0.07657  0.9950 0.11068  0.586    
# Residuals     120    9.2340 0.07695         0.32557           
# Total         329   28.3624                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###Now I would like to plot the centroids of the Data to get a better idea of the microbial community distrbution
#we are running centriods. We are using tapply, creates a new data fame based on the pervious 
#Transfroming the ordinational data into other states

#First we will perfrom a correspondence analysis on the matrix
cca.16S <- vegan::cca(NIL.16S.spec)
cca.16S
ordiplot(cca.16S, type="n")

COL.treat=(c(rainbow(50)))
points(cca.16S$CA$u, 
       col=COL.treat[NIL.16S.treatment$Pedigree],
       pch=16)
str(cca.16S)
cca.16S$CA$u[,1:2]

env.16S <- envfit(cca.16S,NIL.16S.treatment)

# cca.16SnoOL <- vegan::cca(MM.16S.spec[-c(161,302),])
# cca.16SnoOL
# ordiplot(cca.16SnoOL, type="n")
# points(cca.16SnoOL$CA$u, 
#        col=COL.treat[MM.16S.treatment[-c(161,302),]$Pedigree],
#        pch=16)

cca.16Spoints <- cbind(NIL.16S.treatment,cca.16S$CA$u[,1:2])

#All of the steps before this were to just place the points in ordinatioal space

cca.16S.means <- cbind(
  tapply(cca.16Spoints$CA1,cca.16Spoints$Pedigree, 
         FUN=mean),
  tapply(cca.16Spoints$CA2,cca.16Spoints$Pedigree, 
         FUN=mean)
)

#Function finds the distrbution of the points
sem <- function(x) sd(x)/sqrt(length(x))

cca.16S.sem <- cbind(
  tapply(cca.16Spoints$CA1,cca.16Spoints$Pedigree, 
         FUN=sem),
  tapply(cca.16Spoints$CA2,cca.16Spoints$Pedigree, 
         FUN=sem)
)


ordiplot(cca.16S, type="n",
         xlim=c(-2,1),ylim=c(-.5,1.5))
points(cca.16S.means, 
       col=COL.treat,
       pch=16)
arrows(cca.16S.means[,1], cca.16S.means[,2], 
       cca.16S.means[,1]+cca.16S.sem[,1], cca.16S.means[,2],
       lwd=1, length=0.025, angle=90)
arrows(cca.16S.means[,1], cca.16S.means[,2], 
       cca.16S.means[,1]-cca.16S.sem[,1], cca.16S.means[,2],
       lwd=1, length=0.025, angle=90)
arrows(cca.16S.means[,1], cca.16S.means[,2], 
       cca.16S.means[,1], cca.16S.means[,2]+cca.16S.sem[,2],
       lwd=1, length=0.025, angle=90)
arrows(cca.16S.means[,1], cca.16S.means[,2], 
       cca.16S.means[,1], cca.16S.means[,2]-cca.16S.sem[,2],
       lwd=1, length=0.025, angle=90)

with(row.names(cca.16S.means),text,cca.16S.means)

cbind(row.names(cca.16S.means),COL.treat)



NIL.16S.Filtered=rbind(
  subset(NIL.16S, Pedigree=="B73"),
  subset(NIL.16S, Pedigree=="PI384071"),
  subset(NIL.16S, Pedigree=="Z031E0047"),
  subset(NIL.16S,Pedigree=="Z031E0021"),
  subset(NIL.16S,Pedigree=="Z031E0560"),
  subset(NIL.16S, Pedigree=="Z031E0545")
)

NIL.16S.spec <- NIL.16S.Filtered[,8:140]
NIL.16S.treatment <- NIL.16S.Filtered[,1:7]

MDS <- metaMDS(NIL.16S.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
points(MDS$points,  
       col=COL.treat[NIL.16S.treatment$Pedigree],
       pch=18, cex=2)

MDS.points <- cbind(NIL.16S.treatment,MDS$points[,1:2])
boxplot(MDS.points$MDS1 ~(MDS.points$Pedigree))
boxplot(MDS.points$MDS2 ~(MDS.points$Pedigree))

#Screet plot; Ploting the varaince. 
#I've optained the ordniationl means of the vaules. on one axis the hybrid, inbred and B73 are very similar. 
#It seems like Ill need to pull out the OTUs that are driving the primiary differences in these classes. 
#It may also be useful to 
#Something that may be useful will be to fgo back into the breeding book and prima for ecological genetics and redefine the questiosn 
#Right now it is clear to me that I will need to use a more complicated 
#########

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.005 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}

###nifH NonSign in ANOVA ====
############
#nihH
library(vegan)

NIL.nifh <- read.csv("otu_table_nifH_500.csv")
NIL.nifh<-read.csv("otu_nifH.csv")
#This Function below should be used to get rid of sample s
NIL.nifh<-NIL.nifh[-c(226,229),]
#Remove all of the low abundance samples
#Outliers are
#T2.090
#T2.078
#T1.131
#T2.144

# All = read.csv(file ="otu_nifH.csv")
# row.names(All) <- All$X
# NIL.nifh=rbind(
#   subset(All, Pedigree=="B73xPI384071"),
#   subset(All, Pedigree=="B73"),
#   subset(All, Pedigree=="PI384071"))

NIL.nifh<-subset(NIL.nifh, Pedigree!="B73xPI384071")
NIL.nifh<-subset(NIL.nifh, Pedigree!="PI384071")
NIL.nifh<-subset(NIL.nifh,Pedigree!="Bulk")
NIL.nifh<-subset(NIL.nifh,Pedigree!="Bulk ")


NIL.nifh.spec <- NIL.nifh[,7:1685]
NIL.nifh.treatment <- NIL.nifh[,1:6]
library(tidyverse)
NIL.nifh.treatment<-NIL.nifh.treatment %>% unite(RangeRow, c(Range, Row), remove = FALSE)


##Testing different ANOVA Models
adonis(NIL.nifh.spec ~ Time, NIL.nifh.treatment)
adonis(NIL.nifh.spec ~ Pedigree, NIL.nifh.treatment)
adonis(NIL.nifh.spec ~ Range, NIL.nifh.treatment)
adonis(NIL.nifh.spec ~ Row, NIL.nifh.treatment)

adonis(formula= NIL.nifh.spec~ Time+Pedigree,NIL.nifh.treatment)
adonis(formula= NIL.nifh.spec~ Pedigree+Time+RangeRow,strata=NIL.nifh.treatment$RangeRow,NIL.nifh.treatment)
adonis(formula= NIL.nifh.spec~ Pedigree+Time+Range,strata=NIL.nifh.treatment$Range,NIL.nifh.treatment)
adonis(formula= NIL.nifh.spec~ Pedigree+Time+Row+Range,strata=NIL.nifh.treatment$Row,NIL.nifh.treatment)
adonis(formula= NIL.nifh.spec~ Pedigree+Time+Row+Range,strata=NIL.nifh.treatment$Range,NIL.nifh.treatment)

adonis(formula= NIL.nifh.spec~ Pedigree+Time+Row+Range+RangeRow,strata=NIL.nifh.treatment$RangeRow,NIL.nifh.treatment)
#anosim test====
anosim(NIL.nifh.spec, NIL.nifh.treatment$Pedigree)

#Testing it under filering
X.filtered = filter(NIL.nifh.spec)
rrna.rel <- decostand(X.filtered, method="hellinger")
adonis(rrna.rel ~ Pedigree, NIL.nifh.treatment)
adonis(formula= X.filtered~ Pedigree+Time+Row+Range,strata=NIL.nifh.treatment$Range,NIL.nifh.treatment)

##Here is the filtered data 
adonis(formula= X.filtered~ Pedigree+Time+Row+Range+RangeRow,strata=NIL.nifh.treatment$RangeRow,NIL.nifh.treatment)
##Here is the transfromed fata
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range+RangeRow,
       strata=NIL.nifh.treatment$RangeRow,NIL.nifh.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   45    19.551 0.43447  1.0246 0.13011  0.244  
# Time        1     0.490 0.49024  1.1561 0.00326  0.164  
# Row        23    10.971 0.47698  1.1248 0.07301  0.079 .
# Range       7     3.627 0.51812  1.2219 0.02414  0.277  
# RangeRow  104    45.234 0.43495  1.0257 0.30103  0.440  
# Residuals 166    70.391 0.42405         0.46845         
# Total     346   150.265                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Model Results of introgression
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range+RangeRow,
       strata=NIL.nifh.treatment$RangeRow,NIL.nifh.treatment)

#NifH introgression model =====
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   41    17.954 0.43790  1.0360 0.12677  0.101  
# Time        1     0.545 0.54519  1.2898 0.00385  0.054 .
# Row        21     9.795 0.46642  1.1034 0.06916  0.077 .
# Range       7     3.707 0.52952  1.2527 0.02617  0.214  
# RangeRow   98    42.833 0.43707  1.0340 0.30245  0.234  
# Residuals 158    66.787 0.42270         0.47159         
# Total     326   141.620                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Here is the unfiltered transformed data
NIL.nifh.spec.Trans <- decostand(NIL.nifh.spec, method="hellinger")
adonis(formula= NIL.nifh.spec.Trans~ Pedigree+Time+Row+Range+RangeRow,
       strata=NIL.nifh.treatment$RangeRow,
       NIL.nifh.treatment)


#Time 1
NILT1<-subset(NIL.nifh, Time=="T1")
NILT2<-subset(NIL.nifh, Time=="T2")


T1.nifh.spec <- NILT1[,7:1685]
T1.nifh.treatment <- NILT1[,1:6]

##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T1.nifh.spec~ Pedigree+Row+Range,T1.nifh.treatment)

#Time 2
T2.nifh.spec <- NILT2[,7:1685]
T2.nifh.treatment <- NILT2[,1:6]
T2.nifh.treatment<-T2.nifh.treatment %>% unite(RangeRow, c(Range, Row), remove = FALSE)

##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T2.nifh.spec~ Pedigree+Row+Range,T2.nifh.treatment)

##Wow the second time point acutally had an effect
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#   Pedigree   45    21.772 0.48383  1.0706 0.26789  0.004 **
#   Row        23    11.053 0.48055  1.0633 0.13599  0.031 * 
#   Range       7     3.709 0.52981  1.1724 0.04563  0.005 **
#   Residuals  99    44.741 0.45192         0.55049          
# Total     174    81.274                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# T2.nifh.spec.T <- decostand(T2.nifh.spec, method="hellinger")
# adonis(formula= T2.nifh.spec.T~ Pedigree+Row+Range+RangeRow,strata=T2.nifh.treatment$RangeRow,T2.nifh.treatment)
# 


#The filtering seems to result in some samples have a matrix full og zeros 

COL.treat <- c(rainbow(50))
MDS <- metaMDS(NIL.nifh.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
points(MDS$points,  
       col=COL.treat[NIL.nifh.treatment$Pedigree],
       pch=18, cex=2)

MDS.nifhpoints <- cbind(NIL.nifh.treatment,MDS$points[,1:2])

MDS.nifh.means <- cbind(
  tapply(MDS.nifhpoints$MDS1,MDS.nifhpoints$Pedigree, 
         FUN=mean),
  tapply(MDS.nifhpoints$MDS2,MDS.nifhpoints$Pedigree, 
         FUN=mean)
)

#There doesnt appear to be a signifcant effect on the 
#The Data has some serious outliers that need to be deal with to have meaningful 

All = read.csv(file ="otu_nifH.csv")
row.names(All) <- All$X
X=rbind(
  subset(All, Pedigree=="B73xPI384071"),
  subset(All, Pedigree=="B73"),
  subset(All, Pedigree=="PI384071"))

NIL.Set.spec <- X[,7:1799]
NIL.Set.treatment <- X[,1:6]

adonis(formula= NIL.Set.spec~ Pedigree+Time+Row+Range,strata=NIL.Set.treatment$Range,NIL.Set.treatment)
COL.treat<-c("Red","blue","Green","Green")

MDS <- metaMDS(NIL.Set.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
points(MDS$points,  
       col=COL.treat[NIL.Set.treatment$Pedigree],
       pch=18, cex=2)
##########

###ITS Signifcant in perMANOVA ====
#########
#

library(vegan)

NIL.ITS <- read.csv("otu_table_ITS_2302.csv")
###Below is to remove specifc samples from dataset
NIL.ITS<-NIL.ITS[-c(39),]

# All = read.csv(file ="otu_table_ITS_2302.csv")
# row.names(All) <- All$X
# NIL.ITS=rbind(
#   subset(All, Pedigree=="B73xPI384071"),
#   subset(All, Pedigree=="B73"),
#   subset(All, Pedigree=="PI384071"))
NIL.ITS<-subset(NIL.ITS, Pedigree!="B73xPI384071")
NIL.ITS<-subset(NIL.ITS, Pedigree!="PI384071")
NIL.ITS<-subset(NIL.ITS,Pedigree!="Bulk")



NIL.ITS.spec <- NIL.ITS[,8:2253]
NIL.ITS.treatment <- NIL.ITS[,1:7]
adonis(formula= NIL.ITS.spec~ Pedigree+Time+Row+Range+RangeRow,strata=NIL.ITS.treatment$Range,NIL.ITS.treatment)


filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.05 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.ITS.spec)


rrna.rel <- decostand(X.filtered, method="hellinger")
adonis(rrna.rel ~ Pedigree, NIL.ITS.treatment)
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range,strata=NIL.ITS.treatment$Time,NIL.ITS.treatment)
# This is when we treat time as a random effect 
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   44    15.443 0.35098  1.0424 0.12637  0.041 *  
#   Time        1     2.142 2.14161  6.3604 0.01752  0.002 ** 
#   Row        23     8.176 0.35549  1.0558 0.06691  0.059 .  
# Range       7     3.849 0.54991  1.6332 0.03150  0.001 ***
#   Residuals 275    92.595 0.33671         0.75770           
# Total     350   122.205                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(formula= rrna.rel~ Pedigree+Time+Row+Range+RangeRow,strata=NIL.ITS.treatment$RangeRow,NIL.ITS.treatment)
##Fungal ITS Introgression Model=====
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   41    13.138 0.32043  1.0334 0.12399  0.001 ***
#   Time        1     2.285 2.28518  7.3697 0.02157  0.001 ***
#   Row        21     6.805 0.32405  1.0451 0.06422  0.001 ***
#   Range       7     3.730 0.53279  1.7182 0.03520  0.001 ***
#   RangeRow   98    29.766 0.30374  0.9795 0.28093  0.001 ***
#   Residuals 162    50.233 0.31008         0.47409           
# Total     330   105.956                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Time 1
NILT1<-subset(NIL.ITS, Time=="T1")
NILT2<-subset(NIL.ITS, Time=="T2")


T1.ITS.spec <- NIL.ITS[,7:2252]
T1.ITS.treatment <- NIL.ITS[,1:6]

##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T1.ITS.spec~ Pedigree+Row+Range,T1.ITS.treatment)
#non-sigificant. 

#Time 2
T2.ITS.spec <- NIL.ITS[,7:2252]
T2.ITS.treatment <- NIL.ITS[,1:6]

##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T2.ITS.spec~ Pedigree+Row+Range,T2.ITS.treatment)
#Non-sig



COL.treat <- c(rainbow(5))
COL.treat<-c("Red","blue")
MDS <- metaMDS(NIL.ITS.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
points(MDS$points,  
       col=COL.treat[NIL.ITS.treatment$Pedigree],
       pch=18, cex=2)
MDS.ITSpoints <- cbind(NIL.ITS.treatment,MDS$points[,1:2])

#Try to work with Data that is not rarefyed 
###nosZ Nonsingicant in ANOVA ====
##########

library(vegan)

NIL.nosZ <- read.csv("otu_table_nosZ_294.csv")
###Below is to remove specifc samples from dataset and cleaning up the data
NIL.nosZ<-NIL.nosZ[-c(201,225),]

##Here we subset the samples for the introgressopm
NIL.nosZ<-subset(NIL.nosZ, Pedigree!="B73xPI384071")
NIL.nosZ<-subset(NIL.nosZ, Pedigree!="PI384071")
NIL.nosZ<-subset(NIL.nosZ,Pedigree!="Bulk")
NIL.nosZ<-subset(NIL.nosZ,Pedigree!="Bulk ")

NIL.nosZ.spec <- NIL.nosZ[,7:1661]
NIL.nosZ.treatment <- NIL.nosZ[,1:6]

adonis(formula= NIL.nosZ.spec~ Pedigree+Time+Row+Range,strata=NIL.nosZ.treatment$Time,NIL.nosZ.treatment)

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.005 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.nosZ.spec)

#Here we are making the range row data file- This acts as the block 
NIL.nosZ.treatment<-NIL.nosZ.treatment %>% unite(RangeRow, c(Range, Row), remove = FALSE)

rrna.rel <- decostand(X.filtered, method="hellinger")
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range,strata=NIL.nosZ.treatment$Time,NIL.nosZ.treatment)

##This is the introgression model 
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range+RangeRow,
       strata=NIL.nosZ.treatment$RangeRow,
       NIL.nosZ.treatment)

##It seems like this is the model that holds the most value. I need to control for spatial effects 
##These regional starting microbes play a major role in shaping th microbes
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# Pedigree   44    19.394 0.44077  1.0060 0.12447  0.012 * 
#   Time        1     0.651 0.65072  1.4852 0.00418  0.007 **
#   Row        24    10.807 0.45027  1.0277 0.06936  0.048 * 
#   Range       7     3.587 0.51246  1.1697 0.02302  0.142   
# RangeRow  104    46.887 0.45084  1.0290 0.30093  0.031 * 
#   Residuals 170    74.483 0.43813         0.47804          
# Total     350   155.809                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#NosZ Model Statisitics PERMANOVA Model====
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   41    18.163 0.44300  1.0116 0.12399  0.046 *
#   Time        1     0.609 0.60897  1.3906 0.00416  0.021 *
#   Row        21     9.394 0.44736  1.0216 0.06413  0.087 .
# Range       7     3.649 0.52129  1.1904 0.02491  0.229  
# RangeRow   98    44.170 0.45071  1.0292 0.30152  0.106  
# Residuals 161    70.504 0.43791         0.48129         
# Total     329   146.489                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Set Colors
COL.treat <- c(rainbow(50))
COL.treat<-c("Red","blue")

MDS <- metaMDS(NIL.nosZ.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n")
points(MDS$points,  
       col=COL.treat[NIL.nosZ.treatment$Pedigree],
       pch=18, cex=2)
MDS.nosZpoints <- cbind(NIL.nosZ.treatment,MDS$points[,1:2])

#Seems like there may be a potential statically important values
#Potentially it may be useful to look at everysingle otu seperatly 
###########
#Time 1
NILT1<-subset(NIL.nosZ, Time=="T1")
NILT2<-subset(NIL.nosZ, Time=="T2")


T1.nosZ.spec <- NILT1[,7:1661]
T1.nosZ.treatment <- NILT1[,1:6]

##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T1.nosZ.spec~ Pedigree+Row+Range,T1.nosZ.treatment)
#non-sigificant. 

#Time 2
T2.nosZ.spec <- NILT2[,7:1661]
T2.nosZ.treatment <- NILT2[,1:6]

##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T2.nosZ.spec~ Pedigree+Row+Range,T2.nosZ.treatment)
#Non-sig




###nirK:: Signifcantly infelucned by the NILs=====
##########
library(vegan)

NIL.nirK <- read.csv("otu_table_nirK_5004.csv")
###Below is to remove specifc samples from dataset
NIL.nirK<-NIL.nirK[-c(201,225),]
NIL.nirK<-subset(NIL.nirK, Pedigree!="B73xPI384071")
NIL.nirK<-subset(NIL.nirK, Pedigree!="PI384071")
NIL.nirK<-subset(NIL.nirK,Pedigree!="Bulk")
NIL.nirK<-subset(NIL.nirK,Pedigree!="Bulk ")


NIL.nirK.spec <- NIL.nirK[,7:4743]
NIL.nirK.treatment <- NIL.nirK[,1:6]
adonis(formula= NIL.nirK.spec~ Pedigree+Time+Row+Range,strata=NIL.nirK.treatment$Time,NIL.nirK.treatment)

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.005 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.nirK.spec)

NIL.nirK.treatment<-NIL.nirK.treatment %>% unite(RangeRow, c(Range, Row), remove = FALSE)


rrna.rel <- decostand(X.filtered, method="hellinger")
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range,strata=NIL.nirK.treatment$Time,NIL.nirK.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   45    17.341 0.38537  1.0419 0.13089  0.063 .  
# Time        1     0.510 0.51011  1.3792 0.00385  0.068 .  
# Row        23     9.268 0.40295  1.0895 0.06995  0.010 ** 
#   Range       7     3.293 0.47040  1.2719 0.02485  0.001 ***
#   Residuals 276   102.080 0.36985         0.77046           
# Total     352   132.492                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#This is the model that Im interested in using-- it FUCKING works 
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range+RangeRow,
       strata=NIL.nirK.treatment$RangeRow,
       NIL.nirK.treatment)
##The introgression model would have a DF of 42|| This model has bulk soil and 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   45    17.286 0.38414  1.0614 0.13135  0.014 *
#   Time        1     0.513 0.51267  1.4165 0.00390  0.016 *
#   Row        23     9.252 0.40224  1.1114 0.07030  0.012 *
#   Range       7     3.247 0.46379  1.2814 0.02467  0.091 .
# RangeRow  104    39.783 0.38253  1.0569 0.30229  0.025 *
#   Residuals 170    61.528 0.36193         0.46751         
# Total     350   131.609                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##NirK MANOVA Results-----

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#   Pedigree   41    15.609 0.38070  1.0523 0.12643  0.022 *
#   Time        1     0.506 0.50563  1.3977 0.00410  0.020 *
#   Row        21     8.414 0.40067  1.1075 0.06815  0.018 *
#   Range       7     3.199 0.45702  1.2633 0.02591  0.140  
#   RangeRow   98    37.489 0.38254  1.0574 0.30365  0.031 *
#   Residuals 161    58.245 0.36177         0.47177         
# Total     329   123.461                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Time 1
NILT1<-subset(NIL.nirK, Time=="T1")
NILT2<-subset(NIL.nirK, Time=="T2")


T1.nirK.spec <- NILT1[,7:4743]
T1.nirK.treatment <- NILT1[,1:6]

X.filtered = filter(T1.nirK.spec)
T1.nirK.spec <- decostand(X.filtered, method="hellinger")
##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T1.nirK.spec~ Pedigree+Row+Range,T1.nirK.treatment)
#non-sigificant. 

#Time 2
T2.nirK.spec <- NIL.nirK[,7:4743]
T2.nirK.treatment <- NIL.nirK[,1:6]

X.filtered = filter(T2.nirK.spec)
T2.nirK.spec <- decostand(X.filtered, method="hellinger")
adonis(formula= T2.nirK.spec~ Pedigree+Row+Range,T2.nirK.treatment)


# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   45    17.286 0.38414  1.0374 0.13135  0.083 .  
# Row        23     9.247 0.40203  1.0858 0.07026  0.014 *  
#   Range       7     3.249 0.46417  1.2536 0.02469  0.001 ***
#   Residuals 275   101.827 0.37028         0.77371           
# Total     350   131.609                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Set Colors
COL.treat <- c(rainbow(50))
COL.treat<-c("Red","blue")

MDS <- metaMDS(NIL.nirK.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n")
points(MDS$points,  
       col=COL.treat[NIL.nirK.treatment$Pedigree],
       pch=18, cex=2)
MDS.nirKpoints <- cbind(NIL.nirK.treatment,MDS$points[,1:2])


##########

###nirS =====
##########
library(vegan)

NIL.nirS <- read.csv("otu_table_nirS_707.csv")
###Below is to remove specifc samples from dataset
NIL.nirS<-NIL.nirS[-c(201,225),]
NIL.nirS<-subset(NIL.nirS, Pedigree!="B73xPI384071")
NIL.nirS<-subset(NIL.nirS, Pedigree!="PI384071")
NIL.nirS<-subset(NIL.nirS,Pedigree!="Bulk")
NIL.nirS<-subset(NIL.nirS,Pedigree!="Bulk ")

NIL.nirS.spec <- NIL.nirS[,7:1193]
NIL.nirS.treatment <- NIL.nirS[,1:6]
adonis(formula= NIL.nirS.spec~ Pedigree+Time+Row+Range,strata=NIL.nirS.treatment$Time,NIL.nirS.treatment)

filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.005 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.nirS.spec)


rrna.rel <- decostand(X.filtered, method="hellinger")
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range,strata=NIL.nirS.treatment$Range,NIL.nirS.treatment)

##Working Model 
NIL.nirS.treatment<-NIL.nirS.treatment %>% unite(RangeRow, c(Range, Row), remove = FALSE)
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range+RangeRow,
       strata=NIL.nirS.treatment$RangeRow,
       NIL.nirS.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   45    12.823 0.28496  1.0950 0.13023  0.220  
# Time        1     0.299 0.29936  1.1503 0.00304  0.234  
# Row        23     7.840 0.34087  1.3098 0.07962  0.078 .
# Range       7     2.865 0.40927  1.5727 0.02910  0.298  
# RangeRow  104    30.659 0.29480  1.1328 0.31136  0.291  
# Residuals 169    43.981 0.26024         0.44666         
# Total     349    98.467                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#NirS PERMANOVA Model=====
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Pedigree   41    11.621 0.28345  1.0947 0.12665  0.273
# Time        1     0.310 0.31033  1.1985 0.00338  0.200
# Row        21     7.162 0.34103  1.3171 0.07805  0.209
# Range       7     2.829 0.40410  1.5607 0.03083  0.265
# RangeRow   98    28.926 0.29517  1.1400 0.31524  0.441
# Residuals 158    40.910 0.25892         0.44584       
# Total     326    91.758                 1.00000       


#Time 1
NILT1<-subset(NIL.nirS, Time=="T1")
NILT2<-subset(NIL.nirS, Time=="T2")


T1.nirS.spec <- NILT1[,7:1193]
T1.nirS.treatment <- NILT1[,1:6]

X.filtered = filter(T1.nirS.spec)
T1.nirS.spec <- decostand(X.filtered, method="hellinger")
##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T1.nirS.spec~ Pedigree+Row+Range,T1.nirS.treatment)
#non-sigificant. 

#Time 2
T2.nirS.spec <- NILT2[,7:1193]
T2.nirS.treatment <- NILT2[,1:6]

X.filtered = filter(T2.nirS.spec)
T2.nirS.spec <- decostand(X.filtered, method="hellinger")
adonis(formula= T2.nirS.spec~ Pedigree+Row+Range,T2.nirS.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   45    12.698 0.28217  1.0057 0.25591  0.449  
# Row        23     6.829 0.29691  1.0583 0.13764  0.138  
# Range       7     2.314 0.33062  1.1784 0.04664  0.030 *
#   Residuals  99    27.776 0.28056         0.55981         
# Total     174    49.617                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




#Set Colors
COL.treat <- c(rainbow(50))
COL.treat<-c("Red","blue")

MDS <- metaMDS(NIL.nirS.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n")
points(MDS$points,  
       col=COL.treat[NIL.nirS.treatment$Pedigree],
       pch=18, cex=2)
MDS.nirSpoints <- cbind(NIL.nirS.treatment,MDS$points[,1:2])

boxplot(MDS.nirSpoints$MDS1 ~(MDS.nirSpoints$Pedigree))
boxplot(MDS.nirSpoints$MDS2 ~(MDS.nirSpoints$Pedigree))
############

##amoA Bacteria =====
##########
library(vegan)

NIL.amoA <- read.csv("otu_table_amoA_144.csv")
###Below is to remove specifc samples from dataset
NIL.amoA<-NIL.amoA[-c(12),]
NIL.amoA<-subset(NIL.amoA, Pedigree!="B73xPI384071")
NIL.amoA<-subset(NIL.amoA, Pedigree!="PI384071")
NIL.amoA<-subset(NIL.amoA,Pedigree!="Bulk")
NIL.amoA<-subset(NIL.amoA,Pedigree!="Bulk ")

NIL.amoA.spec <- NIL.amoA[,7:61]
NIL.amoA.treatment <- NIL.amoA[,1:6]
adonis(formula= NIL.amoA.spec~ Pedigree+Time+Row+Range,strata=NIL.amoA.treatment$Range,NIL.amoA.treatment)
anosim(NIL.amoA.spec, NIL.amoA.treatment$Time)


# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   44    13.212 0.30028  1.0132 0.13076  0.618  
# Time        1     0.583 0.58283  1.9667 0.00577  0.059 .
# Row        23     6.938 0.30166  1.0179 0.06867  0.526  
# Range       7     4.144 0.59198  1.9975 0.04101  0.248  
# Residuals 257    76.163 0.29636         0.75379         
# Total     332   101.040                 1.00000         
# ---
#   
#   
filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.005 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.amoA.spec)


rrna.rel <- decostand(X.filtered, method="hellinger")
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range,strata=NIL.amoA.treatment$Range,NIL.amoA.treatment)

##Here we are combing the range row and using it as a block factor
NIL.amoA.treatment<-NIL.amoA.treatment %>% unite(RangeRow, c(Range, Row), remove = FALSE)
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range+RangeRow,
       strata=NIL.amoA.treatment$RangeRow,
       NIL.amoA.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Pedigree   44    11.587 0.26333 0.98430 0.12766  0.016 *
#   Time        1     0.617 0.61665 2.30494 0.00679  0.022 *
#   Row        23     6.560 0.28521 1.06606 0.07227  0.037 *
#   Range       7     3.816 0.54518 2.03779 0.04205  0.032 *
#   RangeRow  104    27.250 0.26202 0.97938 0.30023  0.015 *
#   Residuals 153    40.932 0.26753         0.45099         
# Total     332    90.761                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Bact amoA PERMANOVA=====
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# Pedigree   41    10.688 0.26068 0.97963 0.12605  0.008 **
#   Time        1     0.595 0.59455 2.23433 0.00701  0.015 * 
#   Row        21     5.455 0.25978 0.97626 0.06434  0.083 . 
# Range       7     3.622 0.51742 1.94445 0.04272  0.026 * 
#   RangeRow   98    25.843 0.26371 0.99102 0.30480  0.005 **
#   Residuals 145    38.584 0.26610         0.45507          
# Total     313    84.787                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Time 1
NILT1<-subset(NIL.amoA, Time=="T1")
NILT2<-subset(NIL.amoA, Time=="T2")


T1.amoA.spec <- NILT1[,7:61]
T1.amoA.treatment <- NILT1[,1:6]

X.filtered = filter(T1.amoA.spec)
T1.amoA.spec <- decostand(X.filtered, method="hellinger")
##The adonis model is giving me an issue with the number of permuations needed. Basically I dont have enough data for the complex model
adonis(formula= T1.amoA.spec~ Pedigree+Row+Range,T1.amoA.treatment)
#non-sigificant. 

#Time 2
T2.amoA.spec <- NILT2[,7:61]
T2.amoA.treatment <- NILT2[,1:6]

X.filtered = filter(T2.amoA.spec)
T2.amoA.spec <- decostand(X.filtered, method="hellinger")
adonis(formula= T2.amoA.spec~ Pedigree+Row+Range,T2.amoA.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Pedigree   44    13.165 0.29921  1.1718 0.28442  0.050 *  
#   Row        23     6.865 0.29848  1.1689 0.14831  0.102    
#   Range       7     3.532 0.50452  1.9759 0.07630  0.001 ***
#   Residuals  89    22.726 0.25534         0.49096           
# Total     163    46.288                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Set Colors
COL.treat <- c(rainbow(50))
COL.treat<-c("Red","blue")

MDS <- metaMDS(T2.amoA.spec, trymax=50, perm=1e6)
ordiplot(MDS, type="n")
points(MDS$points,  
       col=COL.treat[T2.amoA.treatment$Pedigree],
       pch=18, cex=2)
MDS.amoApoints <- cbind(T2.amoA.treatment,MDS$points[,1:2])

##Plotting NMDS Axis
ggplot(MDS.amoApoints, aes(x = fct_reorder(Pedigree, MDS1, .fun = median, .desc=TRUE), y=MDS1, fill=Pedigree))+
  geom_boxplot()+labs(y="NMDS 1", title="Differences in T2 amoA community", x="Pedigree")+
  theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(MDS.amoApoints, aes(x = fct_reorder(Pedigree, MDS2, .fun = median, .desc=TRUE), y=MDS2, fill=Pedigree))+
  geom_boxplot()+labs(y="NMDS 2", title="Differences in T2 amoA community", x="Pedigree")+
  theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

boxplot(MDS.amoApoints$MDS1 ~(MDS.amoApoints$Pedigree))
boxplot(MDS.amoApoints$MDS2 ~(MDS.amoApoints$Pedigree))

#Im going to subset our NILS of interest. 
#This step allows me to bind together different NIls to have a better view of the 

NIL.amoA.Filtered=rbind(
  subset(NIL.amoA, Pedigree=="B73"),
  subset(NIL.amoA, Pedigree=="PI384071"),
  subset(NIL.amoA, Pedigree=="Z031E0047"),
  subset(NIL.amoA,Pedigree=="Z031E0021"),
  subset(NIL.amoA,Pedigree=="Z031E0560"),
  subset(NIL.amoA, Pedigree=="Z031E0545")
)

NIL.amoA.data <- NIL.amoA.Filtered[,7:61]
NIL.amoA.treat <- NIL.amoA.Filtered[,,1:6]

adonis(NIL.amoA.data~ Pedigree, NIL.amoA.treat)

COL.treat <- c(rainbow(50))
COL.treat<-c("Red","blue")

MDS <- metaMDS(NIL.amoA.data, trymax=50, perm=1e6)
ordiplot(MDS, type="n")
points(MDS$points,  
       col=COL.treat[NIL.amoA.treat$Pedigree],
       pch=18, cex=2)
MDS.amoApoints <- cbind(NIL.amoA.treat,MDS$points[,1:2])

boxplot(MDS.amoApoints$MDS1 ~(MDS.amoApoints$Pedigree))
boxplot(MDS.amoApoints$MDS2 ~(MDS.amoApoints$Pedigree))
###Seems like the time points is having a strong effect on our results. 
#Lets seperate it out babay

NIL.amoA.T1<-subset(NIL.amoA, Time=="T1")
NIL.amoA.T2<-subset(NIL.amoA, Time=="T2")

NIL.amoA.Filtered=rbind(
  subset(NIL.amoA.T2, Pedigree=="B73"),
  subset(NIL.amoA.T2, Pedigree=="PI384071"),
  subset(NIL.amoA.T2, Pedigree=="Z031E0047"),
  subset(NIL.amoA.T2,Pedigree=="Z031E0021"),
  subset(NIL.amoA.T2,Pedigree=="Z031E0560"),
  subset(NIL.amoA.T2, Pedigree=="Z031E0545")
)

NIL.amoA.data <- NIL.amoA.Filtered[,7:20]
NIL.amoA.treat <- NIL.amoA.Filtered[,1:6]

adonis(NIL.amoA.data~ Pedigree, NIL.amoA.treat)

MDS <- metaMDS(NIL.amoA.data, trymax=50, perm=1e6)
ordiplot(MDS, type="n")
points(MDS$points,  
       col=COL.treat[NIL.amoA.treat$Pedigree],
       pch=18, cex=2)
MDS.amoApoints <- cbind(NIL.amoA.treat,MDS$points[,1:2])

boxplot(MDS.amoApoints$MDS1 ~(MDS.amoApoints$Pedigree))
boxplot(MDS.amoApoints$MDS2 ~(MDS.amoApoints$Pedigree))

object<-simper(NIL.amoA.data, NIL.amoA.treat$Pedigree, permutations = 0, trace = FALSE, parallel = getOption("mc.cores"), )
## S3 method for class 'simper'
summary(object, ordered = TRUE, digits = max(3,getOption("digits") - 3), )

##Arch amoA ====
library(tidyverse)
NIL.arcamoA <- read.csv("arch_amoA_149.csv")

NIL.arcamoA<-subset(NIL.arcamoA, Pedigree!="B73xPI384071")
NIL.arcamoA<-subset(NIL.arcamoA, Pedigree!="PI384071")
NIL.arcamoA<-subset(NIL.arcamoA,Pedigree!="Bulk")
NIL.arcamoA<-subset(NIL.arcamoA,Pedigree!="Bulk ")

NIL.arcamoA.spec <- NIL.arcamoA[,7:80]
NIL.archamoA.treatment <- NIL.arcamoA[,1:6]



adonis(formula= NIL.arcamoA.spec~ Pedigree+Time+Row+Range,strata=NIL.archamoA.treatment$Range,NIL.archamoA.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   44    4.8491 0.11021  1.4269 0.15868  0.001 ***
#   Time        1    0.5871 0.58708  7.6015 0.01921  0.001 ***
#   Row        24    2.6474 0.11031  1.4283 0.08663  0.012 *  
#   Range       7    1.3139 0.18770  2.4304 0.04300  0.101    
# Residuals 274   21.1617 0.07723         0.69248           
# Total     350   30.5592                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##Filter the Data 
filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.005 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
X.filtered = filter(NIL.arcamoA.spec)


rrna.rel <- decostand(X.filtered, method="hellinger")
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range,strata=NIL.archamoA.treatment$Range,NIL.archamoA.treatment)

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   44    4.1598 0.094540  1.3388 0.15271  0.001 ***
#   Time        1    0.2933 0.293321  4.1539 0.01077  0.002 ** 
#   Row        24    2.3367 0.097361  1.3788 0.08578  0.010 ** 
#   Range       7    1.1012 0.157308  2.2277 0.04043  0.094 .  
# Residuals 274   19.3481 0.070613         0.71031           
# Total     350   27.2390                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Here we are combing the range row and using it as a block factor
NIL.archamoA.treatment<-NIL.archamoA.treatment %>% unite(RangeRow, c(Range, Row), remove = FALSE)
adonis(formula= rrna.rel~ Pedigree+Time+Row+Range+RangeRow,
       strata=NIL.archamoA.treatment$RangeRow,
       NIL.archamoA.treatment)

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#   Pedigree   44    4.1598 0.094540  1.5940 0.15271  0.001 ***
#   Time        1    0.2933 0.293321  4.9456 0.01077  0.001 ***
#   Row        24    2.3367 0.097361  1.6416 0.08578  0.001 ***
#   Range       7    1.1012 0.157308  2.6523 0.04043  0.001 ***
#   RangeRow  104    9.2654 0.089090  1.5021 0.34015  0.001 ***
#   Residuals 170   10.0827 0.059310         0.37016           
# Total     350   27.2390                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Arch amoA PERMANOVA Model====

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Pedigree   41    3.8991 0.095101  1.6246 0.15176  0.001 ***
#   Time        1    0.2996 0.299628  5.1185 0.01166  0.001 ***
#   Row        21    1.9779 0.094187  1.6090 0.07698  0.001 ***
#   Range       7    1.0956 0.156511  2.6737 0.04264  0.001 ***
#   RangeRow   98    8.9960 0.091796  1.5681 0.35014  0.001 ***
#   Residuals 161    9.4246 0.058538         0.36682           
# Total     329   25.6929                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Here Im going to try env fit on the community
#nmds code
set.seed(123)
nmds = metaMDS(rrna.rel, distance = "bray")
nmds
en = envfit(nmds, NIL.archamoA.treatment, permutations = 999, na.rm = TRUE)
en
plot(nmds)
plot(en)

dis <- vegdist(rrna.rel)        # Bray-Curtis distances between samples
groups <- NIL.archamoA.treatment$Pedigree

mod <- betadisper(dis, groups)  # Calculate multivariate dispersions
plot(mod)
plot(mod, hull = FALSE, ellipse = TRUE)
plot.betadisper(mod)

anova(mod)
##This test is determining whether group is a signifanct factor in dispersal from centeriod
# Response: Distances
# Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups     44 0.26703 0.0060689  1.0489 0.3946
# Residuals 306 1.77045 0.0057858     

## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)


####Looking at Rhizosphere samples from NOVA sq Run 9/14/2021=====

#Importing Biom File
biomOTU = import_biom('Phyloseq/Rhiz_otu_table_rarefied_132893.biom')
#Importaning Mapping
NILMap = import_qiime_sample_data("NIL2018.RhizoMapping.txt")
#Functional Bins are present
#NILMap = import_qiime_sample_data("Sample_map.txt")
#NILMap = import_qiime_sample_data("Sample_map_ActualHighLow.txt")

#Setting themes 
theme_set(theme_bw())

##Loading and Editing Data ====
#Combining Biom and Mapping file to Make the phyloseq file
NILOTU = merge_phyloseq(biomOTU, NILMap)
NILOTU <- subset_samples(NILOTU, X.SampleID != "NIL.144")

plot_richness(NILOTU, x = "Pedigree", measures = c("Observed", "Chao1", "Simpson", "Shannon"), sortby = "Shannon") + geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 90))+labs( x="Treatment", y="Alpha Diversity")+ggtitle("16S rRNA Diversity")


NitrifersOTU=subset_taxa(NILOTU, Rank2 == "D_1__Nitrospirae"|Rank4 == "D_3__Nitrosomonadales"|Rank6=="D_5__Candidatus Nitrososphaera")


# FunctionNiBin <- subset_samples(NitrifersOTU, Nitrification == "Low"|Nitrification =="High")
# FunctionNitrosoBin <- subset_samples(NitrosoOTU, Nitrification == "Low"|Nitrification =="High")

plot_bar(NitrifersOTU, "Pedigree", fill = "Rank3", facet_grid = ~Rank2)+geom_bar(stat="identity")+labs(title = "nitrifier taxa")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()
plot_bar(NitrifersOTU, "Pedigree", fill = "Rank4")+geom_bar(stat="identity")+labs(title = "nitrifier taxa")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()


OTU_merge = merge_samples(NILOTU, "Pedigree", fun = mean) #averages each OTU in all samples belonging to each habitat class
#This line changed in code ====
sample_data(OTU_merge)$Pedigree <- levels(sample_names(NILOTU)$Pedigree)
OTU_rel_merge = transform_sample_counts(OTU_merge, function(x) 100 * x/sum(x) )
NitrifersOTU=subset_taxa(OTU_rel_merge, Rank2 == "D_1__Nitrospirae"|Rank4 == "D_3__Nitrosomonadales"|Rank6=="D_5__Candidatus Nitrososphaera")

GP.ord <- ordinate(NILOTU, "NMDS", "bray",)
p1 = plot_ordination(NILOTU, GP.ord, type="Samples", color="Pedigree", title="Bacterial 16S Ordination")+geom_point(size=3)+stat_ellipse()
print(p1)


plot_bar(OTU_rel_merge, "Pedigree", fill = "Rank3")+geom_bar(stat="identity")+theme(legend.position = "none")
plot_bar(NitrifersOTU, "Pedigree", fill = "Rank4")+geom_bar(stat="identity")+labs(title = "nitrifier taxa")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()

GP.ord <- ordinate(OTU_rel_merge, "RDA", "bray",)
p1 = plot_ordination(OTU_rel_merge, GP.ord, type="Samples", color="Pedigree", title="Means")+geom_point(size=3)
print(p1)


GP.ord <- ordinate(NitrifersOTU, "RDA", "bray",)
p1 = plot_ordination(NitrifersOTU, GP.ord, type="Samples", color="Pedigree", title="Nitrfiers")+geom_point(size=3)
print(p1)


MM20156x610
GP.ord <- ordinate(MM20156x610, "NMDS", "bray",)
p1 = plot_ordination(MM20156x610, GP.ord, type="Samples", color="Pedigree", title="Sig Taxa")+geom_point(size=3)
print(p1)

##From the rhizosphere data we can see that there is a difference in the B73 microbiome vs NIL
#Need to determine list of OTUs different from the RHizosphere microbiome of Maize and teosinte in the field 

#Test Example
NILOTUB73.Teo <- subset_samples(NILOTU,Pedigree== "B73"| Pedigree=="PI384071")


library("DESeq2")

diagdds = phyloseq_to_deseq2(NILOTUB73.Teo, ~ Pedigree)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(NILOTUB73.Teo)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
# [1] 424  13

ggplot(sigtab, aes(x=Rank3, y=log2FoldChange, color=Rank4)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

TopNOTUs<-c(rownames(sigtab))

#Here we are getting the enriched list of taxa to subset them- them 
MM20156x610 = prune_taxa(TopNOTUs, NILOTU)
##This commented out code will be used if I want to plot the Raw abudannce.
#MM20156x610 = prune_taxa(TopNOTUs, TvM.OTU)

MM20156x610 = subset_taxa(MM20156x610, Rank3 != "NA")


##I orginally plotted them in abudnace now I will plot them in Relvativee abudnance
##Here I am pulling the signifcant OTUS and using them in a stack plot
plot_bar(MM20156x610, "Pedigree", fill = "Rank3")+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(legend.position = "none",strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


plot_bar(MM20156x610, "Pedigree", fill = "Rank3")+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


NitrifersOTU=subset_taxa(MM20156x610, Rank2 == "D_1__Nitrospirae"|Rank4 == "D_3__Nitrosomonadales"|Rank6=="D_5__Candidatus Nitrososphaera")
NILOTUB73.Teo <- subset_samples(NitrifersOTU,Pedigree== "B73"| Pedigree=="PI384071" | Pedigree=="Z031E0021" | Pedigree=="Z031E0021" |Pedigree=="B73xPI384071"| Pedigree=="Z031E0047")


plot_bar(NILOTUB73.Teo, "Pedigree", fill = "Rank4", title="Signficant Nitrfier OTUs")+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance (rarefied reads)", x="Treatment")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##The rhizosphere data is pretty clear cut, I think if we need to show a specifc rhizosphere relationship we should use rhizosphere data
##The rhizopsphere data also shows that B73 has Nitrfier OTUs comparted to 