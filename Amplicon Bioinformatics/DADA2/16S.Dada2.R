#Mcrobial Profiling Workshop ###
# These codes are for a 16S paired-end data set sequenced on an Illumina MiSeq platform. 
# Reads are 300nt, and from the V4 region


#### Load packages ####
library(dada2); packageVersion("dada2")
library(ggplot2)


#### Set the correct working directory ####
setwd("/home/a-m/favela3/Research/NIL2018/V4_515F_New_V4_806R_New")


#### Import files into dada2 format ####

# Confirm that the FASTQ files are in the correct location
path <- "/home/a-m/favela3/Research/NIL2018/V4_515F_New_V4_806R_New"
list.files(path)

message(" Sort paired-end reads")

# Remember, forward and reverse fastq filenames have format SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnF <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnR <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_R1.fastq
sample.names <- sub("_R1.fastq", "", basename(fnF))



#### Quality Control ####

# Inspect read quality profiles using dada2.
# Note: We already did this using FASTQC and MultiQC, but let's see what the dada2 version looks like
#plotQualityProfile(fnF, aggregate = TRUE)
#ggsave("forward-quality-profile.pdf")

#plotQualityProfile(fnR, aggregate = TRUE)
#ggsave("Reverse-quality-profile.pdf")


# What happened to the plot?
#?plotQualityProfile

# We can use ggsave() to save it to a specific place and filename
# ggsave() acts on the last ggplot made
#ggsave("forward-quality-profile.pdf")


#### Exercise 1 ####
# Create a quality profile image of the reverse reads using dada2, and save it using ggsave().

#ggsave("forward-quality-profile.pdf")


message("### Filtering the sequences ####")

# Create filtered filenames and paths (in filtered/ subdirectory)
dir.create("results/dada2/filtered")
filtF <- file.path("results/dada2/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtR <- file.path("results/dada2/filtered", paste0(sample.names, "_R_filt.fastq.gz"))

message("# Assign names to these file paths")
names(filtF) <- sample.names
names(filtR) <- sample.names

# Perform filtering:
# This is the trickiest step because you need to consider multiple attributes of your data.
# (1) You'll want to trim from the left of your sequences the length of your primers. See
# your sequencing report for this information. Then (2) from the right of your reads you'll want
# to trim low quality data plotted in the previous visualization (or FASTQC/MultiQC) and (3) anything 
# in excess of your targeted regions length, because beyond that can only be primers. In our 
# case, the V4 region is 252nt. It helps to draw a visualization of this.

# Let's go through each part of this command before we start it. 
## First, go back to the slide deck and I'll explain these ideas in more detail.
## Second, run the following help page in Rstudio:

library(dada2) #If you haven't run it in Rstudio yet
#?filterAndTrim

# NOTE: The following line of code does not work in R 4.0.3 because of bug that was fixed in the
# more recent version of R and Bioconductor. This will work on the Desktop computers because they
# have R 4.0.5. Biocluster will be updated in a week when the new Bioconductor comes out.
message("######filerAndTrim#######")
##Need to run this first 
out <- filterAndTrim(fnF, filtF, 
                      fnR, filtR, 
                      trimLeft=c(19,20), 
                      truncLen=c(240,230),
                      maxN=0, 
                      maxEE=c(2,2),
                      rm.phix=TRUE,
                     multithread=6) 
# Instead of running the above code. Load in the output from the above code here:
#out <- readRDS("results/dada2/filter-and-trim-output.RDS")
#head(out)

saveRDS(out, file = "results/dada2/filter-and-trim-output.RDS")

# NOTE: TruncLen is not desirable for ITS reads, since they vary in length.





###Learning Error Rate to use in model for correction=====

message("#### Learn Error Rates ####")
# This step learns the error rates of your sequences for more accurate sample inference
# in the next step.
# This step is the second most computationally intensive. This may not run on your 
# personal computer. 
# NOTE: On Windows set multithread=FALSE.
errF <- learnErrors(filtF, multithread=6)
errR <- learnErrors(filtR, multithread=6)

##This script takes a while so I saved it just in case 
saveRDS(errF, file = "results/dada2/errF.RDS")
saveRDS(errR, file = "results/dada2/errR.RDS")


# Plot errors to make sure everything looks as expected. 
# If you used the NovaSeq platform,
# and the error distribution has large dips in comparison to the theoretical red line,
# then check the recommendations for this here:  

plotErrors(errF, nominalQ=TRUE)
ggsave("results/dada2/forward-error-distribution.pdf")
plotErrors(errR, nominalQ=TRUE)
ggsave("results/dada2/reverse-error-distribution.pdf")








#### Dereplicate samples ####
# Dereplication combines all identical reads into into “unique sequences” with a 
# corresponding “abundance”: the number of reads with that unique sequence. This reduces 
# computation time by eliminating redundant comparisons. This set could require a lot of 
# RAM.


#### Exercise 2 ####
# Look at the manual for derepFastq() and try to dereplicate the filtF and filtR objects 
# (these contain file pathes to our filtered fastqs) and use the verbose mode.
derepF <- derepFastq(filtF, verbose = TRUE)
derepR <- derepFastq(filtR, verbose = TRUE)
# If you had PacBio reads, what change would you need to make to the above code?

  
  
message('#### Sample Inference ####')
# This is the core algorithm of dada2 that will infer real sequence variants from our
# unique sequences using the error models we generated to guide its decisions. See the 
# dada2 publication for more details: https://doi.org/10.1038/nmeth.3869
dadaF <- dada(derepF, err=errF, multithread=6)
dadaR <- dada(derepR, err=errR, multithread=6)






#### Merge paired reads ####
# This is where we'll finally merge the paired-end reads together into one read, which will
# further improve accuracy. Note that any unmerged reads will be lost. It is possible to rescue 
# these unmerged reads (see the returnRejects option), but proceed with caution. Also note that 
# it is possible to artificially merge  reads that don't overlap using the justConcatenate option.

## See slides for visuals and explanations of PE merging. ##

#### Exercise 3 ####
# Set the minimum overlap of our PE reads to 200 in the command below
mergers <- mergePairs(
  dadaF, filtF, 
  dadaR, filtR, 
  verbose=TRUE
)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Make a sequence table.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths (should be around size of V4 region)
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, file = "results/dada2/seqtab-withchim.RDS")


#### Remove chimeras. ####
# Chimeras are a mistake that happens during PCR where to different sequences merge 
# together to form a hybrid sequence. This are easier to identify after all filtering has 
# happened.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=6, verbose=TRUE)
dim(seqtab.nochim)

# The number of chimeras made up 21% of the inferred sequence variants. This is not unusual,
# and sometimes they might make up a majority. Let's check the percentage of sequences that
# this makes up...
sum(seqtab.nochim)/sum(seqtab)

# This was about 11% of total sequences. Not terrible, but we wouldn't expect chimeras to
# take up too much more than this. If too many are removed, the most common reason is that 
# you haven't properly removed primer sequences in the filtering step.


#### Summarize results ####
# Create a summary table of sequences lost at each filtering step. Confirm that there 
# was a gradual loss of sequencing reads at each step. If there was a drastic loss,
# consider adjusting parameters at that step.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# We did lose anywhere from 25-40% of sequencing reads depending on the sample. This
# would be more concerning if this weren't an older dataset (2014) using less accurate
# lab protocols and sequencing.
# If this were a modern data set, we'd want to play around with the filtering and merging
# steps. See the slides for possible reasons and solutions for too much data loss.


#### Save final data ####
saveRDS(track, file = "results/dada2/summary-of-dada2.RDS")
saveRDS(seqtab.nochim, file = "results/dada2/seqtab_final.RDS")


