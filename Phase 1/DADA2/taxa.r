

#### Load packages ####
library(dada2); packageVersion("dada2")
library(ggplot2)

setwd("/home/a-m/favela3/Research/NIL2018/V4_515F_New_V4_806R_New")

seqtab.nochim <- readRDS("results/dada2/seqtab_final.RDS")


### Assign taxonomy ####
# We assign taxonomy in two steps using the Silva 138.1 database. This is the most
# computationally intensive step, and requires a decent amount of memory/RAM. There's
# good chance this won't run on a laptop for a full data set.

# Step 1: Assign taxonomy through genus
taxa <- assignTaxonomy(seqtab.nochim, 
                       "/home/a-m/favela3/Research/SlivaDB/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=6
                       )

# Step 2: Assign species if 100% match
taxa <- addSpecies(taxa, 
                   "/home/a-m/favela3/Research/SlivaDB/silva_species_assignment_v138.1.fa.gz"
                   )

# Preview of results
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, file = "results/dada2/taxonomy_final.RDS")


