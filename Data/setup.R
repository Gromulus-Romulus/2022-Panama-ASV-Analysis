# Setup Script for Panama ASV analysis.
#
# Output: .rds file containing phyloseq object, microbial taxa,
# ASV table, sample metadata, and guild classifications.
#
# Side Effects:
#   Will rename the first column
#   of the metadata and ASV tables to "ASV_ID".
#
# Author: Nathan Malamud.

library(phyloseq)
library(dplyr)
library(reshape2)
library(spaa)
library(tidyverse)

# Define source files here
ASV.F = "./PanamaMicSig_ITS_fwdASVtabVT.txt" # ASV table
META.F = "./Panama_micsig_metadata.csv" # Sample metadata
GUILD.F = "./PanamaMicSig_fwd_funguild.guilds.txt" # FunGUILD output


# TODO - unknown columns error
build.phyloseq <- function (metadata, asvtab) {
  # Handles all of the grunt work for building a phyloseq object.
  #
  # Parameters: metadata - sample metadata, asvtab - asv table
  # Returns: physeq.obj - a phyloseq object
  
  # Necessary data wrangling
  map = asvtab[,-c(1:8)] 
  
  # TODO: Unknown columns in metadata but not asvtab
  # Remove these outside of the function.
  # `HITR_038CL`, `HITR_190HL`, `SIAM_026CL`, `SIAM_027CS`, `TEPA_005S`
  
  map2 = map %>% select(one_of(as.character(metadata$Sample_ID)))
  names = colnames(map2)
  
  map2 = cbind(asvtab[,1:8], map2) 
  map3 = map2[,9:249]
  
  map3.t = data.frame(t(map3))
  
  # Remove duplicate Sample IDs
  filter_meta = metadata %>% distinct(Sample_ID, .keep_all = TRUE)
  
  
  # Set our rownames and get the sample names
  # into the same order as the ASV table
  rownames(filter_meta) = filter_meta$Sample_ID
  imeta = filter_meta[rownames(map3.t),]
  
  # Convert to a matrix of microbial taxa
  # so that we can build our phyloseq object
  tax.matrix=as.matrix(map2[,1:8])
  
  # otu_table is equivalent to asv_table here, don't worry about it
  ASV = otu_table(map3, taxa_are_rows = TRUE)
  TAX = tax_table(tax.matrix)
  
  physeq.obj = phyloseq(ASV,TAX)
  samples = sample_data(as.data.frame(imeta,
                                      row.names=sample_names(physeq.obj),
                                      stringsAsFactors=FALSE))
  
  # Almost done, just need to prune it, and tada!
  physeq.obj = phyloseq(ASV,TAX,samples)
  physeq.obj <- prune_samples(sample_sums(physeq.obj) > 0, physeq.obj)
  
  return(physeq.obj)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - #
#   Driver Code   #
# - - - - - - - - #

# Import asvtable, metadata, and FunGUILD output
asvtab = read.delim(ASV.F)
metadata <- read.delim(META.F, sep=',')
guilds <- read.delim(GUILD.F)

# Reset column names, originally called "X"
colnames(guilds)[1] <- "ASV_ID"
colnames(asvtab)[1] <- "ASV_ID"

# Build phyloseq object for ITS sequences
itsphyseq <- build.phyloseq(metadata, asvtab)

# Package in separate .rds files - to be accessed in other scripts
saveRDS(asvtab, './asvtab.rds')
saveRDS(metadata, './metadata.rds')
saveRDS(guilds, './guilds.rds')
saveRDS(itsphyseq, './itsphyseq.rds')


