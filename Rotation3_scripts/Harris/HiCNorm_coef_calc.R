# R script to calculate 
# Pearson and Spearman correlation
# coefficients based on the normalized
# Lieberman-Aiden data from the paper:
# Hu et. al "HiCNorm: removing biases in Hi-C data via Poisson regression"
# Bioinformatics 28, 23, 2012, pp 3131-3133
# Data available on http://www.people.fas.harvard.edu/~junliu/HiCNorm
# Date: Mar 12 2014
# Harris A. Lazaris

# Import the normalized HindIII and NcoI Lieberman-Aiden data (normalized 
# using HiCNorm).

setwd("/Users/chlazaris/Desktop/Sackler_PhD/ROTATIONS/ROTATION_3/Data/20140312_HiCNorm_data/Lieberman-Aiden")

HindIII_norm <- as.vector(read.delim("normalized_whole_genome_contact_map_HindIII.txt", header = TRUE,  sep = "\t",  stringsAsFactors = FALSE)) 
NcoI_norm <- as.vector(read.delim("normalized_whole_genome_contact_map_NcoI.txt", header = TRUE,  sep = "\t",  stringsAsFactors = FALSE))

# Question: What is this extra bin that appears before bin 1?

A <- data.matrix(HindIII_norm, rownames.force = NA)
B <- data.matrix(NcoI_norm, rownames.force = NA)

pearson_coef <- cor(A, B, use="everything",method=c("pearson"))
spearman_coef <- cor(A, B, use="everything",method=c("spearman"))

# In this case we get one coefficient for one bin. What we want though
# is one coefficient for each chromosome.