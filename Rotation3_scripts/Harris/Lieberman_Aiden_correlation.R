setwd("~/Desktop/20140314_Dekker_cor_analysis/Lieberman-Aiden")

HindIII_norm_feat <- read.delim("local_genomic_features_HindIII.txt", header=T, sep="\t")
NcoI_norm_feat    <- read.delim("local_genomic_features_NcoI.txt", header=T, sep="\t")

# Given that the matrix files downloaded from Liu, had additional tabs, vi was used
# to get rid of the extra tabs
# vi + file_name and then: :%s/\t\t/\t/g
# Then and by using row.names=1 the matrices appear in the right format

HindIII_norm_contact_map <- read.delim("normalized_whole_genome_contact_map_HindIII.txt", header=T, row.names=1, sep="\t")
NcoI_norm_contact_map    <- read.delim("normalized_whole_genome_contact_map_NcoI.txt", header=T, row.names=1, sep="\t")

# Now add the chromosome column from the local_genomic_features_NcoI and local_genomic_features_HindIII
# in order to split based on the chromosome number

#test_HindIII <- cbind(HindIII_norm_contact_map,HindIII_norm_feat$chr)
#test_NcoI    <- cbind(NcoI_norm_contact_map,NcoI_norm_feat$chr)

# Then split those based on the chrom numbers and create the corresponding lists
#test_HindIII_split <- split(test_HindIII, HindIII_norm_feat$chr)
#test_NcoI_split    <- split(test_NcoI, NcoI_norm_feat$chr)

# Given that the two lists (test_HindIII_split and test_NcoI_split)
# are of the same size, loop over and calculate the Spearman correlation
# for each one of the chromosomes

# Chromosome information
chrom <- HindIII_norm_feat$chr

# Create a vector to store the correlations
chrom_num <- length(unique(chrom))
cor_vector <- rep(0,chrom_num)

# Set the number of chromosome in a variable:
c <- 1:chrom_num

for (i in c) {
	I <- chrom==i
	cor_vector[i] <- cor(as.vector(as.matrix(HindIII_norm_contact_map[I,I])),as.vector(as.matrix(NcoI_norm_contact_map[I,I])),method='spearman')
}

#Create the boxplot to summarise the results
boxplot(cor_vector, xlab = "HiCNorm", ylab = "Spearman correlation coefficients")
# Do also the summary
summary(cor_vector)

# To calculate correlation for the whole HindIII and NcoI genomes
# whole_genome_cor <- cor(as.vector(as.matrix(HindIII_norm_contact_map)),as.vector(as.matrix(NcoI_norm_contact_map)), method="spearman")
# Result
# 0.07944521
# Check the datasets



