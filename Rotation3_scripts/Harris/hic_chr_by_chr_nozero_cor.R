#!/usr/local/bin/Rscript

# This is an R script that is used
# for comparison (find correlation
# between Hi-C contact matrices). This
# script outputs only the summary for 
# the resulting correlation vector to
# the standard output (console). Also 
# the lines and columns with zeros are 
# removed from the input matrices.

args <- commandArgs(trailingOnly=T);
if (length(args)!=6) {
  cat("USAGE: hicseq_chr_by_chr_nozero_cor MATRIX1 MATRIX2 CHROMOSOME-VECTOR VALUE METHOD={pearson,spearman} HEADER={true,false}\n");
  quit(save="no");
}

m1_file <- args[1];
m2_file <- args[2];
fc <- args[3];
value <- args[4];
cmethod <- args[5];
header <- args[6];

if (header=="true") {
	m1 <- as.matrix(read.table(m1_file, header=T, row.names=1));
	m2 <- as.matrix(read.table(m2_file, header=T, row.names=1)); 
} else if (header=="false") {
	m1 <- as.matrix(read.table(m1_file, header=F));
    m2 <- as.matrix(read.table(m2_file, header=F));
} else {
	cat ("Error: Use true/false for header\n");
	quit(save="no");
}
 
# Get rid of the rows and columns that sum to zero
a <- apply(m1,1,sum) > value
b <- apply(m2,1,sum) > value
w <- a&b

# Get the new matrices that have
# no rows or columns that sum to zero
m1 <- m1[w,w]
m2 <- m2[w,w]

# Get the unique chromosomes
chrom <- as.matrix(read.table(fc));
chrom <- chrom[w]

chr_num <- max(chrom);

# Create a vector to store the correlations
cor_vector <- rep(0, chr_num)

# Set the chromosome number in a variable and
# calculate Hi-C contact correlation for each
# chromosome and all chromosomes.
# Set the number of chromosome in a variable:
c <- 1:chr_num

for (i in c) {
	I <- chrom==i
	cor_vector[i] <- cor(as.vector(m1[I,I]),as.vector(m2[I,I]),method=cmethod)
}

# Report the summary of the correlation 
# vector (quantiles, mean, median).

# Report also the summary
cat("Correlation summary")
cat("\n")
print(summary(cor_vector))
cat("\n")

