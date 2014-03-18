#/usr/local/bin/Rscript

# This is an R script that is used
# for comparison (find correlation
# between Hi-C contact matrices
# and output boxplot with the results)

args <- commandArgs(trailingOnly=T);
if (length(args)!=5) {
  cat("USAGE: hicseq_matrix_compare MATRIX1 MATRIX2 CHROMOSOME-VECTOR METHOD={pearson,spearman} HEADER={true,false}\n");
  quit(save="no");
}

m1_file <- args[1];
m2_file <- args[2];
fc <- args[3];
cmethod <- args[4];
header <- args[5];


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
 

# Get the unique chromosomes
chrom <- as.matrix(read.table(fc));

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

#Create the boxplot to summarise the results
pdf(paste(args[1],args[2],".pdf",sep="_"))
boxplot(cor_vector, xlab = "HiCNorm", ylab = paste(cmethod,"correlation coefficients"))
dev.off()
# Report also the summary
cat (mean(cor_vector))
cat ("\n")