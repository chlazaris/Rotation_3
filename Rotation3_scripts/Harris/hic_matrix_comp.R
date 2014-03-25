#!/usr/local/bin/Rscript

# This is an R script that is used
# for comparison (find correlation
# between Hi-C contact matrices
# and output boxplot with the results)

args <- commandArgs(trailingOnly=T);
if (length(args)!=5) {
  cat("USAGE: hicseq_matrix_comp MATRIX1 MATRIX2 CHROMOSOME-VECTOR METHOD={pearson,spearman} HEADER={true,false}\n");
  quit(save="no");
}

m1_file <- args[1];
m2_file <- args[2];
fc <- args[3];
cmethod <- args[4];
header <- args[5];
#of <- args[6] # This is the output file with name specified by user


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
#chrom <- as.matrix(read.table(fc));

#chr_num <- max(chrom);

# Create a vector to store the correlations
#cor_vector <- rep(0, chr_num)

# Set the chromosome number in a variable and
# calculate Hi-C contact correlation for each
# chromosome and all chromosomes.
# Set the number of chromosome in a variable:
#c i <- 1:chr_num

#for (i in c) {
#	I <- chrom==i
#	cor_vector[i] <- as.vector(m1[I,I])
#}

#Create the boxplot to summarise the results
#pdf(paste(args[1],args[2],".pdf",sep="_"))
#boxplot(cor_vector, xlab = "HiCNorm", ylab = paste(cmethod,"correlation coefficients"))
#dev.off()

#The output will be a file containing a matrix with three columns
#chrom number, correlation value, output file
#Given that there are 23 chromosomes, the matrix will be 23*3

# Report also the summary
#cat (mean(cor_vector))
#cat ("\n")
#filename_vector <- rep(of, chr_num) 
#result_matrix <- cbind(c, cor_vector, filename_vector)

# Write the resulting matrix to a file
#write.table(result_matrix, file=of, sep="\t", quote=F, row.names=F, 
#	col.names=F)

#cat("Mean correlation:\n")
#mean_cor <- mean(cor_vector)
#print(mean_cor)
#cat("Median correlation:\n")
#median_cor <- median(cor_vector)
#print(median_cor)

total_cor <- cor(as.vector(m1),as.vector(m2),method=cmethod)
cat("The total correlation of the matrices is:")
print(total_cor)


