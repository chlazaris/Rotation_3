# R code to import all .txt files
# containing the cis-contact data 
# from Lieberman-Aiden and combine
# them into one. Then the script
# generates the matrix for the whole
# genome (populating the rest of the
# matrix (off-diagonal) with zeros).
# The script also outputs the chromosome
# vector to be used for the creation of 
# boxplots.

# Author: Harris Lazaris
# Date: March 28 2014

# Read in the files in the
# right order
temp <- c(list.files(pattern="*chr[0-9].txt"),list.files(pattern="*chr[0-9][0-9].txt"))
myfiles <- lapply(temp, read.table)

# Chromosome vector
c <- 1:length(myfiles)

# Now access each element of the myfiles list
# which correspond to each chromosome and calculate
# its dimensions. Write these dimensions to a matrix
dim_matrix <- matrix(rep(0,max(c)*3),nrow=max(c),ncol=2)

# Generate the chromosome vector
chrom_vector <- vector()

for (i in c) {
	dim_matrix[i,1] <- i
	dim_matrix[i,2] <- dim(myfiles[[i]])[1]
	dim_matrix[i,3] <- dim(myfiles[[i]])[2]
}

# Populate the chromosome vector
for (i in 1:nrow(dim_matrix)) {
	v <- c(v,rep(dim_matrix[i,1],dim_matrix[i,2]))
}


genome_matrix_dim <- colSums(dim_matrix)

# Get the number of rows and columns
# for the genome matrix
genmat_rows <- genome_matrix_dim[2] 
genmat_cols <- genome_matrix_dim[3]

# Construct the matrix
# Populate with zeros
genome_matrix <- matrix(rep(0, genmat_rows*genmat_cols), nrow=genmat_rows, ncol=genmat_cols)

# Go through each chromosome and store
# the cis contact data to the right place
# in the genome matrix

r1 <- 1
c1 <- 1

for (i in 1) {
	# Get the matrix that corresponds
	# to the chromosome and find dimensions
    chrom_rows <- dim(myfiles[[i]])[1]
    chrom_cols <- dim(myfiles[[i]])[2]
	# Columns and rows have equal 
	# dimensions
    r2 <- r1 + chrom_rows - 1
    c2 <- c1 + chrom_cols - 1

    # Put the chromosome data in the matrix
    genome_matrix[r1:r2,c1:c2] <- as.matrix(myfiles[[i]])

    # Move to the next one
    r1 <- r2 + 1
    c1 <- c2 + 1
}

# Now write both the genome matrix and chromosome vector 
# to files
write.table(chrom_vector, file="chrom_vector.txt", row.names=F, col.names=F)
write.table(genome_matrix, file="genome_matrix.txt", row.names=F, col.names=F)

