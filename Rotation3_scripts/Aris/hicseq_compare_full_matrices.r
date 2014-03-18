#!/usr/local/bin/Rscript

##
## USAGE: hicseq_compare_matrices MATRIX-SET-1 MATRIX-SET-2 CHROMOSOME-VECTOR METHOD={pearson,spearman} 
##

args <- commandArgs(trailingOnly=T);
if (length(args)!=4) {
  cat("USAGE: hicseq_compare_matrices MATRIX-SET-1 MATRIX-SET-2 CHROMOSOME-VECTOR METHOD={pearson,spearman}\n");
  quit(save="no");
}

fa <- strsplit(args[1],' ')[[1]];
fb <- strsplit(args[2],' ')[[1]];
fc <- args[3];
cmethod <- args[4];

n_a <- length(fa);
n_b <- length(fb);
n <- n_a + n_b;

chr <- as.matrix(read.table(fc));
nn <- nrow(chr);

M <- NULL;
Asum <- matrix(0,nn,nn);
Bsum <- matrix(0,nn,nn);
for (i in 1:n_a) { cat('reading',fa[i],'...\n'); M[[i]] <- as.matrix(read.table(fa[i])); Asum <- Asum + M[[i]]; }
for (i in 1:n_b) { cat('reading',fb[i],'...\n'); M[[n_a+i]] <- as.matrix(read.table(fb[i])); Bsum <- Bsum + M[[n_a+i]]; }

# correlations
mask <- matrix(FALSE,nn,nn);
for (i in 1:nn) { for (j in i:nn) { if (chr[i]==chr[j]) { mask[i,j] <- mask[j,i] <- TRUE; }; }; }

sum_global_cor <- cor(as.vector(Asum),as.vector(Bsum),method=cmethod);
sum_intra_cor <- cor(Asum[mask],Bsum[mask],method=cmethod);
sum_inter_cor <- cor(Asum[!mask],Bsum[!mask],method=cmethod);

intra_cor <- matrix(0,n,n);
inter_cor <- matrix(0,n,n);
for (i in 1:(n-1)) {
  intra_cor[i,i] <- 1.0;
  inter_cor[i,i] <- 1.0;
  for (j in (i+1):n) {
    intra_cor[i,j] <- intra_cor[j,i] <- cor(M[[i]][mask],M[[j]][mask],method=cmethod);
    inter_cor[i,j] <- inter_cor[j,i] <- cor(M[[i]][!mask],M[[j]][!mask],method=cmethod);
  }
}
intra_cor[n,n] <- 1.0;
inter_cor[n,n] <- 1.0;

cat("Global correlation (sum):\n");
round(sum_global_cor,3)
cat("Intra-chromosomal correlation (sum):\n");
round(sum_intra_cor,3)
cat("Inter-chromosomal correlation (sum):\n");
round(sum_inter_cor,3)

cat("Intra-chromosomal correlation:\n");
round(intra_cor,3)
cat("Inter-chromosomal correlation:\n");
round(inter_cor,3)

#d <- row(a)-col(a); cor(a[d>0],b[d>0]);


