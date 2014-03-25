#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly=T);
if (length(args)!=3) {
  cat("USAGE: hicseq_cor_boxplot INPUT METHOD={Spearman|Pearson} OUTPUT\n");
  quit(save="no");
}

in_file <- args[1]
cmethod <- args[2]
out_file <- args[3]

HiC_cor_data <- read.table(in_file, header=F, sep="\t")

# Set the boxplots to appear in the order you want
HiC_1 <- "HindIII_1_2_flt"
HiC_2 <- "NcoI_1_2_flt"
HiC_3 <- "HindIII_NcoI_1_flt"
HiC_4 <- "HindIII_NcoI_2_flt"
HiC_5 <- "HindIII_1_2_cor"
HiC_6 <- "NcoI_1_2_cor"
HiC_7 <- "HindIII_NcoI_1_cor"
HiC_8 <- "HindIII_NcoI_2_cor"

labels <- paste(c(HiC_1,HiC_2,HiC_3,HiC_4,HiC_5,HiC_6,HiC_7,HiC_8))
HiC_cor_data$V3 = factor(HiC_cor_data$V3, c(HiC_1,HiC_2,HiC_3,HiC_4,HiC_5,HiC_6,HiC_7,HiC_8))
colors=c(rep("dark grey",4),rep("light grey",4))

#Output the resulting combined boxplot in .pdf format
pdf(paste(out_file,".pdf",sep=""))
boxplot(HiC_cor_data$V2 ~ HiC_cor_data$V3, col=colors, las=2, ylim=c(0,1), cex.axis=0.6, ylab=paste(cmethod,"correlation coefficient"), xlab="")
legend(x=6.5,y=0.2,legend=c("filtered","corrected"),fill=c("dark grey","light grey"))
dev.off()
