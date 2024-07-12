#! /usr/bin/env Rscript

# This wrapper allows easy command line access to removing batch effects
# from counts data using the ComBat_seq function from the Surrogate
# Variable Analysis (sva) package.

# The expected format for the input counts matrix is samples x genes
# (columns x rows) with names for each sample.

# Also taken as input is a sample annotation with batch and condition listed.
# The samples listed in both the counts matrix and the annotation file must
# match.

# Version 1.0, 18 May, 2021
# Troy Whitfield - Bioinformatics and Research Computing, Whitehead Institute

argv <- commandArgs()
if('--args' %in% argv){
   nignore <- which (argv == "--args")
   argv <- argv[-(1:nignore)]
}
if((length (argv) <= 4) | (argv[1] == "/usr/lib/R/bin/exec/R")) {
   message("Convert a raw counts matrix with batch and condition")
   message("information to a batch-adjusted counts matrix using")
   message("ComBat-seq.\n")
   message("USAGE: ./remove_batch_effects_counts.R counts samples covariates adjCounts\n")
   message("Example: ./remove_batch_effects_counts.R counts.txt samples.txt default adjCounts.txt\n")
   message("The two input files should be tab delimited text.")
   message("The counts matrix file should be samples x genes (columns x rows).")
   message("The annotation file should have samples listed in colunm 1,")
   message("batches listed in column 2 and biological condition listed in column 3.")
   message("The samples listed in the counts matrix should match those from the")
   message("sample annotations.\n")
   message("The 'covariates' keyword controls whether the condition annotations")
   message("from the annotation file should be used to explain some of the")
   message("between-sample variance in the counts matrix.  This keyword can be")
   message("set to 'true' or 'false', but the default is to include such features.\n")
# stop("Follow the expected input format above.")
   q()
}

message("\nAdjusting input counts matrix for batch effects using ComBat-seq.")
message("If you use these results in a publication, please cite")
message("Y. Zhang et al., NAR Genomics and Bioinformatics 2, lqaa078 (2020).\n")

suppressMessages(require(sva))

# Assign variables from input.
countsfile<-argv[1]
samplefile<-argv[2]
covariates<-argv[3]
outfile<-argv[4]

covKeys<-c("true","false","default")
if (covariates == "default"){
   covariates<-"true"
} else if (!(covariates %in% covKeys)){
   covariates<-"true"
}

# Read counts and annotations.
counts<-read.table(countsfile,header=TRUE,sep="\t",check.names=FALSE,row.names=1)
counts<-as.matrix(counts)
sample<-read.table(samplefile,header=TRUE,sep="\t")

# Make sure that sample names match between annotations and counts.
nsamps<-length(colnames(counts))
if (abs(length(which(colnames(counts) %in% sample$sample))-nsamps) != 0){
   message("Exiting. Please ensure that sample annotations match the counts matrix.")
   q()
}   

Batch<-sample[,2]

# Adjust for batch effects.
if(covariates == "true"){
    Group<-sample[,3]
    adjcounts<-ComBat_seq(counts, batch=Batch, group=Group, full_mod=TRUE)
} else {
    adjcounts<-ComBat_seq(counts, batch=Batch, group=NULL, full_mod=FALSE)
}

# Write the results to file.
out.table = cbind(rownames(adjcounts),adjcounts)
colnames(out.table)[1] = "feature.name"
write.table(out.table, outfile, row.names=FALSE, quote=FALSE,sep="\t")

message("Batch adjustment complete.\n")
