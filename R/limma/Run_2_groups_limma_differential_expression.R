#!/usr/bin/env Rscript

###
###  George Bell
###  Bioinformatics and Research Computing, Whitehead Institute
###
###  Automate 2-group comparisons using limma.
###  Inputs: un-log-transformed matrix and design file
###
###  Version 1.0: 16 November 2020
###  Version 1.1: 6 February 2023 -- Fix correction of contrast matrix [GB]
###

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./Run_2_groups_limma_differential_expression.R" }

# Statistics require 'limma' package
# See limma Users Guide for more info
# http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
library(limma)

matrix.of.values = commandArgs()[6 - offset]
design.file = commandArgs()[7 - offset]
pseudocounts = as.numeric(commandArgs()[8 - offset])
output.file = commandArgs()[9 - offset]

if ( ! is.na(output.file) ) {
	message(paste("Using matrix file *", matrix.of.values, "* for input.", sep=""))
	message(paste("Design file: *", design.file, "*", sep=""))
	message(paste("Pseudocounts for log-transformation: *", pseudocounts, "*", sep=""))
	message(paste("Output file: *", output.file, "*", sep=""))
} else {
	message("\nDo 2-group limma analysis,")
	message("  starting with a matrix of pre-normalized values (NOT log-transformed).")
	message("\nUSAGE: ", this.script, " inputMatrix designFile pseudocounts outputFile\n")
	quit()
}

# Read matrix of values
values = read.delim(matrix.of.values, row.names=1)
message("\nAdding pseudocounts (", pseudocounts, ") and log2-transforming matrix ...")
log2.values = log(values + pseudocounts, 2)

# Read design as a file
design = read.delim(design.file, row.names=1)
design.samples = rownames(design)
# design.samples
# colnames(values)

# Check that columns in the matrix file match rows in the design file
if (! identical(colnames(values), design.samples))
{
	message("\nUh oh -- the column names of your matrix (after \"cleaning up\" by R)")
	print(colnames(values))
	message("don't seem to match the sample names in the first column of your design file")
	print(design.samples)
	message()
	quit()
}

# Convert design into a design matrix
design = design[,1]
design.matrix = model.matrix(~0+design)
rownames(design.matrix) = design.samples
# Make nicer design matrix group names
colnames(design.matrix) = gsub("design", "", colnames(design.matrix))

# Fit linear model for each feature given a series of samples (normalized)
fit = lmFit(log2.values, design.matrix)

groups = unique(as.vector(design))
num.groups = length(groups)
comparisons = vector(length = 2)
comparison.num = 1
for (i in 1:2)
{
	for (j in 1:2)
	{
		if (i != j)
		{
			comparisons[comparison.num] = paste(groups[i], "-", groups[j])
			comparison.num = comparison.num + 1
		}
	}
}

# GB -- 7 Nov 2022
message("Comparisons: ")
cat(paste(comparisons, collapse="\n"), "\n")

contrast.matrix = matrix(data=NA, nrow=2, ncol=2)
rownames(contrast.matrix) = colnames(fit$coefficients)
colnames(contrast.matrix) = comparisons
message("Contrast matrix (between groups, in both directions):")
for (i in 1:length(comparisons))
{
	# This code chunk from Saroj Mohapatra (https://support.bioconductor.org/p/27900/)
	prestr="makeContrasts("
	# poststr=",levels=design)"
	# poststr=",levels=unique(design))"
	# Need this sorting to work -- added 6 Feb 2023
	poststr=",levels=sort(unique(design)))"
	commandstr=paste(prestr,comparisons[i],poststr,sep="")
	contrast.matrix[,i] = eval(parse(text=commandstr))
}
contrast.matrix
message("WARNING: Check the correctness of the contrast matrix !!!")

fit.contrasts = contrasts.fit(fit, contrast.matrix)

# Compute moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
fit.ebayes = eBayes(fit.contrasts)

# topTable(fit.ebayes)
write.fit(fit.ebayes, file=output.file, digits=8, adjust="fdr")

message("\nAll done (!) -- see ", output.file, "\n")

message("Output file contents (where you can treat each \"-\" as \"/\"):
A             mean (log2-transformed) of all samples
Coef          log2 fold change
t             t-statistic [can be ignored]
p.value       raw p-value for moderated t-test
p.value.adj   adjusted p-value (FDR)
F             F statistic [can be ignored]
F.p.value     p-value reflecting whether the feature shows a difference in any comparison\n")

