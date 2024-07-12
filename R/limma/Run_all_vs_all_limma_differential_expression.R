#!/usr/bin/env Rscript

###
###  George Bell
###  Bioinformatics and Research Computing, Whitehead Institute
###
###  Automate all vs.all comparisons using limma.
###  Inputs: un-log-transformed matrix and design file
###
###  Version 1.0: 4 October 2019
###  Version 1.1: 20 January 2023 -- Fix correction of contrast matrix [GB]
###  Version 1.2: 23 January 2023 -- Fix making-unique of feature IDs [GB]
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
if (is.na(this.script)) { this.script = "./Run_all_vs_all_limma_differential_expression.R" }

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
	message("Do an all vs. all limma analysis,")
	message("  starting with a matrix of pre-normalized values (NOT log-transformed).")
	message("\nUSAGE: ", this.script, " inputMatrix designFile pseudocounts outputFile\n")
	quit()
}

# Added make.unique.3 [GB -- 16 Dec 2022] and corrected it [23 Jan 2023]

make.unique.3 <- function(x){
    # From https://rdrr.io/github/vladpetyuk/vp.misc/src/R/make.unique.2.R
    # but then further modified by GB, Whitehead BaRC
    if(!any(duplicated(x)))
        return(x)
    x.dups = unique(x[duplicated(x)])
    # xu <- unique(x)
    x2 <- c(x.dups, x)
    x3 <- make.unique(x2)
    return(x3[-seq_along(x.dups)])
}

values = read.delim(matrix.of.values, header=T, sep="\t", check.names=FALSE)

message("Done reading input matrix.  Now processing ....")

# Get total values per gene (for sorting redundant rows)
values.sum.by.feature = apply(values[,2:ncol(values)], 1, sum)
# Order matrix by decreasing values.sum.by.gene, adding a numeric suffix to each symbol we've already seen.
values = values[order(values.sum.by.feature, decreasing=TRUE),]
feature.symbols.original = as.character(values[,1])
# feature.symbols.unique = .Internal(make.names(feature.symbols.original, unique=TRUE))
feature.symbols.unique = make.unique.3(feature.symbols.original)
o = order(feature.symbols.original != feature.symbols.unique, values.sum.by.feature, decreasing=TRUE)
feature.symbols.unique[o] = make.unique(feature.symbols.unique[o])
# Make row names from unique gene symbols.
rownames(values) = feature.symbols.unique
# Drop the gene symbol (first column) since we have them as row names
values = values[,-1]

values.new.file = gsub(".txt", ".unique.names.txt", matrix.of.values)
message("Printing out non-redundant matrix: ", values.new.file)
write.table(values, values.new.file, sep="\t", quote=F)

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

# Use just the first column of the design file (second column really, since first column should be sample IDs/names)
design = design[,1]
design.matrix = model.matrix(~0+design)
colnames(design.matrix) = gsub("design", "", colnames(design.matrix))
rownames(design.matrix) = design.samples
message("\nExperimental design:")
design.matrix

# Get all possible comparisons
groups = unique(as.vector(design))
num.groups = length(groups)
comparisons = vector(length = num.groups * num.groups - num.groups)
comparison.num = 1
for (i in 1:length(groups))
{
	for (j in 1:length(groups))
	{
		if (i != j)
		{
			comparisons[comparison.num] = paste(groups[i], "-", groups[j])
			comparison.num = comparison.num + 1
		}
	}
}

# Fit linear model for each gene given a series of arrays (normalized)
# Agilent 2-color arrays (with offset of 0) have been normalized by loess and then Aquantile.
fit = lmFit(log2.values, design.matrix)

# Optional: Ignore control spots in linear model
# Make a list of weights so negative and positive control spots have a weight of 0 and experimental spots have a weight of 1
# Get vector of control types (0 (experimental), -1 (neg control), 1 (pos control))
# controlType = MA.loess.q.0$genes[,"ControlType"]
# wt.noControls = as.numeric(controlType == 0)
# fit = lmFit(MA.loess.q.0, design, weights = wt.noControls)

# Describe the comparisons you want to make
# To find out:
# 1 - What is the effect of genotype (mutant A vs wt and mutant b vs wt)
#     for both untreated and drugged cells? (within array comparison)
# 2 - How does drug treatment change the effect of genotype? (between array comparison)
#

contrast.matrix = matrix(data=NA, nrow=num.groups, ncol=length(comparisons))
rownames(contrast.matrix) = colnames(fit$coefficients)
colnames(contrast.matrix) = comparisons
message("Contrast matrix (all possible combinations):")
for (i in 1:length(comparisons))
{
	# This code chunk from Saroj Mohapatra (https://support.bioconductor.org/p/27900/)
	prestr="makeContrasts("
	# poststr=",levels=design)"
	# Need this to work -- added 20 Jan 2023
	poststr=",levels=sort(unique(design)))"
	commandstr=paste(prestr,comparisons[i],poststr,sep="")
	contrast.matrix[,i] = eval(parse(text=commandstr))
}
contrast.matrix

# Compute Contrasts from Linear Model Fit
fit2 = contrasts.fit(fit, contrast.matrix)

# Compute moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
fit2.ebayes = eBayes(fit2)

# Write a microarray linear model fit to a file
# 'Coef' refers to log2 ratios
# Output will contain both raw p-values and FDR-corrected p-values
write.fit(fit2.ebayes, file=output.file, digits=8, adjust="fdr")

message("\nAll done (!) -- see ", output.file, "\n")

message("Output file contents (where you can treat each \"-\" as \"/\"):
A             mean (log2-transformed) of all samples
Coef          log2 fold change
t             t-statistic [can be ignored]
p.value       raw p-value for moderated t-test
p.value.adj   adjusted p-value (FDR)
F             F statistic [can be ignored]
F.p.value     p-value reflecting whether the feature shows a difference in any comparison\n")
