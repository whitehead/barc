#!/usr/bin/env Rscript

###
###  Given a matrix of feature counts, calculate correlation to a reference feature.
###  The reference feature should be the first feature in the inputCounts file.
###  The inputCounts table is log2-transformed (after adding a pseudocount) and median-centered
###    before correlations are calculated.
###
###  George Bell - BaRC
###
###  Input file: expression matrix (counts with one row per gene/transcript and a column per sample)
###  First column should be feature IDs
###
###  Version 0.1 - first version - 28 January 2020
###  Version 0.2 - 22 April 2020: Allow reference feature to be anywhere in the file.
###  Version 0.3 - 30 April 2020 - Add 'use="pairwise.complete.obs"' to cor.test() and add status messages
###

# For user interface choice
offset = 0
# For log-transformation
pseudocounts = 1
# For output table
significant.figures = 4

if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./calculate_correlation_to_reference_feature.R" }

if (length(commandArgs()) < (6 - offset))
{
	message("\nGiven a matrix of (positive, non-log-transformed) feature counts, calculate correlation to a reference feature.")
	message(paste("USAGE:   ", this.script, "inputCounts referenceID OutputFile"))
	message(paste("Example: ", this.script, "MySample.counts.txt GATA4 MySample.gene_correlations.txt\n"))
	message("The reference feature should be the first feature in the inputCounts file.")
	message("The inputCounts table is log2-transformed (after adding a pseudocount) and median-centered")
	message("  before correlations are calculated.\n")
	q()
}

# File of values, with first column of feature IDs
# The first feature is the reference feature
values.file = commandArgs()[6]
reference.feature = commandArgs()[7]
output.file = commandArgs()[8]

# Read the input matrix, with feature names as first column
message("Reading values.file ...")
values = read.delim(values.file, row.names=1)
message("Done reading.  Now log-transforming matrix ...")

# Log-transform
values.log2 = log(values + pseudocounts, 2)

# Get median value for each feature
values.log2.feature.median = apply(values.log2, 1, median)
# Median-normalize each feature
message("Now median-centering matrix ...")
values.median.centered = values.log2 - values.log2.feature.median

# Create output matrix
correlation.table = matrix(data=NA, ncol=4, nrow=nrow(values))
colnames(correlation.table) = c("Feature.ID", "Pearson correlation", "Pearson p-value", "Pearson FDR")
correlation.table[,1] = rownames(values)

# Get reference profile
if (reference.feature %in% rownames(values.median.centered))
{
	values.median.centered.reference = values.median.centered[reference.feature,]
	message("Found reference feature.")
} else {
	message("\nError: Cannot find reference feature (", reference.feature, ") in your data table.\n")
	q()
}

message("Now calculating correlations relative to reference feature ...")

for (i in 1:nrow(values.median.centered))
{
	# If all values are 0, leave the metric and statistic as NA
	if (max(as.numeric(values.median.centered[i,]), na.rm=T) != min(as.numeric(values.median.centered[i,]), na.rm=T))
	{
		# this.corr.test = cor.test(as.numeric(values.median.centered[i,]), as.numeric(values.median.centered.reference), alternative="two.sided", method = "pearson")
		this.corr.test = cor.test(as.numeric(values.median.centered[i,]), as.numeric(values.median.centered.reference), alternative="two.sided", method = "pearson", use="pairwise.complete.obs")
		this.corr = this.corr.test$estimate
		this.p.value = this.corr.test$p.value
		correlation.table[i,2:3] = c(signif(this.corr, significant.figures), signif(this.p.value, significant.figures))
		# In case this is useful info
		# numNAs.this.row = sum(is.na(values[i,]))
	}
}

# Now correct p-values for FDR
FDR.P.value = p.adjust(as.numeric(correlation.table[,3]), "fdr")
# Add FDR values to matrix
correlation.table[,4] = signif(FDR.P.value, significant.figures)

correlation.table.sorted = correlation.table[order(as.numeric(correlation.table[,2]), na.last = TRUE, decreasing=TRUE),]

write.table(correlation.table.sorted, file=output.file, sep="\t", quote=F, row.names=F)

message(paste("\nAll done!  For output, see", output.file, "\n"))
