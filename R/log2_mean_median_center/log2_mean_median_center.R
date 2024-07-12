#!/usr/bin/env Rscript

# Mean or median center a matrix of values:
# Input a matrix of values
# Add pseudocounts to all values
# Log2-transform
# subtract mean or median for each row
#
# For input file: assumes first row is header and first column has row IDs.
# George W Bell, BaRC

version = "Version 1.0 -- 29 July 2016"
# How many significant digits?
num.rounded.digits = 4
# How should we print NAs?
na = ""

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./mean_median_center.R" }

matrix.filename = commandArgs()[6 - offset]
summary.method = commandArgs()[7 - offset]
pseudocounts = as.numeric(commandArgs()[8 - offset])
# Should we add an extra row ID column to view output in Java TreeView?
jtv = commandArgs()[9 - offset]
if (! exists(jtv)) { jtv = "no" }

if ( ! is.na(pseudocounts) ) {
	message(paste("\nGetting raw values from *", matrix.filename, "*", sep=""))
	message(paste("Summarizing by *", summary.method, "*", sep=""))
	message(paste("Adding *", pseudocounts, "* pseudocounts before log-transformation", sep=""))
	if (jtv == "jtv") {
		message("Adding another column of row IDs to make output ready for Java TreeView viewing")
	}
} else {
	message("\nMean or median transform matrix of values (with row names in first column)")
	message("USAGE: ", this.script, " inputFile [mean|median] pseudocounts [jtv] > out.txt\n")
	# message(paste(version, "\n"))
	quit()
}

# Input a matrix of values
# Add pseudocounts to all values
# Log2-transform and subtract mean or median for each row
center = function(my.Matrix, summary.method, pseudocounts)
{
		values.plus = my.Matrix + pseudocounts
		log2.values.plus = log2(values.plus)
		log2.values.featureSummary = apply(log2.values.plus, 1, summary.method, na.rm=TRUE)
		values.centered = log2.values.plus - log2.values.featureSummary
		signif(values.centered, num.rounded.digits)
}

# Read values, assuming first column has row IDs and first row is header
values = read.delim(matrix.filename, row.names=1)
# Center values
if (summary.method == "mean" || summary.method == "median") {
	centered.values = center(values, summary.method, pseudocounts)
	message("Log-transforming matrix, summarizing rows, and centering ...")
} else {
	message("\nSummary method must be \'mean\' or \'median\'\n")
}

if (jtv == "jtv") {
	out.matrix = cbind(rownames(centered.values), rownames(centered.values), centered.values)
	colnames(out.matrix) = c("Position", "Position", colnames(values))	
} else {
	out.matrix = cbind(rownames(centered.values), centered.values)
	colnames(out.matrix) = c("Position", colnames(values))

}
write.table(out.matrix, sep="\t", quote=F, row.names=F, na=na)

message("\nAll done!\n")
