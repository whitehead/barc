#!/usr/bin/env Rscript

###
###  Correct batch effects with sva's ComBat() and limma's removeBatchEffect()
###
###  George Bell - BaRC
###  Version 1.0: 12 July 2019 
###  Version 1.1: 30 July 2019 -- Check that column names of expresssion matrix match row names of sample info file.  If not, exit.
###                               Remove all rows with all 0 values and then rows with at least half 0s.
###  Version 1.2: 24 November 2020 -- Allow choice to include log-transformation or not
###
###  Input files: 
###		Matrix of log2-transformed expression values
###		Sample information (including batch info and (optionally) covariate info to use for batch effect removal
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
if (is.na(this.script)) { this.script = "./remove_batch_effects.R" }

expression.file = commandArgs()[6 - offset]
sample.info.file = commandArgs()[7 - offset]
output.prefix = commandArgs()[8 - offset]

if (is.na(output.prefix))
{
	message("\nCorrect batch effects with sva's ComBat() and limma's removeBatchEffect().")
	message("USAGE: ", this.script, " log2Matrix sampleInfo OutfilePrefix [logTransform(yes|no)] [pseudocounts]")
	message("Ex:    ", this.script, " Expression.log2_values.txt Sample_information.txt MyBatchAdjustment yes 1\n")
	message("Notes:")
	message("  The rows of sampleInfo should match the columns of log2Matrix, and")
	message("  sampleInfo needs a column named \'batch\' and an optional column named \'covariate\'.\n")
	quit()
}

# Optional arguments for log-transformation
log.transform = commandArgs()[9 - offset]
optLT<-c("yes","no")
if (is.na(log.transform)){
   log.transform<-"no"
} else if (!(log.transform %in% optLT)){
   log.transform<-"no"
}
pseudocounts = as.numeric(commandArgs()[10 - offset])

message("Loading sva, limma, and other required modules ...")
###  See https://www.bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
suppressMessages(library(sva))
###  See https://bioconductor.org/packages/release/bioc/html/limma.html
suppressMessages(library(limma))
###  For PCA
suppressMessages(library(made4))

# Read expression matrix (log2 values)
# First column has gene IDs
message("Loading input matrix ...")
expression = read.delim(expression.file, check.names=F, row.names=1)

if (log.transform == "yes")
{
	if (is.na(pseudocounts))
	{
		message("\nIf log-transformation is desired, please provide pseudocounts to add before transforming.")
		message("USAGE: ", this.script, " log2Matrix sampleInfo OutfilePrefix [logTransform(yes|no)] [pseudocounts]")
		message("Ex:    ", this.script, " Expression.log2_values.txt Sample_information.txt MyBatchAdjustment yes 1\n")
		quit()
	} else {
		message("Adding ", pseudocounts, " pseudocounts and log2-transforming ...")
		expression = log(expression + pseudocounts, 2)
	}
} else {
	message("\nSkipping log-transformation, so assuming that input matrix is already log2-transformed.\n")
}

# Removing all rows with all 0s [to prevent "Error in while (change > conv) { : missing value where TRUE/FALSE needed" error]
message(paste("Removing all rows with all 0s, starting with", nrow(expression), "rows ..."))
row.sums = apply(expression, 1, sum)
expression = expression[row.sums > 5, ]
message(paste("Continuing with", nrow(expression), "non-0 rows ..."))

# Removing all rows with at least half 0s [to prevent "Error in while (change > conv) { : missing value where TRUE/FALSE needed" error]
message(paste("Removing all rows with all at least half 0s, starting with", nrow(expression), "rows ..."))
num.0s.in.row = apply(expression == 0, 1, sum)
expression = expression[num.0s.in.row < (ncol(expression)/2), ]
message(paste("Continuing with", nrow(expression), "rows ..."))

# Read sample data
# Rownames (Column 1 in this example) contain sample IDs that match column IDs in expression file
# A column is named 'batch' (Column 2 in this example)
sample.info = read.delim(sample.info.file, check.names=F)
# Make rownames from the first column
rownames(sample.info) = sample.info[,1]
# Change (if needed) the name of the second column to "batch"
colnames(sample.info)[2] = "batch"

# Make sure we have the same order as the expression matrix
sample.info = sample.info[colnames(expression),]

if (length(setdiff(colnames(expression), rownames(sample.info))) > 0)
{
	message("ERROR!: Column names of expresssion matrix don't match row names of sample info file.")
	quit()
}

batch = sample.info$batch
# Include covariates if provided in sample.info file
# Column needs to be called "covariate"
if (is.factor(sample.info$covariate))
{
	message("Found covariates in sample info file.")
	modcombat = model.matrix(~as.factor(covariate), data=sample.info)
} else
{
	message("No covariates found in sample info file.")
	modcombat = model.matrix(~1, data=sample.info)
}

message("\nUsing sva's ComBat to try to correct for batch ...")
expression.after.ComBat = ComBat(dat=as.matrix(expression), batch=batch, mod=modcombat, par.prior=TRUE)
message("\nUsing limma's removeBatchEffect to try to correct for batch ...")
expression.after.limma = removeBatchEffect(expression, batch=batch, covariates=modcombat)

# Output expression matrices after batch correction
ComBat.corrected.outfile = paste(output.prefix, "ComBat-corrected.log2.txt", sep=".")
limma.corrected.outfile = paste(output.prefix, "limma-corrected.log2.txt", sep=".")
unlogged.ComBat.corrected.outfile = paste(output.prefix, "ComBat-corrected.txt", sep=".")
unlogged.limma.corrected.outfile = paste(output.prefix, "limma-corrected.txt", sep=".")

write.table(round(expression.after.ComBat, 4), ComBat.corrected.outfile, sep="\t", quote=F)
write.table(round(expression.after.limma, 4), limma.corrected.outfile, sep="\t", quote=F)

write.table(round(2**expression.after.ComBat, 4), unlogged.ComBat.corrected.outfile, sep="\t", quote=F)
write.table(round(2**expression.after.limma, 4), unlogged.limma.corrected.outfile, sep="\t", quote=F)



message("\nAll done!\n")

# This works but the output matrix has a completely different mean and SD as the input matrix (as with zero-centering)
# library(pamr)
# expression.list = list(x=as.matrix(expression),y=factor(colnames(expression)),batchlabels=factor(as.numeric(batch)))
# expression.after.pamrBatchAdjust = pamr.batchadjust(expression.list)
