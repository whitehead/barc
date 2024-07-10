#!/usr/bin/env Rscript

###
###  George Bell
###  Bioinformatics and Research Computing, Whitehead Institute
###
###  Automate all vs.all comparisons using DESeq2.
###  Inputs: raw counts matrix and design file
###  Outputs: for each comparison, print a results table and a MA plot
###
###  Version 1.0: 8 October 2109
###  Version 1.1: 5 February 2020 -- Add lfcShrink with 'type="ashr"' (since we set the contrast)
###  Version 1.2: 28 April 2021   -- Floor input matrix to 0 and convert to integers (to allow bending the rules a bit)
###  Version 1.3: 18 November 2021   -- Print original gene symbols in output files (rather than R-corrected versions)
###  Version 1.4: 31 January 2022   -- Add volcano plots
###  Version 1.5: 20 July 2022   -- Change order of groups for creation of contrast matrix to correct occassional bug
###  Version 1.6: 25 January 2023  -- In rnk file (for GSEA) exclude unexpressed genes (instead of listing them with log2FC=0)
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
if (is.na(this.script)) { this.script = "./Run_all_vs_all_DESeq2_differential_expression.R" }

matrix.of.counts = commandArgs()[6 - offset]
design.file = commandArgs()[7 - offset]
output.dir = commandArgs()[8 - offset]

# For MA plots
fdr.threshold = 1e-05
log2FC.threshold = 1
# For normalized counts
significant.figures.after.decimal = 4
# How do we log-transform?
# pcaLogMethod = "rlog"
pcaLogMethod = "vst"
# pcaLogMethod = "log2"
# Do lfcShrink? If not, comment this out.  [5 February 2020]
do.lfcShrink = "yes"


# Sample command: ./Run_DESeq2_differential_expression.R Three_groups.counts.txt Three_groups.design.txt Three_groups.counts.DESeq2_comparisons_v2
# setwd("/nfs/BaRC_Public/BaRC_code/R/DESeq2/complex_design_analysis")
# matrix.of.counts = "Three_groups.counts.txt"
# design.file = "Three_groups.design.txt"
# output.dir = "Three_groups.counts.DESeq2_comparisons_v3"

if ( ! is.na(output.dir) ) {
	message(paste("Using matrix file *", matrix.of.counts, "* for input.", sep=""))
	message(paste("Design file: *", design.file, "*", sep=""))
	message(paste("Output directory: *", output.dir, "*", sep=""))
} else {
	message("Do an all vs. all DESeq2 analysis,")
	message("  starting with a matrix of raw counts (NOT pre-normalized).")
	message("\nUSAGE: ", this.script, " inputMatrix designFile outputDirectory\n")
	quit()
}

# Load library (quietly)
message("\nLoading DESeq2 and other required R packages ....")
library("limma")
suppressMessages(library(DESeq2))
library("vsn")

# Read matrix of counts
# Permit gene matrices with duplicated gene symbols [GB, 7 Sep 2018]	Added "check.names=FALSE" on 28 Jan 2021
message("Reading counts matrix: ", matrix.of.counts)
counts = read.delim(matrix.of.counts, header=T, sep="\t", check.names=FALSE)
# Get total counts per gene (for sorting redundant rows)
counts.sum.by.gene = apply(counts[,2:ncol(counts)], 1, sum)
# Order matrix by decreasing counts.sum.by.gene, adding a numeric suffix to each symbol we've already seen.
gene.symbols.original = as.character(counts[,1])
# gene.symbols.unique = .Internal(make.names(gene.symbols.original, unique=TRUE))
gene.symbols.unique = make.unique(gene.symbols.original)
o = order(gene.symbols.original != gene.symbols.unique, counts.sum.by.gene, decreasing=TRUE)
gene.symbols.unique[o] = make.unique(gene.symbols.unique[o])
# Make row names from unique gene symbols.
rownames(counts) = gene.symbols.unique
# Drop the gene symbol (first column) since we have them as row names
counts = counts[,-1]

# Read design as a file
message("Reading design file: ", design.file)
design = read.delim(design.file, row.names=1, check.names=FALSE)
design.samples = rownames(design)
# design.samples
# colnames(values)

# Check that columns in the matrix file match rows in the design file
if (! identical(colnames(counts), design.samples))
{
	message("\nUh oh -- the column names of your matrix (after \"cleaning up\" by R)")
	print(colnames(counts))
	message("don't seem to match the sample names in the first column of your design file")
	print(design.samples)
	message()
	quit()
}

# In case we want to bend the rules and not start with an integer matrix,
# see whether we have integers.  If not, convert them.
message("Flooring \'counts\' to 0 and converting to integers (if they aren't already)...")
counts.pos = apply(counts, c(1, 2), function(x) ifelse(x<0,0,x))
counts = round(counts.pos, 0)

# Use just the first column of the design file (second column really, since first column should be sample IDs/names)
design = design[,1]
design.matrix = model.matrix(~0+design)
colnames(design.matrix) = gsub("design", "", colnames(design.matrix))
rownames(design.matrix) = design.samples
message("\nExperimental design:")
design.matrix

# Get all possible comparisons
# groups = unique(as.vector(design))
groups = sort(unique(as.vector(design)))  # Explicitly sort [7/20/2022]
num.groups = length(groups)
comparisons = vector(length = num.groups * num.groups - num.groups)
comparisons.text = vector(length = num.groups * num.groups - num.groups)
comparison.num = 1
for (i in 1:length(groups))
{
	for (j in 1:length(groups))
	{
		if (i != j)
		{
			comparisons[comparison.num] = paste(groups[i], "-", groups[j])
			comparisons.text[comparison.num] = paste(groups[i], "vs", groups[j], sep=".")
			comparison.num = comparison.num + 1
		}
	}
}

# Make a DESeqDataSet
dds = DESeqDataSetFromMatrix(countData = counts, colData = DataFrame(condition=design), design = design.matrix)

###  MAJOR COMMAND:
# Estimate size factors and dispersions and fit generalized linear model
message("Running DESeq() with 'betaPrior=FALSE' (no fold-change shrinkage) ...")
# dds = DESeq(dds, betaPrior=FALSE)
dds = DESeq(dds, betaPrior=FALSE)

# Print size factors
sizeFactors(dds)

# Create directory (if it doesn't already exist) for output files
if (! file.exists(output.dir))
{
	message(paste("Creating directory ", output.dir))
	dir.create(output.dir)
}

# Add normalized counts for this version (GB - 19 Mar 2013)
counts.normalized = round(t(t(counts(dds))/sizeFactors(dds)), significant.figures.after.decimal)
# Print output only of normalized counts
colnames(counts.normalized) = paste(colnames(counts), "norm", sep=".")
# output.table = cbind(rownames(dds), counts.normalized)
# Replace gene symbols with original ones in input file
output.table = cbind(gene.symbols.original, counts.normalized)
colnames(output.table)[1] = "Feature.ID"
write.table(output.table, file=paste(output.dir, paste(gsub(".txt", "", matrix.of.counts), "DESeq2_normalized_counts.txt", sep="."), sep="/"), sep="\t", quote=F, row.names=F)


###############  Create PCA

if (pcaLogMethod == "log2")
{
	message("Log-transforming normalized expression matrix by adding 1 pseudocount and then calculating log2....")
	dds.log = log2(counts.normalized + 1)
}
if (pcaLogMethod == "rlog")
{
	message("Log-transforming expression matrix with rlog....")
	dds.log = rlog(dds, blind=FALSE)
}
if (pcaLogMethod == "vst")
{
	message("Log-transforming expression matrix with vst....")
	dds.log = varianceStabilizingTransformation(dds, blind=FALSE, fitType="mean")
}

DESeq2_plotPCA_BaRC = 
function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
    library(ggplot2)
    if (pcaLogMethod == "log2")
    {
    	colData = DataFrame(condition=design)
    	object = SummarizedExperiment(assays=SimpleList(counts=object), colData=colData)
    }
    x = assay(object)
    rv = rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
    	stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
    	factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
    	colData(object)[[intgroup]]
    }

    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
    point.size = 8
    # plot.title = paste(gsub(".txt", "", matrix.of.counts), " (", ntop, " features)", sep="")
    plot.title = paste("Using ", ntop, " features", sep="")
    print(ggplot(data = d, aes_string(x = "PC1", y = "PC2", fill = "group")) + 
        geom_point(shape=21, size = point.size, color="gray") + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed() + theme_bw() + ggtitle(plot.title))
    # Add variation of this plot that fills the page
    print(ggplot(data = d, aes_string(x = "PC1", y = "PC2", fill = "group")) + 
        geom_point(shape=21, size = point.size, color="gray") + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + theme_bw() + ggtitle(plot.title))
}

pdf(file=paste(output.dir, paste(gsub(".txt", "", matrix.of.counts), "PCA_etc.pdf", sep="."), sep="/"), useDingbats=FALSE)
# PCA plot

# ntops = c(100, 300, 500, 1000, 5000)
# Add more to this list (18 Jan 2021)
ntops = c(100, 300, 500, 1000, 2500, 5000, 10000, 20000, 40000, 50000)
ntops = sort(c(ntops , nrow(dds.log)), decreasing=F)

for (i in 1:length(ntops))
{
	DESeq2_plotPCA_BaRC(dds.log, ntop=ntops[i], intgroup=c("condition"))
}

if (pcaLogMethod == "log2")
{
	meanSdPlot(dds.log)
	par(las=2, mai=c(2,1,1,1))
	boxplot(dds.log, main=paste(gsub(".txt", "", matrix.of.counts), " (rlog)", sep=""))
} else {
	meanSdPlot(assay(dds.log))
	par(las=2, mai=c(2,1,1,1))
	boxplot(assay(dds.log), main=paste(gsub(".txt", "", matrix.of.counts), " (rlog)", sep=""))
}

# sd vs ranked mean
foo = dev.off()

contrast.matrix = matrix(data=NA, nrow=num.groups, ncol=length(comparisons))
rownames(contrast.matrix) = resultsNames(dds)
colnames(contrast.matrix) = comparisons
message("Contrast matrix (all possible combinations):")
for (i in 1:length(comparisons))
{
	# This code chunk from Saroj Mohapatra (https://support.bioconductor.org/p/27900/)
	prestr="makeContrasts("
	# poststr=",levels=design)"	# 18 Nov 2021 => change to levels=unique(design)
	# poststr=",levels=unique(design))"
	poststr=",levels=sort(unique(design)))"  # Explicitly sort
	commandstr=paste(prestr,comparisons[i],poststr,sep="")
	contrast.matrix[,i] = eval(parse(text=commandstr))
}
contrast.matrix

# Go to the output directory so we can create MA plots there
setwd(output.dir)
# Make all of the desired comparisons
for (i in 1:length(comparisons))
{
	# Do stats
	#NOTE: by default, results function has independentFiltering=TRUE => in general, genes with low counts (based on mean normalized count) are removed (ie. adj p-val is NA for these)
	#use results(dds, independentFiltering=FALSE) to get adj p-val for all genes, however, you may want to do extra filtering (eg. remove genes with low counts).
	#see: https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt
	# res = results(dds)
	res = results(dds, contrast=contrast.matrix[,i], independentFiltering=FALSE)
	
	# Do we need lfcShrink()?
	# The default shrinkage method fails when supplying a contrast matrix:
	# res = lfcShrink(dds, contrast=contrast.matrix[,i])
	# We need to use 'type="ashr"' if we want to set the contrast [5 February 2020]
	if (do.lfcShrink == "yes")
	{
		res = lfcShrink(dds, contrast=contrast.matrix[,i], type="ashr")
	}
	
	res.file = paste(paste("DESeq2_out", comparisons.text[i], "txt", sep="."))
	feature.ID = rownames(res)
	# write.table(cbind(feature.ID, res), file=res.file, sep="\t", quote=F, row.names=F)
	# Replace with correct (original) symbols to match those in input file
	write.table(cbind(gene.symbols.original, res), file=res.file, sep="\t", quote=F, row.names=F)

	# Create rnk file for GSEA
	# Exclude all genes with expression (baseMean) of 0 [25 January 2023 -- GB]
	rnk.file = paste(paste("For_GSEA", comparisons.text[i], "rnk", sep="."))
	rnk.table = cbind(feature.ID, round(res[,2], 5))[res$baseMean>0, ]
	write.table(rnk.table, file=rnk.file, sep="\t", quote=F, row.names=F, col.names=F)

	samples = unlist(strsplit(comparisons[i], " - "))
	MA.plot.command =      paste("/nfs/BaRC_Public/BaRC_code/R/DESeq2/draw_MA_plot_from_DESeq2_analysis.R", res.file, samples[1], samples[2], fdr.threshold, log2FC.threshold)
	message(MA.plot.command)
	system(MA.plot.command)

	volcano.plot.command = paste("/nfs/BaRC_Public/BaRC_code/R/DESeq2/draw_volcano_plot_from_DESeq2_analysis.R", res.file, samples[1], samples[2], fdr.threshold, log2FC.threshold)
	message(volcano.plot.command)
	system(volcano.plot.command)

}

############

# Save details of your session
sink(paste("sessionInfo", Sys.Date(), "txt", sep="."))
date()
print(sessionInfo(), locale=FALSE)
sink()

###################

message("\nAll done (!) -- see files in directory ", output.dir, "\n")
