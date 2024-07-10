#!/usr/bin/env Rscript

### 
###  Run DESeq2 from the command line
###  George Bell - Bioinformatics and Research Computing, Whitehead Institute
###
###  USAGE: DESeq2_normalize_only.R inputCounts OutputFile 
###  EX: DESeq2_normalize_only.R InputDEseq2.txt DESeq2_output.txt
###
###
###  12 August 2016 [GWB]:
###    Modify DESeq2 code for Rscript syntax 
###      and make PCA figures with similar code as DESeq_normalize_only.R
###  6 Dec 2017: In varianceStabilizingTransformation, change 'blind' to FALSE
###  15 Oct 2018: Add 'check.names=FALSE' to read.delim()
###
###  30 Jan 2019
###  Added optional sampleColor file allowing users to define PCA colors
###  in plots.  
###
###  7 July 2019 -- Change variance-stabilizing transformation from vsn to rlog; in DESeq(), change fitType from 'mean' to default.
###  26 August 2019 -- Add legend to each PCA plot (so need to change format of colorFile to include group and color).
###                    Add choice of log method (log2|rlog|vst) where log2 includes first adding 1 pseudocount
###                    Add integer check for input counts => convert if they're not already integers
###  11 October 2019 -- Allow duplicate row names (1st column) which will be made unique [GB]
###  16 October 2019 -- Modify geom_text_repel() options and names of output table files to include log-transform method. [GB]
###  21 October 2019 -- Check colors file; clarify normalized output counts [GB]
###  2 September 2020 -- Adjust handling of colors file so order doesn't need to match that of matrix file [GB]
###  30 April 2021 -- Get colors more dependably match groups [see added commands starting on line 266] [GB]
###  3 December 2021 -- Add gray circles around points [GB]

# For user interface choice
offset = 0
# For output table
significant.figures.after.decimal = 2
# For PCA figures
point.size = 8
# For legend on PCA figures
legend.point.size = 16

if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./DESeq2_normalize_only.R" }

if (length(commandArgs()) < (8 - offset))
{
	message("\nNormalize a matrix of raw counts with DESeq2 BUT don't adjust by size factors (with optional PCA plot group colors).")
	message(paste("USAGE:   ", this.script, "inputCounts OutputFile pcaLogMethod[log2|rlog|vst] [colorFile]"))
	message(paste("Example: ", this.script, "Input_DESeq2.txt DESeq2_output.norm_only.txt rlog sampleGroupColors.txt\n"))
	q()
}

# Load library without messages
message("\nLoading DESeq2 and other required R packages ....")
suppressMessages(library(DESeq2))
suppressMessages(library(genefilter))
suppressMessages(library(RColorBrewer))
suppressMessages(library(lattice))
suppressMessages(library(ggrepel))
suppressMessages(library(reshape2))
suppressMessages(library(vsn))

message(paste("Running DESeq2 version", packageDescription("DESeq2")$Version))

# First argument is input file
input.filename = commandArgs()[6 - offset]
# Second argument is output file
output.filename = commandArgs()[7 - offset]
# Third argument is pcaLogMethod; needs to be log2|rlog|vst
pcaLogMethod = commandArgs()[8 - offset]
if (pcaLogMethod != "log2" && pcaLogMethod != "rlog" && pcaLogMethod != "vst")
{
	message("ERROR: pcaLogMethod [third argument] needs to be log2, rlog, or vst.")
	quit()
}

# Fourth argument is colorGroup file
groupColor.filename = commandArgs()[9 - offset]

message(paste("Reading input file: ", input.filename))

# Read matrix of counts -- allow duplicate row names (which will be made unique)
# Permit gene matrices with duplicated gene symbols [GB, 7 Sep 2018]
counts = read.delim(input.filename, header=T, sep="\t", check.names=FALSE)

message("Done reading input matrix.  Now processing ....")

# Get total counts per gene (for sorting redundant rows)
counts.sum.by.feature = apply(counts[,2:ncol(counts)], 1, sum)
# Order matrix by decreasing counts.sum.by.gene, adding a numeric suffix to each symbol we've already seen.
feature.symbols.original = as.character(counts[,1])
# feature.symbols.unique = .Internal(make.names(feature.symbols.original, unique=TRUE))
feature.symbols.unique = make.unique(feature.symbols.original)
o = order(feature.symbols.original != feature.symbols.unique, counts.sum.by.feature, decreasing=TRUE)
feature.symbols.unique[o] = make.unique(feature.symbols.unique[o])
# Make row names from unique gene symbols.
rownames(counts) = feature.symbols.unique
# Drop the gene symbol (first column) since we have them as row names
counts = counts[,-1]

# Order sample names alphabetically
original.sample.order = colnames(counts)
alphabetical.sample.order = sort(original.sample.order)
counts = counts[,alphabetical.sample.order]

# Check for file containing sample color values
if (length(commandArgs()) > 8)
{
	message("Reading sample colors....")
	group.color = read.delim(groupColor.filename, row.names=1, check.names=FALSE)
	
	# Check that samples in colors file matches those in counts matrix [GB -- 2 Sep 2020]
	if (all(sort(colnames(counts)) != sort(as.vector(rownames(group.color)))))
	{
		message("\nERROR: First column of the colors file should include the same entries as the first row of the counts matrix.\n")
		message(colnames(counts))
		message(as.vector(group.color[,1]))
		quit()
	} else {
		message("OK: First column of the colors file includes the same entries as the first row of the counts matrix.\n")
	}
	# Order colors in same order as input data matrix (to get colors to match counts samples)
	group.color = group.color[colnames(counts),]

	# Ignore given names; set column names to what we want
	colnames(group.color) = c("group", "color")
	sample.colors = group.color[,2]
	
	print(group.color)
	
	colors.provided = 1
	group.color.unique = unique(group.color)
	# Check to see if we have one color per group (rather than per sample)
	if (length(group.color.unique[,1]) != length(unique(group.color.unique[,1])))
	{
		message("\nERROR: Colors file should have one color per group.\n")
		quit()
	}
	rownames(group.color.unique) = group.color.unique[,1]
} else {
	message("No sample colors file is being used.")
	sample.colors = "black"
	colors.provided = 0
}

# In case we want to bend the rules and not start with an integer matrix,
# see whether we have integers.  If not, convert them.
message("Flooring \'counts\' to 0 and converting to integers (if they aren't already)...")
counts.pos = apply(counts, c(1, 2), function(x) ifelse(x<0,0,x))
counts = round(counts.pos, 0)

# Make a DESeqDataSet
samples = factor(colnames(counts))
# Use a design stating that there are no replicates
# If there is replication, replace 'design' with the real design
message("Treating all samples as independent (so assuming no replication)....\n")
dds = DESeqDataSetFromMatrix(countData=counts, colData=DataFrame(samples), design=~1)

###  MAJOR COMMAND:
# Estimate size factors and dispersions and fit generalized linear model
# For designs without replication, fitType has no influence on size factors
# dds = DESeq(dds, fitType="local")

dds = DESeq(dds)

# This version only!!!
# Set size factors to 1 (so don't do estimateSizeFactor scaling)
message("\n!!! WARNING: NOT performing any scaling by library size !!!\n")
sizeFactors(dds) = rep(1, ncol(dds))

# Print size factors
message("\nSize factors:")
sizeFactors(dds)
message()

# Add normalized counts for this version (GB - 19 Mar 2013)
counts.normalized = round(t(t(counts(dds))/sizeFactors(dds)), significant.figures.after.decimal)

# Print output (including norm counts)
# output.table = cbind(rownames(dds), counts, counts.normalized)
# Print output only of normalized counts
colnames(counts.normalized) = colnames(counts)

# Put back into original order
counts.normalized = counts.normalized[,original.sample.order]

output.table = cbind(rownames(dds), counts.normalized)
colnames(output.table)[1] = "Feature.ID"
write.table(output.table, file=output.filename, sep="\t", quote=F, row.names=F)

ggplot.pca.noColors = function (x, y, fac, pca.x.name, pca.y.name, ntop=500)
{
	# Set up colors and shapes of points
	ggplot() +
	aes(x, y) +
	geom_point(shape=21, color = sample.colors, size = point.size, fill = sample.colors) + 
	geom_point(color = "gray", size = point.size, shape=1) + 
	labs(x = pca.x.name) +
	labs(y = pca.y.name) + 
	theme_classic(base_size = 16) +
	ggtitle(paste("Using top", ntop, "features")) + 
	geom_text_repel(aes(x, y, label = fac), point.padding=0.5, size=5) + 
	# coord_fixed(ratio = 1) + 					# Add this to make x- and y-axes have the same scale
	theme_bw()	# To add axis lines
	# coord_equal(ratio = 1)	# Make isometric axes (so each scale is represented by the same distance)
	# geom_text(aes(x, y), label=fac, hjust=0, nudge=nudge, check_overlap = TRUE)
}

ggplot.pca = function (xy, group.color, pca.x.name, pca.y.name, ntop=500)
{
	group.color.unique = unique(group.color)
	rownames(group.color.unique) = group.color.unique[,1]
	group.color.unique = group.color.unique[order(group.color.unique[,1]),]
	my.groups = group.color.unique[,1]
	my.group.colors = group.color.unique[my.groups,2]
	my.group.colors = as.vector(my.group.colors)
	my.groups = as.vector(my.groups)
	
	if (length(my.group.colors) != length(my.groups))
	{
		print(paste("ERROR: number of group colors:", paste(unique(my.group.colors), collapse=";"), "must match number of groups:", paste(unique(my.groups), collapse=";")))
	}
	
	groups = as.vector(group.color[,2])

	# Set up colors and shapes of points
	ggplot(xy, aes(x, y), group=as.vector(group.color[,1])) +
	geom_point(aes(color = as.vector(group.color[,2])), size = point.size, fill="gray") + 
	geom_point(color = "gray", size = point.size, shape=1) + 
	scale_color_manual(values = my.group.colors, breaks = my.group.colors, labels = my.groups, name="group") + 
	theme_bw() + 
	theme(legend.text=element_text(size=legend.point.size), legend.title = element_text(size=legend.point.size)) + 
	labs(x = pca.x.name) +
	labs(y = pca.y.name) + 
	labs(fill = "Groups") + 
	ggtitle(paste("Using top", ntop, "features")) + 
	geom_text_repel(aes(x, y, label = rownames(xy)), point.padding=0.5, size=5, max.overlaps=100)
}

plotPCA_BaRC_v3 = function (x, ntop = 500)
{
	# Need to convert SummarizedExperiment into a plain old matrix
	if (pcaLogMethod != "log2")
	{
		x = assay(x)
	}
    	rv = rowVars(x)
	rv = genefilter::rowVars(x)
	select = order(rv, decreasing = TRUE)[seq_len(ntop)]
	pca = prcomp(t(x[select, ]))
	
	#add this line to get the weights of the genes on each component
	write.table(pca$rotation, file=paste(gsub(".txt", "", input.filename), pcaLogMethod, "pca_genes.using.top", ntop, "txt", sep="."),sep="\t", quote=F)
	# Get the weights of the samples on each component
	write.table(pca$x, file=paste(gsub(".txt", "", input.filename), pcaLogMethod, "pca_samples.using.top", ntop, "txt", sep="."),sep="\t", quote=F)

	# Eigenvalues
	eig = (pca$sdev)^2
	variances = round(eig*100/sum(eig), 2)
	names(variances) = colnames(pca$x)
	
	xy.12 = data.frame(pca$x[,c("PC1","PC2")])
	xy.13 = data.frame(pca$x[,c("PC1","PC3")])
	xy.23 = data.frame(pca$x[,c("PC2","PC3")])
	colnames(xy.12) = c("x", "y")
	colnames(xy.13) = c("x", "y")
	colnames(xy.23) = c("x", "y")

	# Needed to get colors to correctly match groups?!?
	xy.12 = xy.12[colnames(counts),]
	xy.13 = xy.13[colnames(counts),]
	xy.23 = xy.23[colnames(counts),]
	
	par(las=1)
	if (colors.provided > 0)
	{
		print(ggplot.pca(xy.12, group.color, paste("PC1:   ", variances[1], "% variance", sep=""), paste("PC2:   ", variances[2], "% variance", sep=""), ntop))
		print(ggplot.pca(xy.13, group.color, paste("PC1:   ", variances[1], "% variance", sep=""), paste("PC3:   ", variances[3], "% variance", sep=""), ntop))
		print(ggplot.pca(xy.23, group.color, paste("PC2:   ", variances[2], "% variance", sep=""), paste("PC3:   ", variances[3], "% variance", sep=""), ntop))
	} else {
		print(ggplot.pca.noColors(pca$x[,"PC1"], pca$x[,"PC2"], rownames(pca$x), paste("PC1:   ", variances[1], "% variance", sep=""), paste("PC2:   ", variances[2], "% variance", sep=""), ntop))
		print(ggplot.pca.noColors(pca$x[,"PC1"], pca$x[,"PC3"], rownames(pca$x), paste("PC1:   ", variances[1], "% variance", sep=""), paste("PC3:   ", variances[3], "% variance", sep=""), ntop))
		print(ggplot.pca.noColors(pca$x[,"PC2"], pca$x[,"PC3"], rownames(pca$x), paste("PC2:   ", variances[2], "% variance", sep=""), paste("PC3:   ", variances[3], "% variance", sep=""), ntop))
	}
	message(paste("Variances (as percentage) for top", ntop, "features"))
	print(variances)
	# barplot(eig, ylab="eigenvalues", main=paste("Eigenvalues for top", ntop, "features"))
	barplot(variances, ylab="variance (as percentage)", main=paste("Variances (as percentage) for top", ntop, "features"))
}


### Do VST transformation (required to  for PCA plot to remove the dependence of the variance on the mean)
# From manual:
# blind=TRUE should be used for comparing samples in an manner unbiased by prior information on samples, for example to perform sample QA (quality assurance). 
# blind=FALSE should be used for transforming data for downstream analysis, where the full use of the design information should be made.
# dds.vst = varianceStabilizingTransformation(dds, blind=FALSE, fitType="mean")
# Change to rlog (July 17, 2019)
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

### Run Principal Components Analysis
pca.figure.filename = paste(gsub(".txt", "", input.filename), pcaLogMethod, "PCA_etc.pdf", sep=".")
pdf(pca.figure.filename, w=8.5, h=8.5, useDingbats=FALSE)
# Fill the page
# plotPCA_BaRC_DESeq2(dds.rlog, aspect="fill")
# Make both axes on the same scale
# plotPCA_BaRC_DESeq2(dds.rlog, aspect="iso")

# Add more to this list (18 Jan 2021)
ntops = c(100, 300, 500, 1000, 2500, 5000, 10000, 20000, 40000, 50000)
ntops = sort(c(ntops , nrow(dds.log)), decreasing=F)

for (i in 1:length(ntops))
{
	if (nrow(dds.log) >= ntops[i])
	{
		plotPCA_BaRC_v3(dds.log, ntop=ntops[i])
	}
}

# Draw plot of SD vs mean
message("\nDrawing meanSdPlot for quality control (using log-transformed matrix without 0-count features) ....")
if (pcaLogMethod == "log2")
{
	dds.log.matrix = dds.log
} else {
	dds.log.matrix = assay(dds.log)
}
gene.total.counts = apply(dds.log.matrix, 1, sum)
dds.log.matrix.no0 = dds.log.matrix[gene.total.counts > 0,]
meanSdPlot(dds.log.matrix.no0)

# check dispersion distribution
# If fitted curve is off, you can try with local dispersion in stead of default parametric dispersion.
# You can estimated the dispersion with residue as mentioned in the post by Michael Love:
# https://support.bioconductor.org/p/81094/
message(paste("Create dispersion distribution figure","\n","If fitted curve is off, you can try with local dispersion in stead of default parametric dispersion","\n"))
plotDispEsts(dds)

num.samples = length(samples)
# Increases file size too much; not worth it
# par(las=2, mai=c(2,1,1,1))
# boxplot(log(cbind(counts, counts.normalized) + 1, 2), ylab="log2 ( counts + 1 )", col="wheat")

# Close output file
foo = dev.off()

message("\nAll done!")
message(paste("  See", pca.figure.filename, "for PCA plots."))
message(paste("  See", output.filename, "for size-factor normalized counts, not log-transformed.\n"))


###################

# save.image()
