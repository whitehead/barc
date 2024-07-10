#!/usr/bin/env Rscript

###
###  Draw a MA plot for a DESeq2 analysis
###
###  George Bell - BaRC
###  Version 1.0: 2 Nov 2018
###  Input file: Output file from DESeq2's results() object
###    with gene symbols in the first column
###  Draw MA plot, coloring genes that meet threshold, 
###    and labeling up to num.genes.to.label (set to 100) genes.
###  Plot points with NA for FDR but just don't label them (instead of not plotting them at all) [20 March 2019, GB] 
###

# version = "Version 2.0 -- 2 November 2018"
# version = "Version 2.1 -- 8 February 2019"
# version = "Version 2.2 -- 20 March 2019"
# version = "Version 2.3 -- 2 May 2022"	# Keep gene IDs as is instead of "correcting" them
version = "Version 2.4 -- 27 January 2023"	# Allow for adjustment of y.axis.min and y.axis.max, and max.overlaps


suppressMessages(library(ggplot2))
library(ggrepel)

# Maximum number of genes to label with gene names
num.genes.to.label = 100
# Color for up genes
up.color = "#FFAAAA"
# Color for down genes
down.color = "#AAAAFF"
# Offset to add before calculating log2(baseMean)
counts.mean.offset = 1
# Y-axis limits
# y.axis.max = 6
# y.axis.min = -6
### Should we label selected points or just color them?
# Label them
max.overlaps = 1000
# Do not label them
# max.overlaps = 0


offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./draw_MA_plot_from_DESeq_analysis.R" }

input.filename = commandArgs()[6 - offset]
sampleA.name = commandArgs()[7 - offset]
sampleB.name = commandArgs()[8 - offset]
fdr.threshold = as.numeric(commandArgs()[9 - offset])
log2FC.threshold = as.numeric(commandArgs()[10 - offset])
genes.to.label.filename = commandArgs()[11 - offset]

if ( ! is.na(genes.to.label.filename) ) {
	message(paste("\nUsing DESeq2 output file *", input.filename, "* for input", sep=""))
	message(paste("Group A name is *", sampleA.name, "*", sep=""))
	message(paste("Group B name is *", sampleB.name, "*", sep=""))
	message(paste("FDR threshold is *", fdr.threshold, "*", sep=""))
	message(paste("Log2 fold change threshold is *", log2FC.threshold, "*", sep=""))
} else {
	message("\nDraw a MA plot from DESeq2 output,")
	message("  coloring selcted points that meet FDR and fold change thresholds,")
	message("  and labeling user-selected points.")
	message("USAGE: ", this.script, " inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold selectedGenesFile")
	message(paste(version, "\n"))
	quit()
}

# Read DESeq output file, and make the first column row names
# DESeq2.stats = read.delim(input.filename, row.names=1)
###  Modified to permit gene matrices with duplicated gene symbols [GB, 7 Sep 2018]
# Input gene matrix (no row names)
DESeq2.stats = read.delim(input.filename, header=T, sep="\t", check.names=F)
# Order matrix by decreasing baseMean, adding a numeric suffix to each symbol we've already seen.
gene.symbols.original = as.character(DESeq2.stats[,1])
# gene.symbols.unique = .Internal(make.names(gene.symbols.original, unique=TRUE))
gene.symbols.unique = make.unique(gene.symbols.original)
o = order(gene.symbols.original != gene.symbols.unique, DESeq2.stats$baseMean, decreasing=TRUE)
gene.symbols.unique[o] = make.unique(gene.symbols.unique[o])
# Make row names from unique gene symbols.
rownames(DESeq2.stats) = gene.symbols.unique
# Drop the gene symbol (first column) since we have them as row names
DESeq2.stats = DESeq2.stats[,-1]

# Drop rows with a FDR of NA [Keep these - GB - 3/20/2019]
# DESeq2.stats = DESeq2.stats[! is.na(DESeq2.stats[,6]),]

# Order genes by increasing p-values
DESeq2.stats = DESeq2.stats[order(DESeq2.stats[,5], decreasing=FALSE), ]

# For plot labels
comparison = paste("log2 ( ", sampleA.name, " / ", sampleB.name, " )", sep="")
figure.title = paste("MA plot:", sampleA.name, "vs.", sampleB.name)

# If we limit the axis, use the limits as the floor and/or the ceiling
if (exists("y.axis.min"))
{
	DESeq2.stats[,2] = ifelse(DESeq2.stats[,2] < y.axis.min, y.axis.min, DESeq2.stats[,2])
	message("\nSetting floor of log2FC to ", y.axis.min)
}
if (exists("y.axis.max"))
{
	DESeq2.stats[,2] = ifelse(DESeq2.stats[,2] > y.axis.max, y.axis.max, DESeq2.stats[,2])
	message("Setting ceiling of log2FC to ", y.axis.max)
}

# Get IDs
gene.symbols = rownames(DESeq2.stats)
p.value = DESeq2.stats[,4]
# Get the FDR and then do -log10
fdr = DESeq2.stats[,5]
minus.log10.fdr = -log(fdr, 10)
# Get the log2 fold change
log2FC = DESeq2.stats[,2]
# log2 mean
log2.counts.mean = log(DESeq2.stats[,1] + counts.mean.offset, 2)

# Filter by FDR and fold change (for coloring)
# Filter genes with no FDR (FDR=NA) here - GB - 3/20/2019
filtered.rows =      ((fdr < fdr.threshold) & (abs(log2FC) >= log2FC.threshold))
filtered.rows =      (((!is.na(fdr)) & (fdr <= fdr.threshold) & (abs(log2FC) >= log2FC.threshold)))
filtered.rows.up =   ((fdr <= fdr.threshold) & (abs(log2FC) >= log2FC.threshold) & log2FC > 0)
filtered.rows.down = ((fdr <= fdr.threshold) & (abs(log2FC) >= log2FC.threshold) & log2FC < 0)

# Use p-value threshold of 100th gene in this list
DESeq2.stats.filtered = DESeq2.stats[filtered.rows,]

# Look for genes to label in expression matrix
genes.to.label = unlist(read.delim(genes.to.label.filename, header=T, sep="\t", check.names=F))
num.original.genes.to.label = length(genes.to.label)
genes.to.label = intersect(genes.to.label, rownames(DESeq2.stats))
genes.to.label.in.matrix = length(genes.to.label)
message(paste("\nOf", num.original.genes.to.label, "genes in your list for labeling, the expression matrix contains", genes.to.label.in.matrix, "of these."))

filtered.rows.labels = (rownames(DESeq2.stats) %in% genes.to.label) & (DESeq2.stats$baseMean > 0)
num.filtered.rows.labels = length(log2.counts.mean[filtered.rows.labels])
message(paste("Of these,", num.filtered.rows.labels, "are expressed (with a level greater than 0)."))

filtered.rows.up.labels =   ((filtered.rows.labels) & log2FC >= 0)
filtered.rows.down.labels = ((filtered.rows.labels) & log2FC < 0)
num.filtered.rows.up.labels = length(sort(log2.counts.mean[filtered.rows.up.labels]))
num.filtered.rows.down.labels = length(sort(log2.counts.mean[filtered.rows.down.labels]))

message(paste("Of the expressed genes in your list and in the matrix,\n",
  " ", num.filtered.rows.up.labels, "have a logFC >= 0, and\n", 
  " ", num.filtered.rows.down.labels, "have a logFC < 0."))

num.points.to.label = length(log2FC[filtered.rows])
message(paste(num.points.to.label, "points meet FDR and FC thresholds."))
message(paste(length(log2FC[filtered.rows.labels]), "points will be labeled (user-selected gene symbols)."))

comparison = paste(sampleB.name, "vs",  sampleA.name)
log2FC.label = paste("log2 ( ", sampleA.name, " / ",  sampleB.name, " )", sep="")

# Create plot with ggplot2

pdf.filename = paste("MA_plot", sampleA.name, "vs",  sampleB.name, fdr.threshold, log2FC.threshold, gsub(".txt", "", genes.to.label.filename), "pdf", sep=".")
pdf(pdf.filename, w=11, h=8.5, useDingbats=FALSE)

ggplot() +
# xlim(0, 20) + 
# ylim(y.axis.min, y.axis.max) + 
geom_point(aes(log2.counts.mean, log2FC), color="darkgray", size=1) + 
geom_point(aes(log2.counts.mean[filtered.rows.up], log2FC[filtered.rows.up]), color=up.color, size=1) + 
geom_point(aes(log2.counts.mean[filtered.rows.down], log2FC[filtered.rows.down]), color=down.color, size=1) + 
geom_point(aes(log2.counts.mean[filtered.rows.up.labels], log2FC[filtered.rows.up.labels]), color="red", size=3) + 
geom_point(aes(log2.counts.mean[filtered.rows.down.labels], log2FC[filtered.rows.down.labels]), color="blue", size=3) + 
geom_hline(yintercept = 0, color="grey80") + 
labs(x = "log2 mean mRNA level") +
labs(y = log2FC.label) + 
theme_classic(base_size = 16) +
ggtitle(figure.title) +
theme(plot.title = element_text(hjust = 0.5)) +
geom_text_repel(aes(log2.counts.mean[filtered.rows.labels], log2FC[filtered.rows.labels]), label = gene.symbols[filtered.rows.labels], size=4, max.overlaps=max.overlaps)

foo = dev.off()

message(paste("\nCreated plot!  See", pdf.filename, "\n"))

# Save details of your session
sink("sessionInfo.txt")
date()
print(sessionInfo(), locale=FALSE)
sink()
