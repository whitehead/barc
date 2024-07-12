#!/usr/bin/env Rscript

###
###  Draw a MA plot for a one-comparison limma analysis
###
###  George Bell - BaRC
###  Version 1.0: 2 Nov 2018
###  Input file: Output file from limma's results() object
###    with gene symbols in the first column
###  Draw MA plot, coloring genes that meet threshold, 
###    and labeling up to num.genes.to.label (set to 100) genes.
###  Derived from draw_MA_plot_from_DESeq2_analysis.R
###

version = "Version 1.0 -- 22 December 2021"

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
# y.axis.min = -10
# y.axis.max = 10

######  Input file structure (column number for each important field)  ######

log2mean.colNum = 1
log2FC.colNum = 2
tStat.colNum = 3
pvalue.colNum = 4
fdr.colNum = 5

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
gene.symbols.filename = commandArgs()[11 - offset]

# Get path to input.filename (so we can save the output file in the same dir)
input.dir = dirname(input.filename)

if ( ! is.na(log2FC.threshold) ) {
	message(paste("\nUsing limma output file *", input.filename, "* for input", sep=""))
	message(paste("Group A name is *", sampleA.name, "*", sep=""))
	message(paste("Group B name is *", sampleB.name, "*", sep=""))
	message(paste("FDR threshold is *", fdr.threshold, "*", sep=""))
	message(paste("Log2 fold change threshold is *", log2FC.threshold, "*", sep=""))
} else {
	message("\nDraw a MA plot from limma output,")
	message("  coloring points that meet FDR and fold change thresholds,")
	message("  and labeling top 100 points based on p-value + t-statistic.")
	message("USAGE: ", this.script, " inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold")
	message(paste(version, "\n"))
	quit()
}

# Read limma output file, and make the first column row names
# limma.stats = read.delim(input.filename, row.names=1)

# Input gene matrix
limma.stats = read.delim(input.filename, header=T, sep="\t", check.names=F)

# Drop rows with a FDR of NA [Keep these - GB - 3/20/2019]
# limma.stats = limma.stats[! is.na(limma.stats[,6]),]

# Order genes by increasing p-values
limma.stats = limma.stats[order(limma.stats[,4], -1*abs(limma.stats[,3]), decreasing=FALSE), ]

# For plot labels
comparison = paste("log2 ( ", sampleA.name, " / ", sampleB.name, " )", sep="")
figure.title = paste("MA plot:", sampleA.name, "vs.", sampleB.name)

# Get IDs
gene.symbols = rownames(limma.stats)
absT = abs(limma.stats[,tStat.colNum])
p.value = limma.stats[,pvalue.colNum]
# Get the FDR and then do -log10
fdr = limma.stats[,fdr.colNum]
minus.log10.fdr = -log(fdr, 10)
# Get the log2 fold change
log2FC = limma.stats[,log2FC.colNum]
# log2 mean
log2.counts.mean = limma.stats[,log2mean.colNum]

# Filter by FDR and fold change (for coloring)
# Filter genes with no FDR (FDR=NA) here - GB - 3/20/2019
filtered.rows =      ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold))
filtered.rows =      (((!is.na(fdr)) & (fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold)))
filtered.rows.up =   ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold) & log2FC > 0)
filtered.rows.down = ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold) & log2FC < 0)

# Use p-value threshold of 100th gene in this list
limma.stats.filtered = limma.stats[filtered.rows,]
# message(paste("limma.stats.filtered has", nrow(limma.stats.filtered), "rows."))
# Select the top num.genes.to.label genes (or all of them, if limma.stats.filtered is smaller)
# Fixed 2 Feb 2019
pvalue.threshold.labels = limma.stats.filtered[min(num.genes.to.label,nrow(limma.stats.filtered)),pvalue.colNum]
absT.threshold.labels = abs(limma.stats.filtered[min(num.genes.to.label,nrow(limma.stats.filtered)),tStat.colNum])

# Testing
message(paste("Raw p-value threshold for labels (top", num.genes.to.label, "genes by p-value) = ", pvalue.threshold.labels))
message(paste("abs(t) threshold for labels (top", num.genes.to.label, "genes by p-value) = ", absT.threshold.labels))

filtered.rows.labels =      ((!is.na(fdr)) & (p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold) & (absT > absT.threshold.labels))
filtered.rows.up.labels =   ((!is.na(fdr)) & (p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold) & (absT > absT.threshold.labels) & log2FC > 0)
filtered.rows.down.labels = ((!is.na(fdr)) & (p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold) & (absT > absT.threshold.labels) & log2FC < 0)

num.points.to.label = length(log2FC[filtered.rows])
message(paste(num.points.to.label, "points meet FDR and FC thresholds."))
message(paste(length(log2FC[filtered.rows.labels]), "points will be colored (top ~", num.genes.to.label, "by p-value, after filtering by FDR and FC)."))
gene.symbols[filtered.rows.labels]

comparison = paste(sampleB.name, "vs",  sampleA.name)
log2FC.label = paste("log2 ( ", sampleA.name, " / ",  sampleB.name, " )", sep="")

# Create plot with ggplot2

pdf.filename = paste("MA_plot", sampleA.name, "vs",  sampleB.name, fdr.threshold, log2FC.threshold, "pdf", sep=".")
pdf(paste(input.dir, pdf.filename, sep="/"), w=11, h=8.5, useDingbats=FALSE)

if (num.points.to.label < num.genes.to.label)
{
	ggplot() +
	# xlim(0, 20) + 
	# ylim(y.axis.min, y.axis.max) + 
	geom_point(aes(log2.counts.mean, log2FC), color="darkgray", size=1) + 
	geom_point(aes(log2.counts.mean[filtered.rows.up.labels], log2FC[filtered.rows.up.labels]), color="red", size=3) + 
	geom_point(aes(log2.counts.mean[filtered.rows.down.labels], log2FC[filtered.rows.down.labels]), color="blue", size=3) + 
	geom_point(aes(log2.counts.mean[filtered.rows.up], log2FC[filtered.rows.up]), color=up.color, size=1) + 
	geom_point(aes(log2.counts.mean[filtered.rows.down], log2FC[filtered.rows.down]), color=down.color, size=1) + 
	geom_hline(yintercept = 0, color="grey80") + 
	labs(x = "log2 mean mRNA level") +
	labs(y = log2FC.label) + 
	theme_classic(base_size = 16) +
	ggtitle(figure.title) +
	theme(plot.title = element_text(hjust = 0.5)) +
	geom_text_repel(aes(log2.counts.mean[filtered.rows.labels], log2FC[filtered.rows.labels]), label = gene.symbols[filtered.rows.labels], size=4, max.overlaps=100)
} else {
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
	geom_text_repel(aes(log2.counts.mean[filtered.rows.labels], log2FC[filtered.rows.labels]), label = gene.symbols[filtered.rows.labels], size=4, max.overlaps=100)
}

foo = dev.off()

message(paste("Created plot!  See", pdf.filename, "\n"))

# Save details of your session
sink("sessionInfo.txt")
date()
print(sessionInfo(), locale=FALSE)
sink()
