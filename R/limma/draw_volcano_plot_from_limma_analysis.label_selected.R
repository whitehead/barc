#!/usr/bin/env Rscript

###
###  Draw a volcano plot for a limma analysis
###
###  George Bell - BaRC
###  Version 1.0: 2 Nov 2018
###  Input file: Output file from limma's results() object
###    with gene symbols in the first column
###  Draw MA plot, coloring genes that meet threshold, 
###    and labeling up to num.genes.to.label (set to 100) genes.
###  Plot points with NA for FDR but just don't label them (instead of not plotting them at all) [20 March 2019, GB] 
###

# version = "Version 2.0 -- 2 November 2018"
# version = "Version 2.1 -- 8 February 2019"
# version = "Version 2.2 -- 20 March 2019"
version = "Version 1.0 -- 25 August 2023"

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
### Should we label selected points or just color them?
# Label them
max.overlaps = 1000
# Do not label them
# max.overlaps = 0

# Y-axis limits
# y.axis.min = -2.5
# y.axis.max = 2.5

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./draw_volcano_plot_from_limma_analysis.label_selected.R" }

input.filename = commandArgs()[6 - offset]
sampleA.name = commandArgs()[7 - offset]
sampleB.name = commandArgs()[8 - offset]
fdr.threshold = as.numeric(commandArgs()[9 - offset])
log2FC.threshold = as.numeric(commandArgs()[10 - offset])
genes.to.label.filename = commandArgs()[11 - offset]

if ( ! is.na(genes.to.label.filename) ) {
	message(paste("\nUsing limma output file *", input.filename, "* for input", sep=""))
	message(paste("Group A name is *", sampleA.name, "*", sep=""))
	message(paste("Group B name is *", sampleB.name, "*", sep=""))
	message(paste("FDR threshold is *", fdr.threshold, "*", sep=""))
	message(paste("Log2 fold change threshold is *", log2FC.threshold, "*", sep=""))
} else {
	message("\nDraw a volcano plot from limma output,")
	message("  coloring selected points that meet FDR and fold change thresholds,")
	message("  and labeling user-selected points.")
	message("USAGE: ", this.script, " inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold selectedGenesFile")
	message(paste(version, "\n"))
	quit()
}

# Read DESeq output file, and make the first column row names
# limma.stats = read.delim(input.filename, row.names=1)
###  Modified to permit gene matrices with duplicated gene symbols [GB, 7 Sep 2018]
# Input gene matrix (no row names)
limma.stats = read.delim(input.filename, header=T, sep="\t", check.names=F)

# Drop rows with a FDR of NA [Keep these - GB - 3/20/2019]
# limma.stats = limma.stats[! is.na(limma.stats[,6]),]

# Order genes by increasing p-values
limma.stats = limma.stats[order(limma.stats[,5], decreasing=FALSE), ]

# For plot labels
comparison = paste("log2 ( ", sampleA.name, " / ", sampleB.name, " )", sep="")
figure.title = paste("MA plot:", sampleA.name, "vs.", sampleB.name)

# Get IDs
gene.symbols = rownames(limma.stats)
p.value = limma.stats[,4]
# Get the FDR and then do -log10
fdr = limma.stats[,5]
minus.log10.fdr = -log(fdr, 10)
negLog10FDR = minus.log10.fdr

# Get the log2 fold change
log2FC = limma.stats[,2]
# log2 mean
log2.counts.mean = limma.stats[,1]

# Filter by FDR and fold change (for coloring)
# Filter genes with no FDR (FDR=NA) here - GB - 3/20/2019
filtered.rows =      ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold))
filtered.rows =      (((!is.na(fdr)) & (fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold)))
filtered.rows.up =   ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold) & log2FC > 0)
filtered.rows.down = ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold) & log2FC < 0)

# Use p-value threshold of 100th gene in this list
limma.stats.filtered = limma.stats[filtered.rows,]

# Look for genes to label in expression matrix
genes.to.label = unlist(read.delim(genes.to.label.filename, header=T, sep="\t", check.names=F))
num.original.genes.to.label = length(genes.to.label)
genes.to.label = intersect(genes.to.label, rownames(limma.stats))
genes.to.label.in.matrix = length(genes.to.label)
message(paste("\nOf", num.original.genes.to.label, "genes in your list for labeling, the expression matrix contains", genes.to.label.in.matrix, "of these."))

filtered.rows.labels = (rownames(limma.stats) %in% genes.to.label) & (log2.counts.mean > 0)
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

pdf.filename = paste("Volcano_plot", sampleA.name, "vs",  sampleB.name, fdr.threshold, log2FC.threshold, gsub(".txt", "", genes.to.label.filename), "pdf", sep=".")
pdf(pdf.filename, w=11, h=8.5, useDingbats=FALSE)

ggplot() +
# xlim(0, 20) + 
# ylim(y.axis.min, y.axis.max) + 
geom_point(aes(log2FC, negLog10FDR), color="darkgray", size=1) + 
geom_point(aes(log2FC[filtered.rows.up], negLog10FDR[filtered.rows.up]), color=up.color, size=1) + 
geom_point(aes(log2FC[filtered.rows.down], negLog10FDR[filtered.rows.down]), color=down.color, size=1) + 
geom_point(aes(log2FC[filtered.rows.up.labels], negLog10FDR[filtered.rows.up.labels]), color="red", size=3) + 
geom_point(aes(log2FC[filtered.rows.down.labels], negLog10FDR[filtered.rows.down.labels]), color="blue", size=3) + 
geom_vline(xintercept = 0, color="grey80") + 
labs(x = log2FC.label) +
labs(y = "-log10 (FDR)") + 
theme_classic(base_size = 16) +
ggtitle(figure.title) +
theme(plot.title = element_text(hjust = 0.5)) +
geom_text_repel(aes(log2FC[filtered.rows.labels], negLog10FDR[filtered.rows.labels]), label = gene.symbols[filtered.rows.labels], size=4, max.overlaps=100)


foo = dev.off()

message(paste("\nCreated plot!  See", pdf.filename, "\n"))

# Save details of your session
sink("sessionInfo.txt")
date()
print(sessionInfo(), locale=FALSE)
sink()
