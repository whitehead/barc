#!/usr/bin/env Rscript

###
###  Use DESeq2 to calculate FPKM and normalized FPKM from raw counts
###  George Bell - Bioinformatics and Research Computing, Whitehead Institute
###
###  USAGE: counts_to_FPKM_using_DESeq2.R inputCounts geneLengths OutputFile
###  EX: counts_to_FPKM_using_DESeq2.R Mouse_RNAseq_counts.txt Mouse.medianGeneLength.txt Mouse_RNAseq.FPKM.txt
###
###
###  Version 1.0: 16 June 2017 [GWB]
###
###  Version 1.1: Updated 26 Feb 2018
###    Round output values to significant.figures [initially set to 6]
###  Version 1.2: Updated 23 Oct 2019
###    Add TPM to output values
###  Version 1.21: Updated 1 June 2020 -- Change TPM calculation so it works even if some genes have 'NA' FPKMs (like if they're missing lengths)
###  Add "check.names=F" to reading of counts matrix
###  Version 1.22: Fix bug from fpkm() -- 28 Apr 2022
###  Version 1.3: Allow redundant symbols; make them unique -- 3 Nov 2022

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./counts_to_FPKM_TPM_using_DESeq2.R" }

if (length(commandArgs()) < (7 - offset))
{
	message("\nConvert a matrix of raw counts to RPKMs with DESeq2.")
	message(paste("USAGE:   ", this.script, "inputCounts geneLengths OutputFile"))
	message(paste("Example: ", this.script, "Mouse_RNAseq_counts.txt Mouse.medianGeneLength.txt Mouse_RNAseq.FPKM.TPM.txt\n"))
	q()
}

# Load library without messages
message("\nLoading DESeq2 and other required R packages ....")
suppressMessages(library(DESeq2))
suppressMessages(library(genefilter))
suppressMessages(library(RColorBrewer))
suppressMessages(library(lattice))
suppressMessages(library(ggrepel))

significant.figures = 6

message(paste("Running DESeq2 version", packageDescription("DESeq2")$Version))

# First argument is input file
input.filename = commandArgs()[6 - offset]

# 2-row file (no header) of gene lengths
input.gene.lengths.filename = commandArgs()[7 - offset]

# Output file (of FPKMs)
output.filename = commandArgs()[8 - offset]

message(paste("Raw counts read from", input.filename))
message(paste("Gene lengths read from", input.gene.lengths.filename))

# Read matrix of counts
# Change method so duplicate row names will be tolerated.  28 January 2021
# counts.matrix = read.delim(data.matrix.file, row.names=1, check.names=F)
# Permit gene matrices with duplicated gene symbols [GB, 7 Sep 2018]	Added "check.names=FALSE" on 28 Jan 2021
counts.matrix = read.delim(input.filename, header=T, sep="\t", check.names=FALSE)
# Get total counts per gene (for sorting redundant rows)
counts.sum.by.gene = apply(counts.matrix[,2:ncol(counts.matrix)], 1, sum)
# Order matrix by decreasing counts.sum.by.gene, adding a numeric suffix to each symbol we've already seen.
gene.symbols.original = as.character(counts.matrix[,1])
# gene.symbols.unique = .Internal(make.names(gene.symbols.original, unique=TRUE))
gene.symbols.unique = make.unique(gene.symbols.original)
o = order(gene.symbols.original != gene.symbols.unique, counts.sum.by.gene, decreasing=TRUE)
gene.symbols.unique[o] = make.unique(gene.symbols.unique[o])
# Make row names from unique gene symbols.
rownames(counts.matrix) = gene.symbols.unique
# Drop the gene symbol (first column) since we have them as row names
counts = counts.matrix[,-1]
# Force counts values into integers [GB -- 3 Nov 2022]
counts = round(counts, 0)

# Load gene lengths (which can have header or not, but let's assume not)
gene.lengths = read.delim(input.gene.lengths.filename, row.names=1, header=F)

# Make a DESeqDataSet
samples = factor(colnames(counts))
# Use a design stating that there are no replicates
# If there is replication, replace 'design' with the real design
dds = DESeqDataSetFromMatrix(countData=counts, colData=DataFrame(samples), design=~1)

# These commands make no difference in the output
# dds = DESeq(dds, fitType="mean")
# dds = DESeq(dds, fitType="parametric")

# Add the gene lengths to the data structure
gene.lengths.ordered = gene.lengths[rownames(dds),]
# mcols(dds)$basepairs = gene.lengths.ordered	# Fixed on 28 Apr 2022
mcols(dds)$basepairs = as.numeric(gene.lengths.ordered)

# MAIN COMMANDS: calculate FPKMs
# Calculate typical FPKMs
fpkm.raw = fpkm(dds, robust = FALSE)
# Calculate normalized FPKMs
# (use size factors to normalize rather than taking the column sums of the raw counts)
fpkm.norm = fpkm(dds, robust = TRUE)

# Calculate TPM from raw FPKM
# Get total FPKM for each sample (row sums)
# Keep going even if we have NAs
# fpkm.raw.sums = apply(fpkm.raw, 2, sum)
fpkm.raw.sums = apply(fpkm.raw, 2, sum, na.rm=T)
# Apply FPKM => TPM formula
tpm = t(t(fpkm.raw * 1e6) / fpkm.raw.sums)

colnames(fpkm.raw) = paste(colnames(dds), "FPKM.raw", sep=".")
colnames(fpkm.norm) = paste(colnames(dds), "FPKM.norm", sep=".")
colnames(tpm) = paste(colnames(dds), "TPM", sep=".")

# Print output (FPKMs and TPMs)
message(paste("Rounding FPKMs and TPMs to", significant.figures, "significant figures....\n"))
output.table = cbind(rownames(dds), signif(fpkm.raw, significant.figures), signif(fpkm.norm, significant.figures), signif(tpm, significant.figures))
colnames(output.table)[1] = "Feature.ID"
write.table(output.table, file=output.filename, sep="\t", quote=F, row.names=F)

message(paste("\nAll done!  See", output.filename, "for output.\n"))
