#!/usr/bin/env Rscript

# Author: Bingbing Yuan
# July 1 2021

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./draw_scRNAseq_3D_UMAP.R" }


if (length(commandArgs()) < (8 - offset))
{
	message("Calculate the fragment size distribution for ATAC-seq\n.")
	message(paste("USAGE:   ", this.script, "input.bam output.pdf title_of_the_figure"))
	message(paste("Example: ", this.script, "GL1.bam sample_fragment_sizes.pdf Sample \n"))
	q()
}


input.bam = commandArgs()[6 - offset]
out.pdf = commandArgs()[7 - offset]
title = commandArgs()[8 - offset]


##############  end of user-defined variables  ##############

suppressMessages(library("ATACseqQC"))
pdf(out.pdf, w=11, h=8.5, useDingbats=FALSE)
fragSizeDist(input.bam, title)
dev.off()

message("Done")
