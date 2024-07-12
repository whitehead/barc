#!/usr/bin/env Rscript

###
###  George Bell
###  Bioinformatics and Research Computing, Whitehead Institute
###
###  Split a multiple-comparison limma output into separate files by comparison
###  Input: Output from Run_all_vs_all_limma_differential_expression.R
###
###  Version 1.0: 22 December 2021
###  Version 1.1: 20 January 2023 -- Add input for log2FCthreshold
###

offset = 0
this.script = "./split_limma_output_by_comparison.R"
limma.all.comparisons.file = commandArgs()[6 - offset]
fdr.threshold.for.coloring = as.numeric(commandArgs()[7 - offset])
output.dir = "limma_out_by_comparison"
log2FC.threshold.for.coloring = as.numeric(commandArgs()[8 - offset])

if ( ! is.na(fdr.threshold.for.coloring) ) {
	message(paste("\nProcessing limma output file *", limma.all.comparisons.file, "* for input ...", sep=""))
} else {
	message("\nSplit a multiple-comparison limma output into separate files by comparison")
	message("  and draw MA and volcano plots for each comparison.")
	message("USAGE: ", this.script, " inputFile fdrThreshold log2FCthreshold\n")
	quit()
}

if (! dir.exists(output.dir))	{ dir.create(output.dir) }

# Read matrix of values [Don't "fix" names -- 20 Dec 2021]
limma.all.comparisons = read.delim(limma.all.comparisons.file, row.names=1, check.names=F)

num.comparisons = (ncol(limma.all.comparisons) - 3) / 4

A.colNum = 1
F.colNum = ncol(limma.all.comparisons) - 1
F.p.value.colNum = ncol(limma.all.comparisons)

comparisons = gsub(" - ", ".vs.", colnames(limma.all.comparisons)[2:(num.comparisons+1)])
comparisons = gsub("Coef.", "", comparisons)

for (i in 1:num.comparisons)
{
	comparison = comparisons[i]
	Coef.colNum = i + 1
	t.colNum = num.comparisons*1 + i + 1
	p.value.colNum = num.comparisons*2 + i + 1
	p.value.adj.colNum = num.comparisons*3 + i + 1

	columns.this.comparison = colnames(limma.all.comparisons)[c(A.colNum, Coef.colNum, t.colNum, p.value.colNum, p.value.adj.colNum, F.colNum, F.p.value.colNum)]
	limma.this.comparison = limma.all.comparisons[,c(A.colNum, Coef.colNum, t.colNum, p.value.colNum, p.value.adj.colNum, F.colNum, F.p.value.colNum)]
	
	# Sort by increasing p-value and then by abs(t)
	limma.this.comparison.sorted = limma.this.comparison[order(limma.this.comparison[,4], -1*abs(limma.this.comparison[,3])),]
	
	output.file = paste(output.dir, paste(comparison, "limma.txt", sep="."), sep="/")
	write.table(limma.this.comparison.sorted, sep="\t", quote=F, file=output.file)
	
	print(message("Colored points in plots have FDR <", fdr.threshold.for.coloring))
	
	first.second = unlist(strsplit(comparison, ".vs.", fixed=T))
	message(paste("Comparing", first.second[1], "and", first.second[2], "..." ))
	MA.cmd = "/nfs/BaRC_Public/BaRC_code/R/limma/draw_MA_plot_from_limma_analysis.R"
	volcano.cmd = "/nfs/BaRC_Public/BaRC_code/R/limma/draw_volcano_plot_from_limma_analysis.R"
	arguments = paste(output.file, first.second[1], first.second[2], fdr.threshold.for.coloring, log2FC.threshold.for.coloring, sep=" ")
	message(paste(MA.cmd, arguments, sep=" "))	
	system(paste(MA.cmd, arguments, sep=" "))
	message(paste(volcano.cmd, arguments, sep=" "))	
	system(paste(volcano.cmd, arguments, sep=" "))
}
