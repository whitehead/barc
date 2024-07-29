#!/usr/bin/env Rscript

###
###  Draw an enrichment plot from a fgsea output file
###
###  George Bell - Whitehead Institute BaRC
###  Version 1.0: 9 April 2020
###  Input file: Output Excel file from fgsea
###    with these header IDs (in any order): pathway, padj, NES
###  Draw two enrichment plots with 
###    enriched terms printed starting with most extreme (negative or positive) NES
###    gene set size reflected in point size, and 
###    adjusted p-value (padj) reflected by point color.
###
###  12.2.2021. not using default pdf font because Illustrator converts circles to squares: (added by BY) 
### 		https://stackoverflow.com/questions/9992275/ggplot2-pdf-import-in-adobe-illustrator-missing-font-adobepistd
###  allow choose by FDR cutoff
###  change 'aes_string(x = "NES", y = "pathway")' to 'aes(x = NES, y = pathway)' [to prevent warning; 18 Dec 2023]

# What's the largest size point we want to print?
max.point.cex = 6
# Set size of PDF output file (in inches)
pdf.width.inches = 11
pdf.height.inches = 8.5
# Drop initial word in pathway?
drop.pathway.prefix = FALSE

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./draw_enrichment_figures.R" }

input.filename = commandArgs()[6 - offset]
input.file.type = commandArgs()[7 - offset]
fdr.cutoff = as.numeric(commandArgs()[8 - offset])
num.terms.to.plot = as.numeric(commandArgs()[9 - offset])
max.term.length = as.numeric(commandArgs()[10 - offset])

# input.filename = "Res_DE1_7_Hallmarksigpath.txt"
# input.file.type = "fgsea"
# fdr.cutoff = 0.01
# num.terms.to.plot = 10
# max.term.length = 50

# input.filename = "paired_normShrink_logfc.GseaPreranked.1683652035539/gsea_report_for_na_pos_1683652035539.tsv"
# input.file.type = "GSEA"
# fdr.cutoff = 0.15
# num.terms.to.plot = 4
# max.term.length = 50


if ( ! is.na(max.term.length) ) {
	message(paste("\nDraw enrichment figures using *", input.filename, "* for input", sep=""))
	message(paste("Terms with FDR less than *", fdr.cutoff, "* will be plotted.\n", sep=""))
} else {
	message("\nDraw enrichment figures from a GSEA/fgsea output file")
	message("USAGE: ", this.script, " inputFile fileType[GSEA|fgsea] fdr.cutoff num.terms.to.plot max.term.length [max point cex, default=10]")
	message("Note: For GSEA files, provide one file (pos or neg) as an argument, but both files are needed.\n")
	quit()
}

if (length(commandArgs()) >= 11 - offset)
{
	max.point.cex = as.numeric(commandArgs()[11 - offset])
	message(paste("Max value for point size in R is now set to *", max.point.cex, "*.", sep=""))
}

# For plotting
library(ggplot2)

# assume no terms to be printed out unless finding gene sets fitting the cutoffs
num.pos.NES = 0
num.neg.NES = 0

if (input.file.type == "fgsea")
{
	# Convert xlsx file into tibble
	go = read.delim(input.filename)
	
	# Sort by decreasing NES and get top rows (and then by p-value, in case of ties)
	go.pos.NES = go[order(go$NES, go$pval, decreasing=c(TRUE, FALSE)),]
	# Sort by increasing NES and get top rows (and then by p-value, in case of ties)
	go.neg.NES = go[order(go$NES, go$pval, decreasing=c(FALSE, FALSE)),]
	
	go.pos.NES.filtered = go.pos.NES[which ((go.pos.NES$NES > 0) & (go.pos.NES$padj < fdr.cutoff)),]
	go.neg.NES.filtered = go.neg.NES[which ((go.neg.NES$NES < 0) & (go.neg.NES$padj < fdr.cutoff)),]
	
	num.pos.NES = dim(go.pos.NES.filtered)[1]
	num.neg.NES = dim(go.neg.NES.filtered)[1]
	
	if (num.pos.NES <= 0)
	{
	  message(paste("NOTE: No term for positive NES values with FDR less than", fdr.cutoff))
	} else {
	  
	  if (num.pos.NES < num.terms.to.plot)
	  {
	    message(paste("NOTE: num.terms.to.plot is being changed to *", num.pos.NES, "* for positive NES values."))
	    num.terms.to.plot.pos = num.pos.NES
	  } else {
	    num.terms.to.plot.pos = num.terms.to.plot
	  }
	  # Sort by decreasing NES and get top rows (and then by p-value, in case of ties)
	  go.pos.NES.filtered2 = head(go.pos.NES.filtered[order(go.pos.NES.filtered$NES, go.pos.NES.filtered$pval, decreasing=c(TRUE, FALSE)),], num.terms.to.plot.pos)	  
	}
	  
	if (num.neg.NES <= 0)
	{
	  message(paste("NOTE: No term for negative NES values with FDR less than ", fdr.cutoff))
	}
	else {
	  if (num.neg.NES < num.terms.to.plot)
	  {
	    message(paste("NOTE: num.terms.to.plot is being changed to *", num.neg.NES, "* for negative NES values."))
	    num.terms.to.plot.neg = num.neg.NES
	  } else {
	    num.terms.to.plot.neg = num.terms.to.plot
	  }
	  # Sort by increasing NES and get top rows (and then by p-value, in case of ties)
	  go.neg.NES.filtered2 = head(go.neg.NES.filtered[order(go.neg.NES.filtered$NES, go.neg.NES.filtered$pval, decreasing=c(FALSE, FALSE)),], num.terms.to.plot.neg)
	}



} else if (input.file.type == "GSEA") {
	if (grepl('_pos_', input.filename))
	{
		input.filename.pos = input.filename
		input.filename.neg = gsub("_pos_", "_neg_", input.filename)
	} else {
		input.filename.neg = input.filename
		input.filename.pos = gsub("_neg_", "_pos_", input.filename)	
	}

	# Read the xls file (which is really tab-delimited text)
	go.pos = read.delim(input.filename.pos)
	message(paste("Reading GSEA file ", input.filename.pos))
	go.neg = read.delim(input.filename.neg)
	message(paste("Reading GSEA file ", input.filename.neg))
	# NAME	GS<br> follow link to MSigDB	GS DETAILS	SIZE	ES	NES	NOM p-val	FDR q-val	FWER p-val	RANK AT MAX	LEADING EDGE
	# Adjust column names to match those of fgsea file format
	colnames(go.pos)[1] = "pathway"
	colnames(go.pos)[4] = "size"
	colnames(go.pos)[7] = "pval"
	colnames(go.pos)[8] = "padj"
	colnames(go.neg)[1] = "pathway"
	colnames(go.neg)[4] = "size"
	colnames(go.neg)[7] = "pval"
	colnames(go.neg)[8] = "padj"

	# Sort by decreasing NES and get top rows (and then by p-value, in case of ties)
	go.pos.NES = go.pos[order(go.pos$NES, go.pos$pval, decreasing=c(TRUE, FALSE)),]
	# Sort by increasing NES and get top rows (and then by p-value, in case of ties)
	go.neg.NES = go.neg[order(go.neg$NES, go.neg$pval, decreasing=c(FALSE, FALSE)),]

	go.pos.NES.filtered = go.pos.NES[which ((go.pos.NES$NES > 0) & (go.pos.NES$padj < fdr.cutoff)),]
  go.neg.NES.filtered = go.neg.NES[which ((go.neg.NES$NES < 0) & (go.neg.NES$padj < fdr.cutoff)),]
  
  num.pos.NES = dim(go.pos.NES.filtered)[1]
  num.neg.NES = dim(go.neg.NES.filtered)[1]
  

  if (num.pos.NES <= 0)
  {
    message(paste("NOTE: No term for positive NES values with FDR less than", fdr.cutoff))
  } else {
    if (num.pos.NES < num.terms.to.plot)
    {
      message(paste("NOTE: num.terms.to.plot is being changed to *", num.pos.NES, "* for positive NES values."))
      num.terms.to.plot.pos = num.pos.NES
    } else {
      num.terms.to.plot.pos = num.terms.to.plot
    }
    # Sort by decreasing NES and get top rows (and then by p-value, in case of ties)
    go.pos.NES.filtered2 = head(go.pos.NES.filtered[order(go.pos.NES.filtered$NES, go.pos.NES.filtered$pval, decreasing=c(TRUE, FALSE)),], num.terms.to.plot)

  }
  
  if (num.neg.NES <= 0)
  {
    message(paste("NOTE: No terms for negative NES values with FDR less than ", fdr.cutoff))
    
  }
  else {
    if (num.neg.NES < num.terms.to.plot)
    {
      message(paste("NOTE: num.terms.to.plot is being changed to *", num.neg.NES, "* for negative NES values."))
      num.terms.to.plot.neg = num.neg.NES
    } else {
      num.terms.to.plot.neg = num.terms.to.plot
    }
    # Sort by increasing NES and get top rows (and then by p-value, in case of ties)
    go.neg.NES.filtered2 = head(go.neg.NES.filtered[order(go.neg.NES.filtered$NES, go.neg.NES.filtered$pval, decreasing=c(FALSE, FALSE)),], num.terms.to.plot)
  }
  
  
} else {
	message("ERROR: Input fileType needs to GSEA or fgsea !")
	exit()
}

# Function to reverse order of rows in a matrix, data frame, or tibble
reverseRowOrder = function(df) { df[nrow(df):1,] }
# Make a color scale by p-value (or something else)
map2color = function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)] }

gave.truncation.warning = 0
gave.prefix.trimming.warning = 0

# Make the plots
fileBase = gsub(".xls.?", "", input.filename)
fileBase = gsub(".tsv$", "", fileBase)
fileBase = gsub(".txt$", "", fileBase)

if (input.file.type == "GSEA" && grepl('_pos_', fileBase)) {
	fileBase = gsub("_pos_", "_neg+plus_", fileBase)
} else if (input.file.type == "GSEA") {
	fileBase = gsub("_neg_", "_neg+plus_", fileBase)
}

drawEnrichmentPlot = function(enrichment.table, max.cex = 5, low.color="red", high.color="blue", log.point.size=FALSE, num.datasets=1, sort.type="as.is")
{
	# Required fields of enrichment.table
	#   pathway, pval, padj, NES
	
	my.palette = colorRampPalette( c(low.color, high.color) )( 50 )
	go.term.colors = map2color(enrichment.table$padj, my.palette)
	
	if (num.datasets > 1)
	{ max.cex = max.cex / 2 }

	# Clean up selected categories
	# Get first word in pathway
	go.terms = tolower(enrichment.table$pathway)
	pathway.type = sub("_.+", "", go.terms[1])
	
	# Drop initial word in pathway?
	if (drop.pathway.prefix == TRUE)
	{
		go.terms = gsub(paste(pathway.type, "_", sep=""), "", go.terms)
	}
	
	# Convert underscores to spaces
	go.terms = gsub("_", " ", go.terms)

	if (gave.prefix.trimming.warning == 0 & drop.pathway.prefix) {
		message(paste("Trimming '", pathway.type, "_' from the beginning of pathway names.", sep=""))
		# This syntax changes the global variable
		gave.prefix.trimming.warning <<- 1
	}

	if (max(nchar(go.terms)) > max.term.length & gave.truncation.warning == 0) {
		message(paste("!!! Terms that are longer than", max.term.length, "are being truncated !!!"))
		# This syntax changes the global variable
		gave.truncation.warning <<- 1
	}

	# go.terms = substr(go.terms, 0, max.term.length)
	
	for (i in 1:length(go.terms))
	{
		add.dots = FALSE
		if (nchar(go.terms[i]) > max.term.length)
		{ add.dots = TRUE }
		go.terms[i] = substr(go.terms[i], 0, max.term.length)
		if (add.dots)
		{ go.terms[i] = paste(go.terms[i], "...", sep="") }
	}
	
	# Put the modified terms back into the matrix
	enrichment.table$pathway = go.terms
	
	# Sort so most extreme is plotted first)
	# enrichment.table = enrichment.table[order(abs(enrichment.table$NES), decreasing=TRUE),]
	# Need to order pathway factors so plot is printed in correct order (and not sorted by pathway)
	if (sort.type == "as.is") {
		# print(1:nrow(enrichment.table))	
		enrichment.table$pathway = factor(enrichment.table$pathway, levels = enrichment.table$pathway[order(1:nrow(enrichment.table))])	
	} else if (sort.type == "dec") {
		enrichment.table$pathway = factor(enrichment.table$pathway, levels = enrichment.table$pathway[order(abs(enrichment.table$NES), decreasing=FALSE)])
	} else if (sort.type == "inc") {
		
		enrichment.table$pathway = factor(enrichment.table$pathway, levels = enrichment.table$pathway[order(abs(enrichment.table$NES), decreasing=TRUE)])	
	}

	if (log.point.size == FALSE) {
		gene.set.length = max.cex * max.cex*enrichment.table$size/max(enrichment.table$size)
	} else {
		gene.set.length = max.cex * log2(max.cex*enrichment.table$size)/max(log2(enrichment.table$size))
	}

	FDR = signif(enrichment.table$padj, 3)
	
	# calculate the size of gene size to minimize the distortion
	size.min = min(enrichment.table$size)
	size.max = max(enrichment.table$size)
	size.ratio = as.integer(size.max/size.min)

	enrichment.plot = 
	ggplot(data = data.frame(enrichment.table), aes(x = NES, y = pathway)) + 
	geom_point(aes(size = size, color = FDR)) + 
	# use scale_size_area is better to control point size with count as input than scale_size_continuous
	scale_size_area(name="gene\nset\nsize", max_size=max.cex) + 
	scale_shape_identity() +
	scale_colour_gradientn(colours=my.palette, name="FDR") + 
	guides(size = guide_legend(order = 1)) + 	# print this legend first
	xlab("NES") + 
	ylab(paste(toupper(pathway.type), "terms")) + 
	ggtitle("Most enriched terms") +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 12))

	print(enrichment.plot)
}

#message(paste("Number of terms positively enriched is", num.pos.NES, sep=" "))
#message(paste("Number of terms negatively enriched is", num.neg.NES, sep=" "))

if(num.pos.NES <= 0 & num.neg.NES <=0 ) {
  message(paste ("No terms with FDR less than ", fdr.cutoff, sep=""))
} else {
  pdf.file = paste(fileBase, ".", "fdr", fdr.cutoff, ".", num.terms.to.plot, ".", max.term.length, ".enrichment.plots.pdf", sep="")
  pdf(pdf.file, w=pdf.width.inches, h=pdf.height.inches, useDingbats=FALSE)

  if(num.pos.NES>0) {
    drawEnrichmentPlot(go.pos.NES.filtered2, max.cex=max.point.cex, sort.type="dec")
  }
  if(num.neg.NES>0) {
    drawEnrichmentPlot(go.neg.NES.filtered2, max.cex=max.point.cex, sort.type="dec")
  }
  # Put both together
  if(num.pos.NES>=1 & num.neg.NES >=1) {
    drawEnrichmentPlot(rbind(go.pos.NES.filtered2, go.neg.NES.filtered2), max.cex=max.point.cex, sort.type="as.is")
  }
  # Finish figure
  foo = dev.off()

  message(paste("All done!  See ", pdf.file, " for output.\n", sep=""))
}