#! /usr/bin/env Rscript

# Make a volcano plot to highlight genes within a specific pathway after
# doing GSEA.  Requires output from DESeq2 and GSEA.

# Version 1.0, 22 Nov., 2022
# Troy Whitfield - Bioinformatics and Research Computing, Whitehead Institute

argv <- commandArgs()
if('--args' %in% argv){
   nignore <- which (argv == "--args")
   argv <- argv[-(1:nignore)]
}
if((length (argv) != 5) | (argv[1] == "/usr/lib/R/bin/exec/R")) {
   message("Draw a gene-set level volcano plot")
   message("USAGE: ./drawGSEAvolcano.R DESeq2.out GSEA.out MSigDB.in contrastLabel output.pdf")
   message("Example: ./drawGSEAvolcano.R DESeq2.csv GSEA.tsv MSigDB.txt fold-change GSEAvolcano.pdf\n")

   # stop("Follow the expected input format above.")
   q()
}

message("Drawing a volcano plot to highlight genes in a specific pathway.")

suppressMessages(require("ggplot2"))
suppressMessages(require("ggrepel"))
suppressMessages(require("Cairo")) # Prettier fonts for plotting.

# Assign variables from input.
inres<-argv[1]
inGSEA<-argv[2]
ingset<-argv[3]
contrastLab<-as.character(argv[4])
outfile<-argv[5]

# Read in DESeq2 output.
res<-read.table(inres,sep=",",header=TRUE)

# Set limits for the volcano plot.
xpad<-1
ypad<-1
values<-c(0.25,0.5,1,2,3,4,5,10,15,20,25,30,50,100)
yupper<--log(min(res$pvalue,na.rm=TRUE),base=10)+ypad
ylower<-0
xupper<-ceiling(max(abs(res$log2FoldChange))+xpad)
xlower<--xupper
dx<-round((xupper-xlower)/5)
xstep<-values[which(abs(values-dx) == min(abs(values-dx)))]
xmax<-values[which(abs(values-xupper) == min(abs(values-xupper)))]
xmin<--xmax

# Read in MSigDB input.
gset<-read.table(ingset,header=FALSE)

# Select genes to highlight at the colour-only level and with explicit symbols/names too.
resCrd<-as.data.frame(cbind(res$log2FoldChange,-log10(res$pvalue)))
colnames(resCrd)<-c("log2FoldChange","log10pvalue")
pathCrd<-as.data.frame(cbind(res$log2FoldChange[which(res$symbol %in% gset$V1)],-log10(res$pvalue[which(res$symbol %in% gset$V1)])))
pathSymbol<-res$symbol[which(res$symbol %in% gset$V1)]
colnames(pathCrd)<-c("log2FoldChange","log10pvalue")
pathHits<-as.data.frame(cbind(pathCrd,res$symbol[which(res$symbol %in% gset$V1)]))
colnames(pathHits)<-c("log2FoldChange","log10pvalue","symbol")
# Start with the GSEA "leading edge" or "core enrichment" genes.
leadingEdge<-read.table(inGSEA,header=TRUE,sep="\t")
fcHits<-leadingEdge$SYMBOL[which(leadingEdge$CORE.ENRICHMENT == "Yes")]
fcHits<-match(fcHits,pathSymbol)
if (mean(pathCrd[fcHits,1]) > 0){
   fcHits<-fcHits[which(pathCrd[fcHits,1] > 1)] # Selects only leading edge with FC > 2.
}else{
   fcHits<-fcHits[which((pathCrd[fcHits,1] < 0) & (abs(pathCrd[fcHits,1]) > 1))]
}
pathHits<-pathHits[fcHits,]

sink("/dev/null") # Suppress default reporting.
cairo_pdf(outfile,width = 8.5, height = 8.0, onefile = TRUE, family = "Helvetica", symbolfamily="NimbusSans", bg="white", pointsize=14)
ggplot(data = resCrd, # Original data  
       aes(x = log2FoldChange, y = log10pvalue)) + 
  geom_point(colour = "grey", alpha = 0.5) +
  geom_point(data = pathCrd,
             size = 2,
             shape = 21,
             fill = "red",
             colour = "black") +
  geom_text_repel(data = pathHits,   
                   aes(label = symbol),
                   force = 20,
                   nudge_x = 0.5,
		   nudge_y = 0.5,
		   max.overlaps = 20,
		   size = 2.5) +     
  scale_x_continuous(breaks = seq(xmin, xmax, xstep),
             limits = c(xlower, xupper)) +
  scale_y_continuous(breaks = c(seq(0, yupper, 5)),     
             limits = c(0, yupper)) +
  geom_hline(yintercept = -log10(0.00001),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +	     
  labs(title = "",
             x = bquote(log[2](.(contrastLab))),
             y = expression(-log[10](p))) +
  theme(panel.grid.major = element_blank(),
  	     panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
	     axis.text.y = element_text(size=14),
             axis.text.x = element_text(size=14),
             axis.title.y = element_text(size=14),
             axis.title.x = element_text(size=14))
dev.off()
sink()
