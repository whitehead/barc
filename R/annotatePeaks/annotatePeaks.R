#!/usr/bin/env Rscript

###
### Annotate regions (such as ChIP-seq peaks) by promoter, exon, etc. using 'ChIPseeker' (from Bioconductor)
### 
### 
### Added code to export the annotation details to text files 


offset=0
if (length(commandArgs()) < (9 - offset))
{
        message("\nAnnotate ChIPSeq Peaks.")
        message(paste("USAGE: ./annotatePeaks.R inputFile_Peaks distanceFromTSS_inbp genome[hg38|hg19|mm10|mm9] output"))
		message(paste("e.g. ./annotatePeaks.R MACS_peaks.txt 2000 hg38 MACS_peaks_annotated.pdf"))
	    message(paste("Note: This script/package may not correctly annotate for (very) broad peaks, e.g. enhancers, spanning several kb or longer."))
        q()
}

message("Loading ChIPseeker and other required packages ...")
suppressMessages(library(ChIPseeker))

regions.file = commandArgs()[6 - offset]
distFromTSS = as.numeric(commandArgs()[7 - offset])
genome = commandArgs()[8 - offset]
output = commandArgs()[9 - offset]

peak <- readPeakFile(regions.file)
txdb <-NULL;
peakAnno<-NULL;

if(genome == "hg38"){
	message("Loading TxDb.Hsapiens.UCSC.hg38.knownGene ...")
	suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	peakAnno <- annotatePeak(regions.file, tssRegion=c(-distFromTSS,distFromTSS),TxDb=txdb, annoDb="org.Hs.eg.db")

}else if (genome == "hg19"){
	message("Loading TxDb.Hsapiens.UCSC.hg19.knownGene ...")
	suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
	peakAnno <- annotatePeak(regions.file, tssRegion=c(-distFromTSS,distFromTSS),TxDb=txdb, annoDb="org.Hs.eg.db")

}else if (genome == "mm10"){
	message("Loading TxDb.Mmusculus.UCSC.mm10.knownGene ...")
	suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
	txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
	peakAnno <- annotatePeak(regions.file, tssRegion=c(-distFromTSS,distFromTSS),TxDb=txdb, annoDb="org.Mm.eg.db")

}else if (genome == "mm9"){
	message("Loading TxDb.Mmusculus.UCSC.mm9.knownGene ...")
	suppressMessages(library(TxDb.Mmusculus.UCSC.mm9.knownGene))
	txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
	peakAnno <- annotatePeak(regions.file, tssRegion=c(-distFromTSS,distFromTSS),TxDb=txdb, annoDb="org.Mm.eg.db")

}else{
	message(paste("Did not find",genome,sep=" "))
	q()
}



#export the annotation details
#Uncomment this if you would like to explore the Annotation obeject
#saveRDS(peakAnno, file = "peakAnno.RDS") 
#peakAnno@anno
output2 = gsub(".pdf", ".AnnotOut.txt",output)
write.table(as.data.frame(peakAnno@anno), file=output2, quote = F, row.names = F, sep = "\t")

#peakAnno@annoStat
output3 = gsub(".pdf", ".StatsAnnot.txt",output)
write.table(as.data.frame(peakAnno@annoStat), file=output3, quote = F, row.names = F, sep = "\t")

#peakAnno@detailGenomicAnnotation
#output4 = gsub(".pdf", ".DetailAnnot.txt",output)
#write.table(as.data.frame(peakAnno@detailGenomicAnnotation), file=output4, quote = F, row.names = F, sep = "\t")



pdf(output, w=11, h=8.5, useDingbats=FALSE)
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")
foo = dev.off()
