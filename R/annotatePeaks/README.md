# Script to annotate regions (such as ChIP-seq peaks) by promoter, exon, etc. using 'ChIPseeker' (from Bioconductor)

#### Required: ChIPseeker
####  Supports hg38, hg19, mm10, and mm9 coordinates
####  For more info about the package, see 
####  http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
#### Note: This script/package may not correctly annotate for (very) broad peaks, e.g. enhancers, spanning several kb or longer.
## Basic command

```
#USAGE: ./annotateChIPSeqPeaks.R inputFile_Peaks distanceFromTSStoDefinePromoters genome[hg38|hg19|mm10|mm9] output.pdf

./annotatePeaks.R MACS_peaks.txt 3000 hg38 MACS_peaks_annotated.pdf

```