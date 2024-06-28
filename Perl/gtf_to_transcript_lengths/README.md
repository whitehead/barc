# From a GTF file, extract the length of each transcript, creating a file of gene ID, transcript ID, and length.
## Then summarize transcript lengths for each gene, creating output files of mean (*.meanGeneLength.txt) and median (*.medianGeneLength.txt) lengths.

Required: Statistics::Lite, File::Basename


### USAGE: ./gtf_to_transcript_lengths.pl genes.gtf > transcript_lengths.txt
```
./gtf_to_transcript_lengths.pl My_genes.gtf >| My_genes.transcript_lengths.txt
```

Note that these summarized gene lengths can be used for calculations like FPKM and RPKM.
