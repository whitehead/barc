# Calculate the fragment size distribution with ATAC-seq alignment file
 The alignment file .bam should have been indexed (ie. .bam.bai should be in the same folder as .bam)
 sample input was downloaded from ATACseqQC package, a subset of ATAC-seq data from the original Greenleaf paper (Nature methods 10, 1213–1218 (2013))

### Requirements: the R package 'ATACseqQC'

### USAGE:    ./calculate_ATACseq_fragment_size_distribution.R input.bam output.pdf title_of_the_figure
```
./calculate_ATACseq_fragment_size_distribution.R GL1.bam sample_fragment_sizes.pdf Sample
```
