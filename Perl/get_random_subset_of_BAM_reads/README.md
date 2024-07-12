# Extract a random subset of mapped reads from a SAM or BAM file

```
# USAGE: ./get_random_subset_of_BAM_reads.pl SAMfile numReadsWanted > out.sam

# Filter BAM file (SAM output)
./get_random_subset_of_BAM_reads.pl mapped_reads.bam 10 > mapped_reads.10.sam

# Filter BAM file (BAM output)
./get_random_subset_of_BAM_reads.pl mapped_reads.bam 10 | samtools view -bS > mapped_reads.10.bam

# Filter SAM file (SAM output)
./get_random_subset_of_BAM_reads.pl mapped_reads.sam 25 > mapped_reads.25.sam

# Filter SAM file (BAM output)
./get_random_subset_of_BAM_reads.pl mapped_reads.sam 25 | samtools view -bS > mapped_reads.25.bam
```
