# Read compressed fastq file and count guides


### USAGE: ./read_count_CRISPR_guides.pl guidesFile Reads.fq.gz output


```
read_count_CRISPR_guides.pl My_guides_list.txt Reads_sample.fq.gz counts

```

In addition to calculating the counts for each guide, this script generates an additional file, _other_abundant_guides.txt, which contains sequences that weren't listed in the guide list. This can be useful for quality control purposes.

