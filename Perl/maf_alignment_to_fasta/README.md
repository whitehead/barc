# Convert MAF format multiz alignment to fasta format

The output fasta include gaps so can be directly inspected with an alignment viewer (like Jalview).

This works for multiz 100-way alignments as of June 2021.

```
# Usage:   ./maf_alignment_to_fasta.hg38.pl mafFile genomeBuild[UCSC] geneName strand > out.fa
# Usage:   ./maf_alignment_to_fasta.hg19.pl mafFile genomeBuild[UCSC] geneName strand > out.fa

# Reference species is hg38
./maf_alignment_to_fasta.hg38.pl My_region.hg38.multiz.maf hg38 Region_ABC - > Region_ABC.hg38.multiz.fa

# Reference species is hg19
./maf_alignment_to_fasta.hg19.pl My_region.hg19.multiz.maf hg19 Region_ABC + > Region_ABC.hg19.multiz.fa
```
