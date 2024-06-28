# Translate each sequence in a multiple-sequence file in six frames, keep the longest ORF for each gene sequence, and print to STDOUT.

Required: BioPerl (Bio::SeqIO) and the EMBOSS suite (sixplack) 

### USAGE: ./sixpack_iterate.pl inputGenesFile > outputProteinsFile
```
./sixpack_iterate.pl DNA.fa > DNA.sixpack_output.fa
```
