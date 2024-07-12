# Convert bed file to GTF/GFF file

### Convert BED to GFF
./bed2gff.pl Input.bed Ensembl gene_body > Input.gff

### Convert BED to GTF
Note that 4th (name) field needs to include gene ID and transcript ID (semicolon-delimited, in that order) to include different gene and transcript IDs
```
./bed2gtf.pl Regions_for_GTF.bed Ensembl gene_body > Regions_for_GTF.gtf
```
### Convert BED to GTF
If the 4th field contains something else, this will be used to get identical gene and transcript IDs.
```
./bed2gtf.pl Input.bed Ensembl gene_body > Input.gtf
```
