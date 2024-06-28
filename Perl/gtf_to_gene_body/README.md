# Convert a GTF file to a file of gene bodies as features (GTF or BED output)

### Create GTF output, with one gene body per gene
```
./gtf_to_gene_body_gtf.pl MyGenes.gtf gene | sort -k1,1 -k4,4n >| MyGenes.gene.gene_bodies.gtf
```

### Create GTF output, with one gene body per transcript
```
./gtf_to_gene_body_gtf.pl MyGenes.gtf transcript | sort -k1,1 -k4,4n >| MyGenes.transcript.gene_bodies.gtf
```

### Create BED output, with one gene body per gene
```
./gtf_to_gene_body_bed.pl MyGenes.gtf gene | sort -k1,1 -k4,4n >| MyGenes.gene.gene_bodies.bed
```

### Create BED output, with one gene body per transcript
```
./gtf_to_gene_body_bed.pl MyGenes.gtf transcript | sort -k1,1 -k4,4n >| MyGenes.transcript.gene_bodies.bed
```