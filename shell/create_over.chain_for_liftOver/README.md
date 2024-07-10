# Align two homologous genome regions and create a liftOver file to convert coordinates from one assembly/genome to another.  All important steps use code from UCSC Bioinformatics.


### Main command
```
./create_over.chain_for_liftOver.sh oldGenome newGenome chunkLength portNum
```

### Create over.chain for liftOver to convert from hg19 to hg38
```
./create_over.chain_for_liftOver.sh hg19.chr21.fa hg38.chr21.fa 3000 1234
```

### Create over.chain for liftOver to convert from hg38 to hg19
```
./create_over.chain_for_liftOver.sh hg38.chr21.fa hg19.chr21.fa 3000 4321
```

### Repo does not includes input files (hg19.chr21.fa hg38.chr21.fa) or most output files, which should look like those below.  The key file for 'liftOver' is the *over.chain one (which is present).

```
hg19.chr21_to_hg38.chr21_liftOver/
hg19.chr21_to_hg38.chr21_liftOver/hg19.chr21.fa
hg19.chr21_to_hg38.chr21_liftOver/hg38.chr21.fa
hg19.chr21_to_hg38.chr21_liftOver/hg19.chr21.2bit
hg19.chr21_to_hg38.chr21_liftOver/hg19.chr21.chrom.sizes
hg19.chr21_to_hg38.chr21_liftOver/hg38.chr21.2bit
hg19.chr21_to_hg38.chr21_liftOver/hg38.chr21.chrom.sizes
hg19.chr21_to_hg38.chr21_liftOver/hg38.chr21.split.fa
hg19.chr21_to_hg38.chr21_liftOver/hg38.chr21.lft
hg19.chr21_to_hg38.chr21_liftOver/hg19.chr21_to_hg38.chr21.split.psl
hg19.chr21_to_hg38.chr21_liftOver/hg19.chr21_to_hg38.chr21.psl
hg19.chr21_to_hg38.chr21_liftOver/hg19.chr21_to_hg38.chr21.chain
hg19.chr21_to_hg38.chr21_liftOver/chain
hg19.chr21_to_hg38.chr21_liftOver/chain/meta.tmp
hg19.chr21_to_hg38.chr21_liftOver/chain/chr21.chain
hg19.chr21_to_hg38.chr21_liftOver/chr.net
hg19.chr21_to_hg38.chr21_liftOver/hg19.chr21_to_hg38.chr21.over.chain
```
