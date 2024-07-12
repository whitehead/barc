## Adjust for batch effects in RNA-seq data using ComBat-seq

Convert a raw counts matrix with batch and condition
information to a batch-adjusted counts matrix using
ComBat-seq.

USAGE: `./remove_batch_effects_counts.R counts samples covariates adjCounts`

Example: `./remove_batch_effects_counts.R counts.txt samples.txt default adjCounts.txt`

The two input files should be tab delimited text.
The counts matrix file should be samples x genes (columns x rows).
The annotation file should have samples listed in colunm 1,
batches listed in column 2 and biological condition listed in column 3.
The samples listed in the counts matrix should match those from the
sample annotations.

The 'covariates' keyword controls whether the condition annotations
from the annotation file should be used to explain some of the
between-sample variance in the counts matrix.  This keyword can be
set to 'true' or 'false', but the default is to include such features.

If you use these results in a publication, please cite
Y. Zhang et al., _NAR Genomics and Bioinformatics_ **2**, lqaa078 (2020)

### Required packages
sva (3.46.0)