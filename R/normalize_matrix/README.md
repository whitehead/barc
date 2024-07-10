# Normalize by 75th percentile

Required packages: limma, lattice

### Run quantile normalization
```
./normalize_matrix.R Sample_input.txt quantile Sample_input.quantileNorm.txt nothing
```

### Normalize to the 75th percentile
```
./normalize_matrix.R Sample_input.txt percentile=0.75 Sample_input.norm.75thpctile.txt nothing
```

### Run quantile normalization but using different code (quantile_normalize.R) wiht somewhat difference interface
```
./quantile_normalize.R Sample_input.txt Sample_input.quantileNorm2.txt
```

### Quantile normalize matrix (with quantile_normalize.R) but add noise to get truly matching quantiles (for cases where there might be a lot of 0 values).  Matrices with varying numbers of 0s may not otherwise end up with truly equal distributions.
```
./quantile_normalize.R Sample_input.txt Sample_input.quantileNorm.noise.txt 1e-05
```
