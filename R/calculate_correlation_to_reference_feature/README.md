# Calculate correlation for any feature to a reference feature

The inputCounts table is log2-transformed (after adding a pseudocount) and median-centered before correlations are calculated.

Reference feature can be anywhere in file

### USAGE:
./calculate_correlation_to_reference_feature.R inputCounts referenceID OutputFile
```
./calculate_correlation_to_reference_feature.R MySample.counts.txt ENSG00000223972.4 MySample.gene_correlations.txt
```

Run a variant of the code that skips all pre-processing of the matrix (log-transformation and median-centering) and uses the input values for the correlation calculation
```
./calculate_correlation_to_reference_feature.as_is.R MySample.counts.txt ENSG00000223972.4 MySample.gene_correlations_as_is.txt
```
