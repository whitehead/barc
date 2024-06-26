# Correct batch effects with sva's ComBat() and limma's removeBatchEffect()

Required R packages: sva, limma, made4

Input files:
- 1 - Log2-transformed matrix with row IDs in first column
- 2 - Sample info file: rows are identical to columns of data matrix file
   -   One column is named "batch" 
   -   One column can be names "covariate" [optional]
   -   Other columns are ignored

Output files (2) will be batch-corrected log2 matrix files.


### No covariates (sample differences are only due to batch)
```
./remove_batch_effects.R Expression.log2_values.txt Sample_information.no_covariate.txt MyBatchAdjustment
```

### Covariates (sample differences are due to batch and covariate)
```
./remove_batch_effects.R Expression.log2_values.txt Sample_information.txt MyBatchCovariateAdjustment
```

Caveat: See "Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses"
[PubMed: 26272994](https://www.ncbi.nlm.nih.gov/pubmed/26272994)
