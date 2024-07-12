# This script takes a matrix of values, adds pseudocounts to all values, log2-transforms and subtracts the mean or median for each row
## Additionally it can add an extra row ID column to view output in Java TreeView
### Input file format: tab delimited first row is header and first column has row IDs.
### Run mean-centering
```
./log2_mean_median_center.R  MyValues.txt mean 0.001 > MyValues.mean_centered.txt
```

### Run median-centering
```
./log2_mean_median_center.R  MyValues.txt median 0.001 > MyValues.median_centered.txt
```


### Run median-centering and add extra column of row IDs to make output ready for Java TreeView viewing
```
./log2_mean_median_center.R  MyValues.txt median 0.001  jtv  > MyValues.median_centered_forJTV.txt
```