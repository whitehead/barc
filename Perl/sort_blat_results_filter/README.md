# Sort a BLAT file by choosing the hit with the best score for each query

The score is defined as "matches - 3 * mismatches".

In the case of >1 hits with the same score, all are selected.

If filter is set as 'T', aligned bases must be at least 50% of length of query.  If they aren't, they'll be removed and appear in the "*.poor.blat" output file.

```
# USAGE:  ./sort_blat_results_filter.pl BLAT_file sorted_BLAT_results Filter[T/F]
./sort_blat_results_filter.pl BLAT_output.psl BLAT_output.sorted_filtered.psl T
```
