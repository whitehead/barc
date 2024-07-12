# calculate Fraction of Reads in Peaks (FRiP)
Note that a BED or a 'narrowPeak' file (as created by MACS) can be used to identify peaks.

### Typical usage:
```
./calculate_FRiP_score.py SampleA.bam SampleA_peaks.narrowPeak
./calculate_FRiP_score.py SampleA.bam SampleA_peaks.bed
```


### To include multimapped reads with the counted reads:
```
./calculate_FRiP_score.py SampleA.bam SampleA_peaks.narrowPeak -M
```
