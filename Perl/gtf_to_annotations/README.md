# Extract feature annotations from the GTF (or GTF.gz) file into a tab-delimited text table

```
USAGE: ./gtf_to_annotations.pl features.gtf [featureType] > features.txt

### Extract all features
./gtf_to_annotations.pl Features.gtf.gz > Features.gtf.txt

### Extract only CDS features 
./gtf_to_annotations.pl Features.gtf.gz CDS > CDS_features.gtf.txt
```
