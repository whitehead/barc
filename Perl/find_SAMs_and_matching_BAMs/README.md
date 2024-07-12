# Recursively search a desired directory for SAM files and matching BAM files, pairing Foo.sam with Foo.bam or Foo.sorted.bam

This method is designed to reduce redundancy in file storage and archiving.  If a SAM file alone is identified, it is recommended to convert it to BAM and delete the SAM.

### Search in this directory
```
./find_SAMs_and_matching_BAMs.pl . >| SAM_BAM.fileList+commands.thisDir.txt
```

### Look for SAM+BAM files up one directory
```
./find_SAMs_and_matching_BAMs.pl .. >| SAM_BAM.fileList+commands.txt
```
