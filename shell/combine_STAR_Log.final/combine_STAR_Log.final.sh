#!/bin/bash

###
###  Combine *Log.final.out files from multple STAR runs.
###  Bingbing Yuan and George Bell, Whitehead BaRC
###  Version 1.0: 28 March 2018
###

if [ "$#" -ne 1 ]
then
echo "
Merge *Log.final.out files from multple STAR runs.
USAGE: ./combine_STAR_Log.final.sh <STAR output directory> > Merged_STAR_Log.final.txt
"

exit
fi

# Get argument (directory with STAR output)
directory=$1

# Combine content:
paste $directory/*Log.final.out | awk 'BEGIN{FS="\t"; OFS="\t"} {out=""; for(i=2;i<=NF;i=i+2){out=out"\t"$i}; print $1"\t"out}' | tail --lines=+4 | cut -f1,3- | perl -pe 's/^ +//' >| $directory.stat.table.tmp

# Get column names (merged filenames):
/bin/ls $directory/*Log.final.out | perl -pe 's/_Log.final.out//g' | perl -pe 's/\n/\t/g' | perl -pe 's/\t$/\n/' | awk 'BEGIN{FS="\t"; OFS="\t"} { if (NR==1) $1="Statistic\t"$1; print $0 }' > $directory.header.tmp

# Add colnames on top of content (and print to STDOUT)
cat $directory.header.tmp $directory.stat.table.tmp 

# Remove temporary files
rm -f $directory.header.tmp $directory.stat.table.tmp 
