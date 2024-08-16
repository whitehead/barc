
# Join two or more files using matching IDs in the first column of each file.

### Usage: 
    
    ./joinFilesByFirstColumn.pl input1 input2 input3 > output


### Note:

Input files must be tab-delimited.

The first column in each file should contain the IDs to be matched.

IDs are case-sensitive and should be unique within each file.

The order of rows in the merged file will match the order of IDs from the first input file.

The order of columns in the merged file will follow the sequence specified in the command line.

If an ID from the first file is not found in a file, the corresponding record will have 'UNMATCHED' in place of missing data.


