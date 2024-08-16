
# Merge two files allowing duplicated keys

### Usage:
    ./merge_two_files_allowing_dupl_key.pl file.txt anno_file 1st_file_id_col out_col field_seperator

### Example:
    ./merge_two_files_allowing_dupl_key.pl input.tsv anno.tsv 3 4 \, > output.tsv

### Note:
Assume both input text files are tab-delimited.

The first column in the 2nd file must be unique and will be used as the key (ID) to match with the ID in the 1st file, specified as the 3rd argument

The corresponding value from the user-defined output column of the 2nd file, the 4th argument,  will be appended to the 1st file.

If a record in the 1st file contains multiple IDs, separated by a user-defined delimiter (specified as the last argument), each ID will be used for matching against the 2nd file.

The resulting merged file will maintain the same order of IDs as in the 1st file.


