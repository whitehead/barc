# Convert standard NCBI BLAST output to a tab-delimited table of all fields

Requirements: BioPerl (Bio::SearchIO)

###  USAGE: ./parse_BLAST.pl blastFile > parsed_output

### Typical BLAST parsing to save all fields
```
./parse_BLAST.pl sample_BLAST_output.txt > sample_BLAST_output_parsed.txt
```

### Concise BLAST parsing to save selected fields
```
./parse_BLAST_to_m8.pl sample_BLAST_output.txt > sample_BLAST_output.m8.txt
```

The input BLAST report must be in text (not HTML).
	
Each HSP (alignment; "High-scoring Segment Pair") is represented
as one line of output.
Since a query-hit pair of sequences can have multiple HSPs,
the output for a single hit can be spread over several lines.
	
The output file is tab-delimited text, which can be opened in Excel.

Because Excel tries too hard to be smart, some cells in the "hspConsensus"
column may appear as "#NAME?".  If that happens, you can try formatting 
that column as text and repasting in that column.
