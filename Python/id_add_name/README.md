## Append gene symbols to an matrix of values with ENSEMBL ID row names

`id_add_name.py` is intended to append names to an input matrix
using a given (input) mapping between the ID keys and names. The
example here maps ENSEMBL IDs to gene symbols.

It is assumed that the input matrix and mapping file are both tab
delimited.  Furthermore, the matrix IDs that will be mapped are
assumed to be in column 1.

For the mapping file, the name of the column containing gene names
should be specified as input.

Output is written as tab delimited text to the terminal, with
gene names appearing in the last column of the output.  ID keys
with no matching gene name will have empty entries in the "Gene
name" column.

The ordering and formatting of the output can be controled by the
-c flag.  Using this flag replaces IDs in column 1 of the input
and appends the original IDs as a different column.  Otherwise,
gene symbols will appear in a separate column, with IDs retained in
column 1 (default).

Finally, if the mapping file contains columns in addition to ID
keys and gene names (e.g. descriptions of genes), these "extra"
annotations can optionally be included in the output by adding the
-e flag to the command.

Example command: `./id_add_name.py -r ref_annot.txt -q expr.txt -k 'Gene name' -c -e`

### Required packages
pandas, argparse
