## Draw a Sankey diagram

Sankey diagrams can illustrate flux in dynamic processes.  In this application
a Sankey diagram is generated to depict the relationships between two sets
of class labels.

A good way to make Sankey diagrams with R is to use the riverplot
package, which requires nodes and edges as input.  The wrapper function in this
directory generates nodes and edges from an input matrix.  This input matrix is
expected to be a 3-column data frame.  The makeRiver function is used
below to convert nodes and edges into a riverplot object before plotting.

This is the expected format for matrix input:

Column 1: object identifier.
Column 2: class labels 1.
Column 3: class labels 2.

Note that missing class labels will be converted to NA and will not
contribute to edge widths.  To capture the birth and death of objects,
define and label the appropriate null states in the input file (i.e. fill in
the missing values).

Finally, the colours for the sources and sinks in the Sankey diagram with
two sets of class labels are given by the input arguments col1 and col2,
which can be assigned to valid R colour names or RGB codes.

### Required packages
riverplot (0.10), Cairo (1.6-2)
