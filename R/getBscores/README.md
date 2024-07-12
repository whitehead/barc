## Compute B-scores for high-throughput plate-based assays

This script allows easy command line access to computing b-scores
(see C. Brideau et al. "Improved Statistical Methods for Hit Selection
in High-Throughput Screening", _Journal of Biomolecular Screening_ **8**,
634-647 (2003)) from raw single plate matrices for high throughput
screens.

Convert a raw input matrix (e.g. from high throughput
screening experiments) into a matrix of B-scores, thereby
correcting for plate effects.

USAGE: `./getBscores.R input scale output`

Example: `./getBscores.R input.csv additive output.csv`

The input and output files are each a matrix in comma
separated variable format.

The 'scale' keyword can be either 'additive' (default)
or 'multiplicative'.  Multiplicative data (e.g. ratios)
will be (base 2) logaritmically transformed before b-scores
are computed.
