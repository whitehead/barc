## Draw a volcano plot for differential gene expression, highlighting an enriched gene-set/pathway.

`drawGSEAvolcano.R` can be used to draw a volcano plot that highlights genes
belonging to specific pathway/gene-set.  It is assumed that the estimated
effect sizes (reported as base 2 log fold-changes) and p-values drawn in
the volcano plot are reported by DESeq2 (or edgeR etc., as long as the
file format from the example is adhered to).  Likewise, users are expected
to have run GSEA already, so that genes within a specific gene-set are
known, along with whether they belong to the 'leading edge' of the enrichment.

All genes from within the GSEA-tested pathway are coloured in red, as compared
with genes not belonging to the pathway, which are coloured grey.  Pathway
highlighted genes that belong to the 'leading edge' and exhibit a greater than
two-fold effect are indicated by name in the plot.

Example command: ./drawGSEAvolcano.R DESeq2.csv GSEA.tsv MSigDB.txt fold-change GSEAvolcano.pdf

### Required packages
ggplot (3.5.1), ggrepel (0.9.4), Cairo (1.6-2)
