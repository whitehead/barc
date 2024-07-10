# Use RunDESeq2.R to get differentially expressed genes and normalize counts

```
./RunDESeq2.R Input_DESeq2.txt DESeq2_output.txt UHR UHR brain brain
```

### Run DESeq2 to normalize counts and make PCA and other plots.  Colors file should list the same color for each member of a group
```
./DESeq2_normalize_only.R Input_DESeq2.txt DESeq2_output.norm_only.txt rlog sampleGroupColors.txt
./DESeq2_normalize_only.R Input_DESeq2.txt DESeq2_output.norm_only.txt vst
```

### Calculate FPKM values from raw counts and gene lengths
### To get "gene length", one way is to use Perl/gtf_to_transcript_lengths/gtf_to_transcript_lengths.pl GENE_MODELS.gtf > GENE_MODELS.transcript_lengths.txt
### which also makes files of mean and median gene lengths. 
```
./counts_to_FPKM_TPM_using_DESeq2.R Mouse_RNAseq_counts.txt Mouse.medianGeneLength.txt Mouse_RNAseq.FPKM.TPM.txt
```

### Draw plots, using a previously performed DESeq2 analysis
```
./draw_MA_plot_from_DESeq2_analysis.R DESeq2_output.txt brain UHR 1e-10 5
./draw_volcano_plot_from_DESeq2_analysis.R DESeq2_output.txt brain UHR 1e-10 5
./draw_volcano_plot_from_DESeq2_analysis.label_selected.R DESeq2_output.txt brain UHR 1e-10 5 Genes_to_label.txt
```

### Draw a plot and label user-selected genes
```
# USAGE: ./draw_MA_plot_from_DESeq2_analysis.label_selected.R inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold selectedGenesFile
draw_MA_plot_from_DESeq2_analysis.label_selected.R DESeq2_output.txt brain UHR 1e-10 5 Genes_to_label.txt
```

### Do a differential expression analysis comparing each group to each other group
```
# USAGE: ./Run_all_vs_all_DESeq2_differential_expression.R inputMatrix designFile outputDirectory
./Run_all_vs_all_DESeq2_differential_expression.R Three_groups.counts.txt Three_groups.design.txt Three_groups.counts.DESeq2_comparisons
```

### Draw a volcano plot, labeling user-selected genes
./draw_volcano_plot_from_DESeq2_analysis.label_selected.R Volcano_plot_input.txt brain UHR 1e-10 5 Volcano_plot_genes.txt

### Draw a MA plot, labeling user-selected genes
```
./draw_MA_plot_from_DESeq2_analysis.label_selected.R Volcano_plot_input.txt brain UHR 1e-10 5 Volcano_plot_genes.txt
```

### Use RunDESeq2.R to get differentially expressed genes and normalize counts, but in the process, try different lfcShink methods ()
```
./RunDESeq2_compare_shrinkage_methods.R Input_DESeq2.txt DESeq2_lfcShrink_comparison.out.txt UHR UHR brain brain
```

### Normalize a matrix of raw counts with DESeq2 BUT don't adjust by size factors (with optional PCA plot group colors).
```
./DESeq2_normalize_only.no_size_factors.R Input_DESeq2.txt DESeq2_output.norm_only.txt rlog sampleGroupColors.txt
```

### Use the design file (sampleGroupColors.txt) to create a DESeq object with a design
```
./DESeq2_normalize_only.with_design.R Input_DESeq2.d.txt DESeq2_output.d.norm_only.txt rlog sampleGroupColors.txt
```

### Run DESeq2 to normalize counts and make PCA and other plots *BUT* use custom size factors
```
RunDESeq2.custom_size_factors.R Input_DESeq2.txt DESeq2_output.customSF.txt UHR UHR brain brain customSizeFactors.txt
```
