# Use RunDESeq2.R to get differentially expressed genes and normalize counts

Required: DESeq2, RColorBrewer, genefilter, ggplot2, ggrepel, lattice, limma, reshape2, vsn

### Basic command to do 2-way comparison:
```
# USAGE: ./RunDESeq2.R inputCounts OutputFile group1 group2 [...]
./RunDESeq2.R Input_DESeq2.txt DESeq2_output.txt UHR UHR brain brain
```

### Run DESeq2 to normalize counts and make PCA and other plots.  Colors file should list the same color for each member of a group
```
# USAGE: ./DESeq2_normalize_only.R inputCounts OutputFile pcaLogMethod[log2|rlog|vst] [colorFile]
./DESeq2_normalize_only.R Input_DESeq2.txt DESeq2_output.norm_only.txt rlog sampleGroupColors.txt
./DESeq2_normalize_only.R Input_DESeq2.txt DESeq2_output.norm_only.txt vst
```

### Calculate FPKM values from raw counts and gene lengths
### To get "gene length", one way is to use Perl/gtf_to_transcript_lengths/gtf_to_transcript_lengths.pl GENE_MODELS.gtf > GENE_MODELS.transcript_lengths.txt
### which also makes files of mean and median gene lengths. 
```
# USAGE: ./counts_to_FPKM_TPM_using_DESeq2.R inputCounts geneLengths OutputFile
./counts_to_FPKM_TPM_using_DESeq2.R Mouse_RNAseq_counts.txt Mouse.medianGeneLength.txt Mouse_RNAseq.FPKM.TPM.txt
```

### Draw plots, using a previously performed DESeq2 analysis
```
# USAGE: ./draw_MA_plot_from_DESeq2_analysis.R inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold
./draw_MA_plot_from_DESeq2_analysis.R DESeq2_output.txt brain UHR 1e-10 5
# USAGE: ./draw_volcano_plot_from_DESeq2_analysis.R inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold
./draw_volcano_plot_from_DESeq2_analysis.R DESeq2_output.txt brain UHR 1e-10 5
# USAGE: ./draw_volcano_plot_from_DESeq2_analysis.label_selected.R inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold genesFile
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

### Use RunDESeq2.R to get differentially expressed genes and normalize counts, but in the process, try different lfcShink methods ()
```
# USAGE: ./RunDESeq2_compare_shrinkage_methods.R inputCounts OutputFile group1 group2 [...]
./RunDESeq2_compare_shrinkage_methods.R Input_DESeq2.txt DESeq2_lfcShrink_comparison.out.txt UHR UHR brain brain
```

### Normalize a matrix of raw counts with DESeq2 BUT don't adjust by size factors (with optional PCA plot group colors).
```
# USAGE: ./DESeq2_normalize_only.no_size_factors.R inputCounts OutputFile pcaLogMethod[log2|rlog|vst] [colorFile]
./DESeq2_normalize_only.no_size_factors.R Input_DESeq2.txt DESeq2_output.norm_only.txt rlog sampleGroupColors.txt
```

### Use the design file (sampleGroupColors.txt) to create a DESeq object with a design (from the colors file)
```
# USAGE:  ./DESeq2_normalize_only.with_design.R inputCounts OutputFile pcaLogMethod[log2|rlog|vst] [colorFile]
./DESeq2_normalize_only.with_design.R Input_DESeq2.d.txt DESeq2_output.d.norm_only.txt rlog sampleGroupColors.txt
```

### Run DESeq2 to normalize counts and make PCA and other plots *BUT* use custom size factors
```
# USAGE: ./RunDESeq2.custom_size_factors.R inputCounts OutputFile group1 group2 [...] sizeFactorsFile
RunDESeq2.custom_size_factors.R Input_DESeq2.txt DESeq2_output.customSF.txt UHR UHR brain brain customSizeFactors.txt
```
