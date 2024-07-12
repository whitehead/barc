# Perform differential expression analysis on microarray experiments or other types of continuous-value data.  


### Inputs:

* Matrix file: first column contains feature names/IDs
   * first row contains sample names/IDs (and can't contain special characters, according to R)
   * Values are expected NOT to be log-transformed; the code will do that (which is required before limmma analysis).
 * Design file: first column contains sample names/IDs
   * second column lists groups; any other columns are ignored

### Note that this code does not include a normalization step, which is a typical prerequisite for limma analysis. We have separate code to do that: R/normalize_matrix/normalize_matrix.R

### Run a simple 2-group comparison

Design file: first column contains sample names/IDs; second column lists groups; any other columns are ignored

```
./Run_2_group_limma_differential_expression.R Experiment_B.values.txt Experiment_B.design.txt 0.1 Experiment_B.limma_out.txt
```

### Simple 2-group comparison with batches or a paired design

Design file: first column contains sample names/IDs; second column lists groups; third column lists batches

```
./Run_2_groups_with_batches_limma_differential_expression.R Experiment_C.values.txt Experiment_C.design.txt 0.1 Experiment_C.limma_out.txt
```

### Run an all vs. all comparison (for more than 2 groups)

```
# Arguments: matrixFile designFile pseudocounts outputFile
./Run_all_vs_all_limma_differential_expression.R Experiment_A.values.txt Experiment_A.design.txt 0.1 Experiment_A.limma_all_vs_all.txt
```

### Split all vs. all comparison limma output file into single-comparison files
```
./split_limma_output_by_comparison.R Experiment_A.limma_all_vs_all.txt
```

Output files are in limma_out_by_comparison

###  Given a limma output file, create a MA plot

```
# USAGE: ./draw_MA_plot_from_limma_analysis.R inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold
./draw_MA_plot_from_limma_analysis.R Experiment_D.limma_out.txt Sample.1 Sample.2 1e-05 1
```

###  Given a limma output file, create a volcano plot

```
# USAGE: ./draw_volcano_plot_from_limma_analysis.R inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold
./draw_volcano_plot_from_limma_analysis.R Experiment_D.limma_out.txt Sample.1 Sample.2 1e-05 1
```

### Given a limma output file, create a MA plot but label only genes in input list

```
./draw_MA_plot_from_limma_analysis.label_selected.R Experiment_D.limma_out.txt Sample.1 Sample.2 1e-05 1 Genes_to_label.txt
```

### Given a limma output file, create a volcano plot but label only genes in input list

```
./draw_volcano_plot_from_limma_analysis.label_selected.R Experiment_D.limma_out.txt Sample.1 Sample.2 1e-05 1 Genes_to_label.txt
```