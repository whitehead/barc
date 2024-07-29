# filer and summarize gene set enrichment results from either fgsea or GSEA 



### Using fgsea output file


Plot top and bottom 10 terms (and truncate any term that's longer than 50 characters) with FDR less than 0.05
but make points with a cex (circle size) limit that's 10 greater than the default (6) 

```
./draw_enrichment_figures.R Res_DE1_7_Hallmarksigpath.txt fgsea 0.05 10 50 10 
```

Plot top and bottom 20 terms (and truncate any term that's longer than 100 characters)

```
./draw_enrichment_figures.R Res_DE1_7_GOsigpath.txt fgsea 1 20 100
```

Plot top and bottom 30 terms (and truncate any term that's longer than 100 characters)
but make points with a cex limit that's less than the default (6)

```
./draw_enrichment_figures.R Res_DE1_7_GOsigpath.txt fgsea 1 30 100 3
```


### Using GSEA output file


Same output with either 'pos' or 'neg' file, but both 'neg' and 'pos' files are read/required
Plot top and bottom 20 terms (and truncate any term that's longer than 100 characters) with FDR less than 0.01
Note: The number of terms displayed is also depended on the number of terms below the FDR cutoff.
In this example, only 13 instead of 20 terms are plotted for positive NES values since only 13 terms are with FDR less than 0.01 


```
./draw_enrichment_figures.R gsea_report_for_na_neg_1586451599999.xls GSEA 0.01 20 100
```

The above command creates the same figure as the one below.

```
./draw_enrichment_figures.R gsea_report_for_na_pos_1586451599999.xls GSEA 0.01 20 100
```
The 'gene set size' displayed on the figure is the number of genes in the gene set found in the input gene+rank file.

Increase maximum point size to 20 (on ggplot scale)

```
./draw_enrichment_figures.R gsea_report_for_na_pos_1586451599999.xls GSEA 0.01 20 100 20
```
