**Starting with actual values,**
- **log-transform, mean-center, and run hierarchical clustering on samples (columns) and on features (rows), all with Cluster 3.0**
- **draw dendrograms of samples and features, and draw a heatmap (with dendrograms)**

Required Linux programs: 'cluster' ([Cluster 3.0](http://bonsai.hgc.jp/~mdehoon/software/cluster/software))

Required R packages: pheatmap



Typical simplest command    **./cluster_draw_pheatmap.R matrix outputPDF**
```
./cluster_draw_pheatmap.R Expression_values.txt Cluster3_pheatmap.pdf
```

Change line 22 to 'display_numbers = TRUE' to get log2 fold changes printed on each cell
```
./cluster_draw_pheatmap.R Expression_values.txt Cluster3_pheatmap.larger_matrix.display_numbers.pdf
```

Include sample annotations 
```
./cluster_draw_pheatmap.R Expression_values.txt Cluster3_pheatmap.s.pdf Sample_annotations.txt
```

Include feature (gene) annotations 
```
./cluster_draw_pheatmap.R Expression_values.txt Cluster3_pheatmap.f.pdf NA Feature_annotations.txt
```

Include sample and feature (gene) annotations 
```
./cluster_draw_pheatmap.R Expression_values.txt Cluster3_pheatmap.sf.pdf Sample_annotations.txt Feature_annotations.txt
```

Further modifications by changing the actual code:
- Change line 22 to 'display_numbers = TRUE' to get log2 fold changes printed on each cell
- ./cluster_draw_pheatmap.R Expression_values.txt Cluster3_pheatmap.sf.display_numbers.pdf Sample_annotations.txt Feature_annotations.txt
- Change line 24 to 'breaks = seq(-3, 3, length.out = color.range.length+1)' to set the desired color range
- ./cluster_draw_pheatmap.R Expression_values.txt Cluster3_pheatmap.fixed_range.pdf Sample_annotations.txt Feature_annotations.txt
