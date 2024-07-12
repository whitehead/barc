# Select rows and/or columns from a matrix using lists of row and/or column names


```
# USAGE: /home/gbell/BaRC_code/submatrix_selector.pl matrixFile rowIDfile columnIDfile > submatrixFile[output]

# To select only by row IDs:
./submatrix_selector.pl matrixFile rowIDfile > submatrixFile[output]

# To select only by columns IDs:
./submatrix_selector.pl matrixFile 0 columnIDfile > submatrixFile[output]

./submatrix_selector.pl matrix.txt Rows_to_extract.txt Columns_to_extract.txt > matrix_selected.txt
```