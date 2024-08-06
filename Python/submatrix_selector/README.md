# Select rows and/or columns from a matrix using lists of row and/or column names


```
# USAGE: submatrix_selector.py matrixFile rowIDfile columnIDfile > submatrixFile[output]

# To select only by row IDs:
./submatrix_selector.py matrix.txt Rows_to_extract.txt >| matrix_selected_rows.txt

# To select only by column IDs:
./submatrix_selector.py matrix.txt 0 Columns_to_extract.txt >| matrix_selected_columns.txt

# To select by row and column IDs:
./submatrix_selector.py matrix.txt Rows_to_extract.txt Columns_to_extract.txt >| matrix_selected.txt
```