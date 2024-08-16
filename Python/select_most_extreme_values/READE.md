# select most extreme values

### This program takes a tab delimited file, such as genome-wise gene correlation table, sorts the value for each row, and picks up the most associated/anti-associated items (such as genes) with extreme values.

Usage:

    select_most_extreme_values.py correlation_table.txt number top|bottom

Examples:

    ./select_most_extreme_values.py input_tab.txt 5  top >| top_5_output.txt
    
    ./select_most_extreme_values.py input_tab.txt 5  bottom >| bottom_5_output.txt

