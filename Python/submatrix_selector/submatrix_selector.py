#! /usr/bin/env python3

#######################################################################
#
# Copyright(c) 2022 Whitehead Institute for Biomedical Research.
#			  All Rights Reserved
#
# Given a tab-delimited file matrix, select desired rows and columns.
# Adapted from submatrix_selector.pl
# Version 1.0: 7 March 2022 
#
# George Bell, BaRC
#
#######################################################################

import sys
import csv

# Ensure that at least two arguments were entered
try:
	matrixFile = sys.argv[1]
	rowFile = sys.argv[2]
except IndexError:
	sys.stderr.write("\nSubmatrix selector (select rows and/or columns from a matrix)\n")
	sys.stderr.write("USAGE: submatrix_selector.py matrixFile rowIDfile columnIDfile > submatrixFile[output]\n\n")
	sys.stderr.write("To select only by row IDs:\n")
	sys.stderr.write("       submatrix_selector.py matrixFile rowIDfile > submatrixFile[output]\n\n")
	sys.stderr.write("To select only by columns IDs:\n")
	sys.stderr.write("       submatrix_selector.py matrixFile 0 columnIDfile > submatrixFile[output]\n\n")
	exit()
# Get this (optional) argument
try:
	colFile = sys.argv[3]
except IndexError:
	sys.stderr.write("No column file given.  Printing all columns ....\n")
	colFile = "0"

if (rowFile == "0"):
	sys.stderr.write("No row file given.  Printing all rows ....\n")

# Initialize lists of rows and columns we want
rowsToGet = []
colsToGet = []

# sys.stderr.write("Col file: " + colFile)

# Read the list of row IDs to print
if (rowFile and rowFile != "0"):
	with open(rowFile) as rows:
		reader = csv.reader(rows, delimiter="\t")
		for row in reader:
			# print(row)
			rowsToGet.append(row[0])
row_dictionary = dict.fromkeys(rowsToGet, "get")

# Read the list of column IDs to print
if (colFile and colFile != "0"):
	# sys.stderr.write("Reading " + colFile + " to get columns to print ....\n")
	with open(colFile) as cols:
		reader = csv.reader(cols, delimiter="\t")
		for col in reader:
			colsToGet.append(col[0])
col_dictionary = dict.fromkeys(colsToGet, "get")

# Read through big matrix file
with open(matrixFile) as matrixData:
	reader = csv.reader(matrixData, delimiter="\t")
	for matrixRow in reader:
		fieldsToPrint = []
		if (reader.line_num == 1):
			matrixHeaderRow = matrixRow
			for item in matrixHeaderRow:
				if item in col_dictionary.keys() or colFile == "0":
					fieldsToPrint.append(item)
			print("\t".join(fieldsToPrint))
		else:
			# Iterate through this row (and header row, again)
			for headerItem, item in zip(matrixHeaderRow, matrixRow):
				if headerItem in col_dictionary.keys() or colFile == "0":
					fieldsToPrint.append(item)
			# Is the first field of this row in our "rows to get"?
			if matrixRow[0] in row_dictionary.keys() or rowFile == "0":
				print("\t".join(fieldsToPrint))

# End
