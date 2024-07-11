#!/usr/bin/env python3

import sys
import argparse as ap
import pandas as pd

# Read a reference table of keys and values, use the values to update
# keys and values in query table.

# Version 1.0, 1 Mar., 2023
# Copyright, Troy Whitfield
# Bioinformatics and Research Computing
# Whitehead Institute

if __name__ == '__main__':
    
    parser = ap.ArgumentParser()
    parser.add_argument("-r","--ref", dest="REF", help="Input reference table.", default=ap.SUPPRESS)
    parser.add_argument("-q","--query", dest="QUERY", help="Input query table.", default=ap.SUPPRESS)
    parser.add_argument("-k","--geneKey", dest="GK", help="Key for gene IDs in reference table.", default=ap.SUPPRESS)
    parser.add_argument("-c","--convertIDs", action="store_true", dest="CV", default=False, help="When possible, replace IDs with gene symbols (default is false).")
    parser.add_argument("-e","--extraCols", action="store_true", dest="EC", default=False, help="Add extra referance columns to query table (default is false).")

    try:
        options = parser.parse_args()
        REF = pd.read_csv(options.REF, sep='\t')
        QUERY = pd.read_csv(options.QUERY, sep='\t')
        GK = options.GK
        REF.rename(columns = {list(REF)[0]:'UID'}, inplace=True)
        QUERY.rename(columns = {list(QUERY)[0]:'UID'}, inplace=True)
        
    except:
        print('id_add_name.py is intended to append gene names to an input matrix')
        print('using a given (input) mapping between the ID keys and gene names.')
        print('')
        print('It is assumed that the input matrix and mapping file are both tab')
        print('delimited.  Furthermore, the matrix IDs that will be mapped are')
        print('assumed to be in column 1.')
        print('')
        print('For the mapping file, the name of the column containing gene names')
        print('should be specified as input.')
        print('')
        print('Output is written as tab delimited text to the terminal, with') 
        print('gene names appearing in the last column of the output.  ID keys')
        print('with no matching gene name will have empty entries in the "Gene')
        print('name" column.')
        print('')
        print('The ordering and formatting of the output can be controled by the')
        print('-c flag.  Using this flag replaces IDs in column 1 of the input')
        print('and appends the original IDs as a different column.  Otherwise,')
        print('gene symbols will appear in a separate column, with IDs retained in')
        print('column 1 (default).')
        print('')
        print('Finally, if the mapping file contains columns in addition to ID')
        print('keys and gene names (e.g. descriptions of genes), these "extra"')
        print('annotations can optionally be included in the output by adding the')
        print('-e flag to the command.')
        print('')
        print("Example command: ./id_add_name.py -r ref_annot.txt -q expr.txt -k 'Gene name' -c -e")
        print('')
        sys.exit()       
    
    refName = REF.pop(GK); REF.insert(0,GK,refName)
    RG = REF[['UID',GK]]
    genes = dict(zip(RG['UID'],RG[GK]))
    if options.EC:
        QUERY = pd.merge(QUERY,REF,on=['UID'],how='left')
    else:
        QUERY = pd.merge(QUERY,RG,on=['UID'],how='left')
    if options.CV:
        QUERY[GK] = QUERY['UID']
        QUERY['UID'] = QUERY['UID'].map(genes).fillna(QUERY['UID'])    
        QUERY.rename(columns = {GK:'Gene.ID'}, inplace=True)
        QUERY.rename(columns = {'UID':GK}, inplace=True)
    else:
        QUERY.rename(columns = {'UID':'Gene.ID'}, inplace=True)
    QUERY.to_csv(sys.stdout, sep='\t', index=False)
