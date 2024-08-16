#!/usr/bin/python3
"""Add pseudo count to a table
Copyright (c) 2020
@author: Bingbing Yuan
"""


import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        usage='This script is for adding a value to a table.',
        formatter_class=argparse.RawTextHelpFormatter,
        description=("""
        Sample command:
        %(prog)s input.txt value > output.txt
        """)
    )

    parser.add_argument('infile', type=str, help='a tab delimited input file')
    parser.add_argument('pseudoCount', type=float, help="a number, could be an integer or float")

    return parser.parse_args()


def add_pseudo_count(file, number):
    content = pd.read_csv(file, sep="\t", index_col=0)
    print(content + number)


if __name__ == '__main__':
    args = parse_arguments()
    add_pseudo_count(args.infile, args.pseudoCount)

