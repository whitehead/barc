#!/usr/bin/python3
"""Select the top ranked name based on order of columns
Copyright (c) 2020
@author: Bingbing Yuan
"""

import argparse
import textwrap


def parse_arguments():
    parser = argparse.ArgumentParser(
        usage="\nThis program takes a tab delimited file, such as genome-wise gene correlation table,\n" +
        "sorts the value for each row, and \n" + 
        "picks up the most associated/anti-associated items (such as genes) with extreme values.\n",
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent("""
        Sample commands: 
        %(prog)s input_tab.txt 5 top > top_5_output.txt
        %(prog)s input_tab.txt 5 bottom >| bottom_5_output.txt
        """) + '\n\n'
    )

    parser.add_argument('infile', type=str, help='a tab delimited input file, such as correlation table')
    parser.add_argument('number', type=int, help="number of items (such as genes) with extreme values")
    parser.add_argument('type', metavar='type', type=str.lower, help='most associate (top) or anti-associated (bottom) items',
                        choices=['top', 'bottom'])

    return parser.parse_args()


def select_items(file, number, type):
    f = open(file, 'r')
    genes = ()
    row_num = 0

    for line in f:
        row_num = row_num + 1
        if row_num == 1:
            #       A1BG  A1CF  A2M
            genes = line.rstrip().split("\t")
            # skip the first one
            genes = genes[1:]
            # print (genes)

        else:
            # A1BG  1.0 -0.07  0.08
            cor = line.rstrip().split("\t")
            # first item is targeted gene
            target_gene = cor[0]
            cor = cor[1:]
            # convert string to float
            cor = [float(idx) for idx in cor]

            # create a dictionary
            # one coefficient value could potentially associates with multiple genes
            # re-create dictionary whose values are list of genes

            res = {}
            for i in range(len(cor)):
                if cor[i] not in res.keys():
                    res[cor[i]] = [genes[i]]
                else:
                    res[cor[i]].append(genes[i])
                    # print (res[cor[i]])

            # print(res)

            # print results
            if type.upper() == "TOP":
                # up items
                n = 0
                for i in reversed(sorted(res.keys())):
                    n = n + 1
                    if n > number:
                        break
                    print(target_gene + "\t" + str(i) + "\t" + ','.join(res[i]))

            elif type.upper() == "BOTTOM":

                # bottom items
                n = 0
                for i in sorted(res.keys()):
                    n = n + 1
                    if n > number:
                        break
                    print(target_gene + "\t" + str(i) + "\t" + ','.join(res[i]))


if __name__ == '__main__':
    args = parse_arguments()
    select_items(args.infile, args.number, args.type)
