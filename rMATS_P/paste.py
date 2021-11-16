#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse


USAGE = '''usage: %(prog)s [options] arg1 arg2'''


def get_args():
    """TODO: Docstring for get_args.
    :returns: TODO

    """
    parser = argparse.ArgumentParser(usage=USAGE)

    parser.add_argument('-i', action='store', default='', required=True,
                        help='input filename.', dest='i')
    parser.add_argument('--o1', action='store', default='', required=True,
                        help='output1 filename.', dest='o1')
    parser.add_argument('--o2', action='store', default='', required=True,
                        help='output2 filename.', dest='o2')

    args = parser.parse_args()

    return args


def main():
    """TODO: Docstring for main.
    :returns: TODO

    """
    args = get_args()

    with open(args.i, 'r') as ifp:
        with open(args.o1, 'w') as ofp1:
            with open(args.o2, 'w') as ofp2:
                for line in ifp:
                    columns = line.strip().split('\t')
                    ofp1.write(columns[0]+'\n')
                    ofp2.write('\t'.join(columns[1:])+'\n')

    return


if __name__ == "__main__":
    main()
