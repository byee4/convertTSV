#!/usr/bin/env python

"""
basic module to use as an example.
"""
import argparse
import pandas as pd
import collections
import tables
import scipy.io
import scipy.sparse
import csv
import h5py

def get_matrix_from_h5(filename, genome):
    """
    Returns matrix and attributes using cellranger's h5 accessor method.

    :param filename:
    :param genome:
    :return:
    """
    GeneBCMatrix = collections.namedtuple('GeneBCMatrix',
                                          ['gene_ids', 'gene_names',
                                           'barcodes', 'matrix'])

    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print "That genome does not exist in this file."
            return None
        gene_ids = getattr(group, 'genes').read()
        gene_names = getattr(group, 'gene_names').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)
        return GeneBCMatrix(gene_ids, gene_names, barcodes, matrix)


def read_mtx_as_dataframe(mtx_file, columns_file, rows_file):
    """
    Reads a mtx file and returns a pandas dataframe
    :param mtx_file: sparse matrix
    :return df: Pandas.DataFrame()
    """
    mat = scipy.io.mmread(mtx_file)
    columns = [
        row[0] for row in csv.reader(open(columns_file), delimiter="\t")
    ]
    rows = [
        row[0] for row in csv.reader(open(rows_file), delimiter="\t")
    ]
    df = pd.DataFrame(mat.todense(), columns=columns, index=rows)
    return df


def read_sv_as_dataframe(fn, sep='\t'):
    """
    Assumes tabbed file includes rownames (first column) and colnames (first row).

    :param fn: string
        filename
    :return:
    """
    df = pd.read_table(fn, index_col=0, sep=sep)
    return df, df.columns, df.index


def write_dataframe_as_matrix(df, output_file):
    """
    Writes the pandas dataframe as a matrix.

    :param df:
    :param output_file:
    :return:
    """
    mtx = scipy.sparse.csr_matrix(df.values)
    scipy.io.mmwrite(output_file, mtx)


def write_dataframe_as_sv(df, output_file, sep='\t'):
    """
    Writes a {sep}-separated-values (sv) dataframe to file

    :param df:
    :param output_file:
    :param sep:
    :return:
    """
    df.to_csv(output_file, sep=sep)


def convert(input_file, input_type, output_file, output_type,
            columns_file, rows_file, genome):
    """
    Runs the conversion between input_type to output_type

    :param input_file:
    :param input_type:
    :param output_file:
    :param output_type:
    :param columns_file:
    :param rows_file:
    :return:
    """
    if input_type == 'csv':
        if output_type == 'mtx':
            csv2mtx(input_file, output_file)
    elif input_type == 'tsv':
        if output_type == 'mtx':
            tsv2mtx(input_file, output_file)
    elif input_type == 'mtx':
        if output_type == 'tsv':
            mtx2tsv(input_file, rows_file, columns_file, output_file)
        elif output_type == 'csv':
            mtx2csv(input_file, rows_file, columns_file, output_file)
    elif input_type == 'h5':
        if output_type == 'tsv':
            h52tsv(input_file, output_file, genome)

def csv2mtx(input_file, output_file):
    """
    Converts a comma separated file into an mtx file

    :param input_file:
    :param output_file:
    :return:
    """
    df, cols, rows = read_sv_as_dataframe(input_file, ',')
    write_dataframe_as_matrix(df, output_file)

    o = open(output_file + '.columns', 'w')
    for column in cols:
        o.write(column + '\n')
    o.close()

    o = open(output_file + '.rows', 'w')
    for row in rows:
        o.write(row + '\n')
    o.close()


def tsv2mtx(input_file, output_file):
    """
    Converts a tab separated file into an mtx file

    :param input_file:
    :param output_file:
    :return:
    """
    df, cols, rows = read_sv_as_dataframe(input_file, '\t')
    write_dataframe_as_matrix(df, output_file)

    o = open(output_file + '.columns', 'w')
    for column in cols:
        o.write(column + '\n')
    o.close()

    o = open(output_file + '.rows', 'w')
    for row in rows:
        o.write(row + '\n')
    o.close()


def mtx2tsv(input_file, rows_file, columns_file, output_file):
    """
    Converts mtx file + rows + columns into a tab-separated file.

    :param input_file:
    :param rows_file:
    :param columns_file:
    :param output_file:
    :return:
    """
    df = read_mtx_as_dataframe(
        input_file, columns_file, rows_file
    )
    write_dataframe_as_sv(df, output_file, sep="\t")


def mtx2csv(input_file, rows_file, columns_file, output_file):
    """
    Writes a csv from a mtx file

    :param input_file:
    :param rows_file:
    :param columns_file:
    :param output_file:
    :return:
    """
    df = read_mtx_as_dataframe(
        input_file, columns_file, rows_file
    )
    write_dataframe_as_sv(df, output_file, sep=",")


def tsv2h5(input_file, output_file, genome):
    pass


def h52tsv(input_file, output_file, genome):
    """
    Converts an h5 input file into a tab separated file.

    :param input_file:
    :param output_file:
    :param genome:
    :return:
    """
    gene_ids, gene_names, barcodes, mat = get_matrix_from_h5(input_file, genome)
    df = pd.DataFrame(mat.todense(), columns=barcodes, index=gene_ids)
    write_dataframe_as_sv(df, output_file, '\t')

def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        required=True,
    )
    parser.add_argument(
        "--input_type",
        required=True,
        default='tsv'
    )
    parser.add_argument(
        "--output",
        required=True,
    )
    parser.add_argument(
        "--output_type",
        required=True,
        default='mtx'
    )
    parser.add_argument(
        "--columns",
        required=False,
        default=None,
        help="barcodes.tsv"
    )
    parser.add_argument(
        "--rows",
        required=False,
        default=None,
        help="genes.tsv"
    )
    parser.add_argument(
        "--genome",
        required=False,
        default=None,
        help="genome (hg19, mm10, etc.)"
    )
    args = parser.parse_args()

    input_file = args.input
    input_type = args.input_type
    output_file = args.output
    output_type = args.output_type
    columns_file = args.columns
    rows_file = args.rows
    genome = args.genome

    convert(
        input_file, input_type, output_file, output_type,
        columns_file, rows_file, genome
    )

if __name__ == "__main__":
    main()
