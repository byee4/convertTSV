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
import numpy as np

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
            print("That genome does not exist in this file.")
            return None
        gene_ids = getattr(group, 'genes').read()
        gene_names = getattr(group, 'gene_names').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        # print('data', data)
        indices = getattr(group, 'indices').read()
        # print('indices', indices)
        indptr = getattr(group, 'indptr').read()
        # print('indptr', indptr)
        shape = getattr(group, 'shape').read()
        # print('shape', shape)
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


def write_dataframe_as_h5(df, output_file, genome, gene_ids, gene_names,
                          barcodes):
    h5_mtx = scipy.sparse.csc_matrix(df.values)


    flt = tables.Filters(complevel=1)
    with tables.open_file(output_file, 'w', filters=flt) as f:
        f.set_node_attr(f.root, "chemistry_description", "Single Cell 3\' V2")
        f.set_node_attr(f.root, "filetype", "matrix")
        # f.set_node_attr(f.root, "library_ids", ['id1'])
        # f.set_node_attr(f.root, "original_gem_groups", [1])
        f.root._f_setattr('original_gem_groups', np.array([1]))
        f.root._f_setattr('library_ids', np.array(['id1']))
        try:
            group = f.create_group(f.root, genome)
            # for attribute in ('indices', 'indptr'):
            #     arr = np.ndarray(getattr(h5_mtx, attribute))
            #     f.create_carray(group, attribute, obj=arr)
            f.create_carray(group, 'data', obj=np.asarray(h5_mtx.data, dtype=np.dtype('int32')))
            f.create_carray(group, 'genes', obj=np.asarray(gene_ids, dtype=np.dtype('S18')))
            f.create_carray(group, 'gene_names', obj=np.asarray(gene_names, dtype=np.dtype('S18')))
            f.create_carray(group, 'barcodes', obj=np.asarray(barcodes, dtype=np.dtype('S18')))
            f.create_carray(group, 'indices', obj=np.asarray(h5_mtx.indices, dtype=np.dtype('uint32')))
            f.create_carray(group, 'indptr', obj=np.asarray(h5_mtx.indptr, dtype=np.dtype('uint32')))
            f.create_carray(group, 'shape', obj=np.array(getattr(h5_mtx, 'shape'), dtype=np.dtype('int32')))
        except Exception as e:
            print('cannot write h5', e)

def write_dataframe_as_h5_h5py(df, output_file, genome, gene_ids, gene_names,
                          barcodes):
    """
    Writes a pandas dataframe as an h5 file.

    :param df:
    :param output_file:
    :param genome:
    :param gene_ids:
    :param gene_names:
    :param barcodes:
    :return:
    """
    f = h5py.File(output_file, "w")
    h5_mtx = scipy.sparse.csc_matrix(df.values)
    h5_group = f.create_group(genome)

    h5_group.attrs['CLASS'] = 'GROUP'
    h5_group.attrs['FILTERS'] = 65793
    h5_group.attrs['TITLE'] = '.'
    h5_group.attrs['VERSION'] = '1.0'

    f.attrs['CLASS'] = 'GROUP'
    f.attrs['FILTERS'] = 65793
    f.attrs['PYTABLES_FORMAT_VERSION'] = '2.1'
    f.attrs['TITLE'] = '.'
    f.attrs['VERSION'] = '1.0'
    f.attrs['chemistry_description'] = "Single Cell 3\' v2"
    f.attrs['filetype'] = "matrix"
    f.attrs['library_ids'] = ["expression_csv"]
    f.attrs['original_gem_groups'] = [1]

    for attribute in ('data', 'shape'):
        arr = np.array(getattr(h5_mtx, attribute))
        ds = h5_group.create_dataset(
            name=attribute, data=arr, dtype='int32',
            compression="gzip", compression_opts=1,
            shuffle=True,
        )
        ds.attrs['CLASS'] = 'CARRAY'
        ds.attrs['TITLE'] = '.'
        ds.attrs['VERSION'] = '1.1'

    arr = np.array(getattr(h5_mtx, 'indices'))
    h5_indices = h5_group.create_dataset(
        name='indices', data=arr, dtype='int64',
        compression="gzip", compression_opts=1,
        shuffle=True,
    )
    h5_indices.attrs['CLASS'] = 'CARRAY'
    h5_indices.attrs['TITLE'] = '.'
    h5_indices.attrs['VERSION'] = '1.1'

    arr = np.array(getattr(h5_mtx, 'indptr'))
    h5_indptr = h5_group.create_dataset(
        name='indptr', data=arr, dtype='int64',
        compression="gzip", compression_opts=1,
        shuffle=True,
    )
    h5_indptr.attrs['CLASS'] = 'CARRAY'
    h5_indptr.attrs['TITLE'] = '.'
    h5_indptr.attrs['VERSION'] = '1.1'


    h5_gene_ids = h5_group.create_dataset(
        name='genes', data=list(gene_ids),
        compression="gzip", compression_opts=1,
        shuffle=True,
    )
    h5_gene_ids.attrs['CLASS'] = 'CARRAY'
    h5_gene_ids.attrs['TITLE'] = '.'
    h5_gene_ids.attrs['VERSION'] = '1.1'

    h5_gene_names = h5_group.create_dataset(
        name='gene_names', data=list(gene_names),
        compression="gzip", compression_opts=1,
        shuffle=True,
    )
    h5_gene_names.attrs['CLASS'] = 'CARRAY'
    h5_gene_names.attrs['TITLE'] = '.'
    h5_gene_names.attrs['VERSION'] = '1.1'

    h5_barcodes = h5_group.create_dataset(
        name='barcodes', data=list(barcodes),
        compression="gzip", compression_opts=1,
        shuffle=True,
    )
    h5_barcodes.attrs['CLASS'] = 'CARRAY'
    h5_barcodes.attrs['TITLE'] = '.'
    h5_barcodes.attrs['VERSION'] = '1.1'
    print("barcodes", list(barcodes))
    f.close()


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
            print('warning: gene name+id may not be present in h5 file')
            sv2mtx(input_file, output_file, ',')
        elif output_type == 'h5':
            print('warning: gene name+id may not be present in h5 file')
            sv2h5(input_file, output_file, genome, ',')
    elif input_type == 'tsv':
        if output_type == 'mtx':
            print('warning: gene name+id may not be present in h5 file')
            sv2mtx(input_file, output_file, '\t')
        elif output_type == 'h5':
            print('warning: gene name+id may not be present in h5 file')
            sv2h5(input_file, output_file, genome, '\t')
    elif input_type == 'mtx':
        if output_type == 'tsv':
            mtx2sv(input_file, rows_file, columns_file, output_file, sep='\t')
        elif output_type == 'csv':
            mtx2sv(input_file, rows_file, columns_file, output_file, sep=',')
    elif input_type == 'h5':
        if output_type == 'tsv':
            h52sv(input_file, output_file, genome, '\t')
        elif output_type == 'csv':
            h52sv(input_file, output_file, genome, ',')


def sv2mtx(input_file, output_file, sep):
    """
    Converts a tab/comma/sep separated file into an mtx file

    :param input_file:
    :param output_file:
    :return:
    """
    df, cols, rows = read_sv_as_dataframe(input_file, sep)
    write_dataframe_as_matrix(df, output_file)

    o = open(output_file + '.columns', 'w')
    for column in cols:
        o.write(column + '\n')
    o.close()

    o = open(output_file + '.rows', 'w')
    for row in rows:
        o.write(row + '\n')
    o.close()


def mtx2sv(input_file, rows_file, columns_file, output_file, sep):
    """
    Converts mtx file + rows + columns into a tab/comma-separated file.

    :param input_file:
    :param rows_file:
    :param columns_file:
    :param output_file:
    :return:
    """
    df = read_mtx_as_dataframe(
        input_file, columns_file, rows_file
    )
    write_dataframe_as_sv(df, output_file, sep=sep)


def mtx2h5(input_file, output_file, columns_file, rows_file, genome):
    """
    Converts mtx/mtx associated files into h5 format.

    :param input_file:
    :param output_file:
    :param columns_file:
    :param rows_file:
    :param genome:
    :return:
    """
    df = read_mtx_as_dataframe(
        input_file, columns_file, rows_file
    )
    columns = [
        row[0] for row in csv.reader(open(columns_file), delimiter="\t")
        ]
    rows = [
        row[0] for row in csv.reader(open(rows_file), delimiter="\t")
        ]
    write_dataframe_as_h5(df, output_file, genome, rows, rows, columns)


def h52mtx(input_file, output_file, genome):
    """
    Converts from h5 file to mtx

    :param input_file:
    :param output_file:
    :param genome:
    :return:
    """
    gene_ids, gene_names, barcodes, matrix = get_matrix_from_h5(
        input_file, genome
    )

    o = open(output_file + '.rows', 'w')
    for i in range(0, len(gene_ids)):
        o.write('{}\t{}\n'.format(gene_ids[i], gene_names[i]))
    o.close()

    o = open(output_file + '.columns', 'w')
    for column in barcodes:
        o.write(column + '\n')
    o.close()

    mtx = scipy.sparse.csr_matrix(matrix)
    scipy.io.mmwrite(output_file, mtx)


def sv2h5(input_file, output_file, genome, sep='\t'):
    """
    Convert tab or csv separated files to h5.

    :param input_file:
    :param output_file:
    :param genome:
    :param sep:
    :return:
    """
    df, barcodes, gene_ids = read_sv_as_dataframe(input_file, sep=sep)
    write_dataframe_as_h5(df, output_file, genome, gene_ids, gene_ids,
                          barcodes)


def h52sv(input_file, output_file, genome, sep='\t'):
    """
    Converts an h5 input file into a tab separated file.

    :param input_file:
    :param output_file:
    :param genome:
    :return:
    """
    gene_ids, gene_names, barcodes, mat = get_matrix_from_h5(input_file, genome)
    df = pd.DataFrame(mat.todense(), columns=barcodes, index=gene_ids)
    write_dataframe_as_sv(df, output_file, sep)


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
