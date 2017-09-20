from convertTSV import convertTSV as c
from pandas.util.testing import assert_frame_equal
import pandas as pd
import scipy.io
import numpy as np
import csv
import scipy.sparse

def get_reference_h5():
    return '../../data/filtered_gene_bc_matrices_h5.h5'

def get_reference_mtx():
    mtx = '../../data/filtered_gene_bc_matrices/mm10/matrix.mtx'
    barcodes = '../../data/filtered_gene_bc_matrices/mm10/barcodes.tsv'
    genes = '../../data/filtered_gene_bc_matrices/mm10/genes.tsv'

    mat = scipy.io.mmread(mtx)
    barcodes = [
        row[0] for row in csv.reader(open(barcodes), delimiter="\t")
        ]
    genes = [
        row[0] for row in csv.reader(open(genes), delimiter="\t")
        ]

    barcodes_ndarray = np.asarray(barcodes)
    genes_ndarray = np.asarray(genes)

    csc_matrix = scipy.sparse.csc_matrix(mat)

    return csc_matrix, barcodes_ndarray, genes_ndarray


def get_reference_csv():
    return '../../data/mat2csv_conversion.csv'


def test_get_matrix_from_h5():
    # First convert mtx to tsv
    genome = 'mm10'
    input_file = get_reference_h5()
    gene_ids, gene_names, barcodes, mat = c.get_matrix_from_h5(
        input_file, genome
    )
    ref_mtx, ref_barcodes, ref_genes = get_reference_mtx()
    print('from mtx', type(mat), mat.shape)
    print('from h5', type(ref_mtx), ref_mtx.shape)
    assert np.array_equal(barcodes, ref_barcodes)
    assert np.array_equal(gene_ids, ref_genes)
    assert (mat!=ref_mtx).nnz==0


def test_tsv2h5():
    genome = 'mm10'
    input_file = get_reference_h5()
    input_csv = get_reference_csv()
    output_file = 'data/tmp.h5'



    # get dataframe, barcodes, gene_ids from mat2csv
    df, barcodes, gene_ids = c.read_sv_as_dataframe(
        input_csv, sep=','
    )

    # save the h5 to temp file
    c.write_dataframe_as_h5(df, output_file, genome, gene_ids, gene_ids, barcodes)

    # open the h5 file using get_matrix_from_h5
    test_gene_ids, test_gene_names, test_barcodes, test_mtx = c.get_matrix_from_h5(
        output_file, genome
    )

    # open the (reference) h5
    ref_gene_ids, ref_gene_names, ref_barcodes, ref_mtx = c.get_matrix_from_h5(
        input_file, genome
    )

    # compare the reference h5 to the h5 created from tsv2h5
    assert list(ref_gene_ids) == list(test_gene_ids)
    assert list(ref_barcodes) == list(test_barcodes)
    assert (test_mtx != ref_mtx).nnz==0


def test_h52sv():
    genome = 'mm10'
    input_file = get_reference_h5()
    input_csv = get_reference_csv()
    output_file = 'data/tmp.csv'

    # convert h5 to tsv
    c.h52sv(input_file, output_file, genome, ',')

    # get dataframe, barcodes, gene_ids from mat2csv to get reference csv
    ref_df, ref_barcodes, ref_gene_ids = c.read_sv_as_dataframe(
        input_csv, sep=','
    )

    # read in the tmp csv file that h52sv made
    test_df = pd.read_table(output_file, sep=',', index_col=0)

    # assert that the two dataframes are equal.
    assert_frame_equal(ref_df, test_df)


def test_h52mtx():
    print("Converts reference h5 to mtx and compares the resulting "
          "mtx/genes/barcodes files to the reference mtx/genes/barcodes. "
          "Compares the two mtx files by converting each into dataframes.")
    genome = 'mm10'
    input_file = get_reference_h5()

    ref_mtx = '../../data/filtered_gene_bc_matrices/mm10/matrix.mtx'
    ref_barcodes = '../../data/filtered_gene_bc_matrices/mm10/barcodes.tsv'
    ref_genes = '../../data/filtered_gene_bc_matrices/mm10/genes.tsv'

    output_file = 'data/tmp.mtx'
    output_cols = 'data/tmp.mtx.columns'
    output_rows = 'data/tmp.mtx.rows'

    # convert h5 to mtx/genes/barcodes
    c.h52mtx(input_file, output_file, genome)

    # get reference mtx files
    ref_df = c.read_mtx_as_dataframe(ref_mtx, ref_barcodes, ref_genes)

    # get test mtx files
    test_df = c.read_mtx_as_dataframe(output_file, output_cols, output_rows)

    assert_frame_equal(ref_df, test_df)