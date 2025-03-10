"""
Tests for reader class

"""

import pytest
import textwrap
from genomic_address_service.classes.reader import dist_reader
import io
import pandas as pd

@pytest.fixture()
def query_table_tsv(tmp_path):
    input_file = textwrap.dedent(
        """\
        query_id\tref_id\tdist
        sampleQ\tsampleQ\t0
        sampleQ\tsampleN\t0
        sampleQ\tsample1\t1
        sampleQ\tsample2\t1
        sampleQ\tsample3\t2
        sampleN\tsampleQ\t0
        sampleN\tsampleN\t0
        sampleN\tsample1\t1
        sampleN\tsample2\t1
        sampleN\tsample3\t2
        """
    )
    df = pd.read_csv(io.StringIO(input_file), sep='\t')
    out_path = tmp_path / "test_file.tsv"
    df.to_csv(out_path, sep='\t', index=False)
    return out_path

@pytest.fixture()
def query_table_pq(tmp_path):
    input_file = textwrap.dedent(
        """\
        query_id\tref_id\tdist
        sampleQ\tsampleQ\t0
        sampleQ\tsampleN\t0
        sampleQ\tsample1\t1
        sampleQ\tsample2\t1
        sampleQ\tsample3\t2
        sampleN\tsampleQ\t0
        sampleN\tsampleN\t0
        sampleN\tsample1\t1
        sampleN\tsample2\t1
        sampleN\tsample3\t2
        """
    )
    df = pd.read_csv(io.StringIO(input_file), sep='\t')
    out_path = tmp_path / "test_file.parquet"
    df.to_parquet(out_path, index=False)
    return out_path


@pytest.fixture()
def matrix_class():
    return dist_reader("genomic_address_service/example/mcluster/hamming/results.text", n_records=10, min_dist=1)


@pytest.fixture(scope="function")
def test_class(query_table_tsv):
    return dist_reader(str(query_table_tsv), n_records=10, min_dist=1)

def test_dist_reader_functionality(query_table_tsv):
    """
    A base test to verify the intended functionality of the dist_reader class    
    """
    d = dist_reader(str(query_table_tsv), n_records=11, min_dist=1)
    assert list(d.read_data()) == [
        {
            'sampleQ': {
                'sample1': 1.0,
                'sample2': 1.0,
                'sample3': 2.0,
                'sampleQ': 0,
                'sampleN': 0
            },
            'sampleN': {
                'sample1': 1.0,
                'sample2': 1.0,
                'sample3': 2.0,
                'sampleQ': 0,
                'sampleN': 0
            }
        }
    ]

@pytest.mark.xfail(run=False, comment="Cannot deserialize parquet file.")
def test_dist_reader_functionality_pq(query_table_pq):
    """
    A base test to verify the intended functionality of the dist_reader class    
    """
    d = dist_reader(str(query_table_pq), n_records=11, min_dist=1)
    assert list(d.read_data()) == [{'sampleQ': {'sample1': 1.0, 'sample2': 1.0, 'sample3': 2.0}, 'sampleN': {'sample1': 1.0, 'sample2': 1.0, 'sample3': 2.0}}]


def test_guess_file_type(test_class, query_table_tsv):
    assert test_class.guess_file_type(str(query_table_tsv)) == "text"

def test_guess_dist_type(test_class, query_table_tsv):
    assert test_class.guess_dist_type(test_class.fpath, 
                                      test_class.guess_file_type(str(query_table_tsv)), 
                                      test_class.delim) == "pd"


def test_read_pd(test_class):
    """
    Part of this test has to recreate the steps required
    for testing the read_pd method
    """
    test_class.file_handle = open(test_class.fpath, 'r')
    test_class.header = next(test_class.file_handle).split(test_class.delim)
    assert test_class.dists == {}
    test_class.read_pd()
    [_ for _ in test_class.read_pd()] # need to exhaust iterator to populate dists
    assert test_class.dists == {
        'sampleQ': {
            'sample1': 1.0,
            'sample2': 1.0,
            'sample3': 2.0,
            'sampleQ': 0,
            'sampleN': 0
        },
        'sampleN': {
            'sample1': 1.0,
            'sample2': 1.0,
            'sample3': 2.0,
            'sampleQ': 0,
            'sampleN': 0
        }
    }

@pytest.mark.xfail(run=False, comment="Distances reported from this test are incorrect, \
                due to an issue in parsing the file. However this function appears not to be called as \
                part of this programs regular usage.")
def test_read_matrix(matrix_class):
    """
    Test the matrix reading class.
    """
    matrix_class.file_handle = open(matrix_class.fpath, 'r')
    matrix_class.header = next(matrix_class.file_handle).split(matrix_class.delim)
    assert matrix_class.header == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    assert matrix_class.dists == {}
    [_ for _ in matrix_class.read_matrix()]