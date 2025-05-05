"""
Tests for reader class

"""

import pytest
import textwrap
from genomic_address_service.classes.reader import dist_reader
import io
import pandas as pd
from os import path

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

@pytest.fixture(scope="function")
def test_class(query_table_tsv):
    return dist_reader(str(query_table_tsv), n_records=10)

def get_path(location):
    directory = path.dirname(path.abspath(__file__))
    return path.join(directory, location)

def test_dist_reader_functionality(query_table_tsv):
    """
    A base test to verify the intended functionality of the dist_reader class    
    """
    d = dist_reader(str(query_table_tsv), n_records=11)
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

def test_reader_basic():

    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    chunk_size = 5

    distance_reader = dist_reader(pairwise_distances_path, n_records=chunk_size)
    chunks = distance_reader.read_data()

    chunk = next(chunks)
    assert chunk == {
        'E': {
            'E': 0.0,
            'B': 0.0,
            'A': 1.0,
            'C': 4.0,
            'D': 5.0,
            'F': 5.0},
        'F': {
            'F': 0.0,
            'D': 0.0,
            'A': 3.0,
            'B': 5.0,
            'E': 5.0,
            'C': 6.0}
        }

def test_reader_chunk_size_same():
    # number of items is 2: {E, F}
    # chunk_size = 2

    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    chunk_size = 2

    distance_reader = dist_reader(pairwise_distances_path, n_records=chunk_size)
    chunks = distance_reader.read_data()

    chunk = next(chunks)
    assert chunk == {
        'E': {
            'E': 0.0,
            'B': 0.0,
            'A': 1.0,
            'C': 4.0,
            'D': 5.0,
            'F': 5.0},
        'F': {
            'F': 0.0,
            'D': 0.0,
            'A': 3.0,
            'B': 5.0,
            'E': 5.0,
            'C': 6.0}
        }

def test_reader_chunk_size_smaller():
    # number of items is 2: {E, F}
    # chunk_size = 1

    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    chunk_size = 1

    distance_reader = dist_reader(pairwise_distances_path, n_records=chunk_size)
    chunks = distance_reader.read_data()

    chunk = next(chunks)
    assert chunk == {
        'E': {
            'E': 0.0,
            'B': 0.0,
            'A': 1.0,
            'C': 4.0,
            'D': 5.0,
            'F': 5.0}
        }

    chunk = next(chunks)
    assert chunk == {
        'F': {
            'F': 0.0,
            'D': 0.0,
            'A': 3.0,
            'B': 5.0,
            'E': 5.0,
            'C': 6.0}
        }
