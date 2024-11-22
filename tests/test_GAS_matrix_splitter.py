import pytest
import os
from tempfile import NamedTemporaryFile, TemporaryDirectory
from genomic_address_service.classes.matrix_splitter import matrix_splitter

@pytest.fixture
def mock_matrix_file():
    content = "Header\nRow1\t0.1\t0.2\nRow2\t0.2\t0.1\n"
    with NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
    yield tmp.name
    os.unlink(tmp.name)

@pytest.fixture
def big_mock_matrix_file():
    content = "Header\nRow1\t0.1\t0.2\t0.3\t0.4\t0.5\nRow2\t0.6\t0.7\t0.8\t0.9\t1.0\nRow3\t1.1\t1.2\t1.3\t1.4\t1.5\nRow4\t1.6\t1.7\t1.8\t1.9\t2.0\nRow5\t2.1\t2.2\t2.3\t2.4\t2.5\n"
    with NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
    yield tmp.name
    os.unlink(tmp.name)

@pytest.fixture
def output_directory():
    with TemporaryDirectory() as tmpdir:
        yield tmpdir


def test_initialization(mock_matrix_file, output_directory):
    batch_size = 7
    ms = matrix_splitter(mock_matrix_file, output_directory, batch_size)
    assert ms.file_path == mock_matrix_file
    assert ms.out_path == output_directory
    assert ms.batch_size == batch_size
    assert ms.is_ok == True  # Assuming the mock file and output directory meet the requirements

def test_get_file_length(mock_matrix_file):
    dummy_out_path = "dummy_out_path"
    dummy_batch_size = 1
    ms = matrix_splitter(file_path=mock_matrix_file, 
                         out_path=dummy_out_path, 
                         batch_size=dummy_batch_size)
    expected_line_count = 3
    assert ms.get_file_length() == expected_line_count, "get_file_length did not return the expected number of lines"

def test_prep_batch_ranges(big_mock_matrix_file, output_directory):
    batch_size = 1
    ms = matrix_splitter(big_mock_matrix_file, output_directory, batch_size)
    ms.prep_batch_ranges()
    assert ms.num_batches == 5 
    assert len(ms.ranges) == ms.num_batches
    assert ms.ranges == [(i, i+1) for i in range(0, ms.num_batches)]

@pytest.mark.parametrize("method_name", ["parse_distance_matrix_bins", "parse_distance_matrix_partitions"])
def test_parse_methods(big_mock_matrix_file, output_directory, method_name):
    batch_size = 1
    ms = matrix_splitter(big_mock_matrix_file, output_directory, batch_size)
    ms.prep_batch_ranges()
    parse_method = getattr(ms, method_name)
    parse_method()
    for i in range(ms.num_batches):
        assert os.path.exists(os.path.join(ms.out_path, f"{ms.prefix}-{i}.matrix"))
