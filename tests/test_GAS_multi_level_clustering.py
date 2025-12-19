import pytest
from genomic_address_service.classes.multi_level_clustering import multi_level_clustering
import tempfile
import textwrap
import os
import re

@pytest.fixture
def sample_distance_matrix():
    content = textwrap.dedent(
        """\
        Header\tLabel1\tLabel2\tLabel3
        Label1\t0.0\t0.1\t0.2
        Label2\t0.1\t0.0\t0.3
        Label3\t0.2\t0.3\t0.0
        """
    )
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def unsymmetrical_sample_distance_matrix():
    content = textwrap.dedent(
        """\
        Header\tLabel3\tLabel2\tLabel1
        Label1\t0.2\t0.1\t0.0
        Label2\t0.3\t0.0\t0.1
        Label3\t0.0\t0.3\t0.2
        """
    )
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def unsymmetrical_sample_distance_matrix_2():
    content = textwrap.dedent(
        """\
        Header\tLabel3\tLabel2
        Label1\t0.2\t0.1
        Label2\t0.3\t0.0
        Label3\t0.0\t0.3
        """
    )
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)
@pytest.fixture
def unsymmetrical_sample_distance_matrix_3():
    content = textwrap.dedent(
        """\
        Header\tLabel1\tLabel2\tLabel3
        Label1\t0.0\t0.1\t0.2
        Label2\t0.1\t0.0\t0.4
        Label3\t0.2\t0.3\t0.0
        """
    )
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)


def test_initialization(sample_distance_matrix):
    thresholds = [0.15]
    mlc = multi_level_clustering(dist_mat_file=sample_distance_matrix, thresholds=thresholds, method="single", sort_matrix=False)
    assert len(mlc.labels) == 3  # Expecting 3 labels based on the sample matrix
    assert mlc.linkage is not None  # Linkage matrix should be created
    assert 'Label1' in mlc.cluster_memberships  # Initial membership should be populated

def test_assign_clusters(sample_distance_matrix):
    thresholds = [0.15]
    mlc = multi_level_clustering(dist_mat_file=sample_distance_matrix, thresholds=thresholds, method="single", sort_matrix=False)
    assert all(len(clusters) == 1 for clusters in mlc.cluster_memberships.values())

def test_newick_string(sample_distance_matrix):
    thresholds = [0.15]
    mlc = multi_level_clustering(dist_mat_file=sample_distance_matrix, thresholds=thresholds, method="single",  sort_matrix=False)
    assert mlc.newick.endswith(";")  # Newick strings should end with a semicolon

def test_assymetical_matrix_error_diagnonal(unsymmetrical_sample_distance_matrix):
    ### Confirms that matrix with diagnonal not starting at [0,0] raises an error
    ### Example:
    ### Header Label3  Label2  Label1                        
    ### Label1     0.2     0.1     0.0
    ### Label2     0.3     0.0     0.1
    ### Label3     0.0     0.3     0.2
    ### Rows and columns are in different orders, and diagonal does not start at [0,0]
    thresholds = [0.15]
    expected_msg = "Incorrect Distance Matrix Format: --matrix must have (n x n) dimensions, 0 diagonal starting at position [0,0] and rows/columns must in the same order."
    with pytest.raises(ValueError, match=re.escape(expected_msg)):
        mlc = multi_level_clustering(
            dist_mat_file=unsymmetrical_sample_distance_matrix, 
            thresholds=thresholds, 
            method="single",  
            sort_matrix=False
        )

def test_assymetical_matrix_error_shape(unsymmetrical_sample_distance_matrix_2):
    ## Confirms that a 2x3 matrix raises an error
    thresholds = [0.15]
    expected_msg = "Incorrect Distance Matrix Format: --matrix must have (n x n) dimensions, 0 diagonal starting at position [0,0] and rows/columns must in the same order."
    with pytest.raises(ValueError, match=re.escape(expected_msg)):
        mlc = multi_level_clustering(
            dist_mat_file=unsymmetrical_sample_distance_matrix_2, 
            thresholds=thresholds, 
            method="single",  
            sort_matrix=False
        )

def test_assymetical_matrix_error_values(unsymmetrical_sample_distance_matrix_3):
    ## Confirms that both upper and lower triangles do not match raises an error
    ## Header   Label1  Label2  Label3                        
    ## Label1     0.0     0.1     0.2
    ## Label2     0.1     0.0     0.4
    ## Label3     0.2     0.3     0.0

    ## In this case, position [2,1] = 0.4 in upper triangle does not match position [1,2] = 0.3 in lower triangle
    thresholds = [0.15]
    expected_msg = 'Distance matrix has non-symmetrical values'
    with pytest.raises(ValueError, match=re.escape(expected_msg)):
        mlc = multi_level_clustering(
            dist_mat_file=unsymmetrical_sample_distance_matrix_3, 
            thresholds=thresholds, 
            method="single",  
            sort_matrix=True
        )