import pytest
from genomic_address_service.classes.multi_level_clustering import multi_level_clustering
import tempfile
import textwrap
import os
import re

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
    
@pytest.fixture
def sample_distance_matrix_sample_int():
    content = textwrap.dedent(
        """\
        Header\t1\t2\t3
        1\t0.0\t0.1\t0.2
        2\t0.1\t0.0\t0.3
        3\t0.2\t0.3\t0.0
        """
    )
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def sample_distance_matrix_sample_int_float():
    content = textwrap.dedent(
        """\
        Header\t1\t2.0\t3
        1\t0.0\t0.1\t0.2
        2.0\t0.1\t0.0\t0.3
        3\t0.2\t0.3\t0.0
        """
    )
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def sample_distance_matrix_string_values():
    content = textwrap.dedent(
        """\
        Header\tLabel1\tLabel2\tLabel3
        Label1\t0.0\t0.1\ttest
        Label2\t0.1\t0.0\t0.3
        Label3\ttest\t0.3\t0.0
        """
    )
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def sample_distance_matrix_nullvalues():
    content = textwrap.dedent(
        """\
        Header\tLabel1\tLabel2\tLabel3
        Label1\t0.0\t0.1\tNaN
        Label2\t0.1\t0.0\t0.3
        Label3\tNaN\t0.3\t0.0
        """
    )
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)
    
### Tests for read_distance_matrix function in multi_level_clustering class:

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
    ## Confirms that a 2x3 matrix raises an error (matrix should be n x n)
    ### Example:
    ### Header Label3  Label2                        
    ### Label1     0.2     0.1
    ### Label2     0.3     0.0
    ### Label3     0.0     0.3
    ### Missing column for Label1
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

def test_integer_labels(sample_distance_matrix_sample_int):
    ### Confirms that integer labels are correctly read as strings and converted back to float
    ## Header   1  2  3                        
    ## 1     0.0     0.1     0.2
    ## 2     0.1     0.0     0.3
    ## 3     0.2     0.3     0.0
    thresholds = [0.15]
    mlc = multi_level_clustering(dist_mat_file=sample_distance_matrix_sample_int, thresholds=thresholds, method="single", sort_matrix=False)
    assert len(mlc.labels) == 3  # Expecting 3 labels based on the sample matrix
    assert '1' in mlc.cluster_memberships  # Initial membership should be populated with string '1'
    assert '2' in mlc.cluster_memberships  # Initial membership should be populated with string '2.0'
    assert '3' in mlc.cluster_memberships  # Initial membership should be populated with string '3' 

def test_integer_float_labels(sample_distance_matrix_sample_int_float):
    ### Confirms that integer and float labels are correctly read as strings
    ## Header   1  2.0  3                        
    ## 1     0.0     0.1     0.2
    ## 2.0     0.1     0.0     0.3
    ## 3     0.2     0.3     0.0
    thresholds = [0.15]
    mlc = multi_level_clustering(dist_mat_file=sample_distance_matrix_sample_int_float, thresholds=thresholds, method="single", sort_matrix=False)
    assert len(mlc.labels) == 3  # Expecting 3 labels based on the sample matrix
    assert '1' in mlc.cluster_memberships  # Initial membership should be populated with string '1'
    assert '2.0' in mlc.cluster_memberships  # Initial membership should be populated with string '2.0'
    assert '3' in mlc.cluster_memberships  # Initial membership should be populated with string '3'

def test_string_values_error(sample_distance_matrix_string_values):
    ### Confirms that a distance matrix with string values raises an error
    ## Header   Label1  Label2  Label3                        
    ## Label1     0.0     0.1     test
    ## Label2     0.1     0.0     0.3
    ## Label3     test     0.3     0.0
    
    thresholds = [0.15]
    expected_msg = "Input matrix should only contain numerical values"
    with pytest.raises(ValueError, match=re.escape(expected_msg)):
        mlc = multi_level_clustering(
            dist_mat_file=sample_distance_matrix_string_values, 
            thresholds=thresholds, 
            method="single",  
            sort_matrix=False
        )

def test_null_values_error(sample_distance_matrix_nullvalues):
    ### Confirms that a distance matrix with null values raises an error
    ## Header   Label1  Label2  Label3                        
    ## Label1     0.0     0.1     NaN
    ## Label2     0.1     0.0     0.3
    ## Label3     NaN     0.3     0.0
    thresholds = [0.15]
    expected_msg = "Distance matrix contains NaN, null or NA values."
    with pytest.raises(ValueError, match=re.escape(expected_msg)):
        mlc = multi_level_clustering(
            dist_mat_file=sample_distance_matrix_nullvalues, 
            thresholds=thresholds, 
            method="single",  
            sort_matrix=False
        )   