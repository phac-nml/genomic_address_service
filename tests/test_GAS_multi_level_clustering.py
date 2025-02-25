import pytest
from genomic_address_service.classes.multi_level_clustering import multi_level_clustering
import tempfile
import textwrap
import os

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

def test_initialization(sample_distance_matrix):
    thresholds = [0.15]
    mlc = multi_level_clustering(dist_mat_file=sample_distance_matrix, thresholds=thresholds, method="single")
    assert len(mlc.labels) == 3  # Expecting 3 labels based on the sample matrix
    assert mlc.linkage is not None  # Linkage matrix should be created
    assert 'Label1' in mlc.cluster_memberships  # Initial membership should be populated

def test_assign_clusters(sample_distance_matrix):
    thresholds = [0.15]
    mlc = multi_level_clustering(dist_mat_file=sample_distance_matrix, thresholds=thresholds, method="single")
    mlc.assign_clusters()
    # This assertion may need to be adjusted based on expected cluster assignments
    # I gussed 2?
    assert all(len(clusters) == 2 for clusters in mlc.cluster_memberships.values())

def test_newick_string(sample_distance_matrix):
    thresholds = [0.15]
    mlc = multi_level_clustering(dist_mat_file=sample_distance_matrix, thresholds=thresholds, method="single")
    assert mlc.newick.endswith(";")  # Newick strings should end with a semicolon