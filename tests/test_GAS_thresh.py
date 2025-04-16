import pytest
import pandas as pd
from tempfile import NamedTemporaryFile
import os
import textwrap
from genomic_address_service.classes.assign import assign

@pytest.fixture
def results_dist():
    content = textwrap.dedent(
        """\
        query_id\tref_id\tdist
        sampleQ\tsampleQ\t0
        sampleQ\tsample1\t1
        sampleQ\tsample2\t1
        sampleQ\tsample3\t2
        """
    )
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def reference_clusters():
    content = textwrap.dedent(
        """\
        id\taddress\tlevel_1
        sample1\t1\t1
        sample2\t1\t1
        sample3\t1\t1
        """
    )
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

def test_check_thresh(results_dist, reference_clusters):
    #Test the threshold of being included
    threshold_map = {0: 1.0}
    z = assign(dist_file=results_dist, membership_file=reference_clusters, threshold_map=threshold_map, linkage_method='single', address_col = "address", sample_col='id', batch_size=100, delimiter=".")
    assert z.memberships_dict == {'sample1': '1', 'sample2': '1', 'sample3': '1', 'sampleQ': '1'}
    #Test the threshold of being excluded
    threshold_map = {0: 0.0}
    z = assign(dist_file=results_dist, membership_file=reference_clusters, threshold_map=threshold_map, linkage_method='single', address_col = "address", sample_col='id', batch_size=100, delimiter=".")
    assert z.memberships_dict == {'sample1': '1', 'sample2': '1', 'sample3': '1', 'sampleQ': '2'}


@pytest.fixture
def reference_clusters():
    content = textwrap.dedent(
        """\
        id\taddress
        sample1\t1.1
        sample2\t1.1
        sample3\t1.1
        """
    )
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

def test_check_thresh(results_dist, reference_clusters):
    #Test the threshold of being excluded
    threshold_map = {0: 1.0, 1: 0.0}
    z = assign(dist_file=results_dist, membership_file=reference_clusters, threshold_map=threshold_map, linkage_method='single', address_col = "address", sample_col='id', batch_size=100, delimiter=".")
    assert z.memberships_dict == {'sample1': '1.1', 'sample2': '1.1', 'sample3': '1.1', 'sampleQ': '1.2'}
    #Test the threshold of being included
    threshold_map = {0: 1.0, 1: 1.0}
    z = assign(dist_file=results_dist, membership_file=reference_clusters, threshold_map=threshold_map, linkage_method='single', address_col = "address", sample_col='id', batch_size=100, delimiter=".")
    assert z.memberships_dict == {'sample1': '1.1', 'sample2': '1.1', 'sample3': '1.1', 'sampleQ': '1.1'}
    threshold_map = {0: 1.0, 1: 2.0}
    z = assign(dist_file=results_dist, membership_file=reference_clusters, threshold_map=threshold_map, linkage_method='single', address_col = "address", sample_col='id', batch_size=100, delimiter=".")
    assert z.memberships_dict == {'sample1': '1.1', 'sample2': '1.1', 'sample3': '1.1', 'sampleQ': '1.1'}
    threshold_map = {0: 1.0, 1: 3.0}
    z = assign(dist_file=results_dist, membership_file=reference_clusters, threshold_map=threshold_map, linkage_method='single', address_col = "address", sample_col='id', batch_size=100, delimiter=".")
    assert z.memberships_dict == {'sample1': '1.1', 'sample2': '1.1', 'sample3': '1.1', 'sampleQ': '1.1'}