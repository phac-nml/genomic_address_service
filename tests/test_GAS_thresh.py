import pytest
import pandas as pd
from tempfile import NamedTemporaryFile
import os
from genomic_address_service.classes.assign import assign

@pytest.fixture
def results_dist():
    content = """query_id\tref_id\tdist
sampleQ\tsampleQ\t0
sampleQ\tsample1\t1
sampleQ\tsample2\t1
sampleQ\tsample3\t2
"""
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def reference_clusters():
    content = """id\taddress\tlevel_1
sample1\t1\t1
sample2\t1\t1
sample3\t1\t1
"""
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

def test_check_thresh(results_dist, reference_clusters):
    threshold_map = {0: 1.0}
    z = assign(dist_file=results_dist, membership_file=reference_clusters, threshold_map=threshold_map, linkage_method='average', address_col = "address", sample_col='id', batch_size=100)
    assert z.memberships_dict == {'sample1': '1', 'sample2': '1', 'sample3': '1', 'sampleQ': '1'}