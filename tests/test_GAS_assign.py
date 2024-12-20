import pytest
import pandas as pd
from tempfile import NamedTemporaryFile
import os
from genomic_address_service.classes.assign import assign

@pytest.fixture
def mock_dist_file():
    content = """query_id\tref_id\tdist
q1\tq1\t0.0
q1\tr1\t0.1
q1\tr2\t0.2
"""
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def mock_membership_file():
    content = """id\taddress_levels_notsplit
r1\t1.1
r2\t2.1
"""
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

def test_initialization(mock_dist_file, mock_membership_file):
    threshold_map = {"level_0": 0.1, "level_1": 0.2}
    a = assign(dist_file=mock_dist_file, membership_file=mock_membership_file, threshold_map=threshold_map, linkage_method='single', sample_col='id', address_col='address_levels_notsplit', batch_size=100)
    assert a.status, "Initialization failed, check error_msgs for details"
    assert not a.error_msgs, f"Unexpected errors during initialization: {a.error_msgs}"
    assert isinstance(a.memberships_df, pd.DataFrame), "Memberships DataFrame"

def test_check_membership_columns(mock_dist_file, mock_membership_file):
    threshold_map = {"level_0": 0.1, "level_1": 0.2}
    a = assign(dist_file=mock_dist_file, membership_file=mock_membership_file, threshold_map=threshold_map, linkage_method='single', sample_col='id', address_col='address_levels_notsplit', batch_size=100)
    cols = ['level_0', 'level_1']
    assert a.check_membership_columns(cols), "Membership column check failed for valid columns"