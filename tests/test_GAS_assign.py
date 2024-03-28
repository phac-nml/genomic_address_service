import pytest
import pandas as pd
from tempfile import NamedTemporaryFile
import os
from genomic_address_service.classes.assign import assign

@pytest.fixture
def mock_dist_file():
    content = """query_id\tref_id\tdist
q1\tr1\t0.1
q2\tr2\t0.2
"""
    with NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def mock_membership_file():
    content = """id\tthreshold_0.1\tthreshold_0.2
r1\t1\t1
r2\t2\t1
"""
    with NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

def test_initialization(mock_dist_file, mock_membership_file):
    threshold_map = {"threshold_0.1": 0.1, "threshold_0.2": 0.2}
    a = assign(dist_file=mock_dist_file, membership_file=mock_membership_file, threshold_map=threshold_map, linkage_method='single')
    
    assert a.status, "Initialization failed, check error_msgs for details"
    assert not a.error_msgs, f"Unexpected errors during initialization: {a.error_msgs}"
    assert isinstance(a.query_df, pd.DataFrame), "Query DataFrame not initialized properly"
    assert isinstance(a.memberships_df, pd.DataFrame), "Memberships DataFrame"

def test_check_membership_columns(mock_dist_file, mock_membership_file):
    threshold_map = {"threshold_0.1": 0.1, "threshold_0.2": 0.2}
    a = assign(dist_file=mock_dist_file, membership_file=mock_membership_file, threshold_map=threshold_map, linkage_method='single')
    cols = ['threshold_0.1', 'threshold_0.2']
    assert a.check_membership_columns(cols), "Membership column check failed for valid columns"