import pytest
import pandas as pd
from tempfile import NamedTemporaryFile
import os
import textwrap
from genomic_address_service.classes.assign import assign

@pytest.fixture
def mock_dist_file():
    content = textwrap.dedent(
        """\
        query_id\tref_id\tdist
        q1\tq1\t0.0
        q1\tr1\t0.1
        q1\tr2\t0.2
        """
    )
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

@pytest.fixture
def mock_membership_file():
    content = textwrap.dedent(
        """\
        id\taddress_levels_notsplit
        r1\t1.1
        r2\t2.1
        """
    )
    with NamedTemporaryFile('w+', suffix='.tsv', delete=False) as tmp:
        tmp.write(content)
        tmp.flush()
        yield tmp.name
        os.unlink(tmp.name)

def get_path(location):
    directory = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(directory, location)

def test_initialization(mock_dist_file, mock_membership_file):
    threshold_map = {"level_0": 0.1, "level_1": 0.2}
    a = assign(dist_file=mock_dist_file, membership_file=mock_membership_file, threshold_map=threshold_map, linkage_method='single', sample_col='id', address_col='address_levels_notsplit', batch_size=100, delimiter=".")
    assert a.status, "Initialization failed, check error_msgs for details"
    assert not a.error_msgs, f"Unexpected errors during initialization: {a.error_msgs}"
    assert isinstance(a.memberships_df, pd.DataFrame), "Memberships DataFrame"

def test_check_membership_columns(mock_dist_file, mock_membership_file):
    threshold_map = {"level_0": 0.1, "level_1": 0.2}
    a = assign(dist_file=mock_dist_file, membership_file=mock_membership_file, threshold_map=threshold_map, linkage_method='single', sample_col='id', address_col='address_levels_notsplit', batch_size=100, delimiter=".")
    cols = ['level_0', 'level_1']
    assert a.check_membership_columns(cols), "Membership column check failed for valid columns"

def test_check_file_type():
    threshold_map = {"level_0": 5.0, "level_1": 3.0, "level_2": 0.0}
    assignment = assign(dist_file=get_path("data/pairwise_distances/basic.tsv"),
                        membership_file=get_path("data/clusters/basic.tsv"),
                        threshold_map=threshold_map,
                        linkage_method='single',
                        sample_col='id',
                        address_col='address',
                        batch_size=100, delimiter=".")
    
    # Test check_file_type (tests file types of CLUSTER files):

    # .txt
    cluster_path = get_path("data/clusters/basic.txt")
    assignment.check_file_type(cluster_path)
    # no exception

    # .tsv
    cluster_path = get_path("data/clusters/basic.tsv")
    assignment.check_file_type(cluster_path)
    # no exception

    # .mat
    cluster_path = get_path("data/clusters/basic.mat")
    assignment.check_file_type(cluster_path)
    # no exception

    # .text
    cluster_path = get_path("data/clusters/basic.text")
    assignment.check_file_type(cluster_path)
    # no exception

    # .csv
    cluster_path = get_path("data/clusters/csv.csv")
    with pytest.raises(Exception) as exception:
        assignment.check_file_type(cluster_path)

    assert exception.type == Exception
    assert str(exception.value) == f"{cluster_path} does not have a valid extension (.csv): ['.txt', '.tsv', '.mat', '.text']"

def test_add_memberships_lookup_add_similar():
    threshold_map = {"level_0": 5.0, "level_1": 3.0, "level_2": 0.0}
    assignment = assign(dist_file=get_path("data/pairwise_distances/basic.tsv"),
                        membership_file=get_path("data/clusters/basic.tsv"),
                        threshold_map=threshold_map,
                        linkage_method='single',
                        sample_col='id',
                        address_col='address',
                        batch_size=100, delimiter=".")

    # memberships_dict will have:
    # A: 1.1.1
    # B: 1.1.2
    # C: 1.1.3
    # D: 1.1.4
    # E: 1.1.2
    # F: 1.1.4

    # memberships_lookup will have:
    # {'1': ['A', 'B', 'C', 'D', 'E', 'F'],
    # '1.1': ['A', 'B', 'C', 'D', 'E', 'F'],
    # '1.1.1': ['A'],
    # '1.1.2': ['B', 'E'],
    # '1.1.3': ['C'],
    # '1.1.4': ['D', 'F']}

    assignment.add_memberships_lookup("G", [1, 1, 5])
    assert assignment.memberships_dict["G"] == "1.1.5"
    assert assignment.memberships_lookup["1"] == ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    assert assignment.memberships_lookup["1.1"] == ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    assert assignment.memberships_lookup["1.1.5"] == ['G']

def test_add_memberships_lookup_add_different():
    threshold_map = {"level_0": 5.0, "level_1": 3.0, "level_2": 0.0}
    assignment = assign(dist_file=get_path("data/pairwise_distances/basic.tsv"),
                        membership_file=get_path("data/clusters/basic.tsv"),
                        threshold_map=threshold_map,
                        linkage_method='single',
                        sample_col='id',
                        address_col='address',
                        batch_size=100, delimiter=".")

    # memberships_dict will have:
    # A: 1.1.1
    # B: 1.1.2
    # C: 1.1.3
    # D: 1.1.4
    # E: 1.1.2
    # F: 1.1.4

    # memberships_lookup will have:
    # {'1': ['A', 'B', 'C', 'D', 'E', 'F'],
    # '1.1': ['A', 'B', 'C', 'D', 'E', 'F'],
    # '1.1.1': ['A'],
    # '1.1.2': ['B', 'E'],
    # '1.1.3': ['C'],
    # '1.1.4': ['D', 'F']}

    assignment.add_memberships_lookup("G", [2, 1, 1])
    assert assignment.memberships_dict["G"] == "2.1.1"
    assert assignment.memberships_lookup["2"] == ['G']
    assert assignment.memberships_lookup["2.1"] == ['G']
    assert assignment.memberships_lookup["2.1.1"] == ['G']

def test_add_memberships_lookup_add_same():
    threshold_map = {"level_0": 5.0, "level_1": 3.0, "level_2": 0.0}
    assignment = assign(dist_file=get_path("data/pairwise_distances/basic.tsv"),
                        membership_file=get_path("data/clusters/basic.tsv"),
                        threshold_map=threshold_map,
                        linkage_method='single',
                        sample_col='id',
                        address_col='address',
                        batch_size=100, delimiter=".")

    # memberships_dict will have:
    # A: 1.1.1
    # B: 1.1.2
    # C: 1.1.3
    # D: 1.1.4
    # E: 1.1.2
    # F: 1.1.4

    # memberships_lookup will have:
    # {'1': ['A', 'B', 'C', 'D', 'E', 'F'],
    # '1.1': ['A', 'B', 'C', 'D', 'E', 'F'],
    # '1.1.1': ['A'],
    # '1.1.2': ['B', 'E'],
    # '1.1.3': ['C'],
    # '1.1.4': ['D', 'F']}

    assignment.add_memberships_lookup("G", [1, 1, 1])
    assert assignment.memberships_dict["G"] == "1.1.1"
    assert assignment.memberships_lookup["1"] == ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    assert assignment.memberships_lookup["1.1"] == ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    assert assignment.memberships_lookup["1.1.1"] == ['A', 'G']

def test_add_memberships_lookup_add_double_digit():
    threshold_map = {"level_0": 5.0, "level_1": 3.0, "level_2": 0.0}
    assignment = assign(dist_file=get_path("data/pairwise_distances/double_digit.tsv"),
                        membership_file=get_path("data/clusters/double_digit.tsv"),
                        threshold_map=threshold_map,
                        linkage_method='single',
                        sample_col='id',
                        address_col='address',
                        batch_size=100, delimiter=".")

    # memberships_dict will have:
    # A: 1.1.1
    # B: 2.2.2
    # C: 3.3.3
    # ...
    # N: 14.14.14

    assignment.add_memberships_lookup("O", [15, 15, 15])
    assert assignment.memberships_dict["O"] == "15.15.15"
    assert assignment.memberships_lookup["15"] == ['O']
    assert assignment.memberships_lookup["15.15"] == ['O']
    assert assignment.memberships_lookup["15.15.15"] == ['O']
