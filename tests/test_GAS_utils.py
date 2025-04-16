import pytest
import os
import json
import pandas as pd
import tempfile
from genomic_address_service.utils import (
    get_file_length, get_file_header, get_file_footer,
    is_matrix_valid, is_file_ok, format_threshold_map,
    write_threshold_map, write_cluster_assignments,
    init_threshold_map
)

def test_get_file_length():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmpfile:
        tmpfile.write("Line 1\nLine 2\nLine 3\n")
    assert get_file_length(tmpfile.name) == 3
    os.unlink(tmpfile.name)

def test_get_file_header():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmpfile:
        tmpfile.write("Header\nLine 1\nLine 2")
    assert get_file_header(tmpfile.name).strip() == "Header"
    os.unlink(tmpfile.name)

def test_get_file_footer():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmpfile:
        tmpfile.write("Line 1\nLine 2\nFooter")
    assert get_file_footer(tmpfile.name).strip() == "Footer"
    os.unlink(tmpfile.name)

def test_is_matrix_valid():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmpfile:
        tmpfile.write("Header1\tHeader2\nValue1\tValue2\nValue3\tValue4")
    assert is_matrix_valid(tmpfile.name) == True
    os.unlink(tmpfile.name)

def test_is_file_ok():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmpfile:
        tmpfile.write("Some content")
    assert is_file_ok(tmpfile.name) == False
    os.unlink(tmpfile.name)

def test_format_threshold_map():
    thresholds = [0.1, 0.2, 0.3]
    expected_output = {'level_1': 0.1, 'level_2': 0.2, 'level_3': 0.3}
    assert format_threshold_map(thresholds) == expected_output

def test_write_threshold_map():
    data = {'level_1': 0.1, 'level_2': 0.2, 'level_3': 0.3}
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmpfile:
        write_threshold_map(data, tmpfile.name)
        tmpfile.seek(0)
        content = json.load(tmpfile)
    assert content == data
    os.unlink(tmpfile.name)

def test_write_cluster_assignments():
    memberships = {'1': '1.2.3', '2': '2.5.6'}
    threshold_map = {'level_1': '10.0', 'level_2': '5.0', 'level_3': '0.0'}
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as tmpfile:
        write_cluster_assignments(tmpfile.name, memberships, threshold_map)

        tmpfile.seek(0)
        df = pd.read_csv(tmpfile.name, sep="\t")
        assert 'id' in df.columns and 'address' in df.columns
        assert df.iloc[0]['address'] == '1.2.3'
        assert df.iloc[1]['address'] == '2.5.6'

    os.unlink(tmpfile.name)

def test_init_threshold_map():
    thresholds = [0.1, 0.2, 0.3]
    result = init_threshold_map(thresholds)
    expected_result = {0: 0.1, 1: 0.2, 2: 0.3}
    assert result == expected_result
    