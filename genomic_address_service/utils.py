import os.path
import shutil
import sys
import time
import psutil
import pandas as pd
import numpy as np
import fastparquet as fp
import tables
from numba import jit
from numba.typed import List
import pyarrow.parquet as pq
import re
import json

from genomic_address_service.constants import MIN_FILE_SIZE


def get_file_length(f):
    return int(os.popen(f'wc -l {f}').read().split()[0])

def get_file_header(f):
    return str(os.popen(f'head -n1 {f}').read())

def get_file_footer(f):
    return str(os.popen(f'tail -n1 {f}').read())

def is_matrix_valid(f):
    num_lines = get_file_length(f)
    footer = get_file_footer(f).split("\t")
    if num_lines == len(footer):
        return True
    return False

def is_file_ok(f):
    status = True
    if not os.path.isfile(f):
        status = False
    elif get_file_length(f) < 2:
        status = False
    elif os.path.getsize(f) < MIN_FILE_SIZE:
        status = False

    return status

def format_threshold_map(thresholds):
    data = {}
    for i,value in enumerate(thresholds):
        data[f'level_{i+1}'] = value
    return data

def write_threshold_map(data,file):
    with open(file,'w') as fh:
        fh.write(json.dumps(data, indent=4))

    fh.close()

def write_cluster_assignments(file ,memberships, threshold_map, delimiter=".", sample_col='id', address_col='address'):
    results = {}
    threshold_keys = list(threshold_map.keys())
    for id in memberships:
        address = memberships[id]
        results[id] = {'id':id,'address':address}
        for idx,value in enumerate(address.split(delimiter)):
            results[id][threshold_keys[idx]] = value
    df = pd.DataFrame.from_dict(results,orient='index')
    df = df[[sample_col,address_col]]
    df.to_csv(file,header=True,sep="\t",index=False)


def init_threshold_map(thresholds):
    thresh_map = {}
    for idx,value in enumerate(thresholds):
        thresh_map[idx] = value

    return thresh_map

def process_thresholds(thresholds):

    try:
        processed = [float(x) for x in thresholds]
    except ValueError:
        message = f'thresholds {thresholds} must all be integers or floats'
        raise Exception(message)

    # Thresholds must be strictly decreasing:
    if not all(processed[i] > processed[i+1] for i in range(len(processed)-1)):
        message = f'thresholds {thresholds} must be in decreasing order'
        raise Exception(message)

    return processed

def has_valid_header_matrix(file_path):
    """
    This file can contain a variable number of delimiters,
    but the minimum should be 1 (2 tokens):

    dists   A
    A   0
    """
    MIN_TOKENS = 2

    with open(file_path) as tsv_file:
        header = tsv_file.readline()

        valid = len(header.split("\t")) >= MIN_TOKENS
        return valid

def has_valid_header_pairwise_distances(file_path):
    """
    This file must contain 2 delimiters (3 tokens):

    query_id    ref_id    dist
    A    A    0
    """
    MIN_TOKENS = 3

    with open(file_path) as tsv_file:
        header = tsv_file.readline()

        valid = len(header.split("\t")) == MIN_TOKENS
        return valid

def has_valid_header_cluster(file_path):
    """
    This file can contain a variable number of delimiters,
    but the minimum should be 1 (2 tokens):

    id    address
    A    1.1.1
    """
    MIN_TOKENS = 2

    with open(file_path) as tsv_file:
        header = tsv_file.readline()

        valid = len(header.split("\t")) >= MIN_TOKENS
        return valid
