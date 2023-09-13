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
import pyarrow.parquet as pq

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

def write_cluster_assignments(file,memberships,threshold_map,outfmt,delimeter="."):
    results = {}
    threshold_keys = list(threshold_map.keys())
    for id in memberships:
        address = memberships[id]
        results[id] = {'id':id,'address':address}
        for idx,value in enumerate(address.split(delimeter)):
            results[id][threshold_keys[idx]] = value
    df = pd.DataFrame.from_dict(results,orient='index')
    if outfmt == 'text':
        df.to_csv(file,header=True,sep="\t",index=False)
    else:
        fp.write(file, df, compression='GZIP')

def init_threshold_map(file,thresholds):
    thresh_map = {}
    df = pd.read_csv(file,header=0,sep="\t")
    cols = df.columns.values.tolist()
    valid_cols = []
    for idx,value in enumerate(cols):
        if value in ['sample_id','id','ST','nomenclature','address']:
            continue
        valid_cols.append(value)

    if len(valid_cols) != len(thresholds):
        return {}

    for idx,value in enumerate(valid_cols):
        thresh_map[value] = thresholds[idx]

    return thresh_map

