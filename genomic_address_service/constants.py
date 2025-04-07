from genomic_address_service.version import __version__

TEXT = 'text'
EXTENSIONS = {}
EXTENSIONS.update(dict.fromkeys(['.txt', '.tsv', '.mat', '.text'], TEXT))
# '.txt' -> TEXT, '.tsv' -> TEXT, etc.

PD_HEADER = [
    'query_id',
    'ref_id',
    'dist'
]

MIN_FILE_SIZE = 32
CLUSTER_METHODS = ['average','complete','single']

def build_mc_run_data():
    run_data = {
        'genomic address service: de novo clustering': f'version: {__version__}',
        'analysis_start_time':'',
        'analysis_end_time':'',
        'parameters':{},
        'threshold_map':{},
        'result_file':''
    }

    return run_data

def build_call_run_data():
    run_data = {
        'genomic address service: cluster assignment': f'version: {__version__}',
        'analysis_start_time':'',
        'analysis_end_time':'',
        'parameters':{},
        'threshold_map':{},
        'result_file':''
    }

    return run_data
