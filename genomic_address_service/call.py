import os
import sys
import json
from datetime import datetime
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from genomic_address_service.version import __version__
from genomic_address_service.constants import EXTENSIONS, CLUSTER_METHODS, build_call_run_data
from genomic_address_service.utils import is_file_ok, write_threshold_map, write_cluster_assignments, \
init_threshold_map, process_thresholds, has_valid_header_pairwise_distances, has_valid_header_cluster
from genomic_address_service.classes.assign import assign

def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Genomic Address Service: Assignment of samples to existing groupings",
        formatter_class=CustomFormatter)
    parser.add_argument('-d','--dists', type=str, required=True,help='Three column file [query_id,ref_id,dist] in TSV format')
    parser.add_argument('-r', '--rclusters', type=str, required=True, help='Existing cluster file in TSV format')
    parser.add_argument('-m', '--method', type=str, required=False, help='cluster method [single, complete, average]',
                        default='average')
    parser.add_argument('-j', '--thresh_map', type=str, required=False, help='Json file of colname:threshold',
                        default=None)
    parser.add_argument('-s', '--sample_col', type=str, required=False, help='Column name for sample id',
                        default='id')
    parser.add_argument('-c', '--address_col', type=str, required=False, help='Column name for genomic address',
                        default='address')
    parser.add_argument('-t', '--thresholds', type=str, required=False, help='thresholds delimited by , columns will be treated in sequential order')
    parser.add_argument('-o','--outdir', type=str, required=True, help='Output directory to put cluster results')
    parser.add_argument('-l', '--delimiter', type=str, required=False, help='The delimiter used within addresses in the input cluster file, as well as the delimiter to use for addresses in the output. The delimiter must not be a tab or newline character.', default=".")
    parser.add_argument('-b', '--batch_size', type=int, required=False, help='Number of records to process at a time',default=100)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    return parser.parse_args()

def call(config):
    dist_file = config['dists']
    membership_file = config['rclusters']
    thresh_map_file = config['thresh_map']
    outdir = config['outdir']
    linkage_method = config['method']
    thresholds = config['thresholds']
    if thresholds is not None:
        thresholds = process_thresholds(config["thresholds"].split(','))
    delimiter = config['delimiter']
    force = config['force']
    address_col = config['address_col']
    sample_col = config['sample_col']
    run_data = build_call_run_data()
    batch_size = config['batch_size']

    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = config

    if len(delimiter) > 1 or delimiter == "\t" or delimiter == "\n":
        message = f'please specify a different delimiter {delimiter} ie. ,|.|\\||-'
        raise Exception(message)

    if thresholds is None and thresh_map_file is None:
        message = f'you must specify --thresholds or --threshold_map'
        raise Exception(message)

    valid_extensions = list(EXTENSIONS.keys())

    extension = os.path.splitext(dist_file)[1]
    if not extension in valid_extensions:
        message = f'{dist_file} does not have a valid extension {valid_extensions}'
        raise Exception(message)

    extension = os.path.splitext(membership_file)[1]
    if not extension in valid_extensions:
        message = f'{membership_file} does not have a valid extension {valid_extensions}'
        raise Exception(message)

    if not is_file_ok(dist_file):
        message = f'{dist_file} does not exist or is empty'
        raise Exception(message)

    if not is_file_ok(membership_file):
        message = f'{membership_file} does not exist or is empty'
        raise Exception(message)

    if thresh_map_file is not None and not is_file_ok(thresh_map_file ):
        message = f'{thresh_map_file} does not exist or is empty'
        raise Exception(message)

    if not has_valid_header_cluster(membership_file):
        message = f'{membership_file} does not appear to be a properly TSV-formatted file'
        raise Exception(message)

    if not has_valid_header_pairwise_distances(dist_file):
        message = f'{dist_file} does not appear to be a properly TSV-formatted file'
        raise Exception(message)

    if not linkage_method in CLUSTER_METHODS:
        message = f'{linkage_method} is not one of the accepeted methods {CLUSTER_METHODS}'
        raise Exception(message)

    if batch_size < 1:
        message = f'batch size ({batch_size}) must be >=1'
        raise Exception(message)

    if os.path.isdir(outdir) and not force:
        message = f'{outdir} exists, if you would like to overwrite, then specify --force'
        raise Exception(message)

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    if thresh_map_file is not None:
        threshold_map = json.load(open(thresh_map_file,'r'))
    else:
        threshold_map = init_threshold_map(thresholds)


    run_data['threshold_map'] = threshold_map
    write_threshold_map(threshold_map, os.path.join(outdir, "thresholds.json"))

    assignment = assign(dist_file,membership_file,threshold_map,linkage_method,address_col,sample_col,batch_size, delimiter)

    if assignment.status == False:
        exception_message = "something went wrong with cluster assignment"
        exception_message += "\ndistance file: " + str(dist_file)
        exception_message += "\nmembership file: " + str(membership_file)
        exception_message += "\nthreshold map: " + str(threshold_map)
        exception_message += "\nlinkage method: " + str(linkage_method)
        exception_message += "\ndelimiter: " + str(delimiter)
        exception_message += "\n\nCheck error messages:\n" + "\n".join(assignment.error_msgs)

        raise Exception(exception_message)

    cluster_assignments = assignment.memberships_dict

    run_data['result_file'] = os.path.join(outdir, "results.text")

    write_cluster_assignments(run_data['result_file'], cluster_assignments, threshold_map, delimiter, sample_col, address_col)

    with open(os.path.join(outdir,"run.json"),'w') as fh:
        fh.write(json.dumps(run_data, indent=4))

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

def run():
    cmd_args = parse_args()

    try:
        call(vars(cmd_args))

    except Exception as exception:
        print("Exception: " + str(exception))
        sys.exit(1)

# call main function
if __name__ == '__main__':
    run()
