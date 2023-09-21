import os
import sys
import json
from datetime import datetime
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from genomic_address_service.version import __version__
from genomic_address_service.constants import CLUSTER_METHODS, CALL_RUN_DATA
from genomic_address_service.utils import is_file_ok, write_threshold_map, write_cluster_assignments, \
init_threshold_map
from genomic_address_service.classes.assign import assign


def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Genomic Address Service: Assignment of samples to existing groupings",
        formatter_class=CustomFormatter)
    parser.add_argument('-d','--dists', type=str, required=True,help='Three column file [query_id,ref_id,dist] in TSV or parquet format')
    parser.add_argument('-r', '--rclusters', type=str, required=True, help='Existing cluster file in TSV or parquet format')
    parser.add_argument('-m', '--method', type=str, required=False, help='cluster method [single, complete, average]',
                        default='average')
    parser.add_argument('-j', '--thresh_map', type=str, required=False, help='Json file of colname:threshold',
                        default=None)
    parser.add_argument('-t', '--thresholds', type=str, required=False, help='thresholds delimited by , columns will be treated in sequential order')
    parser.add_argument('-o','--outdir', type=str, required=True, help='Output directory to put cluster results')
    parser.add_argument('-u', '--outfmt', type=str, required=False, help='Output format for assignments [text, parquet]',default='text')
    parser.add_argument('-l', '--delimeter', type=str, required=False, help='delimeter desired for nomenclature code',default=".")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')

    return parser.parse_args()



def run():
    cmd_args = parse_args()
    dist_file = cmd_args.dists
    membership_file = cmd_args.rclusters
    thresh_map_file = cmd_args.thresh_map
    outdir = cmd_args.outdir
    linkage_method = cmd_args.method
    thresholds = cmd_args.thresholds
    if thresholds is not None:
        thresholds = [float(x) for x in cmd_args.thresholds.split(',')]
    delimeter = cmd_args.delimeter
    force = cmd_args.force
    outfmt = cmd_args.outfmt

    run_data = CALL_RUN_DATA

    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    if outfmt not in ['text','parquet']:
        print(f'Error please specify a either text or parquet as the output format you specified: {outfmt}')
        sys.exit()


    if len(delimeter) > 1 or delimeter == "\t" or delimeter == "\n":
        print(f'Error please specify a different delimeter {delimeter} ie. ,|.|\||-')
        sys.exit()

    if thresholds is None and thresh_map_file is None:
        print(f'Error you must specify --thresholds or --threshold_map')
        sys.exit()

    if not is_file_ok(dist_file):
        print(f'Error {dist_file} does not exist or is empty')
        sys.exit()

    if not is_file_ok(membership_file ):
        print(f'Error {membership_file } does not exist or is empty')
        sys.exit()

    if thresh_map_file is not None and not is_file_ok(thresh_map_file ):
        print(f'Error {thresh_map_file } does not exist or is empty')
        sys.exit()

    if not linkage_method in CLUSTER_METHODS:
        print(f'Error {linkage_method} is not one of the accepeted methods {CLUSTER_METHODS}')
        sys.exit()


    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    if thresh_map_file is not None:
        threshold_map = json.load(open(thresh_map_file,'r'))
    else:
        threshold_map = init_threshold_map(membership_file,thresholds)

    if len(threshold_map) == 0 and thresh_map_file is not None:
        print(f'Error {thresh_map_file } is malformed or empty')
        sys.exit()
    elif  len(threshold_map) == 0 :
        print(f'Error {thresholds }  supplied thresholds do not equal the number of threshold columns in {membership_file}')
        sys.exit()

    run_data['threshold_map'] = threshold_map
    write_threshold_map(threshold_map, os.path.join(outdir, "thresholds.json"))

    obj = assign(dist_file,membership_file,threshold_map,linkage_method)

    if obj.status == False:
        print(f'Error something went wrong with cluster assignment. check error messages {obj.error_msgs}')
        sys.exit()


    cluster_assignments = obj.memberships_dict

    run_data['result_file'] = os.path.join(outdir, "results.{}".format(outfmt))

    write_cluster_assignments(run_data['result_file'], cluster_assignments, threshold_map, outfmt,delimeter)

    with open(os.path.join(outdir,"run.json"),'w') as fh:
        fh.write(json.dumps(run_data, indent=4))

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

# call main function
if __name__ == '__main__':
    run()
