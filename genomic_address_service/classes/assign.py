import copy
import sys

import pandas as pd
from genomic_address_service.constants import EXTENSIONS, PD_HEADER
from genomic_address_service.utils import is_file_ok

class assign:
    avail_methods = ['average','complete','single']
    threshold_map = {}
    thresholds = []
    linkage_method = 'single'
    sample_labels = []
    query_df = None
    query_labels = set()
    memberships_df = None
    memberships_dict = {}
    memberships_lookup = {}
    ref_labels = set()
    dist_type = 'matrix'
    error_msgs = []
    status = True
    assignments = {}
    nomenclature_cluster_tracker = {}


    def __init__(self,dist_file,membership_file,threshold_map,linkage_method,address_col, sample_col):
        file_type = None
        
        if not linkage_method in self.avail_methods:
            self.status = False
            self.error_msgs.append(f'Provided {linkage_method} is not one of the accepted {self.avail_methods}')

        if is_file_ok(dist_file):
            (file_type, self.query_df) = self.read_data( dist_file)
            self.dist_type = self.guess_dist_type()
        else:
            self.error_msgs.append(f'Provided {dist_file} file does not exist or is empty')
            self.status = False

        if file_type is None:
            self.status = False

        if not self.status:
            return

        if is_file_ok(membership_file):
            (_, self.memberships_df) = self.read_data(membership_file)
        else:
            self.error_msgs.append(f'Provided {membership_file} file does not exist or is empty')
            self.status = False

        if file_type is None:
            self.status = False

        self.address_col = address_col
        self.sample_col = sample_col

        #Important sort distances smallest to largest
        self.threshold_map = {k: v for k, v in sorted(threshold_map.items(), key=lambda item: item[1])}
        self.thresholds = sorted(list(self.threshold_map.values()))

        columns = self.memberships_df.columns.values.tolist()
        if not sample_col in columns:
            self.status = False
            self.error_msgs.append(f'Could not find sample column: {sample_col} in the file {membership_file}: columns: {columns}')
        if not address_col in columns:
            self.status = False
            self.error_msgs.append(f'Could not find address column: {address_col} in the file {membership_file}: columns: {columns}')            
        if not self.status:
            return

        self.memberships_df = self.memberships_df[[sample_col,address_col]]
        self.memberships_df = self.format_df(self.memberships_df.set_index(sample_col).to_dict()[address_col])
        if len(self.error_samples) > 0:
            self.status = False
            self.error_msgs.append(f'Genomic address too short for samples: {self.error_samples} based on {self.threshold_map}')
        if not self.status:
            return      
         
        self.query_labels = set(self.query_df.columns.values.tolist())
        self.process_memberships()
        self.ref_labels = set(self.memberships_dict.keys())
        self.init_nomenclature_tracker()
        self.assign()

    def format_df(self,data,delim='.'):
        num_thresholds = len(self.thresholds)
        self.error_samples = []
        membership = {}
        for sample_id in data:
            address = data[sample_id].split(delim)
            if len(address) != num_thresholds:
                self.error_samples.append(sample_id)
                continue
            membership[sample_id] = {}
            for idx,value in enumerate(address):
                l = f'level_{idx}'
                try:
                    value = int(value)
                except Exception:
                    self.error_samples.append(sample_id)
                    del membership[sample_id]
                    break                 
                membership[sample_id][l] = int(value)
        return pd.DataFrame.from_dict(membership,orient='index')


    def check_membership_columns(self,cols):
        is_ok = True
        for c in cols:
            if not c in self.threshold_map:
                is_ok = False
                break
        return  is_ok


    def guess_file_type(self,f):
        file_type = None
        for t in EXTENSIONS:
            for e in EXTENSIONS[t]:
                if e in f:
                    file_type = t
                    break
        return file_type

    def guess_dist_type(self):
        cols = self.query_df.columns.values.tolist()
        is_pd_fmt = True
        for c in PD_HEADER:
            if not c in cols:
                is_pd_fmt = False
                break
        if is_pd_fmt:
            return 'pd'
        else:
            return 'matrix'

    def read_data(self,f):
        df = pd.DataFrame()
        file_type = self.guess_file_type(f)
        if file_type == 'text':
            df = pd.read_csv(f, header=0, sep="\t", low_memory=False)
        elif file_type == 'parquet':
            df = pd.read_parquet(
                f,
                engine='auto',
                columns=None,
                storage_options=None,)

        return (file_type, df)

    def init_nomenclature_tracker(self):
        nomenclature_cluster_tracker = self.memberships_df.max().to_frame().T.to_dict()
        for col in nomenclature_cluster_tracker:
            nomenclature_cluster_tracker[col] = nomenclature_cluster_tracker[col][0] + 1
        self.nomenclature_cluster_tracker = nomenclature_cluster_tracker


    def process_memberships(self):
        lookup = {}

        for row in self.memberships_df.itertuples():
            values = [str(x) for x in list(row)]
            id = values[0]
            values = values[1:]
            self.memberships_dict[id] = ".".join([str(x) for x in values])
            for idx,value in enumerate(values):
                code = ".".join(values[0:idx+1])
                if not code in lookup:
                    lookup[code] = list()
                lookup[code].append(id)
        self.memberships_lookup = lookup

    def get_dist_summary(self,q_id,r_id):
        filt_df = self.query_df[self.query_df['query_id'] == q_id]
        address = self.memberships_dict[r_id]
        filt_df = filt_df[filt_df['ref_id'].isin(self.memberships_lookup[address])]
        min_dist = filt_df['dist'].min()
        ave_dist = filt_df['dist'].mean()
        max_dist = filt_df['dist'].max()
        return {'min':min_dist,'mean':ave_dist,'max':max_dist}

    def lookup_col_by_dist(self,d):
        for col in self.threshold_map:
            if self.threshold_map[col] == d:
                return col
        return None


    def assign(self):
        subset = self.query_df
        unassigned_ids = set(subset['query_id'].unique())
        num_ranks = len(self.thresholds)
        rank_ids = list(self.nomenclature_cluster_tracker.keys())

        for idx,thresh in enumerate(self.thresholds):
            if len(unassigned_ids) == 0:
                break
            subset = subset[subset['query_id'].isin(list(unassigned_ids))]
            filt_df = subset[subset['dist'] <= thresh].sort_values(by=['dist'])
            query_ids = filt_df['query_id'].unique()
            for q_id in query_ids:
                a = [None] * num_ranks
                if q_id in self.memberships_dict:
                    unassigned_ids = unassigned_ids - set(q_id)
                    continue
                tmp_df = filt_df[filt_df['query_id'] == q_id]
                ref_ids = tmp_df['ref_id'].unique()
                if len(ref_ids) == 0:
                    continue
                checked_addresses = set()

                for r_id in ref_ids:
                    if r_id == q_id or not r_id in self.memberships_dict:
                        continue
                    a = copy.deepcopy(self.memberships_dict[r_id])
                    if a in checked_addresses:
                        continue
                    summary = self.get_dist_summary(q_id, r_id)
                    checked_addresses.add(a)
                    is_eligible = True
                    if self.linkage_method == 'complete' and summary['max'] > thresh:
                        is_eligible = False
                    elif self.linkage_method == 'average' and summary['mean'] > thresh:
                        is_eligible = False
                    if is_eligible:
                        break

                #Sample is too distant to be assigned at this threshold
                if r_id == q_id:
                    continue

                unassigned_ids = unassigned_ids - set([q_id])

                #remove the last n codes from the address based on the threshold
                #Pad the code out with None
                if isinstance(a,str):
                    a = a.split('.')

                #Blank the digits of the code which are past the rank
                a.reverse()
                for i in range(0,idx):
                    a[i] = None
                a.reverse()

                for i,value in enumerate(a):
                    if value is not None:
                        continue
                    a[i] = self.nomenclature_cluster_tracker[rank_ids[i]]
                    self.nomenclature_cluster_tracker[rank_ids[i]]+=1
                self.memberships_dict[q_id] = ".".join([str(x) for x in a])

                for i, value in enumerate(a):
                    code = ".".join([str(x) for x in a[0:i + 1]])
                    if not code in self.memberships_lookup:
                        self.memberships_lookup[code] = list()
                    self.memberships_lookup[code].append(q_id)


        for q_id in unassigned_ids:
            a = [None] * num_ranks
            for i, value in enumerate(a):
                if value is not None:
                    continue
                a[i] = self.nomenclature_cluster_tracker[rank_ids[i]]
                self.nomenclature_cluster_tracker[rank_ids[i]] += 1
            self.memberships_dict[q_id] = ".".join([str(x) for x in a])