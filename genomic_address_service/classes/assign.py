import copy
import sys
from statistics import mean
import pandas as pd
from genomic_address_service.constants import EXTENSIONS, PD_HEADER
from genomic_address_service.utils import is_file_ok
from genomic_address_service.classes.reader import dist_reader
class assign:
    avail_methods = ['average','complete','single']
    threshold_map = {}
    thresholds = []
    linkage_method = 'single'
    sample_labels = []
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


    def __init__(self,dist_file,membership_file,threshold_map,linkage_method,address_col, sample_col, batch_size):
        self.dist_file = dist_file
        self.batch_size = batch_size
        file_type = None
        
        if not linkage_method in self.avail_methods:
            self.status = False
            self.error_msgs.append(f'Provided {linkage_method} is not one of the accepted {self.avail_methods}')

        if not is_file_ok(dist_file):
            self.error_msgs.append(f'Provided {dist_file} file does not exist or is empty')
            self.status = False

        if not self.status:
            return

        if is_file_ok(membership_file):
            (_, self.memberships_df) = self.read_data(membership_file)
        else:
            self.error_msgs.append(f'Provided {membership_file} file does not exist or is empty')
            self.status = False
        
        if not self.status:
            return

        self.address_col = address_col
        self.sample_col = sample_col

        #Important sort distances smallest to largest
        self.threshold_map = threshold_map
        self.thresholds = list(self.threshold_map.values())

        columns = self.memberships_df.columns.values.tolist()
        if sample_col not in columns:
            self.status = False
            self.error_msgs.append(f'Could not find sample column: {sample_col} in the file {membership_file}: columns: {columns}')
        if address_col not in columns:
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
         
        self.process_memberships()
        self.ref_labels = set(self.memberships_dict.keys())
        self.init_nomenclature_tracker()
        self.assign(n_records=batch_size)

    def format_df(self,data,delim='.'):
        num_thresholds = len(self.thresholds)
        self.error_samples = []
        membership = {}
        for sample_id in data:
            address = str(data[sample_id]).split(delim)
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

    def get_dist_summary(self,dists):
        min_dist = min(dists)
        ave_dist = mean(dists)
        max_dist = max(dists)
        return {'min':min_dist,'mean':ave_dist,'max':max_dist}

    def get_threshold_idx(self,dist):

        for i in reversed(range(0,len(self.thresholds))):
            if dist <= self.thresholds[i]:
                return i
        return 0
    
    def guess_file_type(self,f):
        file_type = None
        for t in EXTENSIONS:
            for e in EXTENSIONS[t]:
                if e in f:
                    file_type = t
                    break
        return file_type
    
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
    
    def assign(self, n_records=1000,delim="\t"):
        min_dist = min(self.thresholds)
        reader_obj = dist_reader(f=self.dist_file,min_dist=min_dist, max_dist=None, n_records=n_records,delim=delim)
        query_ids = set()
        rank_ids = list(self.nomenclature_cluster_tracker.keys())
        num_ranks = len(self.thresholds)
        for dists in reader_obj.read_data():
            query_ids = query_ids | set(dists.keys())
            for qid in dists:
                is_eligible = False
                self.query_labels.add(qid)
                query_addr = [None] * num_ranks
                if qid in self.memberships_dict:
                    continue

                for rid in dists[qid]:
                    if rid == qid or rid not in self.memberships_dict:
                        continue
                    pairwise_dist = dists[qid][rid]
                    thresh_idx = self.get_threshold_idx(pairwise_dist)
                    thresh_value = self.thresholds[thresh_idx]
                    
                    #save unnecessary work
                    if thresh_value > pairwise_dist:
                        ref_address = self.memberships_dict[rid].split('.')[0:thresh_idx+1]
                        alen = len(ref_address)
                        for i in range(0,len(ref_address)):
                            addr = ".".join(ref_address[0:alen-i])
                            if addr not in self.memberships_lookup:
                                continue
                            addr_members = self.memberships_lookup[addr]
                            addr_dists = []
                            for id in addr_members:
                                addr_dists.append(dists[qid][id])
                            if len(addr_dists) == 0:
                                continue
                            summary = self.get_dist_summary(addr_dists)
                            
                            is_eligible = True
                            if self.linkage_method == 'complete' and summary['max'] > thresh_value:
                                is_eligible = False
                            elif self.linkage_method == 'average' and summary['mean'] > thresh_value:
                                is_eligible = False
                            
                            if is_eligible:
                                for idx,value in enumerate(addr.split('.')):
                                    query_addr[idx] = value
                                break
                        
                    for idx,value in enumerate(query_addr):
                        if value is None:
                            query_addr[idx] = self.nomenclature_cluster_tracker[rank_ids[idx]]
                            self.nomenclature_cluster_tracker[rank_ids[idx]]+=1

                    break
        
                self.memberships_dict[qid] = ".".join([str(x) for x in query_addr])