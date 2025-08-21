import copy
import sys
import os
from statistics import mean
import pandas as pd
from genomic_address_service.constants import EXTENSIONS, TEXT
from genomic_address_service.utils import is_file_ok
from genomic_address_service.classes.reader import dist_reader

class assign:
    ERROR_MISSING_DELIMITER = "delimiter was not found"
    ERROR_LENGTH = "genomic address length is incorrect"
    ERROR_NON_INTEGER = "address could not be converted to an integer"

    AVAILABLE_METHODS = ["average", "complete", "single"]

    def __init__(self,dist_file,membership_file,threshold_map,linkage_method,address_col, sample_col, batch_size, delimiter):
        self.dist_file = dist_file
        self.batch_size = batch_size
        file_type = None
        self.threshold_map = threshold_map
        self.thresholds = []
        self.linkage_method = linkage_method
        self.sample_labels = []
        self.query_labels = set()
        self.memberships_df = None
        self.memberships_dict = {}
        self.memberships_lookup = {}
        self.ref_labels = set()
        self.dist_type = 'matrix'
        self.error_msgs = []
        self.status = True
        self.assignments = {}
        self.nomenclature_cluster_tracker = {}
        self.query_ids = set()
        self.delimiter = delimiter

        self.error_samples = {
            self.ERROR_MISSING_DELIMITER: [],
            self.ERROR_LENGTH: [],
            self.ERROR_NON_INTEGER: []
        } # message -> list of IDs

        if not linkage_method in self.AVAILABLE_METHODS:
            self.status = False
            self.error_msgs.append(f'Provided {linkage_method} is not one of the accepted {self.AVAILABLE_METHODS}')

        if not is_file_ok(dist_file):
            self.error_msgs.append(f'Provided {dist_file} file does not exist or is empty')
            self.status = False

        if not self.status:
            return

        if is_file_ok(membership_file):
            self.memberships_df = self.read_data(membership_file)
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
        self.memberships_df = self.format_df(self.memberships_df.set_index(sample_col).to_dict()[address_col], self.delimiter)

        if len(self.error_samples[self.ERROR_MISSING_DELIMITER]) > 0:
            self.status = False
            self.error_msgs.append(f'Error: {self.ERROR_MISSING_DELIMITER} for samples {self.error_samples[self.ERROR_MISSING_DELIMITER]}.')

        if len(self.error_samples[self.ERROR_LENGTH]) > 0:
            self.status = False
            self.error_msgs.append(f'Error: {self.ERROR_LENGTH} for samples {self.error_samples[self.ERROR_LENGTH]}; expected length ({len(self.thresholds)}) based on thresholds {self.threshold_map}.')

        if len(self.error_samples[self.ERROR_NON_INTEGER]) > 0:
            self.status = False
            self.error_msgs.append(f'Error: {self.ERROR_NON_INTEGER} for samples {self.error_samples[self.ERROR_NON_INTEGER]}.')

        if not self.status:
            return

        self.process_memberships()
        self.ref_labels = set(self.memberships_dict.keys())
        self.init_nomenclature_tracker()
        self.assign(n_records=batch_size)

    def format_df(self,data, delim='.'):
        num_thresholds = len(self.thresholds)
        membership = {}

        for sample_id in data:
            address = str(data[sample_id]).split(delim)

            # Did we get the wrong number of address components?
            if len(address) != num_thresholds:

                # No delimiter:
                if (not delim in data[sample_id]) and num_thresholds > 1:
                    self.error_samples[self.ERROR_MISSING_DELIMITER].append(sample_id)
                # Length problem:
                else:
                    self.error_samples[self.ERROR_LENGTH].append(sample_id)

                continue
            membership[sample_id] = {}
            for idx,value in enumerate(address):
                l = f'level_{idx}'
                try:
                    value = int(value)
                except Exception:
                    self.error_samples[self.ERROR_NON_INTEGER].append(sample_id)
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
        self.nomenclature_cluster_tracker = self.memberships_df.max().to_frame().T.to_dict()
        for col in self.nomenclature_cluster_tracker:
            self.nomenclature_cluster_tracker[col] = self.nomenclature_cluster_tracker[col][0] + 1

    def process_memberships(self):
        lookup = {}
        for row in self.memberships_df.itertuples():
            values = [str(x) for x in list(row)]
            id = values[0]
            values = values[1:]
            self.memberships_dict[id] = self.delimiter.join([str(x) for x in values])
            for idx,value in enumerate(values):
                code = self.delimiter.join(values[0:idx+1])
                if not code in lookup:
                    lookup[code] = list()
                lookup[code].append(id)
        self.memberships_lookup = lookup

    def add_memberships_lookup(self,sample_id, address):
        self.memberships_dict[sample_id] = self.delimiter.join([str(x) for x in address])
        for idx in range(0,len(address)):
            code = self.delimiter.join([str(x) for x in address[0:idx+1]])
            if not code in self.memberships_lookup:
                self.memberships_lookup[code] = list()
            self.memberships_lookup[code].append(sample_id)


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

    def check_file_type(self,f):
        extension = os.path.splitext(f)[1]
        valid_extensions = list(EXTENSIONS.keys())

        if not extension in valid_extensions:
            message = f'{f} does not have a valid extension ({extension}): {valid_extensions}'
            raise Exception(message)

    def read_data(self, f):
        self.check_file_type(f)
        df = pd.read_csv(f, header=0, sep="\t", low_memory=False)

        return df

    def assign(self, n_records=1000,delim="\t"):
        reader_obj = dist_reader(f=self.dist_file, n_records=n_records, delim=delim)
        self.query_ids = set()
        rank_ids = list(self.nomenclature_cluster_tracker.keys())
        num_ranks = len(self.thresholds)
        for dists in reader_obj.read_data():
            self.query_ids = self.query_ids | set(dists.keys())
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
                    if thresh_value >= pairwise_dist:
                        ref_address = self.memberships_dict[rid].split(self.delimiter)[0:thresh_idx+1]
                        alen = len(ref_address)
                        for i in range(0,len(ref_address)):
                            addr = self.delimiter.join(ref_address[0:alen-i])
                            
                            if addr not in self.memberships_lookup:
                                continue
                            addr_members = self.memberships_lookup[addr]
                            addr_dists = []
                            for id in addr_members:
                                if id in dists[qid]:
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
                                for idx,value in enumerate(addr.split(self.delimiter)):
                                    query_addr[idx] = value
                                break
                            thresh_value = self.thresholds[thresh_idx-(i+1)]
                            
                    for idx,value in enumerate(query_addr):
                        if value is None:
                            query_addr[idx] = self.nomenclature_cluster_tracker[rank_ids[idx]]
                            self.nomenclature_cluster_tracker[rank_ids[idx]]+=1
                    break

                self.add_memberships_lookup(qid, query_addr)
