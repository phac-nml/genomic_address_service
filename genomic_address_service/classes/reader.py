from genomic_address_service.constants import EXTENSIONS
from pyarrow.parquet import ParquetFile
import pyarrow as pa 
from genomic_address_service.utils import  get_file_length, get_file_header

class dist_reader:
    record_ids = set()
    dists = {}
    file_handle = None
    row_number = 0

    def __init__(self, f, min_dist=None, max_dist=None, n_records=1000,delim="\t") -> None:
        self.fpath = f
        self.delim = delim
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.n_records = n_records
        self.filter = True
        if min_dist is None and max_dist is None:
            self.filter = False

    
    def guess_file_type(self,f):
        file_type = None
        for t in EXTENSIONS:
            for e in EXTENSIONS[t]:
                if e in f:
                    file_type = t
                    break
        return file_type
    
    def guess_dist_type(self, fpath, ftype, delim="\t"):
        if ftype == 'text':
            header = get_file_header(fpath).split(delim)
            num_rows = get_file_length(fpath)
        elif ftype == 'parquet':
            pf = ParquetFile(fpath) 
            first_rows = next(pf.iter_batches(batch_size = self.n_records))
            df = pa.Table.from_batches([first_rows]).to_pandas() 
            header = list(df.columns)
            num_rows = pf.num_row_groups
        
        if len(header) == 3 and num_rows != 3:
            return 'pd'
        else:
            return 'matrix'
    
    def read_pd(self):
        for line in self.file_handle:
            self.row_number+=1
            line = line.rstrip().split(self.delim)
            if len(line) < 3:
                continue
            qid = line[0]
            rid = line[1]
            d = float(line[2])
            if qid not in self.record_ids and len(self.dists) >= self.n_records:
                self.sort_distances()
                yield self.dists
                self.dists = {}

            if qid not in self.record_ids:  
                self.record_ids.add(qid)
                self.dists[qid] = {}
            if self.filter:
                if self.min_dist is not None:
                    if d < self.min_dist:
                        continue
                if self.max_dist is not None:
                     if d > self.max_dist:
                        continue                   
            self.dists[qid][rid] = d
        self.sort_distances()

    
    def sort_distances(self):
        for qid in self.dists:
            self.dists[qid] = {k: v for k, v in sorted(self.dists[qid].items(), key=lambda item: item[1])}

    def read_matrix(self):
        for line in self.file_handle:
            line = line.rstrip().split(self.delim)
            self.row_number+=1
            qid = self.header[self.row_number]
            if qid not in self.record_ids and len(self.dists) >= self.n_records:
                self.sort_distances()
                yield self.dists
                self.dists = {}
            
            if qid not in self.record_ids:  
                self.record_ids.add(qid)
                self.dists[qid] = {}
                
            values = list(map(float, line[1:]))
            for i in range(0,len(values)):
                rid = self.header[i]
                d = values[i]
                if self.filter:
                    if self.min_dist is not None:
                        if d < self.min_dist:
                            continue
                    if self.max_dist is not None:
                        if d > self.max_dist:
                            continue     
                self.dists[qid][rid] = d
        self.sort_distances()

           
    def read_data(self):
        ftype = self.guess_file_type(self.fpath)
        dist_type = self.guess_dist_type(self.fpath, ftype, self.delim)
        chunk = None
        if ftype == 'text':
            self.file_handle = open(self.fpath,'r')
            self.header = next(self.file_handle).split(self.delim)
        elif ftype == 'parquet':
            self.file_handle = ParquetFile(self.fpath)
        if ftype == 'text' and dist_type == 'pd':
            for chunk in self.read_pd():
                if chunk is not None:
                    yield chunk
            if chunk is None:
                chunk = self.dists
        elif ftype == 'text' and dist_type == 'matrix':
            for chunk in self.read_matrix():
                if chunk is not None:
                    yield chunk
            if chunk is None:
                chunk = self.dists
        yield chunk



    

    
