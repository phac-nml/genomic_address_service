class dist_reader:

    def __init__(self, f, n_records=1000, delim="\t") -> None:
        self.record_ids = set()
        self.dists = {}
        self.file_handle = None
        self.row_number = 0
        self.fpath = f
        self.delim = delim
        self.n_records = n_records

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
            self.dists[qid][rid] = d
        self.sort_distances()
       
        yield self.dists


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
                self.dists[qid][rid] = d
        self.sort_distances()

    def read_data(self):
        self.file_handle = open(self.fpath,'r')
        self.header = next(self.file_handle).split(self.delim)

        for chunk in self.read_pd():
            if chunk is not None:
                yield chunk
        if chunk is None:
            chunk = self.dists

        self.file_handle.close()
        return chunk
