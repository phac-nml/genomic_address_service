import numpy as np
import scipy

class multi_level_clustering:

    def __init__(self,dist_mat_file,thresholds,method):
        self.labels, matrix = self.read_distance_matrix(dist_mat_file)
        self.thresholds = thresholds
        self.linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric='precomputed')
        self.init_membership()
        self.assign_clusters()
        self.linkage_to_newick()

    def init_membership(self):
        self.cluster_memberships = {}

        for label in self.labels:
            self.cluster_memberships[label] = []

    def read_distance_matrix(self,file_path, delim="\t"):
        labels = []
        v = []
        with open(file_path, 'r') as f:
            header = next(f)
            offset = 2
            for line in f:
                line_split = line.strip().split(delim)
                if len(line_split) < 1:
                    continue
                labels.append(line_split[0])
                v += [float(x) for x in line_split[offset:]]
                offset += 1
        return (labels, np.array(v))

    def assign_clusters(self):
        for idx,dist in enumerate(self.thresholds):
            cid = 0
            clusters = scipy.cluster.hierarchy.fcluster(self.linkage, dist, criterion='distance')
            for label in self.labels:
                self.cluster_memberships[label].append(f'{clusters[cid]}')
                cid+=1

    def linkage_to_newick(self):
        tree = scipy.cluster.hierarchy.to_tree(self.linkage)
        self.newick = self.buildNewick(tree, "", tree.dist, self.labels)

    def buildNewick(self,node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = f"):{parentdist - node.dist}{newick}"
            else:
                newick = ");"
            newick = self.buildNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = self.buildNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick

    def get_memberships(self):
        return self.cluster_memberships
