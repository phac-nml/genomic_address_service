import pandas as pd
import scipy


class multi_level_clustering:
    cluster_memberships = {}
    thresholds = []
    labels = []
    linkage = None
    newick = None


    def __init__(self,dist_mat_file,thresholds,method):
        df = self.read_matrix(dist_mat_file).astype(float)
        self.labels = df.columns.values.tolist()
        matrix = scipy.spatial.distance.squareform(df.values)
        del df
        self.thresholds = thresholds
        self.linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric='precomputed')
        self.init_membership()
        self.assign_clusters()
        self.linkage_to_newick()

    def init_membership(self):
        for label in self.labels:
            self.cluster_memberships[label] = []

    def read_matrix(self,f):
        return pd.read_csv(f,header=0,sep="\t",index_col=0)

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
                newick = f"):{(parentdist - node.dist) / 2}{newick}"
            else:
                newick = ");"
            newick = self.buildNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = self.buildNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick

    def get_memberships(self):
        return self.cluster_memberships




