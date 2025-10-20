import numpy as np
import scipy
import skbio.tree

class multi_level_clustering:
    """
    Perform hierarchical clustering with multiple threshold levels.

    This class wraps around SciPy's hierarchical clustering tools to:
    - Read a precomputed distance matrix.
    - Generate a linkage matrix.
    - Assign cluster memberships at multiple distance thresholds.
    - Export the resulting tree in Newick format.

    Attributes
    ----------
    cluster_memberships : dict[str, list[str]]
        Mapping from each label to a list of cluster IDs across thresholds.
    thresholds : list of float
        Distance thresholds at which the dendrogram is cut to form clusters.
    labels : list of str
        Labels corresponding to the original observations (tips).
    linkage : np.ndarray
        The hierarchical clustering linkage matrix of shape (n-1, 4).
    newick : str
        The Newick-formatted string representation of the tree.

    Notes
    -----
    - Cluster memberships are stored as strings, one per threshold.
    - The linkage matrix is built using SciPy's `scipy.cluster.hierarchy.linkage`.
    - Newick export relies on scikit-bio's `TreeNode`.
    """

    def __init__(self,dist_mat_file,thresholds,method):
        """
        Initialize the clustering object.

        Parameters
        ----------
        dist_mat_file : str
            Path to a tab-delimited distance matrix file.
        thresholds : list of float
            Distance thresholds at which to assign clusters.
        method : str
            Linkage method to use (e.g., 'single', 'complete', 'average', 'ward').
        """

        #init class attributes
        self.thresholds = thresholds
        self.labels = []
        self.linkage = None
        self.newick = None
        self.cluster_memberships = {}

        #perform clustering
        self.labels, matrix = self.read_distance_matrix(dist_mat_file)  
        self.linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric='precomputed')
        self.init_membership()
        self.assign_clusters()
        self.linkage_to_newick()

    def init_membership(self):
        """
        Initialize the cluster membership dictionary.

        For each label in `self.labels`, create an empty list to store
        its cluster IDs across thresholds.
        """
        self.cluster_memberships = {}

        for label in self.labels:
            self.cluster_memberships[label] = []

    def read_distance_matrix(self,file_path, delim="\t"):
        """
        Read a precomputed distance matrix from file.

        Parameters
        ----------
        file_path : str
            Path to the distance matrix file.
        delim : str, optional (default="\\t")
            Delimiter used in the file (default: tab).

        Returns
        -------
        labels : list of str
            Observation labels parsed from the first column.
        np.ndarray
            Flattened upper-triangular distance values as a 1D NumPy array.

        Notes
        -----
        - The function assumes the first line is a header and skips it.
        - Each subsequent line should start with a label followed by distances.
        - Distances are extracted row by row, omitting redundant lower-triangle entries.
        """
        labels: list[str] = []
        values: list[float] = []

        with open(file_path, 'r',encoding='utf-8') as f:
            next(f, None)  # skip header
            for i, raw in enumerate(f):  # i = 0 for the first data row
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(delim)
                if not parts:
                    continue

                labels.append(parts[0])

                # For row i, skip: 1 (label) + (i + 1) entries up to and including the diagonal
                start = 1 + (i + 1)
                if start < len(parts):
                    try:
                        values.extend(float(x) for x in parts[start:] if x != "")
                    except ValueError as e:
                        raise ValueError(
                            f"Non-numeric value on line {i + 2} (after header): {parts[start:]}"
                        ) from e
        self.validate_distance_matrix(len(labels), len(values))   
        return (labels, np.array(values))

    def validate_distance_matrix(self, num_labels, num_values):
        n = num_labels
        expected = n * (n - 1) // 2
        if num_values != expected:
            raise ValueError(
                f"Expected {expected} upper-triangular distances for n={n}, got {num_values}. "
                "Check file formatting and delimiter."
            )

    def assign_clusters(self):
        """
        Assign cluster memberships for each threshold distance.


        For each threshold in ``self.thresholds``, this method applies 
        hierarchical clustering (using SciPy's ``fcluster``) to the stored 
        linkage matrix. Each label in ``self.labels`` is assigned a cluster 
        identifier at that threshold, and the result is appended to 
        ``self.cluster_memberships``.

        Attributes Used
        ---------------
        self.linkage : np.ndarray
            The hierarchical clustering linkage matrix (n-1, 4).
        self.thresholds : list of float
            A sequence of distance thresholds at which to cut the dendrogram.
        self.labels : list of str
            Labels corresponding to the original observations (tip labels).
        self.cluster_memberships : dict[str, list[str]]
            A mapping of labels to a list of cluster IDs. After execution, 
            each list will contain one entry per threshold.

        Notes
        -----
        - Cluster IDs are integers as returned by ``scipy.cluster.hierarchy.fcluster``.
        - Memberships are stored as strings (converted with ``f'{clusters[cid]}'``).
        - The order of assignments in ``self.cluster_memberships`` corresponds 
        to the order of ``self.thresholds``.
        """
        for dist in self.thresholds:
            clusters = scipy.cluster.hierarchy.fcluster(
                self.linkage, dist, criterion="distance"
            )
            for label, cluster_id in zip(self.labels, clusters):
                self.cluster_memberships[label].append(str(cluster_id))

    def linkage_to_newick(self,cophenetic=True):
        """
        Convert a SciPy linkage matrix into a Newick-formatted tree string.

        This function optionally rescales the distance (cophenetic) values in the
        linkage matrix before building a tree. The tree is constructed using
        scikit-bio's `TreeNode.from_linkage_matrix` and then exported as a Newick string.

        Parameters
        ----------
        cophenetic : bool, optional (default=True)
            If True, multiply the distance values (third column) of the linkage
            matrix by 2 before tree construction. If False, use the original
            distances as-is.

        Attributes Updated
        ------------------
        self.newick : str
            The Newick-formatted string representation of the hierarchical tree.

        Notes
        -----
        - The function assumes `self.linkage` is a valid SciPy linkage matrix of
        shape (n-1, 4).
        - The tip labels are taken from `self.labels`.
        - Single quotes in the Newick string are removed for consistency.
        """
        lmat = self.linkage
        if cophenetic:
            for i in range(lmat.shape[0]):
                lmat[i, 2] *= 2
        
        skb_tree = skbio.tree.TreeNode.from_linkage_matrix(lmat, id_list=self.labels)
        self.newick = str(skb_tree).strip().replace("'", "")



    def get_memberships(self):
        """
        Get the cluster membership mapping.

        Returns
        -------
        dict[str, list[str]]
            Dictionary mapping each label to a list of cluster IDs
            across thresholds.
        """
        return self.cluster_memberships 

