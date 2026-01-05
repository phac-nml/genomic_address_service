import numpy as np
import pandas as pd
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
    VALID_TREE_DISTANCES = ['patristic', 'cophenetic']

    def __init__(self, dist_mat_file, thresholds, method, sort_matrix, tree_distances='patristic'):
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
        tree_distances : string, optional (default=patristic)
            Defines how distances in the original distance matrix are represented in the tree (Newick file).
            If pastristic, the distances in the distance matrix correspond to the patristic distances between nodes or leaves in the tree.
            If 'cophenetic', the distances in the distance matrix correspond to the cophenetic distance (height in tree where leaves first
            share a common ancestor). For ultrametic trees this corresponds to twice the patristic distance.
        """

        #init class attributes
        self.thresholds = thresholds
        self.labels = []
        self.linkage = None
        self.newick = None
        self.cluster_memberships = {}

        #perform clustering
        self.labels, matrix = self.read_distance_matrix(dist_mat_file, sort_matrix=sort_matrix)
        self.linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric='precomputed')
        self._init_membership()
        self._assign_clusters()
        self._linkage_to_newick(tree_distances=tree_distances)

    def _init_membership(self):
        """
        Initialize the cluster membership dictionary.

        For each label in `self.labels`, create an empty list to store
        its cluster IDs across thresholds.
        """
        self.cluster_memberships = {}

        for label in self.labels:
            self.cluster_memberships[label] = []

    def read_distance_matrix(self,file_path, delim="\t", sort_matrix=False):
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

        # Read the distance matrix into a DataFrame
        df = pd.read_csv(file_path, header=0, sep=delim, low_memory=False, dtype=str)

        # There's a bug in some versions of pandas with the read_csv function.
        # If you attempt pd.read_csv(dtype=str, index_col=0),
        # then pandas will not cast the index to the specified type (str).
        # We work around this by not loading an index, and then reshaping the
        # DataFrame to have the index we want, which will be in str format.
        index = df.iloc[:, 0]
        df = df.iloc[:, 1:]
        df = df.set_index(index)
        # Convert all values to float, raising an error if conversion fails
        try:
            df = df.astype(float)
        except ValueError as e:
            raise ValueError("Input matrix should only contain numerical values") from e

        # Check that there are no NaN values in the distance matrix
        if np.isnan(df.values).any():
            raise ValueError("Distance matrix contains NaN, null or NA values.")
        # Check that the distance matrix is square and rows/columns match
        if not df.index.equals(df.columns):
            raise ValueError("Incorrect Distance Matrix Format: --matrix must have (n x n) dimensions, 0 diagonal starting at position [0,0] and rows/columns must in the same order.")
        if sort_matrix:
            df = df.sort_index(ascending=True)
            df = df.sort_index(axis=1, ascending=True)

        # Split matrix into lower and upper triangular matrices
        upper_mask = np.triu(np.ones_like(df, dtype=bool))
        lower_mask = np.tril(np.ones_like(df, dtype=bool))
        
        lower_tri = df.mask(upper_mask) # Mask upper triangle
        upper_tri = df.mask(lower_mask) # Mask lower triangle

        one_dim_tri_lower = lower_tri.values.flatten(order='F') # Reads elements column-wise
        one_dim_tri_upper = upper_tri.values.flatten(order='C') # Reads elements row-wise

        lower_array = one_dim_tri_lower[~np.isnan(one_dim_tri_lower)]
        upper_array = one_dim_tri_upper[~np.isnan(one_dim_tri_upper)]

        # Validate symmetry of distance matrix values. Upper and lower triangles must match.
        if not np.array_equal(lower_array, upper_array):
            raise ValueError("Distance matrix has non-symmetrical values")

        # Extract non-NaN values from one triangle (lower)
        values.extend(lower_array)
        labels.extend(df.index.tolist())
        
        return (labels, np.array(values))


    def _assign_clusters(self):
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

    def _linkage_to_newick(self, tree_distances):
        """
        Convert a SciPy linkage matrix into a Newick-formatted tree string.

        This function optionally rescales the distance (cophenetic) values in the
        linkage matrix before building a tree. The tree is constructed using
        scikit-bio's `TreeNode.from_linkage_matrix` and then exported as a Newick string.

        Parameters
        ----------
        tree_distances : string, optional (default=patristic)
            Defines how distances in the original distance matrix are represented in the tree (Newick file).
            If pastristic, the distances in the distance matrix correspond to the patristic distances between nodes or leaves in the tree.
            If 'cophenetic', the distances in the distance matrix correspond to the cophenetic distance (height in tree where leaves first
            share a common ancestor). For ultrametic trees this corresponds to twice the patristic distance 
            (hence we multiply distances values (third column) in the linkage matrix by 2).

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
        lmat = None
        if tree_distances == 'patristic':
            lmat = self.linkage
        elif tree_distances == 'cophenetic':
            lmat = self.linkage
            for i in range(lmat.shape[0]):
                lmat[i, 2] *= 2
        else:
            raise Exception(f'Invalid tree_distances value [{tree_distances}]. Must be one of {self.VALID_TREE_DISTANCES}')
        
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

