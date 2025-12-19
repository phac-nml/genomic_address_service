import pytest
import pathlib
import csv
import json
import skbio
from skbio.tree import TreeNode
from os import path
from io import StringIO
import pandas as pd

from genomic_address_service.mcluster import mcluster

def get_path(location):
    directory = path.dirname(path.abspath(__file__))
    return path.join(directory, location)

def distance_patristic_from_tree(tree, a, b) -> float:
    """
    Calculates the patristic distance from the tree.
    That is the sum of branch lengths between leaves 'a' and 'b'.
    """
    return tree.find(a).distance(tree.find(b))

def distance_cophenetic_from_tree(tree, a, b) -> float:
    """
    Calculates the cophenetic distance from the tree.
    That is the height of the lowest common ancestor of leaves 'a' and 'b'.
    """
    height, node = tree.lca([a, b]).height(use_length=True)
    return height

def distance_from_matrix(matrix, a, b) -> float:
    return (float)(matrix[a][b])

def test_basic(tmp_path):
    # A basic example with one threshold (0) where every
    # item should only cluster with itself or others with
    # a distance of 0.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # A, B, I, J unique
        # C, D share
        # E, F share
        # G, H share
        assert ["id", "address", "level_1"] in clusters
        assert ["A", "3", "3"] in clusters
        assert ["B", "4", "4"] in clusters
        assert ["C", "2", "2"] in clusters
        assert ["D", "2", "2"] in clusters
        assert ["E", "1", "1"] in clusters
        assert ["F", "1", "1"] in clusters
        assert ["G", "5", "5"] in clusters
        assert ["H", "5", "5"] in clusters
        assert ["I", "6", "6"] in clusters
        assert ["J", "7", "7"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "0"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 1
        assert run_json["threshold_map"]["level_1"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 1
        assert thresholds_json["level_1"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("(((J:1.000000,I:1.000000):2.0,(H:0.000000,G:0.000000):3.0):3.0,(((B:1.000000,A:1.000000):1.0,(D:0.000000,C:0.000000):2.0):1.0,(F:0.000000,E:0.000000):3.0):3.0);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(((E:0.0,F:0.0):3.0,((C:0.0,D:0.0):2.0,(A:1.0,B:1.0):1.0):1.0):3.0,((G:0.0,H:0.0):3.0,(I:1.0,J:1.0):2.0):3.0);"

def test_basic_patristic(tmp_path):
    # A basic example with one threshold (0) where every
    # item should only cluster with itself or others with
    # a distance of 0.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'patristic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # A, B, I, J unique
        # C, D share
        # E, F share
        # G, H share
        assert ["id", "address", "level_1"] in clusters
        assert ["A", "3", "3"] in clusters
        assert ["B", "4", "4"] in clusters
        assert ["C", "2", "2"] in clusters
        assert ["D", "2", "2"] in clusters
        assert ["E", "1", "1"] in clusters
        assert ["F", "1", "1"] in clusters
        assert ["G", "5", "5"] in clusters
        assert ["H", "5", "5"] in clusters
        assert ["I", "6", "6"] in clusters
        assert ["J", "7", "7"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "0"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'patristic'

        assert len(run_json["threshold_map"]) == 1
        assert run_json["threshold_map"]["level_1"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 1
        assert thresholds_json["level_1"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    assert str(actual_tree).strip() == "(((E:0.0,F:0.0):1.5,((C:0.0,D:0.0):1.0,(A:0.5,B:0.5):0.5):0.5):1.5,((G:0.0,H:0.0):1.5,(I:0.5,J:0.5):1.0):1.5);"

def test_compare_tree_patristic_cophenetic(tmp_path):
    args_patristic = {
        "matrix": get_path("data/matrix/basic.tsv"),
        "outdir": path.join(tmp_path, "test_out"),
        "method": "single",
        "thresholds": "0",
        "sort_matrix": False,
        "delimiter": ".",
        "force": False,
        "tree_distances": 'patristic'
    }

    args_cophenetic = {
        "matrix": get_path("data/matrix/basic.tsv"),
        "outdir": path.join(tmp_path, "test_out2"),
        "method": "single",
        "thresholds": "0",
        "sort_matrix": False,
        "delimiter": ".",
        "force": False,
        "tree_distances": 'cophenetic'
    }

    mcluster(args_patristic)
    mcluster(args_cophenetic)

    input_distance_matrix = pd.read_csv(get_path("data/matrix/basic.tsv"), sep='\t', index_col='dists')

    tree_path_patristic = path.join(args_patristic["outdir"], "tree.nwk")
    actual_tree_patristic = skbio.io.registry.read(tree_path_patristic, format='newick', into=TreeNode)

    tree_path_cophenetic = path.join(args_cophenetic["outdir"], "tree.nwk")
    actual_tree_cophenetic = skbio.io.registry.read(tree_path_cophenetic, format='newick', into=TreeNode)

    # The difference in Newick file between patristic and cophenetic is a factor of 2 on branch lengths
    assert str(actual_tree_patristic).strip()  == "(((E:0.0,F:0.0):1.5,((C:0.0,D:0.0):1.0,(A:0.5,B:0.5):0.5):0.5):1.5,((G:0.0,H:0.0):1.5,(I:0.5,J:0.5):1.0):1.5);"
    assert str(actual_tree_cophenetic).strip() == "(((E:0.0,F:0.0):3.0,((C:0.0,D:0.0):2.0,(A:1.0,B:1.0):1.0):1.0):3.0,((G:0.0,H:0.0):3.0,(I:1.0,J:1.0):2.0):3.0);"

    # For context for understanding tests, here is a dendrogram
    # Generated from newick tree string above and using https://github.com/JLSteenwyk/PhyKIT: "phykit print_tree tree.txt"
    #
    #                                                                            , E
    #                                        ____________________________________|
    #                                       |                                    | F
    #                                       |
    #   ____________________________________|                                    , C
    #  |                                    |            ________________________|
    #  |                                    |           |                        | D
    #  |                                    |___________|
    #  |                                                |            ____________ A
    # _|                                                |___________|
    #  |                                                            |____________ B
    #  |
    #  |                                                                         , G
    #  |                                     ____________________________________|
    #  |                                    |                                    | H
    #  |____________________________________|
    #                                       |                        ____________ I
    #                                       |_______________________|
    #                                                               |____________ J
    #
    # And here is the input distance matrix ("data/matrix/basic.tsv")
    # dists   A       B       C       D       E       F       G       H       I       J
    # A       0       1       2       2       5       5       6       6       9       9
    # B       1       0       3       3       6       6       7       7       9       9
    # C       2       3       0       0       3       3       6       6       9       9
    # D       2       3       0       0       3       3       6       6       9       9
    # E       5       6       3       3       0       0       6       6       9       9
    # F       5       6       3       3       0       0       6       6       9       9
    # G       6       7       6       6       6       6       0       0       3       3
    # H       6       7       6       6       6       6       0       0       3       3
    # I       9       9       9       9       9       9       3       3       0       1
    # J       9       9       9       9       9       9       3       3       1       0

    # For the newick tree using "--tree-distances patristic", the patristic distance (sum of branch lengths)
    # between two leaves 'A' and 'B' corresponds to the distance value from the input distance matrix between 'A' and 'B'
    assert distance_patristic_from_tree( actual_tree_patristic,  'A', 'B') == 1.0
    assert distance_from_matrix(         input_distance_matrix,  'A', 'B') == 1.0

    # In this scenario, calculating the cophenetic distance between two leaves 'A' and 'B' (height of lowest common ancestor)
    # Corresponds to a distance of 1/2 the value from the input matrix
    # Hence why the value of the parameter "--tree-distances patristic" is called "patristic"
    assert distance_cophenetic_from_tree(actual_tree_patristic,  'A', 'B') == 0.5
    assert distance_from_matrix(         input_distance_matrix,  'A', 'B') == 1.0

    # For the newick tree using "--tree-distances cophenetic", the patristic distance (sum of branch lengths)
    # between two leaves 'A' and 'B' corresponds to twice the distance value from the input distance matrix between 'A' and 'B'
    assert distance_patristic_from_tree( actual_tree_cophenetic, 'A', 'B') == 2.0
    assert distance_from_matrix(         input_distance_matrix,  'A', 'B') == 1.0
    # However, this means that the cophenetic distance between two leaves 'A' and 'B' (height of lowest common ancestor)
    # corresponds to the distance value from the input distance matrix
    # Hence why the value of the parameter "--tree-distances cophenetic" is called "cophenetic"
    assert distance_cophenetic_from_tree(actual_tree_cophenetic, 'A', 'B') == 1.0
    assert distance_from_matrix(         input_distance_matrix,  'A', 'B') == 1.0


    # Repeats same comparisons as above, but between leaves 'A' and 'C' (which have distance 0 between each other, so all distances are 0)
    # Case patristic tree_distances ("--tree-distances patristic")
    assert distance_patristic_from_tree( actual_tree_patristic,  'A', 'C') == 2.0
    assert distance_from_matrix(         input_distance_matrix,  'A', 'C') == 2.0
    assert distance_cophenetic_from_tree(actual_tree_patristic,  'A', 'C') == 1.0

    # Case cophenetic tree_distances ("--tree-distances cophenetic")
    assert distance_patristic_from_tree( actual_tree_cophenetic, 'A', 'C') == 4.0
    assert distance_from_matrix(         input_distance_matrix,  'A', 'C') == 2.0
    assert distance_cophenetic_from_tree(actual_tree_cophenetic, 'A', 'C') == 2.0


    # Repeats same comparisons as above, but between leaves 'C' and 'D' (which have distance 0 between each other, so all distances are 0)
    # Case patristic tree_distances ("--tree-distances patristic")
    assert distance_patristic_from_tree( actual_tree_patristic,  'C', 'D') == 0.0
    assert distance_from_matrix(         input_distance_matrix,  'C', 'D') == 0.0
    assert distance_cophenetic_from_tree(actual_tree_patristic,  'C', 'D') == 0.0

    # Case cophenetic tree_distances ("--tree-distances cophenetic")
    assert distance_patristic_from_tree( actual_tree_cophenetic, 'C', 'D') == 0.0
    assert distance_from_matrix(         input_distance_matrix,  'C', 'D') == 0.0
    assert distance_cophenetic_from_tree(actual_tree_cophenetic, 'C', 'D') == 0.0


    # Note: Distances from a tree (patristic or cophenetic) don't always correspond to the original input matrix
    assert distance_patristic_from_tree(actual_tree_patristic, 'A', 'J') == 6.0
    assert distance_from_matrix(        input_distance_matrix, 'A', 'J') == 9.0

def test_wikipedia(tmp_path):
    # Ensures mcluster generates the same output as this
    # example on Wikipedia (2025-03-28):
    # https://en.wikipedia.org/wiki/Single-linkage_clustering
    # This mostly tests the linkage part.
    # Note that the article labels the tree with patristic
    # distances, so the thresholds used here are doubled,
    # with respect to that tree.
    args = {"matrix": get_path("data/matrix/wikipedia-single.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "25,18,0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    # t=25
    # a,b,c,e
    # d

    # t=18
    # a,b
    # c
    # e
    # d

    # t=0
    # a
    # b
    # c
    # d
    # e

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address", "level_1", "level_2", "level_3"] in clusters
        assert ["a", "1.1.1", "1", "1", "1"] in clusters
        assert ["b", "1.1.2", "1", "1", "2"] in clusters
        assert ["c", "1.2.3", "1", "2", "3"] in clusters
        assert ["d", "2.4.5", "2", "4", "5"] in clusters
        assert ["e", "1.3.4", "1", "3", "4"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "25,18,0"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["level_1"] == 25.0
        assert run_json["threshold_map"]["level_2"] == 18.0
        assert run_json["threshold_map"]["level_3"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["level_1"] == 25.0
        assert thresholds_json["level_2"] == 18.0
        assert thresholds_json["level_3"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)
    
    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("((((b:17.000000,a:17.000000):4.0,c:21.000000):0.0,e:21.000000):7.0,d:28.000000);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(d:28.0,(e:21.0,(c:21.0,(a:17.0,b:17.0):4.0):0.0):7.0);"

def test_wikipedia_patristic(tmp_path):
    # Tests same wikipedia tree, but where patristic distances on newick correspond to 
    # distances in the distance matrix. Also shows that this parameter change does not effect
    # the clusters identified or cluster addresses
    args = {"matrix": get_path("data/matrix/wikipedia-single.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "25,18,0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'patristic'}

    # t=25
    # a,b,c,e
    # d

    # t=18
    # a,b
    # c
    # e
    # d

    # t=0
    # a
    # b
    # c
    # d
    # e

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address", "level_1", "level_2", "level_3"] in clusters
        assert ["a", "1.1.1", "1", "1", "1"] in clusters
        assert ["b", "1.1.2", "1", "1", "2"] in clusters
        assert ["c", "1.2.3", "1", "2", "3"] in clusters
        assert ["d", "2.4.5", "2", "4", "5"] in clusters
        assert ["e", "1.3.4", "1", "3", "4"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "25,18,0"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'patristic'

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["level_1"] == 25.0
        assert run_json["threshold_map"]["level_2"] == 18.0
        assert run_json["threshold_map"]["level_3"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["level_1"] == 25.0
        assert thresholds_json["level_2"] == 18.0
        assert thresholds_json["level_3"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    assert str(actual_tree).strip() == "(d:14.0,(e:10.5,(c:10.5,(a:8.5,b:8.5):2.0):0.0):3.5);"

def test_compare_tree_patristic_cophenetic_wikipedia(tmp_path):
    args_patristic = {
        "matrix": get_path("data/matrix/wikipedia-single.tsv"),
        "outdir": path.join(tmp_path, "test_out"),
        "method": "single",
        "thresholds": "25,18,0",
        "sort_matrix": False,
        "delimiter": ".",
        "force": False,
        "tree_distances": 'patristic'
    }

    args_cophenetic = {
        "matrix": get_path("data/matrix/wikipedia-single.tsv"),
        "outdir": path.join(tmp_path, "test_out2"),
        "method": "single",
        "thresholds": "25,18,0",
        "sort_matrix": False,
        "delimiter": ".",
        "force": False,
        "tree_distances": 'cophenetic'
    }

    mcluster(args_patristic)
    mcluster(args_cophenetic)

    input_distance_matrix = pd.read_csv(get_path("data/matrix/wikipedia-single.tsv"), sep='\t', index_col='dists')

    tree_path_patristic = path.join(args_patristic["outdir"], "tree.nwk")
    actual_tree_patristic = skbio.io.registry.read(tree_path_patristic, format='newick', into=TreeNode)

    tree_path_cophenetic = path.join(args_cophenetic["outdir"], "tree.nwk")
    actual_tree_cophenetic = skbio.io.registry.read(tree_path_cophenetic, format='newick', into=TreeNode)

    # The difference in Newick file between patristic and cophenetic is a factor of 2 on branch lengths
    assert str(actual_tree_patristic).strip()  == "(d:14.0,(e:10.5,(c:10.5,(a:8.5,b:8.5):2.0):0.0):3.5);"
    assert str(actual_tree_cophenetic).strip() == "(d:28.0,(e:21.0,(c:21.0,(a:17.0,b:17.0):4.0):0.0):7.0);"

    # For context for understanding tests, here is a dendrogram
    # Generated from newick tree string above and using https://github.com/JLSteenwyk/PhyKIT: "phykit print_tree tree.txt"
    #   __________________________________________________________________________ d
    # _|
    #  |                  ________________________________________________________ e
    #  |_________________|
    #                    |________________________________________________________ c
    #                    |
    #                    |           _____________________________________________ a
    #                    |__________|
    #                               |_____________________________________________ b
    #
    # And here is the input distance matrix ("data/matrix/wikipedia-single.tsv")
    # dists   a       b       c       d       e
    # a       0       17      21      31      23
    # b       17      0       30      34      21
    # c       21      30      0       28      39
    # d       31      34      28      0       43
    # e       23      21      39      43      0


    # For the newick tree using "--tree-distances patristic", the patristic distance (sum of branch lengths)
    # between two leaves 'a' and 'b' corresponds to the distance value from the input distance matrix between 'a' and 'b'
    assert distance_patristic_from_tree( actual_tree_patristic,  'a', 'b') == 17.0
    assert distance_from_matrix(         input_distance_matrix,  'a', 'b') == 17.0

    # In this scenario, calculating the cophenetic distance between two leaves 'a' and 'b' (height of lowest common ancestor)
    # Corresponds to a distance of 1/2 the value from the input matrix
    # Hence why the value of the parameter "--tree-distances patristic" is called "patristic"
    assert distance_cophenetic_from_tree(actual_tree_patristic,  'a', 'b') == 8.5
    assert distance_from_matrix(         input_distance_matrix,  'a', 'b') == 17.0

    # For the newick tree using "--tree-distances cophenetic", the patristic distance (sum of branch lengths)
    # between two leaves 'a' and 'b' corresponds to twice the distance value from the input distance matrix between 'a' and 'b'
    assert distance_patristic_from_tree( actual_tree_cophenetic, 'a', 'b') == 34.0
    assert distance_from_matrix(         input_distance_matrix,  'a', 'b') == 17.0
    # However, this means that the cophenetic distance between two leaves 'a' and 'b' (height of lowest common ancestor)
    # corresponds to the distance value from the input distance matrix
    # Hence why the value of the parameter "--tree-distances cophenetic" is called "cophenetic"
    assert distance_cophenetic_from_tree(actual_tree_cophenetic, 'a', 'b') == 17.0
    assert distance_from_matrix(         input_distance_matrix,  'a', 'b') == 17.0


def test_linkage_methods_wikipedia(tmp_path):
    args_single = {
        "matrix": get_path("data/matrix/wikipedia-single.tsv"),
        "outdir": path.join(tmp_path, "test_out_single"),
        "method": "single",
        "thresholds": "25,18,0",
        "sort_matrix": False,
        "delimiter": ".",
        "force": False,
        "tree_distances": 'patristic'
    }

    args_average = {
        "matrix": get_path("data/matrix/wikipedia-single.tsv"),
        "outdir": path.join(tmp_path, "test_out_average"),
        "method": "average",
        "thresholds": "25,18,0",
        "sort_matrix": False,
        "delimiter": ".",
        "force": False,
        "tree_distances": 'patristic'
    }

    args_complete = {
        "matrix": get_path("data/matrix/wikipedia-single.tsv"),
        "outdir": path.join(tmp_path, "test_out_complete"),
        "method": "complete",
        "thresholds": "25,18,0",
        "sort_matrix": False,
        "delimiter": ".",
        "force": False,
        "tree_distances": 'patristic'
    }

    mcluster(args_single)
    mcluster(args_average)
    mcluster(args_complete)

    input_distance_matrix = pd.read_csv(get_path("data/matrix/wikipedia-single.tsv"), sep='\t', index_col='dists')

    tree_path_single = path.join(args_single["outdir"], "tree.nwk")
    actual_tree_single = skbio.io.registry.read(tree_path_single, format='newick', into=TreeNode)

    tree_path_average = path.join(args_average["outdir"], "tree.nwk")
    actual_tree_average = skbio.io.registry.read(tree_path_average, format='newick', into=TreeNode)

    tree_path_complete = path.join(args_complete["outdir"], "tree.nwk")
    actual_tree_complete = skbio.io.registry.read(tree_path_complete, format='newick', into=TreeNode)

    assert str(actual_tree_single).strip()   == "(d:14.0,(e:10.5,(c:10.5,(a:8.5,b:8.5):2.0):0.0):3.5);"
    assert str(actual_tree_average).strip()  == "((e:11.0,(a:8.5,b:8.5):2.5):5.5,(c:14.0,d:14.0):2.5);"
    assert str(actual_tree_complete).strip() == "((e:11.5,(a:8.5,b:8.5):3.0):10.0,(c:14.0,d:14.0):7.5);"

    # Tree single-linkage (generated using https://github.com/JLSteenwyk/PhyKIT)
    #   __________________________________________________________________________ d
    # _|
    #  |                  ________________________________________________________ e
    #  |_________________|
    #                    |________________________________________________________ c
    #                    |
    #                    |           _____________________________________________ a
    #                    |__________|
    #                               |_____________________________________________ b

    # Tree average-linkage #########################################################
    #                            _________________________________________________ e
    #   ________________________|
    #  |                        |           ______________________________________ a
    #  |                        |__________|
    # _|                                   |______________________________________ b
    #  |
    #  |           _______________________________________________________________ c
    #  |__________|
    #             |_______________________________________________________________ d

    # Tree complete-linkage ########################################################
    #                                     ________________________________________ e
    #   _________________________________|
    #  |                                 |           _____________________________ a
    #  |                                 |__________|
    # _|                                            |_____________________________ b
    #  |
    #  |                          ________________________________________________ c
    #  |_________________________|
    #                            |________________________________________________ d

    # And here is the input distance matrix ("data/matrix/wikipedia-single.tsv")
    # dists   a       b       c       d       e
    # a       0       17      21      31      23
    # b       17      0       30      34      21
    # c       21      30      0       28      39
    # d       31      34      28      0       43
    # e       23      21      39      43      0

    # The distance (patristic) from 'a' to 'b' is the same in all three scenarios since these are the
    # closest nodes in the matrix (so get joined into a cluster first by the algorithm)
    assert distance_patristic_from_tree(actual_tree_single,    'a', 'b') == 17.0
    assert distance_patristic_from_tree(actual_tree_average,   'a', 'b') == 17.0
    assert distance_patristic_from_tree(actual_tree_complete,  'a', 'b') == 17.0
    assert distance_from_matrix(        input_distance_matrix, 'a', 'b') == 17.0

    # However, the distance (patristic) from 'e' to 'd' is different since it is
    # the furthest distance from the distance matrix (so the clustering algorithm handles these nodes later)
    # For single-linkage, it's the smallest distance between clusters (d) and (e, c, a, b)
    # See https://en.wikipedia.org/wiki/Single-linkage_clustering#Final_step
    assert distance_patristic_from_tree(actual_tree_single,    'd', 'e') == 28.0
    # For average-linkage, it's the distance from (d) to the internal node for cluster (c, d) [value of 14]
    # plus the distance from cluster (c, d) to the cluster (e) [value of 19]. That is 14 + 19 = 33
    # See https://en.wikipedia.org/wiki/UPGMA#The_UPGMA_dendrogram
    assert distance_patristic_from_tree(actual_tree_average,   'd', 'e') == 33.0
    # For complete-linkage, it's the distance from (d) to the internal node for cluster (c, d) [value fo 14] plus
    # the distance between clusters (c, d) and (e) [value of 29]. That is 14 + 29 = 43
    # See https://en.wikipedia.org/wiki/Complete-linkage_clustering#The_complete-linkage_dendrogram
    assert distance_patristic_from_tree(actual_tree_complete,  'd', 'e') == 43.0
    # This means the distances calculated from a tree won't always correspond to the distance from the
    # input distance matrix (value depends on if using single, average, or complete linkage)
    assert distance_from_matrix(        input_distance_matrix, 'd', 'e') == 43.0


def test_threshold_same(tmp_path):
    # Tests behaviour that similar thresholds create similar clusters.
    args = {"matrix": get_path("data/matrix/wikipedia-single.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "20,19,18",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    """
    The linkage will look like this:

    [ 0.  1. 17.  2.]
    [ 2.  5. 21.  3.]
    [ 4.  6. 21.  4.]
    [ 3.  7. 28.  5.]

    Reminder that the third column is distances, so thresholds of
    20,19,18 will split the linkage at the same place every time:

    {a,b}, {c}, {d}, {e}
    """

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address", "level_1", "level_2", "level_3"] in clusters
        assert ["a", "1.1.1", "1", "1", "1"] in clusters
        assert ["b", "1.1.1", "1", "1", "1"] in clusters
        assert ["c", "2.2.2", "2", "2", "2"] in clusters
        assert ["d", "4.4.4", "4", "4", "4"] in clusters
        assert ["e", "3.3.3", "3", "3", "3"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "20,19,18"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["level_1"] == 20.0
        assert run_json["threshold_map"]["level_2"] == 19.0
        assert run_json["threshold_map"]["level_3"] == 18.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["level_1"] == 20.0
        assert thresholds_json["level_2"] == 19.0
        assert thresholds_json["level_3"] == 18.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("((((b:17.000000,a:17.000000):4.0,c:21.000000):0.0,e:21.000000):7.0,d:28.000000);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(d:28.0,(e:21.0,(c:21.0,(a:17.0,b:17.0):4.0):0.0):7.0);"

def test_thresholds_0_10_0_10(tmp_path):
    # "thresholds": "0,10,0,10"
    # This should fail because thresholds must be decreasing.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0,10,0,10",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['0', '10', '0', '10'] must be in decreasing order"

    assert path.isdir(args["outdir"]) == False

def test_thresholds_0_0(tmp_path):
    # "thresholds": "0,0"
    # This should fail because thresholds must be decreasing.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0,0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['0', '0'] must be in decreasing order"

    assert path.isdir(args["outdir"]) == False

def test_thresholds_1_2_3(tmp_path):
    # "thresholds": "1,2,3"
    # This should fail because thresholds must be decreasing.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,2,3",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['1', '2', '3'] must be in decreasing order"

    assert path.isdir(args["outdir"]) == False

def test_thresholds_string(tmp_path):
    # "thresholds": "cat,dog"
    # This should fail because thresholds must be integers or floats.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "cat,dog",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['cat', 'dog'] must all be integers or floats"

    assert path.isdir(args["outdir"]) == False

def test_no_thresholds(tmp_path):
    # "thresholds": ""
    # This should fail because there are no thresholds.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds [''] must all be integers or floats"

    assert path.isdir(args["outdir"]) == False

def test_delimiter_slash(tmp_path):
    # "delimiter": "/"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,0",
            "sort_matrix": False,
            "delimiter": "/",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # A, B, I, J unique
        # C, D share
        # E, F share
        # G, H share
        assert ["id", "address", "level_1", "level_2"] in clusters
        assert ["A", "3/3", "3", "3"] in clusters
        assert ["B", "3/4", "3", "4"] in clusters
        assert ["C", "2/2", "2", "2"] in clusters
        assert ["D", "2/2", "2", "2"] in clusters
        assert ["E", "1/1", "1", "1"] in clusters
        assert ["F", "1/1", "1", "1"] in clusters
        assert ["G", "4/5", "4", "5"] in clusters
        assert ["H", "4/5", "4", "5"] in clusters
        assert ["I", "5/6", "5", "6"] in clusters
        assert ["J", "5/7", "5", "7"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "1,0"
        assert run_json["parameters"]["delimiter"] == "/"
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 2
        assert run_json["threshold_map"]["level_1"] == 1.0
        assert run_json["threshold_map"]["level_2"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 2
        assert thresholds_json["level_1"] == 1.0
        assert thresholds_json["level_2"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("(((J:1.000000,I:1.000000):2.0,(H:0.000000,G:0.000000):3.0):3.0,(((B:1.000000,A:1.000000):1.0,(D:0.000000,C:0.000000):2.0):1.0,(F:0.000000,E:0.000000):3.0):3.0);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(((E:0.0,F:0.0):3.0,((C:0.0,D:0.0):2.0,(A:1.0,B:1.0):1.0):1.0):3.0,((G:0.0,H:0.0):3.0,(I:1.0,J:1.0):2.0):3.0);"

def test_delimiter_0(tmp_path):
    # "delimiter": "0"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,0",
            "sort_matrix": False,
            "delimiter": "0",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # A, B, I, J unique
        # C, D share
        # E, F share
        # G, H share
        assert ["id", "address", "level_1", "level_2"] in clusters
        assert ["A", "303", "3", "3"] in clusters
        assert ["B", "304", "3", "4"] in clusters
        assert ["C", "202", "2", "2"] in clusters
        assert ["D", "202", "2", "2"] in clusters
        assert ["E", "101", "1", "1"] in clusters
        assert ["F", "101", "1", "1"] in clusters
        assert ["G", "405", "4", "5"] in clusters
        assert ["H", "405", "4", "5"] in clusters
        assert ["I", "506", "5", "6"] in clusters
        assert ["J", "507", "5", "7"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "1,0"
        assert run_json["parameters"]["delimiter"] == "0"
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 2
        assert run_json["threshold_map"]["level_1"] == 1.0
        assert run_json["threshold_map"]["level_2"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 2
        assert thresholds_json["level_1"] == 1.0
        assert thresholds_json["level_2"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("(((J:1.000000,I:1.000000):2.0,(H:0.000000,G:0.000000):3.0):3.0,(((B:1.000000,A:1.000000):1.0,(D:0.000000,C:0.000000):2.0):1.0,(F:0.000000,E:0.000000):3.0):3.0);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(((E:0.0,F:0.0):3.0,((C:0.0,D:0.0):2.0,(A:1.0,B:1.0):1.0):1.0):3.0,((G:0.0,H:0.0):3.0,(I:1.0,J:1.0):2.0):3.0);"

def test_delimiter_1(tmp_path):
    # "delimiter": "1"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,0",
            "sort_matrix": False,
            "delimiter": "1",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # A, B, I, J unique
        # C, D share
        # E, F share
        # G, H share
        assert ["id", "address", "level_1", "level_2"] in clusters
        assert ["A", "313", "3", "3"] in clusters
        assert ["B", "314", "3", "4"] in clusters
        assert ["C", "212", "2", "2"] in clusters
        assert ["D", "212", "2", "2"] in clusters
        assert ["E", "111", "1", "1"] in clusters
        assert ["F", "111", "1", "1"] in clusters
        assert ["G", "415", "4", "5"] in clusters
        assert ["H", "415", "4", "5"] in clusters
        assert ["I", "516", "5", "6"] in clusters
        assert ["J", "517", "5", "7"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "1,0"
        assert run_json["parameters"]["delimiter"] == "1"
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 2
        assert run_json["threshold_map"]["level_1"] == 1.0
        assert run_json["threshold_map"]["level_2"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 2
        assert thresholds_json["level_1"] == 1.0
        assert thresholds_json["level_2"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("(((J:1.000000,I:1.000000):2.0,(H:0.000000,G:0.000000):3.0):3.0,(((B:1.000000,A:1.000000):1.0,(D:0.000000,C:0.000000):2.0):1.0,(F:0.000000,E:0.000000):3.0):3.0);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(((E:0.0,F:0.0):3.0,((C:0.0,D:0.0):2.0,(A:1.0,B:1.0):1.0):1.0):3.0,((G:0.0,H:0.0):3.0,(I:1.0,J:1.0):2.0):3.0);"

def test_delimiter_quote(tmp_path):
    # "delimiter": '"'
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,0",
            "sort_matrix": False,
            "delimiter": '"',
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # A, B, I, J unique
        # C, D share
        # E, F share
        # G, H share
        assert ["id", "address", "level_1", "level_2"] in clusters
        assert ["A", "3\"3", "3", "3"] in clusters
        assert ["B", "3\"4", "3", "4"] in clusters
        assert ["C", "2\"2", "2", "2"] in clusters
        assert ["D", "2\"2", "2", "2"] in clusters
        assert ["E", "1\"1", "1", "1"] in clusters
        assert ["F", "1\"1", "1", "1"] in clusters
        assert ["G", "4\"5", "4", "5"] in clusters
        assert ["H", "4\"5", "4", "5"] in clusters
        assert ["I", "5\"6", "5", "6"] in clusters
        assert ["J", "5\"7", "5", "7"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "1,0"
        assert run_json["parameters"]["delimiter"] == '"'
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 2
        assert run_json["threshold_map"]["level_1"] == 1.0
        assert run_json["threshold_map"]["level_2"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 2
        assert thresholds_json["level_1"] == 1.0
        assert thresholds_json["level_2"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("(((J:1.000000,I:1.000000):2.0,(H:0.000000,G:0.000000):3.0):3.0,(((B:1.000000,A:1.000000):1.0,(D:0.000000,C:0.000000):2.0):1.0,(F:0.000000,E:0.000000):3.0):3.0);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(((E:0.0,F:0.0):3.0,((C:0.0,D:0.0):2.0,(A:1.0,B:1.0):1.0):1.0):3.0,((G:0.0,H:0.0):3.0,(I:1.0,J:1.0):2.0):3.0);"

def test_matrix_missing(tmp_path):
    # Missing input file.
    args = {"matrix": get_path("data/matrix/basic_MISSING.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert "data/matrix/basic_MISSING.tsv does not exist or is empty" in str(exception.value)

    assert path.isdir(args["outdir"]) == False

def test_matrix_empty(tmp_path):
    # Input file is empty.
    args = {"matrix": get_path("data/matrix/empty.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert "data/matrix/empty.tsv does not exist or is empty" in str(exception.value)

    assert path.isdir(args["outdir"]) == False

def test_method_invalid_nope(tmp_path):
    # method = "nope"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "nope",
            "thresholds": "0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "nope is not one of the accepeted methods ['average', 'complete', 'single']"

    assert path.isdir(args["outdir"]) == False

def test_method_invalid_singl(tmp_path):
    # method = "singl"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "singl",
            "thresholds": "0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "singl is not one of the accepeted methods ['average', 'complete', 'single']"

    assert path.isdir(args["outdir"]) == False

def test_method_invalid_1(tmp_path):
    # method = "1"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "1",
            "thresholds": "0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "1 is not one of the accepeted methods ['average', 'complete', 'single']"

    assert path.isdir(args["outdir"]) == False

def test_many_thresholds(tmp_path):
    # Checks the property that so long as thresholds are always
    # decreasing, the correspondingly assigned labels can never
    # decrease. For example:
    # Good: 1.1.1.2.3
    # Bad: 1.1.2.1.1
    args = {"matrix": get_path("data/matrix/wikipedia-single.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "26,24,22,20,18,16,14,12,10,8,6,4,2,0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["a", "1.1.1.1.1.1.1.1.1.1.1.1.1.1", "1","1","1","1","1","1","1","1","1","1","1","1","1","1"] in clusters
        assert ["b", "1.1.1.1.1.2.2.2.2.2.2.2.2.2", "1","1","1","1","1","2","2","2","2","2","2","2","2","2"] in clusters
        assert ["c", "1.1.1.2.2.3.3.3.3.3.3.3.3.3", "1","1","1","2","2","3","3","3","3","3","3","3","3","3"] in clusters
        assert ["d", "2.2.2.4.4.5.5.5.5.5.5.5.5.5", "2","2","2","4","4","5","5","5","5","5","5","5","5","5"] in clusters
        assert ["e", "1.1.1.3.3.4.4.4.4.4.4.4.4.4", "1","1","1","3","3","4","4","4","4","4","4","4","4","4"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "26,24,22,20,18,16,14,12,10,8,6,4,2,0"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 14
        assert run_json["threshold_map"]["level_1"] == 26.0
        assert run_json["threshold_map"]["level_2"] == 24.0
        assert run_json["threshold_map"]["level_3"] == 22.0
        assert run_json["threshold_map"]["level_4"] == 20.0
        assert run_json["threshold_map"]["level_5"] == 18.0
        assert run_json["threshold_map"]["level_6"] == 16.0
        assert run_json["threshold_map"]["level_7"] == 14.0
        assert run_json["threshold_map"]["level_8"] == 12.0
        assert run_json["threshold_map"]["level_9"] == 10.0
        assert run_json["threshold_map"]["level_10"] == 8.0
        assert run_json["threshold_map"]["level_11"] == 6.0
        assert run_json["threshold_map"]["level_12"] == 4.0
        assert run_json["threshold_map"]["level_13"] == 2.0
        assert run_json["threshold_map"]["level_14"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 14
        assert thresholds_json["level_1"] == 26.0
        assert thresholds_json["level_2"] == 24.0
        assert thresholds_json["level_3"] == 22.0
        assert thresholds_json["level_4"] == 20.0
        assert thresholds_json["level_5"] == 18.0
        assert thresholds_json["level_6"] == 16.0
        assert thresholds_json["level_7"] == 14.0
        assert thresholds_json["level_8"] == 12.0
        assert thresholds_json["level_9"] == 10.0
        assert thresholds_json["level_10"] == 8.0
        assert thresholds_json["level_11"] == 6.0
        assert thresholds_json["level_12"] == 4.0
        assert thresholds_json["level_13"] == 2.0
        assert thresholds_json["level_14"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("((((b:17.000000,a:17.000000):4.0,c:21.000000):0.0,e:21.000000):7.0,d:28.000000);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(d:28.0,(e:21.0,(c:21.0,(a:17.0,b:17.0):4.0):0.0):7.0);"

def test_method_single(tmp_path):
    # method = "single"
    # thresholds = "2,1,0"

    """
    The distance matrix at "data/matrix/basic.tsv" will generate
    the following single linkage:

    10: [ 2.  3.  0.  2.] -> (C,D)
    11: [ 4.  5.  0.  2.] -> (E,F)
    12: [ 6.  7.  0.  2.] -> (G,H)
    13: [ 0.  1.  1.  2.] -> (A,B)
    14: [ 8.  9.  1.  2.] -> (I,J)
    15: [10. 13.  2.  4.] -> ((C,D), (A,B))
    16: [11. 15.  3.  6.] -> ((E,F), ((C,D), (A,B)))
    17: [12. 14.  3.  4.] -> ((G,H), (I,J))
    18: [16. 17.  6. 10.] -> (((E,F), ((C,D), (A,B))), ((G,H), (I,J)))

    When threshold=2, the following will be grouped and labelled
    together when flattened:

    - ((C,D), (A,B))
    - (G,H)
    - (I,J)
    - (E,F)

    When threshold=1:

    - (C,D)
    - (A,B)
    - (G,H)
    - (I,J)
    - (E,F)

    When threshold=0:

    - A
    - B
    - (C,D)
    - (G,H)
    - (E,F)
    - I
    - J
    """

    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "2,1,0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address", "level_1", "level_2", "level_3"] in clusters
        assert ["A", "2.3.3", "2", "3", "3"] in clusters
        assert ["B", "2.3.4", "2", "3", "4"] in clusters
        assert ["C", "2.2.2", "2", "2", "2"] in clusters
        assert ["D", "2.2.2", "2", "2", "2"] in clusters
        assert ["E", "1.1.1", "1", "1", "1"] in clusters
        assert ["F", "1.1.1", "1", "1", "1"] in clusters
        assert ["G", "3.4.5", "3", "4", "5"] in clusters
        assert ["H", "3.4.5", "3", "4", "5"] in clusters
        assert ["I", "4.5.6", "4", "5", "6"] in clusters
        assert ["J", "4.5.7", "4", "5", "7"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "2,1,0"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["level_1"] == 2.0
        assert run_json["threshold_map"]["level_2"] == 1.0
        assert run_json["threshold_map"]["level_3"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["level_1"] == 2.0
        assert thresholds_json["level_2"] == 1.0
        assert thresholds_json["level_3"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("(((J:1.000000,I:1.000000):2.0,(H:0.000000,G:0.000000):3.0):3.0,(((B:1.000000,A:1.000000):1.0,(D:0.000000,C:0.000000):2.0):1.0,(F:0.000000,E:0.000000):3.0):3.0);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(((E:0.0,F:0.0):3.0,((C:0.0,D:0.0):2.0,(A:1.0,B:1.0):1.0):1.0):3.0,((G:0.0,H:0.0):3.0,(I:1.0,J:1.0):2.0):3.0);"

def test_method_complete(tmp_path):
    # method = "complete"
    # thresholds = "10,8,5"

    """
    The distance matrix at "data/matrix/small.tsv" will generate
    the following single linkage:

    5: [0. 1. 1. 2.] -> (A,B)
    6: [2. 5. 5. 3.] -> (C,(A,B))
    7: [3. 6. 8. 4.] -> (D,(C,(A,B)))
    8: [4. 7. 10. 5.] -> (E,(D,(C,(A,B))))

    When threshold=10, the following will be grouped and labelled
    together when flattened:

    - {A, B, C, D, E}

    When threshold=8:

    - {A, B, C, D}, {E}

    When threshold=5:

    - {A, B, C}, {D}, {E}
    """

    args = {"matrix": get_path("data/matrix/small.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "complete",
            "thresholds": "10,8,5",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address", "level_1", "level_2", "level_3"] in clusters
        assert ["A", "1.1.1", "1", "1", "1"] in clusters
        assert ["B", "1.1.1", "1", "1", "1"] in clusters
        assert ["C", "1.1.1", "1", "1", "1"] in clusters
        assert ["D", "1.1.2", "1", "1", "2"] in clusters
        assert ["E", "1.2.3", "1", "2", "3"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "complete"
        assert run_json["parameters"]["thresholds"] == "10,8,5"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["level_1"] == 10.0
        assert run_json["threshold_map"]["level_2"] == 8.0
        assert run_json["threshold_map"]["level_3"] == 5.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["level_1"] == 10.0
        assert thresholds_json["level_2"] == 8.0
        assert thresholds_json["level_3"] == 5.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("((((B:1.000000,A:1.000000):4.0,C:5.000000):3.0,D:8.000000):2.0,E:10.000000);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(E:10.0,(D:8.0,(C:5.0,(A:1.0,B:1.0):4.0):3.0):2.0);"

def test_method_average(tmp_path):
    # method = "average"
    # thresholds = "12,8,2"

    """
    The distance matrix at "data/matrix/small-average-linkage.tsv"
    will generate the following single linkage:

    3: [ 0.  1.  4.  2.] -> (A,B)
    4: [ 2.  3. 10.  3.] -> (C, (A,B))

    When threshold=12, the following will be grouped and labelled
    together when flattened:

    - {A, B, C}

    When threshold=8:

    - {A, B}, {C}

    When threshold=2:

    - {A}, {B}, {C}
    """

    args = {"matrix": get_path("data/matrix/small-average-linkage.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "complete",
            "thresholds": "12,8,2",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address", "level_1", "level_2", "level_3"] in clusters
        assert ["A", "1.1.1", "1", "1", "1"] in clusters
        assert ["B", "1.1.2", "1", "1", "2"] in clusters
        assert ["C", "1.2.3", "1", "2", "3"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "complete"
        assert run_json["parameters"]["thresholds"] == "12,8,2"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["level_1"] == 12.0
        assert run_json["threshold_map"]["level_2"] == 8.0
        assert run_json["threshold_map"]["level_3"] == 2.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["level_1"] == 12.0
        assert thresholds_json["level_2"] == 8.0
        assert thresholds_json["level_3"] == 2.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("((B:4.000000,A:4.000000):8.0,C:12.000000);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0
    assert expected_tree.compare_cophenet(actual_tree) == 0
    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(C:12.0,(A:4.0,B:4.0):8.0);"

def test_invalid_header_pairwise_matrix(tmp_path):
    matrix_path = get_path("data/matrix/csv.text")
    args = {"matrix": matrix_path,
            "outdir": path.join(tmp_path, "test_out"),
            "method": "1",
            "thresholds": "0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == f"{matrix_path} does not appear to be a properly TSV-formatted file"

    assert path.isdir(args["outdir"]) == False

def test_double_digit(tmp_path):
    # An example that tests double digit distances in the input
    # and double digit addresses.
    args = {"matrix": get_path("data/matrix/double_digit.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "5,3,0",
            "sort_matrix": False,
            "delimiter": ".",
            "force": False,
            "tree_distances": 'cophenetic'}

    mcluster(args)

    assert path.isdir(args["outdir"])

    clusters_path = path.join(args["outdir"], "clusters.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # all are unique
        assert ["id", "address", "level_1", "level_2", "level_3"] in clusters
        assert ["A", "1.1.1", "1", "1", "1"] in clusters
        assert ["B", "2.2.2", "2", "2", "2"] in clusters
        assert ["C", "3.3.3", "3", "3", "3"] in clusters
        assert ["D", "4.4.4", "4", "4", "4"] in clusters
        assert ["E", "5.5.5", "5", "5", "5"] in clusters
        assert ["F", "6.6.6", "6", "6", "6"] in clusters
        assert ["G", "7.7.7", "7", "7", "7"] in clusters
        assert ["H", "8.8.8", "8", "8", "8"] in clusters
        assert ["I", "9.9.9", "9", "9", "9"] in clusters
        assert ["J", "10.10.10", "10", "10", "10"] in clusters
        assert ["K", "11.11.11", "11", "11", "11"] in clusters
        assert ["L", "12.12.12", "12", "12", "12"] in clusters
        assert ["M", "13.13.13", "13", "13", "13"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "5,3,0"
        assert run_json["parameters"]["delimiter"] == "."
        assert run_json["parameters"]["tree_distances"] == 'cophenetic'

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["level_1"] == 5.0
        assert run_json["threshold_map"]["level_2"] == 3.0
        assert run_json["threshold_map"]["level_3"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["level_1"] == 5.0
        assert thresholds_json["level_2"] == 3.0
        assert thresholds_json["level_3"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

    actual_tree = skbio.io.registry.read(tree_path, format='newick', into=TreeNode)
    expected_tree = skbio.io.registry.read(StringIO("((((((((((((B:25.000000,A:25.000000):0.0,C:25.000000):0.0,D:25.000000):0.0,E:25.000000):0.0,F:25.000000):0.0,G:25.000000):0.0,H:25.000000):0.0,I:25.000000):0.0,J:25.000000):0.0,K:25.000000):0.0,L:25.000000):0.0,M:25.000000);"),
                                           format='newick', into=TreeNode)
    assert sorted([t.name for t in expected_tree.tips()]) == sorted([t.name for t in actual_tree.tips()])
    assert expected_tree.compare_rfd(actual_tree) == 0

    # In this specific example, the compare_cophenet() between the trees results in nan instead of 0
    # so instead I compare the underlying data structure (distance matrix) between each tree
    # Note, even though ids/tree tips are in a different order, all distances between pairs are the same
    # so comparison of the two distance matrices can be done without re-ordering
    #assert expected_tree.compare_cophenet(actual_tree) == 0
    assert sorted(expected_tree.cophenet().ids) == sorted(actual_tree.cophenet().ids)
    assert (expected_tree.cophenet().data == actual_tree.cophenet().data).all()

    assert expected_tree.compare_subsets(actual_tree) == 0

    assert str(actual_tree).strip() == "(M:25.0,(L:25.0,(K:25.0,(J:25.0,(I:25.0,(H:25.0,(G:25.0,(F:25.0,(E:25.0,(D:25.0,(C:25.0,(A:25.0,B:25.0):0.0):0.0):0.0):0.0):0.0):0.0):0.0):0.0):0.0):0.0):0.0);"

def test_sort_matrix(tmp_path):
    # Compare outputs of mcluster on same input matrix
    # With different sample orders
    # Should produce identical outputs
    # when sort_matrix=True
    
    args_sorted = {
        "matrix": get_path("data/matrix/sorted.tsv"),
        "outdir": path.join(tmp_path, "test_out"),
        "method": "single",
        "thresholds": "0",
        "sort_matrix": False,
        "delimiter": ".",
        "force": False,
        "tree_distances": 'cophenetic'
    }

    args_shuffled = {
        "matrix": get_path("data/matrix/shuffled.tsv"),
        "outdir": path.join(tmp_path, "test_out2"),
        "method": "single",
        "thresholds": "0",
        "sort_matrix": True, # This will sort the input matrix to match the sorted.tsv
        "delimiter": ".",
        "force": False,
        "tree_distances": 'cophenetic'
    }

    mcluster(args_sorted)
    mcluster(args_shuffled)
    
    cluster_sorted = path.join(args_sorted["outdir"], "clusters.text")
    cluster_sorted_output = pd.read_csv(cluster_sorted, sep='\t')
    cluster_shuffled = path.join(args_shuffled["outdir"], "clusters.text")
    cluster_shuffled_output = pd.read_csv(cluster_shuffled, sep='\t')

    assert cluster_sorted_output.equals(cluster_shuffled_output)