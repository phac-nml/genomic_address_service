import pytest
import pathlib
import csv
import json

from os import path
from genomic_address_service.mcluster import mcluster

def get_path(location):
    directory = path.dirname(path.abspath(__file__))
    return path.join(directory, location)

def test_basic_0(tmp_path):
    # A basic example with one threshold (0) where every
    # item should only cluster with itself or others with
    # a distance of 0.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0",
            "delimeter": ".",
            "force": False}

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
        assert run_json["parameters"]["delimeter"] == "."

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

def test_basic_wikipedia(tmp_path):
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
            "delimeter": ".",
            "force": False}

    print(args["outdir"])
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
        assert run_json["parameters"]["delimeter"] == "."

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

def test_thresholds_0_10_0_10(tmp_path):
    # "thresholds": "0,10,0,10"
    # This should fail because thresholds must be decreasing.
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0,10,0,10",
            "delimeter": ".",
            "force": False}

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
            "delimeter": ".",
            "force": False}

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
            "delimeter": ".",
            "force": False}

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
            "delimeter": ".",
            "force": False}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['cat', 'dog'] must all be integers or floats"

    assert path.isdir(args["outdir"]) == False

def test_delimeter_slash(tmp_path):
    # "delimeter": "/"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,0",
            "delimeter": "/",
            "force": False}

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
        assert run_json["parameters"]["delimeter"] == "/"

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

def test_delimeter_0(tmp_path):
    # "delimeter": "0"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,0",
            "delimeter": "0",
            "force": False}

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
        assert run_json["parameters"]["delimeter"] == "0"

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

def test_delimeter_1(tmp_path):
    # "delimeter": "1"
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,0",
            "delimeter": "1",
            "force": False}

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
        assert run_json["parameters"]["delimeter"] == "1"

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

def test_delimeter_quote(tmp_path):
    # "delimeter": '"'
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "1,0",
            "delimeter": '"',
            "force": False}

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
        assert run_json["parameters"]["delimeter"] == '"'

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

def test_matrix_missing(tmp_path):
    # Missing input file.
    args = {"matrix": get_path("data/matrix/basic_MISSING.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "single",
            "thresholds": "0",
            "delimeter": ".",
            "force": False}

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
            "delimeter": ".",
            "force": False}

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
            "delimeter": ".",
            "force": False}

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
            "delimeter": ".",
            "force": False}

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
            "delimeter": ".",
            "force": False}

    with pytest.raises(Exception) as exception:
        mcluster(args)

    assert exception.type == Exception
    assert str(exception.value) == "1 is not one of the accepeted methods ['average', 'complete', 'single']"

    assert path.isdir(args["outdir"]) == False

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
            "delimeter": ".",
            "force": False}

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
        assert run_json["parameters"]["delimeter"] == "."

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
            "delimeter": ".",
            "force": False}

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
        assert run_json["parameters"]["delimeter"] == "."

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
