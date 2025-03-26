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
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "average",
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
        assert ["A", "6", "6"] in clusters
        assert ["B", "7", "7"] in clusters
        assert ["C", "5", "5"] in clusters
        assert ["D", "5", "5"] in clusters
        assert ["E", "4", "4"] in clusters
        assert ["F", "4", "4"] in clusters
        assert ["G", "1", "1"] in clusters
        assert ["H", "1", "1"] in clusters
        assert ["I", "2", "2"] in clusters
        assert ["J", "3", "3"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "average"
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

def test_delimeter_slash(tmp_path):
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "average",
            "thresholds": "0,0",
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
        assert ["A", "6/6", "6", "6"] in clusters
        assert ["B", "7/7", "7", "7"] in clusters
        assert ["C", "5/5", "5", "5"] in clusters
        assert ["D", "5/5", "5", "5"] in clusters
        assert ["E", "4/4", "4", "4"] in clusters
        assert ["F", "4/4", "4", "4"] in clusters
        assert ["G", "1/1", "1", "1"] in clusters
        assert ["H", "1/1", "1", "1"] in clusters
        assert ["I", "2/2", "2", "2"] in clusters
        assert ["J", "3/3", "3", "3"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "average"
        assert run_json["parameters"]["thresholds"] == "0,0"
        assert run_json["parameters"]["delimeter"] == "/"

        assert len(run_json["threshold_map"]) == 2
        assert run_json["threshold_map"]["level_1"] == 0.0
        assert run_json["threshold_map"]["level_2"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 2
        assert thresholds_json["level_1"] == 0.0
        assert thresholds_json["level_2"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

def test_delimeter_0(tmp_path):
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "average",
            "thresholds": "0,0",
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
        assert ["A", "606", "6", "6"] in clusters
        assert ["B", "707", "7", "7"] in clusters
        assert ["C", "505", "5", "5"] in clusters
        assert ["D", "505", "5", "5"] in clusters
        assert ["E", "404", "4", "4"] in clusters
        assert ["F", "404", "4", "4"] in clusters
        assert ["G", "101", "1", "1"] in clusters
        assert ["H", "101", "1", "1"] in clusters
        assert ["I", "202", "2", "2"] in clusters
        assert ["J", "303", "3", "3"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "average"
        assert run_json["parameters"]["thresholds"] == "0,0"
        assert run_json["parameters"]["delimeter"] == "0"

        assert len(run_json["threshold_map"]) == 2
        assert run_json["threshold_map"]["level_1"] == 0.0
        assert run_json["threshold_map"]["level_2"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 2
        assert thresholds_json["level_1"] == 0.0
        assert thresholds_json["level_2"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)

def test_delimeter_quote(tmp_path):
    args = {"matrix": get_path("data/matrix/basic.tsv"),
            "outdir": path.join(tmp_path, "test_out"),
            "method": "average",
            "thresholds": "0,0",
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
        assert ["A", '6"6', "6", "6"] in clusters
        assert ["B", '7"7', "7", "7"] in clusters
        assert ["C", '5"5', "5", "5"] in clusters
        assert ["D", '5"5', "5", "5"] in clusters
        assert ["E", '4"4', "4", "4"] in clusters
        assert ["F", '4"4', "4", "4"] in clusters
        assert ["G", '1"1', "1", "1"] in clusters
        assert ["H", '1"1', "1", "1"] in clusters
        assert ["I", '2"2', "2", "2"] in clusters
        assert ["J", '3"3', "3", "3"] in clusters

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "average"
        assert run_json["parameters"]["thresholds"] == "0,0"
        assert run_json["parameters"]["delimeter"] == "\""

        assert len(run_json["threshold_map"]) == 2
        assert run_json["threshold_map"]["level_1"] == 0.0
        assert run_json["threshold_map"]["level_2"] == 0.0

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 2
        assert thresholds_json["level_1"] == 0.0
        assert thresholds_json["level_2"] == 0.0

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)
