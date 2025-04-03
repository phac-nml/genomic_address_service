"""
A simple test for running gas call

"""

import pytest
import pathlib as p
from os import path
import csv
import json

from genomic_address_service.call import call


def get_path(location):
    directory = path.dirname(path.abspath(__file__))
    return path.join(directory, location)

def test_run_call(tmp_path):
    """
    Run call takes in a dictionary of parameters based off of the CLI.
    """
    input_results = p.Path("genomic_address_service/example/call/hamming/gas/results.text")
    config = {}
    config['dists'] = "genomic_address_service/example/call/hamming/results.text"
    config['rclusters'] = "genomic_address_service/example/call/clusters.text"
    config['thresh_map'] = None
    config['outdir'] = str(tmp_path / "test_out")
    config['method'] = "average"
    config['thresholds'] = "10,9,8,7,6,5,4,3,2,1,0"
    config['delimeter'] = "."
    config['force'] = False
    config['outfmt'] = "text"
    config['address_col'] = "address"
    config['batch_size'] = 100
    config['sample_col'] = "id"
    call(config)
    outdir = p.Path(config["outdir"])
    output_file = outdir / "results.text"
    assert output_file.read_text().split("\n").sort() == input_results.read_text().split("\n").sort()

def test_basic(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["outfmt"] = "text"
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimeter"] = "."

    config["batch_size"] = 100

    call(config)

    assert path.isdir(output_path)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # The new E is the same as the existing B (1.1.2)
        # The new F is the same as the existing D (1.1.4)

        assert ["id", "address"] in clusters
        assert ["A", "1.1.1"] in clusters
        assert ["B", "1.1.2"] in clusters
        assert ["C", "1.1.3"] in clusters
        assert ["D", "1.1.4"] in clusters
        assert ["E", "1.1.2"] in clusters
        assert ["F", "1.1.4"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "5,3,0"
        assert run_json["parameters"]["delimeter"] == "."

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["0"] == 5.0
        assert run_json["threshold_map"]["1"] == 3.0
        assert run_json["threshold_map"]["2"] == 0.0

    # Thresholds JSON
    thresholds_path = path.join(output_path, "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["0"] == 5.0
        assert thresholds_json["1"] == 3.0
        assert thresholds_json["2"] == 0.0
