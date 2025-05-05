"""
A simple test for running gas call

"""

import pytest
import pathlib as p
from os import path
import csv
import json

from genomic_address_service.call import call
from genomic_address_service.constants import CLUSTER_METHODS


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
    config['delimiter'] = "."
    config['force'] = False
    config['address_col'] = "address"
    config['batch_size'] = 100
    config['sample_col'] = "id"
    call(config)
    outdir = p.Path(config["outdir"])
    output_file = outdir / "results.text"
    assert output_file.read_text().split("\n").sort() == input_results.read_text().split("\n").sort()

def test_basic(tmp_path):
    # A basic test of GAS call.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

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
        assert run_json["parameters"]["delimiter"] == "."

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

def test_threshold_same(tmp_path):
    # A test GAS call when multiple thresholds cause the same integer
    # label to be assigned.
    config = {}

    clusters_path = get_path("data/clusters/threshold_same.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "10,9,8"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    call(config)

    assert path.isdir(output_path)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # The new E is the same as the existing B (1.1.1)
        # The new F is the same as the existing D (1.1.1)

        assert ["id", "address"] in clusters
        assert ["A", "1.1.1"] in clusters
        assert ["B", "1.1.1"] in clusters
        assert ["C", "1.1.1"] in clusters
        assert ["D", "1.1.1"] in clusters
        assert ["E", "1.1.1"] in clusters
        assert ["F", "1.1.1"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "10,9,8"
        assert run_json["parameters"]["delimiter"] == "."

        assert len(run_json["threshold_map"]) == 3
        assert run_json["threshold_map"]["0"] == 10.0
        assert run_json["threshold_map"]["1"] == 9.0
        assert run_json["threshold_map"]["2"] == 8.0

    # Thresholds JSON
    thresholds_path = path.join(output_path, "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 3
        assert thresholds_json["0"] == 10.0
        assert thresholds_json["1"] == 9.0
        assert thresholds_json["2"] == 8.0

def test_thresholds_0_10_0(tmp_path):
    # "thresholds": "0,10,0"
    # This should fail because thresholds must be decreasing.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "0,10,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['0', '10', '0'] must be in decreasing order"

    assert path.isdir(output_path) == False

def test_thresholds_0_0(tmp_path):
    # "thresholds": "0,0"
    # This should fail because thresholds must be decreasing.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "0,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['0', '0'] must be in decreasing order"

    assert path.isdir(output_path) == False

def test_thresholds_1_2_3(tmp_path):
    # "thresholds": "1,2,3"
    # This should fail because thresholds must be decreasing.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "1,2,3"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['1', '2', '3'] must be in decreasing order"

    assert path.isdir(output_path) == False

def test_thresholds_string(tmp_path):
    # "thresholds": "cat,dog"
    # This should fail because thresholds must be integers or floats.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "cat,dog"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds ['cat', 'dog'] must all be integers or floats"

    assert path.isdir(output_path) == False

def test_no_thresholds(tmp_path):
    # "thresholds": ""
    # This should fail because there are no thresholds.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = ""
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == "thresholds [''] must all be integers or floats"

    assert path.isdir(output_path) == False

def test_thresholds_many(tmp_path):
    # A test where many thresholds are used.
    config = {}

    clusters_path = get_path("data/clusters/thresholds_many.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/thresholds_many.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "26,24,22,20,18,16,14,12,10,8,6,4,2,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    call(config)

    assert path.isdir(output_path)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # The new f is the same as the existing a (1.1.1.1.1.1.1.1.1.1.1.1.1.1)
        # The new g is completelly different (3.3.3.5.5.6.6.6.6.6.6.6.6.6)

        assert ["id", "address"] in clusters
        assert ["a", "1.1.1.1.1.1.1.1.1.1.1.1.1.1"] in clusters
        assert ["b", "1.1.1.1.1.2.2.2.2.2.2.2.2.2"] in clusters
        assert ["c", "1.1.1.2.2.3.3.3.3.3.3.3.3.3"] in clusters
        assert ["d", "2.2.2.4.4.5.5.5.5.5.5.5.5.5"] in clusters
        assert ["e", "1.1.1.3.3.4.4.4.4.4.4.4.4.4"] in clusters
        assert ["f", "1.1.1.1.1.1.1.1.1.1.1.1.1.1"] in clusters
        assert ["g", "3.3.3.5.5.6.6.6.6.6.6.6.6.6"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "26,24,22,20,18,16,14,12,10,8,6,4,2,0"
        assert run_json["parameters"]["delimiter"] == "."

        #WIP
        assert len(run_json["threshold_map"]) == 14
        assert run_json["threshold_map"]["0"] == 26.0
        assert run_json["threshold_map"]["1"] == 24.0
        assert run_json["threshold_map"]["2"] == 22.0
        assert run_json["threshold_map"]["3"] == 20.0
        assert run_json["threshold_map"]["4"] == 18.0
        assert run_json["threshold_map"]["5"] == 16.0
        assert run_json["threshold_map"]["6"] == 14.0
        assert run_json["threshold_map"]["7"] == 12.0
        assert run_json["threshold_map"]["8"] == 10.0
        assert run_json["threshold_map"]["9"] == 8.0
        assert run_json["threshold_map"]["10"] == 6.0
        assert run_json["threshold_map"]["11"] == 4.0
        assert run_json["threshold_map"]["12"] == 2.0
        assert run_json["threshold_map"]["13"] == 0.0

    # Thresholds JSON
    thresholds_path = path.join(output_path, "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 14
        assert thresholds_json["0"] == 26.0
        assert thresholds_json["1"] == 24.0
        assert thresholds_json["2"] == 22.0
        assert thresholds_json["3"] == 20.0
        assert thresholds_json["4"] == 18.0
        assert thresholds_json["5"] == 16.0
        assert thresholds_json["6"] == 14.0
        assert thresholds_json["7"] == 12.0
        assert thresholds_json["8"] == 10.0
        assert thresholds_json["9"] == 8.0
        assert thresholds_json["10"] == 6.0
        assert thresholds_json["11"] == 4.0
        assert thresholds_json["12"] == 2.0
        assert thresholds_json["13"] == 0.0

def test_delimiter_slash(tmp_path):
    # "delimiter": "/"
    config = {}

    clusters_path = get_path("data/clusters/basic_delim_slash.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "/"

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
        assert ["A", "1/1/1"] in clusters
        assert ["B", "1/1/2"] in clusters
        assert ["C", "1/1/3"] in clusters
        assert ["D", "1/1/4"] in clusters
        assert ["E", "1/1/2"] in clusters
        assert ["F", "1/1/4"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "5,3,0"
        assert run_json["parameters"]["delimiter"] == "/"

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

def test_delimiter_r(tmp_path):
    # "delimiter": "r"
    config = {}

    clusters_path = get_path("data/clusters/basic_delim_r.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "r"

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
        assert ["A", "1r1r1"] in clusters
        assert ["B", "1r1r2"] in clusters
        assert ["C", "1r1r3"] in clusters
        assert ["D", "1r1r4"] in clusters
        assert ["E", "1r1r2"] in clusters
        assert ["F", "1r1r4"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "5,3,0"
        assert run_json["parameters"]["delimiter"] == "r"

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

def test_delimiter_mismatch_all(tmp_path):
    # The `--delimiter` doesn't match the delimiter used in
    # the cluster file.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "/"

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception

    expected_message = "something went wrong with cluster assignment"
    expected_message += "\ndistance file: " + str(config["dists"])
    expected_message += "\nmembership file: " + str(config["rclusters"])
    expected_message += "\nthreshold map: {0: 5.0, 1: 3.0, 2: 0.0}"
    expected_message += "\nlinkage method: " + str(config["method"])
    expected_message += "\ndelimiter: " + str(config["delimiter"])
    expected_message += "\n\nCheck error messages:"
    expected_message += "\nError: delimiter was not found for samples ['A', 'B', 'C', 'D']."

    assert str(exception.value) == expected_message

    assert path.isfile(path.join(output_path, "results.text")) == False

def test_delimiter_mismatch_some(tmp_path):
    # The `--delimiter` doesn't match the delimiter used in
    # the cluster file.
    config = {}

    clusters_path = get_path("data/clusters/delimiter_mix.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "/"

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception

    expected_message = "something went wrong with cluster assignment"
    expected_message += "\ndistance file: " + str(config["dists"])
    expected_message += "\nmembership file: " + str(config["rclusters"])
    expected_message += "\nthreshold map: {0: 5.0, 1: 3.0, 2: 0.0}"
    expected_message += "\nlinkage method: " + str(config["method"])
    expected_message += "\ndelimiter: " + str(config["delimiter"])
    expected_message += "\n\nCheck error messages:"
    expected_message += "\nError: delimiter was not found for samples ['A', 'C']."

    assert str(exception.value) == expected_message

    assert path.isfile(path.join(output_path, "results.text")) == False

def test_address_errors(tmp_path):
    # Multiple different address errors:
    # - wrong delimiter
    # - wrong length
    # - non-integer
    config = {}

    clusters_path = get_path("data/clusters/address_errors.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception

    expected_message = "something went wrong with cluster assignment"
    expected_message += "\ndistance file: " + str(config["dists"])
    expected_message += "\nmembership file: " + str(config["rclusters"])
    expected_message += "\nthreshold map: {0: 5.0, 1: 3.0, 2: 0.0}"
    expected_message += "\nlinkage method: " + str(config["method"])
    expected_message += "\ndelimiter: " + str(config["delimiter"])
    expected_message += "\n\nCheck error messages:"
    expected_message += "\nError: delimiter was not found for samples ['A']."
    expected_message += "\nError: genomic address length is incorrect for samples ['B', 'C', 'E']; expected length (3) based on thresholds {0: 5.0, 1: 3.0, 2: 0.0}."
    expected_message += "\nError: address could not be converted to an integer for samples ['D']."


    assert str(exception.value) == expected_message

    assert path.isfile(path.join(output_path, "results.text")) == False

def test_delimiter_too_long(tmp_path):
    # Tests a delimiter that's too long (multiple characters).
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = ".."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'please specify a different delimiter {config["delimiter"]} ie. ,|.|\\||-'

def test_delimiter_tab(tmp_path):
    # Tests a tab delimiter.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "\t"

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'please specify a different delimiter {config["delimiter"]} ie. ,|.|\\||-'

def test_delimiter_newline(tmp_path):
    # Tests a newline delimiter.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "\n"

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'please specify a different delimiter {config["delimiter"]} ie. ,|.|\\||-'

def test_no_thresholds(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = None
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'you must specify --thresholds or --threshold_map'

def test_missing_cluster_file(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/file_does_not_exist123.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'{config["rclusters"]} does not exist or is empty'

def test_missing_pairwise_file(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/file_does_not_exist123.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'{config["dists"]} does not exist or is empty'

def test_missing_threshold_file(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = None
    config["thresh_map"] = get_path("data/does_not_exist123.json")

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'{config["thresh_map"]} does not exist or is empty'

def test_linkage_invalid(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "INVALID"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'{config["method"]} is not one of the accepeted methods {CLUSTER_METHODS}'

def test_extension_bad_cluster(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/basic.badext")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f"{clusters_path} does not have a valid extension ['.txt', '.tsv', '.mat', '.text']"

def test_extension_bad_pairwise(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.badext")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f"{pairwise_distances_path} does not have a valid extension ['.txt', '.tsv', '.mat', '.text']"

def test_extensions(tmp_path):
    # Quickly test that accepted extensions are accepted.
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = True

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    # '.txt'
    config["dists"] = get_path("data/pairwise_distances/basic.txt")
    config["rclusters"] = get_path("data/clusters/basic.txt")
    call(config)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address"] in clusters
        assert ["A", "1.1.1"] in clusters
        assert ["B", "1.1.2"] in clusters
        assert ["C", "1.1.3"] in clusters
        assert ["D", "1.1.4"] in clusters
        assert ["E", "1.1.2"] in clusters
        assert ["F", "1.1.4"] in clusters

    # '.tsv'
    config["dists"] = get_path("data/pairwise_distances/basic.tsv")
    config["rclusters"] = get_path("data/clusters/basic.tsv")
    call(config)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address"] in clusters
        assert ["A", "1.1.1"] in clusters
        assert ["B", "1.1.2"] in clusters
        assert ["C", "1.1.3"] in clusters
        assert ["D", "1.1.4"] in clusters
        assert ["E", "1.1.2"] in clusters
        assert ["F", "1.1.4"] in clusters

    # '.mat'
    config["dists"] = get_path("data/pairwise_distances/basic.mat")
    config["rclusters"] = get_path("data/clusters/basic.mat")
    call(config)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address"] in clusters
        assert ["A", "1.1.1"] in clusters
        assert ["B", "1.1.2"] in clusters
        assert ["C", "1.1.3"] in clusters
        assert ["D", "1.1.4"] in clusters
        assert ["E", "1.1.2"] in clusters
        assert ["F", "1.1.4"] in clusters

    # '.text'
    config["dists"] = get_path("data/pairwise_distances/basic.text")
    config["rclusters"] = get_path("data/clusters/basic.text")
    call(config)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address"] in clusters
        assert ["A", "1.1.1"] in clusters
        assert ["B", "1.1.2"] in clusters
        assert ["C", "1.1.3"] in clusters
        assert ["D", "1.1.4"] in clusters
        assert ["E", "1.1.2"] in clusters
        assert ["F", "1.1.4"] in clusters

# pairwise dists file that triggers this but isn't matrix:
#   if len(header) == 3 and num_rows != 3:
# gas call -t 5,3,0 -o test_call -r tests/data/clusters/small.tsv -d tests/data/pairwise_distances/small.tsv --force --delimiter "."

def test_small_3_by_3(tmp_path):
    # A small input that creates a 3x3 pairwise distances file.
    # Previously, such input was being interpretted as square matrix,
    # when it should not have been.
    config = {}

    clusters_path = get_path("data/clusters/small.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/small.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    call(config)

    assert path.isdir(output_path)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        assert ["id", "address"] in clusters
        assert ["A", "1.1.1"] in clusters
        assert ["B", "1.1.2"] in clusters
        assert ["C", "1.1.2"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "5,3,0"
        assert run_json["parameters"]["delimiter"] == "."

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

def test_method_single(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/wikipedia-single.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/wikipedia-single.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "45,40,35,30,25,20,15,10,5,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    call(config)

    assert path.isdir(output_path)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # The new f is the same as the existing d (1.1.1.1.2.4.5.5.5.5)

        assert ["id", "address"] in clusters
        assert ["a", "1.1.1.1.1.1.1.1.1.1"] in clusters
        assert ["b", "1.1.1.1.1.1.2.2.2.2"] in clusters
        assert ["c", "1.1.1.1.1.2.3.3.3.3"] in clusters
        assert ["d", "1.1.1.1.2.4.5.5.5.5"] in clusters
        assert ["e", "1.1.1.1.1.3.4.4.4.4"] in clusters
        assert ["f", "1.1.1.1.2.4.5.5.5.5"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "45,40,35,30,25,20,15,10,5,0"
        assert run_json["parameters"]["delimiter"] == "."

        assert len(run_json["threshold_map"]) == 10
        assert run_json["threshold_map"]["0"] == 45.0
        assert run_json["threshold_map"]["1"] == 40.0
        assert run_json["threshold_map"]["2"] == 35.0
        assert run_json["threshold_map"]["3"] == 30.0
        assert run_json["threshold_map"]["4"] == 25.0
        assert run_json["threshold_map"]["5"] == 20.0
        assert run_json["threshold_map"]["6"] == 15.0
        assert run_json["threshold_map"]["7"] == 10.0
        assert run_json["threshold_map"]["8"] == 5.0
        assert run_json["threshold_map"]["9"] == 0.0

    # Thresholds JSON
    thresholds_path = path.join(output_path, "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 10
        assert thresholds_json["0"] == 45.0
        assert thresholds_json["1"] == 40.0
        assert thresholds_json["2"] == 35.0
        assert thresholds_json["3"] == 30.0
        assert thresholds_json["4"] == 25.0
        assert thresholds_json["5"] == 20.0
        assert thresholds_json["6"] == 15.0
        assert thresholds_json["7"] == 10.0
        assert thresholds_json["8"] == 5.0
        assert thresholds_json["9"] == 0.0

def test_method_average(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/wikipedia-average.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/wikipedia-average.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "45,40,35,30,25,20,15,10,5,0"
    config["thresh_map"] = None

    config["method"] = "average"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    call(config)

    assert path.isdir(output_path)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # The new f is the same as the existing d (1.1.1.2.3.4.5.5.5.5)

        assert ["id", "address"] in clusters
        assert ["a", "1.1.1.1.1.1.1.1.1.1"] in clusters
        assert ["b", "1.1.1.1.1.1.2.2.2.2"] in clusters
        assert ["c", "1.1.1.2.2.3.4.4.4.4"] in clusters
        assert ["d", "1.1.1.2.3.4.5.5.5.5"] in clusters
        assert ["e", "1.1.1.1.1.2.3.3.3.3"] in clusters
        assert ["f", "1.1.1.2.3.4.5.5.5.5"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "average"
        assert run_json["parameters"]["thresholds"] == "45,40,35,30,25,20,15,10,5,0"
        assert run_json["parameters"]["delimiter"] == "."

        assert len(run_json["threshold_map"]) == 10
        assert run_json["threshold_map"]["0"] == 45.0
        assert run_json["threshold_map"]["1"] == 40.0
        assert run_json["threshold_map"]["2"] == 35.0
        assert run_json["threshold_map"]["3"] == 30.0
        assert run_json["threshold_map"]["4"] == 25.0
        assert run_json["threshold_map"]["5"] == 20.0
        assert run_json["threshold_map"]["6"] == 15.0
        assert run_json["threshold_map"]["7"] == 10.0
        assert run_json["threshold_map"]["8"] == 5.0
        assert run_json["threshold_map"]["9"] == 0.0

    # Thresholds JSON
    thresholds_path = path.join(output_path, "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 10
        assert thresholds_json["0"] == 45.0
        assert thresholds_json["1"] == 40.0
        assert thresholds_json["2"] == 35.0
        assert thresholds_json["3"] == 30.0
        assert thresholds_json["4"] == 25.0
        assert thresholds_json["5"] == 20.0
        assert thresholds_json["6"] == 15.0
        assert thresholds_json["7"] == 10.0
        assert thresholds_json["8"] == 5.0
        assert thresholds_json["9"] == 0.0

def test_method_complete(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/wikipedia-complete.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/wikipedia-complete.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "45,40,35,30,25,20,15,10,5,0"
    config["thresh_map"] = None

    config["method"] = "complete"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    call(config)

    assert path.isdir(output_path)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # The new f is the same as the existing d (1.2.2.2.3.4.5.5.5.5)

        assert ["id", "address"] in clusters
        assert ["a", "1.1.1.1.1.1.1.1.1.1"] in clusters
        assert ["b", "1.1.1.1.1.1.2.2.2.2"] in clusters
        assert ["c", "1.2.2.2.2.3.4.4.4.4"] in clusters
        assert ["d", "1.2.2.2.3.4.5.5.5.5"] in clusters
        assert ["e", "1.1.1.1.1.2.3.3.3.3"] in clusters
        assert ["f", "1.2.2.2.3.4.5.5.5.5"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "complete"
        assert run_json["parameters"]["thresholds"] == "45,40,35,30,25,20,15,10,5,0"
        assert run_json["parameters"]["delimiter"] == "."

        assert len(run_json["threshold_map"]) == 10
        assert run_json["threshold_map"]["0"] == 45.0
        assert run_json["threshold_map"]["1"] == 40.0
        assert run_json["threshold_map"]["2"] == 35.0
        assert run_json["threshold_map"]["3"] == 30.0
        assert run_json["threshold_map"]["4"] == 25.0
        assert run_json["threshold_map"]["5"] == 20.0
        assert run_json["threshold_map"]["6"] == 15.0
        assert run_json["threshold_map"]["7"] == 10.0
        assert run_json["threshold_map"]["8"] == 5.0
        assert run_json["threshold_map"]["9"] == 0.0

    # Thresholds JSON
    thresholds_path = path.join(output_path, "thresholds.json")
    assert path.isfile(thresholds_path)
    with open(thresholds_path) as thresholds_file:
        thresholds_json = json.load(thresholds_file)

        assert len(thresholds_json) == 10
        assert thresholds_json["0"] == 45.0
        assert thresholds_json["1"] == 40.0
        assert thresholds_json["2"] == 35.0
        assert thresholds_json["3"] == 30.0
        assert thresholds_json["4"] == 25.0
        assert thresholds_json["5"] == 20.0
        assert thresholds_json["6"] == 15.0
        assert thresholds_json["7"] == 10.0
        assert thresholds_json["8"] == 5.0
        assert thresholds_json["9"] == 0.0

def test_invalid_header_cluster(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/csv.text")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'{clusters_path} does not appear to be a properly TSV-formatted file'

def test_invalid_header_pairwise_distances(tmp_path):
    config = {}

    clusters_path = get_path("data/clusters/small.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/csv.text")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f'{pairwise_distances_path} does not appear to be a properly TSV-formatted file'

def test_batch_size_n1(tmp_path):
    # batch_size = -1
    # The (pairwise distance) input has 2 elements {E, F}
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = -1

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f"batch size ({config['batch_size']}) must be >=1"

    assert path.isdir(output_path) == False

def test_batch_size_0(tmp_path):
    # batch_size = 0
    # The (pairwise distance) input has 2 elements {E, F}
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 0

    with pytest.raises(Exception) as exception:
        call(config)

    assert exception.type == Exception
    assert str(exception.value) == f"batch size ({config['batch_size']}) must be >=1"

    assert path.isdir(output_path) == False

def test_batch_size_1(tmp_path):
    # batch_size = 1
    # The (pairwise distance) input has 2 elements {E, F}
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 1

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
        assert run_json["parameters"]["delimiter"] == "."

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

def test_batch_size_2(tmp_path):
    # batch_size = 2
    # The (pairwise distance) input has 2 elements {E, F}
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 2

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
        assert run_json["parameters"]["delimiter"] == "."

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

def test_batch_size_3(tmp_path):
    # batch_size = 3
    # The (pairwise distance) input has 2 elements {E, F}
    config = {}

    clusters_path = get_path("data/clusters/basic.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/basic.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 3

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
        assert run_json["parameters"]["delimiter"] == "."

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

def test_batch_double_digit(tmp_path):
    # Tests when working with double digit addresses.
    config = {}

    clusters_path = get_path("data/clusters/double_digit.tsv")
    pairwise_distances_path = get_path("data/pairwise_distances/double_digit.tsv")
    output_path = path.join(tmp_path, "test_out")

    config["dists"] = pairwise_distances_path
    config["rclusters"] = clusters_path
    config["outdir"] = output_path
    config["force"] = False

    config["thresholds"] = "5,3,0"
    config["thresh_map"] = None

    config["method"] = "single"

    config["sample_col"] = "id"
    config["address_col"] = "address"
    config["delimiter"] = "."

    config["batch_size"] = 100

    call(config)

    assert path.isdir(output_path)

    # Clusters
    clusters_path = path.join(output_path, "results.text")
    assert path.isfile(clusters_path)
    with open(clusters_path) as clusters_file:
        clusters = csv.reader(clusters_file, delimiter="\t")

        # The new N is completely different from everything.
        # All are unique:
        assert ["id", "address"] in clusters
        assert ["A", "1.1.1"] in clusters
        assert ["B", "2.2.2"] in clusters
        assert ["C", "3.3.3"] in clusters
        assert ["D", "4.4.4"] in clusters
        assert ["E", "5.5.5"] in clusters
        assert ["F", "6.6.6"] in clusters
        assert ["G", "7.7.7"] in clusters
        assert ["H", "8.8.8"] in clusters
        assert ["I", "9.9.9"] in clusters
        assert ["J", "10.10.10"] in clusters
        assert ["K", "11.11.11"] in clusters
        assert ["L", "12.12.12"] in clusters
        assert ["M", "13.13.13"] in clusters
        assert ["N", "14.14.14"] in clusters

    # Run JSON
    run_path = path.join(output_path, "run.json")
    assert path.isfile(run_path)
    with open(run_path) as run_file:
        run_json = json.load(run_file)

        assert run_json["parameters"]["method"] == "single"
        assert run_json["parameters"]["thresholds"] == "5,3,0"
        assert run_json["parameters"]["delimiter"] == "."

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
