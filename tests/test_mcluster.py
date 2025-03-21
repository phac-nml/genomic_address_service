import pytest
import pathlib
from os import path
from genomic_address_service.mcluster import mcluster, parse_args

def get_path(location):
    directory = path.dirname(path.abspath(__file__))
    return path.join(directory, location)

def test_run_mcluster(tmp_path):
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

    run_path = path.join(args["outdir"], "run.json")
    assert path.isfile(run_path)

    thresholds_path = path.join(args["outdir"], "thresholds.json")
    assert path.isfile(thresholds_path)

    tree_path = path.join(args["outdir"], "tree.nwk")
    assert path.isfile(tree_path)
