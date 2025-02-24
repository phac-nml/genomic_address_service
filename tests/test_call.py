"""
A simple test for running gas call

"""

import pytest
import pathlib as p
from genomic_address_service.call import run_call 




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
    run_call(config)
    outdir = p.Path(config["outdir"])
    output_file = outdir / "results.text"
    assert output_file.read_text().split("\n").sort() == input_results.read_text().split("\n").sort()