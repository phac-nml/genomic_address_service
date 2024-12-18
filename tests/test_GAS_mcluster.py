#!/usr/bin/env python
import os
import tempfile
from genomic_address_service.mcluster import write_clusters  # Adjust the import path based on your project structure

def test_write_clusters():
    # Create mock cluster data
    mock_clusters = {
        '1': ['1', '1', '1'],
        '2': ['1', '1', '2'],
        '3': ['1', '2', '3']
    }
    num_thresholds = 3
    delimiter = "."

    # Create a temporary file
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    try:
        # Write mock clusters to the temporary file
        write_clusters(mock_clusters, num_thresholds, temp_file.name, delimiter)

        # Verify the contents of the file
        with open(temp_file.name, 'r') as file:
            lines = file.readlines()
            # Check the header
            assert lines[0].strip() == "id\taddress\tlevel_1\tlevel_2\tlevel_3"
            # Check the first line of data
            assert lines[1].strip() == "1\t1.1.1\t1\t1\t1"
            assert lines[2].strip() == "2\t1.1.2\t1\t1\t2"
            assert lines[3].strip() == "3\t1.2.3\t1\t2\t3"
    finally:
        # Clean up - delete the temporary file
        os.remove(temp_file.name)