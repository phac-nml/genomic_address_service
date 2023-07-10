[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/genomic_address_service/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/genomic_address_service)
[![Conda](https://img.shields.io/conda/dn/bioconda/profile_dists?color=green)](https://anaconda.org/bioconda/genomic_address_service)
[![License: Apache-2.0](https://img.shields.io/github/license/phac-nml/genomic_address_service)](https://www.apache.org/licenses/LICENSE-2.0)


## Profile Dists

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Quick Start](#quick-start)
- [FAQ](#faq)
- [Citation](#citation)
- [Legal](#legal)
- [Contact](#contact)

## Introduction

Surveillance and outbreak analysis of pathogens has been operationalized by multiple public health laboratories around 
the world using gene-by-gene, SNP and k-mer approaches to produce estimates of genetic distance between sets of samples.
Standard phylogentic or heirarchal clustering approaches group samples based on distances but need to be recalculated 
everytime a new sample is added. The lack of repeatable genetic units or nomenclature between runs of clustering makes
communication between different groups difficult and the slow due to the poor scaling of these approaches to larger
sample sizes. Enterobase implements [HeirCC](https://github.com/zheminzhou/pHierCC) as an algorithm approach to assign 
new cgMLST profiles into existing multi-level cluster nomenclature from a minimum spanning tree, using single linkage in real time, 
along with some tools for evaluating the clusters when run in de novo mode. [SnapperDB](https://github.com/ukhsa-collaboration/snapperdb) 
used by the UKHSA utilizes single linkage clustering based on SNPs to produce a multi-level cluster nomenclature for outbreak and surveillance activities. 
[ReportTree](https://github.com/insapathogenomics/ReporTreer) will perform de novo clustering based on [SciPy](https://scipy.org/)
linkage methods from a sequence alignment, VCF, allele profile, or distance matrix and so provides significant flexibility in inputs
compared to the other two methods. Like Enterobase, ReportTree also provides tools for evaluating regions of cluster stability and
can maintain cluster nomenclature between runs with the caveate that Report tree has the potential to split and merge previous cluster assignments.
This poses a significant issue in terms of Public Health Reporting and it is desirable to maintain existing cluster designations with
occasional updates to the cluster assignments which must then be reported to all clients and partners to make the changes in
their databases as well. HeirCC and be run as a standalone piece of software and applied to new cgMLST schemes and organisms readily, whereas
SnapperDB is a relatively complex pipeline where establishing databases for new organisms is rather complex. HeirCC could
be incorporated into other bioinformatics workflows relatively easily however, SnapperDB cannot.

<br><br>
While HeirCC provides customization and flexibility in terms of thresholds and schemes, it is limited to single-linkage
clustering which is known to have issues with stability and be prone to producing "scraggly" clusters due to chaining.
ReportTree provides customization and reporting functionality which can be of use to public health labs, however the 
potential for changing cluster assignments poses issues for distributed Public Health surveillance and outbreak activities
due to the potential for becoming out of synch with each other. Within our public health partners, there is a need
for a clustering service which can perform de novo clustering based on average, complete, and single linkage which can
then be partitioned into clusters based on multiple thresholds.  Additionally, there is a need to assign new samples into
an existing clustering to provide stable nomenclature for communication between different partners and stakeholders. Any
changes to the clustering assignments needs to be done purposefully with changes clearly communicated to partners and stakeholders.
<br><br>
To address needs of users within our team we have designed an integrated solution for calculating distance matricies and querying of genetically
similar samples within a defined threshold to support outbreak and surveillance activities. We provide the flexibility to have standard text 
based outputs as well as parquet. It is implemented in pure python and currently is only available in a single threaded version but later
refinements may include the support for multiprocessing. To facilitate integration of the tool into larger workflows it will also be implemented 
as a nextflow workflow.

As datasets grow, using text based formats such as CSV, TSV represent significant amounts of runtime in terms
of reading, parsing and writing.  New formats such as [parquet](https://parquet.apache.org/) support compression and are optimized for 
efficiency of storage and retreiveal of data.



## Installation

Install the latest released version from conda:

        conda create -c bioconda -c conda-forge -n genomic_address_service genomic_address_service

Install using pip:

        pip install genomic_address_service

Install the latest master branch version directly from Github:

        pip install git+https://github.com/phac-nml/genomic_address_service.git



## Usage
If you run ``gas``, you should see the following usage statement:

    Usage: gas <command> [options] <required arguments>

    To get minimal usage for a command use:
    gas command

    To get full help for a command use one of:
    gas command -h
    cladeomatic command --help


    Available commands:

    mcluster  De novo nested multi-level clustering
    call      Call genomic address based on existing clusterings
    test      Test functionality on a small dataset

Supported distance matrix formats
=====
**Square**

|  id  |  S1  |  S2  |  S3  |  S4 |  S5 |  S6  | 
| ----------- | ----------- |----------- | ----------- | ----------- |----------- | ----------- | 
|  S1  |	0  |  0  |  3  |  3  |  9  |  9  |
|  S2  |	0  |  0  |  3  |  3  |  9  |  9  | 
|  S3  |	3  |  3  |  0  |  0  |  9  |  9  | 
|  S4  |	3  |  3  |  0  |  0  |  9  |  9  | 
|  S5  |	9  |  9  |  9  |  9  |  0  |  0  | 
|  S6  |	9  |  9 |  9  |  9  |  0  |  0 |  

- Distance matrix units can be of float, or integer type with the constrain that the diagnonal must be 0 and the first
line must be a header with all of the samples

Quick start
=====
**De novo Multi-level Clustering**

Mcluster minimally accepts as input a distance matrix, output directory, and a set of thresholds delimeted by a comma.

<br />

        gas mcluster --matrix distance.matrix.text --outdir results --thresholds 10,5,0

This will produce a cluster file with the following header [id, nomenclature, level_1, level_2, level_3]

<br />

   |  id  |  nomenclature  |  level_1  |  level_2  |  level_3  |
| ----------- | ----------- | ----------- |  ----------- |   ----------- |
|  S1  |  1.1.1  |  1  |  1  |  1  |
|  S2  |  1.1.1  |  1  |  1  |  1  |
|  S3  |  1.1.2  |  1  |  1  |  2  |
|  S4  |  1.1.2  |  1  |  1  |  2  |
|  S5  |  1.2.3  |  1  |  2  |  3  |
|  S6  |  1.2.3  |  1  |  2  |  3  |

**Outputs:**

```
{Output folder name}
├── distances.{text|parquet} - Three column file of [query_id, ref_if, distance]
├── thresholds.json - JSON formated mapping of columns to distance thresholds
├── clusters.{text|parquet} - Either symmetric distance matrix or three column file of [query_id, ref_if, distance]
├── tree.newick - Newick formatted dendrogram of the linkage matrix produced by SciPy
└── run.json - Contains logging information for the run including parameters, newick tree, and threshold mapping info
```

**Cluster assignment**

Coming soon


## Benchmarks

Coming soon

## FAQ

Coming soon

## Citation

Robertson, James, Wells, Matthew, Schonfeld, Justin, Reimer, Aleisha. Genomic Address Service: Convenient package for de novo clustering and sample assignment to existing clusters. 2023. https://github.com/phac-nml/genomic_address_service

## Legal

Copyright Government of Canada 2023

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.


## Contact

**James Robertson**: james.robertson@phac-aspc.gc.ca
