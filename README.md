geo2fastq
=========

Search the GEO database for a specific accession number. With the `-d` tag, the data will be downloaded and converted to fastq.


## Dependencies

Python Modules:

* biopython
* PyYaml
* parallel python

To create fastqs:

* sra tools (fastq-dump and vdb-validate). Available from the [NCBI website](http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software).

For easy bigwig generation:

* [bam2bw](https://bitbucket.org/simonvh/bam2bw)

## Install

Download and run `python setup.py install`. Tests can be run with `python setup.py test`.


## Run

Usage: 

* `geo2fastq <search term>` to search the database.
* `geo2fastq -d <GSE_Accession>` to download and convert data.

