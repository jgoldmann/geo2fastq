geo2fastq & geo2trackhub
========================

Search the GEO database for a specific accession number. With the `-d` tag, the data will be downloaded and converted to fastq. `geo2trackhub` will align the sequencing experiments and generate a trackhub (that can be displayed in a genome browser).


## Dependencies

Python Modules (should be automatically installed during installation):

* biopython
* PyYaml
* parallel python (pp)
* trackhub
* fabric

To create fastqs:

* sra tools (fastq-dump and vdb-validate). Available from the [NCBI website](http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software).

Aligning and Mapping:

* These scripts are specialized to work with the [soladmin](https://bitbucket.org/simonvh/soladmin) experiment management system.

For easy bigwig generation:

* [bam2bw](https://bitbucket.org/simonvh/bam2bw)


## Install

Download and run `python setup.py install`. Tests can be run with `python setup.py test`.


## Run

Usage: 

* `geo2fastq <search term>` to search the database.
* `geo2fastq -d <GSE_Accession>` to download and convert data to `fastq` files.
* `geo2trackhub -d <GSE_Accession>` to do all of the above plus aligning the files and generating a trackhub file.

