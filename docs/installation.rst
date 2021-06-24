============
Installation
============

The recommended installation is using the conda environment manager.
You can find installation instructions for your operating system
`here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.

First, checkout the the code from Github

    .. code-block::

      bash$ git clone https://github.com/tyjo/coptr
      bash$ cd coptr/

You can then create a new environment to run coPTR:

    .. code-block::

      # creates an environment called coptr
      $ conda env create -f coptr.yml

      # you only need to create this once
      # to activate the environment:
      $ conda activate coptr
      $ pip install .

To check if you can run CoPTR:

    .. code-block::

      $ coptr

      usage: coptr.py <command> [options]

      command: index            create a bowtie2 index for a reference database
               map              map reads against a reference database
               merge            merge BAM files from reads mapped to multiple indexes
               extract          compute coverage maps from bam files
               estimate         estimate PTRs from coverage maps
               count            compute read counts for each genome after filtering

      CoPTR (v1.0.0): Compute PTRs from complete reference genomes and assemblies.

      positional arguments:
        command     Command to run.

      optional arguments:
        -h, --help  show this help message and exit


Dependencies
------------
Installation using conda iswill install all required dependencies. However,
CoPTR's dependencies are listed below if you wish to install them manually.

* bowtie2 (>=2.4.1): CoPTR assumes bowtie2 is install and accessible from your PATH variable
* python  (>=3.7.8)
* matplotlib (>=3.3.2)
* numpy (=1.19.1)
* scipy (=1.5.2)
* pysam (=0.16.0.1)