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
      $ conda env create -f coptr-env.yml

      # you only need to create this once
      # to activate the environment:
      $ conda activate coptr

To check if you can run coPTR:

    .. code-block::

      $ python coptr.py -h

      usage: coptr.py <command> [options]

      command: index            create a bowtie2 index for a reference database
               map              map reads against a reference database
               extract          compute coverage maps from bam files
               estimate         estimate PTRs from coverage maps

      Compute PTRs from complete reference genomes and assemblies.

      positional arguments:
        command     Command to run.

      optional arguments:
        -h, --help  show this help message and exit


Dependencies
------------
Installation using conda will install all required dependencies. However,
coPTR's dependencies are listed below if you wish to install them manually.

* bowtie2 (>=2.4.1): coPTR assumes bowtie2 is install and accessible from your PATH variable
* python  (>=3.7.8)
* numpy (=1.19.1)
* scipy (=1.5.2)
* pysam (=0.16.0.1)