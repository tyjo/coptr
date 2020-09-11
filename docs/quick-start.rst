===========
Quick Start
===========

For instructions on installation and environment set up, please see :doc:`Installation <installation>`. 

| Example data:
| ``https://drive.google.com/file/d/1vGAfTDRDKN1wPULRd4EibxzbYFbJofYZ/view?usp=sharing``

    .. code-block::

      git clone https://github.com/tyjo/coptr
      cd coptr

      # Set up the environment:
      conda env create -f coptr-env.yml
      conda activate coptr
      
      # Index the reference database:
      python coptr.py index example-data/ref-db example-data/ref-db/example-index
      
      # Map reads against the reference:
      python coptr.py map example-data/ref-db/example-index example-data/fastq example-data/bam
      
      # Extract read positions:
      python coptr.py extract example-data/bam example-data/coverage-maps
      
      # Estimate ptrs:
      python coptr.py estimate example-data/coverage-maps out --min-reads 2500

      # View the output:
      cat out

See the :doc:`Tutorial <tutorial>` for more details.