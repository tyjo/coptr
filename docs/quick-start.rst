===========
Quick Start
===========

Environment set up and installation: :doc:`Installation <installation>`.

Detailed instructions: :doc:`Tutorial <tutorial>`.

| Example data:
| ``https://dl.dropboxusercontent.com/s/wrsxjimr6l96lcq/example-data.tar.gz``

    .. code-block::

      git clone https://github.com/tyjo/coptr
      cd coptr

      wget https://dl.dropboxusercontent.com/s/wrsxjimr6l96lcq/example-data.tar.gz
      tar -xzvf example-data.tar.gz

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
      # Note the min-reads flag is optional. We recommend the default setting (5000 reads).
      python coptr.py estimate example-data/coverage-maps out.csv --min-reads 2500

      # View the output:
      cat out.csv

      log2(PTR):genome_id/sample_id,ERR969281,ERR969282,ERR969283,ERR969285,ERR969286,ERR969428,ERR969429,ERR969430
      l-gasseri-ref,,,,,,1.1840987863325785,1.1945539660363145,1.2879271469720541
      e-coli-mag,1.2078623467253466,1.0375575947553943,0.9433005522894075,0.759132363901812,0.7846476652840171,,,

