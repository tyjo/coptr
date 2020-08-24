.. coPTR documentation master file, created by
   sphinx-quickstart on Mon Aug 24 10:35:15 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for coPTR
=================================

============
Quick Start
============
1. Index the reference database:

    ``coptr.py index example-data/ref-db example-data/ref-db/example-index``

2. Map reads:

    ``coptr.py map example-data/ref-db/example-index example-data/ example-data/bam``

3. Extract read positions:

    ``coptr.py extract example-data/bam example-data/coverage-maps``

4. Estimate PTRs:

    ``coptr.py estimate coverage-maps``

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules