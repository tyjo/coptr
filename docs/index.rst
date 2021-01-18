.. coPTR documentation master file, created by
   sphinx-quickstart on Mon Aug 24 10:35:15 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================================================
CoPTR: *Co*\ mpute *P*\ eak-to-*T*\ rough *R*\ atio
====================================================
The peak-to-trough ratio is the ratio of sequencing coverage near
the replication origin to the replication terminus in a species of
bacteria. PTRs are a valuable tool for investigating microbiome dynamics
because they reflect a species' replication rate at the time of sampling.

CoPTR is a tool for computing peak-to-trough ratios from metagenomic
sequencing data from complete reference genomes and draft assemblies.
CoPTR takes as input fastq files from multiple metagenomic samples
and a reference database. It outputs the log2 PTR for each species
in a sample.

`View on Github <https://github.com/tyjo/coptr>`_

.. toctree::
   :maxdepth: 1

   quick-start
   installation
   databases
   tutorial
   modules