===================
Reference Databases
===================


IGGdb High-Quality Gut Genomes
------------------------------

* Pre-computed bowtie2 indexes
* Representative genomes for human gut species with a high-quality genome (N=2,935)
* **Citation:** Nayfach, S, et al. "New insights from uncultivated genomes of the global human gut microbiome." Nature 568.7753 (2019): 505-510.
* **Download Link (11.2 GB):** `<https://dl.dropboxusercontent.com/s/9fq2vck1dugmic5/IGG_v1.0_split.tar.gz>`_


IGGdb High-Quality Genomes Across Environments
----------------------------------------------

This database is too large to host bowtie2 indexes. However, you
can download it directly and index with CoPTR. We recommended breaking
the database into several smaller indexes, such that each index contains 
fewer than 4 billion bases. This allows bowtie2 to create
a small index for each subset. You can merge reads mapped to multiple 
indexes with the ``merge`` command.

* Representative genomes for species with a high-quality genome (N=16,136)
* Download link available from `<https://github.com/snayfach/IGGdb>`_
* **Citation:** Nayfach, S, et al. "New insights from uncultivated genomes of the global human gut microbiome." Nature 568.7753 (2019): 505-510.