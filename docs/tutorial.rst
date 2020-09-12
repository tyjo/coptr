==============
coPTR Tutorial
==============

This page contains a worked example for using coPTR. It takes you
through the steps in the :doc:`Quick Start <quick-start>`.
The data used in this example is available here:

``https://drive.google.com/file/d/1vGAfTDRDKN1wPULRd4EibxzbYFbJofYZ/view?usp=sharing``

.. contents::
    :depth: 2


Building the reference database
===============================

The first step for running coPTR is to create a reference database of
fasta files. coPTR is designed to run on both complete reference genomes
and draft assemblies, but there are a few considerations for database 
construction.

Cluster genomes at 95% ANI
--------------------------
coPTR estimates PTRs at the species level, estimates among strains within 95%
ANI are consistent. This means you should cluster genomes by 95% ANI, an
operational defintion of a species, and select one representative genome per
species.

Select only high-quality assemblies or complete genomes
-------------------------------------------------------
Metagenomic assembly is a challenging problem, and assemblies are not guaranteed
to be correct. We set some basic criteria for using an assembly:

* It is estimated to be >90% complete
* It is estimated to have <5% contamination
* It has an N50 > 10kb
* It has an average contig length > 10kb

The first two criteria are similar to the MIMAG criteria for high-quality
assemblies.

When choosing between an assembly, and complete reference genome, preference
should always be given to the complete genome.



Indexing the reference database
===============================

**On the example data:**

.. code-block:: text

    $ python coptr.py index example-data/ref-db example-data/ref-db/example-index
    [INFO] (Sep 10, 2020 18:06:41) ReadMapper: found 2 files totaling 0.006 Gb
    [INFO] (Sep 10, 2020 18:06:41) ReadMapper: copying to fasta files to coptr-fna-1599775601.629647.fna with prepended genome ids (filenames)
    [INFO] (Sep 10, 2020 18:06:41) ReadMapper: writing 2 reference genome ids to example-data/ref-db/example-index.genomes
    [INFO] (Sep 10, 2020 18:06:41) ReadMapper: bowtie2-build coptr-fna-1599775601.629647.fna example-data/ref-db/example-index --noref --threads 1
    ...bowtie2 output...
    [INFO] (Sep 10, 2020 18:06:44) ReadMapper: indexed 2 fasta files for reference database.
    [INFO] (Sep 10, 2020 18:06:44) ReadMapper: cleaning up coptr-fna-1599775601.629647.fna




coPTR provides a few utilities to make indexing the database with bowtie2
easier. It assumes the reference database is provided by a folder of fasta
files, one per reference genome.

Indexing a database of fasta files can be accomplished by calling
``coptr.py index``:

.. code-block:: text

    usage: coptr.py index [-h] [--bt2-bmax BT2_BMAX] [--bt2-dcv BT2_DCV] [--bt2-threads BT2_THREADS] [--bt2-packed] ref-fasta index-out

    positional arguments:
      ref_fasta             File or folder containing fasta to index. If a folder,
                            the extension for each fasta must be one of [.fasta,
                            .fna, .fa]
      index_out             Filepath to store index.

    optional arguments:
      -h, --help            show this help message and exit
      --bt2-bmax BT2_BMAX   Set the --bmax arguement for bowtie2-build. Used to
                            control memory useage.
      --bt2-dcv BT2_DCV     Set the --dcv argument for bowtie2-build. Used to
                            control memory usage.
      --bt2-threads BT2_THREADS
                            Number of threads to pass to bowtie2-build.
      --bt2-packed          Set the --packed flag for bowtie2-build. Used to
                            control memory usage.

It takes two arguments. The first ``ref-fasta`` is either a folder containing
fasta for the database. If it is a folder, coPTR will scan the folder for
all files ending in ``.fasta``, ``.fna``, or ``.fa``. coPTR will combine these
into a single fasta file, prepending the filename to each sequence id in order
to track contigs from the same reference genome. It then calls the ```bowtie2-build```
to index this file.

coPTR additionally outputs an ```index_name.genomes``` file with a list of ids for each
reference genome.

Notes on memory usage
---------------------
For large databases, the options ```--bt2-bmax```, ```--bt2-dcv```, and ```--bt2-packed```
can reduce memory usage. Please see the `bowtie2 manual <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer>`_ for more details.
If these options still do not help, you can split the reference database into multiple
parts, and index each separately. The resulting BAM files can be merged with ```coptr.py merge```.

Mapping reads
=============

**On the example data:**

.. code-block:: text

    $ python coptr.py map example-data/ref-db/example-index example-data/fastq example-data/bam
    [INFO] (Aug 31, 2020 12:12:10) ReadMapper: mapping example-data/fastq/ERR969281.fastq.gz to example-data/bam/ERR969281.sam
    [INFO] (Aug 31, 2020 12:12:10) ReadMapper: bowtie2 -x example-data/ref-db/example-index example-data/fastq/ERR969281.fastq.gz --no-unal -p 1
    10818 reads; of these:
      10818 (100.00%) were unpaired; of these:
        4071 (37.63%) aligned 0 times
        6709 (62.02%) aligned exactly 1 time
        38 (0.35%) aligned >1 times
    62.37% overall alignment rate
    [INFO] (Aug 31, 2020 12:12:11) ReadMapper: converting example-data/bam/ERR969281.sam to example-data/bam/ERR969281.bam
    [INFO] (Aug 31, 2020 12:12:11) ReadMapper: cleaning up example-data/bam/ERR969281.sam
    ....
    [INFO] (Aug 31, 2020 12:12:24) ReadMapper: converting example-data/bam/ERR969285.sam to example-data/bam/ERR969285.bam
    [INFO] (Aug 31, 2020 12:12:24) ReadMapper: cleaning up example-data/bam/ERR969285.sam



Once you have indexed a reference database. You can then map reads against
the database. coPTR provides a wrapper around bowtie2 to make read mapping
convenient:

.. code-block:: text

    usage: coptr.py map [-h] [--threads INT] index input out-folder

    positional arguments:
      index              Name of database index.
      input              File or folder containing fastq reads to map. If a
                         folder, the extension for each fastq must be one of
                         [.fastq, .fq, .fastq.gz, fq.gz]
      out_folder         Folder to save mapped reads. BAM files are output here.

    optional arguments:
      -h, --help         show this help message and exit
      --threads THREADS  Number of threads for bowtie2 mapping.

The name of the database index corresponds to the name used from ``coptr.py index``.

The input can either be a single fastq file, or a folder of fastq files to map.
It also takes an optional ``--threads`` argument that allows bowtie2 to use
multiple threads. Reads are output as ``bam`` files to save space.


Merging mapped reads from multiple indexes
------------------------------------------
For large reference databases, it is sometimes necessary to create several
indexes for subsets of the data and map reads against each index. Results
from each index need to be merged to select reads with the best MAPQ across
indexes. You can use ```coptr.py merge``` to merge multiple bam files.

.. code-block:: text

    usage: coptr.py merge [-h] in-bam1 in-bam2 ... in-bamN out-bam

    positional arguments:
      in-bams     A space separateed list of BAM files to merge. Assumes same
                  reads were mapped against different indexes.
      out-bam     Path to merged BAM.

    optional arguments:
      -h, --help  show this help message and exit


Note on paired end sequencing
-----------------------------
coPTR uses the density of reads along the genome to estimate PTRs. It
uses the starting coordinate at each read to fit its model. Because
mate pairs are not independent, once one read of the mate pair is observed
the second read does not add any additional information.

Therefore, we recommend using only **one mate pair from paired end sequencing**.
The ``map`` command has been designed with this in mind.


Computing coverage maps
=======================
**On the example data:**

.. code-block:: text

    $ python coptr.py extract example-data/bam example-data/coverage-maps
    [INFO] (Aug 31, 2020 12:13:53) BamProcessor: processing example-data/bam/ERR969428.bam
    ...
    [INFO] (Aug 31, 2020 12:13:56) BamProcessor: processing example-data/bam/ERR969285.bam
    [INFO] (Aug 31, 2020 12:13:56) BamProcessor: found 190 reference sequences corresponding to 2 genomes

Once reads have been mapped, the next step is to compute the coverage along
each reference genome. In this step, starting positions of each read are
extracted from each bam file, and reads from different contigs of the same
assembly are collected.

.. code-block:: text

    usage: usage: coptr.py extract [-h] [--ref-genome-regex REF_GENOME_REGEX] [--check-regex]
                    in-folder out-folder

    positional arguments:
      in_folder             Folder with BAM files.
      out_folder            Folder to store coverage maps.

    optional arguments:
      -h, --help            show this help message and exit
      --ref-genome-regex REF_GENOME_REGEX
                            Regular expression extracting a reference genome id
                            from the sequence id in a bam file.
      --check-regex         Check the regular expression by counting reference
                            genomes without processing

The important argument here is the ``--ref-genome-regex``. This is a regular
expression that extracts the reference genome id from a sequence id. The default
argument will work with the index created by ```coptr.py index```, and works by
prepending the name of the fasta file, and special character ```|``` to each
sequence id.


Estimating PTRs
===============

**On the example data:**

.. code-block:: text

    # python coptr.py estimate example-data/coverage-maps out --min-reads 2500
    [INFO] (Aug 31, 2020 13:49:06) CoPTRRef: estimating PTRs for NC_008530.1

The final stage is to estimate PTR ratios from coverage maps. This is accomplished
with the ``estimate`` command.

.. code-block:: text

    usage: coptr.py estimate [-h] [--min-reads MIN_READS] [--min-cov MIN_COV] [--threads THREADS] coverage-map-folder out-file

    positional arguments:
      coverage_map_folder   Folder with coverage maps computed from 'extract'.
      out_file              Filename to store PTR table.

    optional arguments:
      -h, --help            show this help message and exit
      --min-reads MIN_READS
                            Minimum number of reads required to compute a PTR
                            (default 5000).
      --min-cov MIN_COV     Fraction of nonzero 10Kb bins required to compute a
                            PTR (default 0.75).
      --min-samples MIN_SAMPLES
                            CoPTRContig only. Minimum number of samples required
                            to reorder bins (default 5).
      --threads THREADS     Number of threads to use (default 1).

We have tried to set sensible default parameters for PTR estimatation. We set
the minimum number of reads for the example data to 2500 in order to keep the
size of the example data small.

The output is a CSV file where, the rows are reference genomes, and the
columns are samples.