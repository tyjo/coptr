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

Select only high-quality assemblies
-----------------------------------
Metagenomic assembly is a challenging problem, and assemblies are not guaranteed
to be correct. We set some basic criteria for using an assembly:

* It is estimated to be >90% complete
* It is estimated to have <5% contamination
* It has an N50 > 10kb
* It has an average contig length > 10kb
* It has fewer than 400 contigs

The first two criteria are similar to the MIMAG criteria for high-quality
assemblies.

When choosing between an assembly, and complete reference genome, preference
should always be given to the complete genome.

Nomenclature for reference genomes
----------------------------------
coPTR needs to be able to group contigs by their assembly. It is recommended
to choose a consistent naming convention for all contigs. For example,
``REF-ID-FOR-ASSEMBLY-1:CONTIG-ID-1``. When mapping reads against the reference
database (see below), you can provide an argument that maps the id of the contig to
its reference genome.


Indexing the reference database
===============================

**On the example data:**

.. code-block:: text

    $ python coptr.py index example-data/ref-db example-data/ref-db/example-index
    [INFO] (Aug 31, 2020 12:09:50) ReadMapper: bowtie2-build example-data/ref-db/e-coli-mag.fna,example-data/ref-db/l-gasseri-ref.fa, example-data/example-index
    ...bowtie2 output...
    [INFO] (Aug 31, 2020 12:09:53) ReadMapper: indexed 2 fasta files for reference database.



coPTR provides a few utilities to make indexing the database with bowtie2
easier. Indexing a database of fasta files can be accomplished by calling
``coptr.py index``:

.. code-block:: text

    usage: coptr.py index [-h] ref-fasta index-out

    positional arguments:
      ref_fasta   File or folder containing fasta to index. If a folder, the
                  extension for each fasta must be one of [.fasta, .fna, .fa]
      index_out   Filepath to store index.

    optional arguments:
      -h, --help  show this help message and exit

It takes two arguments. The first ``ref-fasta`` is either a folder containing
fasta for the database. If it is a folder, coPTR will scan the folder for
all files ending in ``.fasta``, ``.fna``, or ``.fa``, and call ``bowtie2 index``
on these files. If it is a file, coPTR will index the single file.

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


Note on paired end sequencing
-----------------------------
coPTR uses the density of reads along the genome to estimate PTRs. It
uses the starting coordinate at each read to fit its model. Because
mate pairs are not independent, once one read of the mate pair is observed
the second read adds little information.

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
expression that matches the reference genome id from a sequence id. For example,
consider the sequence id

``REF-ID-FOR-ASSEMBLY-1:CONTIG-ID-1``

The sequence id for the fasta would look something like

``>REF-ID-FOR-ASSEMBLY-1:CONTIG-ID-1 [other data that is not parsed]``

and the reference genome id is ``REF-ID-FOR-ASSEMBLY-1``.

We would want to the regular expression to match ``REF-ID-FOR-ASSEMBLY-1``. This
is straightforward if you use a single special character (such as :) to separate
reference genome ids from contig ids. In this case, the right regex is

``[\w-]+``

which matches one or more word characters or "-" characters in a row. If no
matches are found, coPTR defaults to using the sequence id itself as the reference
genome id.

The optional flag ``--check-regex`` allows you to sanity check your regular
expressions. It outputs the total number of reference genomes found and
their reference ids. On the example data:

.. code-block:: text

    $ python coptr.py extract example-data/bam example-data/coverage-maps --check-regex
    [INFO] (Aug 31, 2020 12:25:32) BamProcessor: found 190 reference sequences corresponding to 2 genomes
    [INFO] (Aug 31, 2020 12:25:32) BamProcessor: reference genome ids:
         ERS235517|65
         NC_008530.1



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