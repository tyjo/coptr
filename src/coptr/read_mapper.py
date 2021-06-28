"""
read_mapper.py
======================
Module to map reads to a reference database of
high-quality assemblies and complete reference genomes.
"""

"""
This file is part of CoPTR.

CoPTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CoPTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CoPTR.  If not, see <https://www.gnu.org/licenses/>.
"""

import logging
import os
import os.path
import subprocess as sub
from datetime import datetime, timezone
from pathlib import Path

import pysam

from .util import get_fastq_name


logger = logging.getLogger(__name__)


class ReadMapper:
    """Wrapper around bowtie2."""

    def index(self, ref_fasta, index_out, bt2_bmax, bt2_dcv, bt2_threads, bt2_packed):
        """Build a bowtie2 index from ref_fasta.

        Parameters
        ----------
            ref_fasta : str
                Fasta file or folder containing fasta files to index.
                Valid extensions include '.fasta', '.fna', '.fa'
            index_out : str
                Path to output the index.
            bt2_bmax : str
                Parameter to pass to bowtie2-build --bmax argument.
            bt2_dcv : str
                Parameter to pass to bowtie2-build --dcv argument.
            bt2_threads : str
                Parameter to pass to bowtie2-build --threads argument.
            bt2_packed : str
                Parameter to pass to bowtie2-build --packed argument.
        """
        files_found = 0
        total_size = 0
        ref_files = []
        if os.path.isfile(ref_fasta):
            ref_files.append(ref_fasta)
            files_found += 1
            total_size += os.stat(ref_fasta).st_size
        elif os.path.isdir(ref_fasta):
            valid_ext = [".fasta", ".fa", ".fna"]
            # for multiple fasta files, bowtie2 takes a comma
            # separated list
            for in_handle in os.listdir(ref_fasta):
                fname, ext = os.path.splitext(in_handle)
                fpath = os.path.join(ref_fasta, in_handle)

                if os.path.isfile(fpath) and ext in valid_ext:
                    ref_files.append(fpath)
                    total_size += os.stat(fpath).st_size

        else:
            logger.error("The index must be either a file or directory.")

        logger.info(
            "Found %d files totaling %.3g GB.", len(ref_files), total_size / (1024 ** 3)
        )

        sequence_collection = Path(
            f"coptr-fna-{datetime.now(timezone.utc).isoformat(timespec='seconds')}.fna"
        )
        genomes = Path(f"{index_out}.genomes")

        logger.info(
            "Copying FASTA files to %s with prepended genome ids (filenames).",
            str(sequence_collection),
        )

        # assume 1 genome per fasta
        # let's set the filename as the genome identifier
        n_genomes = 0
        with sequence_collection.open("w") as out_handle, genomes.open(
            "w"
        ) as genomes_handle:
            for fpath in ref_files:
                fname = os.path.basename(os.path.splitext(fpath)[0])
                genomes_handle.write(os.path.basename(fname) + "\n")
                n_genomes += 1

                with open(fpath, "r") as in_handle:
                    for line in in_handle:
                        if line.startswith(">"):
                            # prepend the filename as the identifier for the genome
                            line = f">{fname}|{line[1:]}"
                        out_handle.write(line)

        logger.info("Writing %d reference genome ids to %s.", n_genomes, str(genomes))
        call = [
            "bowtie2-build",
            str(sequence_collection),
            index_out,
            "--threads",
            bt2_threads,
        ]
        if bt2_bmax is not None:
            call += ["--bmax", bt2_bmax]
        if bt2_dcv is not None:
            call += ["--dcv", bt2_dcv]
        if bt2_packed:
            call += ["--packed"]

        try:
            logger.info("%s", " ".join(call))
            sub.check_call(call)
            logger.info(
                "Indexed %d FASTA files for the reference database.", len(ref_files)
            )
        except Exception as error:
            logger.error("An error occurred while indexing with bowtie2.")
            logger.info("Detailed information:", exc_info=error)
        finally:
            logger.info("Cleaning up %s.", str(sequence_collection))
            sequence_collection.unlink()

    def map(self, index, inputf, outfolder, paired, threads, bt2_k):
        """Map reads from infile against reference database using bowtie2, then
        convert to a bam file.

        Parameters
        ----------
        index : str
            Path of the database index to map against
        inputf : str
            File or folder with the reads to map
        outfolder : str
            Folder to save bam files.
        paired : bool
            True for paired end sequencing.
        threads : int
            Number of threads to use. Passed to the -p argument for bowtie2.
        bt2_k : int
            Number of alignments to report
        """
        bt2_k = str(bt2_k)

        outfolder = Path(outfolder)
        if not outfolder.is_dir():
            logger.error("The output directory %s does not exist.", str(outfolder))

        if os.path.isfile(inputf):
            bn = os.path.basename(inputf)
            bn, ext = os.path.splitext(bn)
            out_sam = outfolder / f"{bn}.sam"
            out_bam = outfolder / f"{bn}.bam"

            # first map to a sam file with bowtie2
            logger.info("Mapping {} to {}".format(inputf, str(out_sam)))
            call = [
                "bowtie2",
                "-x",
                index,
                inputf,
                "--no-unal",
                "-p",
                str(threads),
                "-k",
                bt2_k,
            ]
            logger.info(" ".join(call))
            with out_sam.open("w") as out:
                sub.check_call(call, stdout=out)

            # then convert to a bam file
            logger.info("Converting {} to {}.".format(str(out_sam), str(out_bam)))
            infile = pysam.AlignmentFile(out_sam, "r")
            outfile = pysam.AlignmentFile(out_bam, "wb", template=infile)
            try:
                for s in infile:
                    outfile.write(s)
            finally:
                infile.close()
                outfile.close()

            # now remove sam file
            logger.info("Cleaning up {}.", str(out_sam))
            out_sam.unlink()

        # single end sequencing
        elif os.path.isdir(inputf) and not paired:
            valid_ext = [".fastq", ".fq", ".gz"]
            files_found = 0
            for f in sorted(os.listdir(inputf)):
                fname, ext1 = os.path.splitext(f)
                if ext1 == ".gz":
                    ext2 = fname.split(".")[1]
                else:
                    ext2 = ""

                fpath = os.path.join(inputf, f)
                if ext1 in valid_ext or ext2 in valid_ext:
                    bn = os.path.basename(f)
                    bn, ext = os.path.splitext(bn)
                    out_sam = outfolder / f"{get_fastq_name(bn)}.sam"
                    out_bam = outfolder / f"{get_fastq_name(bn)}.bam"

                    # map reads with bowtie2
                    logger.info("Mapping {} to {}".format(inputf, str(out_sam)))
                    call = [
                        "bowtie2",
                        "-x",
                        index,
                        fpath,
                        "--no-unal",
                        "-p",
                        str(threads),
                        "-k",
                        bt2_k,
                    ]
                    logger.info(" ".join(call))
                    with out_sam.open("w") as out:
                        sub.check_call(call, stdout=out)

                    # then convert to a bam file
                    logger.info(
                        "Converting {} to {}.".format(str(out_sam), str(out_bam))
                    )
                    infile = pysam.AlignmentFile(out_sam, "r")
                    outfile = pysam.AlignmentFile(out_bam, "wb", template=infile)
                    try:
                        for s in infile:
                            outfile.write(s)
                    finally:
                        infile.close()
                        outfile.close()

                    # now remove sam file
                    logger.info("Cleaning up %s.", str(out_sam))
                    out_sam.unlink()

        # paired end sequencing
        elif os.path.isdir(inputf) and paired:
            valid_ext = [".fastq", "fq", ".gz"]
            files_found = 0
            # file prefix -> [pair_1, pair_2]
            read_pairs = {}
            for f in sorted(os.listdir(inputf)):
                fname, ext1 = os.path.splitext(f)
                if ext1 == ".gz":
                    ext2 = fname.split(".")[1]
                else:
                    ext2 = ""

                # the file is not a fastq file
                if not ext1 in valid_ext and not ext2 in valid_ext:
                    continue

                # read pairs are labeled FILENAME_1.fastq or FILENAME_2.fastq
                f_split = f.split("_")
                if f_split[0] not in read_pairs:
                    read_pairs[f_split[0]] = [f]
                elif f_split[0] in read_pairs and f_split[1][0] == "2":
                    read_pairs[f_split[0]].append(f)
                elif f_split[0] in read_pairs and f_split[1][0] == "1":
                    read_pairs[f_split[0]] = [f, read_pairs[f_split[0]][0]]

            # now map paired end reads
            for pair_name in sorted(read_pairs):
                f1 = os.path.join(inputf, read_pairs[pair_name][0])
                f2 = os.path.join(inputf, read_pairs[pair_name][1])
                bn = os.path.basename(f1)
                bn = bn.split("_")[0]
                out_sam = outfolder / f"{get_fastq_name(bn)}.sam"
                out_bam = outfolder / f"{get_fastq_name(bn)}.bam"

                # map reads with bowtie2
                logger.info("Mapping {} to {}".format(inputf, str(out_sam)))
                call = [
                    "bowtie2",
                    "-x",
                    index,
                    "-1",
                    f1,
                    "-2",
                    f2,
                    "--no-unal",
                    "-p",
                    str(threads),
                    "-k",
                    bt2_k,
                ]
                logger.info(" ".join(call))
                with out_sam.open("w") as out:
                    sub.check_call(call, stdout=out)

                # then convert to a bam file
                logger.info(
                    "Converting {} to {}.".format(str(out_sam), str(out_bam))
                )
                infile = pysam.AlignmentFile(out_sam, "r")
                outfile = pysam.AlignmentFile(out_bam, "wb", template=infile)
                try:
                    for s in infile:
                        outfile.write(s)
                finally:
                    infile.close()
                    outfile.close()

                # now remove sam file
                logger.info("Cleaning up %s.", str(out_sam))
                out_sam.unlink()

        else:
            logger.error("Input must be either a file or a folder.")
