"""
read_mapper.py
======================
Module to map reads to a reference database of
high-quality assemblies and complete reference genomes.
"""
import numpy as np
import os
import os.path
import pysam
import subprocess as sub
import sys
import time

from src.print import print_error, print_info
from src.util import get_fastq_name

class ReadMapper:
    """Wrapper around bowtie2.
    """

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
        if os.path.isfile(ref_fasta):
            ref_files = [ref_fasta]
            files_found = 1
            total_size = os.stat(ref_fasta).st_size
        elif os.path.isdir(ref_fasta):
            valid_ext = [".fasta", ".fa", ".fna"]
            # for multiple fasta files, bowtie2 takes a comma
            # separated list
            ref_files = []
            files_found = 0
            total_size = 0
            for f in os.listdir(ref_fasta):
                fname,ext = os.path.splitext(f)
                fpath = os.path.join(ref_fasta,f)

                if os.path.isfile(fpath) and ext in valid_ext:
                    ref_files.append(fpath)
                    total_size += os.stat(fpath).st_size
                    files_found += 1

        else:
            print_error("ReadMapper", "index must either be file or folder.")


        print_info("ReadMapper", "found {} files totaling {:.3f} Gb".format(len(ref_files), total_size / (1024**3)))

        fna_out = open("coptr-fna-{}.fna".format(time.time()), "w")
        genomes_out = open("{}.genomes".format(index_out), "w")
        n_genomes = 0

        print_info("ReadMapper", "copying to fasta files to {} with prepended genome ids (filenames)".format(fna_out.name))

        # assume 1 genome per fasta
        # let's set the filename as the genome identifier
        for fpath in ref_files:
            fname = os.path.basename(os.path.splitext(fpath)[0])
            genomes_out.write(os.path.basename(fname) + "\n")
            n_genomes += 1

            with open(fpath, "r") as f:
                for line in f:
                    if ">" == line[0]:
                        # prepend the filename as the identifier for the genome
                        line = line[0] + fname + "|" + line[1:]
                    fna_out.write(line)

        fna_out.close()
        genomes_out.close()

        print_info("ReadMapper", "writing {} reference genome ids to {}".format(n_genomes, genomes_out.name))
        call = ["bowtie2-build", fna_out.name, index_out, "--threads", bt2_threads]
        if bt2_bmax is not None:
            call += ["--bmax", bt2_bmax]
        if bt2_dcv is not None:
            call += ["--dcv", bt2_dcv]
        if bt2_packed:
            call += ["--packed"]

        try:
            print_info("ReadMapper", " ".join(call))
            sub.check_call(call)
            print_info("ReadMapper", "indexed {} fasta files for reference database.".format(files_found))
        except Exception as e:
            print(e, file=sys.stderr)
        finally:
            print_info("ReadMapper", "cleaning up {}".format(fna_out.name))
            sub.check_call(["rm", fna_out.name])


    def map(self, index, inputf, outfolder, paired, threads):
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
        """
        if not os.path.isdir(outfolder):
            print_error("ReadMapper", "output folder does not exist.")

        if os.path.isfile(inputf):
            bn = os.path.basename(inputf)
            bn,ext = os.path.splitext(bn)
            out_sam = os.path.join(outfolder, bn + ".sam")
            out_bam = os.path.join(outfolder, bn + ".bam")

            # first map to a sam file with bowtie2
            print_info("ReadMapper", "mapping {} to {}".format(inputf, out_sam))
            call = ["bowtie2", "-x", index, inputf, "--no-unal", "-p", str(threads)]
            print_info("ReadMapper", " ".join(call))
            sub.check_call(call, stdout=open(out_sam, "w"))

            # then convert to a bam file
            print_info("ReadMapper", "converting {} to {}".format(out_sam, out_bam))
            infile = pysam.AlignmentFile(out_sam, "r")
            outfile = pysam.AlignmentFile(out_bam, "wb", template=infile)
            for s in infile:
                outfile.write(s)

            # now remove sam file
            print_info("ReadMapper", "cleaning up {}".format(out_sam))
            call = ["rm", out_sam]
            sub.check_call(call)

        # single end sequencing
        elif os.path.isdir(inputf) and not paired:
            valid_ext = [".fastq", ".fq", ".gz"]
            files_found = 0
            for f in os.listdir(inputf):
                fname,ext1 = os.path.splitext(f)
                if ext1 == ".gz":
                    ext2 = fname.split(".")[1]
                else:
                    ext2 = ""

                fpath = os.path.join(inputf,f)
                if ext1 in valid_ext or ext2 in valid_ext:
                    bn = os.path.basename(f)
                    bn,ext = os.path.splitext(bn)
                    out_sam = os.path.join(outfolder, get_fastq_name(bn) + ".sam")
                    out_bam = os.path.join(outfolder, get_fastq_name(bn) + ".bam")

                    # map reads with bowtie2
                    print_info("ReadMapper", "mapping {} to {}".format(fpath, out_sam))
                    call = ["bowtie2", "-x", index, fpath, "--no-unal", "-p", str(threads)]
                    print_info("ReadMapper", " ".join(call))
                    sub.check_call(call, stdout=open(out_sam, "w"))

                    # then convert to a bam file
                    print_info("ReadMapper", "converting {} to {}".format(out_sam, out_bam))
                    infile = pysam.AlignmentFile(out_sam, "r")
                    outfile = pysam.AlignmentFile(out_bam, "wb", template=infile)
                    for s in infile:
                        outfile.write(s)

                    # now remove sam file
                    print_info("ReadMapper", "cleaning up {}".format(out_sam))
                    call = ["rm", out_sam]
                    sub.check_call(call)
        
        # paired end sequencing
        elif os.path.isdir(inputf) and paired:
            valid_ext = [".fastq", "fq", ".gz"]
            files_found = 0
            # file prefix -> [pair_1, pair_2]
            read_pairs = {}
            for f in os.listdir(inputf):
                fname,ext1 = os.path.splitext(f)
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
                out_sam = os.path.join(outfolder, get_fastq_name(bn) + ".sam")
                out_bam = os.path.join(outfolder, get_fastq_name(bn) + ".bam")

                # map reads with bowtie2
                print_info("ReadMapper", "mapping {},{} to {}".format(f1, f2, out_sam))
                call = ["bowtie2", "-x", index, "-1", f1, "-2", f2, "--no-unal", "-p", str(threads)]
                print_info("ReadMapper", " ".join(call))
                sub.check_call(call, stdout=open(out_sam, "w"))

                # then convert to a bam file
                print_info("ReadMapper", "converting {} to {}".format(out_sam, out_bam))
                infile = pysam.AlignmentFile(out_sam, "r")
                outfile = pysam.AlignmentFile(out_bam, "wb", template=infile)
                for s in infile:
                    outfile.write(s)

                # now remove sam file
                print_info("ReadMapper", "cleaning up {}".format(out_sam))
                call = ["rm", out_sam]
                sub.check_call(call)

        else:
            print_error("ReadMapper", "input must either be a file or folder.")
