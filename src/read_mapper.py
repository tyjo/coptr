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

from src.print import print_error, print_info

class ReadMapper:
    """Wrapper around bowtie2.
    """

    def index(self, ref_fasta, index_out):
        """Build a bowtie2 index from ref_fasta.

        Parameters
        ----------
            ref_fasta : str
                Fasta file or folder containing fasta files to index.
                Valid extensions include '.fasta', '.fna', '.fa'
            index_out : str
                Path to output the index.
        """
        if os.path.isfile(ref_fasta):
            print_info("ReadMapper", "found 1 fasta file for reference database.")
            call = ["bowtie2-build", ref_fasta, index_out]
            print_info("ReadMapper", " ".join(call))
            sub.check_call(call)
        elif os.path.isdir(ref_fasta):
            valid_ext = [".fasta", ".fa", ".fna"]
            # for multiple fasta files, bowtie2 takes a comma
            # separated list
            ref_files = ""
            files_found = 0
            for f in os.listdir(ref_fasta):
                fname,ext = os.path.splitext(f)
                fpath = os.path.join(ref_fasta,f)
                if os.path.isfile(fpath) and ext in valid_ext:
                    ref_files += fpath + ","
                    files_found += 1

            call = ["bowtie2-build", ref_files, index_out]
            print_info("ReadMapper", " ".join(call))
            sub.check_call(call)
            print_info("ReadMapper", "indexed {} fasta files for reference database.".format(files_found))
        else:
            print_error("ReadMapper", "index must either be file or folder.")


    def map(self, index, inputf, outfolder, threads):
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

        elif os.path.isdir(inputf):
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
                        out_sam = os.path.join(outfolder, bn + ".sam")
                        out_bam = os.path.join(outfolder, bn + ".bam")

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

        else:
            print_error("ReadMapper", "input must either be a file or folder.")
