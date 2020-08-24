import argparse
import os
import os.path
import sys

from src.bam_processor import BamProcessor, CoverageMapRef, CoverageMapContig
from src.print import print_error, print_info
from src.read_mapper import ReadMapper



class ProgramOptions:

    def __init__(self):
        parser = argparse.ArgumentParser(
            description="Compute PTRs from complete reference genomes and assemblies.",
            usage='''coptr.py <command> [options]

command: index            create a bowtie2 index for a reference database
         map              map reads against a reference database
         extract          extract read positions from bam files
'''
        )

        if len(sys.argv[1:]) < 1:
            parser.print_help()
            exit(1)

        parser.add_argument("command", type=str, help="Command to run.")
        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.command):
            print_error("Main", "Unrecognized command.", exit=False)
            parser.print_help()
        getattr(self, args.command)()


    def index(self):
        parser = argparse.ArgumentParser(usage="coptr.py index [-h] ref-fasta index-out")
        parser.add_argument("ref-fasta", help=
'''File or folder containing fasta to index. If a folder, the extension for each
fasta must be one of [.fasta, .fna, .fa]
'''
        )
        parser.add_argument("index-out", help="Filepath to store index.")

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            exit(1)

        args = parser.parse_args(sys.argv[2:])
        read_mapper = ReadMapper()
        read_mapper.index(args.ref_fasta, args.index_out)


    def map(self):
        parser = argparse.ArgumentParser(usage="coptr.py map [-h] [--threads INT] index input out-folder")
        parser.add_argument("index", help="Name of database index.")
        parser.add_argument("input", help=
'''File or folder containing fastq reads to map. If a folder, the extension for
each fastq must be one of [.fastq, .fq, .fastq.gz, fq.gz]
'''
        )
        parser.add_argument("out-folder",
            help="Folder to save mapped reads. BAM files are output here."
        )
        parser.add_argument("--threads", type=int, default=1, 
            help="Number of threads for bowtie2 mapping."
        )

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            exit(1)


        args = parser.parse_args(sys.argv[2:])
        read_mapper = ReadMapper()
        read_mapper.map(args.index, args.input, args.out_folder, args.threads)


    def extract(self):
        parser = argparse.ArgumentParser(usage=
'''usage: coptr.py extract [-h] [--ref-genome-regex REF_GENOME_REGEX] [--check-regex]
                in-folder
'''
        )
        parser.add_argument("in-folder", help="Folder with BAM files.")
        parser.add_argument("out-folder", help="Folder to store coverage maps.")
        parser.add_argument("--ref-genome-regex", default="\w+\|\d+",
            help="Regular expression extracting a reference genome id from the sequence id in a bam file.",
        )
        parser.add_argument("--check-regex", action="store_true", default=False)

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            exit(1)

        args = parser.parse_args(sys.argv[2:])

        bam_processor = BamProcessor(args.ref_genome_regex)
        ref_sequences = set()
        ref_genomes = set()
        for f in os.listdir(args.in_folder):
            fname, ext = os.path.splitext(f)
            if ext == ".bam":
                fpath = os.path.join(args.in_folder, f)
                seq, gen = bam_processor.get_ref_names(fpath)
                ref_sequences.update(seq)
                ref_genomes.update(gen)

                # don't process the rest of the bam file if we just want to 
                # sanity check the regular expression
                if args.check_regex:
                    continue

                bam_processor.process_bam(fpath)

        print_info("BamProcessor", "found {} reference sequences corresponding to {} genomes".format(len(ref_sequences), len(ref_genomes)))
        #coverage_maps = bam_processorlen




if __name__ == "__main__":
    ProgramOptions()