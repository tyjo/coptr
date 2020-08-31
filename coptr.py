import argparse
import numpy as np
import os
import os.path
import pickle as pkl
import sys

from src.bam_processor import BamProcessor, CoverageMapRef, CoverageMapContig
from src.coptr_contig import estimate_ptrs_coptr_contig
from src.coptr_ref import estimate_ptrs_coptr_ref
from src.print import print_error, print_info
from src.read_mapper import ReadMapper
from src.util import get_fastq_name


class ProgramOptions:

    def __init__(self):
        parser = argparse.ArgumentParser(
            description="Compute PTRs from complete reference genomes and assemblies.",
            usage='''coptr.py <command> [options]

command: index            create a bowtie2 index for a reference database
         map              map reads against a reference database
         extract          compute coverage maps from bam files
         estimate         estimate PTRs from coverage maps
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
        parser.add_argument("ref_fasta", help=
'''File or folder containing fasta to index. If a folder, the extension for each
fasta must be one of [.fasta, .fna, .fa]
'''
        )
        parser.add_argument("index_out", help="Filepath to store index.")

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
        parser.add_argument("out_folder",
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
'''coptr.py extract [-h] [--ref-genome-regex REF_GENOME_REGEX] [--check-regex]
            in-folder out-folder
'''
        )
        parser.add_argument("in_folder", help="Folder with BAM files.")
        parser.add_argument("out_folder", help="Folder to store coverage maps.")
        parser.add_argument("--ref-genome-regex", default="\w+\|\d+",
            help="Regular expression extracting a reference genome id from the sequence id in a bam file.",
        )
        parser.add_argument("--check-regex", action="store_true", default=False,
            help="Check the regular expression by counting reference genomes without processing."
        )

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

                coverage_maps = bam_processor.process_bam(fpath)
                pkl.dump(coverage_maps, open(os.path.join(args.out_folder, get_fastq_name(f) + ".cm.pkl"), "wb"))

        print_info("BamProcessor", "found {} reference sequences corresponding to {} genomes".format(len(ref_sequences), len(ref_genomes)))
        if args.check_regex:
            print_info("BamProcessor", "reference genome ids:")
            for ref in sorted(ref_genomes):
                print("\t", ref, file=sys.stderr)


    def estimate(self):
        parser = argparse.ArgumentParser(usage=
'''usage: coptr.py estimate [-h] [--min-reads MIN_READS] [--min-cov MIN_COV] [--threads THREADS] coverage-map-folder out-file
'''
        )
        parser.add_argument("coverage_map_folder", help="Folder with coverage maps computed from 'extract'.")
        parser.add_argument("out_file", help="Filename to store PTR table.")
        parser.add_argument("--min-reads", type=float, help="Minimum number of reads required to compute a PTR (default 5000).", default=5000)
        parser.add_argument("--min-cov", type=float, help="Fraction of nonzero 10Kb bins required to compute a PTR (default 0.75).", default=0.75)
        parser.add_argument("--min-samples", type=float, help="CoPTRContig only. Minimum number of samples required to reorder bins (default 5).", default=5)
        parser.add_argument("--threads", type=int, help="Number of threads to use (default 1).", default=1)

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            exit(1)

        args = parser.parse_args(sys.argv[2:])
        # reference genome id -> list of coverage maps
        coverage_maps_ref = {}
        coverage_maps_contig = {}
        sample_ids = set()

        for f in os.listdir(args.coverage_map_folder):
            fname, ext = os.path.splitext(f)
            if ext == ".pkl":
                fpath = os.path.join(args.coverage_map_folder, f)
                coverage_maps = pkl.load(open(fpath, "rb"))

                # split into maps from assemblies and compete reference genomes
                for ref_id in coverage_maps:
                    if not coverage_maps[ref_id].is_assembly and ref_id not in coverage_maps_ref:
                        coverage_maps_ref[ref_id] = [coverage_maps[ref_id]]
                    elif not coverage_maps[ref_id].is_assembly:
                        coverage_maps_ref[ref_id].append(coverage_maps[ref_id])
                    elif coverage_maps[ref_id].is_assembly and ref_id not in coverage_maps_contig:
                        coverage_maps_contig[ref_id] = [coverage_maps[ref_id]]
                    else:
                        coverage_maps_contig[ref_id].append(coverage_maps[ref_id])
                    sample_ids.add(coverage_maps[ref_id].sample_id)

        sample_ids = sorted(list(sample_ids))
        results_ref = estimate_ptrs_coptr_ref(coverage_maps_ref, args.min_reads, args.min_cov, threads=args.threads)
        results_contig = estimate_ptrs_coptr_contig(coverage_maps_contig, args.min_reads, args.min_samples, threads=args.threads)

        with open(args.out_file, "w") as f:
            # write the header
            f.write("genome_id/sample_id")
            for sample_id in sample_ids:
                f.write(",{}".format(sample_id))
            f.write("\n")

            for genome_id in sorted(results_ref):
                f.write(genome_id)
                for idx,result in enumerate(sorted(results_ref[genome_id], key=lambda x: x.sample_id)):
                    if result.sample_id != sample_ids[idx]:
                        print_error("Main", "{} is missing from {}".format(sample_id, result.ref_genome))
                    if not np.isnan(result.estimate):
                        f.write(",{}".format(result.estimate))
                    else:
                        f.write(",NA".format(result.estimate))
                f.write("\n")

            for genome_id in sorted(results_contig):
                f.write(genome_id)
                for idx,result in enumerate(sorted(results_contig[genome_id], key=lambda x: x.sample_id)):
                    if result.sample_id != sample_ids[idx]:
                        print_error("Main", "{} is missing from {}".format(sample_id, result.ref_genome))
                    if not np.isnan(result.estimate):
                        f.write(",{}".format(result.estimate))
                    else:
                        f.write(",NA".format(result.estimate))
                f.write("\n")   

if __name__ == "__main__":
    ProgramOptions()