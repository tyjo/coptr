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

import argparse
import logging
import os
import os.path
import pickle as pkl
import sys

import numpy as np

from . import __version__
from .bam_processor import BamProcessor
from .compute_read_counts import compute_read_counts
from .compute_rel_abun import compute_rel_abun
from .coptr_contig import estimate_ptrs_coptr_contig
from .coptr_ref import estimate_ptrs_coptr_ref
from .read_mapper import ReadMapper
from .util import get_fastq_name


logger = logging.getLogger(__name__)


class ProgramOptions:
    def __init__(self):
        parser = argparse.ArgumentParser(
            description="CoPTR (v{}): Compute PTRs from complete reference genomes and assemblies.".format(
                __version__
            ),
            usage="""coptr <command> [options]

command: index            create a bowtie2 index for a reference database
         map              map reads against a reference database
         merge            merge BAM files from reads mapped to multiple indexes
         extract          compute coverage maps from bam files
         estimate         estimate PTRs from coverage maps
         count            compute read counts for each genome after filtering
         rabun            estimate relative abundances for each genomes after filtering
""",
        )
        self.default_bt2_k = 10

        if len(sys.argv[1:]) < 1:
            parser.print_help()
            sys.exit(2)

        parser.add_argument("command", type=str, help="Command to run.")
        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.command):
            logger.critical("Unrecognized command.")
            parser.print_help()
            sys.exit(2)
        getattr(self, args.command)()

    def index(self):
        parser = argparse.ArgumentParser(
            usage="coptr index [-h] [--bt2-bmax BT2_BMAX] [--bt2-dcv BT2_DCV] [--bt2-threads BT2_THREADS] [--bt2-packed] ref-fasta index-out"
        )
        parser.add_argument(
            "ref_fasta",
            help="""File or folder containing fasta to index. If a folder, the extension for each
        fasta must be one of [.fasta, .fna, .fa]
        """,
        )
        parser.add_argument("index_out", help="Filepath to store index.")
        parser.add_argument(
            "--bt2-bmax",
            default=None,
            help="Set the --bmax arguement for bowtie2-build. Used to control memory useage.",
        )
        parser.add_argument(
            "--bt2-dcv",
            default=None,
            help="Set the --dcv argument for bowtie2-build. Used to control memory usage.",
        )
        parser.add_argument(
            "--bt2-threads",
            default="1",
            help="Number of threads to pass to bowtie2-build.",
        )
        parser.add_argument(
            "--bt2-packed",
            action="store_true",
            help="Set the --packed flag for bowtie2-build. Used to control memory usage.",
        )

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            sys.exit(2)

        args = parser.parse_args(sys.argv[2:])
        read_mapper = ReadMapper()
        read_mapper.index(
            args.ref_fasta,
            args.index_out,
            args.bt2_bmax,
            args.bt2_dcv,
            args.bt2_threads,
            args.bt2_packed,
        )

    def map(self):
        parser = argparse.ArgumentParser(
            usage="coptr map [-h] [--threads INT] [--bt2-k INT] [--paired] index input out-folder"
        )
        parser.add_argument("index", help="Name of database index.")
        parser.add_argument(
            "input",
            help="""File or folder containing fastq reads to map. If a folder, the extension for
        each fastq must be one of [.fastq, .fq, .fastq.gz, fq.gz]
        """,
        )
        parser.add_argument(
            "out_folder", help="Folder to save mapped reads. BAM files are output here."
        )
        parser.add_argument(
            "--paired",
            action="store_true",
            help="Set for paired end reads. Assumes fastq files end in _1.* and _2.*",
        )
        parser.add_argument(
            "--threads",
            type=int,
            default=1,
            help="Number of threads for bowtie2 mapping.",
        )
        parser.add_argument(
            "--bt2-k",
            type=int,
            default=self.default_bt2_k,
            help="(Default 10). Number of alignments to report. Passed to -k flag of bowtie2.",
        )

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            sys.exit(2)

        args = parser.parse_args(sys.argv[2:])
        read_mapper = ReadMapper()
        read_mapper.map(
            args.index,
            args.input,
            args.out_folder,
            args.paired,
            args.threads,
            args.bt2_k,
        )

    def merge(self):
        parser = argparse.ArgumentParser(
            usage="coptr merge [-h] in-bam1 in-bam2 ... in-bamN out-bam"
        )
        parser.add_argument(
            "in-bams",
            nargs="+",
            help="A space separated list of BAM files to merge. Assumes same reads were mapped against different indexes. "
            + "Only keeps read 1 of paired end sequencing, since this is used downstream.",
        )
        parser.add_argument("out-bam", help="Path to merged BAM.")

        if len(sys.argv[2:]) < 2:
            parser.print_help()
            sys.exit(2)

        args = vars(parser.parse_args(sys.argv[2:]))
        in_bams = args["in-bams"]
        out_bam = args["out-bam"]
        bam_processor = BamProcessor()
        bam_processor.merge(in_bams, out_bam)

    def extract(self):
        parser = argparse.ArgumentParser(
            usage="coptr extract [-h] [--ref-genome-regex REF_GENOME_REGEX] [--check-regex] in-folder out-folder"
        )
        parser.add_argument("in_folder", help="Folder with BAM files.")
        parser.add_argument("out_folder", help="Folder to store coverage maps.")
        parser.add_argument(
            "--ref-genome-regex",
            default="[^\|]+",
            help="Regular expression extracting a reference genome id from the sequence id in a bam file.",
        )
        parser.add_argument(
            "--check-regex",
            action="store_true",
            default=False,
            help="Check the regular expression by counting reference genomes without processing.",
        )
        parser.add_argument(
            "--bt2-k",
            type=int,
            default=self.default_bt2_k,
            help="Maximum number of alignments.",
        )

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            sys.exit(2)

        args = parser.parse_args(sys.argv[2:])

        bam_processor = BamProcessor(args.ref_genome_regex)
        ref_sequences = set()
        ref_genomes = set()
        for f in sorted(os.listdir(args.in_folder)):
            fname, ext = os.path.splitext(f)
            if ext == ".bam":
                fpath = os.path.join(args.in_folder, f)
                seq, gen = bam_processor.get_ref_names(fpath)
                ref_sequences.update(seq)
                ref_genomes.update(gen)

                if os.path.isfile(
                    os.path.join(args.out_folder, get_fastq_name(f) + ".cm.pkl")
                ):
                    logger.info("Output for %s already found, skipping.", fname)
                    continue

                # don't process the rest of the bam file if we just want to
                # sanity check the regular expression
                if args.check_regex:
                    continue

                coverage_maps = bam_processor.process_bam(fpath, args.bt2_k)
                with open(
                    os.path.join(args.out_folder, get_fastq_name(f) + ".cm.pkl"), "wb"
                ) as f:
                    pkl.dump(coverage_maps, f)

        logger.info(
            "Found %d reference sequences corresponding to %d genomes.",
            len(ref_sequences),
            len(ref_genomes),
        )
        if args.check_regex:
            logger.info("Reference genome ids:")
            for ref in sorted(ref_genomes):
                logger.info("\t%s", ref)

    def estimate(self):
        parser = argparse.ArgumentParser(
            usage="""usage: coptr estimate [-h] [--min-reads MIN_READS] [--min-cov MIN_COV] [--min-samples MIN_SAMPLES] [--threads THREADS] [--plot OUTFOLDER] [--restart] coverage-map-folder out-file
        """
        )
        parser.add_argument(
            "coverage_map_folder",
            help="Folder with coverage maps computed from 'extract'.",
        )
        parser.add_argument("out_file", help="Filename to store PTR table.")
        parser.add_argument(
            "--min-reads",
            type=float,
            help="Minimum number of reads required to compute a PTR (default 5000).",
            default=5000,
        )
        parser.add_argument(
            "--min-cov",
            type=float,
            help="Fraction of nonzero bins required to compute a PTR (default 0.75).",
            default=0.75,
        )
        parser.add_argument(
            "--min-samples",
            type=float,
            help="CoPTRContig only. Minimum number of samples required to reorder bins (default 5).",
            default=5,
        )
        parser.add_argument(
            "--plot", default=None, help="Plot model fit and save the results."
        )
        parser.add_argument(
            "--restart",
            default=False,
            action="store_true",
            help="Restarts the estimation step using the genomes in the coverage-maps-genome folder.",
        )

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            sys.exit(2)

        args = parser.parse_args(sys.argv[2:])
        sample_ids = set()
        ref_genome_ids = set()
        assembly_genome_ids = set()

        grouped_coverage_map_folder = os.path.join(
            args.coverage_map_folder, "coverage-maps-genome"
        )
        if not args.restart and not os.path.isdir(grouped_coverage_map_folder):
            os.mkdir(grouped_coverage_map_folder)

        if args.restart:
            logger.info("Restarting from files in %s.", grouped_coverage_map_folder)
            for cm_file in sorted(os.listdir(grouped_coverage_map_folder)):
                _, ext = os.path.splitext(cm_file)
                if ext != ".pkl":
                    continue

                logger.info("Checking %s.", cm_file)
                with open(
                    os.path.join(grouped_coverage_map_folder, cm_file), "rb"
                ) as f:
                    try:
                        while True:
                            cm = pkl.load(f)
                            if cm.is_assembly:
                                assembly_genome_ids.add(cm.genome_id)
                            else:
                                ref_genome_ids.add(cm.genome_id)
                            sample_ids.add(cm.sample_id)

                    except EOFError:
                        pass

        else:
            logger.info("Grouping reads by reference genome.")
            logger.info("Saving to %s:", grouped_coverage_map_folder)

            # first construct a list of genome_ids for ptr estimates
            for f in sorted(os.listdir(args.coverage_map_folder)):
                fname, ext = os.path.splitext(f)
                if ext != ".pkl":
                    continue
                fpath = os.path.join(args.coverage_map_folder, f)

                logger.info("\t%s", f)
                with open(fpath, "rb") as file:
                    coverage_maps = pkl.load(file)

                    for ref_id in coverage_maps:
                        sample_ids.add(coverage_maps[ref_id].sample_id)

                        # don't load coverage maps of species in reference database without reads
                        if coverage_maps[ref_id].count_reads() < args.min_reads:
                            continue

                        write_mode = (
                            "ab+"
                            if ref_id in ref_genome_ids or ref_id in assembly_genome_ids
                            else "wb"
                        )

                        # append to genome file
                        with open(
                            os.path.join(
                                grouped_coverage_map_folder, ref_id + ".cm.pkl"
                            ),
                            write_mode,
                        ) as tofile:
                            pkl.dump(coverage_maps[ref_id], tofile)

                        if coverage_maps[ref_id].is_assembly:
                            assembly_genome_ids.add(ref_id)
                        else:
                            ref_genome_ids.add(ref_id)

                    del coverage_maps
            logger.info("Grouping by reference genome: Complete.")
            logger.info("The --restart flag can be used to start from here.")

        sample_ids = sorted(list(sample_ids))
        results_ref = estimate_ptrs_coptr_ref(
            ref_genome_ids,
            grouped_coverage_map_folder,
            args.min_reads,
            args.min_cov,
            plot_folder=args.plot,
        )
        results_contig = estimate_ptrs_coptr_contig(
            assembly_genome_ids,
            grouped_coverage_map_folder,
            args.min_reads,
            args.min_samples,
            plot_folder=args.plot,
        )

        out_file = args.out_file
        _, ext = os.path.splitext(out_file)
        if ext != ".csv":
            out_file += ".csv"

        logger.info("Writing %s.", out_file)

        with open(out_file, "w") as f:
            # write the header
            f.write("log2(PTR):genome_id/sample_id")
            for sample_id in sample_ids:
                f.write(",{}".format(sample_id))
            f.write("\n")

            for genome_id in sorted(results_ref):
                # don't write rows without estimates
                estimates = [result.estimate for result in results_ref[genome_id]]
                if np.all(np.isnan(estimates)):
                    continue

                f.write(genome_id + ",")
                row = ["" for s in sample_ids]
                for result in results_ref[genome_id]:
                    if not np.isnan(result.estimate):
                        row[sample_ids.index(result.sample_id)] = str(result.estimate)
                f.write(",".join(row) + "\n")

            for genome_id in sorted(results_contig):
                # don't write rows without estimates
                estimates = [result.estimate for result in results_contig[genome_id]]
                if np.all(np.isnan(estimates)):
                    continue

                f.write(genome_id + ",")
                row = ["" for s in sample_ids]
                for result in results_contig[genome_id]:
                    if not np.isnan(result.estimate):
                        row[sample_ids.index(result.sample_id)] = str(result.estimate)
                f.write(",".join(row) + "\n")

        logger.info("Done.")
        logger.info("You may now remove the folder %s.", grouped_coverage_map_folder)

    def count(self):
        parser = argparse.ArgumentParser(
            usage="""usage: coptr count [-h] [--min-cov MIN_COV] [--min-samples MIN_SAMPLES] coverage-map-folder out-file
        """
        )
        parser.add_argument(
            "coverage_map_folder",
            help="Folder with coverage maps computed from 'extract'.",
        )
        parser.add_argument("out_file", help="Filename to store PTR table.")
        parser.add_argument(
            "--min-cov",
            type=float,
            help="Fraction of nonzero bins required to compute a PTR (default 0.75).",
            default=0.75,
        )
        parser.add_argument(
            "--min-samples",
            type=float,
            help="CoPTRContig only. Minimum number of samples required to reorder bins (default 5).",
            default=5,
        )

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            exit(1)

        args = parser.parse_args(sys.argv[2:])

        logger.info("Computing read counts.")

        counts, genome_ids = compute_read_counts(
            args.coverage_map_folder, args.min_cov, args.min_samples
        )

        out_file = args.out_file
        _, ext = os.path.splitext(out_file)
        if ext != ".csv":
            out_file += ".csv"

        logger.info("Writing %s.", out_file)

        with open(out_file, "w") as f:
            # write the header
            f.write("count:genome_id/sample_id")
            for sample_id in counts:
                f.write(",{}".format(sample_id))
            f.write("\n")

            for genome_id in sorted(genome_ids):
                row = [genome_id]
                for sample_id in counts:
                    if genome_id in counts[sample_id]:
                        row.append(str(counts[sample_id][genome_id]))
                    else:
                        row.append(str(0))
                row = ",".join(row) + "\n"
                f.write(row)

        logger.info("Done.")

    def rabun(self):
        parser = argparse.ArgumentParser(
            usage="""usage: coptr rabun [-h] [--min-reads MIN_READS] [--min-cov MIN_COV] [--min-samples MIN_SAMPLES] coverage-map-folder out-file
        """
        )
        parser.add_argument(
            "coverage_map_folder",
            help="Folder with coverage maps computed from 'extract'.",
        )
        parser.add_argument("out_file", help="Filename to store PTR table.")
        parser.add_argument(
            "--min-reads",
            type=float,
            help="Minimum number of reads required to compute a PTR (default 5000).",
            default=5000,
        )
        parser.add_argument(
            "--min-cov",
            type=float,
            help="Fraction of nonzero bins required to compute a PTR (default 0.75).",
            default=0.75,
        )
        parser.add_argument(
            "--min-samples",
            type=float,
            help="CoPTRContig only. Minimum number of samples required to reorder bins (default 5).",
            default=5,
        )

        if len(sys.argv[2:]) < 1:
            parser.print_help()
            exit(1)

        args = parser.parse_args(sys.argv[2:])

        logger.info("Computing relative abundances.")

        rel_abun, genome_ids = compute_rel_abun(
            args.coverage_map_folder, args.min_reads, args.min_cov, args.min_samples
        )

        out_file = args.out_file
        _, ext = os.path.splitext(out_file)
        if ext != ".csv":
            out_file += ".csv"

        logger.info("Writing %s.", out_file)

        with open(out_file, "w") as f:
            # write the header
            f.write("rel_abun:genome_id/sample_id")
            for sample_id in rel_abun:
                f.write(",{}".format(sample_id))
            f.write("\n")

            for genome_id in sorted(genome_ids):
                row = [genome_id]
                for sample_id in rel_abun:
                    if genome_id in rel_abun[sample_id]:
                        row.append(str(rel_abun[sample_id][genome_id]))
                    else:
                        row.append(str(0))
                row = ",".join(row) + "\n"
                f.write(row)

        logger.info("Done.")


def cli():
    """Serve as an entry point for command line calls."""
    logging.basicConfig(
        level="INFO",
        format="[%(levelname)s] [%(asctime)s] [%(name)s] %(message)s",
        datefmt="%b %d, %Y %H:%M:%S",
    )
    ProgramOptions()
