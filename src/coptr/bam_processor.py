"""
bam_processor.py
======================
Extract coordinates along each reference sequence and group reference
sequences by genome.
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

import array
import bisect
import logging
import math
import os.path
import re

import numpy as np
import pysam
from scipy.sparse import csr_matrix

from .read_assigner import ReadAssigner


logger = logging.getLogger(__name__)


class ReadContainer:
    """Container to store read positions and metadata."""

    def __init__(self):
        # read query_id -> [ids in read_data]
        self.reads = {}
        # index to ref name, read position, alignment score
        # this is more memory efficient than other options
        self.read_data = np.zeros((1000000, 3), dtype=np.uint32)
        # next row to insert data
        self.next_row_idx = 0
        # names of reference sequences
        self.ref_names_index = {}
        self.index_ref_name = {}
        self.next_ref_idx = 0

    def check_add_mapping(self, query_id, ref_name, ref_position, score):
        """Stores a mapping if it has an alignment score greater than or equal
        to the mappings seen so far.

        Parameters
        ----------
            query_id : str
                The name (identifier) of the read.
            ref_name : str
                The name of the reference sequence the read maps to.
            ref_position : int
                The location in the reference sequence the read maps to.
            score : float
                The alignment score
        """
        if query_id not in self.reads:
            self.reads[query_id] = []
            best_score = -np.inf
        else:
            row_ids = list(self.reads[query_id])
            best_score = np.max(self.read_data[row_ids, 2])

        # we've run out of room in the read data matrix,
        if self.next_row_idx == self.read_data.shape[0]:
            new_read_data = np.zeros((2 * self.next_row_idx, 3))
            new_read_data[: self.next_row_idx, :] = self.read_data
            self.read_data = new_read_data

        if ref_name not in self.ref_names_index:
            self.ref_names_index[ref_name] = self.next_ref_idx
            self.index_ref_name[self.next_ref_idx] = ref_name
            self.next_ref_idx += 1

        # don't store mappings if they are worse
        # than the one we've seen so far
        if score < best_score:
            return
        else:
            self.read_data[self.next_row_idx][0] = self.ref_names_index[ref_name]
            self.read_data[self.next_row_idx][1] = ref_position
            self.read_data[self.next_row_idx][2] = score
            self.reads[query_id].append(self.next_row_idx)
            self.next_row_idx += 1

    def get_mappings(self, query_id):
        """Return the highest scoring mappings for a read.

        Parameters
        ----------
            query_id : str
                The name of the read.

        Returns
        -------
            ref_names : list[str]
                A list of the reference sequences the read maps to.
            ref_positions : list[int]
                A list of positions in each reference sequences.
        """
        row_ids = list(self.reads[query_id])
        best_score = np.max(self.read_data[row_ids, 2])
        ref_names = []
        ref_positions = []
        for row in self.read_data[row_ids, :]:
            ref_name_idx = int(row[0])
            ref_pos = int(row[1])
            score = int(row[2])

            if score == best_score:
                ref_names.append(self.index_ref_name[ref_name_idx])
                ref_positions.append(ref_pos)

        return ref_names, ref_positions


class BamProcessor:
    """Extract coordinates from a bam file.

    Parameters
    ----------
        regex : str
            A regular expression that matches a reference genome id from
            a reference sequence name.

    Notes
    -----
        For assemblies, the reference database will contain several contigs
        per assembly, each with a reference sequence name. We need to map
        a reference sequence name to the assembly is came from. This is accomplished
        with a regular expression. The expression matches a reference genome id.

        For example, consider an assembly with contigs named

        | ``ERS235517|65|k99_317``
        | ``ERS235517|65|k99_655``
        | ``ERS235517|65|k99_1708``

        ERS235517|65 gives the reference genome, while the k99_* identify
        each contig. We want a regular expression that matches ERS235517|65.
        The regular expression

        ``\w+\|\\d+``

        matches 1 or more letter, number, or underscore (``\w+``), followed by
        the ``|`` character, followed by 1 or more number (``\d+``). This extracts the
        right reference genome id.
    """

    def __init__(self, regex="[^\|]+", is_bam=True):
        self.regex = regex
        # min percent identity to consider a valid read
        self.min_identity = 0.95

        # bam files without an index will generate a warning on
        # opening. since we don't need an index, setting the
        # verbosity will silence this message
        pysam.set_verbosity(0)

        if is_bam:
            self.read_mode = "rb"
            self.write_mode = "wb"
        else:
            self.read_mode = "r"
            self.write_mode = "w"

    def get_ref_names(self, bam_file):
        """Get the names of reference sequences and reference genomes.

        Parameters
        ----------
            bam_file : str
                Path to bam file.

        Returns
        -------
            ref_seq : set[str]
                A set of reference sequences ids
            ref_genomes : set[str]
                A set of reference genome ids
        """
        ref_seqs = set()
        ref_genomes = set()

        infile = pysam.AlignmentFile(bam_file, self.read_mode)
        for sq in infile.header["SQ"]:
            ref_seq_id = sq["SN"]

            if ref_seq_id not in ref_seqs:
                ref_seqs.add(ref_seq_id)

            # see if we can extract a reference genome id
            match = re.match(self.regex, ref_seq_id)
            if match is None:
                ref_genome = ref_seq_id
            else:
                ref_genome = match.group(0)

            if ref_genome not in ref_genomes:
                ref_genomes.add(ref_genome)

        return ref_seqs, ref_genomes

    def extract_reads(self, bam_file):
        """Get the read positions from sam_file for each reference genome.

        Parameters
        ----------
            bam_file : str
                Path to a bam file.

        Returns
        -------
            read_positions : dict[str -> list]
                A dictionary whose key is a reference sequence id,
                and the value is a list of read coordinates along that
                reference sequence
            lengths : dict[str -> int]
                A dictionary whose key is a reference sequence id,
                and the value is the length of that reference sequence
        """
        lengths = self.get_ref_seq_lengths(bam_file)

        read_container = ReadContainer()

        alignment_file = pysam.AlignmentFile(bam_file, self.read_mode)

        for aln in alignment_file.fetch(until_eof=True):

            if aln.is_unmapped or aln.is_qcfail:
                continue

            # only use one read from each mate pair
            if aln.is_read2:
                continue

            read_name = aln.query_name
            ref_name = aln.reference_name
            pos = aln.get_reference_positions()[0]

            # -6 is the default mismatch penalty for bowtie2
            # This minimum scores means that reads with perfect
            # quality scores that fall below min_identity will
            # be discarded
            min_score = -6 * (1 - self.min_identity) * aln.query_alignment_length
            alignment_score = aln.get_tag("AS")

            if alignment_score < min_score:
                continue

            assert pos < lengths[ref_name]

            read_container.check_add_mapping(read_name, ref_name, pos, alignment_score)

        alignment_file.close()
        return read_container, lengths

    def get_ref_seq_lengths(self, bam_file):
        """Extract the lengths of each reference sequence from a bam file.

        Parameters
        ----------
            bam_file : str
                Path to a bam file.

        Returns
        -------
            ref_seq_lengths : dict[str -> int]
                A dictionary whose key is a reference sequence id,
                and value is the length of that reference sequence.
        """
        ref_seq_lengths = {}
        alignment_file = pysam.AlignmentFile(bam_file, self.read_mode)
        for sq in alignment_file.header["SQ"]:

            if "LN" not in sq:
                continue
            ref_seq_lengths[sq["SN"]] = sq["LN"]

        alignment_file.close()
        return ref_seq_lengths

    def compute_bin_coverage(
        self, genome_ids, read_container, lengths, ref_seq_genome_id
    ):
        """Compute fraction of bins with reads for each reference genome.

        Parameters
        ----------
            genome_ids : set[str]
                A set of genome ids in the sample
            read_container : ReadContainer
                Reads in the sample
            lengths : dict[str -> int]
                Length of each reference sequence

        Returns
        -------
            coverage_frac : dict[str -> float]
                Fraction of nonzero read count bins for each genome_id.
        """
        binned_counts = {}
        for seq_id in ref_seq_genome_id:
            # compute 5000kb bins
            nbins = int(math.ceil(lengths[seq_id] / 5000))
            binned_counts[seq_id] = np.zeros(nbins)

        for read_id in sorted(read_container.reads):
            ref_names, ref_positions = read_container.get_mappings(read_id)

            if len(ref_names) > 1:
                continue
            bin_idx = math.floor(ref_positions[0] / 5000)
            binned_counts[ref_names[0]][bin_idx] += 1

        covered_bins = {}
        total_bins = {}
        for ref_seq in binned_counts:
            genome_id = ref_seq_genome_id[ref_seq]

            covered_bins[genome_id] = (
                covered_bins.get(genome_id, 0) + (binned_counts[ref_seq] != 0).sum()
            )
            total_bins[genome_id] = (
                total_bins.get(genome_id, 0) + binned_counts[ref_seq].size
            )

        # now add contigs without reads
        for ref_seq in ref_seq_genome_id:
            genome_id = ref_seq_genome_id[ref_seq]
            if genome_id in total_bins and ref_seq not in binned_counts:
                nbins = int(math.ceil(lengths[seq_id] / 5000))
                total_bins[genome_id] += nbins

        coverage_frac = {}
        for genome_id in genome_ids:
            if genome_id in covered_bins:
                coverage_frac[genome_id] = (
                    covered_bins[genome_id] / total_bins[genome_id]
                )
            else:
                coverage_frac[genome_id] = 0
        return coverage_frac

    def assign_multimapped_reads(self, read_container, lengths, max_alignments):
        """Assign multiple mapped reads to a single genome.

        Parameters
        ----------
            reads : dict[str] -> Read
                A dictionary whose key is the read query_id and
                value is a Read object

        Returns
        -------
            read_positions : dict[str] -> int
                A dictionary whose key is a sequence id, and
                value are the read positions along that sequence.
        """
        logger.info("Determining reference genomes.")
        genome_ids = set()

        # sequence_id -> genome_id
        ref_seq_genome_id = {}

        for ref_name in read_container.ref_names_index:
            if self.regex == "[^\|]+":
                genome_id = ref_name.split("|")[0]
            else:
                match = re.match(self.regex, ref_name)
                genome_id = match.group(0) if match is not None else ref_name
            ref_seq_genome_id[ref_name] = genome_id
            genome_ids.add(genome_id)

        genome_coverage = self.compute_bin_coverage(
            genome_ids, read_container, lengths, ref_seq_genome_id
        )

        logger.info("Collecting multi-mapped reads.")
        genome_ids = sorted(list(genome_ids))
        # reads that fail filtering criteria
        discarded_reads = set()
        # single-mapped reads can be used to set a prior
        # on the proportios of each genome in the sample
        prior_counts = np.zeros(len(genome_ids))

        # fields to construct a sparse csr_matrix
        data = []
        indices = []
        indptr = [0]
        nreads = 0
        for read_id in sorted(read_container.reads):
            ref_names, ref_positions = read_container.get_mappings(read_id)
            ref_genomes = [ref_seq_genome_id[r] for r in ref_names]

            # the read maps twice to the same genome
            if len(ref_genomes) != np.unique(ref_genomes).size:
                discarded_reads.add(read_id)
                continue

            # discard reads that map to too many genomes
            if len(ref_genomes) >= max_alignments:
                discarded_reads.add(read_id)
                continue

            # if the read maps to only one genome, use it to
            # set the prior for the read assignment model,
            # but only if the coverage is sufficient
            if len(ref_genomes) == 1 and genome_coverage[genome_id] > 0.3:
                prior_counts[bisect.bisect_left(genome_ids, genome_id)] += 1
            # otherwise, store the mappings with indicators
            elif len(ref_genomes) > 1:
                for genome_id in ref_genomes:
                    idx = bisect.bisect_left(genome_ids, genome_id)
                    indices.append(idx)
                    data.append(1)
                nreads += 1
                indptr.append(len(indices))

        if len(data) > 0:
            X = csr_matrix((data, indices, indptr), shape=(nreads, len(genome_ids)))

            logger.info("Assigning multi-mapped reads.")
            read_assigner = ReadAssigner(X, prior_counts)
            assignments = read_assigner.assign_reads()

        # sequence id -> position along sequence
        read_positions = {}
        current_row = 0
        for read_id in sorted(read_container.reads):
            ref_names, ref_positions = read_container.get_mappings(read_id)
            # if there is only one mapping
            if len(ref_names) == 1:
                ref_name = ref_names[0]
                pos = ref_positions[0]

            # skip reads that failed filtering criteria
            elif read_id in discarded_reads:
                continue

            # find the genome assignment for multi-mapped reads
            # use it to set the sequence id (either a contig or
            # full reference genome)
            else:
                mapping_genome_ids = []
                for ref_name in ref_names:
                    if ref_name not in ref_seq_genome_id:
                        match = re.match(self.regex, ref_name)
                        genome_id = match.group(0) if match is not None else ref_name
                        ref_seq_genome_id[ref_name] = genome_id
                    genome_id = ref_seq_genome_id[ref_name]
                    mapping_genome_ids.append(genome_id)

                assignment_id = assignments[current_row]
                assigned_genome = genome_ids[assignment_id]
                sequence_id = mapping_genome_ids.index(assigned_genome)
                ref_name = ref_names[sequence_id]
                pos = ref_positions[sequence_id]
                current_row += 1

            # sanity check
            assert pos < lengths[ref_name]

            if ref_name not in read_positions:
                read_positions[ref_name] = [pos]
            else:
                read_positions[ref_name].append(pos)

        return read_positions, ref_seq_genome_id

    def process_bam(self, bam_file, max_alignments=10):
        """Extract read coordinates along each reference sequence.

        Parameters
        ----------
            bam_file : str
                Path to a bam file.
            max_alignments : int
                Bounds the number alignments to consider for multi-mapped reads.
                Reads greater than or equal to this threshold are discarded.

        Returns
        -------
            coverage_maps : dict[str -> CoverageMap]
                A dictionary whose key is a reference genome id,
                and value is a CoverageMap.
        """
        logger.info("Processing %s.", bam_file)

        # extract read positions
        # reads : dict[query_seq_id] -> Read
        # lengths : dict[ref_seq_id] -> int
        reads, lengths = self.extract_reads(bam_file)

        # assign multiply mapped reads to reference sequences
        # 1. dict[ref_seq_id] -> int (read position)
        # 2. dict[ref_seq_id] -> genome_id
        read_positions, ref_seq_genome_id = self.assign_multimapped_reads(
            reads, lengths, max_alignments
        )

        # fill in read_positions with remaining sequence ids
        for seq_id in lengths:
            if seq_id not in read_positions:
                read_positions[seq_id] = []

        # we now need to group reference sequencies by species
        # assemblies will have multiple contigs, while complete reference
        # genomes will have a single contig
        contig_read_positions = {}
        contig_lengths = {}

        logger.info("Grouping reads by reference genome.")
        # group by reference genome
        # for contigs, this will group all contigs together by species
        for ref in sorted(lengths):
            if ref not in ref_seq_genome_id:
                match = re.match(self.regex, ref)

                if match is not None:
                    genome_id = match.group(0)
                else:
                    genome_id = ref
                ref_seq_genome_id[ref] = genome_id

            genome_id = ref_seq_genome_id[ref]
            if genome_id not in contig_read_positions:
                contig_read_positions[genome_id] = {}
                contig_lengths[genome_id] = {}

            contig_read_positions[genome_id][ref] = read_positions[ref]
            contig_lengths[genome_id][ref] = lengths[ref]

        # now, store coverage maps for each refence genome
        coverage_maps = {}
        for genome_id in contig_read_positions:
            is_assembly = len(contig_read_positions[genome_id]) > 1

            if is_assembly:
                coverage_maps[genome_id] = CoverageMapContig(
                    bam_file,
                    genome_id,
                    contig_read_positions[genome_id],
                    contig_lengths[genome_id],
                )
            else:
                ref_id = list(contig_read_positions[genome_id].keys())[0]
                coverage_maps[genome_id] = CoverageMapRef(
                    bam_file,
                    genome_id,
                    contig_read_positions[genome_id][ref_id],
                    contig_lengths[genome_id][ref_id],
                )

        return coverage_maps

    def merge(self, bam_files, out_bam):
        """Merge many bam files from different indexes into one, taking the
        reads with the highest mapping quality from each bam.

        Parameters
        ----------
            bam_files : list[str]
                A list of bam files to merge.
            out_bam : str
                Location to store the merged bam file.
        """
        logger.info("Merging BAM files %s.", ", ".join(bam_files))
        logger.info("Keeping reads with highest alignment score.")
        # bam header: SN => LN
        seq_len = {}
        # read_id => best_score
        score = {}
        # need to find the bam file with the highest alignment score
        for bam_file in bam_files:
            bf = pysam.AlignmentFile(bam_file, self.read_mode)

            for sq in bf.header["SQ"]:
                seq_len[sq["SN"]] = sq["LN"]

            for aln in bf.fetch(until_eof=True):
                if aln.is_unmapped or aln.is_qcfail:
                    continue

                alignment_score = aln.get_tag("AS")
                if not aln.is_read2 and aln.query_name not in score:
                    score[aln.query_name] = alignment_score
                elif not aln.is_read2 and score[aln.query_name] < alignment_score:
                    score[aln.query_name] = alignment_score

            bf.close()

        # combined header
        header = {
            "HD": {"VN": "1.0", "SO": "unsorted"},
            "SQ": [{"SN": sq, "LN": seq_len[sq]} for sq in sorted(seq_len)],
        }

        seq_names = sorted(seq_len.keys())

        logger.info("Writing merged file %s.", out_bam)
        out = pysam.AlignmentFile(out_bam, self.write_mode, header=header)
        for bam_file in bam_files:
            inf = pysam.AlignmentFile(bam_file, self.read_mode)
            for in_aln in inf.fetch(until_eof=True):
                if in_aln.is_unmapped or in_aln.is_qcfail:
                    continue

                alignment_score = in_aln.get_tag("AS")
                if (
                    in_aln.query_name in score
                    and alignment_score == score[in_aln.query_name]
                    and not in_aln.is_read2
                ):
                    # reference name will be wrong if it is not expicility set
                    out_aln = pysam.AlignedSegment()
                    out_aln.query_name = in_aln.query_name
                    out_aln.query_sequence = in_aln.query_sequence
                    out_aln.flag = in_aln.flag
                    out_aln.reference_id = bisect.bisect_left(
                        seq_names, in_aln.reference_name
                    )
                    out_aln.reference_start = in_aln.reference_start
                    out_aln.cigar = in_aln.cigar
                    out_aln.template_length = in_aln.template_length
                    out_aln.query_qualities = in_aln.query_qualities
                    out_aln.tags = in_aln.tags

                    # check that reference sequence is set correctly
                    assert (
                        seq_names[out_aln.reference_id] == in_aln.reference_name
                    ), "missing reference sequence from {}".format(bam_file)

                    out.write(out_aln)
            inf.close()
        out.close()
        logger.info("Finished writing %s.", out_bam)


class CoverageMap:
    """Data structure to store read positions along reference genomes.

    Parameters
    ----------
        bam_file : str
            Bam file used to construct coverage map
        genome_id : str
            Unique ID of the reference genome
        is_assembly : bool
            If true, the reference genome is an assembly (and has multiple
            contigs). Otherwise is a complete reference genome with a single
            contig.
    """

    __slots__ = ("bam_file", "genome_id", "sample_id", "is_assembly")

    def __init__(self, bam_file, genome_id, is_assembly):
        self.bam_file = bam_file
        self.genome_id = genome_id
        self.sample_id, _ = os.path.splitext(os.path.basename(bam_file))
        self.is_assembly = is_assembly


class CoverageMapRef(CoverageMap):
    """Data structure to store read positions from complete reference genomes.

    Parameters
    ----------
        bam_file : str
            Bam file used to construct coverage map
        genome_id : str
            Unique ID of the reference genome
        read_positions : list[int]
            A list of coordinates, one per read.
    """

    __slots__ = (
        "bam_file",
        "genome_id",
        "sample_id",
        "is_assembly",
        "read_positions",
        "length",
    )

    def __getstate__(self):
        return (
            self.bam_file,
            self.genome_id,
            self.sample_id,
            self.is_assembly,
            self.read_positions,
            self.length,
        )

    def __setstate__(self, state):
        if type(state) == list or type(state) == tuple:
            (
                self.bam_file,
                self.genome_id,
                self.sample_id,
                self.is_assembly,
                self.read_positions,
                self.length,
            ) = state
        else:
            for key in state:
                setattr(self, key, state[key])

    def __init__(self, bam_file, genome_id, read_positions, length):
        super().__init__(bam_file, genome_id, is_assembly=False)
        self.read_positions = np.array(read_positions)
        self.length = length

    def __str__(self):
        return (
            "CoverageMapRef(bam_file={}, genome_id={}, sample_id={}, reads={})".format(
                self.bam_file, self.genome_id, self.sample_id, len(self.read_positions)
            )
        )

    def __repr__(self):
        return self.__str__()

    def get_length(self):
        """Get the lenth of the genome.

        Returns
        -------
            length : int
                The length of the reference genome
        """
        return self.length

    def get_reads(self):
        """Get the read coordinates along the reference genome.

        Returns
        -------
            reads : numpy.array
                A numpy array containing the coordinates of each read
        """
        reads = np.copy(self.read_positions)
        return reads

    def count_reads(self):
        """Return number of mapped reads.

        Returns
        -------
            nreads : int
                Number of reads
        """
        return np.array(self.read_positions).size


class CoverageMapContig(CoverageMap):
    """Data structure to store read positions from assemblies.

    Parameters
    ----------
        bam_file : str
            Bam file used to construct coverage map
        genome_id : str
            Unique ID of the reference genome
        contig_read_positions : dict[str -> list]
            A dictionary whose key is a contig id (sequence id),
            and value is a list of coordinates from the reads
            along that contig.
        contig_lengths: dict[str -> int]
            A dictionary whose key is a contig id (sequence id),
            and value is the length of that contig.
    """

    __slots__ = (
        "bam_file",
        "genome_id",
        "sample_id",
        "is_assembly",
        "contig_ids",
        "contig_read_positions",
        "contig_lengths",
        "bin_size",
        "binned_reads",
        "total_bins",
        "frac_nonzero",
        "reads",
        "passed_qc_flag",
    )

    def __getstate__(self):
        return (
            self.bam_file,
            self.genome_id,
            self.sample_id,
            self.is_assembly,
            self.contig_ids,
            self.contig_read_positions,
            self.contig_lengths,
            self.bin_size,
            self.binned_reads,
            self.total_bins,
            self.frac_nonzero,
            self.reads,
            self.passed_qc_flag,
        )

    def __setstate__(self, state):
        if type(state) == list or type(state) == tuple:
            (
                self.bam_file,
                self.genome_id,
                self.sample_id,
                self.is_assembly,
                self.contig_ids,
                self.contig_read_positions,
                self.contig_lengths,
                self.bin_size,
                self.binned_reads,
                self.total_bins,
                self.frac_nonzero,
                self.reads,
                self.passed_qc_flag,
            ) = state
        else:
            for key in state:
                setattr(self, key, state[key])

        # load old pickle files
        for contig in self.contig_read_positions:
            self.contig_read_positions[contig] = array.array(
                "I", self.contig_read_positions[contig]
            )

    def __init__(self, bam_file, genome_id, contig_read_positions, contig_lengths):
        super().__init__(bam_file, genome_id, is_assembly=True)
        self.contig_ids = np.sort(list(contig_lengths.keys()))
        self.contig_lengths = contig_lengths

        self.contig_read_positions = {}
        for contig in contig_read_positions:
            self.contig_read_positions[contig] = array.array(
                "I", contig_read_positions[contig]
            )

        self.bin_size = None
        self.binned_reads = {}
        # total number of bins
        self.total_bins = None
        # fraction with nonzero read count
        self.frac_nonzero = None
        # total reads prefiltering
        self.reads = None
        # flag for qc check
        self.passed_qc_flag = None

    def __str__(self):
        return "CoverageMapContig(bam_file={}, genome_id={}, ncontigs={}, nreads={}, cov_frac={}, passed_qc={})".format(
            self.bam_file,
            self.genome_id,
            len(self.contig_ids),
            self.reads,
            self.frac_nonzero,
            self.passed_qc_flag,
        )

    def __repr__(self):
        return self.__str__()

    def reset(self):
        self.bin_size = None
        self.binned_reads = {}
        self.total_bins = None
        self.frac_nonzero = None
        self.reads = None
        self.passed_qc_flag = None

    def compute_bin_size(self):
        """Compute bin size for read counts.

        Returns
        -------
            bin_size : int
                Bin size for read counts
        """
        if self.bin_size is not None:
            return self.bin_size

        # we want approximately 500 bins
        target_bins = 500

        total_length = 0
        for contig_id in self.contig_lengths:
            length = self.contig_lengths[contig_id]
            if length >= 11000:
                total_length += length

        # bound bin size below by 1000
        bin_size = np.max((1000, total_length / target_bins))

        # want a number divisible by 100 for downstream steps
        bin_size = bin_size - (bin_size % 100)
        self.bin_size = bin_size

        return bin_size

    def get_length(self, contig_id):
        """Get the length of a contig.

        Parameters
        ----------
            contig_id : str
                Reference sequence id of the contig

        Returns
        -------
            length : int
                The length of the contig
        """
        length = self.contig_lengths[contig_id]
        return length

    def get_reads(self, contig_id):
        """Get the coordinates of all reads for a contig.

        Parameters
        ----------
            contig_id : str
                Reference sequence id of the contig

        Returns
        -------
            reads : numpy.array
                The coordinates of each read along the contig
        """
        reads = np.copy(self.contig_read_positions[contig_id])
        return reads

    def bin_reads(self, contig_id):
        """Count reads in bins.

        Parameters
        ----------
            contig_id : str
                Reference sequence id of the contig

        Returns
        -------
            binned_reads : numpy.array
                Total number of reads in each 10Kb window
        """
        bin_size = self.compute_bin_size()
        contig_length = self.contig_lengths[contig_id]

        if contig_id not in self.binned_reads:
            binned_counts = []

            nbins = int(math.floor(contig_length / bin_size))
            bin_counts = np.zeros(nbins)

            # the contig is too small
            if bin_counts.size == 0:
                return bin_counts

            if contig_id not in self.contig_read_positions:
                self.bin_counts[contig_id].append(bin_counts)
                return self.binned_reads[contig_id]

            for r in self.contig_read_positions[contig_id]:
                rbin = int(math.floor(nbins * r / contig_length))
                bin_counts[rbin] += 1

            self.binned_reads[contig_id] = bin_counts

        return self.binned_reads[contig_id]

    def count_reads(self):
        """Return number of mapped reads.

        Returns
        -------
            nreads : int
                Number of reads
        """
        nreads = 0
        for contig_id in self.contig_ids:
            nreads += np.array(self.contig_read_positions[contig_id]).size
        return nreads

    def passed_qc(self):
        """Run basic quality control checks. Sets the passed_gc_flag attribute."""
        if self.passed_qc_flag is not None:
            return self.passed_qc_flag
        total_bins = 0
        total_reads = 0
        zero_bins = 0
        for cid in self.contig_ids:
            length = self.contig_lengths[cid]
            if length < 11000:
                continue
            bins = self.bin_reads(cid)
            total_bins += bins.size
            total_reads += bins.sum()
            zero_bins += (bins == 0).sum()

        if total_bins == 0 or zero_bins / total_bins > 0.5 or total_bins < 50:
            self.passed_qc_flag = False
        else:
            self.passed_qc_flag = True

        self.frac_nonzero = (
            (total_bins - zero_bins) / total_bins if total_bins > 0 else 0
        )
        self.reads = total_reads
        self.total_bins = total_bins
        return self.passed_qc_flag
